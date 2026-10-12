"""
EO / propylene glycol runaway relief model (stage 2).

Scenario (user decision): the runaway is the NORMAL KOH-catalyzed ethoxylation accelerating because
unreacted EO has built up in the liquid. Start from normal operating T and P mid-batch, raise the
dissolved-EO mole fraction to x_acc, keep the normal cooling duty (or lose it), and vent
(homogeneous disengagement) so the vessel pressure is held at set*(1+accumulation).
The vent mass rate required to hold that pressure is the relieving rate W(t).

Sources (Reference Materials\\...):
  [S96]  Di Serio et al., IECR 1996, 35, 3848  -- k11/k21 Arrhenius, EO Antoine, EO density
  [Z26]  Zahnd et al., OPRD 2026, 30, 628      -- dHrxn, EO liquid cp (eq 5), VLE caveats
  [C59]  Crocco et al., J Chem Phys 1959, 31, 506 -- gas-phase EO decomposition (first order)
  [thesis] Chiu, Swinburne                     -- 6 m3 reactor, 3.62 bar normal, 4.45 bar design
PLACEHOLDERS (replace with Perry's estimates): polyether cp, polyether density, EO latent heat.
Homogeneous venting = every species leaves in proportion to the vessel inventory (CCPS 3.2.8.2).
"""
import numpy as np
from dataclasses import dataclass
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ============================================================================
# Constants
# ============================================================================
R = 8.314462                 # J/mol/K
R_CAL = 1.98720              # cal/mol/K
MW_EO = 44.053e-3            # kg/mol
MW_PG = 76.095e-3            # kg/mol
MW_KOH = 56.106e-3           # kg/mol
MW_N2 = 28.014e-3            # kg/mol
MW_NC = 22.03e-3             # kg/mol, mean of CH4 and CO (decomposition gas)
TC_EO = 469.0                # K   (195.8 C, [Z26])
TB_EO = 283.85               # K
LV_EO_TB = 25.5e3            # J/mol at Tb -- PLACEHOLDER (from memory; verify in Perry's)
DH_RXN = -92.0e3             # J/mol EO   [Z26] lit.; flow calorimetry 88 +/- 4
DH_DEC = -132.0e3            # J/mol EO -> CH4 + CO (own calc from standard dHf; verify)
CP_POLY = 2400.0             # J/kg/K polyether -- PLACEHOLDER (Perry's)
ATM = 101325.0

# Kinetics [S96] Table 5: k = exp(lnA) * exp(-E/(R_cal T)); A in cm3/mol/s, E in cal/mol
K11 = (20.52, 13.0e3)        # EO + primary-OH alkoxide (chain propagation on PG adduct)
K21 = (23.83, 16.8e3)        # EO + secondary-OH alkoxide (PG initiation)
# Decomposition [C59]: k = 1e11 exp(-42000/(R_cal T)) 1/s (gas phase; data 900-1200 K only)
KD = (1.0e11, 42.0e3)


# ============================================================================
# Scenario definition
# ============================================================================
@dataclass
class Scenario:
    V_vessel: float = 6.0        # m3
    m_PG: float = 1000.0         # kg starter charge
    avg_eo_units: float = 3.0    # mean EO units per PG chain already reacted (mid-batch)
    koh_wt: float = 0.003        # kg KOH per kg PG
    T_op: float = 413.15         # K   normal operating temperature (120-180 C, [thesis])
    P_op: float = 3.62e5         # Pa  abs, normal operating pressure [thesis]
    x_ss: float = 0.03           # baseline dissolved-EO mole fraction in normal operation (ASSUMED)
    x_acc: float = 0.05          # accumulated dissolved-EO mole fraction (the swept variable)
    P_set: float = 4.45e5        # Pa  abs relief set = design pressure [thesis]
    acc: float = 0.10            # accumulation fraction (relieving P = P_set*(1+acc))
    gamma_EO: float = 1.0        # Raoult; unvalidated at high x_EO [Z26]
    k_choice: str = "k11"
    k_scale: float = 1.0         # [Z26] found x1.5 needed on literature kinetics
    cooling: str = "normal"      # "normal": UA sized to remove baseline heat; "lost": none
    dT_cool: float = 30.0        # K   T_op - coolant T at normal duty
    phi_vessel_Cp: float = 0.0   # J/K vessel steel thermal mass
    decomp: bool = True
    vent: bool = True
    tau_p: float = 2.0           # s   pressure-hold relaxation


# ============================================================================
# Properties
# ============================================================================
def psat_EO(T):
    """Pa. Antoine [S96] eq 5. Held at its value at Tc above 469 K (diagnostic only beyond)."""
    T = min(max(T, 250.0), TC_EO)
    return np.exp(16.74 - 2568.0 / (T - 29.01)) * 133.322

def rho_EO(T):
    """kg/m3, linear fit to Yen-Woods values in [S96] Table 3."""
    T = min(max(T, 273.15), TC_EO - 5.0)
    return 740.0 - 2.333 * (T - 273.15 - 100.0)

def rho_poly(T):
    """kg/m3 -- PLACEHOLDER for PG/EO adduct (Perry's)."""
    return 1000.0 - 0.7 * (T - 293.15)

def cp_EO_L(T):
    """J/kg/K, [Z26] eq 5 (fitted -40..150 C; OCR'd coefficients -- verify)."""
    T = min(max(T, 233.0), 423.0)
    return 0.11017 * T**2 - 60.153 * T + 10055.2936

def Lv_EO(T):
    """J/mol, Watson scaling from Tb (PLACEHOLDER magnitude)."""
    Tr, Trb = min(T / TC_EO, 0.999), TB_EO / TC_EO
    return LV_EO_TB * ((1 - Tr) / (1 - Trb)) ** 0.38

def k_arr(p, T):
    """m3/mol/s from (lnA [cm3/mol/s], E [cal/mol])."""
    return np.exp(p[0]) * np.exp(-p[1] / (R_CAL * T)) * 1e-6


# ============================================================================
# Vapor-liquid split: inventory z = [n_PG, n_rxn, n_free, n_N2, n_NC, n_KOH]
# ============================================================================
def flash(s: Scenario, T, z):
    n_PG, n_rxn, n_free, n_N2, n_NC = z[0], z[1], z[2], z[3], z[4]
    n_tot = max(n_free, 0.0)
    V_poly = (n_PG * MW_PG + n_rxn * MW_EO) / rho_poly(T)
    Ps = psat_EO(T)
    rE = rho_EO(T)

    def resid(nL):
        V_L = V_poly + nL * MW_EO / rE
        V_g = max(s.V_vessel - V_L, 1e-6)
        x = nL / (nL + n_PG)
        return s.gamma_EO * x * Ps * V_g / (R * T) + nL - n_tot

    if n_tot <= 0.0:
        nL = 0.0
    elif resid(n_tot) <= 0.0:
        nL = n_tot
    else:
        nL = brentq(resid, 0.0, n_tot, xtol=1e-12, rtol=1e-12)
    V_L = V_poly + nL * MW_EO / rE
    V_g = max(s.V_vessel - V_L, 1e-6)
    x = nL / (nL + n_PG)
    P_EO = s.gamma_EO * x * Ps
    P = P_EO + (n_N2 + n_NC) * R * T / V_g
    return dict(n_L=nL, n_V=n_tot - nL, V_L=V_L, V_g=V_g, P=P, P_EO=P_EO, x=x)


def init_state(s: Scenario):
    """Normal operating inventory (x_ss at P_op) then accumulate EO to x_acc at constant T."""
    T = s.T_op
    n_PG = s.m_PG / MW_PG
    n_rxn = s.avg_eo_units * n_PG
    n_KOH = s.koh_wt * s.m_PG / MW_KOH
    V_poly = (n_PG * MW_PG + n_rxn * MW_EO) / rho_poly(T)
    Ps = psat_EO(T)

    def n_free_for(x):
        nL = x / (1 - x) * n_PG
        V_g = s.V_vessel - V_poly - nL * MW_EO / rho_EO(T)
        P_EO = s.gamma_EO * x * Ps
        return nL + P_EO * V_g / (R * T), P_EO, V_g

    n_ss, P_EO_ss, V_g_ss = n_free_for(s.x_ss)
    n_N2 = (s.P_op - P_EO_ss) * V_g_ss / (R * T)
    if n_N2 <= 0:
        raise ValueError("Baseline EO partial pressure exceeds P_op; lower x_ss")
    n_acc, _, _ = n_free_for(s.x_acc)
    z_ss = [n_PG, n_rxn, n_ss, n_N2, 0.0, n_KOH]
    z0 = [n_PG, n_rxn, n_acc, n_N2, 0.0, n_KOH]
    return z_ss, z0


# ============================================================================
# Dynamics
# ============================================================================
def _rate_rx(s, T, f, n_KOH):
    kp = K11 if s.k_choice == "k11" else K21
    return k_arr(kp, T) * s.k_scale * f["n_L"] * n_KOH / f["V_L"]        # mol/s

def _eval(s: Scenario, y, cool):
    """Return (dy/dt, info). y = [T, z(6), n_dec_cum]."""
    T = float(np.clip(y[0], 250.0, 1500.0))
    z = np.maximum(np.array(y[1:7], dtype=float), 0.0)
    f0 = flash(s, T, z)
    r_rx = _rate_rx(s, T, f0, z[5])
    r_dec = KD[0] * np.exp(-KD[1] / (R_CAL * T)) * f0["n_V"] if s.decomp else 0.0
    src = np.array([0.0, r_rx, -r_rx - r_dec, 0.0, 2.0 * r_dec, 0.0])

    C = (z[0] * MW_PG + z[1] * MW_EO) * CP_POLY + f0["n_L"] * MW_EO * cp_EO_L(T) + s.phi_vessel_Cp
    Q_cool = cool["UA"] * (T - cool["T_cw"]) if s.cooling == "normal" else 0.0
    Q = -DH_RXN * r_rx - DH_DEC * r_dec - Q_cool
    Lv = Lv_EO(T)

    # Sensitivities of n_V and P to T and to each inventory component (finite differences)
    dT = 0.01
    fp, fm = flash(s, T + dT, z), flash(s, T - dT, z)
    nV_T = (fp["n_V"] - fm["n_V"]) / (2 * dT)
    P_T = (fp["P"] - fm["P"]) / (2 * dT)
    nV_j = np.zeros(6); P_j = np.zeros(6)
    for j in range(5):                       # KOH (j=5) does not enter the VLE
        eps = max(1e-6 * abs(z[j]), 1e-3)
        zz = z.copy(); zz[j] += eps
        fj = flash(s, T, zz)
        nV_j[j] = (fj["n_V"] - f0["n_V"]) / eps
        P_j[j] = (fj["P"] - f0["P"]) / eps

    D = C + Lv * nV_T
    A = Q - Lv * float(nV_j @ src)
    coefT = (Lv * (float(nV_j @ z) - f0["n_V"]) - f0["P"] * s.V_vessel) / D   # per unit omega
    g0 = P_T * (A / D) + float(P_j @ src)
    g1 = P_T * coefT - float(P_j @ z)

    P_hold = s.P_set * (1.0 + s.acc)
    omega = 0.0
    if s.vent and g1 < 0.0:
        w = min(max((f0["P"] - 0.995 * P_hold) / (0.005 * P_hold), 0.0), 1.0)
        omega = w * max(0.0, (g0 + max(f0["P"] - P_hold, 0.0) / s.tau_p) / (-g1))
    dTdt = A / D + omega * coefT
    dz = src - omega * z
    M_tot = (z[0] * MW_PG + z[1] * MW_EO + z[2] * MW_EO + z[3] * MW_N2 + z[4] * MW_NC + z[5] * MW_KOH)
    info = dict(P=f0["P"], omega=omega, W=omega * M_tot, r_rx=r_rx, r_dec=r_dec, dTdt=dTdt,
                x=f0["x"], n_L=f0["n_L"], n_V=f0["n_V"], Q_cool=Q_cool, M=M_tot)
    return np.concatenate([[dTdt], dz, [r_dec]]), info

def _cool_params(s, z_ss):
    f = flash(s, s.T_op, np.array(z_ss))
    Q_norm = -DH_RXN * _rate_rx(s, s.T_op, f, z_ss[5])
    return dict(UA=Q_norm / s.dT_cool, T_cw=s.T_op - s.dT_cool, Q_norm=Q_norm)


def run(s: Scenario, t_end=3600.0, T_stop=470.0, n_out=400):
    z_ss, z0 = init_state(s)
    cool = _cool_params(s, z_ss)
    y0 = np.array([s.T_op] + z0 + [0.0])
    f_init = flash(s, s.T_op, np.array(z0))

    def rhs(t, y): return _eval(s, y, cool)[0]
    def ev_done(t, y): return y[3] - 1e-3 * z0[2]
    ev_done.terminal = True
    def ev_T(t, y): return y[0] - T_stop
    ev_T.terminal = True

    atol = np.array([1e-4, 1e-2, 1e-2, 1e-2, 1e-4, 1e-3, 1e-3, 1e-4])
    sol = solve_ivp(rhs, (0.0, t_end), y0, method="LSODA", rtol=1e-6, atol=atol,
                    events=[ev_done, ev_T], dense_output=True, max_step=t_end / 300)
    t = np.linspace(0.0, sol.t[-1], n_out)
    Y = sol.sol(t)
    infos = [_eval(s, Y[:, i], cool)[1] for i in range(t.size)]
    g = lambda k: np.array([d[k] for d in infos])
    W, P = g("W"), g("P")
    i_pk = int(np.argmax(W))
    ended = "EO consumed" if sol.t_events[0].size else ("T_stop (>Tc: diagnostic)" if sol.t_events[1].size else "t_end")
    return dict(t=t, T=Y[0], P=P, W=W, r_rx=g("r_rx"), x=g("x"), dTdt=g("dTdt"), n_dec=Y[7],
                P0=f_init["P"], W_peak=W[i_pk], t_peak=t[i_pk], T_peak=Y[0][i_pk],
                T_max=Y[0].max(), P_max=P.max(), dec_frac=Y[7][-1] / z0[2], ended=ended,
                Q_norm=cool["Q_norm"], x0=f_init["x"], z_ss=z_ss, z0=z0,
                eo_excess_kg=(z0[2] - z_ss[2]) * MW_EO)


# ============================================================================
# Sweep of accumulated EO
# ============================================================================
def run_init_pressure(s):
    _, z0 = init_state(s)
    return flash(s, s.T_op, np.array(z0))["P"]

def x_max_below_set(s: Scenario):
    """Largest x_acc whose isothermal starting pressure stays below the set pressure."""
    def f(x):
        return run_init_pressure(Scenario(**{**s.__dict__, "x_acc": x})) - s.P_set
    return brentq(f, s.x_ss, 0.5)


if __name__ == "__main__":
    import sys
    import matplotlib
    if "--save" in sys.argv:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    base = Scenario()
    xmax = x_max_below_set(base)
    print(f"T_op={base.T_op-273.15:.0f} C  P_op={base.P_op/1e5:.2f} bar  P_set={base.P_set/1e5:.2f} bar  "
          f"P_relieve={base.P_set*(1+base.acc)/1e5:.2f} bar  x_ss={base.x_ss}  x_max(P0<P_set)={xmax:.3f}")
    xs = [0.03, 0.04, 0.05, 0.06, round(xmax, 3), 0.08, 0.10, 0.15]
    print(f"{'x_acc':>6} {'coolg':>6} | {'P0 bar':>7} {'EO xs kg':>9} | {'Wpeak kg/s':>10} {'t_pk s':>7} {'T_pk C':>7} "
          f"{'Pmax bar':>8} {'decomp%':>8}  end")
    results = {}
    for cool in ("normal", "lost"):
        for x in xs:
            s = Scenario(x_acc=x, cooling=cool)
            try:
                o = run(s)
            except Exception as e:
                print(f"{x:6.3f} {cool:>6} | failed: {e}"); continue
            results[(cool, x)] = o
            print(f"{x:6.3f} {cool:>6} | {o['P0']/1e5:7.2f} {o['eo_excess_kg']:9.0f} | {o['W_peak']:10.2f} "
                  f"{o['t_peak']:7.0f} {o['T_peak']-273.15:7.1f} {o['P_max']/1e5:8.2f} {o['dec_frac']*100:8.4f}  {o['ended']}")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for (cool, x), o in results.items():
        if cool == "lost":
            ax[0].plot(o["t"], o["W"], label=f"x={x:.3f}")
    ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("required relieving rate W (kg/s)"); ax[0].legend(); ax[0].set_title("cooling lost")
    for cool, m in (("normal", "o-"), ("lost", "s--")):
        xx = [x for (c, x) in results if c == cool]
        ax[1].plot(xx, [results[(cool, x)]["W_peak"] for x in xx], m, label=f"cooling {cool}")
    ax[1].axvline(xmax, c="r", ls=":"); ax[1].set_xlabel("accumulated dissolved-EO mole fraction")
    ax[1].set_ylabel("peak relieving rate (kg/s)"); ax[1].legend()
    fig.tight_layout()
    if "--save" in sys.argv:
        fig.savefig("eo_relief_vs_accumulation.png", dpi=120)
    else:
        plt.show()
