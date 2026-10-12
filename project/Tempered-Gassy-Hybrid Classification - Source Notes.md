# Tempered / Gassy / Hybrid System Classification — Source Notes

Compiled 2026-09-26 for the `Runaway Reaction Venting Requirements.pptx` deck (slide 3,
"System Classification", is the fork in the road). Page numbers are **PDF page indices
(0-based)** unless marked "printed".

## Where each source sits

| Source | File | Location | What it gives you |
|---|---|---|---|
| CCPS *Engineering Design for Process Safety*, 1st ed. (1993), Ch. 14 | `CCPS Engineering Design for Process Safety (1st Ed 1993) - Ch 14 Pressure Relief Systems.pdf` (extract of book printed pp. 483–517) | §14.6.2 (extract p. ~24–25) | Four-way DIERS definitions, in plain language |
| DIERS Project Manual, Ch. VI (Bench-Scale Apparatus) | `Emergency Relief System Design Using DIERS Technology\9EB043C1…pdf` | idx 10–12, **Table VI-2** (idx 12) | Quantitative criterion (dT/dt, dP/dt at set pressure) + which sizing method is valid for each class × disengagement regime |
| DIERS Project Manual, App. VI-A | `Emergency Relief System Design Using DIERS Technology\2C4B61A8…pdf` | idx 36 (VI-A6 Fauske nomograph), idx 44–45 (VI-A7 Leung I/II/III) | Do / Do-Not rules and rationale for each method by class |
| CCPS GPREH (2nd ed.) | `Guidelines for Pressure Relief and Effluent Handling.pdf` | idx 331–332 (§3.4 "Determining System Type"); idx 366–367 (§3.5 "Selection of Applicable Balances", peak loads); idx 369–373 (§3.5.3.2 Leung Omega, tempered/gassy/hybrid) | How to *tell* from calorimetry which class you have; which balances apply |
| CCPS *Engineering Design for Process Safety*, 2nd ed. (2012) | Textbooks `Process Safety\` | — | **Not useful here** — PSV/ERS chapter was dropped; "hybrid" in 2nd ed. means dust + flammable-vapor explosion mixtures (idx 250), a different concept |
| Crosby PRV Engineering Handbook (TP-V300, 1997) | Textbooks `Process Safety\` | — | **No coverage** of reactive classification (0 hits); valve hardware/API sizing only |

## The delineation

GPREH frames it as **three** types by what carries the energy out; DIERS (and CCPS 1st ed.
Ch. 14) split it into **four** by adding the question *does venting at set pressure hold T
constant?* The four-way split is the one that actually decides the sizing method.

| DIERS class | Volatile species? | Non-condensable gas generated? | At set pressure with adequate vent | Governing balances (GPREH §3.5) |
|---|---|---|---|---|
| **Volatile / Tempered** | Yes | No | dT/dt = 0, dP/dt = 0 — latent heat removes reaction heat | Overall mass + **overall energy** balance |
| **Hybrid / Tempered** | Yes | Yes | dT/dt = 0, dP/dt = 0 — vent controls both T and P | Mass + energy + **non-condensable species** balance |
| **Hybrid / Non-tempered** | Yes, but too little / too low ΔHvap | Yes | dT/dt > 0, dP/dt > 0 — vent cannot hold T or P | Mass + energy + non-condensable (behaves gassy-like) |
| **Gassy / Non-tempered** | No (or extremely low volatility) | Yes | dT/dt > 0, dP/dt > 0 — no tempering at any point | Mass + **non-condensable species** balance only |

(Definitions: CCPS 1st ed. §14.6.2; DIERS Ch. VI §4-5 and Table VI-2 footnotes.)

Physical meaning of the fork:
- **Tempered** → P is tied to T by VLE. Stop the temperature rise and you stop the pressure
  rise; relief load is set by how much vapor must boil off to absorb q_rxn. Peak load at
  ~constant T (GPREH Eq. 3.5-19).
- **Gassy** → P is set by moles of gas generated, not T. Venting does not cool; T keeps
  climbing through the event, so reaction rate (and relief load) peak at the max-rate point.
  Relief load set by volumetric gas generation rate (Eq. 3.5-20). Requires much larger vents;
  homogeneous venting is assumed.
- **Hybrid** → both. P and T are still coupled but **P can peak at a different time than T**
  (GPREH idx 373), so you must carry both balances (Eqs. 3.5-21/22).

## How you determine the class (test data, not chemistry intuition)

From GPREH idx 331–332:
1. **Heat-up, pre-onset:** P vs T linear → gas cap only (gassy-like); ln P vs 1/T linear →
   Antoine vapor pressure (tempered-like); in between → hybrid.
2. **Cool-down after reaction (best single clue):** residual pressure well above initial →
   non-condensables formed (gassy). Pressure returns to initial at same T → tempered.
   Some residual *and* pressure still falling with T → hybrid.
3. **Vented closed cell:** if T stays constant while venting is needed to hold P → hybrid.

From DIERS App. VI-A rationale (idx 45): tempering must be proven with **open-system
testing** at the actual set pressure *and* shown to persist for the **entire** venting
incident; a system can start tempered and go non-tempered once the volatile solvent is
depleted. Also: all measured dT/dt, dP/dt must be φ-corrected to 1. Fire caveat (CCPS 1st
ed. §14.6.8.1): tempering by solvent boil-off may be invalid if the fire already boiled the
solvent away.

## Why the fork matters — method validity (DIERS Table VI-2, App. VI-A6/7)

| Method | Volatile/Tempered | Hybrid/Tempered | Hybrid/Non-tempered | Gassy/Non-tempered |
|---|---|---|---|---|
| Leung analytical | **Method I** | **Method II** | **Method III** | **Method III** |
| Fauske FAI nomograph | Preliminary only | Do not use | Do not use | **Do not use** |
| Huff pseudo-steady-state (closed ARC) | Valid | NA | NA | NA |
| Area:Charge scaling, top-vent test | Valid (vapor or homogeneous only, not churn-turbulent/bubbly) | Same | Only if test cell empties | Only if test cell empties |

**Two axes, don't mix them up (added 2026-09-28).** DIERS Leung Methods I/II/III are the
*relief-load* (W) equations sorted by class, and they are the same equations as GPREH §3.5.3.2:
I = tempered Eq. 3.5-26/27; II = hybrid Eqs. 3.5-35–37 when tempered; III = hybrid when
non-tempered, and gassy Eqs. 3.5-28–34. The *mass-flux* (G) half of the omega method is GPREH
Ch. 4 §4.2 (idx ~412–417), sorted by valve-inlet state (saturated, subcooled, with
non-condensables). Web summaries often label those "Method I/II/III" too, but that is not the
DIERS numbering. Area = W/G. In "hybrid/tempered", the first word is the class (what the
chemistry makes) and the second is the result of the venting test, which only splits hybrid.

Leung I is "only valid if dT/dt ≈ 0 at the ERS design pressure"; Leung II requires dT/dt ≈ 0
**and** dP/dt ≈ 0; Leung III needs combined closed- and open-cell data (idx 45). Caution from
DIERS: if two-phase area comes out *smaller* than all-vapor area, the method is misapplied.

## Relevance to the EO case

KOH-catalyzed ethoxylation: EO is volatile (tempering possible) but EO decomposition and
polyol degradation can generate non-condensables → expect **hybrid**, and whether it stays
*tempered* depends on free-EO inventory relative to the reaction heat — i.e., the
loss-of-mixing EO pocket size can push the same reactor from hybrid/tempered to
hybrid/non-tempered. See `EO Case Study\` papers.
