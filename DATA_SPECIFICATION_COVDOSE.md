# BNT211 CLDN6 CAR-T — NONMEM Dataset Specification (Literature-Consistent: Dose as Covariate)

**File:** `nm_bnt211_carT_covdose.csv`
**Format:** comma-delimited, all-numeric, missing = `.`
**Rows:** 2,045 (158 subjects × observations only)
**Columns:** 52

---

## 1. Design Rationale

Published CAR-T cellular-kinetics models (Stein et al. 2019 CPT:PSP 8:285–295; Mueller et al. 2021 CPT 110:1362–1375 [liso-cel]; Singh et al. 2021 CPT:PSP 10:362–376 [ide-cel]) do **not** use traditional PK dosing events. Instead:

- **Initial condition** E(0) or C₀ is a model-estimated parameter, with the infused cell dose as a **covariate** that scales it.
- **There are no EVID=1 rows.** Every record is an observation (EVID=0).
- **Re-dosing** is handled through MTIME-indexed covariates (re-dose times and amounts carried on each observation row) rather than additional bolus events.
- **CARVac (RNA-LPX vaccine)** is represented entirely as time-varying covariates (cumulative count, time since last dose, recent-exposure flag) — no CMT=2 dose events.

This dataset follows those conventions. A companion event-based version (`nm_bnt211_carT.csv`) is available if traditional dosing events are preferred.

---

## 2. NONMEM Control Stream Reference

```
$INPUT
ID TIME TAD DV LNDV MDV EVID CMT AMT RATE
BLQ CENS LLOQ DVID
CARDOSE LCARDOSE DOSELVLN
NRDOSE DOSENUM CYCLE
RDTIM1 RDAMT1 RDTIM2 RDAMT2 RDTIM3 RDAMT3 RDTIM4 RDAMT4
NVAC TPV VACCFL LASTVACAMT TOTVACAMT VACSCHN
AGE SEXN RACEN WEIGHT HEIGHT BMI BSA
ECOGN TUMTYPEN LDLVLN REDLDFL REDOSEFL PROCAUTO
CNTRYN SITEN
PCELL PCELL_LLOQ NDAY

$DATA nm_bnt211_carT_covdose.csv IGNORE=@
```

### Initial Condition Parameterisation (`$DES` / `$PK`)

```
; --- $PK ---
; Initial condition: E(0) scaled by dose as covariate
TVE0  = THETA(1) * (CARDOSE / 1E8)**THETA(2)
E0    = TVE0 * EXP(ETA(1))
A_0(1) = E0            ; effector compartment seed
A_0(2) = 0             ; memory compartment starts at 0

; --- Re-dosing via MTIME (if NRDOSE > 0) ---
MTIME(1) = RDTIM1      ; time of first re-dose
MTIME(2) = RDTIM2      ; time of second re-dose (if applicable)
; etc.

; --- $DES ---
; At re-dose time, add cells to the effector compartment:
IF (MPAST(1).EQ.1 .AND. NRDOSE.GE.1) THEN
  A(1) = A(1) + FR * RDAMT1
ENDIF
IF (MPAST(2).EQ.1 .AND. NRDOSE.GE.2) THEN
  A(1) = A(1) + FR * RDAMT2
ENDIF
```

### CARVac as Covariate on Expansion Rate

```
; Option A: cumulative count
KPROL = KG * (1 + THETA(CARVac) * NVAC)

; Option B: binary recent-exposure flag
KPROL = KG * (1 + THETA(CARVac) * VACCFL)

; Option C: time since last vaccine (decaying effect)
IF (TPV.GT.0) THEN
  VAC_EFF = THETA(CARVac) * EXP(-THETA(KA_VAC) * TPV)
ELSE
  VAC_EFF = 0
ENDIF
KPROL = KG * (1 + VAC_EFF)
```

---

## 3. Column Dictionary

### Block 1 — NONMEM Structural (columns 1–14)

| # | Column | Type | Values | Description |
|---|--------|------|--------|-------------|
| 1 | ID | integer | 1–161 | Numeric subject identifier |
| 2 | TIME | real | 0 to 484 | Days from first CAR-T infusion (continuous across cycles) |
| 3 | TAD | real | 0 to 484 | Days since most recent CAR-T infusion (resets at each infusion) |
| 4 | DV | real | 14–401,540 | Transgene WPRE copies/µg DNA. BLQ records set to LLOQ (M3-ready) |
| 5 | LNDV | real | 2.6–12.9 | ln(DV), quantifiable observations only |
| 6 | MDV | flag | 0 | Always 0 (all rows are observations with measured DV) |
| 7 | EVID | flag | 0 | Always 0 (observations only — no dose events in this version) |
| 8 | CMT | integer | 1 | Always 1 (observation compartment) |
| 9 | AMT | — | . | Always missing (no dose events) |
| 10 | RATE | integer | 0 | Always 0 (placeholder) |
| 11 | BLQ | flag | 0, 1 | 1 = below LLOQ; DV set to LLOQ |
| 12 | CENS | flag | 0, 1 | M3 censoring flag (= BLQ) |
| 13 | LLOQ | real | 14, 27, 78 | Assay LLOQ at this visit |
| 14 | DVID | integer | 1 | Always 1 (transgene observation) |

### Block 2 — CAR-T Dose Covariates / Initial Condition (columns 15–17)

These parameterise E(0) in the model. The dose **does not** enter as a bolus event.

| # | Column | Type | Values | Description |
|---|--------|------|--------|-------------|
| 15 | CARDOSE | real | 1×10⁶–5×10⁸ | First CAR-T infusion: number of transduced CAR+ cells |
| 16 | LCARDOSE | real | 6.0–8.7 | log₁₀(CARDOSE) |
| 17 | DOSELVLN | real | 0–8 | Dose-level cohort (0=DL0, 1=DL1, 1.5=DL1.5, 2=DL2, 3=DL3, 7=DL7, 8=DL8) |

### Block 3 — Re-Dose Covariates / MTIME (columns 18–28)

For 24 re-dosed subjects, these carry the time and amount of each subsequent CAR-T infusion. The model uses MTIME(n)/MPAST(n) to add cells at the appropriate time within `$DES`. Missing (`.`) for subjects without re-doses at that index.

| # | Column | Type | Values | Description |
|---|--------|------|--------|-------------|
| 18 | NRDOSE | integer | 0–4 | Number of CAR-T re-doses (0 = single dose, 85% of subjects) |
| 19 | DOSENUM | integer | 1–5 | Which CAR-T infusion epoch this observation belongs to |
| 20 | CYCLE | integer | 1–5 | Dosing cycle (= DOSENUM, counts from 1) |
| 21 | RDTIM1 | real | 7–525 | TIME (days) of 1st re-dose. Missing if NRDOSE=0 |
| 22 | RDAMT1 | real | 1×10⁷–5×10⁸ | Cell dose of 1st re-dose |
| 23 | RDTIM2 | real | 86–525 | TIME of 2nd re-dose. Missing if NRDOSE<2 |
| 24 | RDAMT2 | real | 1×10⁸ | Cell dose of 2nd re-dose |
| 25 | RDTIM3 | real | — | TIME of 3rd re-dose. Missing if NRDOSE<3 |
| 26 | RDAMT3 | real | — | Cell dose of 3rd re-dose |
| 27 | RDTIM4 | real | — | TIME of 4th re-dose. Missing if NRDOSE<4 |
| 28 | RDAMT4 | real | — | Cell dose of 4th re-dose |

### Block 4 — CARVac Time-Varying Covariates (columns 29–34)

All CARVac information is carried as time-varying covariates, evaluated at each observation time. **No CMT=2 dose events exist in this dataset.**

| # | Column | Type | Values | Description |
|---|--------|------|--------|-------------|
| 29 | NVAC | integer | 0–27 | Cumulative CARVac doses received up to this observation time |
| 30 | TPV | real | 0–466 | Days since most recent CARVac dose. Missing before first vaccine |
| 31 | VACCFL | flag | 0, 1 | **Recent vaccine exposure:** 1 if a CARVac dose was given within 14 days before this observation |
| 32 | LASTVACAMT | real | 25–200 | Dose amount (µg) of the most recent CARVac. Missing before first vaccine |
| 33 | TOTVACAMT | real | 0–2,550 | Cumulative total CARVac dose (µg) to date |
| 34 | VACSCHN | integer | 1, 2 | Vaccination schedule (1 or 2). Missing for monotherapy subjects |

### Block 5 — Demographics / Body Size (columns 35–41)

| # | Column | Type | Values | Description |
|---|--------|------|--------|-------------|
| 35 | AGE | real | 18–76 | Age (years) |
| 36 | SEXN | flag | 0, 1 | 0=male, 1=female |
| 37 | RACEN | integer | 1,2,3,88,99 | 1=White, 2=Asian, 3=Black/AA, 88=Unknown, 99=Not reported |
| 38 | WEIGHT | real | 44–140 | Body weight (kg) |
| 39 | HEIGHT | real | 150–198 | Height (cm) |
| 40 | BMI | real | 17–45 | BMI (kg/m²) |
| 41 | BSA | real | 1–3 | BSA (m²) |

### Block 6 — Disease / Treatment Factors (columns 42–47)

| # | Column | Type | Values | Description |
|---|--------|------|--------|-------------|
| 42 | ECOGN | integer | 0, 1, 2 | Baseline ECOG |
| 43 | TUMTYPEN | integer | 1, 2, 3 | 1=Ovarian, 2=GCT, 3=Others |
| 44 | LDLVLN | real | 50, 75, 100 | LD intensity %. Missing for 55% (all standard LD) |
| 45 | REDLDFL | flag | 0, 1 | 1=reduced/without LD (preferred binary) |
| 46 | REDOSEFL | flag | 0, 1 | 1=received ≥1 re-dose |
| 47 | PROCAUTO | flag | 0, 1 | 1=automated manufacturing |

### Block 7 — Site / Country (columns 48–49)

| # | Column | Type | Values | Description |
|---|--------|------|--------|-------------|
| 48 | CNTRYN | integer | 1–4 | 1=Australia, 2=Germany, 3=Netherlands, 4=Sweden |
| 49 | SITEN | integer | 1–8 | Clinical site |

### Block 8 — Reference (columns 50–52)

| # | Column | Type | Values | Description |
|---|--------|------|--------|-------------|
| 50 | PCELL | real | 0–252 | Paired readout: % WPRE copies per cell (reference only) |
| 51 | PCELL_LLOQ | real | 0.005–0.05 | LLOQ for paired readout |
| 52 | NDAY | real | 0–483 | Nominal study day |

---

## 4. Differences from the Event-Based Version (`nm_bnt211_carT.csv`)

| Feature | Event-based version | This version (covariate-dose) |
|---------|--------------------|-----------------------------|
| Dose events (EVID=1) | 757 rows (191 CAR-T + 566 CARVac) | **None** |
| Total rows | 2,951 | **2,045** (observations only) |
| First CAR-T dose | EVID=1, CMT=1, AMT=cells | **CARDOSE covariate → A_0(1)** |
| Re-doses | EVID=1 at re-dose time | **RDTIM/RDAMT covariates → MTIME** |
| CARVac | EVID=1, CMT=2, AMT=µg | **NVAC, TPV, VACCFL, LASTVACAMT, TOTVACAMT** |
| Pre-dose obs (C=1) | 149 rows, IGNORE(C.EQ.1) | **Excluded** |
| Subjects | 161 | **158** (3 subjects with no post-dose observations dropped) |

---

## 5. BLQ Handling

Identical to the event-based version:
- BLQ records retained with CENS=1, DV=LLOQ (M3-ready)
- 228 / 2,045 observations (11.1%) are BLQ (lower than the 17.2% in the event-based version because pre-dose BLQ observations are excluded here)
- For M1: add `IGNORE(BLQ.EQ.1)` to `$DATA`

---

## 6. Comparison with Published Models

| Published model | DV | Dose handling | Re-dose | Vaccine/stimulation |
|-----------------|-----|---------------|---------|---------------------|
| Stein 2019 (tisa-cel) | copies/µg DNA | Covariate on Cmax/C₀ | N/A (single dose) | N/A |
| Mueller 2021 (liso-cel) | copies/µg DNA | Covariate on initial condition | N/A | N/A |
| Singh 2021 (ide-cel) | copies/µg DNA | Covariate | N/A | N/A |
| **This dataset (BNT211)** | **copies/µg DNA** | **Covariate (CARDOSE)** | **MTIME covariates** | **Time-varying covariates** |

---

## 7. QC Summary

| Check | Result |
|-------|--------|
| All columns numeric | ✓ (52 columns) |
| EVID = 0 for all rows | ✓ |
| No missing TIME | ✓ |
| DV present for all rows | ✓ (BLQ set to LLOQ) |
| Covariates constant within subject (baseline) | ✓ |
| NVAC non-decreasing within subject | ✓ |
| RDTIM values consistent with re-dose count | ✓ |
| VACCFL=1 only when CARVac given within 14d | ✓ |
| Sorted by ID, TIME | ✓ |
