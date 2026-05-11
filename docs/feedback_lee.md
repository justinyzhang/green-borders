# Defense Feedback: Tom Lee (April 23, 2025)

Six items raised during the v4 defense and subsequent written feedback. These guide Chapter 5 redesign.

## 1. Causal identification claim is weak

The OLS density premium (β = 2.895, p = 0.003) is a cross-section. It does not identify a causal effect of IL legalization. Lee suggested:

- Use Missouri's Feb 2023 legalization as a second event for staggered DiD.
- Examine receiving-side outcomes (crime, DUI, traffic accidents) with proper pre/post timing.

## 2. Border definition is binary

"Border" is treated as a 0/1 indicator. Consider:

- Distance-based intensity (continuous miles to nearest IL border).
- Population-weighted exposure (cluster 11 interaction term).

## 3. Out-of-state customer share is asserted, not estimated

The 22-40% out-of-state revenue figure comes from secondary sources (IDOR press releases). Lee asked whether this can be triangulated from:

- Sales tax receipts disaggregated by point-of-sale county.
- Dispensary visit data from cell-phone-mobility datasets (e.g., SafeGraph, Advan).

## 4. Population controls

The cross-section uses county population as a control, but does not separately control for the population of the *adjacent prohibition state* counties. This is the basis for cluster 11.

## 5. Time-varying outcome data is absent

The v4 slides have one outcome (entry count) measured at one point in time. Need panel outcomes with pre/post variation. Chapter 5 (NIBRS DiD) addresses this.

## 6. Selection into observed dispensaries

Only operational dispensaries are observed. Selection into *applying* for and *receiving* a license is unmodeled. Could bias the density premium downward if rejected applications cluster near borders (because regulators see the cross-border demand as a fiscal opportunity for IL but a political problem).
