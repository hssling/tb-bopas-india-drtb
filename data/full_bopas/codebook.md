# TB-BOPAS intervention coding: codebook for one round of AI coding

Each input line is one (record, intervention) pair: `pair_id`, `intervention` (the named intervention),
`title`, `abstract` (empty if `title_only` is true), `year`.

## Code 1: what the record reports about THE INTERVENTION NAMED IN THAT ROW (exactly one; unchanged from the validated instructions)

For each row decide what the record reports about THE INTERVENTION NAMED IN THAT ROW (not about anything else in the study), using exactly one of four codes:
- "favours intervention" — the study evaluates this intervention (or a direct comparison involving it) and reports a benefit for it (better accuracy, yield, outcomes, uptake, cost-effectiveness, etc.).
- "no effect or favours control" — the study evaluates this intervention and reports no benefit, a null result, or harm/worse performance.
- "unclear" — the study evaluates this intervention but the result is mixed, heterogeneous or not reported (e.g. a protocol, a description without results).
- "not an evaluation of this intervention" — the record does not evaluate this intervention at all (it merely mentions related words, studies something else, is a case report/prevalence/risk-factor/KAP study, etc.).
Rules: judge only the stated intervention; descriptive prevalence or risk-factor studies are "not an evaluation"; a diagnostic-accuracy study of the named test counts as an evaluation; a study of a different test or regimen is "not an evaluation" of the named one; if only a title is available, code from the title.

## Scope of each named intervention (added before the second round; apply exactly)

- Shorter all-oral regimen rollout (BPaL/BPaL-M): any shorter (6–9 month) all-oral regimen for rifampicin-resistant TB (BPaL, BPaLM, 9-month bedaquiline-containing all-oral regimens, endTB-type short regimens). Not: longer individualised regimens, drug-susceptible TB, pharmacokinetic-only or animal studies.
- Molecular diagnostic scale-up (Xpert MTB/RIF, Xpert MTB/XDR): WHO-recommended rapid molecular tests (Xpert MTB/RIF, Ultra, MTB/XDR, Truenat) and their implementation, on any specimen. Not: other assays (tNGS, LAMP, ddPCR, LPA, other PCR platforms) or studies using Xpert only as a comparator/reference.
- Whole-genome sequencing for drug-susceptibility prediction: WGS used to predict drug susceptibility. Not: targeted NGS, genomic epidemiology/transmission/lineage or mechanism studies.
- Active case finding: systematic screening of people not seeking care, including household and community contact investigation.
- TB preventive therapy: TPT itself, regimens, timing and delivery models. Not: screening/testing strategies, financial incentives, contact-tracing models.
- Airborne infection prevention and control: measures that reduce airborne TB transmission (building/room ventilation, UVGI, respirators/masks, triage, isolation). Not: mechanical ventilation of patients.
- Nutritional support: food, nutritional supplements or nutrition-linked cash for people with TB or their household contacts. Not: other populations (e.g. HIV only, malnourished children without TB).
- Private-sector engagement: engagement of private providers/pharmacies/labs in TB care (PPM, PPIA, intermediary agencies, provider incentives).
- Integrated TB-HIV-diabetes care: integrated services, screening or treatment linking TB with HIV and/or diabetes care.
- Bovine tuberculosis control: control measures for M. bovis in animals (test-and-slaughter, cattle or wildlife BCG vaccination, movement controls) and diagnostic tests used in bTB control programmes.
- Milk pasteurisation / dairy food chain: pasteurisation, heat treatment or dairy food-chain controls against M. bovis / zoonotic TB. Not: other pathogens.
- Wildlife reservoir surveillance / occupational respiratory protection: surveillance of TB in wildlife reservoirs, or respiratory protection for workers exposed to zoonotic TB.

Narrative reviews, commentaries and editorials are "not an evaluation" unless they report pooled results (then they are systematic reviews).

## Widened intervention families (agreed with the authors after the spot check, 27 Nov 2026; these override the scope rules above where they differ)

- Shorter all-oral regimen node also covers: new or repurposed drugs added to rifampicin-resistant TB regimens (bedaquiline, delamanid, pretomanid, linezolid) and tests evaluated for selecting or guiding these regimens (e.g. Xpert MTB/XDR in people with RR-TB).
- Whole-genome sequencing node also covers targeted next-generation sequencing (tNGS) for drug-susceptibility prediction.
- Programme or implementation reports that give outcome data for the named intervention count as evaluations (design "implementation/programme report (no comparator)"), coded by the direction of the reported outcomes.
- Bovine TB control also covers safety studies of BCG in cattle or wildlife, and performance studies of bTB tests (direction by the result).
- Study protocols (no results) are "not an evaluation" (tier 2), not "unclear".

## Tier 2: relevant evidence that does not evaluate the intervention

Every record coded "not an evaluation of this intervention" also gets `tier2`: true if it is on the topic of the named intervention without testing it (risk-factor or prevalence studies of the target population, transmission/genomic epidemiology relevant to it, qualitative/acceptability/formative studies, narrative reviews or commentaries about it, case reports involving it, protocols, modelling that mentions it only in passing); false if it is off-topic for that intervention (other diseases, other interventions, unrelated clinical topics, word-match only). Tier 2 never counts towards direction. For evaluations, `tier2` is false.

## Additional fields (fill for every row; for "not an evaluation" use "na" where marked)

- `design`: one of "systematic review/meta-analysis", "randomised trial", "non-randomised comparative (cohort, quasi-experimental, before-after, interrupted time series)", "diagnostic accuracy", "economic evaluation", "modelling", "implementation/programme report (no comparator)", "protocol", "other", or "na".
- `setting`: "India", "other low/middle-income", "high-income", "multi-country", "not stated", or "na".
- `outcome`: short phrase naming the main outcome used for the direction (e.g. "treatment success", "sensitivity", "case yield", "incidence"), or "na".
- `confidence`: "high", "medium" or "low".
- `quote`: the shortest exact span copied character-for-character from the title or abstract that supports the code (for "not an evaluation", the span showing what the study actually is). Do not paraphrase; it is checked automatically.

## Output

One JSON object per input row, one per line (JSONL), keys exactly:
`{"pair_id": <int>, "code": "<one of the four codes>", "design": "...", "setting": "...", "outcome": "...", "confidence": "...", "quote": "..."}`
