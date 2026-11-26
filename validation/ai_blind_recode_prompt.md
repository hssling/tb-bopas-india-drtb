# AI-assisted blind recoding: model, inputs and instructions

- Model: Claude Opus 5.5 (Anthropic; model identifier claude-opus-5-5), run on 26 November 2026 as an independent agent.
- Inputs: one CSV with columns item, intervention_node, title, abstract for the 150 sampled pairs. It did not receive the coders' codes, the consensus, the classifier's codes or the scope rules.
- Output: one code per pair (the four codes below), a confidence rating and a one-sentence reason quoting the record.
- Use: pairs where its code differed from the coders' consensus (26) were returned to both coders, who re-read each record and decided the final code. It did not set any code.

## Instructions given (verbatim)

For each row decide what the record reports about THE INTERVENTION NAMED IN THAT ROW (not about anything else in the study), using exactly one of four codes:
- "favours intervention" — the study evaluates this intervention (or a direct comparison involving it) and reports a benefit for it (better accuracy, yield, outcomes, uptake, cost-effectiveness, etc.).
- "no effect or favours control" — the study evaluates this intervention and reports no benefit, a null result, or harm/worse performance.
- "unclear" — the study evaluates this intervention but the result is mixed, heterogeneous or not reported (e.g. a protocol, a description without results).
- "not an evaluation of this intervention" — the record does not evaluate this intervention at all (it merely mentions related words, studies something else, is a case report/prevalence/risk-factor/KAP study, etc.).
Rules: judge only the stated intervention; descriptive prevalence or risk-factor studies are "not an evaluation"; a diagnostic-accuracy study of the named test counts as an evaluation; a study of a different test or regimen is "not an evaluation" of the named one; if only a title is available, code from the title.
