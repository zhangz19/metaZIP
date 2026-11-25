# Real-example data

These CSV files contain observed aggregate counts. The matching R data objects
are available through `data(eus_fna)` and `data(rosiglitazone_mi)`.
No counts were simulated or continuity-corrected.

## EUS-FNA: eus_fna.csv

Source: Chawla A, Maiti T, Pradhan V, Zhang Z, *Meta-analysis of rare events with
application to medical data*, author-supplied manuscript, Table 1
(`ZIPR26.tex`; also present in `ZIPR21.pdf`). Values were parsed directly from
the LaTeX table, reading each column block top to bottom, then left to right.

There are 51 single-arm studies. The six 19G studies contribute 3 events among
215 participants; the 45 22/25G studies contribute 74 events among 5,957
participants. Four and 27 studies, respectively, have zero events.

Columns are `study_id`, `needle`, `x`, `y`, and `ni`. The indicator `x` is 1 for
19G and 0 for 22/25G; `y` and `ni` are the event count and sample size.
`study_id` is an added sequential display identifier, not an original study
identifier. The source does not provide study names. Rows are independent
single-arm study records and must not be paired across groups by row number.

## Rosiglitazone MI: rosiglitazone_mi.csv

Source: Cai T, Parast L, Ryan L (2010). Meta-analysis for rare events.
*Statistics in Medicine*, 29(20), 2078-2089. **Table I, page 2083**.
[DOI: 10.1002/sim.3964](https://doi.org/10.1002/sim.3964).
[PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC2932857/).
The complete table was verified against a
[full-text copy of the original article](https://www.academia.edu/143468083/Meta_analysis_for_rare_events),
because direct PMC retrieval was blocked in the preparation environment.

This dataset contains the **48 paired trials** in that table. It retains the
source study numbers 1-48, including the six additional studies omitted from
the 42-trial Nissen and Wolski (2007) analysis because both MI and cardiovascular
death were absent in both arms. It is distinct from the later **56-trial**
update.

Columns are `study_id`, `y1`, `n1`, `y0`, and `n0`. The first arm is
rosiglitazone; the second is control, which is not uniformly placebo.
**The endpoint is myocardial infarction, not cardiovascular death.** Death
counts in the original table are not included in this file.

Checks against the source:

- 48 consecutive source study IDs, with no missing or added rows.
- Rosiglitazone: 86 MI events among 16,856 participants.
- Control: 72 MI events among 12,962 participants.
- 10 trials with zero MI events in both arms.
- All five output fields in all 48 rows compared programmatically against the
  retrieved source table, with no differences.
- Nonnegative integer counts, each no larger than its arm sample size.

Original clinical source: Nissen SE, Wolski K (2007). Effect of rosiglitazone
on the risk of myocardial infarction and death from cardiovascular causes.
*New England Journal of Medicine*, 356(24), 2457-2471.
[DOI: 10.1056/NEJMoa072761](https://doi.org/10.1056/NEJMoa072761).

## Attribution

The repository's MIT software license retains Copyright (c) 2022 Zhen Zhang.
That license does not relicense the cited articles. The external article is
Copyright (c) 2010 John Wiley & Sons, Ltd. These CSV files contain factual
aggregate counts with source attribution; article prose, figures, and formatted
tables are not bundled.
