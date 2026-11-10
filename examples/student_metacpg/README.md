# metaCpG prognostic signature: student walkthrough

One command runs a complete metaCpG-style analysis on small **synthetic** data:

```bash
git clone https://github.com/dannlbol/pode && cd pode
bash examples/student_metacpg/run_metacpg_demo.sh
```

The first run installs the software into `.venv/` (a few minutes). After that,
the analysis takes under a minute. All outputs go to `runs/metacpg_demo/`.
You need Python 3.10–3.12 and macOS or Linux (on Windows, use WSL).

## The question

Can DNA methylation measured in a tumour predict how long a patient survives?
A **metaCpG signature** answers this by picking a handful of methylation sites
(CpGs) and combining them into one **risk score** per patient.

## What the script does

| Step | Command | In plain terms |
| --- | --- | --- |
| 1 | writes `signature_config.yml` | Says which files hold the data, which columns are survival time/event, and what to try |
| 2 | `pode signatures validate` | Checks sample IDs match, values are valid, etc. — before any modelling |
| 3 | `pode signatures run` | For each patient in turn: hide them, fit a Cox model for every CpG on the others, keep the top *N* CpGs, combine them into a score, then score the hidden patient. Repeated for *N* = 3…15 |
| 4 | small Python join | Puts each patient's clinical data and scores in one table |
| 5 | writes `prognosis_config.yml` | Defines the survival horizon (3 years) and risk-group cut points |
| 6 | `pode prognosis metrics` | Cox model of survival on the score; predicted 3-year survival; risk group per patient |
| 7 | `pode prognosis km/forest/distribution/calibration` | The standard figures |

## What to look at, in order

1. **`signature/python/search_metrics.png`** — performance against number of CpGs.
   The script keeps the size with the best AUC.
2. **`signature/python/selected_features.csv`** — the chosen CpGs. A positive weight
   means higher methylation goes with worse survival.
3. **The numbers printed at the end.** The *C-index* is the probability that, of two
   patients, the one who died first had the higher score (0.5 = guessing, 1 = perfect).
4. **`figures/km_risk_groups.pdf`** — Kaplan–Meier curves: do the Low/Intermediate/High
   groups actually survive differently?
5. **`figures/forest.pdf`**, **`calibration.pdf`**, **`score_by_subgroup.pdf`**.

## Key ideas to discuss

- **Why leave-one-out?** If you choose CpGs using all patients and then test on those
  same patients, the result looks better than it really is (overfitting). `loo_score`
  comes from models that never saw that patient. `full_score` in
  `cohort_with_scores.csv` is the final model applied to everyone. Fit a Cox model to
  each and compare the C-indices.
- **Even LOO is a little optimistic here.** The signature size was picked as the best
  of 13 tries. A truly independent validation cohort is the gold standard.
- **Censoring:** `os_event = 0` means the patient was alive at last contact, not cured.
  Cox models and Kaplan–Meier curves use these patients correctly.

## Exercises

1. In `runs/metacpg_demo/signature_config.yml`, change `selection` to
   `mode: manual` and `feature_count: 5`, then rerun just the signature step:
   `.venv/bin/pode signatures run --config runs/metacpg_demo/signature_config.yml --backend local --engine batched`.
   How do the C-index and the CpG list change?
2. Change `prediction_time` / `horizon` from 3 to 5 years.
3. Change the risk-group `thresholds` in `prognosis_config.yml` and redraw the KM plot.
   Note: editing the generated configs is fine, but rerunning the whole script
   overwrites them.

## Known issue

- `pode prognosis cutoffs` (data-driven cut points) fails with its default
  settings, so this demo uses fixed thresholds.
