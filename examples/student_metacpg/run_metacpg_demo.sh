#!/usr/bin/env bash
# =============================================================================
# metaCpG-style prognostic signature evaluation with PODE: a guided demo
# -----------------------------------------------------------------------------
# Run from anywhere inside the pode repository:
#
#     bash examples/student_metacpg/run_metacpg_demo.sh
#
# It uses only the small SYNTHETIC example data shipped with the repo, so it is
# safe to run on a laptop and takes about a minute (plus a one-off install).
#
# What it does, in plain terms:
#   1. Build a "metaCpG" signature: find the DNA-methylation sites (CpGs) whose
#      levels best predict survival, and combine them into one risk score.
#   2. Check how well that score predicts survival, using leave-one-out
#      cross-validation (each patient is scored by a model that never saw them).
#   3. Turn the score into risk groups and draw the standard survival figures.
#
# Everything is written to runs/metacpg_demo/. Read the README.md next to this
# script for what each output means.
# =============================================================================
set -euo pipefail

# --- go to the repository root (works wherever the script is called from) ----
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"

OUT="$REPO/runs/metacpg_demo"
SIG_OUT="$OUT/signature"
FIGS="$OUT/figures"
mkdir -p "$OUT" "$FIGS"

say() { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }

# --- step 0: software environment (only slow the first time) -----------------
if [ ! -x .venv/bin/pode ]; then
  say "Step 0: installing PODE into .venv (first run only; takes a few minutes)"
  python3 -m venv .venv
  .venv/bin/python -m pip install --quiet --upgrade pip
  .venv/bin/python -m pip install --quiet -r requirements.txt
  .venv/bin/python -m pip install --quiet -e . --no-deps
fi
PODE="$REPO/.venv/bin/pode"
PY="$REPO/.venv/bin/python"

# --- step 1: describe the analysis in a config file --------------------------
# Every PODE analysis is driven by a YAML config. We write one here so you can
# open it afterwards and see each choice (runs/metacpg_demo/signature_config.yml).
say "Step 1: writing the signature config"
cat > "$OUT/signature_config.yml" <<EOF
schema_version: 1
methodology: legacy_exact          # the original metaCpG leave-one-out procedure
analysis_id: metacpg_demo
matrix:                            # methylation data: one row per patient, one column per CpG
  path: "$REPO/examples/signatures/methylation.csv"
  orientation: samples_by_features
  scale: beta                      # beta values (0-1); converted to M-values internally
  feature_id_column: first
phenotype:                         # clinical data: survival time + event per patient
  path: "$REPO/examples/signatures/phenotype.csv"
  sample_id_column: sample_id
endpoints:
  - name: os                       # overall survival
    time: os_time                  # years of follow-up
    status: os_event               # 1 = died, 0 = alive at last follow-up (censored)
    prediction_time: 3             # AUC is measured at 3 years
covariates:
  - column: subgroup               # allow a different baseline risk per subgroup
    type: categorical
    mode: strata
search:                            # try signatures of 3, 4, ... 15 CpGs
  min_features: 3
  max_features: 15
  step: 1
  cross_validation: loo
selection:
  mode: max_auc                    # keep the size with the best cross-validated AUC
  feature_count: null
execution:
  workers: 2
  seed: 1
validation:
  allow_unmatched_samples: true    # some clinical rows have no methylation data
output_dir: "$SIG_OUT"
EOF

# --- step 2: check the inputs before doing any modelling ---------------------
say "Step 2: validating the inputs (sample IDs match, values in range, ...)"
"$PODE" signatures validate --config "$OUT/signature_config.yml"

# --- step 3: build the signature ---------------------------------------------
# For every patient in turn: hold them out, fit one Cox model per CpG on the
# rest, rank the CpGs, build the score, then score the held-out patient.
say "Step 3: searching signature sizes and extracting the final signature"
"$PODE" signatures run --config "$OUT/signature_config.yml" \
    --backend local --engine batched --no-progress

# --- step 4: build one analysis table ----------------------------------------
# Join each patient's clinical data to two versions of the score:
#   loo_score  - from a model that never saw that patient (honest)
#   full_score - from the final model fitted on everyone (optimistic on this data)
say "Step 4: joining the scores to the clinical table"
"$PY" - "$REPO" "$SIG_OUT/python" "$OUT/cohort_with_scores.csv" <<'EOF'
import sys, pandas as pd
repo, sig, out = sys.argv[1:]
clin = pd.read_csv(f"{repo}/examples/signatures/phenotype.csv")
loo = pd.read_csv(f"{sig}/extraction_scores.csv")[["sample_id", "risk_score"]]
full = pd.read_csv(f"{sig}/risk_scores.csv")[["sample_id", "risk_score"]]
df = (clin.merge(loo.rename(columns={"risk_score": "loo_score"}), on="sample_id")
          .merge(full.rename(columns={"risk_score": "full_score"}), on="sample_id"))
df.to_csv(out, index=False)
print(f"{len(df)} patients with methylation + clinical data -> {out}")
EOF

# --- step 5: prognosis config ------------------------------------------------
say "Step 5: writing the prognosis config"
cat > "$OUT/prognosis_config.yml" <<EOF
columns:
  id: sample_id
  time: os_time
  status: os_event
  score: loo_score
  group: risk_group
  covariates: [loo_score, "strata:subgroup"]   # strata = separate baseline risk per subgroup
  categorical: [subgroup, cohort]
survival:
  horizon: 3                       # predict survival at 3 years
  time_unit: years
  landmark_offset: null
  bootstrap: {n: 200, conf: 0.95, seed: 42}
risk_groups:
  # Patients with predicted 3-year survival >= 75% are Low risk, 50-75% are
  # Intermediate and < 50% are High. These are teaching values, not clinical ones.
  thresholds: [0.75, 0.50]
  labels: ["Low", "Intermediate", "High"]
palettes:
  risk: {Low: "#2A9D8F", Intermediate: "#F4A261", High: "#E63946"}
  heatmap: ["#00008B", "#EEE9E9", "#8B0000"]
  categorical: {}
figure:
  page: {width_in: 7.09, height_in: 10.51, base_size: 8}
  device: pdf
  dpi: 300
  outdir: "$FIGS"
EOF

# --- step 6: survival model and risk groups ----------------------------------
say "Step 6: Cox model on the cross-validated score + 3-year predictions"
"$PODE" prognosis metrics --config "$OUT/prognosis_config.yml" \
    --data "$OUT/cohort_with_scores.csv" --out "$OUT/predicted_survival.csv" \
    | tee "$OUT/cox_model_summary.txt"

# --- step 7: figures ---------------------------------------------------------
say "Step 7: drawing figures"
"$PODE" prognosis km           --config "$OUT/prognosis_config.yml" --data "$OUT/cohort_with_scores.csv" \
    --out km_risk_groups.pdf --title "Overall survival by metaCpG risk group"
"$PODE" prognosis forest       --config "$OUT/prognosis_config.yml" --data "$OUT/cohort_with_scores.csv" \
    --covariates loo_score age --out forest.pdf
"$PODE" prognosis distribution --config "$OUT/prognosis_config.yml" --data "$OUT/cohort_with_scores.csv" \
    --group subgroup --kind box --out score_by_subgroup.pdf
"$PODE" prognosis calibration  --config "$OUT/prognosis_config.yml" --data "$OUT/cohort_with_scores.csv" \
    --bins 3 --out calibration.pdf

# --- summary -----------------------------------------------------------------
say "Done. Key numbers:"
"$PY" - "$SIG_OUT/python" "$OUT/cox_model_summary.txt" <<'EOF'
import sys, json, re
sig, cox = sys.argv[1:]
p = json.load(open(f"{sig}/result.json"))["performance"]
print(f"  CpGs in the signature            : {p['feature_count']}")
print(f"  Cross-validated C-index (LOO)    : {p['c_index']:.3f}   (0.5 = coin toss, 1 = perfect)")
print(f"  Cross-validated 3-year AUC (LOO) : {p['auc']:.3f}")
m = re.search(r"Concordance \(in-sample\): ([0-9.]+)", open(cox).read())
if m:
    print(f"  C-index of Cox fit on loo_score  : {m.group(1)}")
EOF
echo
echo "  Signature CpGs + weights : $SIG_OUT/python/selected_features.csv"
echo "  Search curve (pick size) : $SIG_OUT/python/search_metrics.png"
echo "  Analysis table           : $OUT/cohort_with_scores.csv"
echo "  Predicted 3-yr survival  : $OUT/predicted_survival.csv"
echo "  Figures                  : $FIGS/"
