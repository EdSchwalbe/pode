import numpy as np
import prognostk as pk
from conftest import make_toy, toy_cfg


def test_fit_cox_protective_score():
    fit = pk.fit_cox(make_toy(), toy_cfg())
    assert fit.fitter.params_["metaCpG"] < 0


def test_predict_valid_probabilities():
    cfg = toy_cfg(); d = make_toy()
    ps = pk.predict_survival_at_time(pk.fit_cox(d, cfg), d, cfg=cfg)
    assert len(ps) == len(d)
    assert np.all((ps >= 0) & (ps <= 1))


def test_stratified_predict():
    cfg = toy_cfg(); d = make_toy()
    fit = pk.fit_cox(d, cfg, covariates=["metaCpG", "strata:grp"])
    ps = pk.predict_survival_at_time(fit, d, cfg=cfg)
    assert np.all(np.isfinite(ps))


def test_stratified_predict_keeps_row_order():
    # lifelines groups stratified predictions by stratum; each patient must still
    # get their own prediction. Within a stratum survival is then a monotone
    # function of the score (protective here, so rank correlation is +1).
    cfg = toy_cfg(); d = make_toy()
    fit = pk.fit_cox(d, cfg, covariates=["metaCpG", "strata:grp"])
    ps = pk.predict_survival_at_time(fit, d, cfg=cfg)
    for _, g in d.assign(ps=ps).groupby("grp"):
        assert g["metaCpG"].corr(g["ps"], method="spearman") > 0.999


def test_stratified_predict_matches_single_rows():
    # Predicting the cohort at once must equal predicting each patient alone,
    # including with a shuffled, non-default index.
    cfg = toy_cfg(); d = make_toy()
    fit = pk.fit_cox(d, cfg, covariates=["metaCpG", "strata:grp"])
    shuffled = d.sample(frac=1, random_state=3).set_index(np.arange(len(d)) + 500)
    ps = pk.predict_survival_at_time(fit, shuffled, cfg=cfg)
    one_by_one = [pk.predict_survival_at_time(fit, shuffled.iloc[[i]], cfg=cfg)[0]
                  for i in range(len(shuffled))]
    assert np.allclose(ps, one_by_one)


def test_assign_risk_groups_order():
    cfg = toy_cfg()
    g = pk.assign_risk_groups([0.95, 0.8, 0.6, 0.3], cfg)
    assert list(g) == ["Favourable", "Standard", "High", "VeryHigh"]
    assert list(g.categories) == cfg["risk_groups"]["labels"]


def test_bootstrap_ci_columns():
    cfg = toy_cfg(); d = make_toy()
    bs = pk.bootstrap_survival_ci(d, cfg, n=30, seed=1)
    assert set(["id", "pred_original", "pred_bias_corrected",
                "ci_lower", "ci_upper"]).issubset(bs.columns)
    ok = bs["ci_lower"].notna()
    assert (bs.loc[ok, "ci_lower"] <= bs.loc[ok, "ci_upper"] + 1e-9).all()


def test_external_validate_cindex():
    cfg = toy_cfg(); d = make_toy(200)
    fit = pk.fit_cox(d.iloc[:150], cfg)
    ev = pk.external_validate(d.iloc[150:], fit, cfg)
    assert 0.4 < ev["c_index"] <= 1.0


def test_optimise_cutoffs_descending():
    cfg = toy_cfg(); d = make_toy(200)
    fit = pk.fit_cox(d, cfg)
    ps = pk.predict_survival_at_time(fit, d, cfg=cfg)
    oc = pk.optimise_cutoffs(d, ps, cfg, labels=["Low", "Med", "High"],
                             fav_min_surv=0.6, grid=np.round(np.arange(0.3, 0.81, 0.1), 3))
    assert len(oc["best"]) == 2
    assert oc["best"][0] > oc["best"][1]


def test_calibration_eo_ratio():
    cfg = toy_cfg(); d = make_toy()
    fit = pk.fit_cox(d, cfg)
    cal = pk.calibration_table(d, pk.predict_survival_at_time(fit, d, cfg=cfg), cfg, n_bins=4)
    assert np.isfinite(cal["eo_ratio"]) and cal["eo_ratio"] > 0
