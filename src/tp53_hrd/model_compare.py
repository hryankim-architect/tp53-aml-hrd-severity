"""Formal nested Cox PH model comparison.

Three nested models from v0.2:
    M1 = univariate severity_score          (1 parameter)
    M2 = severity_score + hrd_score          (2 parameters)
    M3 = M2 + severity_score:hrd_score       (3 parameters)

This module fits the same models and computes:

    * **LRT** (Likelihood Ratio Test) — chi-square on
      ``-2 * (ll_simple - ll_complex)`` with df = Δ(parameters).
      Tests whether the complex model fits *significantly* better than
      the simple one.
    * **AIC delta** — ``AIC_complex - AIC_simple``. Negative = complex
      model preferred under AIC.
    * **C-index delta** — change in Harrell's C-index (Cox PH
      `concordance_index_`). Positive = complex model has better
      discrimination.

Decision rule (v0.3, deliberately simple):

    "complex model justified" iff (LRT p < 0.05) AND (AIC delta < 0)

The README climax surfaces the three pairwise comparisons (M2-vs-M1,
M3-vs-M2, M3-vs-M1) per run so a reader can see exactly why we did or
did not keep the more complex model.

Caveats (v0.3): with n=15, LRT chi-square approximation is itself
suspect. The decision rule output is reported as a *descriptive* signal
(consistent with / not consistent with the simpler model), not as a
confirmatory hypothesis test.
"""

from __future__ import annotations

import warnings
from dataclasses import asdict, dataclass

import pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import chi2


@dataclass
class ModelComparison:
    """Result of comparing one simple vs one complex nested Cox model."""

    simple: str
    complex: str
    n_params_simple: int
    n_params_complex: int
    ll_simple: float
    ll_complex: float
    lrt_chi2: float
    lrt_df: int
    lrt_p: float
    aic_simple: float
    aic_complex: float
    aic_delta: float  # complex - simple; negative = complex preferred
    cindex_simple: float
    cindex_complex: float
    cindex_delta: float  # complex - simple; positive = complex preferred
    justifies_complex: bool


def _fit_cox(df: pd.DataFrame, covariates: list[str], interaction: bool = False) -> CoxPHFitter:
    """Fit a Cox PH model on the given covariates (+ optional interaction)."""
    fit_df = df[[*covariates, "os_days", "os_event"]].dropna().copy()
    if interaction and len(covariates) == 2:
        a, b = covariates
        fit_df[f"{a}:{b}"] = fit_df[a] * fit_df[b]
    cph = CoxPHFitter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cph.fit(fit_df, duration_col="os_days", event_col="os_event")
    return cph


def compare_nested(
    cox_simple: CoxPHFitter,
    cox_complex: CoxPHFitter,
    name_simple: str,
    name_complex: str,
) -> ModelComparison:
    """Compute LRT + AIC delta + C-index delta for two nested Cox models."""
    ll_s = float(cox_simple.log_likelihood_)
    ll_c = float(cox_complex.log_likelihood_)
    df_s = len(cox_simple.params_)
    df_c = len(cox_complex.params_)
    delta_df = df_c - df_s

    if delta_df <= 0:
        raise ValueError(
            f"complex model must have more parameters than simple: "
            f"got {df_s} -> {df_c}"
        )

    lrt = max(0.0, 2.0 * (ll_c - ll_s))
    p_lrt = float(1.0 - chi2.cdf(lrt, df=delta_df)) if lrt > 0 else 1.0

    aic_s = float(cox_simple.AIC_partial_)
    aic_c = float(cox_complex.AIC_partial_)
    aic_delta = aic_c - aic_s

    c_s = float(cox_simple.concordance_index_)
    c_c = float(cox_complex.concordance_index_)
    c_delta = c_c - c_s

    justifies = (p_lrt < 0.05) and (aic_delta < 0)

    return ModelComparison(
        simple=name_simple,
        complex=name_complex,
        n_params_simple=df_s,
        n_params_complex=df_c,
        ll_simple=ll_s,
        ll_complex=ll_c,
        lrt_chi2=lrt,
        lrt_df=delta_df,
        lrt_p=p_lrt,
        aic_simple=aic_s,
        aic_complex=aic_c,
        aic_delta=aic_delta,
        cindex_simple=c_s,
        cindex_complex=c_c,
        cindex_delta=c_delta,
        justifies_complex=justifies,
    )


def run_nested_model_suite(df: pd.DataFrame) -> dict[str, dict]:
    """Fit M1 (univariate TP53), M2 (bivariate), M3 (+interaction) and
    return all three pairwise comparisons as JSON-serialisable dicts.

    Requires columns: severity_score, hrd_score, os_days, os_event.
    """
    required = {"severity_score", "hrd_score", "os_days", "os_event"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"missing columns: {sorted(missing)}")

    m1 = _fit_cox(df, ["severity_score"])
    m2 = _fit_cox(df, ["severity_score", "hrd_score"])
    m3 = _fit_cox(df, ["severity_score", "hrd_score"], interaction=True)

    return {
        "m2_vs_m1": asdict(compare_nested(m1, m2, "M1_TP53_only", "M2_TP53_plus_HRD")),
        "m3_vs_m2": asdict(compare_nested(m2, m3, "M2_TP53_plus_HRD", "M3_with_interaction")),
        "m3_vs_m1": asdict(compare_nested(m1, m3, "M1_TP53_only", "M3_with_interaction")),
    }
