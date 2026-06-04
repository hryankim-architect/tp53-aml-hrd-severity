"""Causal mediation analysis: TP53 -> HRD -> survival.

Tests the biological hypothesis that TP53 LOF acts as an upstream driver
of HR-repair failure, with HRD-scar burden as a downstream consequence
that then influences survival outcome. The classical Baron-Kenny (1986)
framework is implemented for the survival outcome, with the indirect
effect estimated by:

    indirect = a * b

where:

    a = OLS slope of (HRD_score ~ TP53_severity)
    b = Cox PH log-HR of (HRD_score) when adjusted for TP53_severity

We also estimate a non-parametric bootstrap distribution of the indirect
effect to construct a CI (Sobel's normality assumption is unreliable on
small samples; bootstrap is the modern standard — VanderWeele 2015).

Caveats (v0.3):

    With n=15, a single mediation estimate is descriptive at best. The
    bootstrap CI on the indirect effect will be very wide; the README
    explicitly reports the *direction and consistency* of the estimate
    rather than treating any p-value as confirmatory.

    Causal-mediation interpretation also requires the no-unmeasured-
    confounders assumption (NUCA). In a 15-patient cohort drawn from
    a complex disease (TP53-AML), unmeasured chromosomal events almost
    certainly violate NUCA. The mediation estimate is therefore an
    *association decomposition*, not a strict causal claim.

References:
    Baron & Kenny 1986, J Pers Soc Psychol 51:1173 (classical mediation)
    VanderWeele 2015, Annu Rev Public Health 36:225 (modern causal
        mediation; bootstrap recommendation)
"""

from __future__ import annotations

import warnings
from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from sklearn.linear_model import LinearRegression


@dataclass
class MediationResult:
    """Decomposed effect of TP53 on survival via HRD mediator."""

    n: int
    n_events: int

    # Path a: TP53 -> HRD (OLS)
    a_coef: float
    a_pvalue: float | None  # OLS p; lifelines doesn't expose this, computed manually below

    # Path b: HRD -> survival adjusted for TP53 (Cox)
    b_log_hr: float        # log-HR (Cox coef)
    b_hr: float            # exp(b_log_hr)
    b_pvalue: float

    # Direct effect: TP53 -> survival adjusted for HRD (Cox)
    direct_log_hr: float
    direct_hr: float
    direct_pvalue: float

    # Total effect: TP53 -> survival (univariate Cox)
    total_log_hr: float
    total_hr: float
    total_pvalue: float

    # Indirect effect (Baron-Kenny product)
    indirect_log_hr: float                 # a * b_log_hr (in log-hazard scale)
    indirect_ci_low: float | None          # bootstrap 2.5%-ile
    indirect_ci_high: float | None         # bootstrap 97.5%-ile
    proportion_mediated: float | None      # indirect / total (NaN if total ~ 0)

    n_bootstrap: int


def _path_a_ols(
    df: pd.DataFrame, treatment: str, mediator: str
) -> tuple[float, float]:
    """OLS regression of mediator on treatment. Returns (coef, p-value)."""
    x = df[[treatment]].values
    y = df[mediator].values
    reg = LinearRegression().fit(x, y)
    a = float(reg.coef_[0])

    # Compute p-value manually (sklearn doesn't expose it)
    n = len(x)
    if n <= 2:
        return a, None
    y_pred = reg.predict(x)
    resid = y - y_pred
    rss = float(np.sum(resid ** 2))
    dof = n - 2
    sigma2 = rss / dof if dof > 0 else float("nan")
    x_centered = x[:, 0] - x[:, 0].mean()
    sxx = float(np.sum(x_centered ** 2))
    if sxx <= 0 or not np.isfinite(sigma2):
        return a, None
    se_a = float(np.sqrt(sigma2 / sxx))
    if se_a <= 0:
        return a, None
    t_stat = a / se_a
    from scipy.stats import t as t_dist
    p_val = float(2.0 * (1.0 - t_dist.cdf(abs(t_stat), df=dof)))
    return a, p_val


def _fit_cox(
    df: pd.DataFrame, covariates: list[str]
) -> CoxPHFitter:
    cph = CoxPHFitter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cph.fit(
            df[[*covariates, "os_days", "os_event"]].dropna(),
            duration_col="os_days",
            event_col="os_event",
        )
    return cph


def run_mediation(
    cohort: pd.DataFrame,
    treatment: str = "severity_score",
    mediator: str = "hrd_score",
    n_bootstrap: int = 1000,
    seed: int = 42,
) -> MediationResult:
    """Baron-Kenny mediation analysis with bootstrap CI on indirect effect.

    Args:
        cohort: DataFrame with treatment, mediator, os_days, os_event columns.
        treatment: name of upstream variable (default severity_score = TP53 axis).
        mediator: name of putative mediator (default hrd_score = HRD-scar axis).
        n_bootstrap: number of bootstrap resamples (set to 0 to skip).
        seed: RNG seed for reproducibility.

    Returns:
        MediationResult with all paths + indirect-effect bootstrap CI.
    """
    required = {treatment, mediator, "os_days", "os_event"}
    missing = required - set(cohort.columns)
    if missing:
        raise KeyError(f"mediation needs columns: {sorted(missing)}")

    df = cohort[[*required]].dropna().copy()
    n = len(df)
    n_events = int(df["os_event"].sum())

    # --- Path a: TP53 -> HRD (OLS)
    a, a_p = _path_a_ols(df, treatment, mediator)

    # --- Path b: HRD -> survival adjusted for TP53 (Cox)
    cox_full = _fit_cox(df, [treatment, mediator])
    b_log = float(cox_full.params_.get(mediator, float("nan")))
    b_p = float(cox_full.summary.loc[mediator, "p"])
    b_hr = float(np.exp(b_log)) if np.isfinite(b_log) else float("nan")
    direct_log = float(cox_full.params_.get(treatment, float("nan")))
    direct_p = float(cox_full.summary.loc[treatment, "p"])
    direct_hr = float(np.exp(direct_log)) if np.isfinite(direct_log) else float("nan")

    # --- Total effect: TP53 -> survival (univariate Cox)
    cox_total = _fit_cox(df, [treatment])
    total_log = float(cox_total.params_.get(treatment, float("nan")))
    total_p = float(cox_total.summary.loc[treatment, "p"])
    total_hr = float(np.exp(total_log)) if np.isfinite(total_log) else float("nan")

    # --- Indirect (a * b)
    indirect = a * b_log

    # --- Bootstrap CI on indirect effect
    ci_low: float | None = None
    ci_high: float | None = None
    if n_bootstrap > 0 and n >= 5:
        rng = np.random.default_rng(seed)
        boot_indirect = []
        for _ in range(n_bootstrap):
            idx = rng.integers(0, n, size=n)
            boot_df = df.iloc[idx].reset_index(drop=True)
            # Ensure both event and non-event survive resample
            if boot_df["os_event"].sum() == 0 or boot_df["os_event"].sum() == len(boot_df):
                continue
            try:
                a_b, _ = _path_a_ols(boot_df, treatment, mediator)
                cph_b = _fit_cox(boot_df, [treatment, mediator])
                b_b = float(cph_b.params_.get(mediator, float("nan")))
                if np.isfinite(a_b) and np.isfinite(b_b):
                    boot_indirect.append(a_b * b_b)
            except Exception:  # noqa: BLE001 — bootstrap iterations may fail in small-n
                continue
        if boot_indirect:
            ci_low = float(np.percentile(boot_indirect, 2.5))
            ci_high = float(np.percentile(boot_indirect, 97.5))

    # --- Proportion mediated (indirect / total)
    proportion: float | None = (
        None if abs(total_log) < 1e-9 else float(indirect / total_log)
    )

    return MediationResult(
        n=n,
        n_events=n_events,
        a_coef=float(a),
        a_pvalue=a_p,
        b_log_hr=b_log,
        b_hr=b_hr,
        b_pvalue=b_p,
        direct_log_hr=direct_log,
        direct_hr=direct_hr,
        direct_pvalue=direct_p,
        total_log_hr=total_log,
        total_hr=total_hr,
        total_pvalue=total_p,
        indirect_log_hr=float(indirect),
        indirect_ci_low=ci_low,
        indirect_ci_high=ci_high,
        proportion_mediated=proportion,
        n_bootstrap=n_bootstrap,
    )


def mediation_to_dict(result: MediationResult) -> dict:
    """JSON-serialisable view of the mediation result."""
    return asdict(result)
