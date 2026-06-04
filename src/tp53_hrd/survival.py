"""Survival analysis on the TP53-HRD severity-scored cohort.

The capability claim P3 makes is "*method connects to overall survival*".
This module is where that connection lives:

* **Kaplan-Meier curves** per severity band (low / moderate / high). Returns a
  tidy summary DataFrame (n, n_events, median OS, 1-year and 3-year survival)
  rather than the raw KaplanMeierFitter objects, so the result is JSON-safe
  for the audit ledger and the report.
* **Multivariate log-rank test** across the 3 bands, plus a 2-arm fallback
  (high vs not-high) for when n in any single band is too small for a stable
  3-arm test (a realistic risk at n=15 — see ``docs/what-is-out-of-scope.md``).
* **Cox proportional hazards** with ``severity_score`` as a continuous
  covariate. Returns HR, 95% CI, log-likelihood p-value, and concordance.
* **KM plot** rendered to a static PNG using matplotlib's Agg backend so the
  function works in headless CI without a display server.

n=15 is **below standard statistical power**. The README disclaims this; the
goal is method demonstration. The Cox HR direction (severity ↑ ⇒ hazard ↑)
is informative even when CIs are wide.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")  # headless, must be set before importing pyplot
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from lifelines import CoxPHFitter, KaplanMeierFitter  # noqa: E402
from lifelines.statistics import multivariate_logrank_test  # noqa: E402

# Severity-band display order (low risk first for the KM plot legend)
BAND_ORDER = ("low", "moderate", "high")

# Days lookup for "1-year" and "3-year" survival probabilities
ONE_YEAR_DAYS = 365.0
THREE_YEAR_DAYS = 1095.0


def kaplan_meier_summary(
    cohort: pd.DataFrame, group_col: str = "severity_band"
) -> pd.DataFrame:
    """Per-group KM summary statistics.

    Columns: ``group``, ``n``, ``n_events``, ``median_os_days``,
    ``surv_1yr``, ``surv_3yr``.
    """
    _require_cols(cohort, [group_col, "os_days", "os_event"])

    rows: list[dict[str, Any]] = []
    for group in _ordered_groups(cohort[group_col]):
        sub = cohort[cohort[group_col] == group]
        if sub.empty:
            continue

        kmf = KaplanMeierFitter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            kmf.fit(
                durations=sub["os_days"],
                event_observed=sub["os_event"],
                label=group,
            )

        median = kmf.median_survival_time_
        median_clean = float(median) if pd.notna(median) and median != float("inf") else None

        rows.append(
            {
                "group": group,
                "n": int(len(sub)),
                "n_events": int(sub["os_event"].sum()),
                "median_os_days": median_clean,
                "surv_1yr": _survival_at(kmf, ONE_YEAR_DAYS),
                "surv_3yr": _survival_at(kmf, THREE_YEAR_DAYS),
            }
        )

    return pd.DataFrame(rows)


def multivariate_logrank_p(
    cohort: pd.DataFrame, group_col: str = "severity_band"
) -> float:
    """p-value from the multivariate log-rank test across groups."""
    _require_cols(cohort, [group_col, "os_days", "os_event"])
    result = multivariate_logrank_test(
        cohort["os_days"], cohort[group_col], cohort["os_event"]
    )
    return float(result.p_value)


def two_arm_summary(
    cohort: pd.DataFrame,
    high_band: str = "high",
    group_col: str = "severity_band",
) -> dict[str, Any]:
    """Fallback comparison: ``high`` band vs every other band combined.

    Useful when at least one of the 3 bands is too small for a stable
    3-arm test. Returns KM summary for both arms plus the log-rank p.
    """
    _require_cols(cohort, [group_col, "os_days", "os_event"])
    arm = cohort[group_col].apply(lambda b: "high" if b == high_band else "not-high")
    cohort_2arm = cohort.assign(_arm=arm)

    summary = kaplan_meier_summary(cohort_2arm, group_col="_arm")
    p = multivariate_logrank_p(cohort_2arm, group_col="_arm")
    return {"summary": summary, "logrank_p": p}


def cox_severity(
    cohort: pd.DataFrame, covariate: str = "severity_score"
) -> dict[str, Any]:
    """Cox PH with a single continuous covariate (default: severity_score).

    Returns ``{"hr", "ci_low", "ci_high", "p_value", "concordance", "n", "n_events"}``.
    """
    _require_cols(cohort, [covariate, "os_days", "os_event"])

    df = cohort[[covariate, "os_days", "os_event"]].dropna().copy()
    cph = CoxPHFitter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cph.fit(df, duration_col="os_days", event_col="os_event")

    s = cph.summary.loc[covariate]
    return {
        "covariate": covariate,
        "hr": float(s["exp(coef)"]),
        "ci_low": float(s["exp(coef) lower 95%"]),
        "ci_high": float(s["exp(coef) upper 95%"]),
        "p_value": float(s["p"]),
        "concordance": float(cph.concordance_index_),
        "n": int(len(df)),
        "n_events": int(df["os_event"].sum()),
    }


def cox_bivariate(
    cohort: pd.DataFrame,
    covariates: list[str],
    interaction: bool = False,
) -> dict[str, Any]:
    """Multi-covariate Cox PH on `covariates`, optionally with an interaction term.

    Args:
        cohort: must have all covariates + os_days + os_event.
        covariates: list of column names to include as covariates. For the
            v0.2 TP53+HRD analysis use ``["severity_score", "hrd_score"]``.
        interaction: if True and len(covariates)==2, adds a multiplicative
            interaction term (``covA * covB``) named ``"<a>:<b>"``.

    Returns ``{"model_summary": {covariate: {hr, ci_low, ci_high, p_value}},
                "concordance": float, "n": int, "n_events": int,
                "interaction_term": str | None}``.

    Limitations: this is a multi-parameter Cox; with n<<10/parameter,
    results are descriptive at best. The pipeline emits this regardless
    so the comparison vs the univariate baseline is visible, but the
    README explicitly downgrades the language.
    """
    _require_cols(cohort, [*covariates, "os_days", "os_event"])
    df = cohort[[*covariates, "os_days", "os_event"]].dropna().copy()

    interaction_term: str | None = None
    if interaction and len(covariates) == 2:
        a, b = covariates
        interaction_term = f"{a}:{b}"
        df[interaction_term] = df[a] * df[b]

    cph = CoxPHFitter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cph.fit(df, duration_col="os_days", event_col="os_event")

    out_covs = list(covariates) + ([interaction_term] if interaction_term else [])
    model_summary = {}
    for cov in out_covs:
        if cov not in cph.summary.index:
            continue
        s = cph.summary.loc[cov]
        model_summary[cov] = {
            "hr": float(s["exp(coef)"]),
            "ci_low": float(s["exp(coef) lower 95%"]),
            "ci_high": float(s["exp(coef) upper 95%"]),
            "p_value": float(s["p"]),
        }

    return {
        "model_summary": model_summary,
        "concordance": float(cph.concordance_index_),
        "n": int(len(df)),
        "n_events": int(df["os_event"].sum()),
        "interaction_term": interaction_term,
    }


def make_km_plot(
    cohort: pd.DataFrame,
    out_path: Path,
    group_col: str = "severity_band",
    title: str = "TP53-HRD severity vs overall survival (TCGA-LAML)",
) -> Path:
    """Render a KM curve plot (one curve per group) to ``out_path``.

    Returns the path written.
    """
    _require_cols(cohort, [group_col, "os_days", "os_event"])
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 5))
    groups = _ordered_groups(cohort[group_col])
    for group in groups:
        sub = cohort[cohort[group_col] == group]
        if sub.empty:
            continue
        kmf = KaplanMeierFitter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            kmf.fit(
                durations=sub["os_days"],
                event_observed=sub["os_event"],
                label=f"{group} (n={len(sub)})",
            )
            kmf.plot_survival_function(ax=ax, ci_show=True)

    ax.set_title(title)
    ax.set_xlabel("Days from diagnosis")
    ax.set_ylabel("Overall survival probability")
    ax.set_ylim(0.0, 1.05)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    return out_path


# --- internals --------------------------------------------------------------

def _require_cols(df: pd.DataFrame, cols: list[str]) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"missing required columns: {missing}")


def _ordered_groups(series: pd.Series) -> list[str]:
    """Return unique group labels with ``BAND_ORDER`` priority, then alphabetical."""
    present = list(pd.unique(series.dropna()))
    ordered = [g for g in BAND_ORDER if g in present]
    extras = sorted(g for g in present if g not in BAND_ORDER)
    return ordered + extras


def _survival_at(kmf: KaplanMeierFitter, day: float) -> float | None:
    """Interpolated survival probability at ``day``. Returns None if out of range."""
    sf = kmf.survival_function_
    if sf.empty:
        return None
    max_day = sf.index.max()
    if day > max_day:
        return None
    return float(kmf.predict(day))
