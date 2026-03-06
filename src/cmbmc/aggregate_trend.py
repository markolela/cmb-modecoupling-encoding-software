"""Aggregation, trend fit, detrending, and effect definitions.

Implements binding Sections 2.8, 2.9, 2.10.

Binding scale order.
s_list = [1, 2, 4]

Binding aggregation.
kappa_tilde = median over patches of K[p, t]

Binding trend model.
x = log(s_list)
X = [1, x]
beta = lstsq(X, y)
a = beta[0]
theta_data = beta[1]

Binding detrending.
kappa_theta[t] = kappa_tilde[t] - theta_data * log(s_list[t])

Binding effects.
Delta theta = theta_null_ref - theta_data
theta_null_ref = median_r theta_null[r]
CI95 via numpy.quantile(method="linear") applied to Dtheta = {theta_null[r] - theta_data}

PP = 100 * (max(kappa_theta) - min(kappa_theta)) / max(|mean(kappa_theta)|, eps)
Delta PP = PP_data - PP_null_ref
PP_null_ref = median_r PP_null[r]
CI95(Delta PP) via DPP = {PP_data - PP_null[r]} and quantiles method="linear"
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


S_LIST = [1, 2, 4]
EPS = np.float64(1e-12)


def quantile_linear(x: np.ndarray, q: float) -> np.float64:
    a = np.asarray(x, dtype=np.float64)
    return np.float64(np.quantile(a, q=float(q), method="linear"))


def median_aggregate_kappa(K: np.ndarray) -> np.ndarray:
    """K is patch by scale matrix. Returns kappa_tilde in binding scale order."""
    x = np.asarray(K, dtype=np.float64)
    if x.ndim != 2:
        raise ValueError(f"K must be 2D, got shape {x.shape}")
    if x.shape[1] != 3:
        raise ValueError(f"K must have 3 columns for s_list [1,2,4], got {x.shape}")
    return np.asarray(np.median(x, axis=0), dtype=np.float64, order="C")


def fit_theta(kappa_tilde: np.ndarray) -> tuple[np.float64, np.float64]:
    """Binding OLS fit. Returns (a, theta_data)."""
    y = np.asarray(kappa_tilde, dtype=np.float64)
    if y.shape != (3,):
        raise ValueError(f"kappa_tilde must be shape (3,), got {y.shape}")

    s_arr = np.asarray(S_LIST, dtype=np.float64)
    x = np.log(s_arr)

    X = np.zeros((3, 2), dtype=np.float64)
    X[:, 0] = np.float64(1.0)
    X[:, 1] = x

    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    a = np.float64(beta[0])
    theta = np.float64(beta[1])
    return a, theta


def detrend(kappa_tilde: np.ndarray, theta_data: np.float64) -> np.ndarray:
    """Binding detrending. Returns kappa_theta in s_list order."""
    y = np.asarray(kappa_tilde, dtype=np.float64)
    s_arr = np.asarray(S_LIST, dtype=np.float64)
    out = y - np.float64(theta_data) * np.log(s_arr)
    return np.asarray(out, dtype=np.float64, order="C")


def pp_from_kappa_theta(kappa_theta: np.ndarray) -> np.float64:
    """Binding PP definition."""
    k = np.asarray(kappa_theta, dtype=np.float64)
    mean_theta = np.mean(k, dtype=np.float64)
    denom = np.maximum(np.abs(mean_theta), EPS)
    val = np.float64(100.0) * (np.max(k) - np.min(k)) / denom
    return np.float64(val)


@dataclass(frozen=True)
class Effects:
    theta_data: np.float64
    pp_data: np.float64
    theta_null_ref: np.float64
    pp_null_ref: np.float64
    delta_theta: np.float64
    delta_pp: np.float64
    ci95_delta_theta: tuple[np.float64, np.float64]
    ci95_delta_pp: tuple[np.float64, np.float64]


def effects_with_ci95(
    *,
    theta_data: np.float64,
    pp_data: np.float64,
    theta_null: np.ndarray,
    pp_null: np.ndarray,
) -> Effects:
    """Compute binding effects and CI95 against a given null ensemble."""
    th_null = np.asarray(theta_null, dtype=np.float64)
    pp_n = np.asarray(pp_null, dtype=np.float64)

    if th_null.ndim != 1 or pp_n.ndim != 1:
        raise ValueError("theta_null and pp_null must be 1D")
    if th_null.shape[0] != pp_n.shape[0]:
        raise ValueError("theta_null and pp_null must have same length")

    theta_null_ref = np.float64(np.median(th_null))
    pp_null_ref = np.float64(np.median(pp_n))

    delta_theta = np.float64(theta_null_ref) - np.float64(theta_data)
    dtheta = th_null - np.float64(theta_data)
    ci_theta = (quantile_linear(dtheta, 0.025), quantile_linear(dtheta, 0.975))

    delta_pp = np.float64(pp_data) - np.float64(pp_null_ref)
    dpp = np.float64(pp_data) - pp_n
    ci_pp = (quantile_linear(dpp, 0.025), quantile_linear(dpp, 0.975))

    return Effects(
        theta_data=np.float64(theta_data),
        pp_data=np.float64(pp_data),
        theta_null_ref=theta_null_ref,
        pp_null_ref=pp_null_ref,
        delta_theta=delta_theta,
        delta_pp=delta_pp,
        ci95_delta_theta=ci_theta,
        ci95_delta_pp=ci_pp,
    )
