#!/usr/bin/env python3
"""
plot_beam_h27_comparison.py

Two groups of figures:
  A) Legacy-only: compare beam_h27 meshes (10x2x2 → 80x2x2) within MORFE2.0.
  B) Legacy vs New: matched pairs, legacy MORFE2.0 vs Julia/Ferrite benchmark.

Cumulative-time reconstruction for legacy follows plot_time_cost_data.py:
  total_cumul = monomial_total_time_s.cumsum() + order_step
where order_step is a step function assigning the cumulative order-level cost
(fillrhsG + fillrhsH + fillWf) to every monomial within that order.

Saves all figures to plots_comparison/ next to this script.
"""

import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import seaborn as sns
from pathlib import Path
from numpy.polynomial import Polynomial

# ── Directories ────────────────────────────────────────────────────────────────
OUTPUT_DIR = Path(__file__).parent.resolve()
NEW_ROOT   = (OUTPUT_DIR / "../../../demo/BenchmarkFerrite/benchmark_results").resolve()
PLOTS_DIR  = OUTPUT_DIR / "plots_comparison"

MAX_ORDER = 11

# ── Regex ──────────────────────────────────────────────────────────────────────
LEG_RE = re.compile(r"^beam_h27_(.+?)_o_(\d+)_eps_\d+_Fmode_\d+_.+$")
NEW_RE = re.compile(r"^beam_h27_(.+?)_degree(\d+)_(\w+)$")
SUM_RE = re.compile(r"^(\w+)\s*=\s*(.+)$")

# ── Style constants ────────────────────────────────────────────────────────────
C_LEGACY    = "steelblue"
C_NEW       = "darkorange"
FILL_LEGACY = "lightsteelblue"
FILL_NEW    = "bisque"


# ─────────────────────────────────────────────────────────────────────────────
# Discovery
# ─────────────────────────────────────────────────────────────────────────────
def _discover_legacy():
    runs: dict = {}
    for d in sorted(OUTPUT_DIR.glob("beam_h27_*")):
        if not d.is_dir():
            continue
        m = LEG_RE.match(d.name)
        if not m:
            continue
        mesh, order = m.group(1), int(m.group(2))
        if order < MAX_ORDER:
            continue
        if mesh not in runs or d.name > runs[mesh][1].name:
            runs[mesh] = (order, d)
    return runs


def _discover_new():
    runs: dict = {}
    for d in sorted(NEW_ROOT.glob("beam_h27_*")):
        if not d.is_dir():
            continue
        m = NEW_RE.match(d.name)
        if not m:
            continue
        mesh, order = m.group(1), int(m.group(2))
        if order != MAX_ORDER:
            continue
        if mesh not in runs or d.name > runs[mesh][1].name:
            runs[mesh] = (order, d)
    return runs


# ─────────────────────────────────────────────────────────────────────────────
# Core loaders
# ─────────────────────────────────────────────────────────────────────────────
def _load_summary_new(path: Path) -> dict:
    data: dict = {}
    for line in (path / "summary.txt").read_text().splitlines():
        m = SUM_RE.match(line.strip())
        if m:
            try:
                data[m.group(1)] = float(m.group(2))
            except ValueError:
                data[m.group(1)] = m.group(2).strip()
    return data


def _load_new_order(path: Path) -> pd.DataFrame:
    return pd.read_csv(path / "benchmark_per_order.csv")


def _load_new_mono(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path / "benchmark_per_monomial.csv")
    df["exp_tuple"] = df["exponents"].apply(
        lambda s: tuple(int(x) for x in s.split("_"))
    )
    return df


def _load_leg_order_raw(path: Path) -> pd.DataFrame:
    """Legacy per-order CSV truncated to MAX_ORDER, with recomputed cumulative."""
    df = pd.read_csv(path / "benchmark_per_order.csv")
    df = df[df["order"] <= MAX_ORDER].copy().reset_index(drop=True)
    df["rhs_time_s"]         = df["fillrhsG_time_s"] + df["fillrhsH_time_s"]
    df["solve_time_s"]       = df["fillWf_time_s"]
    df["order_total_time_s"] = df["rhs_time_s"] + df["solve_time_s"]
    df["cumul_time_s"]       = df["order_total_time_s"].cumsum()
    return df


def _leg_total_cumul(path: Path):
    """
    Reconstruct legacy total-cumulative time following plot_time_cost_data.py:
      total = cumsum(monomial_total_time_s) + order_step_function

    Returns
    -------
    global_idx   : 1-D int array  (1-based global monomial index, +4 offset)
    total_s      : 1-D float array (total cumulative time in seconds)
    x_ends       : 1-D float array (global idx of last monomial of each order)
    orders       : 1-D int array  (order labels)
    df_ord       : DataFrame with per-order data
    """
    df_mono = pd.read_csv(path / "benchmark_per_monomial.csv", na_values="NaN")
    df_mono["monomial_total_time_s"] = df_mono["monomial_total_time_s"].fillna(0)
    df_mono = df_mono[df_mono["order"] <= MAX_ORDER].copy().reset_index(drop=True)

    df_ord = pd.read_csv(path / "benchmark_per_order.csv")
    df_ord = df_ord[df_ord["order"] <= MAX_ORDER].copy().reset_index(drop=True)
    df_ord["total_time_order"] = (
        df_ord["fillrhsG_time_s"] + df_ord["fillrhsH_time_s"] + df_ord["fillWf_time_s"]
    )
    cumul_ord = df_ord["total_time_order"].cumsum().values
    df_ord["cumul_monomials"] = df_ord["n_monomials"].cumsum().values
    starts = [1] + (df_ord["cumul_monomials"].iloc[:-1] + 1).tolist()

    cumul_mono = df_mono["monomial_total_time_s"].cumsum().values
    step = np.zeros(len(df_mono))
    for i, s in enumerate(starts):
        e = int(df_ord["cumul_monomials"].iloc[i])
        step[s - 1 : e] = cumul_ord[i]

    total = cumul_mono + step
    # +4 offset: monomials 1-4 are linear modes not in the CSV
    global_idx = np.arange(1, len(df_mono) + 1) + 4
    x_ends     = df_ord["cumul_monomials"].values + 4
    orders     = df_ord["order"].values
    return global_idx, total, x_ends, orders, df_ord


def _load_leg_mono_with_cumul(path: Path) -> pd.DataFrame:
    """
    Legacy per-monomial CSV (filtered to MAX_ORDER) with total_cumul_s column
    computed via the plot_time_cost_data.py reconstruction, and exp_tuple parsed.
    """
    df_mono = pd.read_csv(path / "benchmark_per_monomial.csv", na_values="NaN")
    df_mono["monomial_total_time_s"] = df_mono["monomial_total_time_s"].fillna(0)
    df_mono = df_mono[df_mono["order"] <= MAX_ORDER].copy().reset_index(drop=True)

    df_ord = pd.read_csv(path / "benchmark_per_order.csv")
    df_ord = df_ord[df_ord["order"] <= MAX_ORDER].copy().reset_index(drop=True)
    df_ord["total_time_order"] = (
        df_ord["fillrhsG_time_s"] + df_ord["fillrhsH_time_s"] + df_ord["fillWf_time_s"]
    )
    cumul_ord = df_ord["total_time_order"].cumsum().values
    df_ord["cumul_monomials"] = df_ord["n_monomials"].cumsum().values
    starts = [1] + (df_ord["cumul_monomials"].iloc[:-1] + 1).tolist()

    cumul_mono = df_mono["monomial_total_time_s"].cumsum().values
    step = np.zeros(len(df_mono))
    for i, s in enumerate(starts):
        e = int(df_ord["cumul_monomials"].iloc[i])
        step[s - 1 : e] = cumul_ord[i]

    df_mono["total_cumul_s"] = cumul_mono + step
    df_mono["exp_tuple"] = df_mono["alpha_vector"].apply(
        lambda s: tuple(int(x) for x in s.strip("[]").replace(" ", "").split(","))
    )
    return df_mono


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
def _save(fig: plt.Figure, name: str) -> None:
    PLOTS_DIR.mkdir(exist_ok=True)
    fig.savefig(PLOTS_DIR / name, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}")


def _fom_palette(fom_list) -> dict:
    unique = sorted(set(fom_list))
    colors = sns.color_palette("viridis", len(unique))
    return {fom: colors[i] for i, fom in enumerate(unique)}


def _r2(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - y_true.mean()) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0


def _order_vlines(ax, x_ends, orders, x_min):
    for i, (x_end, order) in enumerate(zip(x_ends, orders)):
        if i < len(x_ends) - 1:
            ax.axvline(x=x_end + 0.5, color="grey", ls="--", lw=0.8, alpha=0.4)
        x_start = x_ends[i - 1] if i > 0 else x_min
        ax.text(
            (x_start + x_end) / 2, 1.01, f"ord {order}",
            fontsize=7, color="grey", ha="center", va="bottom",
            transform=ax.get_xaxis_transform(),
        )


# ─────────────────────────────────────────────────────────────────────────────
# Build run tables
# ─────────────────────────────────────────────────────────────────────────────
leg_runs = _discover_legacy()
new_runs = _discover_new()
common   = sorted(set(leg_runs) & set(new_runs))
print(f"Legacy meshes : {sorted(leg_runs)}")
print(f"New meshes    : {sorted(new_runs)}")
print(f"Matched pairs : {common}")

# Legacy-only list (sorted by FOM via summary lookup where available)
leg_list = []
for mesh, (order, path) in sorted(leg_runs.items()):
    # FOM from matching new run, else None
    fom = None
    if mesh in new_runs:
        try:
            fom = int(_load_summary_new(new_runs[mesh][1]).get("FOM", 0))
        except Exception:
            pass
    leg_list.append({"mesh": mesh, "fom": fom, "path": path})
leg_list.sort(key=lambda r: (r["fom"] or 0, r["mesh"]))

# Matched pairs list
pairs = []
for mesh in common:
    _, leg_path = leg_runs[mesh]
    _, new_path = new_runs[mesh]
    fom = int(_load_summary_new(new_path).get("FOM", 0))
    pairs.append({"mesh": mesh, "fom": fom, "leg": leg_path, "new": new_path})
pairs.sort(key=lambda p: p["fom"])

foms_leg   = [r["fom"] for r in leg_list if r["fom"]]
foms_pairs = [p["fom"] for p in pairs]
leg_palette  = _fom_palette(foms_leg)
pair_palette = _fom_palette(foms_pairs)


# ═════════════════════════════════════════════════════════════════════════════
# A) LEGACY-ONLY PLOTS
# ═════════════════════════════════════════════════════════════════════════════

# ── A1: Cumulative time — all legacy meshes on one figure ─────────────────────
def plot_A1_leg_cumulative():
    fig, ax = plt.subplots(figsize=(12, 6))

    for r in leg_list:
        gidx, total_s, x_ends, orders, _ = _leg_total_cumul(r["path"])
        y_min = total_s / 60

        c   = leg_palette[r["fom"]]
        lbl = f"FOM={r['fom']}  ({r['mesh']})"

        # Polynomial fit at order endpoints
        ye  = total_s[x_ends - 4 - 1] / 60   # values at last monomial of each order
        pol = Polynomial.fit(x_ends.astype(float), ye, deg=3)
        r2  = _r2(ye, pol(x_ends.astype(float)))
        x_sm = np.linspace(gidx[0], gidx[-1], 400)

        ax.fill_between(gidx, y_min, alpha=0.18, color=c)
        ax.plot(gidx, y_min, color=c, lw=2, label=f"{lbl}  R²={r2:.4f}")
        ax.plot(x_sm, pol(x_sm), color=c, lw=1.5, ls="--", alpha=0.85)
        ax.scatter(x_ends, ye, color=c, s=28, zorder=5)

    # Order boundaries from the first (smallest) run
    gidx0, _, x_ends0, orders0, _ = _leg_total_cumul(leg_list[0]["path"])
    _order_vlines(ax, x_ends0, orders0, gidx0[0])

    ax.set_xlabel("Global monomial index")
    ax.set_ylabel("Cumulative time (min)")
    ax.set_title("Legacy MORFE2.0 — cumulative solve time, beam_h27 meshes")
    ax.set_xlim(gidx0[0] - 2, gidx0[-1] + 2)
    ax.set_ylim(0)
    ax.legend(fontsize=9)
    ax.grid(False)
    plt.tight_layout()
    _save(fig, "plot_A1_leg_cumulative.png")


# ── A2: Per-order time — all legacy meshes (log scale) ───────────────────────
def plot_A2_leg_per_order():
    fig, ax = plt.subplots(figsize=(10, 5))

    for r in leg_list:
        _, _, _, orders, df_ord = _leg_total_cumul(r["path"])
        c   = leg_palette[r["fom"]]
        lbl = f"FOM={r['fom']}  ({r['mesh']})"
        ax.plot(orders, df_ord["total_time_order"].values / 60,
                color=c, lw=2, marker="o", ms=5, label=lbl)
        ax.fill_between(orders, df_ord["total_time_order"].values / 60,
                        alpha=0.15, color=c)

    ax.set_yscale("log")
    ax.set_xlabel("Polynomial order")
    ax.set_ylabel("Time per order (min, log scale)")
    ax.set_title("Legacy MORFE2.0 — per-order solve time, beam_h27 meshes")
    ax.set_xticks(range(2, MAX_ORDER + 1))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:.3g}"))
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    _save(fig, "plot_A2_leg_per_order.png")


# ── A3: Scaling — legacy only (log-log, FOM vs total time) ───────────────────
def plot_A3_leg_scaling():
    fig, ax = plt.subplots(figsize=(8, 5))

    fom_arr, t_arr = [], []
    for r in leg_list:
        if not r["fom"]:
            continue
        _, _, _, _, df_ord = _leg_total_cumul(r["path"])
        fom_arr.append(r["fom"])
        t_arr.append(df_ord["cumul_time_s"].iloc[-1] / 60)

    fom_arr = np.array(fom_arr)
    t_arr   = np.array(t_arr)

    lf, lt = np.log(fom_arr), np.log(t_arr)
    alpha_fit, c_log = np.polyfit(lf, lt, 1)[0], np.polyfit(lf, lt, 1)[1]
    C = np.exp(c_log)
    r2 = _r2(lt, alpha_fit * lf + c_log)
    x_sm = np.linspace(fom_arr.min() * 0.9, fom_arr.max() * 1.1, 300)

    ax.scatter(fom_arr, t_arr, color=C_LEGACY, s=60, zorder=5)
    ax.plot(fom_arr, t_arr, color=C_LEGACY, lw=1.5, alpha=0.5)
    ax.plot(x_sm, C * x_sm ** alpha_fit, color=C_LEGACY, ls="--", lw=2,
            label=f"Legacy  {C:.2e}·FOM^{alpha_fit:.2f}  R²={r2:.4f}")

    for fom, t in zip(fom_arr, t_arr):
        ax.annotate(str(fom), (fom, t), textcoords="offset points",
                    xytext=(6, 4), fontsize=8, color=C_LEGACY)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("FOM size (DOFs)")
    ax.set_ylabel("Total solve time at order 11 (min)")
    ax.set_title("Legacy MORFE2.0 — FOM scaling, beam_h27 meshes")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.2)
    plt.tight_layout()
    _save(fig, "plot_A3_leg_scaling.png")


# ── A4: Memory (mem_live_bytes) per order — all legacy meshes ─────────────────
def plot_A4_leg_memory():
    fig, ax = plt.subplots(figsize=(10, 5))

    for r in leg_list:
        _, _, _, orders, df_ord = _leg_total_cumul(r["path"])
        c   = leg_palette[r["fom"]]
        lbl = f"FOM={r['fom']}  ({r['mesh']})"
        mem_gb = df_ord["mem_live_bytes"].values / (1024**3)
        ax.plot(orders, mem_gb, color=c, lw=2, marker="o", ms=5, label=lbl)
        ax.fill_between(orders, mem_gb, alpha=0.15, color=c)

    ax.set_xlabel("Polynomial order")
    ax.set_ylabel("Live memory (GB)")
    ax.set_title("Legacy MORFE2.0 — live memory per order, beam_h27 meshes")
    ax.set_xticks(range(2, MAX_ORDER + 1))
    ax.set_ylim(0)
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    _save(fig, "plot_A4_leg_memory.png")


# ═════════════════════════════════════════════════════════════════════════════
# B) LEGACY vs NEW COMPARISON PLOTS
# ═════════════════════════════════════════════════════════════════════════════

# ── B1: Speedup ratio vs order ────────────────────────────────────────────────
def plot_B1_speedup_ratio():
    fig, ax = plt.subplots(figsize=(9, 5))

    for p in pairs:
        df_leg = _load_leg_order_raw(p["leg"])
        df_new = _load_new_order(p["new"])
        orders = df_leg["order"].values
        ratio  = df_leg["cumul_time_s"].values / df_new["cumul_time_s"].values
        c = pair_palette[p["fom"]]
        ax.plot(orders, ratio, color=c, lw=2, marker="o", ms=6,
                label=f"FOM={p['fom']}  ({p['mesh']})")
        ax.scatter(orders, ratio, color=c, s=40, zorder=5)

    ax.axhline(1.0, color="black", ls="--", lw=1.0, alpha=0.6)
    ax.set_xlabel("Polynomial order")
    ax.set_ylabel("Cumulative-time ratio  legacy / new")
    ax.set_title("Speedup: Legacy MORFE2.0 vs New BenchmarkFerrite")
    ax.set_xticks(range(2, MAX_ORDER + 1))
    ax.set_ylim(0)
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    _save(fig, "plot_B1_speedup_ratio.png")


# ── B2: Per-order time — grouped bars, log scale ─────────────────────────────
def plot_B2_per_order():
    n = len(pairs)
    ncols = min(2, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 4.5 * nrows),
                              squeeze=False)

    for idx, p in enumerate(pairs):
        ax = axes[idx // ncols][idx % ncols]
        df_leg = _load_leg_order_raw(p["leg"])
        df_new = _load_new_order(p["new"])
        orders = df_leg["order"].values
        x = np.arange(len(orders))
        w = 0.35
        ax.bar(x - w / 2, df_leg["order_total_time_s"].values / 60, width=w,
               color=C_LEGACY, alpha=0.8, label="Legacy MORFE2.0")
        ax.bar(x + w / 2, df_new["order_total_time_s"].values / 60, width=w,
               color=C_NEW, alpha=0.8, label="New (Ferrite)")
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(orders)
        ax.set_xlabel("Polynomial order")
        ax.set_ylabel("Time per order (min, log)")
        ax.set_title(f"FOM = {p['fom']}  ({p['mesh']})")
        ax.legend(fontsize=8)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:.3g}"))
        ax.grid(axis="y", alpha=0.2)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle("Per-order solve time: Legacy MORFE2.0 vs New (Ferrite)", y=1.01)
    plt.tight_layout()
    _save(fig, "plot_B2_per_order.png")


# ── B3: Scaling log-log — legacy and new ─────────────────────────────────────
def plot_B3_scaling():
    fig, ax = plt.subplots(figsize=(8, 5))

    fom_a, t_leg_a, t_new_a = [], [], []
    for p in pairs:
        df_leg = _load_leg_order_raw(p["leg"])
        df_new = _load_new_order(p["new"])
        fom_a.append(p["fom"])
        t_leg_a.append(df_leg["cumul_time_s"].iloc[-1] / 60)
        t_new_a.append(df_new["cumul_time_s"].iloc[-1] / 60)

    fom_a = np.array(fom_a)
    t_leg = np.array(t_leg_a)
    t_new = np.array(t_new_a)

    def _pfit(fv, tv):
        lf, lt = np.log(fv), np.log(tv)
        alpha, c = np.polyfit(lf, lt, 1)
        C  = np.exp(c)
        r2 = _r2(lt, alpha * lf + c)
        return C, alpha, r2

    C_l, a_l, r2_l = _pfit(fom_a, t_leg)
    C_n, a_n, r2_n = _pfit(fom_a, t_new)
    x_sm = np.linspace(fom_a.min() * 0.9, fom_a.max() * 1.1, 300)

    ax.scatter(fom_a, t_leg, color=C_LEGACY, s=60, zorder=5)
    ax.plot(fom_a, t_leg, color=C_LEGACY, lw=1.5, alpha=0.5)
    ax.plot(x_sm, C_l * x_sm ** a_l, color=C_LEGACY, ls="--", lw=1.8,
            label=f"Legacy  {C_l:.2e}·FOM^{a_l:.2f}  R²={r2_l:.4f}")

    ax.scatter(fom_a, t_new, color=C_NEW, s=60, zorder=5)
    ax.plot(fom_a, t_new, color=C_NEW, lw=1.5, alpha=0.5)
    ax.plot(x_sm, C_n * x_sm ** a_n, color=C_NEW, ls="--", lw=1.8,
            label=f"New     {C_n:.2e}·FOM^{a_n:.2f}  R²={r2_n:.4f}")

    for fom, tl, tn in zip(fom_a, t_leg, t_new):
        ax.annotate(str(fom), (fom, tl), textcoords="offset points",
                    xytext=(6, 4), fontsize=8, color=C_LEGACY)
        ax.annotate(str(fom), (fom, tn), textcoords="offset points",
                    xytext=(6, -10), fontsize=8, color=C_NEW)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("FOM size (DOFs)")
    ax.set_ylabel("Total solve time at order 11 (min)")
    ax.set_title("FOM scaling at order 11: Legacy MORFE2.0 vs New (Ferrite)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.2)
    plt.tight_layout()
    _save(fig, "plot_B3_scaling.png")


# ── B4: Monomial-level cumulative aligned by exp_tuple ───────────────────────
def plot_B4_monomial_cumulative():
    """
    Align legacy and new per-monomial data by exponent tuple.
    Legacy cumulative uses the plot_time_cost_data.py reconstruction
    (monomial_total_time_s.cumsum + order step), read from total_cumul_s column.
    """
    n = len(pairs)
    ncols = min(2, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 5 * nrows),
                              squeeze=False)

    for idx, p in enumerate(pairs):
        ax = axes[idx // ncols][idx % ncols]

        df_leg = _load_leg_mono_with_cumul(p["leg"])
        df_new = _load_new_mono(p["new"])

        merged = (
            df_new[["monomial_idx", "order", "exp_tuple", "cumul_time_s"]]
            .merge(df_leg[["exp_tuple", "total_cumul_s"]], on="exp_tuple", how="inner")
            .sort_values("monomial_idx")
            .reset_index(drop=True)
        )

        if merged.empty:
            ax.set_title(f"FOM={p['fom']}: no shared monomials")
            continue

        x     = merged["monomial_idx"].values.astype(float)
        y_new = merged["cumul_time_s"].values / 60
        y_leg = merged["total_cumul_s"].values / 60

        order_ends = merged.groupby("order").agg(
            x_end=("monomial_idx", "last"),
            ye_new=("cumul_time_s", "last"),
            ye_leg=("total_cumul_s", "last"),
        ).reset_index()
        x_ends  = order_ends["x_end"].values.astype(float)
        ye_new  = order_ends["ye_new"].values / 60
        ye_leg  = order_ends["ye_leg"].values / 60
        orders  = order_ends["order"].values

        pol_leg = Polynomial.fit(x_ends, ye_leg, deg=3)
        pol_new = Polynomial.fit(x_ends, ye_new, deg=3)
        r2_leg  = _r2(ye_leg, pol_leg(x_ends))
        r2_new  = _r2(ye_new, pol_new(x_ends))
        x_sm    = np.linspace(x.min(), x.max(), 400)

        ax.fill_between(x, y_leg, alpha=0.55, color=FILL_LEGACY)
        ax.plot(x, y_leg, color=C_LEGACY, lw=2,
                label=f"Legacy  (R²={r2_leg:.4f})")
        ax.plot(x_sm, pol_leg(x_sm), color=C_LEGACY, ls="--", lw=1.5)
        ax.scatter(x_ends, ye_leg, color=C_LEGACY, s=30, zorder=5)

        ax.fill_between(x, y_new, alpha=0.75, color=FILL_NEW)
        ax.plot(x, y_new, color=C_NEW, lw=2,
                label=f"New (Ferrite)  (R²={r2_new:.4f})")
        ax.plot(x_sm, pol_new(x_sm), color=C_NEW, ls="--", lw=1.5)
        ax.scatter(x_ends, ye_new, color=C_NEW, s=30, zorder=5)

        _order_vlines(ax, x_ends, orders, x.min())

        ax.set_xlabel("Global monomial index")
        ax.set_ylabel("Cumulative time (min)")
        ax.set_title(f"FOM = {p['fom']}  ({p['mesh']})")
        ax.set_xlim(x.min() - 1, x.max() + 1)
        ax.set_ylim(0)
        ax.legend(fontsize=8)
        ax.grid(False)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle(
        "Monomial-level cumulative time (shared solved monomials): Legacy vs New",
        y=1.01,
    )
    plt.tight_layout()
    _save(fig, "plot_B4_monomial_cumulative.png")


# ── B5: RHS vs solve split ────────────────────────────────────────────────────
def plot_B5_rhs_vs_solve():
    n = len(pairs)
    ncols = min(2, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 4.5 * nrows),
                              squeeze=False)

    for idx, p in enumerate(pairs):
        ax = axes[idx // ncols][idx % ncols]
        df_leg = _load_leg_order_raw(p["leg"])
        df_new = _load_new_order(p["new"])
        orders = df_leg["order"].values

        def _frac(rhs, sol):
            tot = rhs + sol
            return (np.where(tot > 0, rhs / tot, 0),
                    np.where(tot > 0, sol / tot, 0))

        lrhs_f, lsol_f = _frac(df_leg["rhs_time_s"].values, df_leg["solve_time_s"].values)
        nrhs_f, nsol_f = _frac(df_new["rhs_time_s"].values, df_new["solve_time_s"].values)

        x = np.arange(len(orders))
        w = 0.35
        ax.bar(x - w/2, lrhs_f, width=w, color=C_LEGACY, alpha=0.9, label="Legacy – RHS")
        ax.bar(x - w/2, lsol_f, width=w, color=C_LEGACY, alpha=0.4,
               bottom=lrhs_f, label="Legacy – Solve")
        ax.bar(x + w/2, nrhs_f, width=w, color=C_NEW, alpha=0.9, label="New – RHS")
        ax.bar(x + w/2, nsol_f, width=w, color=C_NEW, alpha=0.4,
               bottom=nrhs_f, label="New – Solve")

        ax.set_xticks(x)
        ax.set_xticklabels(orders)
        ax.set_xlabel("Polynomial order")
        ax.set_ylabel("Fraction of order time")
        ax.set_title(f"FOM = {p['fom']}  ({p['mesh']})")
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=7, ncol=2)
        ax.grid(axis="y", alpha=0.2)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle("RHS vs solve fraction per order: Legacy vs New", y=1.01)
    plt.tight_layout()
    _save(fig, "plot_B5_rhs_vs_solve_fraction.png")


# ═════════════════════════════════════════════════════════════════════════════
# C) 10x2x2 order-13 run, truncated to order 11: time + memory
# ═════════════════════════════════════════════════════════════════════════════

LEG_10x2x2_O13 = OUTPUT_DIR / "beam_h27_10x2x2_o_13_eps_13_Fmode_1_2026-05-22T18_33_16"
NEW_10x2x2_O13 = NEW_ROOT / "beam_h27_10x2x2_degree13_20260523T174033"


def _load_leg_mono_full(path: Path, max_order: int = MAX_ORDER):
    """
    Load legacy per-monomial CSV (filtered to max_order) and reconstruct:
      total_cumul_s   : time cumulative (monomial_total_time_s.cumsum + order step)
      total_cumul_gb  : allocation cumulative in GB (monomial_total_alloc_bytes.cumsum + order step)
    Returns df with exp_tuple, total_cumul_s, total_cumul_gb columns.
    """
    df_mono = pd.read_csv(path / "benchmark_per_monomial.csv", na_values="NaN")
    df_mono["monomial_total_time_s"]       = df_mono["monomial_total_time_s"].fillna(0)
    df_mono["monomial_total_alloc_bytes"]  = df_mono["monomial_total_alloc_bytes"].fillna(0)
    df_mono = df_mono[df_mono["order"] <= max_order].copy().reset_index(drop=True)

    df_ord = pd.read_csv(path / "benchmark_per_order.csv")
    df_ord = df_ord[df_ord["order"] <= max_order].copy().reset_index(drop=True)
    df_ord["total_time_order"]  = (df_ord["fillrhsG_time_s"]  + df_ord["fillrhsH_time_s"]  + df_ord["fillWf_time_s"])
    df_ord["total_alloc_order"] = (df_ord["fillrhsG_alloc_bytes"] + df_ord["fillrhsH_alloc_bytes"] + df_ord["fillWf_alloc_bytes"])
    cumul_ord_t = df_ord["total_time_order"].cumsum().values
    cumul_ord_a = df_ord["total_alloc_order"].cumsum().values
    df_ord["cumul_monomials"] = df_ord["n_monomials"].cumsum().values
    starts = [1] + (df_ord["cumul_monomials"].iloc[:-1] + 1).tolist()

    step_t = np.zeros(len(df_mono))
    step_a = np.zeros(len(df_mono))
    for i, s in enumerate(starts):
        e = int(df_ord["cumul_monomials"].iloc[i])
        step_t[s - 1 : e] = cumul_ord_t[i]
        step_a[s - 1 : e] = cumul_ord_a[i]

    df_mono["total_cumul_s"]  = df_mono["monomial_total_time_s"].cumsum().values  + step_t
    df_mono["total_cumul_gb"] = (df_mono["monomial_total_alloc_bytes"].cumsum().values + step_a) / (1024**3)
    df_mono["exp_tuple"] = df_mono["alpha_vector"].apply(
        lambda s: tuple(int(x) for x in s.strip("[]").replace(" ", "").split(","))
    )
    return df_mono


def _load_new_mono_full(path: Path, max_order: int = MAX_ORDER):
    """Load new per-monomial CSV (filtered to max_order) with cumul_gb column."""
    df = pd.read_csv(path / "benchmark_per_monomial.csv")
    df = df[df["order"] <= max_order].copy().reset_index(drop=True)
    df["exp_tuple"] = df["exponents"].apply(
        lambda s: tuple(int(x) for x in s.split("_"))
    )
    df["cumul_gb"] = (df["rhs_alloc_bytes"] + df["solve_alloc_bytes"]).cumsum() / (1024**3)
    return df


def _merge_pair(leg_path, new_path, max_order=MAX_ORDER):
    """Align legacy and new by exp_tuple (inner join), return merged DataFrame."""
    df_leg = _load_leg_mono_full(leg_path, max_order)
    df_new = _load_new_mono_full(new_path, max_order)
    merged = (
        df_new[["monomial_idx", "order", "exp_tuple", "cumul_time_s", "cumul_gb"]]
        .merge(
            df_leg[["exp_tuple", "total_cumul_s", "total_cumul_gb"]],
            on="exp_tuple", how="inner",
        )
        .sort_values("monomial_idx")
        .reset_index(drop=True)
    )
    return merged


def _plot_cumul_pair(ax, merged, y_new_col, y_leg_col, ylabel, title):
    """Shared drawing logic for time or memory cumulative comparison."""
    x     = merged["monomial_idx"].values.astype(float)
    y_new = merged[y_new_col].values
    y_leg = merged[y_leg_col].values

    order_ends = merged.groupby("order").agg(
        x_end  = ("monomial_idx",  "last"),
        ye_new = (y_new_col,       "last"),
        ye_leg = (y_leg_col,       "last"),
    ).reset_index()
    x_ends = order_ends["x_end"].values.astype(float)
    ye_new = order_ends["ye_new"].values
    ye_leg = order_ends["ye_leg"].values
    orders = order_ends["order"].values

    pol_leg = Polynomial.fit(x_ends, ye_leg, deg=3)
    pol_new = Polynomial.fit(x_ends, ye_new, deg=3)
    r2_leg  = _r2(ye_leg, pol_leg(x_ends))
    r2_new  = _r2(ye_new, pol_new(x_ends))
    x_sm    = np.linspace(x.min(), x.max(), 400)

    ax.fill_between(x, y_leg, alpha=0.55, color=FILL_LEGACY)
    ax.plot(x, y_leg, color=C_LEGACY, lw=2, label=f"Legacy  (R²={r2_leg:.4f})")
    ax.plot(x_sm, pol_leg(x_sm), color=C_LEGACY, ls="--", lw=1.5)
    ax.scatter(x_ends, ye_leg, color=C_LEGACY, s=30, zorder=5)

    ax.fill_between(x, y_new, alpha=0.75, color=FILL_NEW)
    ax.plot(x, y_new, color=C_NEW, lw=2, label=f"New (Ferrite)  (R²={r2_new:.4f})")
    ax.plot(x_sm, pol_new(x_sm), color=C_NEW, ls="--", lw=1.5)
    ax.scatter(x_ends, ye_new, color=C_NEW, s=30, zorder=5)

    _order_vlines(ax, x_ends, orders, x.min())

    ax.set_xlabel("Global monomial index")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(x.min() - 1, x.max() + 1)
    ax.set_ylim(0)
    ax.legend(fontsize=9)
    ax.grid(False)


def plot_C1_10x2x2_time():
    merged = _merge_pair(LEG_10x2x2_O13, NEW_10x2x2_O13, max_order=13)
    fig, ax = plt.subplots(figsize=(12, 6))
    _plot_cumul_pair(
        ax, merged,
        y_new_col="cumul_time_s",   y_leg_col="total_cumul_s",
        ylabel="Cumulative time (min)",
        title="Cumulative solve time — 10×2×2 (FOM=1425), orders 2–13\nLegacy MORFE2.0 vs New (Ferrite)",
    )
    # data is in seconds; label in minutes (every 60 min), minor ticks every 15 min
    ax.yaxis.set_major_locator(ticker.MultipleLocator(60 * 60))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(15 * 60))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v/60:.0f}"))
    ax.tick_params(axis="y", which="minor", length=3)
    plt.tight_layout()
    _save(fig, "plot_C1_10x2x2_time.png")


def plot_C3_10x2x2_speedup():
    merged = _merge_pair(LEG_10x2x2_O13, NEW_10x2x2_O13, max_order=13)

    x     = merged["monomial_idx"].values.astype(float)
    y_new = merged["cumul_time_s"].values / 60
    y_leg = merged["total_cumul_s"].values / 60

    order_ends = merged.groupby("order").agg(
        x_end  = ("monomial_idx",  "last"),
        ye_new = ("cumul_time_s",  "last"),
        ye_leg = ("total_cumul_s", "last"),
    ).reset_index()
    x_ends = order_ends["x_end"].values.astype(float)
    ye_new = order_ends["ye_new"].values / 60
    ye_leg = order_ends["ye_leg"].values / 60
    orders = order_ends["order"].values

    # Polynomials fitted in monomial-index space (used for cumul time plots)
    pol_leg = Polynomial.fit(x_ends, ye_leg, deg=3)
    pol_new = Polynomial.fit(x_ends, ye_new, deg=3)

    ratio_actual = ye_leg / np.where(ye_new > 0, ye_new, np.nan)

    # Rational function following plot_comparison.py:
    # N(order) = binom(order+4, 4) - 1  (total non-zero monomials up to that order)
    from math import comb
    def n_mon(order):
        return comb(order + 4, 4) - 1

    orders_ext = np.arange(7, 17)
    N_ext      = np.array([n_mon(o) for o in orders_ext], dtype=float)
    denom_ext  = pol_new(N_ext)
    ratio_ext  = np.where(denom_ext > 0, pol_leg(N_ext) / denom_ext, np.nan)

    print("ratio series =", [f"{r:.2f}" for r in ratio_ext])
    print(f"tends to  {pol_leg.coef[3] / pol_new.coef[3]:.4f}  (ratio of leading cubic coefficients)")

    fig, ax = plt.subplots(figsize=(11, 5))

    ax.scatter(orders, ratio_actual, color="mediumseagreen", s=60, zorder=5,
               label="Speedup ratio at order endpoints  (actual)")
    ax.plot(orders, ratio_actual, color="mediumseagreen", lw=1.5, alpha=0.7)
    ax.plot(orders_ext, ratio_ext, color="mediumseagreen", ls="--", lw=2,
            label=r"Rational function  $p_\mathrm{leg}(N) / p_\mathrm{new}(N)$,"
                  r"  $N = \binom{q+4}{4}-1$")

    asymptote = pol_leg.coef[3] / pol_new.coef[3]
    ax.axhline(asymptote, color="darkgreen", ls="--", lw=1.2, alpha=0.7,
               label=f"Asymptote  $c_{{\\mathrm{{leg}}}}/c_{{\\mathrm{{new}}}} = {asymptote:.2f}$")
    ax.axhline(1.0, color="black", ls=":", lw=1.0, alpha=0.5)

    ax.set_xlabel("Polynomial order")
    ax.set_ylabel("Speedup ratio  legacy / new")
    ax.set_title(
        "Speedup ratio — 10×2×2 (FOM=1425), orders 2–13  (extrapolated to 16)\n"
        "Legacy MORFE2.0 vs New (Ferrite)"
    )
    all_orders = np.arange(orders[0], orders_ext[-1] + 1)
    ax.axvspan(orders[-1] + 0.5, orders_ext[-1] + 0.4, color="grey", alpha=0.10, zorder=0)
    ax.axvline(x=orders[-1] + 0.5, color="grey", ls=":", lw=1.0, alpha=0.6)
    ax.set_xticks(all_orders)
    ax.set_xlim(orders[0] - 0.4, orders_ext[-1] + 0.4)
    ax.set_ylim(0)
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    _save(fig, "plot_C3_10x2x2_speedup.png")


def plot_C4_10x2x2_memory_gain():
    from math import comb
    merged = _merge_pair(LEG_10x2x2_O13, NEW_10x2x2_O13, max_order=13)

    x     = merged["monomial_idx"].values.astype(float)

    order_ends = merged.groupby("order").agg(
        x_end  = ("monomial_idx",  "last"),
        ye_new = ("cumul_gb",      "last"),
        ye_leg = ("total_cumul_gb","last"),
    ).reset_index()
    x_ends = order_ends["x_end"].values.astype(float)
    ye_new = order_ends["ye_new"].values
    ye_leg = order_ends["ye_leg"].values
    orders = order_ends["order"].values

    pol_leg = Polynomial.fit(x_ends, ye_leg, deg=3)
    pol_new = Polynomial.fit(x_ends, ye_new, deg=3)

    ratio_actual = ye_leg / np.where(ye_new > 0, ye_new, np.nan)

    def n_mon(order):
        return comb(order + 4, 4) - 1

    orders_ext = np.arange(7, 17)
    N_ext      = np.array([n_mon(o) for o in orders_ext], dtype=float)
    denom_ext  = pol_new(N_ext)
    ratio_ext  = np.where(denom_ext > 0, pol_leg(N_ext) / denom_ext, np.nan)

    print("memory gain series =", [f"{r:.2f}" for r in ratio_ext])
    print(f"tends to  {pol_leg.coef[3] / pol_new.coef[3]:.4f}  (ratio of leading cubic coefficients)")

    fig, ax = plt.subplots(figsize=(11, 5))

    ax.scatter(orders, ratio_actual, color="steelblue", s=60, zorder=5,
               label="Memory gain at order endpoints  (actual)")
    ax.plot(orders, ratio_actual, color="steelblue", lw=1.5, alpha=0.7)
    ax.plot(orders_ext, ratio_ext, color="steelblue", ls="--", lw=2,
            label=r"Rational function  $p_\mathrm{leg}(N) / p_\mathrm{new}(N)$,"
                  r"  $N = \binom{q+4}{4}-1$")

    asymptote = pol_leg.coef[3] / pol_new.coef[3]
    ax.axhline(asymptote, color="navy", ls="--", lw=1.2, alpha=0.7,
               label=f"Asymptote  $c_{{\\mathrm{{leg}}}}/c_{{\\mathrm{{new}}}} = {asymptote:.2f}$")
    ax.axhline(1.0, color="black", ls=":", lw=1.0, alpha=0.5)

    all_orders = np.arange(orders[0], orders_ext[-1] + 1)
    ax.axvspan(orders[-1] + 0.5, orders_ext[-1] + 0.4, color="grey", alpha=0.10, zorder=0)
    ax.axvline(x=orders[-1] + 0.5, color="grey", ls=":", lw=1.0, alpha=0.6)

    ax.set_xlabel("Polynomial order")
    ax.set_ylabel("Memory gain ratio  legacy / new")
    ax.set_title(
        "Memory allocation gain — 10×2×2 (FOM=1425), orders 2–13  (extrapolated to 16)\n"
        "Legacy MORFE2.0 vs New (Ferrite)"
    )
    ax.set_xticks(all_orders)
    ax.set_xlim(orders[0] - 0.4, orders_ext[-1] + 0.4)
    ax.set_ylim(0)
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    _save(fig, "plot_C4_10x2x2_memory_gain.png")


def plot_C2_10x2x2_memory():
    merged = _merge_pair(LEG_10x2x2_O13, NEW_10x2x2_O13, max_order=13)
    fig, ax = plt.subplots(figsize=(12, 6))
    _plot_cumul_pair(
        ax, merged,
        y_new_col="cumul_gb",       y_leg_col="total_cumul_gb",
        ylabel="Cumulative allocation (GB)",
        title="Cumulative memory allocation — 10×2×2 (FOM=1425), orders 2–13\nLegacy MORFE2.0 vs New (Ferrite)",
    )
    plt.tight_layout()
    _save(fig, "plot_C2_10x2x2_memory.png")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("\n── A) Legacy-only plots ──────────────────────────────────────────")
    plot_A1_leg_cumulative()
    plot_A2_leg_per_order()
    plot_A3_leg_scaling()
    plot_A4_leg_memory()

    print("\n── B) Legacy vs New comparison plots ────────────────────────────")
    plot_B1_speedup_ratio()
    plot_B2_per_order()
    plot_B3_scaling()
    plot_B4_monomial_cumulative()
    plot_B5_rhs_vs_solve()

    print("\n── C) 10x2x2 order-13 ───────────────────────────────────────────")
    plot_C1_10x2x2_time()
    plot_C2_10x2x2_memory()
    plot_C3_10x2x2_speedup()
    plot_C4_10x2x2_memory_gain()

    print(f"\nAll figures saved to {PLOTS_DIR}/")
