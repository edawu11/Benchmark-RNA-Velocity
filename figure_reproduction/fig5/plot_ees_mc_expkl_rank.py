from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MATCH_COLOR = "#2563EB"
MISMATCH_COLOR = "#DC2626"
DATA_DIR = Path(__file__).resolve().parent
METHOD_RENAME_DICT = {
    "velocyto": "Velocyto",
    "scvelo_dyn": "scVelo (dyn)",
    "scvelo_stc": "scVelo (stc)",
    "veloae": "VeloAE",
    "dynamo_m1": "Dynamo (m1)",
    "unitvelo_ind": "UniTVelo (ind)",
    "unitvelo_uni": "UniTVelo (uni)",
    "velovae_std": "VeloVAE (std)",
    "velovae_fullvb": "VeloVAE (Full VB)",
    "deepvelo": "DeepVelo",
    "celldancer": "cellDancer",
    "pyrovelocity_m1": "Pyro-Velocity (m1)",
    "pyrovelocity_m2": "Pyro-Velocity (m2)",
    "kvelo": "  -velo",
    "velovi": "veloVI",
    "cell2fate": "cell2fate",
    "sdevelo": "SDEvelo",
    "svelvetvae": "SvelvetVAE",
    "latentvelo_std": "LatentVelo (std)",
    "tivelo_std": "TIVelo (std)",
    "tivelo_simple": "TIVelo (simple)",
    "sctour_mse": "scTour (MSE)",
    "sctour_nb": "scTour (NB)",
    "sctour_zinb": "scTour (ZINB)",
    "multivelo": "MultiVelo",
    "velvetvae": "VelvetVAE",
    "graphvelo_std": "GraphVelo (std)",
}


def dataset_sort_key(path: Path) -> tuple[int, str]:
    suffix = path.name.replace("data", "")
    return (int(suffix) if suffix.isdigit() else 10**9, path.name)


def discover_datasets(base_dir: Path) -> list[Path]:
    datasets = []
    for path in base_dir.glob("data*"):
        if not path.is_dir():
            continue
        if (path / "EES_df.csv").exists() and (path / "trim_nb_num_df.csv").exists():
            datasets.append(path)
    return sorted(datasets, key=dataset_sort_key)


def read_mean_and_rank(csv_path: Path, value_col: str, rank_col: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    fold_cols = [col for col in df.columns if col != "Method"]
    metric_df = df[["Method"]].copy()
    metric_df[value_col] = df[fold_cols].mean(axis=1)
    metric_df[rank_col] = (
        metric_df[value_col]
        .rank(ascending=False, method="min")
        .astype(int)
    )
    return metric_df


def build_rank_table(dataset_dir: Path) -> pd.DataFrame:
    ees_df = read_mean_and_rank(dataset_dir / "EES_df.csv", "EES_mean", "EES_rank")
    mc_df = read_mean_and_rank(dataset_dir / "trim_nb_num_df.csv", "mc_mean", "mc_rank")

    rank_df = ees_df.merge(mc_df, on="Method", how="inner")
    rank_df["expKL_mean"] = rank_df["EES_mean"] / rank_df["mc_mean"]
    rank_df["expKL_rank"] = (
        rank_df["expKL_mean"]
        .rank(ascending=False, method="min")
        .astype(int)
    )
    rank_df["mc_matches_ees"] = rank_df["mc_rank"] == rank_df["EES_rank"]
    rank_df["point_color"] = np.where(rank_df["mc_matches_ees"], MATCH_COLOR, MISMATCH_COLOR)
    rank_df["Method_display"] = rank_df["Method"].map(METHOD_RENAME_DICT).fillna(rank_df["Method"])
    rank_df = rank_df.sort_values(["EES_rank", "Method"]).reset_index(drop=True)
    return rank_df


def build_xticks(max_rank: int) -> list[int]:
    tick_step = 5 if max_rank > 10 else 1
    positive_ticks = list(range(tick_step, max_rank + 1, tick_step))
    return [-tick for tick in reversed(positive_ticks)] + [0] + positive_ticks


def plot_dataset(rank_df: pd.DataFrame, output_path: Path) -> None:
    max_rank = int(rank_df[["EES_rank", "mc_rank", "expKL_rank"]].to_numpy().max())
    y_values = rank_df["EES_rank"]

    fig_height = max(8, len(rank_df) * 0.38)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    ax.axvline(0, color="#9CA3AF", linestyle="--", linewidth=1.2, zorder=1)

    for row in rank_df.itertuples():
        ax.hlines(
            y=row.EES_rank,
            xmin=-row.mc_rank,
            xmax=0,
            color=row.point_color,
            linewidth=1.2,
            linestyle=":",
            alpha=0.5,
            zorder=2,
        )
        ax.hlines(
            y=row.EES_rank,
            xmin=0,
            xmax=row.expKL_rank,
            color=row.point_color,
            linewidth=1.2,
            linestyle=":",
            alpha=0.5,
            zorder=2,
        )

    ax.scatter(
        -rank_df["mc_rank"],
        y_values,
        s=90,
        c=rank_df["point_color"],
        edgecolors="white",
        linewidths=0.8,
        zorder=3,
    )
    ax.scatter(
        rank_df["expKL_rank"],
        y_values,
        s=90,
        c=rank_df["point_color"],
        edgecolors="white",
        linewidths=0.8,
        zorder=3,
    )

    ax.set_xlim(-(max_rank + 1), max_rank + 1)
    ax.set_ylim(max_rank + 0.5, 0.5)

    xticks = build_xticks(max_rank)
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(abs(tick)) for tick in xticks], fontsize=10)

    ax.set_yticks(rank_df["EES_rank"])
    ax.set_yticklabels(
        [f"{row.Method_display} ({row.EES_rank})" for row in rank_df.itertuples()],
        fontsize=9,
    )

    ax.tick_params(axis="x", which="major", length=4, width=1.0, color="black", labelsize=10)
    ax.tick_params(axis="y", which="major", length=4, width=1.0, color="black", labelsize=9)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_match_summary(summary_tables: list[pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for df in summary_tables:
        dataset_name = df["dataset"].iloc[0]
        same_count = int(df["mc_matches_ees"].sum())
        total_count = int(len(df))
        different_count = total_count - same_count
        rows.append(
            {
                "dataset": dataset_name,
                "same_count": same_count,
                "different_count": different_count,
                "total_count": total_count,
            }
        )
    return pd.DataFrame(rows).sort_values("dataset").reset_index(drop=True)


def format_dataset_label(dataset_name: str) -> str:
    if dataset_name.startswith("data") and dataset_name[4:].isdigit():
        return f"Data {dataset_name[4:]}"
    return dataset_name


def plot_match_summary(match_summary_df: pd.DataFrame, output_path: Path) -> None:
    axis_limit = 20
    y_pos = np.arange(len(match_summary_df))

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.axvline(0, color="#9CA3AF", linestyle="--", linewidth=1.2, zorder=1)

    ax.barh(
        y_pos,
        -match_summary_df["same_count"],
        color=MATCH_COLOR,
        height=0.62,
        zorder=2,
    )
    ax.barh(
        y_pos,
        match_summary_df["different_count"],
        color=MISMATCH_COLOR,
        height=0.62,
        zorder=2,
    )

    for idx, row in match_summary_df.iterrows():
        ax.text(
            -row["same_count"] / 2,
            idx,
            str(int(row["same_count"])),
            ha="center",
            va="center",
            fontsize=10,
            color="black",
            zorder=3,
        )
        ax.text(
            row["different_count"] / 2,
            idx,
            str(int(row["different_count"])),
            ha="center",
            va="center",
            fontsize=10,
            color="black",
            zorder=3,
        )

    ax.set_xlim(-axis_limit, axis_limit)
    xticks = build_xticks(axis_limit)
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(abs(tick)) for tick in xticks], fontsize=10)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([format_dataset_label(name) for name in match_summary_df["dataset"]], fontsize=11)
    ax.invert_yaxis()

    ax.tick_params(axis="x", which="major", length=4, width=1.0, color="black", labelsize=10)
    ax.tick_params(axis="y", which="major", length=4, width=1.0, color="black", labelsize=11)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    ax.grid(False)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    datasets = discover_datasets(DATA_DIR)
    if not datasets:
        raise FileNotFoundError("No dataset directory with both EES_df.csv and trim_nb_num_df.csv was found.")

    summary_tables = []
    for dataset_dir in datasets:
        rank_df = build_rank_table(dataset_dir)
        summary_tables.append(rank_df.assign(dataset=dataset_dir.name))

        csv_output = DATA_DIR / f"{dataset_dir.name}_EES_mc_expKL_rank.csv"
        png_output = DATA_DIR / f"{dataset_dir.name}_EES_mc_expKL_rank.png"
        rank_df.to_csv(csv_output, index=False)
        plot_dataset(rank_df, png_output)

    pd.concat(summary_tables, ignore_index=True).to_csv(
        DATA_DIR / "EES_mc_expKL_rank_summary.csv",
        index=False,
    )
    match_summary_df = build_match_summary(summary_tables)
    match_summary_df.to_csv(DATA_DIR / "EES_sparsity_rank_match_summary.csv", index=False)
    plot_match_summary(match_summary_df, DATA_DIR / "EES_sparsity_rank_match_bar.png")


if __name__ == "__main__":
    main()
