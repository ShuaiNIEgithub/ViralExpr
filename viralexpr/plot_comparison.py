"""
plot_comparison.py
Generate comparison plots for FPKB, FPKM, TPM, and raw counts.

Usage:
    python3 plot_comparison.py --fpkb expr_FPKB.txt --fpkm expr_FPKM.txt \
        --tpm expr_TPM.txt --counts expr_counts.txt --metadata sample_meta.txt \
        --output plot_output
"""

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)


def parse_args():
    parser = argparse.ArgumentParser(description="Viral gene expression comparison plots")
    parser.add_argument("--fpkb", required=True, help="FPKB matrix")
    parser.add_argument("--fpkm", required=True, help="FPKM matrix")
    parser.add_argument("--tpm", required=True, help="TPM matrix")
    parser.add_argument("--counts", required=True, help="Raw count matrix")
    parser.add_argument("--metadata", required=True, help="Sample metadata TSV")
    parser.add_argument("--output", required=True, help="Output prefix")
    parser.add_argument("--width", type=float, default=10)
    parser.add_argument("--height", type=float, default=7)
    parser.add_argument("--format", default="pdf", choices=["pdf", "png", "svg"])
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def load_matrix(filepath, metric_name):
    df = pd.read_csv(filepath, sep="\t", index_col=0)
    long = df.reset_index().melt(id_vars=[df.index.name or "gene_id"], var_name="sample_id", value_name="value")
    long["metric"] = metric_name
    return long


def load_metadata(filepath):
    meta = pd.read_csv(filepath, sep="\t")
    if "sample_id" not in meta.columns or "timepoint" not in meta.columns:
        raise ValueError("Metadata must contain 'sample_id' and 'timepoint' columns")
    return meta


def calc_pvalues(data):
    genes = data["gene_id"].unique()
    timepoints = data["timepoint"].unique()
    if len(timepoints) != 2:
        return pd.DataFrame({"gene_id": genes, "p_value": np.nan, "significant": False})
    t1, t2 = timepoints[0], timepoints[1]
    results = []
    for gene in genes:
        sub = data[data["gene_id"] == gene]
        g1 = sub[sub["timepoint"] == t1]["value"].dropna().values
        g2 = sub[sub["timepoint"] == t2]["value"].dropna().values
        if len(g1) >= 2 and len(g2) >= 2:
            try:
                _, pval = stats.ttest_ind(g1, g2)
                results.append({"gene_id": gene, "p_value": pval, "significant": pval < 0.05})
                continue
            except Exception:
                pass
        results.append({"gene_id": gene, "p_value": np.nan, "significant": False})
    return pd.DataFrame(results)


def plot_comparison(all_data, args):
    metric_order = ["Raw count", "FPKM", "TPM", "FPKB"]
    all_data["metric"] = pd.Categorical(all_data["metric"], categories=metric_order, ordered=True)
    summary = all_data.groupby(["gene_id", "timepoint", "metric"], observed=True)["value"].mean().reset_index()
    g = sns.catplot(data=summary, x="gene_id", y="value", hue="timepoint",
                    col="metric", kind="bar", col_wrap=1,
                    height=args.height * 0.35, aspect=args.width / (args.height * 0.35),
                    palette="Set1", legend_out=True)
    g.set_axis_labels("Viral gene", "Gene expression")
    g.set_titles(col_template="{col_name}")
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=0)
        ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    return g.fig


def plot_fpkb_detailed(fpkb_data, pvals, args):
    summary = fpkb_data.groupby(["gene_id", "timepoint"]).agg(
        mean=("value", "mean"),
        se=("value", lambda x: np.std(x, ddof=1) / np.sqrt(len(x)))
    ).reset_index()
    summary["ymin"] = summary["mean"] - summary["se"]
    summary["ymax"] = summary["mean"] + summary["se"]
    sig_genes = set(pvals[pvals["significant"]]["gene_id"])

    fig, ax = plt.subplots(figsize=(args.width, args.height))
    bar = sns.barplot(data=summary, x="gene_id", y="mean", hue="timepoint",
                      palette="Set1", ax=ax, width=0.7)
    for i, tp in enumerate(summary["timepoint"].unique()):
        sub = summary[summary["timepoint"] == tp].reset_index(drop=True)
        x_offsets = np.arange(len(sub)) + (i - 0.5) * 0.33
        ax.errorbar(x=x_offsets, y=sub["mean"],
                    yerr=[sub["mean"] - sub["ymin"], sub["ymax"] - sub["mean"]],
                    fmt="none", c="black", capsize=3, capthick=1, linewidth=1)
    if sig_genes:
        y_max = summary["ymax"].max()
        for gene in sig_genes:
            try:
                gene_idx = list(summary["gene_id"].unique()).index(gene)
                ax.text(gene_idx, y_max * 1.05, "*", ha="center", va="bottom", fontsize=16, color="red")
            except ValueError:
                pass
    ax.set_xlabel("Viral gene", fontsize=12)
    ax.set_ylabel("FPKB", fontsize=12)
    ax.set_title("Viral gene expression (FPKB normalization)", fontsize=13, fontweight="bold")
    subtitle = f"Significant: {', '.join(sorted(sig_genes))}" if sig_genes else "No significant differences"
    ax.text(0.5, 1.02, subtitle, transform=ax.transAxes, ha="center", va="bottom", fontsize=10, style="italic")
    ax.legend(title="Timepoint", loc="upper right")
    ax.grid(axis="y", alpha=0.3)
    sns.despine(ax=ax)
    plt.tight_layout()
    return fig


def plot_distribution(fpkb_data, args):
    n_genes = fpkb_data["gene_id"].nunique()
    ncols = min(5, n_genes)
    nrows = int(np.ceil(n_genes / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=(args.width * 1.2, args.height * nrows * 0.6),
                             sharey=False)
    axes = np.atleast_1d(axes).flatten()
    genes = sorted(fpkb_data["gene_id"].unique())
    palette = sns.color_palette("Set1", n_colors=fpkb_data["timepoint"].nunique())
    for idx, gene in enumerate(genes):
        ax = axes[idx]
        sub = fpkb_data[fpkb_data["gene_id"] == gene]
        sns.violinplot(data=sub, x="timepoint", y="value", palette=palette, ax=ax, inner=None, alpha=0.6)
        sns.boxplot(data=sub, x="timepoint", y="value", color="white", ax=ax, width=0.2, showcaps=True, showfliers=False)
        sns.stripplot(data=sub, x="timepoint", y="value", color="black", ax=ax, size=4, alpha=0.5, jitter=True)
        ax.set_title(gene, fontsize=11, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("FPKB" if idx % ncols == 0 else "")
        ax.grid(axis="y", alpha=0.3)
    for idx in range(len(genes), len(axes)):
        axes[idx].set_visible(False)
    fig.suptitle("FPKB distribution by gene and timepoint", fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    return fig


def main():
    args = parse_args()
    out_prefix = Path(args.output)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    fpkb = load_matrix(args.fpkb, "FPKB")
    fpkm = load_matrix(args.fpkm, "FPKM")
    tpm = load_matrix(args.tpm, "TPM")
    raw = load_matrix(args.counts, "Raw count")
    meta = load_metadata(args.metadata)

    all_data = pd.concat([fpkb, fpkm, tpm, raw], ignore_index=True)
    all_data = all_data.merge(meta, on="sample_id", how="left")

    fpkb_data = all_data[all_data["metric"] == "FPKB"].copy()
    zero_samples = fpkb_data.groupby("sample_id")["value"].sum()
    zero_ids = zero_samples[zero_samples == 0].index.tolist()
    if zero_ids:
        print(f"Warning: zero expression samples: {', '.join(zero_ids)}", file=sys.stderr)

    pvals = calc_pvalues(fpkb_data)

    fig1 = plot_comparison(all_data, args)
    fig1.savefig(out_prefix.parent / f"{out_prefix.name}_comparison.{args.format}", dpi=args.dpi, bbox_inches="tight", format=args.format)
    plt.close(fig1)

    fig2 = plot_fpkb_detailed(fpkb_data, pvals, args)
    fig2.savefig(out_prefix.parent / f"{out_prefix.name}_FPKB_detailed.{args.format}", dpi=args.dpi, bbox_inches="tight", format=args.format)
    plt.close(fig2)

    fig3 = plot_distribution(fpkb_data, args)
    fig3.savefig(out_prefix.parent / f"{out_prefix.name}_distribution.{args.format}", dpi=args.dpi, bbox_inches="tight", format=args.format)
    plt.close(fig3)

    pvals.to_csv(out_prefix.parent / f"{out_prefix.name}_statistics.txt", sep="\t", index=False)
    print("Plots saved to:", out_prefix.parent)


if __name__ == "__main__":
    main()