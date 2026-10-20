import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.stats import median_abs_deviation

GROUPS= ["Not_AD", "Low", "Intermediate", "High"]
CONTRASTS= [("Not_AD", "Low"), ("Low", "Intermediate"), ("Intermediate", "High"), ("Not_AD", "High")]

BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
OUTFIGDIR = Path(f"{BASE}/output/protein_coding/figures")
OUTDATADIR = Path(f"{BASE}/output/protein_coding/data")
INDATADIR_GENE= Path(f"{OUTDATADIR}/gene_analysis_changes")
OUTDATADIR_GENE= Path(f"{OUTDATADIR}/gene_analysis_changes")

GENE_METRICS = pd.read_csv(INDATADIR_GENE / "gene_metrics_protein_coding.csv")

GENE_STATUS = pd.read_csv(INDATADIR_GENE / "in_out_genes_protein_coding.csv")

GENE_ROLES = pd.read_csv(INDATADIR_GENE / "gene_roles_protein_coding.csv")

STATUS_COLORS = ["#a6cee3", "#1f78b4", "#b2df8a"]

CONTRAST_COLORS = ("#2c7bb6", "#fdae61", "#d7191c", "#abd9e9")
CONTRAST_ORDER = [f"{g1} → {g2}" for g1, g2 in CONTRASTS]
STAGE_COLORS = dict(zip(GROUPS, CONTRAST_COLORS))

QUADRANT_COLORS = {
    "More central, more inter-modular": "#d7191c",
    "More central, more intra-modular": "#fdae61",
    "Less central, more inter-modular": "#abd9e9",
    "Less central, more intra-modular": "#2c7bb6",
}


def plot_status_distribution_contrasts(gene_status, tag="protein_coding"):

    df = gene_status.copy()
    df["contrast"] = df["group1"] + " → " + df["group2"]

    contrast_order = [f"{g1} → {g2}" for g1, g2 in CONTRASTS]
    status_order = ["in", "out", "in_both"]
    status_labels = {
        "in": "In (only in group2)",
        "out": "Out (only in group1)",
        "in_both": "In both networks",
    }
    df["status"] = df["status"].map(status_labels)

    fig, ax = plt.subplots(figsize=(10, 5))

    sns.barplot(
        data=df,
        x="contrast",
        y="count",
        hue="status",
        order=contrast_order,
        hue_order=[status_labels[s] for s in status_order],
        palette=STATUS_COLORS,
        edgecolor="black",
        linewidth=0.5,
        ax=ax,
    )

    for container in ax.containers:
        ax.bar_label(container, fontsize=8, padding=2)

    ax.set_xlabel("Contrast")
    ax.set_ylabel("Number of genes")
    ax.set_title("Gene presence in co-expression networks across stages")
    ax.legend(title="Gene status", frameon=False)
    sns.despine(ax=ax)

    fig.tight_layout()
    fig.savefig(OUTFIGDIR / f"gene_status_distribution_{tag}.png", dpi=300)
    fig.savefig(OUTFIGDIR / f"gene_status_distribution_{tag}.pdf")
    plt.close(fig)

    return fig, ax

def neighbor_retention_table(gene_metrics, min_deg=5):
    """Per-gene neighbor retention, expected retention and excess over chance."""
    df = gene_metrics.copy()
    df["contrast"] = df["group1"] + " → " + df["group2"]

    df = df[
        (df["status"] == "in_both")
        & (df["degree1_shared"] >= min_deg)
        & (df["degree2_shared"] >= min_deg)
    ].copy()

    df["pct_conserved"] = df["n_conserved"] / df["degree1_shared"] * 100
    df["pct_expected"] = df["degree2_shared"] / (df["n_shared"] - 1) * 100
    df["excess"] = df["pct_conserved"] - df["pct_expected"]

    return df


def plot_conserved_neighbors(gene_metrics, min_deg=5, tag="protein_coding"):

    df = neighbor_retention_table(gene_metrics, min_deg)
    colors = dict(zip(CONTRAST_ORDER, CONTRAST_COLORS))
    bins = np.arange(-10, 102.5, 2.5)          # same bins in every panel

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, sharey=True)

    for ax, contrast in zip(axes.flat, CONTRAST_ORDER):
        values = df.loc[df["contrast"] == contrast, "excess"]
        median = values.median()

        ax.hist(values, bins=bins, color=colors[contrast],
                edgecolor="white", linewidth=0.5)
        ax.axvline(0, ls="--", c="0.45", lw=1)
        ax.axvline(median, ls="-", c="black", lw=1.2)

        ax.text(0.97, 0.95, f"n = {len(values)}\nmedian = {median:.1f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=9)
        sns.despine(ax=ax)

    fig.supxlabel("Neighbor retention above chance (percentage points)")
    fig.supylabel("Number of genes")
    fig.suptitle("Neighborhood conservation across stages", fontweight="bold")

    handles = [Patch(facecolor=colors[c], label=c) for c in CONTRAST_ORDER]
    handles += [
        Line2D([0], [0], ls="--", c="0.45", lw=1, label="Chance level (0)"),
        Line2D([0], [0], ls="-", c="black", lw=1.2, label="Median"),
    ]

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
               frameon=False, title="Contrast", alignment="left")

    fig.savefig(OUTFIGDIR / f"neighbor_retention_hist_{tag}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTFIGDIR / f"neighbor_retention_hist_{tag}.pdf", bbox_inches="tight")
    plt.close(fig)

    return fig

def plot_zp_cartography(gene_roles, n_labels=5, tag="protein_coding"):
    """z–P plane per stage, labeling the genes with the highest within-module degree."""

    df = gene_roles.dropna(subset=["z", "P"]).copy()

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, sharey=True)
    ymin, ymax = df["z"].min() - 0.3, df["z"].max() + 0.5

    for ax, group in zip(axes.flat, GROUPS):
        sub = df[df["group"] == group]

        ax.scatter(sub["P"], sub["z"], s=10, alpha=0.6,
                   color=STAGE_COLORS[group], edgecolor="none")

        for _, row in sub.nlargest(n_labels, "z").iterrows():
            ax.annotate(row["gene"], (row["P"], row["z"]), fontsize=7,
                        xytext=(3, 3), textcoords="offset points")

        ax.text(0.97, 0.05, f"n = {len(sub)}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=9)

        ax.set_xlim(-0.02, 1)
        ax.set_ylim(ymin, ymax)
        sns.despine(ax=ax)

    fig.supxlabel("Participation coefficient (P)")
    fig.supylabel("Within-module degree (z)")
    fig.suptitle("Gene roles across stages", fontweight="bold")

    handles = [Patch(facecolor=STAGE_COLORS[g], label=g) for g in GROUPS]

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
               frameon=False, title="Stage", alignment="left")

    fig.savefig(OUTFIGDIR / f"zp_cartography_{tag}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTFIGDIR / f"zp_cartography_{tag}.pdf", bbox_inches="tight")
    plt.close(fig)

def role_change_outliers(gene_metrics, min_deg=5, threshold=3, tag="protein_coding"):
    """Robust outliers of (ΔP, Δz) within each contrast."""

    df = gene_metrics[
        (gene_metrics["status"] == "in_both")
        & (gene_metrics["degree1"] >= min_deg)
        & (gene_metrics["degree2"] >= min_deg)
    ].dropna(subset=["z1", "z2", "P1", "P2"]).copy()

    df["contrast"] = df["group1"] + " → " + df["group2"]
    df["dP"] = df["P2"] - df["P1"]
    df["dz"] = df["z2"] - df["z1"]

    def robust_z(x):
        return (x - x.median()) / median_abs_deviation(x, scale="normal")

    df["rz_dP"] = df.groupby("contrast")["dP"].transform(robust_z)
    df["rz_dz"] = df.groupby("contrast")["dz"].transform(robust_z)
    df["distance"] = np.hypot(df["rz_dP"], df["rz_dz"])
    df["outlier"] = df["distance"] >= threshold

    df["quadrant"] = np.select(
        [(df.dz >= 0) & (df.dP >= 0), (df.dz >= 0) & (df.dP < 0), (df.dz < 0) & (df.dP >= 0)],
        list(QUADRANT_COLORS)[:3],
        default=list(QUADRANT_COLORS)[3],
    )

    out = df[df["outlier"]].copy()
    out["n_contrasts_outlier"] = out.groupby("gene")["contrast"].transform("nunique")

    cols = ["gene", "contrast", "quadrant", "distance", "dP", "dz",
            "P1", "P2", "z1", "z2", "degree1", "degree2", "n_conserved",
            "jaccard", "n_contrasts_outlier"]
    out = out[cols].sort_values(["contrast", "distance"], ascending=[True, False])
    out.to_csv(OUTDATADIR_GENE / f"role_change_outliers_{tag}.csv", index=False)

    return df, out


def plot_role_change_outliers(df, tag="protein_coding"):

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, sharey=True)

    for ax, contrast in zip(axes.flat, CONTRAST_ORDER):
        sub = df[df["contrast"] == contrast]
        bg = sub[~sub["outlier"]]
        hl = sub[sub["outlier"]]

        ax.scatter(bg["dP"], bg["dz"], s=8, color="0.8", edgecolor="none")
        ax.scatter(hl["dP"], hl["dz"], s=18, edgecolor="none",
                   color=hl["quadrant"].map(QUADRANT_COLORS))

        for _, row in hl.iterrows():
            ax.annotate(row["gene"], (row["dP"], row["dz"]), fontsize=6.5,
                        xytext=(3, 3), textcoords="offset points")

        ax.axhline(0, c="0.45", lw=0.8)
        ax.axvline(0, c="0.45", lw=0.8)
        ax.text(0.97, 0.05, f"outliers = {len(hl)} / {len(sub)}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9)
        ax.set_title(contrast, fontsize=10)
        sns.despine(ax=ax)

    fig.supxlabel("ΔP (participation, group2 − group1)")
    fig.supylabel("Δz (within-module degree, group2 − group1)")
    fig.suptitle("Genes with extreme role changes", fontweight="bold")

    handles = [Line2D([0], [0], marker="o", ls="", markersize=6,
                      markerfacecolor=col, markeredgecolor="none", label=q)
               for q, col in QUADRANT_COLORS.items()]
    handles.append(Line2D([0], [0], marker="o", ls="", markersize=6,
                          markerfacecolor="0.8", markeredgecolor="none",
                          label="Not an outlier"))

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
               frameon=False, title="Role change", alignment="left")

    fig.savefig(OUTFIGDIR / f"role_change_outliers_{tag}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTFIGDIR / f"role_change_outliers_{tag}.pdf", bbox_inches="tight")
    plt.close(fig)

def select_trajectory_genes(outliers, min_contrasts=2, max_genes=12):
    """Genes that are role-change outliers in several contrasts, strongest first."""

    rank = (
        outliers.groupby("gene")
                .agg(n_contrasts=("contrast", "nunique"), max_distance=("distance", "max"))
                .sort_values(["n_contrasts", "max_distance"], ascending=False)
    )
    recurrent = rank[rank["n_contrasts"] >= min_contrasts]
    chosen = recurrent if len(recurrent) else rank      # fallback: strongest single outliers
    return chosen.head(max_genes).index.tolist()


def plot_zp_trajectories(gene_roles, genes, ncols=4, tag="protein_coding"):
    """Trajectory of each gene across stages in the z-P plane."""

    df = gene_roles.dropna(subset=["z", "P"])
    nrows = int(np.ceil(len(genes) / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 2.8 * nrows),
                             sharex=True, sharey=True, squeeze=False)
    ymin, ymax = df["z"].min() - 0.3, df["z"].max() + 0.3

    for ax, gene in zip(axes.flat, genes):
        ax.scatter(df["P"], df["z"], s=4, color="0.88", edgecolor="none", zorder=1)

        path = (df[df["gene"] == gene]
                .set_index("group")
                .reindex(GROUPS)
                .dropna(subset=["z", "P"]))

        for (g1, r1), (g2, r2) in zip(path.iterrows(), path.iloc[1:].iterrows()):
            ax.annotate("", xy=(r2["P"], r2["z"]), xytext=(r1["P"], r1["z"]),
                        arrowprops=dict(arrowstyle="->", color="0.3", lw=1.1,
                                        shrinkA=4, shrinkB=4), zorder=2)

        for group, row in path.iterrows():
            ax.scatter(row["P"], row["z"], s=40, color=STAGE_COLORS[group],
                       edgecolor="black", linewidth=0.4, zorder=3)

        ax.set_title(gene, fontsize=10, fontweight="bold")
        ax.set_xlim(-0.02, 1)
        ax.set_ylim(ymin, ymax)
        sns.despine(ax=ax)

    for ax in axes.flat[len(genes):]:
        ax.axis("off")

    fig.supxlabel("Participation coefficient (P)")
    fig.supylabel("Within-module degree (z)")
    fig.suptitle("Role trajectories across stages", fontweight="bold")

    handles = [Line2D([0], [0], marker="o", ls="", markersize=7,
                      markerfacecolor=STAGE_COLORS[g], markeredgecolor="black",
                      markeredgewidth=0.4, label=g) for g in GROUPS]
    handles.append(Line2D([0], [0], marker="o", ls="", markersize=5,
                          markerfacecolor="0.88", markeredgecolor="none",
                          label="All genes, all stages"))

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
               frameon=False, title="Stage", alignment="left")

    fig.savefig(OUTFIGDIR / f"zp_trajectories_{tag}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTFIGDIR / f"zp_trajectories_{tag}.pdf", bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    plot_status_distribution_contrasts(GENE_STATUS)
    plot_conserved_neighbors(GENE_METRICS)
    plot_zp_cartography(GENE_ROLES)

    role_df, outliers = role_change_outliers(GENE_METRICS)
    plot_role_change_outliers(role_df)
    print(outliers.groupby("contrast").size())

    trajectory_genes = select_trajectory_genes(outliers)
    print(trajectory_genes)
    plot_zp_trajectories(GENE_ROLES, trajectory_genes)