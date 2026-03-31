import pandas as pd
import glob
import os
from collections import defaultdict
import itertools
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------
# CONFIG
# -----------------------

ANNOT_FILE = "outputs/annotations_filtered_full.tsv" # annotations regroupées
CLUSTERS_DIR = "outputs/clusters_output"
METRICS_FILE = "outputs/metrics_output/jaccard_metrics.csv"
PLOTS_DIR = "outputs/metrics_output/plots_jaccard_final"

os.makedirs("outputs/metrics_output", exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

VALID_SOURCES = {
    "Pfam", "PANTHER", "Gene3D", "SMART",
    "SUPERFAMILY", "CDD", "PRINTS", "ProSiteProfiles", "ProSitePatterns"
}

MIN_CLUSTER_SIZE = 5  # pour filtrer petits clusters

# -----------------------
# LOAD ANNOTATIONS
# -----------------------

print("Loading annotations...")
annot_df = pd.read_csv(ANNOT_FILE, sep="\t")
annot_df = annot_df[annot_df["source"].isin(VALID_SOURCES)]

annotations = defaultdict(lambda: defaultdict(set))
for _, row in annot_df.iterrows():
    annotations[row["prot_id"]][row["source"]].add(row["annot_id"])

print("Annotated proteins:", len(annotations))

# -----------------------
# READ CLUSTERS
# -----------------------

def read_clusters(file):
    clusters = defaultdict(list)
    with open(file) as f:
        next(f)  # sauter le header
        for line in f:
            prot, cid = line.strip().split("\t")
            clusters[cid].append(prot)
    return clusters


def add_median_labels(ax, df, x_col, y_col, order):
    medians = df.groupby(x_col)[y_col].median()
    for i, val in enumerate(order):
        if val not in medians.index:
            continue
        med = medians.loc[val]
        ax.text(
            i,
            med,
            f"{med:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
            color="black",
            fontweight="bold",
        )


def plot_violin_box(df, x_col, y_col, order, title, xlabel, output_file):
    if df.empty:
        return

    sns.set(style="whitegrid")
    plt.figure(figsize=(11, 6))
    ax = sns.violinplot(data=df, x=x_col, y=y_col, order=order, inner=None)
    sns.boxplot(
        data=df,
        x=x_col,
        y=y_col,
        order=order,
        width=0.2,
        showcaps=True,
        boxprops={"facecolor": "white", "zorder": 3},
        showfliers=False,
    )

    add_median_labels(ax, df, x_col, y_col, order)

    plt.title(title, fontsize=13)
    plt.xlabel(xlabel, fontsize=11)
    plt.ylabel("Jaccard mean", fontsize=11)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

# -----------------------
# JACCARD FUNCTION
# -----------------------

def jaccard(a, b):
    if not a and not b:
        return None
    return len(a & b) / len(a | b) if (a | b) else None

# -----------------------
# CALCUL DES MÉTRIQUES
# -----------------------

results = []

cluster_files = sorted(glob.glob(CLUSTERS_DIR + "/*.tsv"))

for cfile in cluster_files:
    print("Processing", cfile)
    clusters = read_clusters(cfile)
    method = os.path.basename(cfile)
    
    for cid, prots in clusters.items():
        if len(prots) < 2:
            continue  # pas de paires

        for source in VALID_SOURCES:
            values = []
            for p1, p2 in itertools.combinations(prots, 2):
                a = annotations[p1][source] if p1 in annotations else set()
                b = annotations[p2][source] if p2 in annotations else set()
                if not a or not b:
                    continue
                j = jaccard(a, b)
                if j is not None:
                    values.append(j)
            
            if not values:
                continue

            results.append({
                "cluster_file": method,
                "cluster_id": cid,
                "annotation_type": source,
                "cluster_size": len(prots),
                "jaccard_mean": np.mean(values)
            })

df = pd.DataFrame(results)

if df.empty:
    raise RuntimeError("No Jaccard values were computed. Check annotations and clusters.")

df["r"] = df["cluster_file"].str.extract(r"r([0-9.]+)").astype(float)
df["e"] = df["cluster_file"].str.extract(r"e([0-9.]+)").astype(float)
df = df.sort_values(["annotation_type", "e", "r", "cluster_id"])

df.to_csv(METRICS_FILE, index=False)
print("Done — metrics saved:", len(df))

# -----------------------
# PLOTS
# -----------------------

for annot in df["annotation_type"].unique():
    sub = df[df["annotation_type"] == annot]

    # Figure 1: e=2, tous r
    sub_e2 = sub[sub["e"] == 2]
    r_values = sorted(sub_e2["r"].dropna().unique())
    plot_violin_box(
        sub_e2,
        x_col="r",
        y_col="jaccard_mean",
        order=r_values,
        title=f"Jaccard intra-cluster - {annot} - all clusters, e=2",
        xlabel="Inflation r",
        output_file=f"{PLOTS_DIR}/jaccard_{annot}_all_e2.png",
    )

    # Figure 2: e=2, clusters >= MIN_CLUSTER_SIZE
    sub_e2_5 = sub_e2[sub_e2["cluster_size"] >= MIN_CLUSTER_SIZE]
    plot_violin_box(
        sub_e2_5,
        x_col="r",
        y_col="jaccard_mean",
        order=r_values,
        title=f"Jaccard intra-cluster - {annot} - clusters >= {MIN_CLUSTER_SIZE}, e=2",
        xlabel="Inflation r",
        output_file=f"{PLOTS_DIR}/jaccard_{annot}_ge{MIN_CLUSTER_SIZE}_e2.png",
    )

    # Figure 3: r=2, comparaison e
    sub_r2 = sub[sub["r"] == 2]
    e_values = sorted(sub_r2["e"].dropna().unique())
    plot_violin_box(
        sub_r2,
        x_col="e",
        y_col="jaccard_mean",
        order=e_values,
        title=f"Jaccard intra-cluster - {annot} - r=2 (compare e)",
        xlabel="Expansion e",
        output_file=f"{PLOTS_DIR}/jaccard_{annot}_r2_compare_e.png",
    )

print("All Jaccard plots saved in:", PLOTS_DIR)