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

ANNOT_FILE = "outputs/annotations_filtered_full.tsv"
CLUSTERS_DIR = "outputs/clusters_output"
METRICS_FILE = "outputs/metrics_output/jaccard_metrics.csv"
PAIR_FILE = "outputs/metrics_output/jaccard_pair_distribution.csv"
PLOTS_DIR = "outputs/metrics_output/plots_jaccard_new"

os.makedirs("outputs/metrics_output", exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

VALID_SOURCES = {
    "Pfam", "PANTHER", "Gene3D", "SMART",
    "SUPERFAMILY", "CDD", "PRINTS", "ProSiteProfiles", "ProSitePatterns"}

MIN_CLUSTER_SIZE = 5

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
        next(f)
        for line in f:
            prot, cid = line.strip().split("\t")
            clusters[cid].append(prot)
    return clusters

# -----------------------
# PLOT FUNCTIONS
# -----------------------

def add_median_labels(ax, df, x_col, y_col, order):
    medians = df.groupby(x_col)[y_col].median()
    for i, val in enumerate(order):
        if val not in medians.index:
            continue
        med = medians.loc[val]
        ax.text(i, med, f"{med:.2f}", ha="center", va="bottom",
                fontsize=8, color="black", fontweight="bold")


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

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Jaccard mean")
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
pair_stats = []

cluster_files = sorted(glob.glob(CLUSTERS_DIR + "/*clusters.tsv"))

for cfile in cluster_files:
    print("Processing", cfile)
    clusters = read_clusters(cfile)
    method = os.path.basename(cfile)

    for cid, prots in clusters.items():
        if len(prots) < 2:
            continue

        for source in VALID_SOURCES:

            values = []

            bins = {
                "zero": 0,
                "low": 0,
                "mid": 0,
                "high": 0
            }

            nb_total = 0

            for p1, p2 in itertools.combinations(prots, 2):

                a = annotations[p1][source] if p1 in annotations else set()
                b = annotations[p2][source] if p2 in annotations else set()

                if not a or not b:
                    continue

                j = jaccard(a, b)

                if j is None:
                    continue

                values.append(j)
                nb_total += 1

                if j == 0:
                    bins["zero"] += 1
                elif j <= 0.5:
                    bins["low"] += 1
                elif j < 1:
                    bins["mid"] += 1
                else:
                    bins["high"] += 1

            if not values:
                continue

            # métrique classique
            results.append({
                "cluster_file": method,
                "cluster_id": cid,
                "annotation_type": source,
                "cluster_size": len(prots),
                "jaccard_mean": np.mean(values)
            })

            # stats détaillées
            pair_stats.append({
                "cluster_file": method,
                "cluster_id": cid,
                "annotation_type": source,
                "cluster_size": len(prots),
                "nb_pairs": nb_total,
                "zero": bins["zero"],
                "low": bins["low"],
                "mid": bins["mid"],
                "high": bins["high"],
                "ratio_zero": bins["zero"] / nb_total if nb_total > 0 else 0
            })

# -----------------------
# SAVE
# -----------------------

df = pd.DataFrame(results)
pair_df = pd.DataFrame(pair_stats)

if df.empty:
    raise RuntimeError("No Jaccard values computed.")

df["r"] = df["cluster_file"].str.extract(r"r([0-9.]+)").astype(float)
df["e"] = df["cluster_file"].str.extract(r"e([0-9.]+)").astype(float)

df.to_csv(METRICS_FILE, index=False)
pair_df.to_csv(PAIR_FILE, index=False)

print("Metrics saved.")
print("Pair stats saved.")

# -----------------------
# PLOTS
# -----------------------

for annot in df["annotation_type"].unique():
    sub = df[df["annotation_type"] == annot]

    sub_e2 = sub[sub["e"] == 2]
    r_values = sorted(sub_e2["r"].dropna().unique())

    plot_violin_box(
        sub_e2,
        "r",
        "jaccard_mean",
        r_values,
        f"Jaccard - {annot} - e=2",
        "Inflation r",
        f"{PLOTS_DIR}/jaccard_{annot}_e2.png"
    )

# -----------------------
# INCOMPATIBILITY PLOT
# -----------------------

if not pair_df.empty:
    plt.figure(figsize=(10, 6))
    sns.histplot(pair_df["ratio_zero"], bins=50)
    plt.title("Distribution des paires incompatibles (Jaccard = 0)")
    plt.xlabel("Ratio de paires incompatibles")
    plt.ylabel("Nombre de clusters")
    plt.tight_layout()
    plt.savefig(PLOTS_DIR + "/incompatibility_distribution.png", dpi=300)
    plt.close()

    print("Incompatibility plot saved.")