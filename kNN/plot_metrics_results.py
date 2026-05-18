import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("ticks")
sns.set_context("talk")

# Chemins
global_path = "/tempory/21234701"

# BLAST
input_file = f"{global_path}/kNN/metrics_kNN_combined_graph_alpha_1.0.csv"

fig = f"./kNN_metrics.png"

def traitement_graphe(csv_file, fig_name):
    plt.figure(figsize=(9, 6))
    df = pd.read_csv(csv_file)
    
    plt.plot(df["k"], df["Recall"], marker="o", linewidth=2, markersize=8, label="Rappel")
    plt.plot(df["k"], df["Precision"], marker="s", linewidth=2, color="orange", markersize=8, label="Précision")
    plt.plot(df["k"], df["F1"], marker="v", linewidth=2, color="purple", markersize=8, label="F1-score")
             
    plt.xlabel(f"k", fontsize=17)
    plt.ylabel("Score", fontsize=17)
    plt.title(f"Scores en fonction de k", fontsize=14, fontweight="bold")
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(fig_name, bbox_inches="tight")
    
# BLAST & embeddings
traitement_graphe(input_file, fig)
    
print("terminé !")