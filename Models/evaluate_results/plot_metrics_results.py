import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Choix du modèle considéré
model = "leiden"
param_vals = [0.1, 0.5, 1.0, 2.0, 5.0] # gamma
param_name = "gamma"

# model = "mcl"
# param_vals = [1.3, 1.5, 2.0, 2.5, 3.0] # r
# param_name = "inflation"

# Chemins
global_path = "/tempory/21234701"

# BLAST
csv = f"{global_path}/{model}_validation_score/reduced_blast_cluster_validation_scores"

# BLAST & Embeddings
ALPHA_LIST = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
ALPHA_CSV = f"{global_path}/{model}_validation_score/combined_graph_alpha_cluster_validation_scores_{model}"
csvs = [f"{ALPHA_CSV}{a}" for a in ALPHA_LIST]

ALPHA_FIG = f"./alpha_{model}_all_metrics_{param_name}"
alpha_figs = [f"{ALPHA_FIG}{a}.png" for a in ALPHA_LIST]

def traitement_graphe(csv_file, fig_name):
    fig, axs = plt.subplots(3, 2, figsize=(13, 14))
    co_dist_list = []
    ncut_list = []
    nb_clusters_list = []
    df_nodes_list = []
    df_edges_list = []
    
    print(fig_name)

    for param in param_vals:
        df = pd.read_csv(csv_file + f"_{param}.csv")

        # On retire les lignes où le cluster a 0 sommet ou 0 arêtes
        df = df.loc[df["nb_nodes"] >= 10]
        
        nb_clusters_list.append(df["cluster_id"].nunique())
        co_dist_list.append(df["co_distance"].mean())
        ncut_list.append(df["ncut"].mean())
        df_nodes_list.append(df["nb_nodes"])
        df_edges_list.append(df["nb_edges"])

    # Courbe co-distance
    axs[0, 0].plot(param_vals, co_dist_list, marker="o", linewidth=2, markersize=8)
    axs[0, 0].set_xlabel(f"{param_name}", fontsize=12)
    axs[0, 0].set_ylabel("co-distance", fontsize=12)
    axs[0, 0].set_title(f"Co-distance moyenne en fonction de {param_name}", fontsize=14, fontweight="bold")
    axs[0, 0].grid(True, alpha=0.3)

    # Courbe ncut
    axs[0, 1].plot(param_vals, ncut_list, marker="o", color="orange", linewidth=2, markersize=8)
    axs[0, 1].set_xlabel(f"{param_name}", fontsize=12)
    axs[0, 1].set_ylabel("ncut", fontsize=12)
    axs[0, 1].set_title(f"Ncut moyen en fonction de {param_name}", fontsize=14, fontweight="bold")
    axs[0, 1].grid(True, alpha=0.3)

    # Courbe nombre de clusters
    axs[1, 0].plot(param_vals, nb_clusters_list, marker="o", color="purple", linewidth=2, markersize=8)
    axs[1, 0].set_xlabel(f"{param_name}", fontsize=12)
    axs[1, 0].set_ylabel("nombre de clusters", fontsize=12)
    axs[1, 0].set_title(f"Nombre de clusters en fonction de {param_name}", fontsize=14, fontweight="bold")
    axs[1, 0].grid(True, alpha=0.3)
    
    # Histogramme de la taille des clusters (nb nœuds) - changement format long
    nodes_data = []
    for i, param in enumerate(param_vals):
        for value in df_nodes_list[i]:
            nodes_data.append({f"{param_name}": str(param), "Nombre de nœuds": value})
    df_nodes_long = pd.DataFrame(nodes_data)
    
    sns.histplot(data=df_nodes_long, x="Nombre de nœuds", hue=f"{param_name}", element="step", stat="count", ax=axs[1, 1], palette="Set2", fill=True)
    axs[1, 1].set_xlabel(f"Nombre de nœuds", fontsize=12)
    axs[1, 1].set_ylabel("Nombre de clusters", fontsize=12)
    axs[1, 1].set_xscale('log')
    axs[1, 1].set_title(f"Nœuds par cluster en fonction de {param_name}", fontsize=14, fontweight="bold")
    axs[1, 1].grid(True, alpha=0.3, axis="y")
    
    # Violin plot nb_nodes

    sns.violinplot(data=df_nodes_long, log_scale=True, x=f"{param_name}", y="Nombre de nœuds", ax=axs[2, 0], palette="Set2", cut=2)
    axs[2, 0].set_xlabel(f"{param_name}", fontsize=12)
    axs[2, 0].set_ylabel("Nombre de nœuds", fontsize=12)
    axs[2, 0].set_yscale('log')
    axs[2, 0].set_ylim(2, 10**7)
    axs[2, 0].set_title(f"Nœuds par cluster en fonction de {param_name}", fontsize=14, fontweight="bold")
    axs[2, 0].grid(True, alpha=0.3, axis="y")

    # Violin plot nb_edges - changement format long
    edges_data = []
    for i, param in enumerate(param_vals):
        for value in df_edges_list[i]:
            edges_data.append({f"{param_name}": str(param), "Nombre d'arêtes": value})
    df_edges_long = pd.DataFrame(edges_data)

    sns.violinplot(data=df_edges_long, log_scale=True, x=f"{param_name}", y="Nombre d'arêtes", ax=axs[2, 1], palette="Set1", cut=2)
    axs[2, 1].set_xlabel(f"{param_name}", fontsize=12)
    axs[2, 1].set_ylabel("Nombre d'arêtes", fontsize=12)
    axs[2, 1].set_yscale('log')
    axs[2, 1].set_ylim(2, 10**7)
    axs[2, 1].set_title(f"Arêtes par cluster en fonction de {param_name}", fontsize=14, fontweight="bold")
    axs[2, 1].grid(True, alpha=0.3, axis="y")

    # # On masque le dernier subplotplt.show()
    # axs[2, 1].axis("off")
    
    plt.tight_layout()
    plt.savefig(fig_name, bbox_inches="tight")
    
# BLAST
traitement_graphe(csv, f"./Images/blast_{model}_all_metrics_{param_name}.png")

# # BLAST & embeddings
# for csv, fig in zip(csvs, alpha_figs):
#     traitement_graphe(csv, fig)
    
print("terminé !")