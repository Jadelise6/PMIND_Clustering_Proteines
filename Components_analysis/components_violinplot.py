import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Chemins
global_path = "/tempory/21234701"

# BLAST
COMP_FILE = f"{global_path}/blast_components.tsv"   # Fichier contenant les composantes

# BLAST & Embeddings
# Ils ont les mêmes composantes donc on peut tester sur un seul fichier
ALPHA_COMP_FILE = f"{global_path}/components_combined_graph_alpha_0.2.tsv"

comp_nb = {}
with open(COMP_FILE, "r") as f_in:    
    for line in f_in:
        parts = line.split()
        if len(parts) < 4:
            continue
        com_id = parts[3]
        comp_nb[com_id] = comp_nb.get(com_id, 0) + 1

# Séparatin des nombreux nœuds isolés du reste
isolated = sum(1 for v in comp_nb.values() if v == 1)
non_isolated = [v for v in comp_nb.values() if v > 1]

# Statistiques - affichage latex
q75 = sorted(non_isolated)[int(len(non_isolated)*0.75)]
zoomed = [x for x in non_isolated if x <= q75*2]
percentage = (len(zoomed) / len(non_isolated)) * 100

stats = {
    'Statistique': ['Min', 'Q1', 'Médiane', 'Q3', 'Max', 'Moyenne', 'Écart-type'],
    'Valeur': [
        min(non_isolated),
        sorted(non_isolated)[len(non_isolated)//4],
        sorted(non_isolated)[len(non_isolated)//2],
        sorted(non_isolated)[3*len(non_isolated)//4],
        max(non_isolated),
        f"{np.mean(non_isolated):.2f}",
        f"{np.std(non_isolated):.2f}"
    ]
}

df = pd.DataFrame(stats)
print(f"Nœuds isolés: {isolated}")
print(f"Composantes (>1 nœud): {len(non_isolated)}")
print(f"Distribution zoomée: {len(zoomed)} composantes ({percentage:.1f}%)")
print(df.to_latex(index=False))
print("\n")

# Plot juste le zoomed
sns.set_style("whitegrid")
plt.figure(figsize=(8, 4))
sns.violinplot(x=zoomed, color="steelblue", inner="quartile")
plt.margins(y=0.1)
plt.title(f"Distribution des composantes (zoomée) - {percentage:.1f}% des composantes", fontsize=14, fontweight="bold")
plt.ylabel("Nombre de nœuds par composante", fontsize=12)
plt.grid(True, alpha=0.3, axis="y")
plt.tight_layout()
# Changer titre si blast / alpha
plt.savefig("./blast_components.png", bbox_inches="tight", dpi=300)
plt.show()