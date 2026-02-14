import pandas as pd
import matplotlib.pyplot as plt

clusters_file = "graph_mcl.out"
sizes_file = "cluster_sizes.tsv"
hist_file = "cluster_sizes_hist.png"

# --- Partie 1 : calcul des tailles de clusters ---
cluster_sizes = []
with open(clusters_file) as f:
    for i, line in enumerate(f):
        proteins = line.strip().split()
        cluster_sizes.append((i, len(proteins)))

df = pd.DataFrame(cluster_sizes, columns=["cluster_id", "size"])
df.to_csv(sizes_file, sep="\t", index=False)
print(f"Distribution des tailles de clusters enregistrée dans {sizes_file}")

# --- Partie 2 : histogramme ---
df = pd.read_csv(sizes_file, sep="\t")
plt.hist(df["size"], bins=50)
plt.xlabel("Taille du cluster")
plt.ylabel("Nombre de clusters")
plt.title("Distribution des tailles de clusters MCL")
plt.savefig(hist_file)  # Enregistrer directement
plt.close()
print(f"Histogramme enregistré dans {hist_file}")

plt.hist(df["size"], bins=50, log=True)
plt.xlabel("Taille du cluster")
plt.ylabel("Nombre de clusters (log)")
plt.title("Distribution des tailles de clusters MCL")
plt.savefig("cluster_sizes_hist_log.png")
plt.close()

plt.hist(df["size"], bins=50)
plt.xlabel("Taille du cluster")
plt.ylabel("Nombre de clusters")
plt.title("Distribution des tailles de clusters MCL")
plt.xlim(0, 200) # Limite l'axe x
plt.savefig("cluster_sizes_hist_zoom.png")
plt.close()