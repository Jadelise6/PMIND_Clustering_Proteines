import numpy as np
import matplotlib.pyplot as plt

# Définir les bins
num_bins = 100
bins = np.linspace(0, 1, num_bins + 1)
hist = np.zeros(num_bins, dtype=int)

with open("graph_edges.tsv") as f:
    for i, line in enumerate(f):
        parts = line.strip().split()
        if len(parts) == 3:
            w = float(parts[2])
            # Trouver l'indice du bin correspondant
            idx = min(int(w * num_bins), num_bins - 1)
            hist[idx] += 1

        if i % 1_000_000 == 0:
            print(f"Lignes lues : {i}")

# Affichage / sauvegarde
plt.bar(bins[:-1], hist, width=0.01)
plt.xlabel("Poids")
plt.ylabel("Nombre d’arêtes")
plt.title("Distribution des poids")
plt.savefig("histogramme_poids.png", dpi=300)
plt.close()

print("Terminé ! Histogramme enregistré dans histogramme_poids.png")
