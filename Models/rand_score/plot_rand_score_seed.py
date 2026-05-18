import matplotlib.pyplot as plt
import pandas as pd

# Lire le fichier CSV créé avec la bonne structure
df = pd.read_csv("rand_score_seed.csv")

# Récupérer les seed
seed_values = df["seed"].tolist()

# Pour chaque gamma
gamma_cols = [col for col in df.columns if col.startswith("gamma_")]

# Plot pour chaque gamma
plt.figure(figsize=(12, 6))

for gamma_col in gamma_cols:
    gamma_label = gamma_col.replace("gamma_", "$\gamma=$")
    gamma_scores = df[gamma_col].tolist()
    plt.plot(seed_values, gamma_scores, marker="o", linewidth=2, markersize=8, label=gamma_label)

plt.xlabel("Valeurs de seed", fontsize=12)
plt.ylabel("Rand Score (moyenne)", fontsize=12)
plt.title("Rand score pour BLAST en fonction des valeurs de seed et gamma - Leiden", fontsize=14, fontweight="bold")
plt.grid(True, alpha=0.3)
plt.legend()
plt.savefig("./rand_score_seed.png", bbox_inches="tight")
plt.show()