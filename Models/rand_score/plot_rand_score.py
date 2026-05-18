import matplotlib.pyplot as plt
import pandas as pd

# Choix du modèle considéré
#model = "leiden"
#param_vals = [0.1, 0.5, 1.0, 2.0, 5.0] # gamma
#param_name = "gamma"
param_symbol = "\gamma"

model = "mcl"
param_vals = [1.3, 1.5, 2.0, 2.5, 3.0] # r
param_name = "inflation"
param_symbol = "r"


# Lire le fichier CSV créé avec la bonne structure
df = pd.read_csv(f"rand_score_{model}.csv")

# Récupérer les alphas
alpha_values = df["alpha"].tolist()

# Pour chaque gamma
param_cols = [col for col in df.columns if col.startswith(f"{param_name}_")]

# Plot pour chaque gamma
plt.figure(figsize=(12, 6))

for param_col in param_cols:
    param_label = param_col.replace(f"{param_name}_", f"${param_symbol}=$")
    param_scores = df[param_col].tolist()
    plt.plot(alpha_values, param_scores, marker="o", linewidth=2, markersize=8, label=param_label)

plt.xlabel("Valeurs de alpha", fontsize=12)
plt.ylabel("Rand Score (moyenne)", fontsize=12)
plt.title(f"Rand score en fonction des valeurs de alpha et {param_name} - {model}", fontsize=14, fontweight="bold")
plt.grid(True, alpha=0.3)
plt.legend()
plt.savefig(f"./rand_score_{model}.png", bbox_inches="tight")
plt.show()