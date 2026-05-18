import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

df = pd.read_csv('time_leiden.csv')

# Régression linéaire pour approximation de la complexité
X = df["component_size"].to_numpy()
y = df["leiden_time"].to_numpy()
X_log = np.log10(X)
y_log = np.log10(y)
reg = LinearRegression().fit(X_log.reshape(-1, 1), y_log.reshape(-1, 1))
slope = reg.coef_[0][0]
intercept = reg.intercept_[0]
print(f"Pente de la régression linéaire (log-log): {slope:.2f}")

# Prédiction sur les données
y_reg_log = reg.predict(X_log.reshape(-1, 1))
y_reg = 10**y_reg_log.flatten()

# Plot du temps de Leiden en fonction de la taille des composantes
plt.figure(figsize=(10, 6))

# Données
plt.scatter(X, y, label="Données")
# Régression log-log
plt.plot(X, y_reg, label=f"Régression ($O(n^{{{slope:.2f}}})$)", color='red')
# Complexité linéaire O(n) - ancrage sur le plus petit point
plt.plot(X, X * (y[-1] / X[-1]), label="O(n)", color='green', linestyle='--')
# Complexité O(n log(n)) - ancrage sur le plus petit point
plt.plot(X, X * X_log * (y[-1] / X[-1]), label="$O(n \log n)$", color='purple', linestyle='--')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Taille de la composante (nombre de sommets)', fontsize=12)
plt.ylabel('Temps d\'exécution Leiden (secondes)', fontsize=12)
plt.title('Temps d\'exécution de Leiden en fonction de la taille de la composante', fontsize=14, fontweight="bold")
plt.grid(True, alpha=0.3, which='both')
plt.legend(fontsize=11)
plt.savefig('./time_leiden_plot.png', bbox_inches='tight')
plt.show()