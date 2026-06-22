# PMIND_Clustering_Proteines

L’objectif de ce projet est d’explorer l’intégration des embeddings protéiques dans le processus de clustering, en complément des graphes de similarité issus de BLAST.

L’hypothèse sous-jacente est que la combinaison de ces deux sources d’information permettra d’obtenir des clusters plus homogènes du point de vue fonctionnel, et ainsi améliorer la qualité de l’annotation automatique des protéines.

<details>
    <summary>Table des Matières</summary>
    <ol>
        <li><a href="#description-du-projet">Description du projet</a></li>
        <li><a href="#données-utilisées">Données utilisées</a></li>
        <li><a href="#ℹvue-densemble-du-projet">Vue d'ensemble du projet</a></li>
        <li><a href="#structure-du-projet">Structure du projet</a></li>
        <li>
            <a href="#guide-dutilisation">Guide d'utilisation</a>
            <ol>
                <li><a href="#dépendances">Dépendances</a></li>
                <li><a href="#étapes-dexécution">Étapes d'exécution</a></li>
            </ol>
        </li>
        <li><a href="#auteurs">Auteurs</a></li>
    </ol>
</details>

## Description du projet
Lors de ce projet, nous étudions des méthodes de clustering pour résoudre le vaste problème de l'annotation des protéines. Notre objectif est d'améliorer l'homogénéité des clusters obtenus en utilisant ces méthodes, et d'innover en intégrant les embeddings de protéines. L'homogénéité des clusters est évaluée par des mesures classiques, mais également par des mesures fonctionnelles : **l'incompatibilité entre deux annotatons**. En effet, certaines annotations sont incompatibles, et de ce fait, ne doivent pas se retrouver dans le même cluster. Les méthodes que nous employons visent à optimiser cette contrainte en restant très prudent, au risque d'avoir un rappel faible.

### Le kNN comme standard d'annotation fonctionnelle
**Inconvénients du kNN**
- Complexité temporelle : O(M×V×log(V))
avec : 
  - M : nombre de protéines à annoter
  - V : degré de la protéine à annoter
  - log(V) : coût du tri des voisins
- Trop coûteux à grande échelle
- Sensible aux valeurs aberrantes
- Ne peut voir que le connu

### Algorithes de clustering considérés
#### MCL
(Markov Clustering aLgorithm)

**Principe :** marche aléatoire sur le graphe

**Complexité :** O($n^3$)

**Paramètres :**
- Puissance e : contrôle l’expansion de la matrice
- Inflation r : contrôle la granularité des clusters

#### Leiden

**Principe :** optimisation de la modularité

**Complexité :** O(n log(n))

**Paramètres :**
- Résolution gamma : contrôle la granularité des clusters

### Métriques d'évaluation
**Co-distance**
- mesure d'homogénéité
- points proches entre eux au sein du cluster

**Indice de Jaccard**
- mesure de cohérence biologique
- détecte l'incompatibilité entre deux annotations

**Rappel, précision et f1-score**

### Première approche
1. **Combinaison linéaire pondérée**
   - graphe de similarité **BLAST**
   - graphe de la similarité cosinus des **embeddings**
2. **Clustering avec MCL et Leiden**
Recherche des meilleurs paramètres des algorithmes de clustering et du meilleur alpha
3. **Calcul et analyse des métriques de clustering** pour juger de la qualité de la partition 
obtenue (la co-distance pour l'homogénéité et la coupe normalisée pour la séparabilité)
4. **Propagation des annotations** au sein de la partition
5. **Calcul et analyse des mesures de la qualité de la prédiction** des annotations (f1-score, 
rappel, précision)

### Deuxième approche
1. **Clustering avec Leiden et gamma = 1,0** (meilleurs modèle et paramètre trouvés lors de l'approche 1)
   - graphe de similarité BLAST
   - graphe de la similarité cosinus des embeddings
2. **Combinaison des deux partitions obtenues
 Boucle principale**
   - Propagation des annotations dans la partition de BLAST
   - Pour toute protéine, reporter son annotation dans BLAST aux embeddings
   - Propagation des annotations dans la partition des embeddings
   - Pour toute protéine, reporter son annotation dans les embeddings à BLAST
3. **Condition d’arrêt :** pas de nouvelle annotation ajoutée à BLAST ou aux 
embeddings à la fin d’une itération

## Données utilisées
- **Séquences protéiques**
  - environ 7,4 M
  - format FASTA
  - dinoflagellés

- **Annotations fonctionnelles PFAM**
  - 30,9 % de protéines annotées
  - 2,9 M dont 9146 distinctes
  - 1 à 35 par protéine (multi-label)

## Vue d'ensemble du projet

**Pipeline du projet**
- **Étape 1 :** acquisition des Embeddings des protéines - utilisation de proteinBERT (dossier *Embeddings*)
- **Étape 2 :** calcul de la combinaison linéaire pondérée pour différents $\alpha$ entre le graphe de similarité BLAST et le graphe issu de la similarité cosinus des Embeddings (dossier *BLAST_Embeddings*)
- **Étape 3 :** analyse de BLAST et de BLAST & Embeddings (dossier *Components_analysis*)
- **Étape 4 :** étude et utilisation des annotations protéiques (dossier *Annotations*)
- **Étape 5 :** application du kNN (dossier *kNN*)
- **Étape 6 :** clustering avec l'algorithme Leiden et MCL pour la première approche (dossier *Models*)
- **Étape 7 :** propagation des annotations pour la première et la deuxième approche (dossier *Annotation_propagation*)

## Structure du projet
``` 
├───Annotation_propagation/
│   │   algo_2_main_loop.py                # exécute la boucle principale de l'approche 2
│   │   separate_clusters_into_files.py    # sépare les clusters dans différents fichiers pour les indexer - optimisation approche 2
├───Annotations/
│   │   annotations_metrics.py   
│   │   filtrer_annotations.py   
├───BLAST/
│   │   reduce_BLAST.py                    # calcule l'intersection entre BLAST et BLAST & Embeddings, afin de ne garder dans BLAST que les protéines de BLAST & Embeddings
├───BLAST_Embeddings/
│   │   apply_cos_threshold.py             # applique un seuil au graphe avec la similarité cosinus des embeddings
│   │   generate_alpha_variants.py         # génère les graphes résultants de la combinaison linéaire pondérée entre BLAST et la similarité cosinus
├───Components_analysis/
│   │   alpha_components.png               # violinplot de la distribution des composantes de BLAST & Embeddings
│   │   blast_components.png               # violinplot de la distribution des composantes de BLAST
│   │   reduced_blast_components.png       # violinplot de la distribution des composantes de BLAST réduit
│   │   components_violinplot.py           # calcule les violinplot de la distribution des composantes du graphe demandé
│   │   organize_components.py             # crée un fichier associant à chaque protéine l'ID de sa composante
├───Embeddings/
│   │   create_all_darkdino.ipynb          # crée le fichier csv contenant les protéines à partir des génomes .fasta et affiche quelques statistiques
│   │   deduplicate_dat.py                 # supprime les doublons dans le graphe final des embeddings
│   │   esm2_script.py                     # crée les embeddings des protéines en utilisant le modèle ESM_2
│   │   generate_cos_graph.py              # crée le graphe avec la similarité cosinus des embeddings 
|   |   proteinBERT_script.py              # crée les embeddings des protéines en utilisant le modèle ProteinBERT
│   └───Data/
│   |   |   darkdino_cleaned_for_bert.csv  # fichier csv contenant toutes les protéines sans X utilisé pour ProteinBERT
│   |   |   darkdino_cos_graph.tsv         # fichier tsv contenant le graph de la similarité cosinus des embeddings de ProteinBERT
│   |   |   embeddings_esm2.bin            # fichier binaire contenant les embeddings créé par ESM_2
│   |   |   embeddings_proteinbert.dat     # fichier dat contenant les embeddings créé par ProteinBERT
├───kNN/
|   | create_file_kNN.py                   # crée un fichier csv pour le kNN
|   | execute_kNN.py                       # exécute le kNN
|   | plot_metrics_results.py              # affiche les résultats du kNN en créant kNN_metrics.png 
|   | kNN_metrics.png                      # résultats du kNN
├───Models/
│   └───evaluate_results/
│   |   |   analyse_metrics_results.py                  # crée les graphes des métriques
│   |   |   eval_metrics.py                             # calcule les résultats des métriques de clustering
|   |   |   alpha_leiden_all_metrics_résolution*.png    # graphiques des résultats de leiden pour tous les alphas
|   |   |   alpha_mcl_all_metrics_inflation*.png        # graphiques des résultats de mcl pour tous les alphas
│   └───rand_score/
│   |   |   rand_score.py                   # compare BLAST et BLAST & Embeddings avec alpha variant en calculant le score de rand
│   |   |   plot_rand_score.py              # affiche le score de rand et crée rand_score.png
|   |   |   plot_rand_score_seed.py         # affiche le score de rand pour plusieurs graines et crée rand_score_seed.png
│   |   |   rand_score_leiden_seed.png      # rand score de Leiden pour plusieurs graines entre BLAST et BLAST & Embeddings avec alpha variant, crée par plot_rand_score.py
│   |   |   rand_score_leiden.png           # rand score de Leiden entre BLAST et BLAST & Embeddings avec alpha variant, crée par plot_rand_score.py
│   |   |   rand_score_mcl.png              # rand score de mcl entre BLAST et BLAST & Embeddings avec alpha variant, crée par plot_rand_score.py
│   └───time/
│   |   |   calculate_time_leiden.py        # calcule le temps d'exécution de l'agorithme de Leiden sur les 20 première composantes
│   |   |   calculate_time_mcl.py           # calcule le temps d'exécution de l'agorithme MCL sur les 20 première composantes 
│   |   |   plot_time.py                    # affiche la courbe de temps d'exécution et crée time_plot.png
│   |   |   time_leiden_plot.png            # courbe de temps d'exécution de Leinden créée par plot_time.py
│   |   |   time_mcl_plot.png               # courbe de temps d'exécution de mcl créée par plot_time.py
│   |   leiden_clustering.py                # exécute Leiden sur le graphe demandé
│   |   leiden_clustering_seed.py           # exécute Leiden sur le graphe demandé avec plusieurs graines
│   |   mcl_clustering.py                   # exécute mcl sur le graphe demandé
```
*(Certains scripts de Mathilde Gauteur manquent à l'appel, notamment sur les annotations et la fin des approches 1 et 2).*

## Guide d'utilisation

Ce projet peut être exécuté si vous disposez des éléments suivants :
- un ensemble de protéines au format FASTA, avec un identifiant et une séquence d'acides aminés pour chaque entrée ;
- le graphe BLAST correspondant à ces protéines ;
- les annotations fonctionnelles associées aux protéines.

### Dépendances

Outils externes :
- MCL
- proteinBERT

Bibliothèques Python pouvant être installées avec `pip` :
- leidenalg
- igraph
- matplotlib
- numpy
- pandas
- scikit-learn
- seaborn

### Étapes d'exécution

#### 1. Génération des embeddings protéiques
Répertoire : `Embeddings`
1. Lancer `create_all_darkdino.ipynb` pour préparer les séquences.
2. Lancer `proteinBERT_script.py` pour calculer les embeddings.
3. Lancer `generate_cos_graph.py` pour construire le graphe de similarité cosinus.

#### 2. Construction du graphe combiné BLAST + embeddings
Répertoire : `BLAST_Embeddings`
1. Lancer `apply_cos_threshold.py` pour filtrer les arêtes selon le seuil de similarité cosinus.
2. Lancer `generate_alpha_variants.py` pour générer les graphes pondérés pour différentes valeurs de $\alpha$.

#### 3. Analyse des composantes connexes
Répertoires : `BLAST` et `Components_analysis`
1. Lancer `reduce_BLAST.py` pour restreindre le graphe BLAST aux protéines présentes dans le graphe combiné.
2. Lancer `organize_components.py` pour associer à chaque protéine son composant connexe.
3. Lancer `components_violinplot.py` pour visualiser la distribution des tailles de composantes.

#### 4. Exploitation des annotations
Répertoire : `Annotations`
1. Lancer les scripts de filtrage et d'analyse des annotations selon le besoin du pipeline.
2. Utiliser ces résultats pour évaluer l'homogénéité fonctionnelle des clusters.

#### 5. Construction et évaluation du kNN
Répertoire : `kNN`
1. Lancer `create_file_kNN.py` pour préparer le fichier d'entrée.
2. Lancer `execute_kNN.py` pour exécuter l'algorithme.
3. Lancer `plot_metrics_results.py` pour visualiser les scores obtenus.

#### 6. Première approche de clustering
Répertoire : `Models`
1. Lancer `leiden_clustering.py` pour exécuter Leiden.
2. Lancer `mcl_clustering.py` pour exécuter MCL.

##### Évaluation des partitions
Répertoire : `Models/evaluate_results`
1. Lancer `eval_metrics.py` sur les partitions produites par Leiden et MCL.
2. Lancer `analyse_metrics_results.py` pour comparer les métriques d'homogénéité et de séparabilité.

##### Propagation des annotations
Répertoire : `Models/propagation`
1. Lancer `filtrage_clusters.py` pour sélectionner les clusters pertinents.
2. Lancer `propagation_app_1.py` pour propager les annotations dans la première approche.

##### Comparaison des partitions avec l'indice de Rand
Répertoire : `Models/rand_score`
1. Lancer `rand_score.py` pour calculer l'indice de Rand.
2. Lancer `plot_rand_score.py` pour produire les graphiques de comparaison.
3. Lancer `plot_rand_score_seed.py` si vous souhaitez étudier la stabilité selon plusieurs graines.

##### Étude du temps d'exécution
Répertoire : `Models/time`
1. Lancer `calculate_time_leiden.py`.
2. Lancer `calculate_time_mcl.py`.
3. Lancer `plot_time.py` pour comparer les temps d'exécution observés et la complexité théorique attendue.

#### 7. Propagation des annotations dans la seconde approche
Répertoire : `Annotation_propagation`
1. Lancer `separate_clusters_into_files.py` pour répartir les clusters dans des fichiers distincts.
2. Lancer `algo_2_main_loop.py` pour exécuter la boucle principale de l'approche 

## Auteurs
* [Élise THOMAS](https://github.com/Jadelise6)
* [Mathilde GAUTEUR](https://github.com/gauteurmathilde)