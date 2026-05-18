# PMIND_Clustering_Proteines

L’objectif de ce projet est d’explorer l’intégration des embeddings protéiques dans le processus de clustering, en complément des graphes de similarité issus de BLAST. L’hypothèse sous-jacente est que la combinaison de ces deux sources d’information permettra d’obtenir des clusters plus homogènes du point de vue fonctionnel, et ainsi améliorer la qualité de l’annotation automatique des protéines.

## Vue d'ensemble du projet

**Pipeline du projet**
- **Étape 1 :** acquisition des Embeddings des protéines - utilisation de proteinBERT (dossier *Embeddings*)
- **Étape 2 :** calcul de la combinaison linéaire pondérée pour différentes $\alpha$ entre BLAST et la similarité cosinus des Embeddings (dossier *BLAST_Embeddings*)
- **Étape 3 :** analyse de BLAST et de BLAST & Embeddings (dossier *Components_analysis*)
- **Étape 4 :** étude et utilisation des annotations protéiques (dossier *Annotations*)
- **Étape 5 :** application du kNN (dossier *kNN*)
- **Étape 6 :** clustering avec l'algorithme Leiden et MCL pour la première approche (dossier *Models*)
- **Étape 7 :** application de la deuxième approche (dossier *2_app*)

## Structure du projet
``` 
├───2_app/
│   │   algo_2_app.py   # exécute la boucle principale de l'approche 2
│   │   separate_clusters_into_files.py   # sépare les clusters dans différents fichiers pour les indexer - optimisation approche 2
├───Annotations/
│   │   annotations_metrics.py   
│   │   filtrer_annotations.py   
├───BLAST/
│   │   reduce_BLAST.py          # calcule l'intersection entre BLAST et BLAST & Embeddings, afin de ne garder dans BLAST que les protéines de BLAST & Embeddings
├───BLAST_Embeddings/
│   │   apply_cos_threshold.py          # applique un seuil au graphe avec la similarité cosinus des embeddings
│   │   generate_alpha_variants.py      # génère les graphes résultants de la combinaison linéaire pondérée entre BLAST et la similarité cosinus
├───Components_analysis/
│   │   alpha_components.png            # violinplot de la distribution des composantes de BLAST & Embeddings
│   │   blast_components.png            # violinplot de la distribution des composantes de BLAST
│   │   reduced_blast_components.png    # violinplot de la distribution des composantes de BLAST réduit
│   │   components_violinplot.py        # calcule les violinplot de la distribution des composantes du graphe demandé
│   │   organize_components.py          # crée un fichier associant à chaque protéine l'ID de sa composante
├───Embeddings/
│   │   create_all_darkdino.ipynb   # crée le fichier csv contenant les protéines à partir des génomes .fasta et affiche quelques statistiques
│   │   deduplicate_dat.py          # supprime les doublons dans le graphe final des embeddings
│   │   esm2_script.py              # crée les embeddings des protéines en utilisant le modèle ESM_2
│   │   generate_cos_graph.py       # crée le graphe avec la similarité cosinus des embeddings 
|   |   proteinBERT_script.py       # crée les embeddings des protéines en utilisant le modèle ProteinBERT
│   └───Data/
│   |   |   darkdino_cleaned_for_bert.csv   # fichier csv contenant toutes les protéines sans X utilisé pour ProteinBERT
│   |   |   darkdino_cos_graph.tsv          # fichier tsv contenant le graph de la similarité cosinus des embeddings de ProteinBERT
│   |   |   embeddings_esm2.bin             # fichier binaire contenant les embeddings créé par ESM_2
│   |   |   embeddings_proteinbert.dat      # fichier dat contenant les embeddings créé par ProteinBERT
├───kNN/
|   | create_file_kNN.py        # crée un fichier csv pour le kNN
|   | execute_kNN.py            # exécute le kNN
|   | plot_metrics_results.py   # affiche les résultats du kNN en créant kNN_metrics.png 
|   | kNN_metrics.png           # résultats du kNN
├───Models/
│   └───evaluate_results/
│   |   |   analyse_metrics_results.py                  # crée les graphes des métriques
│   |   |   eval_metrics.py                             # calcule les résultats des métriques de clustering
|   |   |   alpha_leiden_all_metrics_résolution*.png    # graphiques des résultats de leiden pour tous les alphas
|   |   |   alpha_mcl_all_metrics_inflation*.png        # graphiques des résultats de mcl pour tous les alphas
│   └───rand_score/
│   |   |   clusters_comparison.py       # compare BLAST et BLAST & Embeddings avec alpha variant en calculant le score de rand
│   |   |   plot_rand_score.py           # affiche le score de rand et crée rand_score.png
|   |   |   plot_rand_score_seed.py           # affiche le score de rand pour plusieurs graines et crée rand_score_seed.png
│   |   |   rand_score_leiden_seed.png               # rand score de Leiden pour plusieurs graines entre BLAST et BLAST & Embeddings avec alpha variant, crée par plot_rand_score.py
│   |   |   rand_score_leiden.png               # rand score de Leiden entre BLAST et BLAST & Embeddings avec alpha variant, crée par plot_rand_score.py
│   |   |   rand_score_mcl.png               # rand score de mcl entre BLAST et BLAST & Embeddings avec alpha variant, crée par plot_rand_score.py
│   └───time/
│   |   |   calculate_time_leiden.py     # calcule le temps d'exécution de l'agorithme de Leiden sur les 20 première composantes
│   |   |   calculate_time_mcl.py     # calcule le temps d'exécution de l'agorithme MCL sur les 20 première composantes 
│   |   |   plot_time.py                 # affiche la courbe de temps d'exécution et crée time_plot.png
│   |   |   time_leiden_plot.png         # courbe de temps d'exécution de Leinden créée par plot_time.py
│   |   |   time_mcl_plot.png         # courbe de temps d'exécution de mcl créée par plot_time.py
│   |   leiden_clustering.py        # exécute Leiden sur le graphe demandé
│   |   leiden_clustering_seed.py        # exécute Leiden sur le graphe demandé avec plusieurs graines
│   |   mcl_clustering.py        # exécute mcl sur le graphe demandé
```

## Guide d'utilisation
Peut être exécuté si vous possédez :
- un dossier contenant les protéines d'un génôme (ID et séquence d'acides aminés)
- le graphe BLAST des dites protéines
- les annotations des protéines

### Étape 1 : acquisition des Embeddings des protéines
Utilisation de proteinBERT (dossier *Embeddings*)
1.  Lancer create_all_darkdino.ipynb
2.  Lancer proteinBERT_script.py
3.  Lancer generate_cos_graph.py

### Étape 2 : calcul de la combinaison linéaire pondérée pour différentes $\alpha$ entre BLAST et la similarité cosinus des Embeddings
(dossier *BLAST_Embeddings*)
1.  Lancer apply_cos_threshold.py
2.  Lancer generate_alpha_variants.py

### Étape 3 : analyse de BLAST et de BLAST & Embeddings
(dossiers *BLAST* et *Components_analyse*)
1.  Lancer reduce_BLAST.py
2.  Lancer organize_components.py en exécutant tour à tour avec BLAST et BLAST & Embeddings
3.  Lancer components_violinplot.py en exécutant tour à tour avec BLAST et BLAST & Embeddings

### Étape 4 : clustering avec l'algorithme Leiden
(dossier *Leiden*)
1.  Lancer leiden_clustering.py en exécutant tour à tour avec BLAST et BLAST & Embeddings
#### Obtention des résultats des métriques d'homogénéité et de sépérabilité des clusters d'une partition
1.  Lancer eval_metrics.py en exécutant tour à tour avec BLAST et BLAST & Embeddings
2.  Lancer analyse_metrics_results.py en exécutant tour à tour avec BLAST et BLAST & Embeddings
#### Comparaison entre les partition à l'aide de l'Index de Rand
1.  Lancer clusters_comparison.py
2.  Lancer plot_rand_score.py
#### Obtention des temps d'exécution de Leiden et de la comparaison avec la complexité temporelle théorique et pratique
1.  Lancer calculate_time_leiden.py
2.  Lancer plot_time.py

### Étape 5 : étude et utilisation des annotations protéiques
(dossier *Annotations*)

