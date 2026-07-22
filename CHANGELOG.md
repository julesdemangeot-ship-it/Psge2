# Historique PSGE II

## Version 2.6.1 — Universe UX Stabilization

- **`psge_engine.py`** : module canonique publié, importable (`import psge_engine as E`)
- **`Universe`** : nouvelle classe de haut niveau pour manipuler un univers discret
  - `Universe(tetrahedra, vertices)` — constructeur simple
  - `run()` / `report` — calcul paresseux et mis en cache
  - `is_closed()` — détecte si la variété est fermée (toutes les charnières intérieures)
  - `total_volume()`, `regge_action()`, `bulk_action()`, `boundary_action()`
  - `topology_distribution()` — répartition des charnières par type topologique
  - `curvature_at(edge)` — angle de déficit en O(1) sur une arête
  - `summary()` — rapport lisible en une ligne
- **`HingeTopology`** : classification topologique des charnières (INTERIOR / BOUNDARY / NON_MANIFOLD / DEGENERATE)
- **`HingeTopologyAnalyzer`** : analyse du link d'une arête (graphe de lien S¹ vs chemin ouvert)
- **`BarycentricDualCell`** : stratégie de mesure duale barycentrique (remplacement de la longueur brute)
- **`ReggeGeometryReport`** : rapport unifié avec `bulk_regge_sum`, `boundary_hinge_sum`, `total_geometric_sum` ; `total_action` conservé comme alias de compatibilité
- **`MetricTensor`** : enrichi de `norm`, `norm2`, `project_orthogonal`, `gram_matrix`, `gram_area`
- **`DeficitComputer`** : calcul du dièdre via projection métrique correcte (paramètre `tensor`)
- Rétrocompatibilité totale avec les scripts existants (`Démonstration psge2`, `PSGE2 test charniere interne`)
- 6 tests de validation mathématique et topologique passent (suite `PSGETestSuite`)

## Version 1.0

- Noyau géométrique validé
- Calcul des volumes
- Matrice de Gram
- Angles dièdres
- Déficit de Regge
- 19 tests validés
