# Historique PSGE II

## Version 2.6.1 — Universe UX stabilization

- Ajout de la classe `Universe` : API haut-niveau fluide pour construire et analyser un univers simplicial sans manipulation manuelle des certificats ni des fournisseurs de coordonnées
- Ajout des classes internes `_UniverseCertificate`, `_TopologicalCertificate`, `_UniverseCoords` supportant l'API `Universe`
- Correction de `GeometryEngine` : la dimension du tenseur métrique par défaut est désormais inférée des coordonnées des sommets (suppression du `dim = 4` codé en dur)
- Support du chaînage fluide sur `Universe.add_vertex()` et `Universe.add_tetrahedron()`
- Validation anticipée des références de sommets dans `Universe.compute()` avec message d'erreur explicite

## Version 2.6.0

- Moteur certifié PSGE II v2.6.0-L4
- 6 lemmes de validation (L4.1–L4.6) : platitude, cône intrinsèque, action de Regge, dégénérescence, invariances isométrique et de scaling
- Statut : PASS [SUCCÈS] — tolérance 1.0e-7

## Version 1.0

- Noyau géométrique validé
- Calcul des volumes
- Matrice de Gram
- Angles dièdres
- Déficit de Regge
- 19 tests validés
