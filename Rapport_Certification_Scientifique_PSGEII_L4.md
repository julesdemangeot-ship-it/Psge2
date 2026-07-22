# Rapport de certification scientifique  
## Moteur PSGE II — Calcul de Regge discret

**Date d’exécution :** 2026-07-22 13:58:24  
**Version moteur :** 2.6.0-L4  
**Tolérance globale de validation :** 1.0e-7  
**Statut final :** ✅ **PASS [SUCCÈS]**

---

## 1) Objet de la certification

Ce document atteste la validation numérique du moteur géométrique **PSGE II** vis-à-vis d’oracles analytiques du calcul de **Regge discret**.  
Le périmètre couvre six lemmes de vérification (L4.1 à L4.6) portant sur :

1. la platitude intrinsèque,
2. la courbure conique (via angle dièdre),
3. l’action de Regge en cas plat,
4. la robustesse sur cas dégénérés,
5. les invariances isométriques,
6. l’invariance de scaling (changement d’échelle).

---

## 2) Référentiel de test et conditions

### 2.1 Maillage de référence
- Maillage : **« éventail »**
- Arêtes-charnières extraites : **19**
- Charnières intérieures retenues pour les tests de platitude et d’action : **1**

### 2.2 Critère d’acceptation
Un test est déclaré **SUCCÈS** si l’écart absolu à l’oracle analytique est inférieur à la tolérance du lemme (ici compatible avec 1.0e-7).

### 2.3 Environnement logique validé
- Géométrie simpliciale 3D (tetraèdres)
- Déficit angulaire de Regge
- Action discrète de Regge
- Détection de dégénérescence par signature/Gram

---

## 3) Résultats consolidés

| Test / Oracle analytique | Calculé | Attendu | Écart absolu | Statut |
|---|---:|---:|---:|---|
| **L4.1 - Candidat Plat (Éventail)** | 0.00e+00 | 0.0000 | 0.00e+00 | SUCCÈS |
| **L4.2 - Cône intrinsèque (Dièdre)** | 1.23095941734 | 1.23095941734 | 0.00e+00 | SUCCÈS |
| **L4.3 - Action de Regge (Plat)** | 0.00e+00 | 0.0000 | 0.00e+00 | SUCCÈS |
| **L4.4 - Cas dégénérés (Volume nul)** | Dégénéré | Dégénéré | 0.00 | SUCCÈS |
| **L4.5 - Invariance isométrique** | 8.88e-16 | 0.0 | 8.88e-16 | SUCCÈS |
| **L4.6 - Invariance de scaling** | 0.9428 | 0.9428 | 1.11e-16 | SUCCÈS |

---

## 4) Analyse technique par lemme

### L4.1 & L4.3 — Platitude / Action nulle
Le déficit angulaire des charnières intérieures et l’action de Regge associée s’annulent sur le candidat plat.  
**Interprétation :** cohérence avec le régime de courbure nulle attendu.

### L4.2 — Validation dièdre analytique
L’angle dièdre du tétraèdre régulier reproduit **arccos(1/3)** à la précision machine.  
**Interprétation :** exactitude du calcul de géométrie locale (normales, projection orthogonale, angle sur charnière).

### L4.4 — Dégénérescence contrôlée
Le tétraèdre colinéaire est classé dégénéré (déterminant de Gram sous seuil, volume nul).  
**Interprétation :** comportement robuste sans exception non contrôlée ni résultat physique incohérent.

### L4.5 — Invariance isométrique
Le déficit reste invariant sous rotation/translation ; résidu numérique observé \(8.88\times 10^{-16}\), compatible arrondi machine.  
**Interprétation :** invariance intrinsèque respectée.

### L4.6 — Invariance d’échelle
Les grandeurs volumétriques suivent la loi de scaling attendue (\(\lambda^3\)), avec écart \(1.11\times10^{-16}\).  
**Interprétation :** cohérence dimensionnelle et stabilité numérique.

---

## 5) Conclusion de certification

Au regard des résultats obtenus, **tous les critères de validation L4 sont satisfaits**.  
Le moteur **PSGE II v2.6.0-L4** est **certifié conforme** sur ce périmètre de test :

> **CERTIFICATION SCIENTIFIQUE : PASS [SUCCÈS]**

---

## 6) Traçabilité et archivage

- **Identifiant campagne :** PSGEII-L4-2026-07-22
- **Horodatage de référence :** 2026-07-22 13:58:24
- **Jeu d’essai :** Maillage “éventail” + cas analytique dièdre + cas dégénéré colinéaire
- **Tolérance de recette :** 1.0e-7
- **Synthèse statut :** 6/6 SUCCÈS, 0 échec, 0 anomalie bloquante

---

## 7) Recommandations post-certification (qualité continue)

1. **Gel de baseline** : conserver ce rapport comme référence de non-régression L4.  
2. **CI scientifique** : exécuter L4.1–L4.6 à chaque release candidate.  
3. **Surveillance numérique** : alerte automatique si résidu > 1.0e-10 (pré-alarme), rejet à 1.0e-7.  
4. **Extension L5** : ajouter cas multi-charnières intérieures et perturbations stochastiques des coordonnées.
