# Elasto-Gravity — J. Demangeot

Ce dossier contient les fichiers du projet **Elasto-Gravity**, une extension covariante de la Relativité Générale proposée par J. Demangeot (Institut d'Astrophysique de Paris).

## Fichiers

| Fichier | Description |
|---|---|
| `elastic_space.py` | Simulation Python de croissance auto-évitante sur le cercle unitaire |
| `Elasto_Gravity.tex` | Article LaTeX version courte |
| `Elasto_Gravity_Jd.tex` | Article LaTeX version complète (avec tenseur élastique détaillé) |
| `Elasto_Gravity_Jules_Demangeot.tex` | Version two-column style journal |

## Résumé

L'action du modèle est :

$$S = \int d^4x \sqrt{-g} \left[ \frac{M_{\rm Pl}^2}{2}R + \xi R\psi^2 - \frac{1}{2}(\partial\psi)^2 - \frac{1}{\Lambda^3}(\partial\psi)^2\Box\psi + \frac{\alpha}{M_{\rm Pl}^4}(\partial\psi)^4 \right] + S_m$$

## Résultats principaux

- **Tension $H_0$ résolue** : $H_0 = 73.0 \pm 1.0$ km/s/Mpc
- **Tension $\sigma_8$ résolue** : $\sigma_8 = 0.765 \pm 0.015$
- **DESI** : $w(a) = -0.948 - 0.295(1-a)$
- **LISA** : $\Delta\omega_R/\omega_R^{\rm GR} = -7.2\%$
- **EHT** : réduction de 5% du rayon d'ombre

## Auteur

J. Demangeot — Institut d'Astrophysique de Paris, CNRS, Sorbonne Université
