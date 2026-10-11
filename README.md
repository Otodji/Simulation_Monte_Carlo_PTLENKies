# Simulation Monte-Carlo de la loi Burr-PTLEN

Simulation, estimation par maximum de vraisemblance et tests d'adéquation pour la loi
**Burr-PTLEN**, une loi continue à 7 paramètres sur $x > 0$. Le notebook part d'un script R et
le traduit en Python. Il le complète par une log-densité numériquement stable, un quantile
explicite, des écarts-types, trois tests d'adéquation et un bootstrap paramétrique.

## La loi Burr-PTLEN

Paramètres : $m, n, p, \lambda, \beta, c, k > 0$. On pose

$$
A = (1 + x^c)^{-k}, \qquad G = \frac{1 - A}{A}, \qquad K(x) = 1 - e^{-\lambda G^{\beta}}, \qquad u = n\,(1 - K^{-p}) \le 0 .
$$

La fonction de répartition s'écrit alors :

$$
F(x) = e^{m u}\,\left(2 - e^{u}\right)^{m}
$$

$A$ est la fonction de survie de la loi de Burr XII.

## Contenu du notebook

| # | Section | Ce qui est fait |
|---|---|---|
| 1 | Bissection | Recherche de racine sur $[a, b]$, utilisée pour inverser $F$ comme dans le script R |
| 2 | Fonction de répartition et densité | $F$ et $f$ vectorisées avec NumPy. **Ajout** : une log-densité calculée terme à terme (`log1p`, `expm1`), pour que la log-vraisemblance reste finie dans les queues de la loi. Vérification que $\int f = 1$ et que $f = F'$ |
| 3 | Simulation par inversion | Méthode de la transformée inverse $X = F^{-1}(U)$ par bissection. **Ajout** : un quantile **explicite** $F^{-1}$, exact et instantané, qui sert à valider la bissection. Simulation de 250 observations, puis comparaison de l'histogramme, de la densité estimée et de la densité théorique |
| 4 | Maximum de vraisemblance | Minimisation de $-\log L$ avec L-BFGS-B, chaque paramètre étant borné dans $[10^{-6}, 10^{3}]$. **Ajout** : écarts-types tirés de l'inverse de la hessienne, l'équivalent de `vcov()` en R |
| 5 | Tests d'adéquation | Kolmogorov-Smirnov, Cramér-von Mises et Anderson-Darling, avec leurs p-values. Comparaison graphique de la fonction de répartition empirique et de la loi ajustée |
| 6 | p-values par bootstrap paramétrique | Les paramètres sont estimés sur les mêmes données que celles testées, donc les p-values classiques sont trop grandes. Le bootstrap paramétrique les corrige : on simule, on ré-estime, puis on recalcule les statistiques |

### Correspondance R → Python

| R | Python |
|---|---|
| `stats4::mle` (L-BFGS-B) | `scipy.optimize.minimize(method="L-BFGS-B")` |
| `runif`, `set.seed(123)` | `numpy.random.default_rng(123).uniform` |
| `tibble`, `dplyr`, `ggplot2` | `numpy`, `matplotlib` |

## Points à retenir

- **La simulation par inversion est validée** : l'histogramme et la densité théorique se superposent,
  et la bissection coïncide avec le quantile exact.
- **Identifiabilité** : avec 7 paramètres et 250 observations, plusieurs jeux de paramètres
  donnent des courbes presque identiques. Les estimations peuvent donc s'éloigner des vraies
  valeurs, ou finir sur une borne, tout en ajustant très bien les données. Les écarts-types
  sont alors grands.
- **Ce sont les p-values bootstrap qu'il faut retenir**, et non les p-values classiques des tests
  KS, CvM et AD.

## Installation et utilisation

```bash
git clone https://github.com/Otodji/simulation-burr-ptlen.git
cd simulation-burr-ptlen
pip install -r requirements.txt
jupyter notebook Simulation_BurrPTLEN.ipynb
```

Le notebook est enregistré **sans résultats d'exécution**. Lance *Run All* pour obtenir les
figures et les estimations. La graine est fixée (`default_rng(123)`), donc les résultats sont
reproductibles. La section 6 (bootstrap) est la plus longue à exécuter.

## Auteur

[OTODJI Théophile](https://github.com/Otodji)
