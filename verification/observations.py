"""Reproduces the numerical observations in ROADMAP.md, section 2.

1. Closed forms at uniform weights (observed, not proved):
     gamma(K_n)     = (n+1)/(2n)
     gamma(K_{m,m}) = (m+1)/(2m-1)
     gamma(C_n)     = 1          (proved: Hessian = (n-2)(J - I))
2. Hypercube values, printed as fractions.
3. Exhaustive scan of all connected graphs on 3..7 vertices (networkx graph
   atlas): minimum of gamma(G, 1) for each vertex count.
"""
from fractions import Fraction

import networkx as nx
import numpy as np
from networkx.generators.atlas import graph_atlas_g

from spectral_gap import hess_normalized


def gamma_uniform(G):
    n = G.number_of_edges()
    ev = np.linalg.eigvalsh(hess_normalized(G, np.ones(n)))
    return (n - 1) * np.abs(ev[:-1]).min() / ev[-1]


def check(name, graphs, formula):
    worst = max(abs(gamma_uniform(G) - formula(p)) for p, G in graphs)
    print(f"{name:12s} max |gamma - formula| over tested sizes: {worst:.1e}")


check("K_n", [(n, nx.complete_graph(n)) for n in range(4, 46)], lambda n: (n + 1) / (2 * n))
check("K_{m,m}", [(m, nx.complete_bipartite_graph(m, m)) for m in range(2, 26)], lambda m: (m + 1) / (2 * m - 1))
check("C_n", [(n, nx.cycle_graph(n)) for n in range(3, 201, 7)], lambda n: 1.0)

ev = np.linalg.eigvalsh(hess_normalized(nx.complete_graph(8), np.ones(28)))
print(f"K_8 at uniform weights: distinct eigenvalues of Hess f / f = {np.unique(np.round(ev, 10))}")

print("\nHypercube Q_d at uniform weights:")
for d in range(2, 10):
    g = gamma_uniform(nx.convert_node_labels_to_integers(nx.hypercube_graph(d)))
    print(f"  d={d}: {g:.10f}  ~ {Fraction(g).limit_denominator(2000)}")

print("\nMinimum gamma(G, 1) over all connected graphs on V vertices:")
best = {}
for G in graph_atlas_g():
    V = G.number_of_nodes()
    if V < 3 or not nx.is_connected(G):
        continue
    g = gamma_uniform(G)
    if g < best.get(V, (np.inf,))[0]:
        best[V] = (g, G.number_of_edges())
for V, (g, E) in sorted(best.items()):
    print(f"  V={V}: {g:.4f}  (a minimizer has {E} edges)")
