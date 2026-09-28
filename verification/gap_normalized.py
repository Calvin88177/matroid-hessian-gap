"""Normalized gap gamma = (n-1) * min|lambda_-| / lambda_+  in (0,1].
Trace(Hess)=0 forces min|lambda_-|/lambda_+ <= 1/(n-1); gamma removes that trivial scaling.
gamma = 1 iff all negative eigenvalues are equal."""
import numpy as np, networkx as nx, csv
from spectral_gap import hess_normalized, families
rng = np.random.default_rng(2)

def gamma(G, x):
    ev = np.linalg.eigvalsh(hess_normalized(G, x))
    lp, neg = ev[-1], ev[:-1]
    assert (neg < 0).all() and lp > 0
    return (len(ev) - 1) * np.abs(neg).min() / lp

def sample(G, K, s=20):
    n = G.number_of_edges()
    return min(gamma(G, np.exp(rng.uniform(-np.log(K), np.log(K), n))) for _ in range(s))

rows = []
for fam, graphs in families.items():
    for p, G in graphs:
        n = G.number_of_edges()
        rows.append((fam, p, n, gamma(G, np.ones(n)), sample(G, 2), sample(G, 10), sample(G, 100)))
with open("gap_normalized.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["family","param","edges","gamma_uniform","gamma_min_ratio2","gamma_min_ratio10","gamma_min_ratio100"]); w.writerows(rows)
print(f"{'family':24s} {'p':>3s} {'|E|':>5s} {'x=1':>7s} {'ratio<=2':>9s} {'<=10':>7s} {'<=100':>8s}")
for r in rows:
    print(f"{r[0]:24s} {r[1]:3d} {r[2]:5d} {r[3]:7.3f} {r[4]:9.3f} {r[5]:7.3f} {r[6]:8.2e}")
