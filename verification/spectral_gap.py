"""Spectral gap of the Kirchhoff-polynomial Hessian along graph families.

Closed form (matrix-tree theorem + Jacobi's formula). With conductances x > 0,
let L be the weighted Laplacian with one row/column removed, b_e the signed
incidence vector of edge e (restricted), and

    T[e,f] = b_e^T L^{-1} b_f          (transfer-current matrix),
    R_e    = T[e,e]                    (effective resistance of e).

Then  d_e f = f * R_e  and  d_e d_f f = f * (R_e R_f - T[e,f]^2),
so    Hess f(x) = f(x) * (R R^T - T∘T)     (∘ = entrywise product).

The diagonal vanishes, as it must for a multiaffine polynomial. The scalar f(x)
does not affect signature or eigenvalue ratios, so we work with
M(x) = R R^T - T∘T.
"""
import sys
import numpy as np
import networkx as nx

rng = np.random.default_rng(1)


def hess_normalized(G, x):
    nodes = list(G.nodes())
    idx = {v: i for i, v in enumerate(nodes)}
    edges = list(G.edges())
    V, n = len(nodes), len(edges)
    B = np.zeros((V, n))
    for k, (u, v) in enumerate(edges):
        B[idx[u], k] = 1.0
        B[idx[v], k] = -1.0
    B = B[1:]                               # ground vertex 0
    L = (B * x) @ B.T
    T = B.T @ np.linalg.solve(L, B)
    R = np.diag(T)
    return np.outer(R, R) - T * T


def spectrum_stats(G, x):
    ev = np.linalg.eigvalsh(hess_normalized(G, x))
    top = np.abs(ev).max()
    pos = int((ev > 1e-10 * top).sum())
    neg = int((ev < -1e-10 * top).sum())
    gap = np.abs(ev).min() / top
    return pos, neg, len(ev), gap


# ---- 1. validate closed form against the symbolic results --------------------
def validate():
    import sympy as sp
    for name, G in [("K4", nx.complete_graph(4)), ("house", nx.house_graph()),
                    ("K2,3", nx.complete_bipartite_graph(2, 3))]:
        edges = list(G.edges())
        xs = sp.symbols(f"x0:{len(edges)}", positive=True)
        nodes = list(G.nodes())
        Lsym = sp.zeros(len(nodes))
        for s, (u, v) in zip(xs, edges):
            i, j = nodes.index(u), nodes.index(v)
            Lsym[i, i] += s; Lsym[j, j] += s; Lsym[i, j] -= s; Lsym[j, i] -= s
        f = sp.expand(Lsym[1:, 1:].det())
        pt = rng.uniform(0.2, 3.0, len(edges))
        sub = dict(zip(xs, pt))
        H_sym = np.array(sp.hessian(f, xs).subs(sub), dtype=float)
        H_num = float(f.subs(sub)) * hess_normalized(G, pt)
        err = np.abs(H_sym - H_num).max() / np.abs(H_sym).max()
        print(f"  validate {name:6s}: rel. error {err:.1e}")


# ---- 2. families ------------------------------------------------------------
def tri_grid(m):
    return nx.convert_node_labels_to_integers(nx.triangular_lattice_graph(m, m))


families = {
    "complete K_n":            [(k, nx.complete_graph(k)) for k in range(3, 23)],
    "bipartite K_{m,m}":       [(m, nx.complete_bipartite_graph(m, m)) for m in range(2, 16)],
    "cycle C_n":               [(k, nx.cycle_graph(k)) for k in range(3, 121, 6)],
    "square grid m x m":       [(m, nx.convert_node_labels_to_integers(nx.grid_2d_graph(m, m))) for m in range(2, 15)],
    "triangular grid (m,m)":   [(m, tri_grid(m)) for m in range(1, 11)],
    "hypercube Q_d":           [(d, nx.convert_node_labels_to_integers(nx.hypercube_graph(d))) for d in range(2, 8)],
}


def run(samples=20):
    rows = []
    for fam, graphs in families.items():
        for p, G in graphs:
            n = G.number_of_edges()
            s_uni = spectrum_stats(G, np.ones(n))
            gaps, sig_ok = [], s_uni[:2] == (1, n - 1)
            for _ in range(samples):
                s = spectrum_stats(G, rng.uniform(0.5, 2.0, n))
                gaps.append(s[3])
                sig_ok &= s[:2] == (1, n - 1)
            rows.append((fam, p, G.number_of_nodes(), n, s_uni[3], min(gaps), sig_ok))
    return rows


if __name__ == "__main__":
    print("Validating closed form against symbolic Hessian:")
    validate()
    rows = run()
    import csv
    with open("spectral_gap.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["family", "param", "vertices", "edges", "gap_uniform", "gap_min_random", "signature_1_n-1"])
        w.writerows(rows)
    print(f"\n{'family':24s} {'p':>3s} {'|V|':>4s} {'|E|':>5s} {'gap(x=1)':>10s} {'min gap rand':>12s} sig")
    for r in rows:
        print(f"{r[0]:24s} {r[1]:3d} {r[2]:4d} {r[3]:5d} {r[4]:10.3e} {r[5]:12.3e} {'ok' if r[6] else 'FAIL'}")
