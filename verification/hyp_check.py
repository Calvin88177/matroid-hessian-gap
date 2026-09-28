"""Hypothesis check: when is the Hessian of the spanning-tree polynomial f_G
non-degenerate with signature (1, n-1) on the positive orthant?

f_G is built with the weighted matrix-tree theorem, so multigraphs and loops
are handled correctly (a loop never lies in a spanning tree).
"""
import numpy as np
import sympy as sp

rng = np.random.default_rng(0)


def tree_poly(num_vertices, edges):
    xs = sp.symbols(f"x0:{len(edges)}", positive=True)
    L = sp.zeros(num_vertices, num_vertices)
    for x, (u, v) in zip(xs, edges):
        if u == v:
            continue  # loop: contributes nothing
        L[u, u] += x; L[v, v] += x
        L[u, v] -= x; L[v, u] -= x
    f = sp.expand(L[1:, 1:].det()) if num_vertices > 1 else sp.Integer(1)
    return xs, f


def signature(num_vertices, edges, samples=200, tol=1e-9):
    xs, f = tree_poly(num_vertices, edges)
    H = sp.hessian(f, xs)
    Hf = sp.lambdify(xs, H, "numpy")
    n = len(edges)
    sigs, min_abs = set(), np.inf
    for _ in range(samples):
        pt = rng.uniform(0.05, 5.0, n)
        ev = np.linalg.eigvalsh(np.array(Hf(*pt), dtype=float))
        scale = max(1.0, np.abs(ev).max())
        pos = int((ev > tol * scale).sum())
        neg = int((ev < -tol * scale).sum())
        sigs.add((pos, neg, n - pos - neg))
        min_abs = min(min_abs, (np.abs(ev) / scale).min())
    return n, sorted(sigs), min_abs


def cycle(k):
    return [(i, (i + 1) % k) for i in range(k)]


def complete(k):
    return [(i, j) for i in range(k) for j in range(i + 1, k)]


cases = {
    # the five graphs already in the README
    "K3": (3, complete(3)),
    "K4": (4, complete(4)),
    "C5": (5, cycle(5)),
    "diamond": (4, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]),
    "house": (5, [(0, 1), (1, 2), (2, 3), (3, 0), (2, 4), (3, 4)]),
    # edge cases
    "single edge (rank 1)": (2, [(0, 1)]),
    "path P3 (tree, 2 bridges)": (3, [(0, 1), (1, 2)]),
    "star K1,3 (tree)": (4, [(0, 1), (0, 2), (0, 3)]),
    "K3 + pendant edge (bridge)": (4, complete(3) + [(2, 3)]),
    "two triangles joined by bridge": (6, complete(3) + [(3, 4), (4, 5), (3, 5), (2, 3)]),
    "bowtie (cut vertex, no bridge)": (5, complete(3) + [(2, 3), (3, 4), (2, 4)]),
    "K3 with a doubled edge (parallel)": (3, complete(3) + [(0, 1)]),
    "K3 with a loop": (3, complete(3) + [(0, 0)]),
    "K5": (5, complete(5)),
    "K2,3": (5, [(i, j) for i in range(2) for j in range(2, 5)]),
}

print(f"{'graph':38s} {'n':>3s}  {'target':>9s}  observed (pos,neg,zero)   min|eig|/max")
for name, (V, E) in cases.items():
    n, sigs, m = signature(V, E, samples=60 if len(E) > 10 else 200)
    ok = sigs == [(1, n - 1, 0)]
    print(f"{name:38s} {n:3d}  {str((1, n-1)):>9s}  {str(sigs):26s} {m:.2e}  {'OK' if ok else 'degenerate (outside hypotheses)'}")
