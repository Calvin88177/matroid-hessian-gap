"""A lattice chemical potential on a Laplacian-type fermion matrix.

    M(mu) = D - W(mu),   W_uv = x_e * exp(mu * t_uv),

with D the ordinary degree matrix (UNCHANGED by mu, as for a lattice chemical potential
in the sense of Hasenfratz-Karsch) and t_uv = +1 for a forward time-like hop, -1 for a
backward one, 0 for a spatial one. (Tutte's version, which also changes the diagonal,
is the driven deformation of driven_deformation.py.)

Kenyon (Ann. Probab. 2011): det M(mu) is a sum over cycle-rooted spanning forests F
(every component contains exactly one cycle) of

    prod_{e in F} x_e * prod_{cycles C of F} (2 - exp(mu * l_C) - exp(-mu * l_C)),

where l_C is the net number of forward-minus-backward time-like steps around C, i.e.
L times its winding number around periodic time. With a mass term, det(m + M) adds
rooted trees, each with weight m * (number of its vertices).

  real mu:       2 - 2 cosh(mu l) = -4 sinh^2(mu l / 2) < 0  -> negative weights
  imaginary mu:  2 - 2 cos(theta l) = 4 sin^2(theta l / 2) >= 0  -> sign-free

This script checks the expansion against the determinant, the Cauchy-Binet form of the
coefficients, the average sign with a mass term, and the Lorentzian property at
imaginary mu (exact Braenden-Huh test).
"""
import itertools

import numpy as np


def periodic_lattice(L, k):
    """C_L x P_k: k periodic time-loops of length L joined by spatial rungs.
    Returns (number of vertices, edges as (u, v, t_uv))."""
    vid = lambda s, i: s * L + i
    edges = [(vid(s, i), vid(s, (i + 1) % L), 1) for s in range(k) for i in range(L)]
    edges += [(vid(s, i), vid(s + 1, i), 0) for s in range(k - 1) for i in range(L)]
    return k * L, edges


def fermion_matrix(V, edges, x, mu, m=0.0):
    M = np.zeros((V, V), dtype=complex) + m * np.eye(V)
    for xe, (u, v, t) in zip(x, edges):
        M[u, u] += xe
        M[v, v] += xe
        M[u, v] -= xe * np.exp(mu * t)
        M[v, u] -= xe * np.exp(-mu * t)
    return M


def forest_weight(V, edges, F, x, mu, m):
    """Weight of the edge set F in the expansion of det(m + M(mu)); 0 if F is not a
    disjoint union of trees and unicyclic components. Union-find carries a time
    potential so that closing a cycle reveals its net time displacement l_C."""
    parent = list(range(V))
    off = [0] * V            # potential of a vertex minus that of its parent
    size = [1] * V
    nedge = [0] * V
    cyc = [None] * V         # net displacement of the component's cycle, if any

    def find(v):
        if parent[v] == v:
            return v, 0
        r, o = find(parent[v])
        parent[v] = r
        off[v] += o
        return r, off[v]

    w = 1.0 + 0j
    for e in F:
        u, v, t = edges[e]
        w *= x[e]
        ru, pu = find(u)
        rv, pv = find(v)
        if ru == rv:
            if cyc[ru] is not None:
                return 0.0
            cyc[ru] = pu + t - pv
            nedge[ru] += 1
        else:
            if cyc[ru] is not None and cyc[rv] is not None:
                return 0.0
            if size[ru] < size[rv]:
                ru, rv, pu, pv, t = rv, ru, pv, pu, -t
            parent[rv] = ru
            off[rv] = t + pu - pv
            size[ru] += size[rv]
            nedge[ru] += nedge[rv] + 1
            if cyc[ru] is None:
                cyc[ru] = cyc[rv]
    for r in range(V):
        if parent[r] == r:
            if cyc[r] is None:
                w *= m * size[r]
            else:
                h = np.exp(mu * cyc[r])
                w *= 2 - h - 1 / h
    return w


def expansion(V, edges, x, mu, m, sizes):
    out = {}
    for r in sizes:
        for F in itertools.combinations(range(len(edges)), r):
            w = forest_weight(V, edges, F, x, mu, m)
            if abs(w) > 1e-12:
                out[frozenset(F)] = w
    return out


def m_convex(J):
    J = set(J)
    for A in J:
        for B in J:
            for a in A - B:
                if not any(((A - {a}) | {b}) in J and ((B - {b}) | {a}) in J for b in B - A):
                    return False
    return True


def max_positive_eigs(terms, d):
    worst = 0
    for S in {frozenset(S) for T in terms for S in itertools.combinations(sorted(T), d - 2)}:
        rest = sorted({e for T in terms if S <= T for e in T - S})
        idx = {e: i for i, e in enumerate(rest)}
        Q = np.zeros((len(rest), len(rest)))
        for T, c in terms.items():
            if S <= T:
                e, f = sorted(T - S)
                Q[idx[e], idx[f]] += c
                Q[idx[f], idx[e]] += c
        ev = np.linalg.eigvalsh(Q)
        worst = max(worst, int((ev > 1e-9 * np.abs(ev).max()).sum()))
    return worst


if __name__ == "__main__":
    rng = np.random.default_rng(0)

    print("1. Kenyon's expansion vs the determinant (massless, random edge weights), ladders C_L x P_2")
    for L in [3, 4]:
        V, E = periodic_lattice(L, 2)
        x = rng.uniform(0.5, 2.0, len(E))
        for label, mu in [("real mu = 0.4", 0.4), ("imaginary mu = 0.4i", 0.4j)]:
            terms = expansion(V, E, x, mu, 0.0, [V])
            det = np.linalg.det(fermion_matrix(V, E, x, mu))
            vals = np.array(list(terms.values()))
            neg = int((vals.real < 0).sum())
            print(f"   L={L} {label:21s} det = {det.real:+11.4f}   forest sum = {sum(vals).real:+11.4f}   "
                  f"nonzero forests = {len(vals)}, negative = {neg}")

    print("\n2. Cauchy-Binet: M = sum_e x_e v_e w_e^T, so each coefficient is det(V_F) det(W_F)")
    V, E = periodic_lattice(3, 2)
    x = np.ones(len(E))
    for label, mu in [("real mu = 0.4", 0.4), ("imaginary mu = 0.4i", 0.4j)]:
        Vm = np.zeros((V, len(E)), dtype=complex)
        Wm = np.zeros((V, len(E)), dtype=complex)
        for k, (u, v, t) in enumerate(E):
            Vm[u, k], Vm[v, k] = 1, -np.exp(-mu * t)
            Wm[u, k], Wm[v, k] = 1, -np.exp(mu * t)
        terms = expansion(V, E, x, mu, 0.0, [V])
        err = 0.0
        for F in itertools.combinations(range(len(E)), V):
            cb = np.linalg.det(Vm[:, F]) * np.linalg.det(Wm[:, F])
            err = max(err, abs(cb - terms.get(frozenset(F), 0.0)))
        conj = np.allclose(Wm, np.conj(Vm))
        print(f"   {label:21s} max |det V_F det W_F - forest weight| = {err:.1e};  W = conj(V): {conj}")

    print("\n3. Average sign  <s> = Z / sum |weights|  with mass m = 0.5, uniform weights")
    for L, k in [(3, 2), (4, 2), (5, 2), (6, 2), (3, 3), (4, 3), (3, 4)]:
        V, E = periodic_lattice(L, k)
        x = np.ones(len(E))
        row = []
        for mu in [0.3, 0.6]:
            terms = expansion(V, E, x, mu, 0.5, range(len(E) + 1))
            Z = sum(terms.values()).real
            det = np.linalg.det(fermion_matrix(V, E, x, mu, 0.5)).real
            assert abs(Z - det) < 1e-8 * max(1.0, abs(det)), (L, k, mu, Z, det)
            row.append(f"mu={mu}: {Z / sum(abs(w) for w in terms.values()):.4f}")
        print(f"   time L={L}, space k={k}:   " + "   ".join(row))

    print("\n4. Imaginary mu = i*theta (massless): is the nonnegative forest polynomial Lorentzian?")
    for L in [3, 4]:
        V, E = periodic_lattice(L, 2)
        x = np.ones(len(E))
        for theta in [0.4, 1.0]:
            terms = {F: w.real for F, w in expansion(V, E, x, 1j * theta, 0.0, [V]).items()}
            assert min(terms.values()) > 0
            print(f"   L={L} theta={theta}: {len(terms)} forests, all weights > 0, "
                  f"support M-convex: {m_convex(terms.keys())}, "
                  f"max # positive eigenvalues over quadratic derivatives: {max_positive_eigs(terms, V)}")
