"""The driven deformation of the spanning-tree polynomial (Tutte's directed Laplacian).

Direction-dependent hopping weights forward hops by e^{+mu} and backward hops by
e^{-mu}. If the diagonal is adjusted so that every column of the matrix sums to zero,
the matrix is (up to sign and transposition) the generator of a biased random walk:
a DRIVEN deformation, not a chemical potential. (A lattice chemical potential leaves
the diagonal alone; see chemical_potential.py.) Tutte's directed matrix-tree theorem
writes its reduced determinant as a POSITIVE sum over arborescences. Directing each
spanning tree T away from a root r:

    f_mu(x) = sum_T  c_T(mu) * prod_{e in T} x_e,
    c_T(mu) = exp(mu * s_r(T)),

where s_r(T) counts edges T traverses "forward" (low label -> high label) minus
"backward". Positivity is automatic. The question is whether f_mu stays LORENTZIAN,
the property behind efficient sampling (Anari et al.).

Exact test (Braenden-Huh): a multiaffine polynomial with nonnegative coefficients and
M-convex support (here: spanning trees, always M-convex) is Lorentzian iff every
quadratic derivative d^S f, |S| = deg - 2, has a Hessian with at most one positive
eigenvalue. Braenden-Huh Theorem 3.14 turns "Lorentzian for every mu >= 0" into
"s_r is M-concave on spanning trees"; the last block of output checks both sides.
"""
import itertools

import networkx as nx
import numpy as np


def spanning_trees(G):
    edges = list(G.edges())
    V = G.number_of_nodes()
    for T in itertools.combinations(range(len(edges)), V - 1):
        H = nx.Graph([edges[i] for i in T])
        if H.number_of_nodes() == V and nx.is_tree(H):
            yield frozenset(T)


def coefficient(G, T, root, mu):
    edges = list(G.edges())
    H = nx.Graph([edges[i] for i in T])
    total = 0
    for parent, child in nx.bfs_edges(H, root):
        total += 1 if parent < child else -1
    return np.exp(mu * total)


def max_positive_eigs(G, root, mu):
    """Largest number of positive eigenvalues over all quadratic derivatives."""
    trees = {T: coefficient(G, T, root, mu) for T in spanning_trees(G)}
    d = G.number_of_nodes() - 1
    worst, witness = 0, None
    subsets = {S for T in trees for S in itertools.combinations(sorted(T), d - 2)}
    for S in subsets:
        S = frozenset(S)
        rest = sorted({e for T in trees if S <= T for e in T - S})
        idx = {e: i for i, e in enumerate(rest)}
        Q = np.zeros((len(rest), len(rest)))
        for T, c in trees.items():
            if S <= T:
                e, f = sorted(T - S)
                Q[idx[e], idx[f]] += c
                Q[idx[f], idx[e]] += c
        ev = np.linalg.eigvalsh(Q)
        pos = int((ev > 1e-9 * max(1.0, np.abs(ev).max())).sum())
        if pos > worst:
            worst, witness = pos, S
    return worst, witness, len(trees)


if __name__ == "__main__":
    graphs = {
        "C4": nx.cycle_graph(4),
        "C5": nx.cycle_graph(5),
        "K4": nx.complete_graph(4),
        "diamond": nx.Graph([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]),
        "K2,3": nx.complete_bipartite_graph(2, 3),
        "2x3 grid": nx.convert_node_labels_to_integers(nx.grid_2d_graph(2, 3)),
        "K5": nx.complete_graph(5),
    }
    def zero_mode_count(G, root=0):
        """Quadratic derivatives of f_0 with a zero eigenvalue (the only places
        where a small mu can create a violation, by continuity)."""
        trees = {T: coefficient(G, T, root, 0.0) for T in spanning_trees(G)}
        d = G.number_of_nodes() - 1
        subsets = {frozenset(S) for T in trees for S in itertools.combinations(sorted(T), d - 2)}
        degenerate = 0
        for S in subsets:
            rest = sorted({e for T in trees if S <= T for e in T - S})
            idx = {e: i for i, e in enumerate(rest)}
            Q = np.zeros((len(rest), len(rest)))
            for T, c in trees.items():
                if S <= T:
                    e, f = sorted(T - S)
                    Q[idx[e], idx[f]] += c
                    Q[idx[f], idx[e]] += c
            ev = np.linalg.eigvalsh(Q)
            degenerate += bool((np.abs(ev) < 1e-9 * np.abs(ev).max()).any())
        return degenerate, len(subsets)

    def tutte_check(G, root, mu, rng):
        """det of the reduced column-sum-zero directed Laplacian vs the tree sum f_mu."""
        nodes = sorted(G.nodes()); idx = {v: i for i, v in enumerate(nodes)}
        E = list(G.edges()); x = rng.uniform(0.5, 2.0, len(E))
        A = np.zeros((len(nodes), len(nodes)))
        for k, (u, v) in enumerate(E):
            A[idx[u], idx[v]] = x[k] * np.exp(mu * (1 if u < v else -1))   # arc u -> v
            A[idx[v], idx[u]] = x[k] * np.exp(mu * (1 if v < u else -1))   # arc v -> u
        L = np.diag(A.sum(axis=0)) - A                                     # columns sum to zero
        keep = [i for i in range(len(nodes)) if i != idx[root]]
        det = np.linalg.det(L[np.ix_(keep, keep)])
        tree_sum = sum(coefficient(G, T, root, mu) * np.prod([x[i] for i in T]) for T in spanning_trees(G))
        return abs(det - tree_sum) / abs(tree_sum)

    rng = np.random.default_rng(0)
    worst = max(tutte_check(G, r, mu, rng) for G in graphs.values() for r in sorted(G.nodes())[:2] for mu in (0.3, 1.0))
    print(f"f_mu equals the reduced determinant of the column-sum-zero directed Laplacian: max rel. error {worst:.1e}")
    print()

    print("Zero modes at mu = 0 (quadratic derivatives with a zero eigenvalue):")
    for name in ["C4", "C5", "K4", "diamond", "2x3 grid"]:
        k, tot = zero_mode_count(graphs[name])
        print(f"  {name:9s} {k} of {tot}")
    for mu in [1e-3, 1e-2]:
        print(f"K4 at mu = {mu}: max # positive eigenvalues = {max_positive_eigs(graphs['K4'], 0, mu)[0]}")
    print(f"K4 at e^mu = 2 (rational coefficients): max # positive eigenvalues = "
          f"{max_positive_eigs(graphs['K4'], 0, np.log(2))[0]}")

    # Zero modes at mu = 0 occur exactly when the contraction G/S has a parallel pair.
    def quadratic_forms(G, root, mu):
        trees = {T: coefficient(G, T, root, mu) for T in spanning_trees(G)}
        d = G.number_of_nodes() - 1
        for S in {frozenset(S) for T in trees for S in itertools.combinations(sorted(T), d - 2)}:
            rest = sorted({e for T in trees if S <= T for e in T - S})
            idx = {e: i for i, e in enumerate(rest)}
            Q = np.zeros((len(rest), len(rest)))
            for T, c in trees.items():
                if S <= T:
                    e, f = sorted(T - S)
                    Q[idx[e], idx[f]] += c
                    Q[idx[f], idx[e]] += c
            yield S, rest, Q

    def has_parallel_pair(G, S, rest):
        E = list(G.edges())
        H = nx.Graph()
        H.add_nodes_from(G.nodes())
        H.add_edges_from(E[i] for i in S)
        comp = {v: i for i, c in enumerate(nx.connected_components(H)) for v in c}
        pairs = [frozenset((comp[E[i][0]], comp[E[i][1]])) for i in rest]
        return len(pairs) != len(set(pairs))

    agree = total = 0
    for G in [graphs["C4"], graphs["C5"], graphs["K4"], graphs["K5"], graphs["diamond"],
              nx.house_graph(), graphs["2x3 grid"], nx.complete_bipartite_graph(3, 3)]:
        for S, rest, Q in quadratic_forms(G, 0, 0.0):
            ev = np.linalg.eigvalsh(Q)
            zero_mode = bool((np.abs(ev) < 1e-9 * np.abs(ev).max()).any())
            agree += zero_mode == has_parallel_pair(G, S, rest)
            total += 1
    print(f"Zero mode at mu = 0 <=> parallel pair in G/S: agrees on {agree} of {total} "
          "quadratic derivatives (C4, C5, K4, K5, diamond, house, 2x3 grid, K3,3)")

    # The K4 witness: the quadratic derivative obtained by contracting edge (2, 3).
    K4 = graphs["K4"]
    E4 = list(K4.edges())
    witness = frozenset({E4.index((2, 3))})
    print("K4, quadratic derivative contracting edge (2, 3), root 0:")
    for mu in [0.0, 0.05, 1.0]:
        Q = next(Q for S, rest, Q in quadratic_forms(K4, 0, mu) if S == witness)
        print(f"  mu = {mu:<5}: eigenvalues {np.round(np.linalg.eigvalsh(Q), 4)}")

    # The full Hessian of f_mu (not just a quadratic derivative) at random positive points.
    trees = {T: coefficient(K4, T, 0, 1.0) for T in spanning_trees(K4)}
    rng = np.random.default_rng(0)
    counts = {}
    for _ in range(2000):
        x = rng.uniform(0.01, 5, len(E4))
        H = np.zeros((len(E4), len(E4)))
        for T, c in trees.items():
            for e, f in itertools.combinations(sorted(T), 2):
                v = c * np.prod([x[g] for g in T if g not in (e, f)])
                H[e, f] += v
                H[f, e] += v
        ev = np.linalg.eigvalsh(H)
        k = int((ev > 1e-9 * np.abs(ev).max()).sum())
        counts[k] = counts.get(k, 0) + 1
    print(f"K4 full Hessian at mu = 1, 2000 random positive points: "
          f"# positive eigenvalues -> # points = {dict(sorted(counts.items()))}")
    print()

    mus = [0.0, 0.1, 0.5, 1.0, 2.0]
    print(f"{'graph':10s} root  " + "  ".join(f"mu={m:<4}" for m in mus) + "   (max # positive eigenvalues; Lorentzian iff all <= 1)")
    for name, G in graphs.items():
        for root in sorted(G.nodes())[:2]:
            row = [max_positive_eigs(G, root, m)[0] for m in mus]
            print(f"{name:10s} {root:4d}  " + "  ".join(f"{r:<7d}" for r in row))

    # Braenden-Huh, Theorem 3.14: sum_a q^{nu(a)} x^a / a! is Lorentzian for all 0 < q <= 1
    # iff nu is M-convex. Here q = e^{-mu} and nu = -s_r, so f_mu is Lorentzian for every
    # mu >= 0 iff s_r is M-concave on spanning trees (a valuated matroid). Check both sides.
    def s_r(G, T, root):
        return round(np.log(coefficient(G, T, root, 1.0)))  # c_T(mu = 1) = e^{s_r(T)}

    def m_concave(G, root):
        """Exchange axiom: for all trees A, B and a in A - B there is b in B - A with
        A - a + b, B - b + a trees and s(A) + s(B) <= s(A - a + b) + s(B - b + a)."""
        trees = list(spanning_trees(G))
        tree_set = set(trees)
        s = {T: s_r(G, T, root) for T in trees}
        for A in trees:
            for B in trees:
                for a in A - B:
                    if not any(
                        ((A - {a}) | {b}) in tree_set and ((B - {b}) | {a}) in tree_set
                        and s[A] + s[B] <= s[(A - {a}) | {b}] + s[(B - {b}) | {a}]
                        for b in B - A
                    ):
                        return False
        return True

    more = dict(graphs)
    more.update({
        "house": nx.house_graph(),
        "bowtie": nx.Graph([(0, 1), (1, 2), (0, 2), (2, 3), (3, 4), (2, 4)]),
        "C6": nx.cycle_graph(6),
        "K3,3": nx.complete_bipartite_graph(3, 3),
        "3x3 grid": nx.convert_node_labels_to_integers(nx.grid_2d_graph(3, 3)),
    })
    test_mus = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0]
    print()
    print(f"{'graph':10s} root  s_r M-concave   Lorentzian at all mu in {test_mus}")
    agree = total = 0
    for name, G in more.items():
        for root in sorted(G.nodes())[:2]:
            mc = m_concave(G, root)
            lor = all(max_positive_eigs(G, root, m)[0] <= 1 for m in test_mus)
            agree += mc == lor
            total += 1
            print(f"{name:10s} {root:4d}  {str(mc):15s} {lor}")
    print(f"M-concavity of s_r agrees with the eigenvalue test in {agree} of {total} cases")
