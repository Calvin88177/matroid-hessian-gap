I want to set down some observations about the spanning-tree polynomial of a graph, together with the questions they have led me to. Some of what follows is proved, some is only computed, and some is speculation. I have tried to say which is which each time. Paragraphs 1–3 recall what is known, paragraph 4 explains why I care, paragraphs 5–10 contain what I think is new, and paragraph 11 is the speculative one. If that last part does not interest you, the rest should stand without it.

---

**1. The object.** Let G be a connected graph with edge weights x_e > 0 on its n edges, and let

$$f_G(x) = \sum_{T} \prod_{e \in T} x_e$$

be the sum over its spanning trees. By Kirchhoff, f_G is the determinant of the weighted Laplacian L with one row and column removed. So it is also a Gaussian Grassmann integral:

$$f_G(x) = \det L_r(x) = \int \prod_i d\bar\psi_i\, d\psi_i\; e^{\bar\psi L_r(x) \psi}$$

(up to sign conventions). In other words, f_G is a free-fermion partition function that happens to be a positive sum of combinatorial objects. Caracciolo, Jacobsen, Saleur, Sokal and Sportiello showed that this extends beyond the Gaussian case: adding a particular four-fermion term produces spanning forests.

**2. What is known about its Hessian.** Brändén and Huh showed that f_G is *Lorentzian*. Its support, the set of spanning trees, is M-convex, and its Hessian has exactly one positive eigenvalue at every point of the positive orthant. Their theorem leaves room for zero eigenvalues. Nagaoka and Yazawa, and then Murai, Nagaoka and Yazawa for all simple matroids of rank at least 2, showed there are none. The Hessian has signature exactly (1, n−1). As I understand their argument, non-degeneracy reduces to the linear independence of the partial derivatives ∂_e f_G, which is where simplicity enters.

The hypotheses are sharp, and the reason is transparent. A loop gives ∂_e f = 0, and a parallel pair gives ∂_e f = ∂_{e'} f. Either one is a linear dependency, and a zero eigenvalue follows. Bridges do no harm. I checked all of this on fifteen small graphs before I understood why. I mention it because the same degeneracy returns in paragraph 7, where it does real work.

**3. The Hessian as a correlation matrix.** Let b_e be the signed incidence vector of edge e. Let T be the transfer-current matrix T_{ef} = b_eᵀ L_r⁻¹ b_f, and let R_e = T_{ee} be the effective resistance across e. Then

$$\nabla^2 f_G = f_G\,\big(RR^{\top} - T \circ T\big),$$

where ∘ is the entrywise product. Equivalently, if 𝒯 is a random spanning tree chosen with probability proportional to its weight,

$$\frac{\partial_e \partial_f f_G}{f_G} = R_e R_f - T_{ef}^2 = \frac{\mathbb{P}(e, f \in \mathcal{T})}{x_e\,x_f}, \qquad e \neq f,$$

which is the Burton–Pemantle transfer-current theorem. Up to scaling, then, the Hessian is the matrix of pairwise inclusion probabilities of a determinantal process: the combinatorial shadow of a free-fermion density correlation. None of this is new, but it makes the Hessian computable with one linear solve. My scripts agree with the symbolic Hessian to about 10⁻¹⁶.

---

**4. Why I care: the sign problem.** A quantum Monte Carlo method for fermions needs a representation of the weight that is both *positive* and *efficiently sampleable*. The sign problem is the failure of the first. Troyer and Wiese showed that it is NP-hard in general. They noted that this does not exclude solving it for a restricted class of systems, but that any such method must be tied to properties of that class. So the productive question is which classes admit both properties, and why.

Spanning trees are the simplest class where both hold visibly. Kirchhoff gives the positivity. Lorentzian structure (log-concavity) is the mechanism behind efficient sampling in the matroid setting: it is what lets random walks on matroid bases mix rapidly (Anari, Liu, Oveis Gharan, Vinzant).

I should be plain about the limits. These fermions never had a sign problem. Even in the deformed setting below, the weights stay positive, the partition function is a determinant, and weighted arborescences can still be sampled exactly by Wilson's algorithm. The physically hard cases, such as doped Hubbard models or QCD at finite baryon density, involve fermion matrices of a different kind, and nothing here reaches them. My interest is structural: *which deformations keep positivity together with log-concavity, and which keep one but not the other?*

**5. The deformation: a chemical potential.** The deformation that matters most physically is a chemical potential. It does not rescale hopping symmetrically. It weights forward and backward hopping differently, by e^{+μ} and e^{−μ}. For a symmetric rescaling x_e ↦ w_e x_e the question is empty: that is a nonnegative linear change of variables, which preserves the Lorentzian property by Brändén–Huh. I had once taken that fact as the answer to the chemical-potential question. It is not.

For the asymmetric deformation, Tutte's directed matrix-tree theorem writes the reduced determinant of the directed Laplacian as a positive sum over arborescences. Directing each spanning tree away from a root r,

$$f_\mu(x) = \sum_T e^{\mu\, s_r(T)} \prod_{e \in T} x_e,$$

where s_r(T) is the number of forward edges minus the number of backward edges. The support is still the set of spanning trees, so it is still M-convex. By Brändén–Huh, f_μ is Lorentzian *if and only if every quadratic derivative* ∂^S f_μ (with |S| = deg − 2) has a Hessian with at most one positive eigenvalue. This makes the question exact and finite for each graph. The quadratic derivative ∂^S f_μ is the tree polynomial of the contraction G/S, a graph on three vertices, with its coefficients inherited from the orientations.

For my first experiments, "forward" means from lower to higher vertex label, on every edge. This is a toy. I return to the faithful version in paragraph 8.

**6. What happens.** I expected the Lorentzian property to survive. It does not, at least not always. Here are the results of the exact test for μ ∈ {0.1, 0.5, 1, 2}, and for K₄ also for μ ∈ {0.001, 0.01}:

| Graph | Root 0 | Root 1 | Triangles |
|---|---|---|---|
| K₄ | **fails** at every μ tested, down to 0.001 | **fails** | 4 |
| K₅ | **fails** | **fails** | 10 |
| Diamond (K₄ minus an edge) | survives | **fails** | 2 |
| C₄ | survives | survives | 0 |
| C₅ | survives | survives | 0 |
| 2×3 grid | survives | survives | 0 |
| K₂,₃ | survives, trivially | survives, trivially | 0 |

On K₄ the failure is not confined to a quadratic derivative. At μ = 1, at 62 of 2,000 random positive points, the full Hessian of f_μ has *two* positive eigenvalues, so the Lorentzian signature itself is lost.

K₂,₃ survives for a reason worth recording, because it shows what a trivial survival looks like. Its labels put one side {0, 1} entirely below the other side {2, 3, 4}, so every edge points from the first side to the second. In a tree directed away from a root on the first side, each vertex on the second side is entered forward, and each non-root vertex on the first side is entered backward. So s_r(T) is the same for every tree, and f_μ is a constant multiple of f_0. The same argument works for a root on the other side, and for any bipartite graph whose labels respect the bipartition in this way.

**7. The mechanism.** Two observations explain the table, one of them a proof.

*(a) Where violations can come from (proved).* Each quadratic derivative is a symmetric matrix depending continuously on μ. If at μ = 0 it has one positive eigenvalue and no zero eigenvalue, it keeps that signature for all small μ. So f_μ is Lorentzian for all sufficiently small μ unless some quadratic derivative of f_0 has a zero eigenvalue. Call these the *zero modes*. The cycles have none (0 of 4 for C₄, 0 of 10 for C₅), so their survival at small μ is automatic. Their survival up to μ = 2 is observed.

*(b) Where the zero modes come from (proved).* At μ = 0 every coefficient is 1, so the Hessian of a quadratic derivative has entry 1 for two edges of G/S that join different pairs of its three vertices, and 0 otherwise. Edges joining the same pair give identical rows, so the rank is at most 3, and there is a zero eigenvalue exactly when G/S has a pair of parallel edges. This is the degeneracy of paragraph 2 again. (It agrees with direct computation on all 204 quadratic derivatives of eight small graphs.) Contracting one edge of a triangle produces a parallel pair; on the grid, contracting two edges of a 4-cycle does. On K₄, every one of the six quadratic derivatives has a zero mode. Take the one obtained by contracting edge {2, 3}:

| μ | Eigenvalues of that quadratic derivative |
|---|---|
| 0 | −2, −1.236, 0, 0, 3.236 |
| 0.05 | −2.005, −1.431, 0, **+0.0025**, 3.434 |
| 1 | −27.21, −3.81, 0, **+0.72**, 30.30 |

A nonzero μ weights the parallel pair differently and pushes one zero eigenvalue positive, at order μ²: 0.0025 = 0.05².

So zero modes are *necessary* for failure at small μ, by (a), but not *sufficient*. The 2×3 grid has zero modes in 27 of its 35 quadratic derivatives, and it survives at every μ tested. What decides is the *direction* in which μ pushes each zero mode. In the language of degenerate perturbation theory: if Q(μ) = Q₀ + μQ₁ + μ²Q₂ + … and P projects onto the kernel of Q₀, the first-order effective matrix is P Q₁ P. When that vanishes, the second-order one is P(Q₂ − Q₁ Q₀⁺ Q₁)P, where Q₀⁺ is the pseudoinverse. The observed μ² scaling on K₄ suggests the first-order term vanishes there. For small μ, a zero mode turns positive when the first non-vanishing effective matrix has a positive eigenvalue. I think this is the right way to attack the problem, but I have not carried it out.

**8. The questions.**

> *(A) For which graphs, roots and orientations does f_μ remain Lorentzian for all μ? Is it enough for the graph to be triangle-free?*

Every failure I have seen is on a graph with triangles, and every triangle-free graph I have tried survives. But triangles do not force failure: the diamond survives for one root and fails for the other.

> *(B) Is there a clean criterion for the sign of the effective perturbation in paragraph 7 on a triangle's zero mode, and on the zero modes that come from longer cycles?*

Answering (B) for triangles would explain the diamond, where the same triangles lead to failure for one root and not for the other. The grid, whose zero modes come from 4-cycles and never turned positive in my tests, suggests longer cycles behave differently.

> *(C) Does the Lorentzian structure survive a physically faithful chemical potential?*

Here the lattice is time-layered, G = C_L × H, with time periodic, and e^{±μ} sits only on the time-like edges. In lattice field theory, the μ-dependence of a fermion determinant enters through paths that wind around periodic time. I would like to know whether the failure of the Lorentzian property tracks winding in the same way. Net direction alone cannot be the whole story: under my toy orientation, C₄ and C₅ have nonzero net direction around the cycle, and both survive.

I have not found (A)–(C) treated in the literature, but I have not yet searched as carefully as I must before calling them open.

---

**9. The size of the margin.** A second, quieter set of questions concerns how non-degenerate the Hessian is. Since the Hessian has zero diagonal, its trace vanishes, and min|λ₋|/λ₊ ≤ 1/(n−1) for every graph. The natural measure is therefore

$$\gamma = (n-1)\,\frac{\min|\lambda_-|}{\lambda_+} \in (0, 1],$$

with γ = 1 exactly when all negative eigenvalues are equal.

At uniform weights the numbers fit clean formulas, to about 10⁻¹²:

$$\gamma(K_n) = \frac{n+1}{2n}, \qquad \gamma(K_{m,m}) = \frac{m+1}{2m-1},$$

both tending to ½. For cycles γ = 1 exactly: every pair of edges lies together in n − 2 spanning trees, so the Hessian is (n − 2)(J − I). For K_n the Hessian commutes with the edge symmetries and lies in the Bose–Mesner algebra of the Johnson scheme J(n, 2). K₈, for instance, has only three distinct eigenvalues (−1/8, −1/32 and 3/2 after scaling). So I expect the K_n formula to fall to a computation, and the K_{m,m} one similarly. For hypercubes I do not see the pattern: 11/14, 31/45, 79/124, 191/315, … for d = 3, 4, 5, 6.

With non-uniform weights, here is the worst γ over twenty random weight vectors with every x_e ∈ [1/K, K], on the largest graph I tried in each family:

| Family | K = 1 | K = 2 | K = 10 | K = 100 |
|---|---|---|---|---|
| K₂₂ | 0.52 | 0.33 | 0.11 | 1e−2 |
| K₁₅,₁₅ | 0.55 | 0.36 | 0.10 | 3e−2 |
| Hypercube Q₇ | 0.59 | 0.35 | 0.06 | 8e−4 |
| Square grid 14 × 14 | 0.64 | 0.27 | 0.009 | 2e−5 |
| Triangular grid | 0.41 | 0.20 | 0.011 | 2e−5 |
| Cycle C₁₁₇ | 1 | 0.18 | 0.001 | 1e−7 |

For bounded K, γ shows no decay across a roughly thirty-fold range of graph sizes in any family. As K grows it collapses, and far faster on lattices and cycles than on dense graphs. That fits the mechanism: a huge weight behaves like contracting an edge, which creates parallel edges. Random sampling overestimates the true minimum, so these numbers are upper bounds on the worst case.

Against this stands a second observation. The minimum of γ(G, 1) over *all* connected graphs on V vertices keeps falling: 1, 0.625, 0.411, 0.320, 0.261 for V = 3, …, 7, computed exhaustively. Several graphs tie at V = 7. Structured families keep their margin while the worst graphs lose it.

**10. The margin questions.**

> *(D) Prove the closed forms for K_n and K_{m,m}, and find the one for hypercubes.*

> *(E) For weight ratios bounded by K, is γ ≥ c(K) > 0 uniformly in the size of the graph, for K_n, K_{m,m} or grids? If not, what are the bad weights?*

> *(F) Which graphs minimize γ(G, 1) on V vertices, and how fast does the minimum decay?*

Known results are qualitative (Murai–Nagaoka–Yazawa). Spectral independence (Anari–Liu–Oveis Gharan) controls the *largest* eigenvalue of the correlation matrix, the opposite end of the spectrum from γ. For (E), the natural tools seem to be perturbation from (D) via Weyl's inequalities for small K, and the negative correlation of spanning-tree measures for large K. Paragraph 7 shows these margins are not idle. For the Lorentzian property under a chemical potential, the relevant margin is the smallest eigenvalue, in absolute value, across *all* quadratic derivatives. It is exactly zero on K₄, on the diamond and on the grid alike; what separates them is the direction in which μ pushes the zero modes.

---

**11. Speculation.** The signature (1, n−1) is formally Lorentzian, and I began this work wondering whether it could be read as an emergent metric: a combinatorial seed of spacetime. I no longer think that reading can be taken for granted. If it is to mean anything physically, it must survive two things. It must survive a large-graph limit, which is what (E) is about. And it must survive finite density, which in general it does not, as paragraph 6 shows. Whether it survives on the graphs that matter, if any do, is question (A). I find it more honest, and more interesting, to let the mathematics decide.

**12. What I intend to do.** I plan to take question (A) to the Caltech Mathathon, with (B) as the route to a proof and (C) as the physical target. I would work on four threads at once:

- *Proof.* Carry out the perturbation analysis of paragraph 7 for a triangle's zero mode, to find exactly when it turns positive (the diamond shows the answer depends on the root). Then either prove the triangle-free case or find the graph that breaks it.
- *Search.* Encode a graph together with its orientation and root as a token sequence, as in Axplorer's built-in square-free-graph environment. Score it by the largest second eigenvalue over all quadratic derivatives, which is positive exactly on a violation. Axplorer alternates a transformer trained on the best examples with classical local search. It can hunt for triangle-free violations beyond the sizes where exhaustive search is possible. Every hit is re-verified exactly by the independent script. The number of quadratic derivatives grows quickly, so this is realistic only for graphs of moderate size. Scaling the scoring function is part of the work. After that come the time-layered lattices of (C).
- *Formalization.* Write a Lean 4 certificate of the K₄ counterexample. Taking e^μ = 2 makes every coefficient rational, and the violation persists there, so the certificate is a finite, exact computation: exhibit a two-dimensional subspace on which one explicit quadratic form is positive definite. When I last checked, Mathlib had matroids but no graphic matroid of a graph and no basis generating polynomial, so the general statements must wait. The K_n spectrum of (D) is the next target.
- *Writing.* A write-up that keeps the proved, the computed and the conjectured visibly apart, and an explanation I can defend step by step.

If the event goes well, the natural sequel is a paper on (A)–(C), with (D)–(E) as supporting results, together with a Lean development and the search data.

**13. What would change my mind.** A triangle-free graph that fails would end the simple form of (A). The perturbation analysis would then have to say which longer cycles matter. A paper answering (A) would mean the right move is to build on it. And if faithful lattices (C) always keep the Lorentzian property while toy orientations lose it, then the failure in paragraph 6 is an artifact of orienting every edge, and the physically interesting case is the robust one. I would count any of these as progress.

Most of this may be naive, and some of it may be known to people I have not read. I would be grateful to be told either.

With best regards,

Calvin Sabastian Tanzil

---

**P.S. on what came before, and on tools.** Everything described above was done before the Mathathon, and none of the questions (A)–(F) is answered here. Every number in this note can be reproduced from `verification/`:

- `hyp_check.py` for paragraph 2;
- `spectral_gap.py`, `gap_normalized.py` and `observations.py` for paragraph 9;
- `chemical_potential.py` for paragraphs 6 and 7.

I used AI assistants for drafting, for algebra, for literature searches and for writing code, and Axplorer is planned for the search in paragraph 12. The computations have an independent check: the closed-form scripts are validated against symbolic computation.

**References.**
- P. Brändén, J. Huh, *Lorentzian polynomials*, Ann. of Math. 192 (2020), 821–891.
- T. Nagaoka, A. Yazawa, *Strict log-concavity of the Kirchhoff polynomial and its applications to the strong Lefschetz property*, J. Algebra 577 (2021), 175–202.
- S. Murai, T. Nagaoka, A. Yazawa, *Strictness of the log-concavity of generating polynomials of matroids*, J. Combin. Theory Ser. A 181 (2021), 105351.
- S. Caracciolo, J. L. Jacobsen, H. Saleur, A. D. Sokal, A. Sportiello, *Fermionic field theory for trees and forests*, Phys. Rev. Lett. 93 (2004), 080601. [arXiv:cond-mat/0403271](https://arxiv.org/abs/cond-mat/0403271)
- M. Troyer, U.-J. Wiese, *Computational complexity and fundamental limitations to fermionic quantum Monte Carlo simulations*, Phys. Rev. Lett. 94 (2005), 170201. [arXiv:cond-mat/0408370](https://arxiv.org/abs/cond-mat/0408370)
- R. Burton, R. Pemantle, *Local characteristics, entropy and limit theorems for spanning trees and domino tilings via transfer-impedances*, Ann. Probab. 21 (1993), 1329–1371.
- N. Anari, K. Liu, S. Oveis Gharan, C. Vinzant, *Log-concave polynomials II*, Ann. of Math. 199 (2024), 259–299.
- N. Anari, K. Liu, S. Oveis Gharan, *Spectral independence in high-dimensional expanders and applications to the hardcore model*, FOCS 2020.

Full citations, with DOIs, are in the README.
