I want to set down some observations about the spanning-tree polynomial of a graph, together with the questions they have led me to. Some of what follows is proved, some is only computed, and some is speculation. I have tried to say which is which each time. Paragraphs 1–3 recall what is known, paragraph 4 explains why I care, paragraphs 5–10 contain what I think is new, paragraph 11 is the speculative one, and paragraph 14 opens a second direction, closer to the sign problem itself. If the speculative part does not interest you, the rest should stand without it.

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

I should be plain about the limits. At zero chemical potential these fermions have no sign problem at all. Under the driven deformation of paragraphs 5–8 the weights stay positive, the partition function is a determinant, and weighted arborescences can still be sampled exactly by Wilson's algorithm. Under a genuine chemical potential some weights do turn negative (paragraph 14), but the partition function is still a single determinant, so this class is a laboratory, not a hard case. The physically hard cases, such as doped Hubbard models or QCD at finite baryon density, involve fermion matrices of a different kind, and nothing here reaches them. My interest is structural: *which deformations keep positivity together with log-concavity, and which keep positivity but lose log-concavity?* (The reverse cannot happen here: a Lorentzian polynomial has nonnegative coefficients by definition.)

**5. The deformation: driven hopping.** The deformations that matter physically make hopping direction-dependent, weighting forward and backward hops by e^{+μ} and e^{−μ}. (For a symmetric rescaling x_e ↦ w_e x_e the question is empty: that is a nonnegative linear change of variables, which preserves the Lorentzian property by Brändén–Huh.) There are two natural ways to put such hopping into a Laplacian, and they differ only in the diagonal. A lattice chemical potential, in the sense of Hasenfratz and Karsch, leaves the diagonal alone; I come to it in paragraph 14. The other way adjusts the diagonal so that every column of the matrix sums to zero. Up to sign and transposition, that matrix is the generator of a biased random walk, the kind of operator that describes driven hopping out of equilibrium. I call it the *driven deformation*, and paragraphs 5–8 are about it. An earlier version of this note called it a chemical potential. It is not one, and the difference changes the answer completely.

For the driven deformation, Tutte's directed matrix-tree theorem writes the reduced determinant as a positive sum over arborescences. Directing each spanning tree away from a root r,

$$f_\mu(x) = \sum_T e^{\mu\, s_r(T)} \prod_{e \in T} x_e,$$

where s_r(T) is the number of forward edges minus the number of backward edges. The support is still the set of spanning trees, so it is still M-convex. By Brändén–Huh, f_μ is Lorentzian *if and only if every quadratic derivative* ∂^S f_μ (with |S| = deg − 2) has a Hessian with at most one positive eigenvalue. This makes the question exact and finite for each graph. The quadratic derivative ∂^S f_μ is the tree polynomial of the contraction G/S, a graph on three vertices, with its coefficients inherited from the orientations.

For my first experiments, "forward" means from lower to higher vertex label, on every edge. This is a toy; paragraph 8 returns to time-layered lattices.

**6. What happens.** I expected the Lorentzian property to survive. It does not, at least not always. Here are the results of the exact test for μ from 0.01 to 5, and for K₄ also for μ = 0.001:

| Graph | Root 0 | Root 1 | Triangles |
|---|---|---|---|
| K₄ | **fails** at every μ tested, down to 0.001 | **fails** | 4 |
| K₅ | **fails** | **fails** | 10 |
| Diamond (K₄ minus an edge) | survives | **fails** | 2 |
| C₄ | survives | survives | 0 |
| C₅ | survives | survives | 0 |
| 2×3 grid | survives | survives | 0 |
| K₂,₃ | survives, trivially | survives, trivially | 0 |
| House | **fails** | **fails** | 1 |
| Bowtie (two triangles at a vertex) | survives | survives | 2 |
| C₆, K₃,₃, 3×3 grid | survive | survive | 0 |

On K₄ the failure is not confined to a quadratic derivative. At μ = 1, at 62 of 2,000 random positive points, the full Hessian of f_μ has *two* positive eigenvalues, so the Lorentzian signature itself is lost.

K₂,₃ survives for a reason worth recording, because it shows what a trivial survival looks like. Its labels put one side {0, 1} entirely below the other side {2, 3, 4}, so every edge points from the first side to the second. In a tree directed away from a root on the first side, each vertex on the second side is entered forward, and each non-root vertex on the first side is entered backward. So s_r(T) is the same for every tree, and f_μ is a constant multiple of f_0. The same argument works for a root on the other side, and for any bipartite graph whose labels respect the bipartition in this way.

**7. The mechanism.** Three observations explain the table. The first two are short proofs; the third is a theorem of Brändén and Huh that turns the whole question into combinatorics.

*(a) Where violations can come from (proved).* Each quadratic derivative is a symmetric matrix depending continuously on μ. If at μ = 0 it has one positive eigenvalue and no zero eigenvalue, it keeps that signature for all small μ. So f_μ is Lorentzian for all sufficiently small μ unless some quadratic derivative of f_0 has a zero eigenvalue. Call these the *zero modes*. The cycles have none (0 of 4 for C₄, 0 of 10 for C₅), so their survival at small μ is automatic. Their survival up to μ = 5 is observed.

*(b) Where the zero modes come from (proved).* At μ = 0 every coefficient is 1, so the Hessian of a quadratic derivative has entry 1 for two edges of G/S that join different pairs of its three vertices, and 0 otherwise. Edges joining the same pair give identical rows, so the rank is at most 3, and there is a zero eigenvalue exactly when G/S has a pair of parallel edges. This is the degeneracy of paragraph 2 again. (It agrees with direct computation on all 204 quadratic derivatives of eight small graphs.) Contracting one edge of a triangle produces a parallel pair; on the grid, contracting two edges of a 4-cycle does. On K₄, every one of the six quadratic derivatives has a zero mode. Take the one obtained by contracting edge {2, 3}:

| μ | Eigenvalues of that quadratic derivative |
|---|---|
| 0 | −2, −1.236, 0, 0, 3.236 |
| 0.05 | −2.005, −1.431, 0, **+0.0025**, 3.434 |
| 1 | −27.21, −3.81, 0, **+0.72**, 30.30 |

A nonzero μ weights the parallel pair differently and pushes one zero eigenvalue positive, at order μ²: 0.0025 = 0.05².

So zero modes are *necessary* for failure at small μ, by (a), but not *sufficient*. The 2×3 grid has zero modes in 27 of its 35 quadratic derivatives, and it survives at every μ tested. What decides is the *direction* in which μ pushes each zero mode. In the language of degenerate perturbation theory: if Q(μ) = Q₀ + μQ₁ + μ²Q₂ + … and P projects onto the kernel of Q₀, the first-order effective matrix is P Q₁ P. When that vanishes, the second-order one is P(Q₂ − Q₁ Q₀⁺ Q₁)P, where Q₀⁺ is the pseudoinverse. The observed μ² scaling on K₄ suggests the first-order term vanishes there. For small μ, a zero mode turns positive when the first non-vanishing effective matrix has a positive eigenvalue. I think this is the right way to attack the problem, but I have not carried it out.

*(c) What decides it for every μ (Brändén–Huh, Theorem 3.14).* Brändén and Huh prove that Σ_α q^{ν(α)} x^α/α! is Lorentzian for all 0 < q ≤ 1 if and only if ν is M-convex. Here the polynomial is multiaffine, q = e^{−μ} and ν = −s_r. So f_μ is Lorentzian for every μ ≥ 0 exactly when T ↦ s_r(T) is M-concave on the spanning trees, that is, when it defines a valuated matroid: for all trees A, B and every a ∈ A∖B there is b ∈ B∖A such that A − a + b and B − b + a are trees and

$$s_r(A) + s_r(B) \le s_r(A-a+b) + s_r(B-b+a).$$

The question is then no longer about eigenvalues at all. As a check, the exchange test and the eigenvalue test agree in all 24 cases I tried: twelve graphs, two roots each, including the house, which fails, and the bowtie, which has triangles and survives for both roots. I suspect the zero-mode mechanism of (b) is the local, two-swap face of this exchange condition, but I have not proved it.

**8. The questions.**

> *(A) For which graphs, roots and orientations is s_r M-concave on the spanning trees, so that f_μ stays Lorentzian for every μ? Is it enough for the graph to be triangle-free?*

Every failure I have seen is on a graph with triangles, and every triangle-free graph I have tried survives. But triangles do not force failure: the diamond survives for one root and fails for the other, and the bowtie survives for both.

> *(B) Is there a clean criterion for the sign of the effective perturbation in paragraph 7 on a triangle's zero mode, and on the zero modes that come from longer cycles?*

Answering (B) for triangles would explain the diamond, where the same triangles lead to failure for one root and not for the other. The grid, whose zero modes come from 4-cycles and never turned positive in my tests, suggests longer cycles behave differently.

> *(C) On time-layered lattices, with e^{±μ} only on the time-like edges, when does the driven deformation stay Lorentzian?*

Here G = C_L × H, with time periodic. Under the genuine chemical potential of paragraph 14, μ enters only through cycles that wind around periodic time, as it does in lattice field theory. I would like to know whether the failure of the Lorentzian property under the driven deformation tracks winding in the same way. Net direction alone cannot be the whole story: under my toy orientation, C₄ and C₅ have nonzero net direction around the cycle, and both survive.

By paragraph 7(c), (A) is equivalent to a question about valuated matroids; that equivalence is Brändén and Huh's, not mine. I have not found the combinatorial question itself, or (B) and (C), treated in the literature, but I have not yet searched as carefully as I must before calling them open.

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

Known results are qualitative (Murai–Nagaoka–Yazawa). Spectral independence (Anari–Liu–Oveis Gharan) controls the *largest* eigenvalue of the correlation matrix, the opposite end of the spectrum from γ. For (E), the natural tools seem to be perturbation from (D) via Weyl's inequalities for small K, and the negative correlation of spanning-tree measures for large K. Paragraph 7 shows these margins are not idle. For the Lorentzian property under the driven deformation, the relevant margin is the smallest eigenvalue, in absolute value, across *all* quadratic derivatives. It is exactly zero on K₄, on the diamond and on the grid alike; what separates them is the direction in which μ pushes the zero modes.

---

**11. Speculation.** The signature (1, n−1) is formally Lorentzian, and I began this work wondering whether it could be read as an emergent metric: a combinatorial seed of spacetime. I no longer think that reading can be taken for granted. First, the signature lives on the n-dimensional space of edge weights, not on a spacetime, and every simple graph has it, so by itself it selects no dimension and no geometry. Second, a model with positive weights is a Euclidean statistical model: Lorentzian physics is reached from one by analytic continuation, and a direct Lorentzian path integral carries oscillating, signed weights. Paragraph 14 shows the same pattern in miniature, with the positive, Lorentzian polynomial at imaginary chemical potential and signs appearing as one continues to real μ. So the sign-free structure this note studies and an emergent Lorentzian spacetime pull in opposite directions. Any spacetime reading would also have to survive a large-graph limit, which is what (E) is about. I find it more honest, and more interesting, to let the mathematics decide.

**12. What I intend to do.** I plan to take question (A) to the Caltech Mathathon, with (B) as the route to a proof and (C) as the lattice case. I would work on four threads at once:

- *Proof.* Characterize the M-concavity of s_r directly from the exchange condition of paragraph 7(c), using the perturbation analysis of paragraph 7 to see which two-tree swaps fail and why (the diamond shows the answer depends on the root). Then either prove the triangle-free case or find the graph that breaks it.
- *Search.* Encode a graph together with its orientation and root as a token sequence, as in Axplorer's built-in square-free-graph environment. Score it by the total violation of the exchange inequality of paragraph 7(c), which is zero exactly when f_μ stays Lorentzian for every μ. Axplorer alternates a transformer trained on the best examples with classical local search. It can hunt for triangle-free violations beyond the sizes where exhaustive search is possible. Every hit is re-verified exactly by the independent script. The number of pairs of trees grows quickly, so this is realistic only for graphs of moderate size. Scaling the scoring function is part of the work. After that come the time-layered lattices of (C).
- *Formalization.* Write a Lean 4 certificate of the K₄ counterexample. The exchange form of paragraph 7(c) makes this purely combinatorial: exhibit two spanning trees A, B and an edge a ∈ A∖B for which every exchange fails the inequality. For a specific μ there is also a linear-algebra certificate: taking e^μ = 2 makes every coefficient rational, and one explicit quadratic form is positive definite on a two-dimensional subspace. When I last checked, Mathlib had matroids but no graphic matroid of a graph and no basis generating polynomial, so the general statements must wait. The K_n spectrum of (D) is the next target.
- *Writing.* A write-up that keeps the proved, the computed and the conjectured visibly apart, and an explanation I can defend step by step.

If the event goes well, the natural sequel is a paper on (A)–(C), with (D)–(E) as supporting results, together with a Lean development and the search data.

**13. What would change my mind.** A triangle-free graph that fails would end the simple form of (A). The perturbation analysis would then have to say which longer cycles matter. A paper answering (A) would mean the right move is to build on it. And if time-layered lattices (C) always keep the Lorentzian property while toy orientations lose it, then the failure in paragraph 6 is an artifact of orienting every edge, and the interesting case is the robust one. I would count any of these as progress.

**14. A second direction: a genuine chemical potential.** A lattice chemical potential, in the sense of Hasenfratz and Karsch, puts e^{±μ} on the time-like hops and leaves the diagonal alone:

$$M(\mu) = D - W(\mu), \qquad W_{uv} = x_e\, e^{\mu t_{uv}},$$

where D is the ordinary degree matrix and t_{uv} is +1, −1 or 0 as the hop goes forward in time, backward, or sideways. Its determinant is no longer a sum over trees. By a theorem of Kenyon it is a sum over *cycle-rooted spanning forests*, subgraphs in which every component contains exactly one cycle:

$$\det M(\mu) = \sum_{F} \prod_{e \in F} x_e \prod_{C \subset F} \big(2 - e^{\mu \ell_C} - e^{-\mu \ell_C}\big),$$

the second product running over the cycles of F. Here ℓ_C is the net number of forward-minus-backward time steps around C, that is, L times the number of times C winds around periodic time. A cycle that does not wind contributes zero, so μ enters only through winding, as it does in lattice field theory. I checked the formula against the determinant on periodic ladders, two time-loops of length L joined by rungs, for real and imaginary μ alike.

The sign of the cycle factor decides everything.

- *Imaginary μ = iθ.* The factor is 4 sin²(θℓ_C/2) ≥ 0, so every weight is nonnegative. More is true. The matrix is then a sum of rank-one positive semidefinite matrices, one per edge, so by Cauchy–Binet each coefficient is |det V_F|² for an explicit matrix V_F, and by Borcea–Brändén the determinant is a real stable polynomial in the edge weights. Homogeneous real stable polynomials with nonnegative coefficients are Lorentzian (Brändén–Huh, Proposition 2.2), and the exact test of paragraph 5 confirms it on the ladders. For generic θ, the forests whose cycles all wind are the bases of a matroid, the frame matroid of the graph with gains e^{iθt} (Zaslavsky), so this is a Lorentzian basis polynomial of a matroid other than the graphic one.
- *Real μ.* The factor is −4 sinh²(μℓ_C/2) < 0, so a forest's weight has sign (−1)^{number of cycles}. On the ladder with L = 3, 50 of the 51 forests with nonzero weight are negative; for L = 4, 192 of 193. In the Cauchy–Binet form the coefficient becomes det V_F · det W_F for two *different* matrices, and positivity is gone.

So the imaginary side has the whole package, positive and log-concave, and the real side has signs. That is the pattern lattice QCD lives with, where one standard route is to simulate at imaginary μ and continue analytically.

Two cautions. This model is Gaussian: its partition function is a single determinant, so the signs exist only in the forest representation, and it is a laboratory, not a hard case. And on the lattices I could enumerate, the sign problem is mild. With a mass term m = 0.5, the average sign (the partition function divided by the sum of the absolute weights) grows with the time extent but falls with the spatial size: at L = 3 and μ = 0.6 it is 0.165, 0.132 and 0.098 for two, three and four spatial sites. That is the direction a sign problem goes, but three points cannot establish the exponential decay in volume that makes one severe.

> *(G) How does the real-stable, Lorentzian structure at imaginary μ break down along the continuation to real μ? At what rate in the spatial volume does the average sign decay, and how does the rate depend on μ?*

> *(H) The interacting case. The arboreal gas, spanning forests with the four-fermion term of Caracciolo et al., is the H^{0|2} sigma model. Bauerschmidt, Crawford, Helmuth and Swan used its OSp(1|2) symmetry to show that trees do not percolate in two dimensions, and Bauerschmidt, Crawford and Helmuth combined it with the renormalization group to prove a percolation transition in three or more dimensions. That model is not exactly solvable, so a chemical potential there would give a sign problem in earnest. Does the μ-deformation preserve OSp(1|2)? If it does, their machinery may apply. If it does not, which part of the positive, log-concave structure survives?*

Troyer and Wiese make it unlikely that a general positive re-expansion exists at real μ. What I would hope for is a special class where one does, and a precise account of why.

Most of this may be naive, and some of it may be known to people I have not read. I would be grateful to be told either.

With best regards,

Calvin Sabastian Tanzil

---

**P.S. on what came before, and on tools.** Everything described above was done before the Mathathon, and none of the questions (A)–(H) is answered here. The reduction in paragraph 7(c) is Brändén and Huh's theorem applied to this family; I checked it numerically but proved nothing new there. Likewise, paragraph 14 combines theorems of Kenyon, Borcea–Brändén, Brändén–Huh and Zaslavsky; what is mine there is the computation and the questions. Every number in this note can be reproduced from `verification/`:

- `hyp_check.py` for paragraph 2;
- `spectral_gap.py`, `gap_normalized.py` and `observations.py` for paragraph 9;
- `driven_deformation.py` for paragraphs 5–7;
- `chemical_potential.py` for paragraph 14.

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
- P. Hasenfratz, F. Karsch, *Chemical potential on the lattice*, Phys. Lett. B 125 (1983), 308–310.
- R. Kenyon, *Spanning forests and the vector bundle Laplacian*, Ann. Probab. 39 (2011). [arXiv:1001.4028](https://arxiv.org/abs/1001.4028)
- J. Borcea, P. Brändén, *Applications of stable polynomials to mixed determinants: Johnson's conjectures, unimodality, and symmetrized Fischer products*, Duke Math. J. 143 (2008), Proposition 2.4. [arXiv:math/0607755](https://arxiv.org/abs/math/0607755)
- T. Zaslavsky, *Biased graphs. II. The three matroids*, J. Combin. Theory Ser. B 51 (1991), 46–72.
- R. Bauerschmidt, N. Crawford, T. Helmuth, A. Swan, *Random spanning forests and hyperbolic symmetry*, Commun. Math. Phys. 381 (2021), 1223–1261. [arXiv:1912.04854](https://arxiv.org/abs/1912.04854)
- R. Bauerschmidt, N. Crawford, T. Helmuth, *Percolation transition for random forests in d ≥ 3*, Invent. Math. 237 (2024), 445–540. [arXiv:2107.01878](https://arxiv.org/abs/2107.01878)

Full citations, with DOIs, are in the README.
