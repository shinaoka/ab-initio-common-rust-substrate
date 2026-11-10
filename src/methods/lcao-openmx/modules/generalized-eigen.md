# Generalized eigenvalue problem and regularization of S

> Status: prototype (awaiting review). The solution flow and the regularization rules were verified as behavior of the existing implementation,
> and the numerical properties were checked with `tools/check_kpoints_eigen.py`.
> For the matrix types (real/complex, dimension \\( N \\)/\\( 2N \\)), see [k points and Bloch sums](kpoints-bloch.md).

## Problem [Public]

As in Eqs. 27–29 of technical note 1, at each k point (only the Γ point for Cluster) and for each spin, the generalized eigenvalue problem

$$
H(\mathbf k)\\, c\_\mu = \varepsilon\_\mu\\, S(\mathbf k)\\, c\_\mu,\qquad c\_\mu^\dagger S(\mathbf k)\\, c\_\nu = \delta\_{\mu\nu}
$$

is solved. \\( H \\) is Hermitian (real symmetric for collinear Cluster), and \\( S \\) is Hermitian positive definite.
However, when the basis functions become nearly linearly dependent, the smallest eigenvalue of \\( S \\) approaches 0 and may even become negative due to rounding errors.
Below, the matrix dimension is denoted \\( n \\) (\\( N \\) for collinear, \\( 2N \\) for noncollinear).

## Solution flow [Convention]

The main paths of the existing implementation (Cluster, Band, the local problems of DC) orthogonalize using the eigendecomposition of \\( S \\), not the Cholesky decomposition.

1. Diagonalize \\( S \\): \\( S = U\Lambda U^\dagger \\), \\( \Lambda = \mathrm{diag}(\lambda\_1 \le \dots \le \lambda\_m) \\).
2. Build \\( \tilde\Lambda \\) from \\( \Lambda \\) by the regularization rule (next section). For rules that drop directions, the corresponding columns of \\( U \\) are also removed.
3. Build the transformation matrix \\( X = U\tilde\Lambda^{-1/2} \\).
4. Build \\( H' = X^\dagger H X \\), and solve the standard eigenvalue problem \\( H' d\_\mu = \varepsilon\_\mu d\_\mu \\) for only the \\( n\_{\mathrm{ev}} \\) lowest eigenpairs.
5. Transform back as \\( c\_\mu = X d\_\mu \\).

Then

$$
C^\dagger S\_{\mathrm{eff}}\\, C = I,\qquad H C = S\_{\mathrm{eff}}\\, C E,\qquad S\_{\mathrm{eff}} = U\tilde\Lambda U^\dagger
$$

hold (\\( C = (c\_1, \dots, c\_{n\_{\mathrm{ev}}}) \\), \\( E = \mathrm{diag}(\varepsilon\_\mu) \\)).
If regularization does not change the eigenvalues, \\( S\_{\mathrm{eff}} = S \\) and the result is exactly the solution of the original problem.

- **Collinear spin-polarized**: \\( S \\) does not depend on spin, so \\( X \\) is built once and used for both spins.
- **Noncollinear**: since \\( S^{(2)} = I\_2\otimes S \\), only the \\( N\times N \\) \\( S \\) is diagonalized, and \\( X^{(2)} = I\_2\otimes X \\).
  \\( H' \\) is the same as transforming each of the three spin blocks as \\( X^\dagger H^{\sigma\sigma'} X \\) and then
  assembling them into a \\( 2N\times 2N \\) matrix (\\( H'^{\downarrow\uparrow} = (H'^{\uparrow\downarrow})^\dagger \\)).
  In noncollinear Cluster calculations, \\( S \\) and \\( X \\) are real and \\( H' \\) is complex.
- **Reuse of \\( X \\)**: since \\( X \\) is determined by \\( S \\) alone, the existing implementation reuses \\( X \\) across SCF iterations
  while the atomic configuration does not change (Cluster, and some parallel configurations of Band).
- How the standard eigenvalue problem is solved (tridiagonalization of the dense matrix and the tridiagonal eigenvalue problem, distributed parallelism, GPU) is an implementation choice that does not affect the numerical results.
  [Public] Technical note 10 explains how to transform a Hermitian matrix into a real symmetric tridiagonal matrix by Householder transformations.
  [Convention] The library can be selected with `scf.eigen.lib` (`lapack`, `elpa1`, `elpa2`, `gpusolver`, `gpusolver2`, etc.),
  and the tridiagonal eigensolver with `scf.lapack.dste` (`dstevx` (default), `dstegr`, `dstedc`, `dsteqr`).
  Note, however, that the choice of `scf.eigen.lib` can change the regularization rule of the next section.

### Relation to Cholesky reduction [Public] (general linear algebra)

Compared with the reduction \\( H'' = L^{-1}HL^{-\dagger} \\) by \\( S = LL^\dagger \\), we have \\( X = L^{-\dagger}Q \\)
(\\( Q = L^\dagger U\Lambda^{-1/2} \\) is unitary), so \\( H' = Q^\dagger H'' Q \\).
As long as regularization does not act, the eigenvalues are the same, and the eigenvectors \\( c \\) are the same up to unitary transformations within degenerate subspaces (and the phase of each vector).
When regularization acts, a problem different from the Cholesky reduction is solved (the problem for \\( S\_{\mathrm{eff}} \\), or a subspace problem).

## Regularization rules [Convention]

The rule for building \\( \tilde\Lambda \\) from \\( \Lambda \\) is one of the following, depending on the path.

| Rule | Content | Dimension | \\( S\_{\mathrm{eff}} \\) |
|---|---|---|---|
| A | Replace only eigenvalues with \\( \lambda < 0 \\) by \\( 10^{-10} \\) | Unchanged | \\( \lVert S\_{\mathrm{eff}} - S\rVert \le 10^{-10} + \lvert\lambda\_{\min}\rvert \\) |
| A′ | Replace only eigenvalues with \\( \lambda < 0 \\) by \\( 10^{-14} \\) | Unchanged | Same as above (\\( 10^{-14} \\)) |
| B | Replace all eigenvalues with \\( \lambda < 10^{-10} \\) by \\( 10^{-10} \\) | Unchanged | \\( \lVert S\_{\mathrm{eff}} - S\rVert \le 10^{-10} + \lvert\lambda\_{\min}\rvert \\) |
| C | Remove eigenvalues with \\( \lambda < 10^{-9} \\) and their eigenvectors | Reduced from \\( n \\) by the number removed | Problem restricted to the subspace |
| D | \\( \tilde\lambda = \lvert\lambda\rvert \\) (no threshold) | Unchanged | \\( U\lvert\Lambda\rvert U^\dagger \\) |

The correspondence with paths is as follows.

| Calculation | Path | Rule |
|---|---|---|
| Cluster, collinear | CPU (`scf.eigen.lib` is `lapack`/`elpa1`/`elpa2`), and the multi-GPU distributed solver (`gpusolver2`) | A |
| Cluster, collinear | Dense solver on GPU (`gpusolver`, when the matrix satisfies the GPU-path conditions) | B |
| Cluster, noncollinear | All paths (CPU and GPU) | B |
| Band, collinear | All paths (CPU and GPU) | A |
| Band, noncollinear | All paths (CPU and GPU) | B |
| DC, collinear/noncollinear | Local problem per atom | C |
| Cluster-LNO | — | A |
| Writing out the density matrix for a specified energy range (post-processing) | Collinear/noncollinear | A/B |
| DC-LNO | Part of the local problems | D |
| Post-processing: band dispersion, density of states | — | A′ ([analysis](analysis.md); B for part of the molecular-orbital output) |
| Post-processing: exchange coupling constants (jx), polarization by Berry phase (polB) | — | No regularization ([analysis](analysis.md)) |

- The thresholds \\( 10^{-10} \\), \\( 10^{-14} \\), and \\( 10^{-9} \\) are all fixed values and **cannot be changed by input keywords**.
- In rule C, the eigenvalues are sorted in ascending order and all eigenvalues below the threshold are removed.
  The existing implementation also removes directions in the Krylov method using the same threshold \\( 10^{-9} \\) on the eigenvalues of the subspace overlap.
- If the "GPU-path conditions" (lower bound on the matrix dimension, whether it fits in memory, etc.) are not satisfied, the CPU path is taken even if `gpusolver` is specified,
  so for Cluster, collinear, A and B switch depending on the system size even for the same input.
- [Public] The difference among the above rules (A, B, C) is also pointed out in a public issue of the GPU version (dc1394/openmx4.0_gpu #1).

### Effect of rule differences on results [Design] (checked)

The "nearly singular \\( S \\)" example in `tools/check_kpoints_eigen.py` (\\( n = 8 \\), two pairs of basis functions nearly linearly dependent,
\\( \lambda\_1 \approx 3\times 10^{-16} \\), \\( \lambda\_2 \approx 7\times 10^{-14} \\), \\( \lambda\_3 \approx 2\times 10^{-2} \\)) shows the following.

- **Rule A** keeps eigenvalues with \\( 0 \le \lambda < 10^{-10} \\) as they are, so the components of \\( X \\) grow up to \\( \lambda^{-1/2} \sim 10^{8} \\),
  and the rounding errors of \\( H' \\) are amplified to about \\( \epsilon\_{\mathrm{mach}}\lVert H\rVert/\lambda \\). In this example \\( \kappa(S\_{\mathrm{eff}}) \approx 3\times 10^{15} \\),
  \\( \lVert C^\dagger S\_{\mathrm{eff}}C - I\rVert \approx 10^{-2} \\), and the lowest eigenvalue deviated by 5% from the values of rules B and C.
  When there are small positive eigenvalues, rule A may produce spurious eigenvalues.
- **Rule B** bounds \\( \kappa(S\_{\mathrm{eff}}) \le \lambda\_{\max}/10^{-10} \\), so the residual and orthogonality with respect to \\( S\_{\mathrm{eff}} \\) stay
  of order \\( \epsilon\_{\mathrm{mach}}\\,\kappa \\), and the physical (low) eigenvalues agreed with rule C to a relative difference of about \\( 10^{-5} \\).
- **Rule C** is the Rayleigh–Ritz method in a subspace, so the residual is 0 in the sense \\( X^\dagger(HC - SCE) = 0 \\), and \\( C^\dagger S C = I \\) holds to rounding error.
  However, the number of eigenpairs becomes smaller than \\( n \\).

## Number of eigenpairs to compute [Convention]

The number of eigenpairs \\( n\_{\mathrm{ev}} \\) computed at each SCF iteration is determined from the number of valence electrons
\\( N\_{\mathrm{el}} = \sum\_I Z\_I^{\mathrm{val}} - Q \\) (\\( Q \\) is `scf.system.charge`) as follows.
\\( \lfloor\cdot\rfloor \\) denotes the floor. In all cases \\( n \\) is the upper bound.

| Calculation | \\( n\_{\mathrm{ev}} \\) |
|---|---|
| Cluster, collinear | \\( n \\) at the first SCF iteration. Thereafter \\( \lfloor N\_{\mathrm{el}}/2\rfloor + \max(400, \lfloor 0.1\\,n\rfloor) \\) |
| Cluster, noncollinear | \\( n \\) at the first SCF iteration. Thereafter \\( n\_{\mathrm{HOMO}} + \max(400, \lfloor 0.2\\,n\rfloor) \\). \\( n\_{\mathrm{HOMO}} \\) is the index of the last level whose occupation exceeded \\( 1/2 \\) in the previous iteration |
| Band, collinear | \\( \lfloor N\_{\mathrm{el}}/2 + \lfloor x\rfloor\rfloor \\), \\( x = 0.2\\,N\_{\mathrm{el}} \\) (but \\( x = 100 \\) if \\( x < 50 \\)) |
| Band, noncollinear | \\( \lfloor N\_{\mathrm{el}} + \lfloor x\rfloor\rfloor \\), \\( x = 0.2\\,N\_{\mathrm{el}} \\) (but \\( x = 200 \\) if \\( x < 60 \\)) |
| Local problems of DC | All eigenpairs of the local problem (dimension reduced by rule C) |

- When Cluster outputs the density of states etc. (at the post-processing stage), the coefficients 0.1 (collinear) and 0.2 (noncollinear) become 0.2 and 0.4, respectively.
  When `partial.charge` is enabled, \\( n\_{\mathrm{ev}} = n \\).
- There is no input keyword that directly specifies \\( n\_{\mathrm{ev}} \\) (`num.HOMOs` and `num.LUMOs` are the numbers of molecular orbitals to output and do not affect this).
- Levels above the \\( n\_{\mathrm{ev}} \\)-th are treated as having occupation 0. Whether this assumption holds is not guaranteed at high electronic temperature or for metals.

[Open] This table gives the values of the version consulted (the GPU version). Whether the coefficients and lower bounds are the same in upstream OpenMX 3.9/4.0 has not been checked.

## Eigenvalue ordering, degeneracy, and return values [Convention]

- Eigenvalues are returned in ascending order (as output by the standard eigensolver). No explicit sorting is done.
- There is no special treatment of degeneracy. The eigenvectors within a degenerate subspace are left to the solver;
  since the occupations are determined by the eigenvalues alone (Fermi distribution), degenerate levels get the same occupation.
  As a result, the density matrix is invariant under unitary transformations within degenerate subspaces.
- For each k point and each spin, the ascending eigenvalues \\( \varepsilon\_1, \dots, \varepsilon\_{n\_{\mathrm{ev}}} \\) and
  the coefficients \\( C \\) (\\( n\times n\_{\mathrm{ev}} \\)) in the original basis satisfying \\( C^\dagger S\_{\mathrm{eff}} C = I \\) are returned.
  In the noncollinear case, the first \\( N \\) rows of the coefficients are spin up and the last \\( N \\) rows are spin down.

## Requirements for implementations [Design]

- **Input**: Hermitian matrices \\( H \\), \\( S \\) (real or complex; only the lower triangle is read), the number of eigenpairs to compute \\( n\_{\mathrm{ev}} \\) (or an upper bound on eigenvalues), and the regularization scheme.
- The **regularization scheme** is made explicit as an enum.
  - `Cholesky`: no regularization. Returns an error if \\( S \\) is not positive definite.
  - `ClampNegative { value }`: rules A and A′.
  - `Floor { threshold }`: rule B.
  - `Drop { threshold }`: rule C.
  - `Abs`: rule D (only for comparison with the existing implementation).
- The **default** is `Floor { threshold: 1e-10 }`. The reasons are that rule A can amplify errors for small positive eigenvalues (the check above),
  and that rule C changes the number of eigenpairs, complicating handling on the user side.
  For comparison with the existing implementation, the scheme is specified explicitly to match the target path.
- **Reporting**: return \\( \lambda\_{\min}(S) \\), \\( \lambda\_{\max}(S) \\), the number of eigenvalues replaced (or removed), \\( \lVert S\_{\mathrm{eff}} - S\rVert \\),
  and the resulting dimension. The caller must be able to know that a replacement occurred.
- **Noncollinear**: exploit the structure \\( S^{(2)} = I\_2\otimes S \\) and use only the \\( N\times N \\) decomposition.
- **Reuse**: return the decomposition of \\( S \\) (\\( U \\), \\( \tilde\Lambda \\), \\( X \\)) as a value so that it can be passed to the next call as long as \\( S \\) does not change.
- **Checks**: verify that the resulting eigenvalues and eigenvectors are finite and that the eigenvalues are in ascending order; return an error otherwise.
- **Consideration for later stages**: when the occupation of the \\( n\_{\mathrm{ev}} \\)-th level is not negligible (e.g. \\( f(\varepsilon\_{n\_{\mathrm{ev}}}) > 10^{-12} \\)),
  or when the \\( n\_{\mathrm{ev}} \\) boundary falls in the middle of degenerate levels, it must be possible to issue a warning or to increase \\( n\_{\mathrm{ev}} \\) and solve again.

## Tests

Residual and orthogonality are judged by the following quantities.

$$
\eta\_\mu = \frac{\lVert H c\_\mu - \varepsilon\_\mu S c\_\mu\rVert\_2}{\left(\lVert H\rVert\_2 + \lvert\varepsilon\_\mu\rvert\\,\lVert S\rVert\_2\right)\lVert c\_\mu\rVert\_2},\qquad
\omega = \lVert C^\dagger S C - I\rVert\_{\max}
$$

When regularization acts, \\( S \\) is replaced by \\( S\_{\mathrm{eff}} \\) in the judgment, and \\( \lVert S\_{\mathrm{eff}} - S\rVert \\) is recorded separately.
For rule C, instead of \\( \eta \\), the residual projected onto the subspace, \\( \lVert X^\dagger(HC - SCE)\rVert/\lVert H\rVert \\), is used.

| Kind | Content | Approximate tolerance |
|---|---|---|
| Invariant | \\( \max\_\mu \eta\_\mu \\) for all schemes with a well-conditioned \\( S \\) (\\( \kappa(S) \lesssim 10^{4} \\)) | \\( 10^{-13} \\) |
| Invariant | \\( \omega \\) with a well-conditioned \\( S \\) | \\( 10^{-12} \\) |
| Comparison | With a well-conditioned \\( S \\), compare eigenvalues with LAPACK `dsygvd` (real symmetric) and `zhegvd` (complex Hermitian) (all schemes) | \\( 10^{-12}\lVert H\rVert \\) |
| Comparison | Likewise, compare eigenvectors via the projector onto the degenerate subspace \\( C\_{\mathcal D} C\_{\mathcal D}^\dagger S \\) | \\( 10^{-10} \\) |
| Invariant | Noncollinear: the blockwise transformation by \\( I\_2\otimes X \\) gives the same eigenvalues as the \\( 2N \\)-dimensional generalized problem (`zhegvd`) | \\( 10^{-11} \\) |
| Invariant | Noncollinear with zero spin-off-diagonal blocks: the eigenvalues are the union of the eigenvalues of the two collinear problems | \\( 10^{-12} \\) |
| Invariant | Collinear spin-polarized: sharing \\( X \\) between both spins agrees with solving each spin independently | Rounding error |
| Invariant | Eigenvalues in ascending order; eigenvalues and eigenvectors finite | Exact |
| Invariant | Eigenpairs computed with \\( n\_{\mathrm{ev}} < n \\) agree with the first \\( n\_{\mathrm{ev}} \\) of all eigenpairs | \\( 10^{-12} \\) |
| Nearly singular | With \\( S \\) whose smallest eigenvalue is \\( 10^{-16} \\)–\\( 10^{-11} \\) (containing pairs of nearly linearly dependent basis functions), \\( \eta \\) of rule B (with respect to \\( S\_{\mathrm{eff}} \\)) | \\( 10^{-12} \\) |
| Nearly singular | Likewise, \\( \omega \\) of rule B (with respect to \\( S\_{\mathrm{eff}} \\)) is at most \\( 10^{3}\epsilon\_{\mathrm{mach}}\\,\kappa(S\_{\mathrm{eff}}) \\) | — |
| Nearly singular | The low eigenvalues (levels unrelated to the nearly linearly dependent directions) of rules B and C agree | \\( 10^{-4} \\) (relative) |
| Nearly singular | The dimension of rule C equals \\( \#\\{\lambda \ge 10^{-9}\\} \\), with \\( C^\dagger S C = I \\) and zero projected residual | \\( 10^{-12} \\) |
| Nearly singular | With \\( S \\) having a negative eigenvalue (around \\( -10^{-12} \\)), rules A and B return finite results and `Cholesky` returns an error | — |
| Nearly singular | \\( S \\) with an eigenvalue exactly 0: rule A does not replace \\( \lambda = 0 \\), so it diverges. Decide whether the reimplementation's `ClampNegative` replaces \\( \lambda \le 0 \\) or returns an error, and test it | — |
| Comparison | With \\( H \\), \\( S \\) written out from the existing implementation (Cluster collinear and noncollinear, Band collinear and noncollinear at general k points), eigenvalues when the same rule is specified | \\( 10^{-10} \\) Hartree |
| Comparison | With the same input, the subspace of occupied levels (density matrix \\( C\_{\mathrm{occ}}C\_{\mathrm{occ}}^\dagger \\)) agrees | \\( 10^{-8} \\) |

## Open issues [Open]

- Details of the treatment of \\( S \\) in the `EGAC` and `Cluster-LNO` paths (the table above covers what could be verified). For the local problems of `Krylov`, `Cluster2` (formerly `ON2`), and `DC-LNO`, see [O(N) methods](on-methods.md); for the treatment of \\( S \\) in NEGF, see [electron transport by NEGF](negf.md).
- The concrete content of the "GPU-path conditions" (dimension threshold, memory check). It is safe to build comparison fixtures from paths whose rule is uniquely determined (Band, noncollinear).
- Differences from upstream OpenMX: the GPU version was consulted; whether Cluster, collinear in the upstream CPU version uses only rule A, and whether the coefficients of \\( n\_{\mathrm{ev}} \\) are the same, has not been checked.
- Whether to include in the comparison the behavior of the existing implementation when \\( \lambda = 0 \\) exactly under rule A (it diverges).
