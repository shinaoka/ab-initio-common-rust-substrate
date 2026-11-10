# k points and Bloch sums

> Status: prototype (awaiting review). The Bloch-sum formulas rewrite technical note 1 (Eqs. 1, 4, 28–29) in the notation of this note;
> the generation and weights of k points were verified as behavior of the existing implementation.
> The numerical properties were checked with `tools/check_kpoints_eigen.py`.

## Scope of this chapter

This chapter defines the conventions for building the matrices \\( H(\mathbf k) \\), \\( S(\mathbf k) \\) at each k point
from the real-space matrices \\( O(\mathbf L) \\) ([sparse block format](../design/sparse-block-format.md)), as well as the set of k points and their weights.
How the eigenvalue problem at each k point is solved is treated in [generalized eigenvalue problem and regularization of S](generalized-eigen.md),
and the details of building the density matrix from eigenvectors in [density matrix and charge mixing](density-mixing.md).

## Lattice and reciprocal lattice [Design]

The unit-cell vectors are \\( \mathbf a\_1, \mathbf a\_2, \mathbf a\_3 \\), and the reciprocal lattice vectors are defined by
\\( \mathbf b\_i\cdot\mathbf a\_j = 2\pi\delta\_{ij} \\).
A lattice vector is represented by an integer triple \\( \mathbf n = (n\_1, n\_2, n\_3) \\) as
\\( \mathbf L = \sum\_i n\_i \mathbf a\_i \\), and a k point by fractional coordinates \\( \mathbf q = (q\_1, q\_2, q\_3) \\) as
\\( \mathbf k = \sum\_i q\_i \mathbf b\_i \\). Then

$$
\mathbf k\cdot\mathbf L = 2\pi\\,(q\_1 n\_1 + q\_2 n\_2 + q\_3 n\_3)
$$

holds. Input, output, and internal representation of k points are based on the fractional coordinates \\( \mathbf q \\).

## Bloch sum [Public] [Convention]

[Public] Eq. 1 of technical note 1 expands the Kohn–Sham orbital with wave vector \\( \mathbf k \\) and band \\( \mu \\) as

$$
\psi\_{\mu\mathbf k}(\mathbf r) = \frac{1}{\sqrt{N\_{\mathrm{cell}}}}\sum\_{\mathbf L} e^{i\mathbf k\cdot\mathbf L}
\sum\_{I\alpha} c\_{I\alpha,\mu}(\mathbf k)\\, \phi\_{I\alpha}(\mathbf r - \mathbf R\_I - \mathbf L)
$$

(the normalization factor is written differently from the original, but is chosen to be consistent with the normalization \\( c^\dagger S c = 1 \\) of \\( c \\)).
The phase \\( e^{i\mathbf k\cdot\mathbf L} \\) contains **only the lattice vector**; the atomic position \\( \mathbf R\_I \\) does not enter.
Substituting this into the Hamiltonian (ibid. Eqs. 28–29) gives

$$
O\_{I\alpha, J\beta}(\mathbf k) = \sum\_{\mathbf L} e^{+i\mathbf k\cdot\mathbf L}\\, O\_{I\alpha, J\beta}(\mathbf L),\qquad
O\_{I\alpha, J\beta}(\mathbf L) = \langle \phi\_{I\alpha}(\mathbf r - \mathbf R\_I) | \hat O | \phi\_{J\beta}(\mathbf r - \mathbf R\_J - \mathbf L) \rangle
$$

where \\( \hat O \\) is the overlap, the kinetic energy, the Kohn–Sham Hamiltonian, etc.
\\( \mathbf L \\) points to the periodic image of the **column-side** (ket-side) atom \\( J \\), and the sign is \\( +i \\).

[Convention] The \\( (I, J) \\) block of \\( O(\mathbf k) \\) in the existing implementation is the sum of \\( e^{2\pi i\\,\mathbf q\cdot\mathbf n}\\, O\_{IJ}(\mathbf L) \\)
over the periodic images of \\( J \\) (positions \\( \mathbf R\_J + \mathbf L \\)) (multiple periodic images are summed).
That is, it follows the formula above, and atomic positions do not enter the phase.

### Properties [Public]

The following follow immediately from the definition and can be used as test invariants.

- Hermiticity: if \\( O\_{J\beta, I\alpha}(-\mathbf L) = \overline{O\_{I\alpha, J\beta}(\mathbf L)} \\) holds in real space, \\( O(\mathbf k) \\) is Hermitian.
- Periodicity: since the phase involves only lattice vectors, \\( O(\mathbf k + \mathbf G) = O(\mathbf k) \\) for a reciprocal lattice vector \\( \mathbf G \\).
  Shifting the fractional coordinates of a k point by 1 does not change the matrix.
- Time reversal: if \\( O(\mathbf L) \\) is a real matrix (\\( H \\) and \\( S \\) in collinear calculations), \\( O(-\mathbf k) = \overline{O(\mathbf k)} \\).
  Hence the eigenvalues satisfy \\( \varepsilon\_\mu(-\mathbf k) = \varepsilon\_\mu(\mathbf k) \\), and the eigenvectors can be chosen as \\( c(-\mathbf k) = \overline{c(\mathbf k)} \\).

### Relation to the convention including atomic positions in the phase [Public]

In the convention where the phase is \\( e^{i\mathbf k\cdot(\mathbf R\_J + \mathbf L - \mathbf R\_I)} \\) (the convention including atomic positions),
using the diagonal matrix \\( P(\mathbf k) = \mathrm{diag}\left(e^{i\mathbf k\cdot\mathbf R\_I}\right) \\) (the same phase for all orbitals of atom \\( I \\)),

$$
\tilde O(\mathbf k) = P(\mathbf k)^\dagger\\, O(\mathbf k)\\, P(\mathbf k),\qquad \tilde c(\mathbf k) = P(\mathbf k)^\dagger c(\mathbf k)
$$

The eigenvalues do not change. The eigenvector components change only by per-atom phases.
Since \\( \tilde O(\mathbf k + \mathbf G) \\) does not coincide with \\( \tilde O(\mathbf k) \\) (it is unitarily equivalent),
post-processing involving k derivatives (velocity operator, Berry phase, etc.) must always state explicitly which convention is used.
This note, like the existing implementation, takes the lattice-vector-only convention as the default.

### Real-space density matrix [Public] [Convention]

[Public] Corresponding to Eq. 4 of technical note 1, the real-space density matrix is built from the k-point eigenvectors as

$$
\rho\_{I\alpha, J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k} \sum\_\mu f\_{\mu\mathbf k}\\,
e^{+i\mathbf k\cdot\mathbf L}\\, \overline{c\_{I\alpha,\mu}(\mathbf k)}\\, c\_{J\beta,\mu}(\mathbf k)
$$

(\\( \sum\_{\mathbf k} w\_{\mathbf k} = 1 \\), \\( f \\) is the occupation). With this phase direction,
the band energy is \\( \sum\_{\mathbf L}\sum\_{I\alpha J\beta} \rho\_{I\alpha,J\beta}(\mathbf L)\\, H\_{I\alpha,J\beta}(\mathbf L) \\),
and the number of electrons is the same expression with \\( H \\) replaced by \\( S \\).

[Convention] The existing implementation for collinear calculations accumulates only the real part of the above expression. Since by time reversal the contributions of \\( \mathbf k \\) and \\( -\mathbf k \\) are complex conjugates of each other,
representing a \\( \pm\mathbf k \\) pair by a single point with weight 2 is exact once the real part is taken ("Reduction of k points" below).

## Calculation types and matrix types [Public] [Convention]

[Public] According to the manual (keyword section), the solver is selected with `scf.EigenvalueSolver`;
options include `Cluster` (isolated systems), `Band` (periodic systems, k-point sampling), `DC`, `DC-LNO`, `Krylov` (O(N) methods), the low-order scaling method, etc.
[Convention] Manual Section 8.2 writes the value for the low-order scaling method as `ON2`, but the value accepted by the existing implementation is `Cluster2` ([O(N) methods](on-methods.md)).
`scf.SpinPolarization` takes `Off` (non-spin-polarized), `On` (collinear spin-polarized), or `NC` (noncollinear),
and `scf.SpinOrbit.Coupling` can be used only in noncollinear calculations.

[Convention] The matrix types solved at each k point (or each local problem) in the existing implementation are as follows.
\\( N \\) is the total number of orbitals in the unit cell (or computational region).

| Solver | Spin | Matrices solved | Type | Count |
|---|---|---|---|---|
| Cluster | `Off` / `On` | \\( H(\Gamma) \\), \\( S(\Gamma) \\) | Real symmetric \\( N\times N \\) | One per spin (one for `Off`) |
| Cluster | `NC` | Spinor version of \\( H(\Gamma) \\), \\( S(\Gamma) \\) | Complex Hermitian \\( 2N\times 2N \\) | One |
| Band | `Off` / `On` | \\( H(\mathbf k) \\), \\( S(\mathbf k) \\) | Complex Hermitian \\( N\times N \\) | k points × spins |
| Band | `NC` | Spinor version | Complex Hermitian \\( 2N\times 2N \\) | One per k point |
| DC | `Off` / `On` | Local matrix per atom | Real symmetric | Atoms × spins |
| DC | `NC` | Spinor version of the local matrix per atom | Complex Hermitian | One per atom |

- \\( H(\Gamma) \\) for Cluster is \\( \sum\_{\mathbf L} O(\mathbf L) \\) (the Bloch sum at \\( \mathbf k = 0 \\)).
  If there are no nearby periodic images, as for an isolated molecule, only \\( \mathbf L = 0 \\) remains. When Cluster is specified for a periodic unit cell,
  the contributions of periodic images are simply summed, giving a Γ-point calculation.
- With `Band` and `scf.Kgrid` equal to `1 1 1`, the existing implementation switches the solver to Cluster (a Γ-point calculation with periodic boundary conditions) and solves with real arithmetic.
- In collinear spin-polarized calculations, the \\( H \\) for spin up and spin down are separate matrices, and \\( S \\) is common.
- In collinear calculations, all \\( O(\mathbf L) \\) are real. In noncollinear calculations, spin-off-diagonal blocks and
  imaginary parts due to spin-orbit coupling, LDA+U, spin constraints, the Zeeman term, etc. are added.

[Convention] The matrix types and handling of the O(N) solvers (`DC`, `Krylov`, `DC-LNO`, `Cluster2`) are defined in [O(N) methods](on-methods.md), and the matrices of `NEGF` (k-space blocks of the central region and the leads) in [electron transport by NEGF](negf.md).
[Open] The matrix type of `Cluster-LNO`.
The existing implementation also retains the values `Recursion`, `Cluster-DIIS`, `Cluster2`, and `EGAC`, which need separate checking together with the manual's description.

## Ordering of noncollinear spinors [Convention]

The \\( 2N \\)-dimensional matrices of noncollinear calculations are ordered with **spin outer and orbital inner**. The row and column index is

$$
p = \sigma N + a,\qquad \sigma \in \\{0\\,(\uparrow),\\, 1\\,(\downarrow)\\},\quad a = 0, \dots, N-1
$$

where \\( a \\) is the global orbital index, the same as in collinear calculations (atoms in input order, and within each atom the order of the [PAO basis](pao-basis.md)).
Therefore

$$
H^{(2)}(\mathbf k) = \begin{pmatrix} H^{\uparrow\uparrow}(\mathbf k) & H^{\uparrow\downarrow}(\mathbf k) \\\\ H^{\downarrow\uparrow}(\mathbf k) & H^{\downarrow\downarrow}(\mathbf k) \end{pmatrix},\qquad
S^{(2)}(\mathbf k) = \begin{pmatrix} S(\mathbf k) & 0 \\\\ 0 & S(\mathbf k) \end{pmatrix} = I\_2 \otimes S(\mathbf k)
$$

The spin quantization axis is the \\( z \\) axis of the global coordinates, and \\( H^{\uparrow\downarrow}\_{ab} = \langle a\uparrow | \hat H | b\downarrow\rangle \\).

- Each spin block is built by the same Bloch sum from its real-space matrix \\( H^{\sigma\sigma'}(\mathbf L) \\).
- The existing implementation builds the three blocks \\( \uparrow\uparrow \\), \\( \downarrow\downarrow \\), \\( \uparrow\downarrow \\) by Bloch sums,
  and fills the remaining one with \\( H^{\downarrow\uparrow}(\mathbf k) = H^{\uparrow\downarrow}(\mathbf k)^\dagger \\).
- In real space, \\( H^{\sigma'\sigma}\_{J\beta, I\alpha}(-\mathbf L) = \overline{H^{\sigma\sigma'}\_{I\alpha, J\beta}(\mathbf L)} \\) holds.
- Since the overlap does not depend on spin, \\( S^{(2)} \\) is block diagonal with \\( S \\) repeated twice.
- For an eigenvector \\( c \\) (length \\( 2N \\)), the first \\( N \\) components are the spin-up coefficients and the last \\( N \\) components are the spin-down coefficients.

[Convention] The real-space density matrix of noncollinear calculations (spin 2×2 components) has the form
\\( \rho^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) = \sum w f e^{i\mathbf k\cdot\mathbf L}\\, \overline{c\_{\sigma I\alpha}}\\, c\_{\sigma' J\beta} \\), naturally derived from the Bloch-sum convention above.
For how the components are stored, see [density matrix and charge mixing](density-mixing.md).

## k-point mesh [Public] [Convention]

[Public] According to manual 4.0 (Chapter 15), `scf.Kgrid` specifies the numbers of divisions \\( (N\_1, N\_2, N\_3) \\) along the reciprocal lattice vectors,
and `scf.Generation.Kpoint` specifies the generation method. The choices are `GCenter` (default) and `Regular`.
For both even and odd numbers of divisions, the `GCenter` mesh contains Γ; it treats \\( \mathbf k \\) and \\( -\mathbf k \\) on an equal footing,
and includes the time-reversal-invariant momenta (TRIM) satisfying \\( \mathbf k \equiv -\mathbf k \pmod{\mathbf G} \\).

[Convention] The default when `scf.Kgrid` is omitted is `4 4 4` (the manual states that it is required for Band).
Keyword values are case-insensitive.

### Coordinates

For each direction \\( i \\) and index \\( l = 0, \dots, N\_i - 1 \\), the fractional coordinates are as follows.

| Method | \\( q\_{i,l} \\) | Remarks |
|---|---|---|
| `GCenter` | \\( l / N\_i \\) | Lie in \\( [0, 1) \\). \\( l = 0 \\) is Γ. If \\( N\_i \\) is even, \\( q = 1/2 \\) (boundary) is also included |
| `Regular` | \\( -\tfrac12 + \tfrac{2l + 1}{2N\_i} \\) | Includes Γ only when all \\( N\_i \\) are odd. Does not include points on the boundary |

- In both methods, \\( q\_i = 0 \\) for directions with \\( N\_i = 1 \\) (the `Regular` formula also gives 0 in this case).
- [Public] The `Regular` formula is identical to the Monkhorst–Pack formula (PRB **13**, 5188 (1976))
  \\( u\_r = (2r - N - 1)/(2N) \\), \\( r = 1, \dots, N \\), with \\( r = l + 1 \\).
  Hence `Regular` corresponds to "a Monkhorst–Pack-type grid (without symmetry reduction)".
- The `GCenter` points lie in \\( [0, 1) \\), but by the periodicity \\( O(\mathbf k + \mathbf G) = O(\mathbf k) \\) they give the same matrices as
  the Γ-centered grid folded into \\( (-\tfrac12, \tfrac12] \\).
  This gives the same set of eigenvalues as a Γ-point calculation of the \\( N\_1\times N\_2\times N\_3 \\) supercell (see Tests).
- [Convention] The existing implementation adds \\( +\delta \\), \\( -\delta \\), \\( +2\delta \\) (\\( \delta = 10^{-12} \\)), in the order of directions 1, 2, 3,
  to the fractional coordinates of directions with 2 or more divisions. The shift is far below any physically relevant scale
  and need not be reproduced in the reimplementation ([bitwise agreement is not required](../../../process/writing-rules.md)).

### Reduction of k points and weights [Convention]

For `Band` in collinear calculations (`scf.SpinPolarization` is `Off` or `On`), k points are reduced using only
\\( \varepsilon(\mathbf k) = \varepsilon(-\mathbf k) \\) from time reversal. The point-group symmetry of the crystal is not used.

- The partner \\( \bar{\mathbf l} \\) of a grid point \\( \mathbf l = (l\_1, l\_2, l\_3) \\) is the grid point corresponding to \\( -\mathbf q \\) (modulo \\( \mathbf G \\)).
  - `Regular`: \\( \bar l\_i = N\_i - 1 - l\_i \\)
  - `GCenter`: \\( \bar l\_i = (N\_i - l\_i) \bmod N\_i \\) (for \\( l\_i = 0 \\) and \\( 2l\_i = N\_i \\), the point is its own partner)
- Points with \\( \bar{\mathbf l} = \mathbf l \\) (TRIM) have weight 1; for the others, one member of each pair is kept with weight 2 and the other is removed.
  The one kept is the one that appears first in lexicographic order with \\( l\_1 \\) outermost and \\( l\_3 \\) innermost.
- The sum of weights is \\( N\_1 N\_2 N\_3 \\). Weights used for the density matrix etc. are divided by this and normalized to sum to 1.
- Number of TRIM: for `GCenter`, the product over directions of 2 if \\( N\_i \\) is even and 1 if odd.
  For `Regular`, only Γ (one point) when all \\( N\_i \\) are odd, and zero otherwise.

For `Band` in noncollinear calculations, no reduction is done, and all \\( N\_1 N\_2 N\_3 \\) points are used with weight 1.
This is because the noncollinear \\( H(\mathbf L) \\) contains complex numbers so that \\( H(-\mathbf k) = \overline{H(\mathbf k)} \\) does not hold,
and time reversal requires an antiunitary transformation accompanied by spin flip, which is also broken in the presence of magnetism [Public] (a general fact).

[Open]

- The existing implementation retains another scheme that reduces the Monkhorst–Pack grid using the symmetry of the atomic configuration (tolerance keyword `scf.MP.criterion`).
  Manual 4.0 documents only `GCenter` and `Regular`, and the value that selects this scheme in the existing implementation is not in the manual either.
  It can be used only with collinear `Band`, not with NEGF or noncollinear calculations. This note treats it as out of scope.
- Conventions for generating band diagrams (k-point paths) (Section 2 of [analysis](analysis.md)).

[Convention] k points used outside SCF are defined in [analysis](analysis.md). The density of states (`Dos.Kgrid`) uses a Γ-centered grid \\( q\_{i,l} = l/N\_i \\), independent of `scf.Generation.Kpoint`, without reduction by time reversal or otherwise.
Exchange coupling constants (jx) use a Monkhorst–Pack-type grid without reduction. Optical conductivity (`CDDF.Kgrid`) uses a Monkhorst–Pack-type grid, merging \\( \pm\mathbf k \\) into a single point with weight 2.

## Requirements for implementations [Design]

- The Bloch sum is a stateless function that multiplies each \\( (I, J, \mathbf n) \\) of the sparse block format by the phase \\( e^{2\pi i\\,\mathbf q\cdot\mathbf n} \\) and accumulates.
  The input is the fractional coordinates \\( \mathbf q \\), and the output is a dense matrix (column-major).
- The phase convention (lattice vectors only / including atomic positions) is made explicit by a type or option, with lattice vectors only as the default.
- For collinear calculations, a path that returns a real matrix at the Γ point is provided, using the fact that \\( O(\mathbf L) \\) is real.
- Spinor matrices of noncollinear calculations use the ordering \\( p = \sigma N + a \\), and \\( S^{(2)} \\) is not built explicitly but can be handled as \\( I\_2\otimes S \\).
- The set of k points is generated as "a sequence of fractional coordinates and normalized weights", with the generation method (`GCenter`/`Regular`) and
  whether to apply time-reversal reduction as arguments. Reduction is allowed only for collinear calculations.
- For comparison with the existing implementation, whether to add the tiny shift \\( \delta \\) may be an option (default: not added).

## Tests [Design]

`tools/check_kpoints_eigen.py` checks many of the following on a small model with random real \\( H(\mathbf L) \\), \\( S(\mathbf L) \\)
(\\( O(-\mathbf L) = O(\mathbf L)^{\mathsf T} \\)).

| Kind | Content | Approximate tolerance |
|---|---|---|
| Invariant | \\( H(\mathbf k) \\), \\( S(\mathbf k) \\) are Hermitian | Rounding error |
| Invariant | \\( O(\mathbf k + \mathbf G) = O(\mathbf k) \\) (shift fractional coordinates by integers) | \\( 10^{-12} \\) |
| Invariant | Collinear: \\( O(-\mathbf k) = \overline{O(\mathbf k)} \\) | Rounding error |
| Invariant | Eigenvalues are invariant when switching to the convention including atomic positions | \\( 10^{-11} \\) |
| Invariant | The sum over reduced k points and weights agrees with the sum over all grid points (for quantities satisfying \\( f(\mathbf k) = f(-\mathbf k) \\), e.g. the sum of occupied eigenvalues) | \\( 10^{-10} \\) |
| Invariant | The sum of weights is \\( N\_1 N\_2 N\_3 \\), and the number of TRIM follows the rule above | Exact |
| Invariant | The set of eigenvalues on the `GCenter` \\( N\_1\times N\_2\times N\_3 \\) grid agrees with the Γ-point eigenvalues of a supercell of the same size | \\( 10^{-10} \\) |
| Invariant | With \\( \rho(\mathbf L) \\) built from the reduced k points, \\( \sum \rho H \\) agrees with the band energy and \\( \sum \rho S \\) with the number of electrons | \\( 10^{-11} \\) |
| Invariant | Noncollinear: \\( H^{(2)} \\) is Hermitian, \\( S^{(2)} = I\_2 \otimes S \\) | Rounding error |
| Invariant | Noncollinear with zero spin-off-diagonal terms: the \\( 2N \\)-dimensional eigenvalues are the union of the collinear eigenvalues for the two spins | \\( 10^{-12} \\) |
| Comparison | Compare the list of fractional coordinates and weights of k points with the output of the existing implementation (`scf.Kgrid` equal to `2 3 4`, `3 3 3`, `4 4 1`, etc., both methods) | \\( 10^{-10} \\) (excluding \\( \delta \\)) |
| Comparison | \\( H(\mathbf L) \\), \\( S(\mathbf L) \\) written out from the existing implementation, and \\( H(\mathbf k) \\), \\( S(\mathbf k) \\) at several k points (collinear and noncollinear) | To be decided |
