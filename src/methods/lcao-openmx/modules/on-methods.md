# O(N) methods: divide-and-conquer, DC-LNO, and the Krylov subspace method

> Status: first draft (awaiting review). The ideas and equations rewrite technical note 3 (Ver. 1.1, Krylov subspace method; hereafter TN3),
> technical note 26 (DC-LNO; hereafter TN26), technical note 16 (low-order scaling method; hereafter TN16), and manual Chapters 27 and 48
> in the notation of this note. Conventions were confirmed as the behavior of the existing implementation (the GPU version of OpenMX 4.0).
> Numerical properties are checked with `tools/check_on_methods.py`.

## Purpose and scope

Instead of solving the generalized eigenvalue problem of the whole system ([Generalized eigenvalue problem](generalized-eigen.md)),
these solvers solve a small "truncated cluster" problem for each atom and build only the rows of the density matrix involving that atom.
The following solvers selected by `scf.EigenvalueSolver` are covered.

- `DC`: divide-and-conquer method. Solves the cluster generalized eigenvalue problem directly.
- `Krylov`: projects the cluster problem onto a Krylov subspace before solving.
- `DC-LNO`: contracts the basis of atoms in the outer part of the cluster to localized natural orbitals (LNOs) before solving.
- `Cluster2` (the low-order scaling method of manual Chapter 48) and other legacy values are treated briefly.

The procedure for building the grid density from the density matrix, mixing, and the expressions for the total energy and forces themselves are
as in the chapters [Density matrix and charge mixing](density-mixing.md) and [Total energy and forces](energy-force.md);
this chapter specifies how to build the density matrix \\( \rho \\) and the energy-weighted density matrix \\( E \\) passed to them.

## Notation [Design]

In addition to the symbols of [Notation](../../../conventions/notation.md), [Sparse block format](../design/sparse-block-format.md), and [Total energy and forces](energy-force.md), the following are used.

| Symbol | Meaning |
|---|---|
| \\( I \\) | Central atom (an atom in the unit cell). One local problem per \\( I \\) |
| \\( s = (J, \mathbf L) \\) | Cluster "site": the periodic image \\( \mathbf R\_J + \mathbf L \\) of atom \\( J \\) |
| \\( d\_{Is} = \lVert \mathbf R\_J + \mathbf L - \mathbf R\_I \rVert \\) | Distance from the central atom to the site |
| \\( r\_{c,J} \\) | PAO cutoff radius of atom \\( J \\) (the largest among its basis functions) |
| \\( R\_{\mathrm{hop}} \\) | `orderN.HoppingRanges` (input in Å, bohr internally) |
| \\( \mathcal F(I) \\) | Set of first-neighbor sites: \\( d\_{Is} \le r\_{c,I} + r\_{c,J} \\) (includes the central atom itself \\( (I, \mathbf 0) \\)) |
| \\( \mathcal S(I) \\) | Set of second-neighbor sites: \\( r\_{c,I} + r\_{c,J} < d\_{Is} \le R\_{\mathrm{hop}} \\) |
| \\( \mathcal C(I) = \mathcal F(I) \cup \mathcal S(I) \\) | Truncated cluster |
| \\( H^{(I)}, S^{(I)} \\) | Hamiltonian and overlap matrix of the cluster (dimension \\( n\_I \\)) |
| \\( \varepsilon^{(I)}\_\mu, c^{(I)}\_\mu \\) | Eigenvalues and eigenvectors of the local problem |
| \\( g\_s \\) | Spin degeneracy. 2 for non-spin-polarized, 1 for collinear spin-polarized and noncollinear |
| \\( N\_{\mathrm e} = \sum\_I Z\_I^{\mathrm{val}} - Q \\) | Number of electrons to occupy ([Density matrix and charge mixing](density-mixing.md)) |
| \\( f(\varepsilon) = 1/(1+e^{(\varepsilon-\mu)/k\_{\mathrm B}T}) \\) | Fermi–Dirac distribution |

In the manual and standard output, the sizes of \\( \mathcal F(I) \\) and \\( \mathcal S(I) \\) are called FNAN and SNAN, respectively.

## Choice of solver [Public] [Convention]

[Public] According to manual Sec. 8.2 and Chapter 27, specifying `DC`, `DC-LNO`, or `Krylov` for `scf.EigenvalueSolver` selects an O(N) method,
and the low-order scaling method of Chapter 48 is selected with `cluster2`. Values are case-insensitive [Convention].

[Convention] The values accepted by the existing implementation and their treatment with respect to O(N) methods are as follows.

| Value | Content | Spin | Remarks |
|---|---|---|---|
| `DC` | Divide-and-conquer method | `Off`, `On`, `NC` | Stops if spin-orbit coupling is specified with collinear spin |
| `Krylov` | Krylov subspace method | `Off`, `On` | `NC` prints "not supported" and stops |
| `DC-LNO` | Divide-and-conquer with LNOs | `Off`, `On`, `NC` | |
| `Cluster2` | Low-order scaling method (manual Chapter 48) | `Off`, `On` | With `NC` returns without computing anything (see "Edge cases" below) |
| `Cluster-LNO` | Cluster calculation contracted with LNOs | | `LNOs.Num` is required. Not treated in this chapter |
| `EGAC` | Solver not described in the public documents | | [Open] |
| `Recursion`, `GDC` | Legacy solvers | | Prints "not supported" and stops |
| `Cluster-DIIS` | Legacy solver | | Accepted, but the eigenvalue problem is not solved in the SCF (effectively unusable) |

- If `scf.EigenvalueSolver` is omitted, the existing implementation treats it as `Recursion` and stops. The keyword is effectively mandatory.
- [Public] Manual Sec. 8.2 gives the value for the low-order scaling method as `ON2`, but the existing implementation does not accept `ON2` (input error);
  `cluster2` of Chapter 48 is the correct value.
- DC, Krylov, and DC-LNO treat real-space periodic images directly as sites, so k-points (`scf.Kgrid`) are not used.
  The procedure is the same for periodic and isolated systems [Convention].

## 1. Truncated cluster

### Definition [Public]

According to TN3 Sec. 2.3 and manual Sec. 27.1, a truncated cluster is built from the atoms inside a sphere of radius \\( R\_{\mathrm{hop}} \\) centered on each atom,
and an independent small problem is solved for each atom. According to manual Sec. 27.4, atoms whose basis functions overlap with those of the central atom are called first neighbors,
and the other atoms of the cluster second neighbors.

Note on [Public]: the English text of manual Sec. 27.4 reads as "atoms for which the sum of cutoff radii is smaller than the interatomic distance are first neighbors", but
there is overlap (i.e. the atom is a first neighbor) when the interatomic distance is smaller than the sum of the cutoff radii. The existing implementation also uses the latter (next section).
In addition, TN3 describes a logical truncation by "number of hoppings" besides truncation by distance, but
the existing implementation only reads `orderN.NumHoppings` (default 2) and does not use it [Convention].

### Site selection [Convention]

- The first neighbors \\( \mathcal F(I) \\) are all sites satisfying \\( d\_{Is} \le r\_{c,I} + r\_{c,J} \\) (including equality).
  This is the same as the range of atom pairs in the [Sparse block format](../design/sparse-block-format.md) and does not depend on \\( R\_{\mathrm{hop}} \\).
- The second neighbors \\( \mathcal S(I) \\) are sites satisfying \\( r\_{c,I} + r\_{c,J} < d\_{Is} \le R\_{\mathrm{hop}} \\).
  If \\( R\_{\mathrm{hop}} \\) is smaller than \\( r\_{c,I} + r\_{c,J} \\), the cluster consists of first neighbors only.
- Sites are periodic images of atoms. If the unit cell is small, multiple images of the same atom (including images of the central atom itself with \\( \mathbf L \ne 0 \\)) enter
  the cluster as separate sites.
- The order of sites in the cluster is: the central atom first, then first neighbors, then second neighbors.
  In DC-LNO, the second neighbors are re-sorted in ascending order of distance (the order for equal distances is [Open]).

### Cluster matrices [Convention]

\\( H^{(I)} \\) and \\( S^{(I)} \\) are the principal submatrices of the whole-system matrices restricted to the cluster. The block for sites \\( s = (J, \mathbf L) \\), \\( s' = (J', \mathbf L') \\) is

$$
H^{(I)}\_{s\alpha, s'\beta} = H\_{J\alpha, J'\beta}(\mathbf L' - \mathbf L),\qquad
S^{(I)}\_{s\alpha, s'\beta} = S\_{J\alpha, J'\beta}(\mathbf L' - \mathbf L)
$$

and is nonzero only if \\( s' \\) is a first neighbor of \\( s \\) (\\( \lVert \mathbf R\_{J'} + \mathbf L' - \mathbf R\_J - \mathbf L\rVert \le r\_{c,J} + r\_{c,J'} \\)); otherwise it is 0.
\\( H \\) is the Kohn–Sham Hamiltonian (per spin) of each SCF iteration; if orbitals are contracted, it is the matrix in the contracted basis.
Connections to outside the cluster are discarded, so \\( H^{(I)} \\) is not a problem of "the cluster cut out into vacuum":
the Hartree and exchange-correlation potentials are built from the density of the whole system (TN3 Sec. 2.3).

### Fixing the number of sites: `orderN.FNAN+SNAN` [Public] [Convention]

[Public] Manual Sec. 27.4: when the lattice constant etc. is changed, \\( \lvert\mathcal C(I)\rvert \\) jumps and the energy curve becomes discontinuous, so
the number of cluster sites can be specified per atom.

[Convention] The behavior of the existing implementation is as follows.

- Each line of the block reads **three integers** \\( (I, a, b) \\), and the target number of sites is \\( a + b \\)
  (lines are in order of atom number; an input error occurs if the first column does not match the atom number).
  The manual's example has two columns (atom number and FNAN+SNAN), but that format is an input error in the existing implementation.
- When \\( a + b \\) is less than or equal to the number of sites determined by distance (\\( \lvert\mathcal F\rvert + \lvert\mathcal S\rvert \\) excluding the central atom),
  the second neighbors are sorted in ascending order of distance and only the first \\( a + b - \lvert\mathcal F(I)\rvert \\) (0 if negative) are kept. First neighbors are never removed.
  When \\( a + b \\) is larger, nothing is done (the cluster is not enlarged).
- This adjustment is made only at the first MD step; from the second step onward the distance-based selection is used again ([Open]: whether this is intended).

## 2. Divide-and-conquer method (DC)

### Local problem [Public]

Following TN3 Eqs. 3–4 and TN26 Eq. 27, for each atom \\( I \\)

$$
H^{(I)} c^{(I)}\_\mu = \varepsilon^{(I)}\_\mu S^{(I)} c^{(I)}\_\mu,\qquad
G^{(I)}(Z) = \sum\_\mu \frac{c^{(I)}\_\mu c^{(I)\dagger}\_\mu}{Z - \varepsilon^{(I)}\_\mu}
$$

is solved, and the density matrix is built using only the rows of the Green's function belonging to the central atom (TN3 calls this "patching").
Writing the residues as

$$
r^{(I)}\_{\mu;\alpha, s\beta} = c^{(I)}\_{(I,\mathbf 0)\alpha,\mu}\\, c^{(I)}\_{s\beta,\mu},\qquad s \in \mathcal F(I)
$$

the integral over the Fermi distribution becomes a sum over residues,

$$
\rho\_{I\alpha, J\beta}(\mathbf L) = \sum\_\mu f(\varepsilon^{(I)}\_\mu)\\, r^{(I)}\_{\mu;\alpha, s\beta},\qquad
E\_{I\alpha, J\beta}(\mathbf L) = \sum\_\mu f(\varepsilon^{(I)}\_\mu)\\, \varepsilon^{(I)}\_\mu\\, r^{(I)}\_{\mu;\alpha, s\beta}
$$

(\\( s = (J, \mathbf L) \\), per spin; TN26 Eqs. 24–25). In collinear calculations the local problem is real symmetric and \\( c \\) is real.

### Treatment of the overlap matrix [Convention]

- Orthogonalize with rule C of [Generalized eigenvalue problem](generalized-eigen.md) (remove eigenvalues of \\( S^{(I)} \\) below \\( 10^{-9} \\) and their eigenvectors), and
  find **all eigenpairs** of the standard eigenvalue problem in the remaining subspace. Eigenvalues corresponding to the removed directions are set to 1000 Hartree and treated as having occupation 0.
- The decomposition of \\( S^{(I)} \\) is built in the first two SCF iterations and the saved one is used thereafter (the result is the same as long as the atomic configuration does not change).

### Electron number and chemical potential [Public] [Convention]

[Public] TN3 Eqs. 29–31: the central atom is given the Mulliken weight

$$
\kappa^{(I)}\_\mu = \sum\_{\alpha}\sum\_{s=(J,\mathbf L)\in\mathcal F(I)}\sum\_\beta r^{(I)}\_{\mu;\alpha,s\beta}\\, S\_{I\alpha, J\beta}(\mathbf L)
$$

and the chemical potential \\( \mu \\), **common** to the whole system, is determined from

$$
N(\mu) = g\_s \sum\_\sigma \sum\_I \sum\_\mu f(\varepsilon^{(I,\sigma)}\_\mu)\\, \kappa^{(I,\sigma)}\_\mu = N\_{\mathrm e}
$$

Electrons can move between clusters, and the electron number of each atom is not fixed. Only the sum of per-atom electron numbers is communicated.
Since it is a sum of the form \\( \sum\_I\sum\_{\mathcal F(I)} \rho S \\), the resulting \\( \rho \\) satisfies exactly the sum rule
\\( g\_s \sum\_\sigma \sum \rho S = N\_{\mathrm e} \\) of [Density matrix and charge mixing](density-mixing.md) (`tools/check_on_methods.py`).

[Convention] \\( \mu \\) is found by bisection.

| Solver | Initial interval (Hartree) | Stopping criterion \\( \lvert N(\mu) - N\_{\mathrm e}\rvert \\) | Clipping of the Fermi-distribution exponent |
|---|---|---|---|
| DC (collinear) | \\( [-30, 30] \\) | \\( < 10^{-13} \\) | \\( x \in [-30, 30] \\) |
| DC (noncollinear) | \\( [-15, 15] \\) | \\( < 10^{-11} \\) | \\( x \in [-50, 50] \\) |
| Krylov | \\( [-10, 10] \\) | \\( < 10^{-11} \\) | \\( x \in [-30, 30] \\) |
| DC-LNO | \\( [-10, 10] \\) in SCF iterations 1 and 2, thereafter the previous \\( \mu \\) \\( \pm 5 \\) | \\( < 10^{-12} \\) | \\( x \in [-30, 30] \\) |

The iteration limit is 1000 in all cases, and \\( x = (\varepsilon - \mu)/k\_{\mathrm B}T \\).
As with Cluster and Band, it is undetermined where in the gap \\( \mu \\) stops for insulators ([Density matrix and charge mixing](density-mixing.md)).

### Band energy [Convention]

\\( U\_{\mathrm{ele}} = g\_s\sum\_\sigma\sum\_I\sum\_\mu f(\varepsilon^{(I,\sigma)}\_\mu)\\, \varepsilon^{(I,\sigma)}\_\mu\\, \kappa^{(I,\sigma)}\_\mu \\) is computed, and
\\( \sum\_{\mathrm{pairs}} \rho H \\) is also computed separately. Unlike full diagonalization, the two do not agree in O(N) methods (because the weight \\( \kappa \\) uses \\( S \\), not \\( H \\)).

### Noncollinear [Convention]

In the noncollinear case (`scf.SpinPolarization NC`), the local problem is a \\( 2n\_I \\)-dimensional complex Hermitian matrix with overlap \\( I\_2\otimes S^{(I)} \\).
The treatment of \\( S^{(I)} \\) (rule C), residues, Mulliken weights, and the determination of \\( \mu \\) have the same form as in the collinear case, and \\( \rho \\) and \\( E \\) are built per 2×2 spin block
(for how spin components are stored, see the noncollinear section of [Density matrix and charge mixing](density-mixing.md)).

### Asymmetry of the density matrix [Convention]

\\( \rho\_{I\alpha,J\beta}(\mathbf L) \\) is built from the cluster of atom \\( I \\), and \\( \rho\_{J\beta,I\alpha}(-\mathbf L) \\) from the cluster of atom \\( J \\), so
in general \\( \rho\_{I\alpha,J\beta}(\mathbf L) \ne \rho\_{J\beta,I\alpha}(-\mathbf L) \\) (they agree only in the limit where the cluster is extended to the whole system).
The existing implementation does not symmetrize; it keeps both blocks as they are and uses them for the density, energy, and forces.
The grid density \\( \sum \rho\_{I\alpha,J\beta}\phi\_{I\alpha}\phi\_{J\beta} \\) effectively uses only the symmetric part.

## 3. Krylov subspace method

### Idea [Public]

TN3 Eqs. 5–8: from \\( A = (S^{(I)})^{-1} H^{(I)} \\) and a set of starting vectors \\( W\_0 \\) (\\( n\_I \times m \\)) consisting of orbitals near the center,

$$
\mathcal K\_q = \mathrm{span}\\{ W\_0, A W\_0, A^2 W\_0, \dots, A^{q} W\_0 \\}
$$

is built, and the standard eigenvalue problem projected onto a basis \\( U \\) orthonormal in the \\( S^{(I)} \\) inner product (\\( U^{\mathsf T} S^{(I)} U = 1 \\)),

$$
H\_K b\_\mu = \varepsilon\_\mu b\_\mu,\qquad H\_K = U^{\mathsf T} H^{(I)} U,\qquad c\_\mu = U b\_\mu
$$

is solved. The projected Green's function reproduces exactly the moments \\( W\_0^{\mathsf T} S A^{p} W\_0 \\) sandwiched by \\( W\_0 \\) up to \\( p \le 2q+1 \\).
[Design] (checked) With \\( t = U^{\mathsf T} S W\_0 \\), \\( W\_0^{\mathsf T} S A^{p} W\_0 = t^{\mathsf T} H\_K^{p} t \\) holds (\\( p \le 2q+1 \\)),
and in general does not hold for \\( p = 2q+2 \\) (`tools/check_on_methods.py`; with \\( q = 3 \\), relative error \\( 10^{-15} \\) for \\( p \le 7 \\) and about \\( 10^{-11} \\) for \\( p = 8 \\)).

TN3 reports that the Krylov dimension needed for convergence of the density matrix is about 40–50% of the number of basis functions of the cluster for insulators.

### Starting vectors (core) [Public] [Convention]

[Public] Manual Sec. 27.3: if `orderN.Expand.Core` is `on` (default), the atoms within a sphere of 1.2 times the nearest-neighbor distance \\( r\_{\min} \\),
and if `off`, only the central atom, form the "core", whose orbitals are used as starting vectors.

[Convention] The rules of the existing implementation are as follows; the factor differs slightly from the manual's description (1.2).

- \\( r\_{\min} \\) is the minimum distance among the first neighbors (excluding the central atom).
- For `on`, first the sites with \\( d\_{Is} < 1.10\\, r\_{\min} \\) (including the central atom, chosen from the first neighbors) form the core,
  and if the number of their orbitals is 30 or less, the selection is redone with the factor multiplied by 1.2, \\( d\_{Is} < 1.32\\, r\_{\min} \\) (redone only once).
- For `off`, the factor is 0.80 and the core is only the central atom.
- The starting vectors \\( W\_0 \\) are the unit vectors of the core orbitals (\\( m \\) = number of core orbitals), \\( S \\)-orthonormalized by the method of the next subsection.

### Krylov dimension [Public] [Convention]

[Public] Manual Sec. 27.3: `orderN.KrylovH.order` is the dimension of the Krylov subspace of the Hamiltonian,
and `orderN.KrylovS.order` (default 4 times `orderN.KrylovH.order`) is the dimension of the Krylov subspace of the overlap matrix.

[Convention] The dimension is determined by the number of blocks. With \\( n\_I \\) the number of basis functions of the cluster,

$$
q + 1 = \min\left( \left\lceil \frac{K\_H}{m} \right\rceil,\ \left\lceil \frac{n\_I}{m} \right\rceil \right),\qquad
\dim U = (q+1)\\, m
$$

(\\( K\_H \\) is `orderN.KrylovH.order`, default 400). The Krylov subspace of \\( S \\) is determined from \\( K\_S \\) by the same rule.

- If \\( \dim U < n\_I \\), the Krylov subspace of the next subsection is built.
- If \\( \dim U \ge n\_I \\), no Krylov subspace is built; instead \\( U \\) is the \\( S \\)-orthonormalized set of unit vectors of all cluster orbitals
  (directions dropped with the same threshold as rule C). In this case the result is the same as DC (except for the update of the embedding matrix).

### Basis generation and orthogonalization [Public] [Convention]

[Public] The block Gram–Schmidt method in the \\( S \\) inner product of TN3 Eqs. 9–16.

[Convention] The existing implementation builds \\( U \\) in the following order.

1. \\( S \\)-orthonormalize \\( W\_0 \\): diagonalize \\( G = W^{\mathsf T} S W \\), set directions with eigenvalue \\( g \le 10^{-9} \\) to 0, and \\( W \leftarrow W X g^{-1/2} \\).
2. For \\( k = 0, \dots, q-1 \\), form \\( V = A W\_k \\) (how \\( A \\) is formed is in the next subsection), apply against all previous blocks
   \\( V \leftarrow V - W\_j (W\_j^{\mathsf T} S V) \\) (classical block Gram–Schmidt, once), and \\( S \\)-orthonormalize by the same method as step 1 to obtain \\( W\_{k+1} \\).
3. For \\( U\_0 = (W\_0, \dots, W\_q) \\), all blocks concatenated, diagonalize \\( F = U\_0^{\mathsf T} S U\_0 \\),
   drop directions with eigenvalue \\( \le 10^{-9} \\), and set \\( U = U\_0 X F^{-1/2} \\) (the final orthogonal transformation of TN3 Eq. 9).
   Columns of the dropped directions are left as 0, and their number \\( z \\) is recorded.
4. \\( U \\) is built with the Hamiltonian of the **first SCF iteration** of the Krylov subspace and is not rebuilt in subsequent iterations with the same atomic configuration
   (it is rebuilt in the first iteration of each MD or geometry-optimization step). A separate \\( U \\) is kept per spin.

The threshold \\( 10^{-9} \\) is the same fixed value as rule C of [Generalized eigenvalue problem](generalized-eigen.md).

### Inverse of the overlap matrix [Public] [Convention]

[Public] Manual Sec. 27.3: if `orderN.Exact.Inverse.S` is `on` (default), \\( (S^{(I)})^{-1} \\) is computed exactly;
if `off`, it is approximated in the Krylov subspace of \\( S \\) (dimension `orderN.KrylovS.order`). According to TN3 Eqs. 17–26, the approximate inverse is

$$
Q = V (V^{\mathsf T} S V)^{-1} V^{\mathsf T},\qquad V = \text{span}\\{ V\_0, S V\_0, \dots, S^{q'} V\_0 \\}
$$

and \\( V \\) is orthonormalized in the ordinary (Euclidean) inner product.

[Convention]

- `on`: the inverse is built from the Cholesky decomposition of \\( S^{(I)} \\). No regularization is applied (if not positive definite, an error is printed and the calculation continues).
  The inverse is built only once, in the first SCF iteration.
- `off`: \\( V\_0 \\) consists of the same core unit vectors as \\( W\_0 \\); \\( V \\) is generated in the same way as steps 1–2 above with the Euclidean inner product;
  \\( (V^{\mathsf T} S V)^{-1} \\) is a pseudo-inverse dropping directions with eigenvalue \\( \le 10^{-9} \\).
- The internal flag `orderN.Inverse.S` is overwritten with the opposite value of `orderN.Exact.Inverse.S`, so specifying it in the input has no effect.
- [Public] Manual Sec. 27.3 states that convergence may improve by setting `orderN.Exact.Inverse.S` and `orderN.Expand.Core` to `on` for covalent systems,
  and both to `off` for simple metals.

### Projected Hamiltonian and embedding [Public] [Convention]

[Public] TN3 Eqs. 27–28: dividing the cluster into a "core" part and a "buffer" part, with the corresponding rows of \\( U \\) denoted \\( u\_c \\) and \\( u\_b \\),

$$
H\_K = u\_c^{\mathsf T} H\_{cc} u\_c + \left( u\_c^{\mathsf T} H\_{cb} u\_b + u\_b^{\mathsf T} H\_{bc} u\_c + u\_b^{\mathsf T} H\_{bb} u\_b \right)
$$

and the part in parentheses (the long-range part) can be computed in the first SCF iteration and kept fixed.

[Convention]

- Here the "core" is the first neighbors \\( \mathcal F(I) \\) (including the central atom) and the "buffer" is the second neighbors \\( \mathcal S(I) \\)
  (a different concept from the "core" of the starting vectors).
- If `orderN.Recalc.Buffer` is `on` (default), the part in parentheses is rebuilt every iteration. Then \\( H\_K = U^{\mathsf T} H^{(I)} U \\) holds exactly
  (`tools/check_on_methods.py`).
- If `off`, it is computed in the first iteration and fixed thereafter. However, it is rebuilt once in the iteration where the square root of the residual norm
  (\\( \lVert R\rVert \\) of [Density matrix and charge mixing](density-mixing.md)) first falls below \\( 0.2 + 0.1\\, N\_{\mathrm{atom}} \\).
- The \\( z \\) dropped directions (zero columns of \\( U \\)) are handled as unoccupied eigenvalues by setting the corresponding diagonal elements of \\( H\_K \\) to 1000 Hartree.

### Eigenvalues, density matrix, and chemical potential [Convention]

- All eigenpairs of \\( H\_K \\) are found, and only the core (\\( \mathcal F(I) \\)) rows of \\( c\_\mu = U b\_\mu \\) are used to build the same residues and Mulliken weights as in DC.
- To reduce cost, the eigenvalues are split into three groups by the energy window \\( [\mu' - \Delta, \mu' + \Delta] \\) (\\( \Delta = 10 \\) eV \\( = 0.367493245 \\) Hartree).
  \\( \mu' \\) is the chemical potential at that time (the value of the previous iteration; in the first iteration `orderN.InitChemPot`, default \\( -3.0 \\) eV).
  - States with eigenvalues below the lower end are added to \\( \rho \\), \\( E \\), and the electron number with occupation 1 (fixed while searching for \\( \mu \\)).
  - For states inside the window (up to the first eigenvalue exceeding the upper end), residues are stored and weighted with the Fermi distribution.
  - States above that are discarded (occupation 0).
- \\( \mu \\) is found by bisection of the same equation as DC (intervals etc. in the table above).
- If in the first iteration the true \\( \mu \\) is 10 eV or more away from `orderN.InitChemPot`, the occupations in that iteration are badly wrong. From the second iteration onward the window follows.

## 4. Divide-and-conquer with localized natural orbitals (DC-LNO)

### Definition of LNOs [Public]

TN26 Eqs. 17–22: for atom \\( I \\), form the diagonal block of the product of the density matrix and the overlap matrix,

$$
\Lambda\_I = \sum\_{(J,\mathbf L)\in\mathcal F(I)} \rho\_{I, J}(\mathbf L)\\, S\_{J, I}(-\mathbf L)
$$

(\\( \rho\_{I,J}(\mathbf L) \\) is the block of orbitals of \\( I \\) × orbitals of \\( J \\)); then \\( g\_s\sum\_\sigma\sum\_I \mathrm{tr}\\,\Lambda\_I = N\_{\mathrm e} \\).
Since \\( \Lambda\_I \\) is not necessarily symmetric, the nonsymmetric eigenvalue problem \\( \Lambda\_I v\_{I\gamma} = \lambda\_{I\gamma} v\_{I\gamma} \\) is solved, and
the eigenvectors whose eigenvalues \\( \lambda \\) (corresponding to orbital occupations) are at or above a threshold are the LNOs.
An LNO is a linear combination of PAOs on atom \\( I \\) only, and gives a local low-rank approximation spanning the subspace of occupied states.
In the TN26 benchmarks, with threshold 0.1, C and Si have 4 per atom, Al 4 (9 to reproduce up to the conduction band), and Li 2 (5 likewise).

Note on [Public]: in the TN26 text, there is a place where the value of \\( \kappa \\), the fraction of second neighbors using PAOs, reads as "310", but
from the domain (0 to 1) and the numerical example of Eq. 28 it means 3/10.

### Conventions for computing LNOs [Convention]

- \\( \Lambda\_I \\) is rebuilt every iteration from the **input** density matrix of that SCF iteration (\\( \rho \\) built in the previous iteration and passed through mixing).
  In collinear spin-polarized calculations, a separate \\( \Lambda\_I \\) and separate LNOs are kept per spin.
  In noncollinear calculations, \\( \Lambda\_I \\) is built from \\( \rho^{\uparrow\uparrow} + \rho^{\downarrow\downarrow} \\), and the same LNOs are used for both spin components.
- The eigenvalue problem is solved as a general real nonsymmetric matrix, and eigenvalues are sorted in descending order of their **real part**. Imaginary parts are discarded; there is no special treatment of complex-conjugate pairs.
- Each eigenvector is normalized to unit Euclidean norm (not \\( S \\)-normalized).
- The number of LNOs \\( n\_{\mathrm L}(I) \\) is the count of \\( \lambda \\), sorted in descending order, up to the position where it first falls below `orderN.LNO.Occ.Cutoff` (default 0.1).
  - In collinear spin-polarized calculations, the larger of the counts of the two spins is used (the number of LNOs is common to both spins).
  - In noncollinear calculations the threshold is doubled (0.2 by default).
  - [Public] Manual Sec. 27.2: specifying the `LNOs.Num` block (number per species) fixes the number regardless of the threshold.
    This is recommended for geometry optimization and MD so that the number of LNOs does not change midway.
- [Public] `scf.dclno.threading` of manual Sec. 27.2 (switch for OpenMP parallelization) is only read and not used in the existing implementation (the version consulted).

### PAO region and LNO region [Public] [Convention]

[Public] TN26 Eq. 26 and manual Sec. 27.2: the inner part of the cluster (the first neighbors and some of the second neighbors close to the center) is kept in PAOs,
and the remaining second neighbors are represented by LNOs. `orderN.LNO.Buffer` (default 0.2) is a guide to the fraction of second-neighbor atoms using PAOs:
0 means (almost) all LNO, and 1 means all PAO (the same as DC). The manual recommends 0.1–0.2.

[Convention] The number \\( n\_P \\) of second neighbors using PAOs is determined as follows. Let \\( r\_0 \\) be the maximum distance of the first neighbors, \\( b \\) be `orderN.LNO.Buffer`,
and \\( \lvert\mathcal S\rvert \\) the number of second neighbors.

$$
r(\phi) = r\_0 + \phi\\,(R\_{\mathrm{hop}} - r\_0),\qquad n(\phi) = \\#\\{ s\in\mathcal S(I) : d\_{Is} \le r(\phi) \\}
$$

Checking \\( \phi = 0.1 \times 1.05^{k} \\) (\\( k = 0, 1, \dots, 399 \\)) in order, the first \\( n(\phi) \\) with \\( n(\phi) > \lfloor b \lvert\mathcal S\rvert \rfloor \\) is taken as \\( n\_P \\)
(if not satisfied in 400 steps, the last \\( n(\phi) \\)). The first \\( n\_P \\) second neighbors in ascending order of distance use PAOs, and the rest use LNOs.

- Since \\( \phi \\) is increased in 5% steps, \\( n\_P \\) is a number at least \\( \lfloor b \lvert\mathcal S\rvert \rfloor + 1 \\) determined by the boundaries of distance shells.
- Even with \\( b = 0 \\), at least one second neighbor (all of them within \\( r(0.1) \\), if any) uses PAOs.
- With \\( b = 1 \\) the condition cannot be satisfied, and \\( \phi \\) becomes large enough that all use PAOs.
- This classification is made when the atomic configuration is determined (each MD step) and is not changed during the SCF.

### Contracted local problem [Public] [Convention]

[Public] TN26 Eq. 26: the cluster matrices are built divided into PAO–PAO, PAO–LNO, LNO–PAO, and LNO–LNO blocks.

[Convention] With per-site transformation matrices \\( T\_s \\), the identity for PAO sites and, for an LNO site \\( s = (J, \mathbf L) \\),
the \\( M\_J \times n\_{\mathrm L}(J) \\) matrix of the LNO coefficients of atom \\( J \\), and \\( T = \bigoplus\_s T\_s \\),

$$
\tilde H^{(I)} = T^{\mathsf T} H^{(I)} T,\qquad \tilde S^{(I)} = T^{\mathsf T} S^{(I)} T
$$

is solved (in spin-polarized calculations, per-spin LNOs are used, so \\( \tilde S \\) also differs per spin).

- Orthogonalization uses rule D of [Generalized eigenvalue problem](generalized-eigen.md) (eigenvalues of \\( \tilde S \\) replaced by their absolute values, no threshold).
  Nearly linearly dependent LNOs make \\( \lvert\lambda\rvert^{-1/2} \\) large.
- In the first two iterations all eigenpairs are computed; thereafter only as many as the index of the upper end of the previous iteration's window plus 100 (capped at \\( \tilde n\_I \\)).
- The energy window is \\( \pm 10 \\) eV as in the Krylov method: occupation 1 below the lower end, Fermi distribution inside the window.
  However, 30 is added to the index of the upper end, capped at the number of computed eigenpairs.
- The central atom and the first neighbors are always in the PAO region, so all indices of the residues \\( r^{(I)}\_{\mu;\alpha,s\beta} \\) are PAOs, and
  the expressions for \\( \rho \\), \\( E \\), and the Mulliken weights are the same as in DC.
- \\( \mu \\) is determined by the same equation as DC, with intervals etc. as in the table above.
- [Design] (checked) If all LNOs are used (\\( n\_{\mathrm L} = M\_J \\)), \\( T \\) is nonsingular and the result agrees with DC (`tools/check_on_methods.py`).

[Public] TN26 Eq. 28 estimates the ratio of computation time to DC as \\( [M\_P(N\_F + \kappa N\_S) + M\_L(1-\kappa)N\_S]^3 / [M\_P^3 (N\_F+N\_S)^3] \\)
(\\( M\_P \\), \\( M\_L \\) are the numbers of PAOs and LNOs per atom, and \\( N\_F \\), \\( N\_S \\) the numbers of first and second neighbors).

## 5. Other solvers

### Low-order scaling method (`Cluster2`) [Public] [Convention]

[Public] TN16 and PRB 82, 075131 (2010), manual Chapter 48: a numerically exact method that, without truncation as in O(N) methods, computes only the elements of the whole-system \\( G(Z) = (ZS - H)^{-1} \\)
where \\( H \\) and \\( S \\) are nonzero (the elements needed for the density matrix).

- The Fermi distribution is represented as a sum over poles obtained from a continued-fraction expansion, and \\( \rho \\) is expressed as a sum of \\( G(\alpha\_p) \\) at finitely many complex energies \\( \alpha\_p = \mu + i z\_p k\_{\mathrm B} T \\) plus
  a moment term evaluated at a large imaginary energy (TN16 Eqs. 7–8). The energy-weighted density matrix is analogous (Eqs. 10–17).
- The needed elements of \\( G(\alpha\_p) \\) are obtained by ordering the basis by recursive spatial bisection (nested dissection), block \\( LDL^{\mathsf T} \\) factorization, and
  recurrences for Schur complements (Eqs. 18–39). The cost is \\( O(N(\log N)^2) \\) in 1D, \\( O(N^2) \\) in 2D, and \\( O(N^{7/3}) \\) in 3D.
- The chemical potential is determined by combining bisection, linear interpolation, Muller's method, and extrapolation from previous iterations (Sec. 2.3.3).
- The number of poles needed does not depend on system size; it is determined by the spectral radius and the temperature. The manual says about 100 poles suffice at 300 K or above.
- Only collinear spin, only isolated systems or the Γ point of periodic systems, and no output of densities of states or wave functions are supported (manual Chapter 48).

[Convention] The number of poles is specified by `scf.Npoles.ON2` (default 100). The energy-weighted density matrix is built after SCF convergence in a separate stage for the force calculation.
In the noncollinear case it returns without computing anything.

### Other values [Convention]

- `EGAC`: reads `scf.Npoles.EGAC` (default 80), `scf.DIIS.History.EGAC` (2), `scf.AC.flag.EGAC` (off), and `scf.GF.EGAC` (off), but
  is not described in the manual and is not treated in this chapter ([Open]).
- `Cluster-LNO`: full diagonalization with the basis contracted by LNOs. `LNOs.Num` is required, and regularization is rule A.
- The keywords of legacy O(N) methods (recursion method etc.) `orderN.RecursiveLevels` (10), `orderN.TerminatorType`, `orderN.InverseSolver`,
  `orderN.InvRecursiveLevels` (30), `orderN.ChargeDeviation` (\\( 10^{-7} \\)), `orderN.AvNumTerminater` (1), `orderN.EC.Sub.Dim` (400),
  `orderN.NumHoppings` (2) are read but do not affect the solvers above. `orderN.NumPoles` (10) is used for another purpose ([Open]).

## 6. Energy and forces [Public] [Convention]

[Public] TN3 Sec. 4.3 and TN26: forces are computed by regarding the eigenstates of the local problems as exact Kohn–Sham states and inserting \\( \rho \\) and \\( E \\) into the same expressions as Cluster and Band
(Hellmann–Feynman terms and Pulay corrections; [Total energy and forces](energy-force.md)).
As long as the cluster is finite, these forces are not the exact derivatives of the total energy (not variational).
TN3 considers the error acceptable if geometry optimization is stopped at about \\( 10^{-4} \\) Hartree/bohr.

[Convention]

- Each term of the total energy uses the same expression as Cluster and Band, simply replacing \\( \rho \\) and \\( E \\) by the results of the O(N) method.
- Since \\( \rho \\) is asymmetric, sums over atom pairs (and over pairs of three-center integrals) are not halved by symmetry but taken over **all ordered pairs**.
  In Cluster and Band the \\( (J, K) \\) and \\( (K, J) \\) terms are combined with a factor of 2, but in DC, Krylov, and DC-LNO each ordered pair is added with factor 1
  (the factor 2 for non-spin-polarized remains).
- Whether stress (derivatives with respect to lattice vectors) can be used with O(N) methods is [Open].

## Parameters and defaults [Public] [Convention]

| Keyword | Default | Applies to | Meaning |
|---|---|---|---|
| `orderN.HoppingRanges` | 7.0 (Å) | DC, Krylov, DC-LNO | Radius \\( R\_{\mathrm{hop}} \\) of the truncated cluster |
| `orderN.FNAN+SNAN` | None | Same | Per-atom specification of the number of sites (lines of 3 integers) |
| `orderN.KrylovH.order` | 400 | Krylov | Target Krylov dimension \\( K\_H \\) for the Hamiltonian |
| `orderN.KrylovS.order` | \\( 4K\_H \\) | Krylov | Krylov dimension for \\( S \\) (with `orderN.Exact.Inverse.S off`) |
| `orderN.Exact.Inverse.S` | on | Krylov | Whether to build \\( (S^{(I)})^{-1} \\) exactly |
| `orderN.Recalc.Buffer` | on | Krylov | Whether to rebuild the embedding (buffer) part every iteration |
| `orderN.Expand.Core` | on | Krylov | Whether to include neighboring atoms in the starting vectors |
| `orderN.InitChemPot` | \\( -3.0 \\) (eV) | Krylov | Center of the energy window in the first iteration |
| `orderN.LNO.Occ.Cutoff` | 0.1 | DC-LNO | Occupation threshold for adopting LNOs (doubled for noncollinear) |
| `orderN.LNO.Buffer` | 0.2 | DC-LNO | Guide to the fraction of second neighbors using PAOs |
| `LNOs.Num` | None | DC-LNO (required for `Cluster-LNO`) | Number of LNOs per species |
| `scf.Npoles.ON2` | 100 | `Cluster2` | Number of poles |
| `scf.ElectronicTemperature` | 300 (K) | All | Temperature of the Fermi distribution |

- [Public] The manual does not give the default of `orderN.HoppingRanges`. It recommends 6.0–7.0 Å for molecular systems (DC), and for DC-LNO a radius containing 200–300 atoms.
  Comments in some distributed input examples say "default=5.0", but the default of the existing implementation is 7.0 Å [Convention].
- [Public] Manual Sec. 27.1: the error of the DC total energy is on the order of millihartree per atom; for gapped systems it decreases almost exponentially
  with cluster size, while for metals it converges slowly. The appendix of TN26 shows with a 1D model that the error can oscillate for metals and that
  raising the electronic temperature suppresses it.

## Edge cases [Convention]

- **Small unit cell**: multiple images of the same atom enter as separate sites, including images of the central atom. The result is handled correctly as a real-space problem distinguishing periodic images.
- **Isolated atom / atom without first neighbors**: \\( r\_{\min} \\) of Krylov is undefined, so the core is only the central atom.
- **Krylov dimension at least the number of basis functions of the cluster**: no Krylov subspace is built and all orbitals are used (equivalent to DC).
- **Nearly singular \\( S^{(I)} \\)**: DC and Krylov drop directions at \\( 10^{-9} \\), and the dropped directions become unoccupied states at 1000 Hartree.
  The exact inverse (Cholesky) of Krylov is not regularized; if not positive definite, an error is printed and the calculation continues. DC-LNO uses rule D and so can diverge.
- **Complex eigenvalues of \\( \Lambda\_I \\)**: the imaginary parts are discarded. Nothing special is done for the eigenvectors either ([Open]: effect on the result).
- **First iteration**: Krylov uses `orderN.InitChemPot` as the center of the energy window. DC-LNO builds LNOs from the input density matrix even in the first iteration
  (what that density matrix is in the first iteration without restart is [Open]. The distributed DC-LNO input examples specify `LNOs.Num`).
- **Unsupported combinations**: Krylov with noncollinear, and spin-orbit coupling with collinear spin (DC, Krylov, `Cluster2`) stop. `Cluster2` with noncollinear does nothing.

## Requirements for implementations [Design]

- Make **cluster construction** an independent module. Inputs are atomic positions, lattice, cutoff radii, and \\( R\_{\mathrm{hop}} \\); outputs are the sequence of sites (atom index and integer lattice-vector coordinates),
  the distinction between first and second neighbors, and the distance of each site. The order is defined as "center, first neighbors, second neighbors (ascending distance; ties broken lexicographically by atom index and lattice vector)".
  The default for specifications equivalent to `orderN.FNAN+SNAN` is to apply them at every MD step, with "first step only" selectable for comparison with the existing implementation.
- **Extraction of local matrices** is a function assembling principal submatrices from the [Sparse block format](../design/sparse-block-format.md); the structure (correspondence between site pairs and blocks of H(R)) is built only once per atomic configuration.
- Switch **local solvers** through a common interface (local problem → eigenvalues and residues for the rows of the central atom and the columns of the first neighbors).
  - `Dense` (DC): regularization specified by the enumeration of schemes of [Generalized eigenvalue problem](generalized-eigen.md), default `Drop { threshold: 1e-9 }`.
  - `Krylov`: make the core selection, number of blocks, inverse-\\( S \\) scheme (exact / Krylov approximation), and embedding update policy explicit settings.
    The timing of building the basis \\( U \\) (first iteration only / every iteration) is also a setting, defaulting to "first iteration only" as in the existing implementation.
  - `Lno`: how to determine the number of LNOs (threshold / fixed number per species), how to determine the PAO region (rule above), and regularization (default `Floor { threshold: 1e-10 }`, `Abs` for comparison).
- The **chemical potential** is a function taking the pairs \\( (\varepsilon, \kappa) \\) of all atoms (and the electron number of fixed occupations below the window); with MPI only per-atom partial sums are gathered.
  The bisection interval and stopping criterion are given by settings, defaulting to \\( [-30, 30] \\) Hartree and \\( 10^{-12} \\).
- The **energy window** is an optional optimization; by default no window is used and the Fermi distribution is applied to all states. For comparison with the existing implementation, the window (\\( \pm 10 \\) eV) is enabled.
- Represent the **asymmetry of the density matrix** by types. The output of O(N) methods is a set of blocks "not necessarily symmetric", holding \\( (I, J, \mathbf L) \\) and \\( (J, I, -\mathbf L) \\) independently.
  The force and energy modules have a path that sums over ordered pairs in this case. Symmetrization (\\( (\rho + \rho^{\mathsf T})/2 \\)) is an optional choice, not done by default.
- **Reporting**: record per atom the number of cluster sites and basis functions, the Krylov dimension and the number of dropped directions, the number of LNOs and the minimum adopted occupation, and the minimum eigenvalue of \\( S \\).

## Tests

### Analytic checks and invariants [Design]

Confirmed with `tools/check_on_methods.py` (1D ring, 2 orbitals per atom, non-orthogonal, up to second neighbors).

| Content | Tolerance |
|---|---|
| When cluster = whole system, \\( \rho \\) and \\( E \\) of DC agree with full diagonalization | \\( 10^{-10} \\) |
| \\( g\_s\sum\_I\sum\_{\mathcal F(I)}\rho S = N\_{\mathrm e} \\) with common \\( \mu \\) | \\( 10^{-9} \\) (electron number) |
| For a gapped system, the error of DC \\( \rho \\) decreases exponentially with \\( R\_{\mathrm{hop}} \\) (by a factor \\( 10^{-6} \\) or less for \\( R = 3 \to 11 \\) sites) | — |
| \\( \rho\_{IJ} \ne \rho\_{JI}^{\mathsf T} \\) for a finite cluster (the asymmetry actually occurs) | — |
| \\( S \\)-orthonormality of the Krylov basis \\( \lVert U^{\mathsf T} S U - 1\rVert\_{\max} \\) | \\( 10^{-11} \\) |
| The sum of \\( H\_K \\) split into core and buffer agrees with \\( U^{\mathsf T} H U \\) | \\( 10^{-12} \\) |
| Moments \\( W\_0^{\mathsf T} S A^p W\_0 \\) reproduced for \\( p \le 2q+1 \\), not for \\( p = 2q+2 \\) | \\( 10^{-10} \\) (relative) |
| With Krylov dimension \\( \ge n\_I \\), eigenvalues agree with the DC cluster | \\( 10^{-9} \\) |
| Krylov-approximate \\( S^{-1} \\) agrees with exact \\( S^{-1} \\) at full dimension (on the starting block) | \\( 10^{-9} \\) |
| \\( g\_s\sum\_I \mathrm{tr}\\,\Lambda\_I = N\_{\mathrm e} \\) | \\( 10^{-9} \\) |
| DC-LNO with all LNOs agrees with DC | \\( 10^{-10} \\) |
| PAO region rule: all second neighbors are PAO with \\( b = 1 \\); at least one is PAO even with \\( b = 0 \\) | Exact |

Additional invariants to place on the implementation side:

- The difference between the \\( \mathrm{tr} \\) sum rule \\( g\_s\sum \rho H \\) and \\( U\_{\mathrm{ele}} \\) approaches 0 as \\( R\_{\mathrm{hop}} \to \infty \\).
- In collinear spin-polarized calculations with equal \\( H \\) for both spins, the result agrees with the non-spin-polarized result (\\( 10^{-10} \\)).
- In noncollinear DC with zero spin-off-diagonal blocks, the result agrees with collinear DC (\\( 10^{-10} \\)).
- Under relabeling (permutation) of atoms, the per-atom blocks of \\( \rho \\) move by the same permutation (with the order of equidistant sites defined; rounding error).

### Comparison with the existing implementation [Design]

Fixtures built with the existing implementation (input, and total energy, chemical potential, per-atom Mulliken charges, forces).

| System | Solver | Comparison | Tolerance (guide) |
|---|---|---|---|
| `DIA8_DC.dat` of manual Sec. 27.1 (diamond 8 atoms, \\( R\_{\mathrm{hop}} = 6.0 \\) Å) | DC | `Utot`, `Chemical potential (Hartree)`, forces | \\( 10^{-6} \\) Hartree, \\( 10^{-5} \\) Hartree/bohr |
| Isolated molecule (C60 etc.) with \\( R\_{\mathrm{hop}} \\) larger than the molecular diameter | DC and Cluster | Total energies agree (cluster = whole system) | \\( 10^{-8} \\) Hartree |
| `Si8-LNO.dat` of manual Sec. 27.2 | DC-LNO | `Utot`, number of LNOs | \\( 10^{-6} \\) Hartree |
| Same system with `LNOs.Num` set to the number of PAOs | DC-LNO and DC | Total energies agree | \\( 10^{-8} \\) Hartree |
| Diamond of about 64 atoms with `orderN.KrylovH.order` at least the number of basis functions of the cluster | Krylov and DC | Total energies agree (`orderN.Recalc.Buffer on`) | \\( 10^{-8} \\) Hartree |
| Same system with \\( K\_H \\) varied over 20%–60% | Krylov | Curve of total energy vs. \\( K\_H \\) | \\( 10^{-6} \\) Hartree |
| `C60_LO.dat` of manual Chapter 48 (600 K, 90 poles) | `Cluster2` and Cluster | Total energies agree | \\( 10^{-9} \\) Hartree |
| Large-scale examples in the distribution (`-runtestL`): Krylov (diamond 512 atoms, `Utot` = −2964.574331104677) and DC-LNO (650 atoms, `Utot` = −4111.573022279569) | Krylov, DC-LNO | `Utot` | \\( 10^{-5} \\) Hartree |

- Since forces are not variational, only check that the difference from finite differences of the total energy decreases with \\( R\_{\mathrm{hop}} \\) (e.g. decreases monotonically for \\( R\_{\mathrm{hop}} = 6, 8, 10 \\) Å).
  For comparison with the forces of the existing implementation themselves, use \\( 10^{-5} \\) Hartree/bohr as a guide.
- For comparison, match the interval and stopping criterion for \\( \mu \\), the energy window, and the scope of `orderN.FNAN+SNAN` to the existing implementation.
- The reference values of the large-scale examples are outputs included in the distribution, and the version that produced them is unconfirmed. Recompute with the existing implementation of the same version before using them for comparison.

## Open issues [Open]

- The line format of `orderN.FNAN+SNAN` (discrepancy between the manual's 2 columns and the existing implementation's 3 integers), and whether applying it only at the first MD step is intended.
- The content of the density matrix used to build LNOs in the first SCF iteration of DC-LNO (no restart, no `LNOs.Num`), and the resulting number of LNOs.
- The handling of complex eigenvalues of \\( \Lambda\_I \\) (whether there are systems where it actually occurs, and its effect when it does).
- The order of equidistant sites when sorting second neighbors by distance (may affect the boundary of the PAO region in DC-LNO).
- Whether stress can be computed with O(N) methods, and its formula.
- The content of `EGAC` (no public documents), and the purpose of `orderN.NumPoles`.
- Details of the handling of the imaginary part of \\( \rho \\) (spin-orbit coupling etc.) in noncollinear DC and DC-LNO.
- The consulted version is the GPU version; whether the above thresholds, intervals, and defaults are the same in upstream OpenMX 3.9/4.0 has not been confirmed.
