# Electron transport by NEGF

> Status: first draft (awaiting review). The equations from technical note 13 (Eq. 6–56), technical note 4 (Eq. 1–18), technical note 20 (Eq. 7–72),
> T. Ozaki, PRB **75**, 035123 (2007), T. Ozaki, K. Nishio, H. Kino, PRB **81**, 035116 (2010), and Section 46 of the manual have been
> rewritten in the notation of this note; conventions and defaults not stated in the public documents have been supplemented as [Convention] from the behavior of the existing implementation.
> The numerical properties of the pole expansion, the surface Green's function, the self-energy, the transmission, and the nonequilibrium density matrix are checked with `tools/check_negf.py`.

## Scope of this chapter

This chapter covers the self-consistent solution, by the nonequilibrium Green's function (NEGF) method, of a central region sandwiched between left and right semi-infinite leads (`scf.EigenvalueSolver NEGF`),
and the post-processing that computes the transmission, current, conductance, eigenchannels, and current density from its result.

- The conventions for real-space matrices and Bloch sums follow [H(R)/S(R) sparse block format](../design/sparse-block-format.md) and [k-points and Bloch sums](kpoints-bloch.md).
- The density, Hartree potential, and XC on the grid presuppose [Real-space grid, Poisson, and exchange-correlation](grid-poisson-xc.md);
  this chapter describes only the NEGF-specific parts (supplementing the density from the leads, and the Poisson equation with boundary conditions).
- For the general theory of charge mixing and convergence criteria, see [Density matrix and charge mixing](density-mixing.md).
- The lead calculation (Step 1) is an ordinary `Band` calculation. This chapter describes only how its output is passed on.

## Notation [Design]

| Symbol | Meaning |
|---|---|
| \\( \mathrm L, \mathrm R \\) | Left lead, right lead. \\( s \in \\{\mathrm L, \mathrm R\\} \\) denotes a generic lead |
| \\( \mathrm L\_0, \mathrm C\_0, \mathrm R\_0 \\) | One cell of the left lead, the central part, one cell of the right lead. Extended central region \\( \mathrm C = \mathrm L\_0\\,\vert\\,\mathrm C\_0\\,\vert\\,\mathrm R\_0 \\) |
| \\( \mathrm L\_j, \mathrm R\_j \\) (\\( j \ge 1 \\)) | Lead cells continuing outside C. \\( \mathrm L\_1, \mathrm R\_1 \\) are adjacent to C |
| \\( \mathbf a\_1 \\) | Transport direction (a axis). \\( \mathbf a\_2, \mathbf a\_3 \\) are the in-plane (bc-plane) lattice vectors, common to the leads and C |
| \\( \mathbf a^{(s)}\_1 \\) | a-axis vector of the unit cell of lead \\( s \\) |
| \\( \mathbf q = (q\_2, q\_3) \\) | In-plane k-point (fractional coordinates with respect to \\( \tilde{\mathbf b}, \tilde{\mathbf c} \\)) |
| \\( N\_{\mathrm C}, N\_s \\) | Number of orbitals in C, number of orbitals in one cell of lead \\( s \\) (each doubled in the noncollinear case) |
| \\( H\_{\mathrm C}, S\_{\mathrm C} \\) | \\( N\_{\mathrm C}\times N\_{\mathrm C} \\) matrices of C (per k-point) |
| \\( H^{(s)}\_{00}, S^{(s)}\_{00} \\) | Blocks within one cell of lead \\( s \\) |
| \\( H^{(s)}\_{\mathrm{out}}, S^{(s)}\_{\mathrm{out}} \\) | Coupling blocks from a cell of lead \\( s \\) to the neighboring cell on the side away from C |
| \\( H\_{\mathrm Cs}, S\_{\mathrm Cs} \\) | Coupling blocks between C and \\( s\_1 \\) (\\( \mathrm L\_1 \\) or \\( \mathrm R\_1 \\)) (\\( N\_{\mathrm C}\times N\_s \\)) |
| \\( z \\) | Complex energy. On the real axis \\( z = E + i\eta \\) (\\( \eta > 0 \\)) |
| \\( g\_s(z) \\) | Surface Green's function of lead \\( s \\) (block of the \\( s\_1 \\) cell) |
| \\( \Sigma\_s(z), \Gamma\_s(z) \\) | Self-energy, broadening matrix |
| \\( G(z) \\) | Retarded Green's function of C |
| \\( \mu\_s \\) | Chemical potential of lead \\( s \\) (including bias) |
| \\( f(x) = 1/(1+e^{x}) \\), \\( \beta = 1/(k\_{\mathrm B}T) \\) | Fermi distribution and inverse temperature (\\( T \\) is `scf.ElectronicTemperature`) |
| \\( V\_{\mathrm b} \\) | Source-drain bias (`NEGF.bias.voltage`, in energy units) |

The Hermitian conjugate of a matrix \\( X \\) is written \\( X^\dagger \\), and the "imaginary part" of a matrix as \\( \mathrm{Im}\_{\mathrm m} X = (X - X^\dagger)/(2i) \\).
Energies are in Hartree internally; energies in the `NEGF.*` input are in eV (\\( 1\\,\mathrm{Hartree} = 27.2113845\\,\mathrm{eV} \\) is used [Convention]).

## 1. System setup

### 1.1 Three regions [Public]

According to technical note 13 (Section 2.1) and manual Section 46.1, the system is a one-dimensional chain arranged along the a axis as \\( \cdots\vert\mathrm L\_2\vert\mathrm L\_1\vert\mathrm C\vert\mathrm R\_1\vert\mathrm R\_2\vert\cdots \\),
periodic within the bc plane. For each in-plane Bloch wave vector \\( \mathbf q \\) the problem becomes an independent one-dimensional problem.
Since the localized basis is truncated at a finite radius, the matrices are block tridiagonal in the cells if the cells are chosen to satisfy the following two conditions.

1. The basis functions of \\( \mathrm C\_0 \\) may overlap those of \\( \mathrm L\_0, \mathrm R\_0 \\), but not those of \\( \mathrm L\_1, \mathrm R\_1 \\).
2. The basis functions of each lead cell do not overlap those of cells farther away than the neighboring cell.

\\( \mathrm L\_0, \mathrm R\_0 \\) are included in C in order to incorporate self-consistently the relaxation of the electronic states near the interface between the leads and the central part.

### 1.2 Construction of the extended central region [Public] [Convention]

[Public] According to manual Section 46.3 A, \\( \mathrm C\_0 \\) is given by `Atoms.Number` and `Atoms.SpeciesAndCoordinates`,
\\( \mathrm L\_0 \\) by `LeftLeadAtoms.Number` and `LeftLeadAtoms.SpeciesAndCoordinates`,
and \\( \mathrm R\_0 \\) by `RightLeadAtoms.Number` and `RightLeadAtoms.SpeciesAndCoordinates`.
The coordinate unit is restricted to `Ang` or `AU`; `FRAC` cannot be used. The unit cell of C is not input.

[Convention] The existing implementation assembles C as follows.

- **Atom order**: the atoms of C are numbered consecutively in the order \\( \mathrm L\_0 \\), \\( \mathrm C\_0 \\), \\( \mathrm R\_0 \\).
- **Matching with the leads**: the \\( i \\)-th atom of \\( \mathrm L\_0 \\) must coincide, within \\( 10^{-12} \\) bohr in each coordinate, with the \\( i \\)-th atom of the unit cell of the lead calculation (Step 1)
  translated so that the first atoms coincide
  (likewise for \\( \mathrm R\_0 \\)). Otherwise the program stops. That is, the number and order of atoms in \\( \mathrm L\_0 \\) are the same as in the lead unit cell.
- **Unit cell**: \\( \mathbf a\_2, \mathbf a\_3 \\) of C are taken as is from the left lead (the program stops if the b and c axes of the left and right leads differ by \\( 10^{-10} \\) or more).
  The a axis of C is the vector from the point \\( \mathbf O\_{\mathrm L} \\), obtained by translating the origin of the left-lead grid to match \\( \mathrm L\_0 \\),
  to the point obtained by adding \\( \mathbf a^{(\mathrm R)}\_1 \\) to the point \\( \mathbf O\_{\mathrm R} \\), obtained likewise by matching the origin of the right-lead grid to \\( \mathrm R\_0 \\):
  \\( \mathbf a\_1 = \mathbf O\_{\mathrm R} + \mathbf a^{(\mathrm R)}\_1 - \mathbf O\_{\mathrm L} \\). The origin of the grid of C is fixed at \\( \mathbf O\_{\mathrm L} \\).
  Hence the two end planes of the grid of C in the a direction coincide with the start of the left-lead cell and the end of the right-lead cell.
- **Number of grid points**: \\( N\_2, N\_3 \\) of C are forced to the lead values (the program stops if they differ between the left and right leads).
  \\( N\_1 \\) is determined by the usual rule ([Real-space grid](grid-poisson-xc.md) Section 1.2); there is no step aligning it with the a-direction grid spacing of the leads.
- **Region labels**: a per-atom label is printed to standard output. `1` is \\( \mathrm C\_0 \\), `2`/`3` are \\( \mathrm L\_0 \\)/\\( \mathrm R\_0 \\),
  and `12`/`13` are atoms of \\( \mathrm L\_0 \\)/\\( \mathrm R\_0 \\) whose basis functions overlap those of atoms of \\( \mathrm L\_1 \\)/\\( \mathrm R\_1 \\)
  (interatomic distance smaller than the sum of the cutoff radii of the two atoms).
- **Overlap check**: the program stops if the distance between an atom of \\( \mathrm C\_0 \\) and an atom of \\( \mathrm L\_1, \mathrm R\_1 \\) is smaller than the sum of the cutoff radii (condition 1 of Section 1.1).
  Condition 2 within the leads is checked likewise.

### 1.3 Requirements on the lead calculation (Step 1) [Public] [Convention]

[Public] Manual Section 46.2: Step 1 is an ordinary band calculation with `NEGF.output_hks on` and `NEGF.filename.hks` added.
The transport direction is the a axis; this is specified by permuting the order of the lattice vectors (in which case `scf.Kgrid` is permuted as well).
To make the semi-infiniteness consistent with Step 2, Step 1 uses a sufficiently large number of k-points along a, and the number of bc-plane k-points is the same in Step 1 and Step 2.

[Convention] When `NEGF.output_hks on`, the a axis of the lead must be parallel to one of the \\( x, y, z \\) axes and point in the positive direction
(\\( (a\_x, 0, 0) \\), \\( (0, a\_y, 0) \\), \\( (0, 0, a\_z) \\) with a positive component). Otherwise the program stops.
The gate voltage, the linear correction of the Poisson equation, and the mixing of \\( H \\) at the boundary described below use the Cartesian coordinate along this axis as the transport-direction coordinate \\( x \\).

## 2. Passing lead data

### 2.1 Contents of the hks file [Public] [Convention]

[Public] According to manual Section 46.2, the file `NEGF.filename.hks` contains the Hamiltonian matrix elements, the electron distribution, and the difference Hartree potential.
In Step 2, the files of the left and right leads are specified by `NEGF.filename.hks.l` and `NEGF.filename.hks.r` (they may be the same file).

[Convention] The existing implementation writes the following quantities of the lead after SCF convergence (binary; the byte layout need not be reproduced by a reimplementation, only the contents are specified).

| Contents | Use |
|---|---|
| Version, spin type, number of atoms, number of species, number of grid points \\( (N^{(s)}\_1, N\_2, N\_3) \\) | Consistency checks |
| Unit-cell vectors, reciprocal-lattice vectors, grid origin, atomic coordinates | Assembly of C (Section 1.2) |
| Chemical potential \\( \mu^{(0)}\_s \\) (Fermi level of Step 1) | Chemical potential of the lead |
| Number of orbitals per species, neighbor atoms and periodic-image indices per atom | Structure of the sparse block format |
| Overlap \\( S(\mathbf L) \\), \\( H(\mathbf L) \\) (per spin), and in the noncollinear case the imaginary part of \\( H \\) (spin-orbit etc.) | Surface Green's function, self-energy, boundary mixing |
| Real-space density matrix \\( \rho(\mathbf L) \\) (per spin) | Supplementing the density near the boundary (Section 5.5), supplementing Mulliken charges |
| Difference density on the grid and difference Hartree potential \\( \delta V\_{\mathrm H} \\) | Boundary conditions of the Poisson equation (Chapter 7) |

The program stops if the spin type (`scf.SpinPolarization`) differs between Step 1 and Step 2.

### 2.2 k-space blocks of the leads [Public] [Convention]

[Public] As in Eq. 3–4 of technical note 13, the per-k-point matrices are built by an in-plane Bloch sum.

[Convention] From the real-space blocks \\( O^{(s)}(n\_1, n\_2, n\_3) \\) of lead \\( s \\) (\\( \mathbf n \\) are the integer coordinates of the lattice vector in the lead unit cell),

$$
O^{(s)}\_{00}(\mathbf q) = \sum\_{n\_2, n\_3} e^{2\pi i (q\_2 n\_2 + q\_3 n\_3)}\\, O^{(s)}(0, n\_2, n\_3),\qquad
O^{(s)}\_{\mathrm{out}}(\mathbf q) = \sum\_{n\_2, n\_3} e^{2\pi i (q\_2 n\_2 + q\_3 n\_3)}\\, O^{(s)}(d\_s, n\_2, n\_3)
$$

(\\( O = H, S \\)). \\( d\_{\mathrm L} = -1 \\), \\( d\_{\mathrm R} = +1 \\); both are couplings to the neighboring cell "in the direction away from C".
As in [sparse block format](../design/sparse-block-format.md), the phase is attached to the lattice vector only. The coupling in the opposite direction is
\\( O^{(s)}\_{\mathrm{in}} = (O^{(s)}\_{\mathrm{out}})^\dagger \\).

The coupling \\( H\_{\mathrm Cs}, S\_{\mathrm Cs} \\) between C and \\( s\_1 \\) is also built from the lead data. Since \\( \mathrm L\_0 \\) is identical to a lead cell,
\\( H^{(\mathrm L)}\_{\mathrm{out}} \\) is placed as is in the \\( \mathrm L\_0 \\) rows of \\( H\_{\mathrm C\mathrm L} \\), and the \\( \mathrm C\_0, \mathrm R\_0 \\) rows are set to 0
(likewise on the right). Matrix elements computed on the C side are not used for this coupling.

### 2.3 Matrices of C and mixing at the boundary [Convention]

\\( H\_{\mathrm C}(\mathbf q), S\_{\mathrm C}(\mathbf q) \\) of C are Bloch sums over only those terms of the real-space matrices obtained from the SCF of C that **do not involve periodic images along a**
(\\( n\_1 = 0 \\)). The coupling to the neighboring images along a is replaced by the self-energy.

Furthermore, the blocks between atoms labeled `12` (atoms of \\( \mathrm L\_0 \\) that overlap \\( \mathrm L\_1 \\)) are replaced/mixed with the lead values as

$$
S\_{\mathrm C} \leftarrow S^{(\mathrm L)}\_{00},\qquad
H\_{\mathrm C} \leftarrow c\\, H\_{\mathrm C} + (1 - c)\\, H^{(\mathrm L)}\_{00}
$$

(likewise for `13` and the right lead). The weight \\( c \\) is the constant `NEGF.Mixing.H.Boundary` if it is 0 or larger;
by default (\\( -1 \\)) it is computed from the midpoint \\( \bar x \\) of the transport-direction coordinates of the atom pair and the edge \\( x\_{\mathrm e} \\) of C (the grid origin on the left, origin \\( + N\_1 h\_1 \\) on the right) as

$$
c(\bar x) = \tfrac12\left[\tanh\\!\left(\tfrac12\left(\lvert\bar x - x\_{\mathrm e}\rvert - 2\right)\right) + 1\right]\qquad(\bar x\\ \text{in bohr})
$$

\\( c = 1/2 \\) at 2 bohr from the edge, \\( c \approx 0.12 \\) at the edge, and \\( c \approx 0.998 \\) at 8 bohr from the edge.
The lead \\( H^{(s)}\_{00} \\) is the value after applying the bias (Section 6.1).

## 3. Surface Green's function [Public] [Convention]

[Public] Technical note 13 (Section 2.2, Eq. 9–21) obtains the surface Green's function \\( g\_s(z) \\) by the renormalization (decimation) method.
For lead \\( s \\), with initial values

$$
\varepsilon\_{\mathrm s}^{(0)} = \varepsilon^{(0)} = zS^{(s)}\_{00} - H^{(s)}\_{00},\qquad
\alpha^{(0)} = -\left(zS^{(s)}\_{\mathrm{out}} - H^{(s)}\_{\mathrm{out}}\right),\qquad
\beta^{(0)} = -\left(zS^{(s)}\_{\mathrm{in}} - H^{(s)}\_{\mathrm{in}}\right)
$$

and iterating, for \\( i = 0, 1, \dots \\),

$$
\begin{aligned}
\varepsilon\_{\mathrm s}^{(i+1)} &= \varepsilon\_{\mathrm s}^{(i)} - \alpha^{(i)} (\varepsilon^{(i)})^{-1} \beta^{(i)}, &
\varepsilon^{(i+1)} &= \varepsilon^{(i)} - \alpha^{(i)} (\varepsilon^{(i)})^{-1} \beta^{(i)} - \beta^{(i)} (\varepsilon^{(i)})^{-1} \alpha^{(i)},\\\\
\alpha^{(i+1)} &= \alpha^{(i)} (\varepsilon^{(i)})^{-1} \alpha^{(i)}, &
\beta^{(i+1)} &= \beta^{(i)} (\varepsilon^{(i)})^{-1} \beta^{(i)}
\end{aligned}
$$

one obtains \\( g\_s = \lim\_{i\to\infty} (\varepsilon\_{\mathrm s}^{(i)})^{-1} \\). The \\( i \\)-th renormalization corresponds to eliminating \\( 2^i \\) cells, and
\\( \lVert\alpha^{(i)}\rVert, \lVert\beta^{(i)}\rVert \\) approach 0 doubly exponentially if \\( \mathrm{Im}\\,z > 0 \\).
The resulting \\( g\_s \\) satisfies the Dyson equation

$$
g\_s = \left[zS^{(s)}\_{00} - H^{(s)}\_{00} - \left(zS^{(s)}\_{\mathrm{out}} - H^{(s)}\_{\mathrm{out}}\right) g\_s \left(zS^{(s)}\_{\mathrm{in}} - H^{(s)}\_{\mathrm{in}}\right)\right]^{-1}
$$

(residual \\( < 10^{-13} \\) confirmed with `tools/check_negf.py`).

[Convention]

- Convergence is reached when the maximum absolute element-wise difference between two successive \\( (\varepsilon\_{\mathrm s}^{(i)})^{-1} \\) falls below
  `NEGF.Surfgreen.convergeeps` (default \\( 10^{-12} \\)). The maximum number of iterations is `NEGF.Surfgreen.iterationmax` (default 600);
  reaching it does not stop the program: an error line is printed and the last value is used.
- At every stage of the density matrix and the transmission, \\( z \\) has either \\( \mathrm{Im}\\,z > 0 \\) or \\( \mathrm{Im}\\,z < 0 \\); it is never evaluated directly on the real axis.
  When the advanced side (\\( \mathrm{Im}\\,z < 0 \\)) is needed, the same renormalization is performed with \\( z^\* \\) (a second computation is done).
- Besides the renormalization method, an implementation of a method based on an eigenvalue problem (the method of Sanvito et al.) remains in the code but is not used.

## 4. Self-energy and the Green's function of the center [Public]

Rewriting Eq. 6–8 of technical note 13,

$$
\Sigma\_s(z) = \left(zS\_{\mathrm Cs} - H\_{\mathrm Cs}\right) g\_s(z) \left(zS\_{s\mathrm C} - H\_{s\mathrm C}\right),\qquad
S\_{s\mathrm C} = S\_{\mathrm Cs}^\dagger,\\ H\_{s\mathrm C} = H\_{\mathrm Cs}^\dagger
$$

$$
G(z) = \left[zS\_{\mathrm C} - H\_{\mathrm C} - \Sigma\_{\mathrm L}(z) - \Sigma\_{\mathrm R}(z)\right]^{-1}
$$

Because the basis is nonorthogonal, the self-energy depends explicitly on \\( z \\) through the coupling overlap \\( S\_{\mathrm Cs} \\).
\\( G(z) \\) is analytic in \\( z \\) in the upper and lower half-planes, and \\( G(z^\*) = G(z)^\dagger \\). The broadening matrix is

$$
\Gamma\_s(E) = i\left[\Sigma\_s(E + i\eta) - \Sigma\_s(E + i\eta)^\dagger\right]
$$

(technical note 13 Eq. 45, technical note 20 Eq. 7).

[Convention] \\( \Sigma\_s \\) and \\( G \\) are built as dense matrices, and \\( G \\) is obtained as an inverse by LU decomposition. In the collinear spin-polarized case this is done per spin;
in the noncollinear case in \\( 2N \\) dimensions (Chapter 12).

## 5. Density matrix

### 5.1 Equilibrium part: pole expansion of the Fermi function [Public]

Technical note 4 (Eq. 1–10) and PRB **75**, 035123 truncate the continued-fraction representation of the Fermi function at a finite depth \\( M = 2N\_{\mathrm p} \\) and decompose it into partial fractions.
The derivation is written here in the notation of this note. The truncated continued fraction can be written as

$$
f(x) \simeq \frac12 - \frac{x}{4}\left[(D + ixB)^{-1}\right]\_{11},\qquad
D = \mathrm{diag}(1, 3, 5, \dots, 2M - 1),\quad B\_{p,p\pm 1} = \tfrac12
$$

(\\( B \\) is tridiagonal with zero diagonal). Solving the generalized eigenvalue problem \\( Bv\_p = \lambda\_p D v\_p \\), \\( v\_p^{\mathsf T}Dv\_p = 1 \\), gives
\\( (D + ixB)^{-1} = \sum\_p v\_p v\_p^{\mathsf T}/(1 + ix\lambda\_p) \\); the eigenvalues appear in pairs \\( \pm\lambda \\), and \\( v\_{1p}^2 \\) are equal within a pair.
Combining the pairs,

$$
f(x) \simeq \frac12 + \sum\_{p=1}^{N\_{\mathrm p}} R\_p\left[\frac{1}{x - iz\_p} + \frac{1}{x + iz\_p}\right],\qquad
z\_p = \frac{1}{\lambda\_p}\\ (\lambda\_p > 0),\quad R\_p = -\frac{v\_{1p}^2\\, z\_p^2}{4}
$$

(technical note 4 Eq. 3). \\( z\_p \\) are real; the smallest roughly 60% coincide with the Matsubara poles \\( (2p - 1)\pi \\) and their \\( R\_p \\) with \\( -1 \\),
while the remaining poles rapidly become sparse and extend far out (`tools/check_negf.py`: already at \\( N\_{\mathrm p} = 20 \\) the error for \\( \lvert x\rvert \le 60 \\) is \\( 10^{-15} \\)).

With \\( x = \beta(E - \mu) \\) and \\( \alpha\_p = \mu + iz\_p/\beta \\), when the poles of \\( G(z) \\) lie on (or near) the real axis,

$$
\rho\_{\mathrm{eq}}(\mu) = -\frac{1}{\pi}\mathrm{Im}\_{\mathrm m}\\!\int\_{-\infty}^{\infty} dE\\, G(E + i0^+) f\\big(\beta(E - \mu)\big)
{}= \frac12\\, m^{(0)} - \frac{1}{\beta}\sum\_{p=1}^{N\_{\mathrm p}} R\_p\left[G(\alpha\_p) + G(\alpha\_p)^\dagger\right]
$$

(per spin). Since \\( R\_p < 0 \\), the weight \\( -R\_p/\beta \\) of each term is positive.
Technical note 4 Eq. 13 is written including the spin factor 2 (\\( \mu^{(0)} + \mathrm{Im}[-4i\beta^{-1}\sum G(\alpha\_p)R\_p] \\)), i.e., twice the above expression.
For an isolated system (no self-energy) this agrees with \\( \rho\_{\mathrm{eq}} = \sum\_n c\_n f(\varepsilon\_n) c\_n^\dagger \\) (to \\( 10^{-15} \\) with `tools/check_negf.py`).

### 5.2 Zeroth moment [Public] [Convention]

[Public] \\( m^{(0)} = -\pi^{-1}\mathrm{Im}\_{\mathrm m}\int G(E+i0^+)\\,dE \\) is the zeroth moment of the Green's function,
evaluated for large \\( R \\) as \\( m^{(0)} \simeq iR\\,G(iR) \\) (technical note 4 Eq. 14–18, technical note 13 Eq. 26).

[Convention] The existing implementation uses \\( R = 10^{12} \\) Hartree and computes \\( G(iR) \\) including the self-energy (with the surface Green's function also evaluated at \\( z = iR \\)).
When \\( S\_{\mathrm Cs} \ne 0 \\), \\( \Sigma\_s(z) \\) grows in proportion to \\( z \\), so \\( m^{(0)} \\) does not coincide with \\( S\_{\mathrm C}^{-1} \\)
(difference \\( 7\times 10^{-3} \\) in the model of `tools/check_negf.py`). A reimplementation must not substitute \\( S\_{\mathrm C}^{-1} \\) either.

[Convention] The number of points at which \\( G \\) is evaluated for the equilibrium part is \\( N\_{\mathrm p} + 1 \\) (\\( N\_{\mathrm p} \\) is `NEGF.Num.Poles`, default 150).

### 5.3 Finite bias: nonequilibrium part [Public] [Convention]

[Public] According to technical note 13 (Section 3.1, Eq. 42–46) and manual Section 46.3 B, the density matrix under bias is the sum of
an equilibrium part at the chemical potential of one lead and a nonequilibrium correction on the real axis. Since the nonequilibrium part is not analytic in \\( z \\),
it is integrated with an equally spaced rectangle rule on a line displaced from the real axis by a small imaginary part \\( \eta\_{\mathrm{neq}} \\).

[Convention] The existing implementation assembles it as follows. Let \\( \mu\_< = \min(\mu\_{\mathrm L}, \mu\_{\mathrm R}) \\), \\( \mu\_> = \max(\mu\_{\mathrm L}, \mu\_{\mathrm R}) \\),
and let \\( s\_> \\) be the lead with \\( \mu\_> \\).

$$
\rho = \rho\_{\mathrm{eq}}(\mu\_<) + \frac{\Delta E}{2\pi}\sum\_{i=0}^{n\_{\mathrm{neq}}-1}
\left[f\\big(\beta(E\_i - \mu\_>)\big) - f\\big(\beta(E\_i - \mu\_<)\big)\right] G(z\_i)\\,\Gamma\_{s\_>}(z\_i)\\,G(z\_i)^\dagger,\qquad z\_i = E\_i + i\eta\_{\mathrm{neq}}
$$

Here \\( \Gamma\_{s\_>}(z) = i[\Sigma\_{s\_>}(z) - \Sigma\_{s\_>}(z)^\dagger] \\),
\\( \Delta E \\) is `NEGF.bias.neq.energy.step` (default 0.02 eV), and \\( \eta\_{\mathrm{neq}} \\) is `NEGF.bias.neq.im.energy` (default 0.01 eV).
The poles of the equilibrium part are centered at \\( \mu\_< \\).

The integration interval \\( [E\_{\min}, E\_{\max}] \\) is symmetric about \\( \bar\mu = (\mu\_< + \mu\_>)/2 \\). \\( E\_{\max} \\) is
the first \\( E \\), on a grid going upward from \\( \bar\mu \\) in steps of \\( \Delta E/10 \\), at which \\( \lvert f(\beta(E - \mu\_>)) - f(\beta(E - \mu\_<))\rvert \\) falls below
`NEGF.bias.neq.cutoff` (default \\( 10^{-8} \\)), and \\( E\_{\min} = 2\bar\mu - E\_{\max} \\).
The number of points is \\( n\_{\mathrm{neq}} = \lfloor (E\_{\max} - E\_{\min})/\Delta E\rfloor \\), and the points are \\( E\_i = E\_{\min} + i\\,\Delta E \\) (\\( i = 0, \dots, n\_{\mathrm{neq}} - 1 \\)).

At \\( V\_{\mathrm b} = 0 \\) the difference is 0 already at the first point, so \\( n\_{\mathrm{neq}} = 0 \\) and only the equilibrium part remains.
The total number of points at which \\( G \\) is evaluated is \\( N\_{\mathrm p} + 1 + n\_{\mathrm{neq}} \\), and the computation time is proportional to it (manual Section 46.3).
Example: in the output example of manual Section 46.3 (carbon chain, \\( V\_{\mathrm b} = 0.5 \\) V), the interval is \\( E\_{\max} - \mu\_> = 0.954 \\) eV and \\( n\_{\mathrm{neq}} = 120 \\), corresponding to \\( 18.4\\,k\_{\mathrm B}T \\) at the electronic temperature of 600 K of the example input (\\( e^{-18.4} \approx 10^{-8} \\)).

[Design] The choice of reference lead is arbitrary; strictly, this equals \\( \rho\_{\mathrm{eq}}(\mu\_>) - \int G\Gamma\_{s\_<}G^\dagger\\,\Delta f\\,dE/2\pi \\)
(since \\( G(\Gamma\_{\mathrm L} + \Gamma\_{\mathrm R})G^\dagger = i(G - G^\dagger) \\), in the absence of bound states).
`tools/check_negf.py` confirms that the two agree within the errors due to \\( \eta\_{\mathrm{neq}} \\) and the step size (\\( 10^{-4} \\)).
\\( \eta\_{\mathrm{neq}} \\) is also a regularization that captures the contribution of bound states within the interval as Lorentzians of width \\( \eta \\); unless \\( \Delta E \lesssim \eta\_{\mathrm{neq}} \\) or so, the integration becomes coarse.

### 5.4 k-points and the real-space density matrix [Public] [Convention]

[Public] Technical note 13 Eq. 22, 43 builds the real-space density matrix by integrating over the in-plane Brillouin zone (that note writes the phase as \\( e^{-i\mathbf k\cdot\mathbf R} \\)).
Manual Section 46.3: the number of in-plane k-points is `NEGF.scf.Kgrid` (default `1 1`); none is specified along a.

[Convention]

- The k-points are \\( q\_j = -\tfrac12 + \tfrac{2l + 1}{2N\_j} \\) (\\( j = 2, 3 \\), \\( l = 0, \dots, N\_j - 1 \\)), the same formula as `Regular` in [k-points](kpoints-bloch.md)
  (with small shifts of \\( +10^{-12} \\) added to \\( q\_2 \\) and \\( -10^{-12} \\) to \\( q\_3 \\)). Monkhorst–Pack reduction (the `scf.MP.criterion` scheme) cannot be used with NEGF; specifying it stops the program.
- In collinear calculations, \\( (l\_2, l\_3) \\) and \\( (N\_2 - 1 - l\_2, N\_3 - 1 - l\_3) \\) are treated as a time-reversal pair and only one of them is computed with weight 2; a point that is its own partner has weight 1.
  No reduction is made in noncollinear calculations.
- The real-space density matrix is obtained from the per-k-point \\( \rho(\mathbf q) \\) (the matrices of Sections 5.1–5.3) as
  \\( \rho\_{I\alpha,J\beta}(\mathbf n) = \frac{1}{N\_2N\_3}\sum\_{\mathbf q} w\_{\mathbf q}\\,\mathrm{Re}\left[\rho\_{I\alpha,J\beta}(\mathbf q)\\,e^{-2\pi i(q\_2n\_2 + q\_3n\_3)}\right] \\)
  (the imaginary part is kept separately with the same phase). Since \\( \rho(\mathbf q)\_{I\alpha,J\beta} \\) has the form \\( c\_{I\alpha}\overline{c\_{J\beta}} \\),
  this is the same convention as \\( e^{+i\mathbf k\cdot\mathbf L}\\,\overline{c\_{I\alpha}}c\_{J\beta} \\) in [sparse block format](../design/sparse-block-format.md).
- Density-matrix elements between atoms of C and their periodic images along a (\\( n\_1 = \pm 1 \\)) are set to 0.
- No search for the chemical potential (electron-number constraint) is performed in the NEGF stage. The number of electrons in C is a result determined by the chemical potentials of the leads.

### 5.5 Density on the grid and supplementation from the leads [Public] [Convention]

[Public] According to technical note 13 (Section 2.4, Eq. 27–31), the density on the grid of C is the sum of terms between basis functions of C, cross terms between basis functions of C and of the leads, and terms between basis functions of the leads.

[Convention] The existing implementation builds it as follows.

- The contribution from the density matrix of C (Section 5.4) is the same as the usual density on the grid ([Real-space grid](grid-poisson-xc.md) Section 2.1),
  except that contributions from the periodic images of orbitals along a (\\( n\_1 \ne 0 \\)) are not included. The same applies to the superposition of atomic densities (the reference of \\( \delta n \\)).
- In slab-shaped regions at the two boundaries, the contribution of the tails of the basis functions of the atoms of \\( \mathrm L\_1 \\) (\\( \mathrm R\_1 \\)) is added using **the density matrix of the Step 1 lead**:
  the terms between orbitals of \\( \mathrm L\_1 \\), and twice the cross terms between orbitals of \\( \mathrm L\_1 \\) and \\( \mathrm L\_0 \\).
  In-plane, images up to \\( -1, 0, +1 \\) along \\( \mathbf a\_2, \mathbf a\_3 \\) are summed.
  The boundary region on the left extends from grid plane 0 along a to "the plane of the rightmost edge of the cutoff spheres of the atoms of \\( \mathrm L\_0 \\), shifted by \\( -\mathbf a^{(\mathrm L)}\_1 \\), \\( + 1 \\)"; the right is defined symmetrically.
  The contributions of the atoms of \\( \mathrm L\_1, \mathrm R\_1 \\) are likewise added to the superposition of atomic densities in the same regions.
- To the Mulliken charges of the orbitals of the atoms of \\( \mathrm L\_0 \\) (\\( \mathrm R\_0 \\)), the term
  \\( \sum\_\beta S\_{\alpha\beta}\rho\_{\alpha\beta} \\) with the orbitals of \\( \mathrm L\_1 \\) (\\( \mathrm R\_1 \\)) is added, using \\( S \\) and \\( \rho \\) of the Step 1 lead.

### 5.6 Old integration contour (`NEGF.Integration OLD`) [Convention]

The default (`CF`) is the pole expansion of Section 5.1. Specifying `OLD` uses a contour consisting of the Matsubara poles within imaginary part \\( \Delta = 1.2\\,\mathrm{eV} \\) from \\( \mu \\),
a line segment at imaginary part \\( \Delta \\) (\\( \mu \pm 40k\_{\mathrm B}T \\), 400-point Gauss–Legendre), and an arc returning to the lower end \\( -1 \\) Hartree (40 points).
It can be used only at zero bias; specifying it with a finite bias stops the program. The contour in the appendix of technical note 13 (\\( \Delta = 1.2 \\) eV, \\( \gamma = 40k\_{\mathrm B}T \\), 40-point arc) corresponds to this.
It is out of scope for a reimplementation [Design].

## 6. Bias and gate

### 6.1 Bias voltage [Public] [Convention]

[Public] Technical note 13 Eq. 47–49: the bias on the right lead is applied by \\( H\_{\mathrm R} \to H\_{\mathrm R} + V\_{\mathrm b}S\_{\mathrm R} \\), \\( \mu\_{\mathrm R} \to \mu\_{\mathrm R} + V\_{\mathrm b} \\),
and by shifting the boundary difference Hartree potential by the same amount. A uniform shift of the electrostatic potential of a lead is equivalent, in a localized basis, to adding \\( V S \\).
Manual Section 46.3: `NEGF.bias.voltage` is given in eV, and a chemical-potential difference of 1 eV is read as 1 V.

[Convention] The existing implementation splits the voltage equally between left and right.

$$
V\_{\mathrm L} = -\tfrac12 V\_{\mathrm b},\quad V\_{\mathrm R} = +\tfrac12 V\_{\mathrm b};\qquad
H^{(s)}(\mathbf L) \to H^{(s)}(\mathbf L) + V\_s S^{(s)}(\mathbf L),\quad
\mu\_s = \mu^{(0)}\_s + V\_s,\quad
\delta V^{(s)}\_{\mathrm H}(\mathbf r) \to \delta V^{(s)}\_{\mathrm H}(\mathbf r) + V\_s
$$

(added as an electron potential energy; a positive bias raises the levels of the right lead). In the collinear case it is added to all spin components;
in the noncollinear case only to the \\( \uparrow\uparrow, \downarrow\downarrow \\) components (\\( V\_s\\,I\_2\otimes S \\)).
It is added to all blocks of the lead (\\( H\_{00} \\), \\( H\_{\mathrm{out}} \\), the coupling to C).

> **Note on manual Section 46.3**: the output example in that section corresponds to adding "0 on the left and \\( V\_{\mathrm b} \\) on the right", whereas the implementation consulted uses \\( \pm V\_{\mathrm b}/2 \\).
> \\( \mu\_{\mathrm R} - \mu\_{\mathrm L} \\) and the width of the nonequilibrium window are the same; only the overall energy origin shifts by \\( V\_{\mathrm b}/2 \\).
> However, the reference of the electrostatic potential of C (the boundary values of Chapter 7) changes accordingly, so the splitting must be matched when comparing results.

Setting `NEGF.bias.apply` (default `on`) to `off` skips the entire bias treatment (this path has a question listed under [Open] at the end).

### 6.2 Gate voltage [Public] [Convention]

[Public] Technical note 13 Eq. 50, manual Section 46.3 D:

$$
V\_{\mathrm g}(x) = V\_{\mathrm g}^{(0)}\exp\\!\left[-\left(\frac{x - x\_{\mathrm c}}{d}\right)^{8}\right]
$$

is added. \\( V\_{\mathrm g}^{(0)} \\) is `NEGF.gate.voltage` (eV, default 0).

[Convention] \\( x \\) is the transport-direction Cartesian coordinate of a grid plane (plane index \\( \times\\,\lvert\mathbf h\_1\rvert \\) + origin),
\\( x\_{\mathrm c} \\) is the mean \\( x \\) coordinate of the atoms of \\( \mathrm C\_0 \\), and \\( d = \lvert a\_{1,x}\rvert - \lvert a^{(\mathrm L)}\_{1,x}\rvert - \lvert a^{(\mathrm R)}\_{1,x}\rvert \\)
(the length of the \\( \mathrm C\_0 \\) part). \\( V\_{\mathrm g} \\) is added in every SCF iteration to the difference Hartree potential (as an electron potential energy) after solving the Poisson equation.
Hence \\( V\_{\mathrm g}^{(0)} > 0 \\) pushes the levels of \\( \mathrm C\_0 \\) up.

### 6.3 Non-SCF bias (`NEGF.NonSCF.bias`) [Convention]

When `on` (default `off`), a linear function evaluated at the midpoint \\( \bar x \\) of an atom pair,
\\( V(\bar x) = V\_{\mathrm L} + (V\_{\mathrm R} - V\_{\mathrm L})(\bar x - x\_0)/(x\_1 - x\_0) \\), is added to \\( H\_{\mathrm C} \\) of C as \\( V(\bar x)\\,S\_{\mathrm C} \\) (\\( x\_0, x\_1 \\) are the two ends of the grid of C).

## 7. Boundary conditions of the Poisson equation

### 7.1 Two-dimensional FFT + one-dimensional finite differences (`FD`, default) [Public] [Convention]

[Public] Technical note 13 (Section 2.5, Eq. 32–37): the Poisson equation for the difference density \\( \delta n \\) is Fourier-transformed in the bc plane
and solved by finite differences along a. With the in-plane reciprocal-lattice vectors \\( \mathbf G\_\parallel \\),

$$
\left(\frac{d^2}{dx^2} - \lvert\mathbf G\_\parallel\rvert^2\right)\delta\tilde V\_{\mathrm H}(x, \mathbf G\_\parallel) = -4\pi\\,\delta\tilde n(x, \mathbf G\_\parallel)
$$

with Dirichlet conditions fixing the values at both ends to \\( \delta V\_{\mathrm H} \\) of Step 1 of the leads. Since the electrostatic potential at the interfaces is continuous with the leads,
this is theoretically the more appropriate choice (manual Section 46.3).

[Convention]

- With grid planes \\( k = 0, \dots, N\_1 - 1 \\) and spacing \\( h\_1 = \lvert\mathbf h\_1\rvert \\),
  \\( V\_{k-1} - (2 + h\_1^2\lvert\mathbf G\_\parallel\rvert^2)V\_k + V\_{k+1} = -4\pi h\_1^2\\,\delta\tilde n\_k \\)
  is solved for each \\( \mathbf G\_\parallel \\) as a tridiagonal linear system (\\( h\_1 \\) is the length of the a-axis grid vector; this implicitly assumes that the a axis is orthogonal to the bc plane).
- The boundary value \\( V\_{-1} \\) is the bc-plane Fourier transform of the **last** a plane of \\( \delta V\_{\mathrm H} \\) of the left lead (including the bias \\( V\_{\mathrm L} \\)), and \\( V\_{N\_1} \\) that of the **first** a plane of the right lead.
  Since the two ends of the grid of C coincide with lead cell boundaries (Section 1.2), these are the planes adjacent to C on the outside.
- `NEGF.Poisson_Gparazero.scaling` (default 1.0) is a factor multiplying \\( \delta\tilde n \\) of the \\( \mathbf G\_\parallel = 0 \\) component.
- After solving, the gate voltage of Section 6.2 is added.
- Charge mixing uses an NEGF-specific variant performed in this representation (bc-plane Fourier components and a-direction grid planes) (details under [Open]).

### 7.2 Three-dimensional FFT (`FFT`) [Public] [Convention]

[Public] Manual Section 46.3: solving by a three-dimensional FFT (periodic boundaries) can also be selected. For nonpolar systems it gives nearly the same results as FD, and SCF tends to converge more easily.

[Convention] To the usual three-dimensional FFT solution (with the \\( \mathbf G = 0 \\) component set to 0), a linear function along a, \\( v(x) = v\_0 + (v\_1 - v\_0)(x - x\_0)/(x\_1 - x\_0) \\), is added.
\\( v\_0 \\) (\\( v\_1 \\)) is the plane average of the first a plane of \\( \delta V\_{\mathrm H} \\) of the left (right) lead minus the plane average of the first (last) a plane of C in the FFT solution, and
\\( x\_0, x\_1 \\) are the coordinates of the two ends of the grid of C. Only the plane averages (\\( \mathbf G\_\parallel = 0 \\) components) match at the boundaries.

In addition, there is an undocumented scheme `FDG` (two-dimensional finite differences with two boundary values for a gate, giving two values to `NEGF.gate.voltage`) [Convention]. It is out of scope for a reimplementation [Design].

### 7.3 Hartree potential when switching from Band to NEGF [Convention]

When `NEGF.Band2NEGF.Hartree on` (default `off`), and only for FD, in the first 5 iterations after switching to NEGF
the solved \\( \delta\tilde V\_{\mathrm H} \\) is mixed with the linear interpolation between the two boundary values (the straight line connecting planes \\( -1 \\) and \\( N\_1 \\)).
The weights of the solution in iterations 1–5 are 0.05, 0.15, 0.25, 0.50, 0.75, respectively.

## 8. SCF flow [Public] [Convention]

[Public] According to manual Sections 46.1 and 46.3, (i) the SCF of the leads (Step 1), then (ii) the equilibrium Green's function of C, at finite bias (iii) the nonequilibrium Green's function,
and (iv) density update and mixing are repeated until convergence, after which the transmission etc. are computed (Step 3).
Running the first few iterations with ordinary diagonalization (Band), assuming periodicity along a as well, speeds up the final convergence (`NEGF.scf.Iter.Band`).
The NEGF convergence criterion applies `scf.criterion` to the residual norm `NormRD` (`dUele` for other solvers).

[Convention]

- While the SCF iteration number is at most `NEGF.scf.Iter.Band`, the density matrix is built by a Band calculation (only the Γ point along a, `NEGF.scf.Kgrid` in-plane);
  thereafter it is built by NEGF (Chapter 5). The default is **10** (different from 6 in the manual). When restarting from a density-matrix file, NEGF is used from the start.
- The chemical potential in the Band stage is fixed at \\( (\mu\_{\mathrm L} + \mu\_{\mathrm R})/2 \\) for `NEGF.ChemP.Band fix` (default), and determined from the electron count for `free`.
- The Poisson equation is solved with the scheme of Chapter 7 from the Band stage onward.
- Convergence is reached when the square root of the residual norm \\( \sqrt{\lvert\mathrm{NormRD}\rvert} \\) is below `scf.criterion`.
- In the NEGF stage the band-energy term (the quantity corresponding to \\( E\_{\mathrm{band}} \\)) is set to 0, and the energy-weighted density matrix is not computed (Section 9.5).

## 9. Transmission, current, conductance (Step 3)

### 9.1 Transmission [Public] [Convention]

[Public] Technical note 13 Eq. 51–52 (Caroli formula):

$$
T\_\sigma(E; \mathbf q) = \mathrm{Tr}\left[\Gamma\_{\mathrm L}(E)\\,G(E + i\eta)\\,\Gamma\_{\mathrm R}(E)\\,G(E - i\eta)\right],\qquad
T\_\sigma(E) = \frac{1}{N\_{\mathbf q}}\sum\_{\mathbf q} T\_\sigma(E; \mathbf q)
$$

The k-points are given by `NEGF.tran.Kgrid` (default: the value of `NEGF.scf.Kgrid`); using more than in the SCF is recommended (manual Section 46.4.1).

[Convention]

- The energy grid is \\( E\_i = \mu\_{\mathrm L} + E\_{\mathrm{lo}} + (E\_{\mathrm{hi}} - E\_{\mathrm{lo}})\\,i/(n\_{\mathrm{div}} - 1) \\) (\\( i = 0, \dots, n\_{\mathrm{div}} - 1 \\)),
  where the three values of `NEGF.tran.energyrange` are \\( E\_{\mathrm{lo}}, E\_{\mathrm{hi}}, \eta \\) (default \\( -10, 10, 10^{-3} \\) eV), and \\( n\_{\mathrm{div}} \\) is `NEGF.tran.energydiv` (default 200).
  The energy origin is \\( \mu\_{\mathrm L} \\) including the bias.
- \\( G(E - i\eta) \\) and \\( \Sigma\_s(E - i\eta) \\) are recomputed at \\( z^\* \\), and \\( \Gamma\_s = i[\Sigma\_s(E + i\eta) - \Sigma\_s(E - i\eta)] \\) (mathematically the same as \\( \Sigma^\dagger \\)). The real part of the trace is output.
- The k-point grid and the time-reversal reduction are the same as in Section 5.4; the value at a reduced point is copied to its partner point.
- A finite \\( \eta \\) acts as absorption in C, so even for a perfect chain \\( 1 - T \\) grows in proportion to \\( \eta \\)
  (\\( 1 - T = 4\times 10^{-4} \\) for \\( \eta = 10^{-5} \\) Hartree in the 5-cell chain of `tools/check_negf.py`). Comparisons at the default \\( 10^{-3} \\) eV (\\( 3.7\times 10^{-5} \\) Hartree) must allow for this systematic error.

### 9.2 Current [Public] [Convention]

[Public] Technical note 13 Eq. 53 (Landauer–Büttiker):

$$
I\_\sigma = \frac{e}{h}\int dE\\, T\_\sigma(E)\left[f\\big(\beta(E - \mu\_{\mathrm R})\big) - f\\big(\beta(E - \mu\_{\mathrm L})\big)\right]
$$

(per spin; the sign is positive for \\( \mu\_{\mathrm R} > \mu\_{\mathrm L} \\)).

[Convention]

- In atomic units (\\( e = \hbar = 1 \\), \\( e/h = 1/2\pi \\)), \\( I\_\sigma = -\frac{\Delta E}{2\pi}\sum\_i [f\_{\mathrm L}(E\_i) - f\_{\mathrm R}(E\_i)]\\,T\_\sigma(E\_i + i\eta\_{\mathrm I}) \\) is computed
  and multiplied by \\( eE\_{\mathrm h}/\hbar = 6.6236178\times 10^{-3} \\) A to convert to amperes (CODATA 2018 gives \\( 6.6236182\times 10^{-3} \\) A; relative difference \\( 7\times 10^{-8} \\)).
- The integration interval and points follow the same rule as Section 5.3, with step \\( \Delta E \\) = `NEGF.current.energy.step` (default 0.01 eV), cutoff `NEGF.current.cutoff` (default \\( 10^{-8} \\)),
  and imaginary part \\( \eta\_{\mathrm I} \\) = `NEGF.current.im.energy` (default: the third value of `NEGF.tran.energyrange`). If the number of points is less than 10, a warning that the step is too coarse is issued.
- The temperature is `scf.ElectronicTemperature` (default 300 K).
- The current is computed per k-point, and the per-k-point values and their average are written to `System.Name.current`. Without spin polarization the same value (the per-spin value) is written in both spin columns, so the total current is twice that.

### 9.3 Conductance [Public] [Convention]

[Public] Manual Section 46.4.1: in the low-bias limit, \\( G = (e^2/h)\\,T(\mu\_{\mathrm L}) \\) is output per spin in units of \\( G\_0 = e^2/h \\).

[Convention] The implementation obtains the transmission at \\( \bar\mu = (\mu\_{\mathrm L} + \mu\_{\mathrm R})/2 \\) by linear interpolation between two adjacent points of the transmission energy grid
(identical to the manual's \\( \mu\_{\mathrm L} \\) at zero bias). The per-k-point values and their average are written to `System.Name.conductance`.

### 9.4 Interpolation of bias effects [Public] [Convention]

[Public] Technical note 13 Eq. 54–56, manual Section 46.7: after finishing SCF at two biases, \\( H\_{\mathrm C} \\), \\( H\_{\mathrm R} \\), \\( \mu\_{\mathrm R} \\) are
linearly interpolated with \\( \lambda \\) and \\( 1 - \lambda \\) to obtain the transmission and current (`NEGF.tran.interpolate`, `NEGF.tran.interpolate.file1`, `file2`, `NEGF.tran.interpolate.coes`).

[Convention] The implementation interpolates all of \\( H \\) of C, \\( H \\) of both the left and right leads, the left and right chemical potentials, and the chemical potential of C with the same coefficients
(because \\( V\_{\mathrm b} \\) is split between left and right, the left lead also changes).

### 9.5 Forces, stress, total energy [Convention]

In the NEGF stage, the energy-weighted density matrix is not computed and the band-energy term is set to 0. The automatic tests of manual Section 46.11 compare total energies and forces,
but the meaning of the forces and total energy obtained by NEGF (variationality in an open system) is not described in the public documents. This chapter leaves them out of scope and lists them under [Open].

## 10. Eigenchannels [Public] [Convention]

[Public] Technical note 20 (Eq. 8–33) applies the method of Paulsson–Brandbyge (PRB **76**, 115117 (2007)) to a nonorthogonal basis.
Under Löwdin orthogonalization \\( S\_{\mathrm C} = S^{1/2}S^{1/2\dagger} \\),

$$
\bar A\_{\mathrm L} = S^{1/2\dagger}\\, G\Gamma\_{\mathrm L}G^\dagger\\, S^{1/2},\qquad
\bar\Gamma\_{\mathrm R} = S^{-1/2}\\,\Gamma\_{\mathrm R}\\,S^{-1/2\dagger},\qquad
\hat T = \bar A\_{\mathrm L}\bar\Gamma\_{\mathrm R}
$$

From \\( \bar A\_{\mathrm L} = \sum\_n a\_n\\, u\_n u\_n^\dagger \\), \\( \bar A\_{\mathrm L}^{1/2} = (\sqrt{a\_1}u\_1, \sqrt{a\_2}u\_2, \dots) \\) is formed, and
the Hermitian matrix \\( \bar A\_{\mathrm L}^{1/2\dagger}\bar\Gamma\_{\mathrm R}\bar A\_{\mathrm L}^{1/2} \\) is diagonalized to obtain the transmission eigenvalues \\( t\_n \\) and the orthogonal \\( \tilde t\_n \\).
The eigenchannels are plotted in real space using \\( t\_n' = S^{-1/2\dagger}\bar A\_{\mathrm L}^{1/2}\tilde t\_n \\) as coefficients in the original basis.
\\( \sum\_n t\_n = T \\) and \\( 0 \le t\_n \le 1 \\) (confirmed with `tools/check_negf.py`).

[Public] Manual Section 46.4.3: `NEGF.tran.Channel` (default `on`); k-points via `NEGF.Channel.Nkpoint` and `<NEGF.Channel.kpoint` (in-plane fractional coordinates);
energies via `NEGF.Channel.Nenergy` and `<NEGF.Channel.energy` (eV relative to \\( \mu\_{\mathrm L} \\)); the number plotted in real space via `NEGF.Channel.Num` (5 collinear, 10 noncollinear).
Transmission eigenvalues are output in descending order.

[Convention]

- \\( S^{1/2} \\) is built from the eigendecomposition of \\( S\_{\mathrm C} \\), removing directions with eigenvalues below \\( 10^{-6} \\) (`N_{ort} / N_{nonort}` in standard output gives the retained and original dimensions).
  Channels in the removed directions are filled in with eigenvalue 0.
- The imaginary part is the same \\( \eta \\) as for the transmission (the third value of `NEGF.tran.energyrange`).
- The overall phase of each channel is chosen so that the sum of coefficients (with orientation aligned by the sign of the imaginary part) is real. The phase is a convention that does not affect physical quantities; comparisons are made with \\( \lvert t'\_n\rvert \\) or the density \\( \lvert T\_n(\mathbf r)\rvert^2 \\) [Design].

## 11. Current density [Public] [Convention]

[Public] Technical note 20 (Eq. 34–58) splits the current density into a local and a nonlocal term.
The local term, using \\( D\_{\mathrm R} = G\Gamma\_{\mathrm R}G^\dagger \\) (a matrix in the nonorthogonal basis), is

$$
\mathbf J\_{\mathrm{Loc}}(\mathbf r) = \int\frac{dE}{2\pi}\left[f(E - \mu\_{\mathrm L}) - f(E - \mu\_{\mathrm R})\right] i \sum\_{\mathbf q}\sum\_{ij} D^{\mathbf q}\_{ij}
\left\\{2i\mathbf k\\,\chi\_j^\*\chi\_i + \chi\_j^\*\nabla\chi\_i - \chi\_i\nabla\chi\_j^\*\right\\}
$$

(technical note 20 Eq. 36–40, 64), and the nonlocal term due to the nonlocal pseudopotential is obtained as the solution of the Poisson equation \\( \nabla^2\varphi = -\rho\_{\mathrm{NL}} \\),
\\( \mathbf J\_{\mathrm{NL}} = -\nabla\varphi \\), determined by charge conservation (Eq. 41–47). At the boundaries, Neumann conditions are given from the Landauer–Büttiker current (Eq. 48–55).
The charge current is \\( \mathbf J\_\uparrow + \mathbf J\_\downarrow \\) and the spin current is \\( \mathbf J\_\uparrow - \mathbf J\_\downarrow \\) (Eq. 57–58).
In the noncollinear case they are defined for each 2×2 spin component (Eq. 63–72).

[Public] Manual Section 46.4.2: `NEGF.tran.CurrentDensity` (default `on`); the noncollinear spin off-diagonal components via `NEGF.OffDiagonalCurrent` (default `off`).

[Convention]

- The energy integration uses the same grid and imaginary part as the current in Section 9.2. It is not computed when \\( \lvert\mu\_{\mathrm L} - \mu\_{\mathrm R}\rvert < 10^{-6} \\) Hartree.
- \\( S^{-1} \\) is the pseudo-inverse from the eigendecomposition of \\( S\_{\mathrm C} \\), removing directions with eigenvalues below \\( 10^{-6} \\).
- Averages over the Voronoi polyhedron of each atom are written to the `_atom` files.

The detailed discretization (evaluation of gradients, implementation of the Neumann conditions) is not described in this version ([Open]).

## 12. Noncollinear NEGF [Public] [Convention]

[Public] Manual Section 46.9: NEGF can be combined with noncollinear DFT and used together with spin-orbit coupling, DFT+U, and constraints on the directions of spin and orbital magnetic moments.
The procedure is the same as in the collinear case.

[Convention]

- The matrices of C and the leads are \\( 2N \\)-dimensional with spin as the outer index (the noncollinear spinor ordering of [k-points](kpoints-bloch.md)); the overlap is \\( I\_2\otimes S \\).
  The \\( \downarrow\uparrow \\) block is filled with the Hermitian conjugate of \\( \uparrow\downarrow \\).
- No time-reversal reduction of k-points is made (in both the SCF and Step 3).
- The bias adds \\( V\_sS \\) only to \\( \uparrow\uparrow \\) and \\( \downarrow\downarrow \\) (Section 6.1).
- The transmission file does not separate spins and has a single column for the total transmission (real and imaginary parts). The default of `NEGF.Channel.Num` is 10.
- In supplementing the density near the boundary (Section 5.5), the imaginary part of the \\( \uparrow\downarrow \\) component of the density matrix is added with its sign reversed.

## 13. Other features [Public] [Convention]

- **Recomputing Step 3 only** (manual Section 46.5): with `NEGF.tran.SCF.skip on`, `System.Name.tranb` is read and only Step 3 is performed.
  To stop after Step 2 only, set `NEGF.tran.SCF.skip off`, `NEGF.tran.Analysis off`, `NEGF.tran.Channel off`.
  [Convention] If any of `NEGF.tran.Channel`, `NEGF.tran.CurrentDensity`, `NEGF.OffDiagonalCurrent` is `on`, `NEGF.tran.Analysis` is automatically set to `on`.
  `NEGF.tran.SCF.skip` also accepts the value `Periodic` (for the three-dimensionally periodic systems of the next item).
- **Zero-bias periodic systems** (manual Section 46.6): adding `NEGF.Output.for.TranMain on` to Step 1 creates `*.tranb` from \\( H \\) of an ordinary band calculation
  of a system that is periodic also along a, so that the transmission can be computed without Step 2.
- **Density of states** (manual Section 46.3 E): when `Dos.fileout on`, the density of states is built from \\( -\pi^{-1}\mathrm{Im}\\,G \\) with `NEGF.Dos.energyrange` (lower bound, upper bound, imaginary part; default \\( -10, 10, 5\times 10^{-3} \\) eV),
  `NEGF.Dos.energy.div` (default 200), and `NEGF.Dos.Kgrid` (default `1 1`), and written to
  `System.Name.Dos.val` and `System.Name.Dos.vec`. The energy origin is \\( \mu\_{\mathrm L} \\).
- **Parallelization** (manual Section 46.8): combinations of energy points, k-points, and spins are distributed with MPI, and matrix operations are parallelized with OpenMP. This does not affect the results.
- **Incompatible features** [Convention]: cannot be used together with basis optimization (`orbitalOpt.Method`); the program stops.

## 14. Output files [Public] [Convention]

| File | Contents | Source |
|---|---|---|
| Specified by `NEGF.filename.hks` (Step 1) | Lead data (Section 2.1) | Manual 46.2 |
| `System.Name.tranb` (Step 2) | Matrices of C and the leads, chemical potentials, region labels, etc. needed for Step 3 | Manual 46.1, 46.5 |
| `System.Name.paranegf` | Chemical potentials of the leads (eigenvalue-based and after bias), net bias, interval, step, and number of points of the nonequilibrium integration. Merged into `System.Name.out` [Convention] | [Convention] |
| `System.Name.tran#_%` | Transmission per k-point \\( (\\#, \\%) \\). Columns: index, \\( E - \mu\_{\mathrm L} \\) and imaginary part (Hartree), same (eV), real and imaginary parts of \\( T\_\uparrow \\), real and imaginary parts of \\( T\_\downarrow \\) | Manual 46.4.1, [Convention] |
| `System.Name.current` | Current per k-point and average (A, per spin) | Manual 46.4.1 |
| `System.Name.conductance` | Conductance per k-point and average (in units of \\( G\_0 = e^2/h \\), per spin) | Manual 46.4.1 |
| `System.Name.traneval#k_#E_#s`, `tranevec#k_#E_#s` | Transmission eigenvalues, LCAO coefficients of the eigenchannels | Manual 46.4.3 |
| `System.Name.tranec#k_#E_#s_#c_r.cube`/`_i.cube` | Real and imaginary parts of the eigenchannels (Gaussian cube) | Manual 46.4.3 |
| `System.Name.curden.xsf`, `scurden.xsf` | Charge and spin current densities (vector fields on the grid) | Manual 46.4.2 |
| `System.Name.curden1.cube`, `scurden1.cube` | Same, component along a | Manual 46.4.2 |
| `System.Name.curden_atom.xsf`, `scurden_atom.xsf` | Same, averaged over the Voronoi polyhedra of the atoms | Manual 46.4.2 |
| `System.Name.odcurden*_r`/`_i` | Noncollinear spin off-diagonal components (real and imaginary parts) | Manual 46.4.2 |

\\( \\# \\) and \\( \\% \\) are the indices of the in-plane k-point (0-based).

## 15. Input keywords and defaults (summary)

[Public] marks those described in manual Section 46; [Convention] marks those not described or whose default differs. Keywords are case-insensitive [Convention].

| Keyword | Default | Meaning | Category |
|---|---|---|---|
| `scf.EigenvalueSolver` | — | `NEGF` for Step 2 | Public |
| `NEGF.output_hks`, `NEGF.filename.hks` | `off`, `NEGF.hks` | Write lead data in Step 1 | Public |
| `NEGF.filename.hks.l`, `.r` | `NEGF.hks.l`, `NEGF.hks.r` | Left and right lead data | Public |
| `LeftLeadAtoms.*`, `RightLeadAtoms.*` | — | Atoms of \\( \mathrm L\_0 \\), \\( \mathrm R\_0 \\) | Public |
| `NEGF.Num.Poles` | 150 | Number of poles \\( N\_{\mathrm p} \\) | Public |
| `NEGF.scf.Kgrid` | `1 1` | In-plane k-points for SCF | Public |
| `NEGF.scf.Iter.Band` | 10 | Number of initial iterations run with Band (6 in the manual) | Convention |
| `NEGF.ChemP.Band` | `fix` | Chemical potential in the Band stage | Convention |
| `NEGF.Band2NEGF.Hartree` | `off` | Hartree mixing at switching (Section 7.3) | Convention |
| `NEGF.bias.apply` | `on` | Whether to apply bias treatment | Convention |
| `NEGF.bias.voltage` | 0.0 eV | \\( V\_{\mathrm b} \\) | Public |
| `NEGF.bias.neq.im.energy` | 0.01 eV | \\( \eta\_{\mathrm{neq}} \\) (stops if negative) | Public |
| `NEGF.bias.neq.energy.step` | 0.02 eV | \\( \Delta E \\) (stops if negative) | Public |
| `NEGF.bias.neq.cutoff` | \\( 10^{-8} \\) | Cutoff of the nonequilibrium window | Convention |
| `NEGF.NonSCF.bias` | `off` | Section 6.3 | Convention |
| `NEGF.gate.voltage` | 0.0 eV | \\( V\_{\mathrm g}^{(0)} \\) | Public |
| `NEGF.Poisson.Solver` | `FD` | `FD`/`FFT` (/`FDG`) | Public (`FDG` is Convention) |
| `NEGF.Poisson_Gparazero.scaling` | 1.0 | Section 7.1 | Convention |
| `NEGF.Mixing.H.Boundary` | \\( -1 \\) (tanh weight) | Section 2.3 | Convention |
| `NEGF.Integration` | `CF` | `CF`/`OLD` | Convention |
| `NEGF.Surfgreen.iterationmax` | 600 | Maximum number of renormalization steps | Convention |
| `NEGF.Surfgreen.convergeeps` | \\( 10^{-12} \\) | Convergence criterion of the renormalization | Convention |
| `NEGF.tran.Analysis` | `on` | Transmission and current in Step 3 | Public |
| `NEGF.tran.SCF.skip` | `off` | Perform Step 3 only | Public |
| `NEGF.tran.energyrange` | \\( -10, 10, 10^{-3} \\) eV | Range and imaginary part for the transmission | Public |
| `NEGF.tran.energydiv` | 200 | Number of transmission points | Public |
| `NEGF.tran.Kgrid` | value of `NEGF.scf.Kgrid` | In-plane k-points for Step 3 | Public (default is Convention) |
| `NEGF.current.energy.step` | 0.01 eV | Step of the current integration | Convention |
| `NEGF.current.im.energy` | imaginary part of `NEGF.tran.energyrange` | Imaginary part of the current integration | Convention |
| `NEGF.current.cutoff` | \\( 10^{-8} \\) | Cutoff of the current window | Convention |
| `NEGF.tran.Channel` | `on` | Eigenchannels | Public |
| `NEGF.Channel.Nkpoint`, `<NEGF.Channel.kpoint` | 1, \\( (0, 0) \\) | k-points for eigenchannels | Public |
| `NEGF.Channel.Nenergy`, `<NEGF.Channel.energy` | 1, 0 eV | Energies for eigenchannels | Public |
| `NEGF.Channel.Num` | 5 (10 noncollinear) | Number of channels plotted in real space | Public |
| `NEGF.tran.CurrentDensity` | `on` | Current density | Public |
| `NEGF.OffDiagonalCurrent` | `off` | Noncollinear spin off-diagonal current | Public |
| `NEGF.tran.interpolate` (`.file1`, `.file2`, `.coes`) | `off` (coefficients `1.0 0.0`) | Bias interpolation | Public |
| `NEGF.Output.for.TranMain` | `off` | Create tranb from a three-dimensionally periodic system | Public |
| `NEGF.Dos.energyrange`, `.energy.div`, `.Kgrid` | \\( -10, 10, 5\times 10^{-3} \\) eV, 200, `1 1` | NEGF density of states | Public |
| `NEGF.FFTE.Num.Buffer.Cells` | 1 | Apparently unused ([Open]) | Convention |

## 16. Edge cases and caveats [Convention] [Design]

- **Zero bias**: there are 0 nonequilibrium points and only the equilibrium part remains (Section 5.3). The current is 0, and the current density is not computed.
- **Sign of the bias**: for \\( V\_{\mathrm b} < 0 \\), \\( \mu\_<, \mu\_> \\) are swapped and the same formula is used. The reference lead is always the one with the lower chemical potential.
- **Non-convergence of the surface Green's function**: the calculation continues even if it does not converge within 600 iterations. This can happen very close to band edges when \\( \eta \\) is small.
  A reimplementation counts non-convergence as a warning and reports it [Design].
- **\\( \mathrm C\_0 \\) identical to the leads**: in a perfect crystal, the transmission equals the number of channels (an integer) at each in-plane k-point (see Tests).
- **Coarse interval step**: if \\( E\_{\max} - E\_{\min} < \Delta E \\), a warning is issued and the nonequilibrium part has 0 points.
- **Grid spacing along a**: since \\( N\_1 \\) of C is determined independently of the lead grid spacing, the FD boundary values (values on lead grid planes) and the grid-plane spacing of C do not match exactly.

## Requirements for implementations [Design]

- **Module decomposition**: split into (1) a type for lead data (the equivalent of hks: real-space \\( H, S, \rho \\), \\( \delta V\_{\mathrm H} \\) on the grid, density, chemical potential),
  (2) assembly of the extended central region (atom correspondence, region labels, boundary mixing), (3) generation of per-k-point blocks (\\( H\_{00}, H\_{\mathrm{out}}, H\_{\mathrm Cs}, H\_{\mathrm C} \\)),
  (4) the surface Green's function, (5) the self-energy and \\( G \\), (6) generation of integration points and weights (pole expansion, moment, nonequilibrium rectangle rule),
  (7) assembly of the density matrix (sum over k-points), (8) a Poisson solver with boundary conditions, and (9) Step 3 analysis (transmission, current, eigenchannels, current density).
- **Representation of integration points**: generate integration points as a sequence of "complex energy \\( z \\), weight \\( w \\), kind (pole, moment, nonequilibrium)",
  and accumulate the density matrix as \\( \sum w\\,X(z) \\) (\\( X \\) is one of \\( G + G^\dagger \\), \\( G \\), \\( G\Gamma G^\dagger \\)).
  Since each point is independent, parallelization over combinations of points, k-points, and spins is possible.
- **Pole expansion**: make the function returning \\( z\_p, R\_p \\) from \\( N\_{\mathrm p} \\) independent so that it can be tested with the formula of Section 5.1.
- **Bias splitting**: for comparison with the existing implementation, use \\( (V\_{\mathrm L}, V\_{\mathrm R}) = (-V\_{\mathrm b}/2, +V\_{\mathrm b}/2) \\) as the default, changeable through an argument.
- **Surface Green's function**: use the renormalization method as the default, and return the convergence threshold, maximum number of iterations, actual number of iterations, and final difference. Make it possible to compute the Dyson residual for testing.
- **Reuse of the self-energy**: \\( g\_s(z) \\) depends only on the lead, \\( z \\), and \\( \mathbf q \\), not on the SCF of C, so it can be cached across SCF iterations
  (the existing implementation recomputes it every time; the results are unchanged).
- **Zeroth moment**: compute \\( iR\\,G(iR) \\) including the self-energy. \\( R \\) is an argument (default \\( 10^{12} \\) Hartree).
- **k-points**: use the [k-points](kpoints-bloch.md) functions in two dimensions for the in-plane `Regular` grid and time-reversal reduction.
- **Units**: convert input eV to Hartree for storage, and make \\( eE\_{\mathrm h}/\hbar \\) an explicit constant for the output current (allowing a choice between the existing implementation's value and the CODATA value).

## Tests [Design]

`tools/check_negf.py` confirms many of the following on a one-dimensional model with a two-orbital nonorthogonal lead cell.

| Kind | Contents | Approximate tolerance |
|---|---|---|
| Analytic | Maximum error of the pole expansion (\\( N\_{\mathrm p} = 150 \\)) for \\( \lvert x\rvert \le 60 \\) | \\( 10^{-12} \\) |
| Analytic | The first 20 \\( z\_p \\) equal \\( (2p-1)\pi \\), \\( R\_p = -1 \\) | \\( 10^{-10} \\) (relative) |
| Analytic | Surface Green's function of a one-orbital chain (\\( \varepsilon\_0, t \\)) \\( g = [w - \sqrt{w - 2t}\sqrt{w + 2t}]/(2t^2) \\), \\( w = z - \varepsilon\_0 \\) | \\( 10^{-7} \\) |
| Invariant | The renormalization result satisfies the Dyson equation of the surface Green's function | \\( 10^{-9} \\) |
| Invariant | For a perfect chain (\\( \mathrm C\_0 \\) = lead), \\( T(E) \\) = number of channels counted from the lead bands (\\( \eta = 10^{-8} \\)) | \\( 10^{-5} \\) |
| Invariant | \\( T\_{\mathrm{LR}} = \mathrm{Tr}[\Gamma\_{\mathrm L}G\Gamma\_{\mathrm R}G^\dagger] = T\_{\mathrm{RL}} \\) | \\( 10^{-10} \\) |
| Invariant | Sum of transmission eigenvalues = \\( T \\), \\( 0 \le t\_n \le 1 \\) | \\( 10^{-10} \\) |
| Invariant | For an isolated system, the pole-expansion density matrix = \\( \sum\_n c\_n f(\varepsilon\_n) c\_n^\dagger \\) | \\( 10^{-12} \\) |
| Invariant | For a system with leads, the density matrices for \\( N\_{\mathrm p} = 150 \\) and 300 agree, and the density matrix is Hermitian | \\( 10^{-11} \\) |
| Comparison (independent computation) | Difference between the pole-expansion density matrix and a real-axis integral (\\( \eta = 2\times 10^{-3} \\)) | \\( 5\times 10^{-3} \\) (\\( O(\eta) \\)) |
| Invariant | The nonequilibrium density matrix is independent of the reference lead (\\( \mu\_< \\) vs. \\( \mu\_> \\)) | \\( 2\times 10^{-3} \\) (depends on \\( \eta\_{\mathrm{neq}} \\) and the step) |
| Analytic | \\( eE\_{\mathrm h}/\hbar = 6.6236\times 10^{-3} \\) A | \\( 10^{-6} \\) (relative) |
| Invariant | \\( T(E; \mathbf q) = T(E; -\mathbf q) \\) at time-reversal-paired k-points (collinear) | \\( 10^{-10} \\) |
| Invariant | \\( T \\) is invariant for the left-right swapped (mirrored) system, and the sign of the current reverses under \\( V\_{\mathrm b} \to -V\_{\mathrm b} \\) | \\( 10^{-8} \\) |
| Invariant | Adding \\( VS \\) to all lead blocks and shifting \\( \mu \\) by \\( V \\) as well gives \\( T(E + V) \\) equal to the original \\( T(E) \\) (with C shifted by the same amount) | \\( 10^{-10} \\) |
| Comparison | Mulliken charges after convergence of `NormRD` in Step 2 for the manual examples (carbon chain Lead-Chain.dat/NEGF-Chain.dat, graphene, 8-ZGNR at finite bias 0.3 V, noncollinear 8-ZGNR, noncollinear gold chain) | \\( 10^{-6} \\) |
| Comparison | Transmission curves of the same examples (columns 6 and 8 of `System.Name.tran0_0`) | \\( 10^{-6} \\) (matching the existing implementation's \\( \eta \\)) |
| Comparison | Current, conductance, and sum of transmission eigenvalues of the same examples | \\( 10^{-7} \\) (the automatic tests of manual Section 46.11 compare 7 decimal places) |
| Comparison | Reproduce the nonequilibrium integration interval and number of points in the existing implementation's `System.Name.paranegf` with the rule of Section 5.3 | Number of points exact, interval to \\( \Delta E/10 \\) |

Fixtures for comparison are created by running Step 1 and Step 2 of each example with the existing implementation and saving hks (lead data), tranb, standard output, and the Step 3 output files.
Since results may change slightly with the degree of parallelism, the number of MPI processes is recorded.

## [Open]

- The `NEGF.bias.apply off` path: in this path the existing implementation uses the pole weight \\( -2R\_p/\beta \\) (twice that of the default path), while the moment term remains the same \\( \tfrac12 iR\\,G(iR) \\).
  This appears inconsistent with the formula of Section 5.1. Whether it is intentional (combined with a factor elsewhere) or a bug needs to be checked by running it.
- The meaning and computation of forces, stress, and total energy in NEGF (the energy-weighted density matrix is not computed; the appendix of technical note 13 gives formulas for the nonequilibrium energy density matrix).
- Details of the NEGF-specific charge mixing (RMM-DIIS and Kerker in the representation of bc-plane Fourier components and a-direction grid planes).
- Implementation details of the discretization of the current density (gradients, the Poisson equation for the nonlocal term and its Neumann boundary conditions).
- Boundary-current diagnostics in standard output. The manual's output example shows `Left(ideal)`/`Right(ideal)` and `Left(truncated)`/`Right(truncated)`,
  which would be surface integrals of the Landauer–Büttiker boundary current and integrals of the grid current density over the boundary planes,
  their difference indicating the error due to the discretization of the current density and the basis truncation. The implementation consulted
  (GPU version, commit 3871b1f) does not print these labels (only plain left/right values), so which quantities it prints, and whether the upstream
  CPU version prints the labelled ones, must be confirmed by running both.
- The meaning of the `FDG` scheme, the two-valued `NEGF.gate.voltage`, and `NEGF.FFTE.Num.Buffer.Cells`.
- Error estimate of the FD boundary condition when \\( N\_1 \\) of C and the a-direction grid spacing of the leads do not match.
- Details of the file formats of current and conductance in noncollinear NEGF (treatment of spin columns), and the conventions for the output of `NEGF.OffDiagonalCurrent`.
- The byte layouts of `System.Name.tranb` and the hks file (whether to build a compatibility layer that reads the existing implementation's output).
- The implementation consulted was the GPU version of OpenMX 4.0; it has not been checked whether the defaults (especially `NEGF.scf.Iter.Band`) are the same in the upstream OpenMX CPU version.
