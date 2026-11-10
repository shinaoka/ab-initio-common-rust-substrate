# Wannier Functions, Band Unfolding, and Localized-Orbital Analysis

> Status: first version (awaiting review). The formulas for MLWF, CWF, unfolding and orbital optimization are the public documents rewritten in the notation of this note;
> defaults, band selection, normalizations, etc. were confirmed as behavior of the existing implementation.
> Properties of the formulas are checked numerically with `tools/check_wannier_unfolding.py`.

## Scope of this chapter

Among the post-processing and analysis steps that take as input the LCAO eigenstates obtained by the SCF ([k points and Bloch sums](kpoints-bloch.md), [generalized eigenvalue problem](generalized-eigen.md)) and
the real-space matrices \\( H(\mathbf L) \\), \\( S(\mathbf L) \\) ([sparse block format](../design/sparse-block-format.md)),
this chapter covers the following five.

| Section | Feature | Main input keywords | Public documents |
|---|---|---|---|
| 1 | Maximally localized Wannier functions (MLWF) | `Wannier.*`, `Wannier90.fileout` | manual 47.2, Marzari–Vanderbilt PRB **56**, 12847, Souza–Marzari–Vanderbilt PRB **65**, 035109, Weng–Ozaki–Terakura PRB **79**, 235118 |
| 2 | Closest Wannier functions (CWF) | `CWF.*` | manual 47.1, technical note CWF (Ozaki PRB **110**, 125115) |
| 3 | Band unfolding | `Unfolding.*` | manual 55, Lee–Yamada-Takamura–Ozaki JPCM **25**, 345501 (arXiv:1211.4690), technical note unfolding_otw (Lee et al. PRB **102**, 075143) |
| 4 | Natural atomic orbitals (NAO), natural population analysis, and O(N) orbital analysis | `NBO.*`, `NAO.*` | manual 36, Ohwaki–Otani–Ozaki JCP **140**, 244105 |
| 5 | Orbital optimization of the basis (generation of contracted PAOs) | `orbitalOpt.*`, `CntOrb.fileout` | manual 26, technical note 5 (Ozaki PRB **67**, 155108) |

Here "manual" refers to section numbers of the OpenMX 4.0 manual (https://www.openmx-square.org/openmx_man4.0/index.html).
The SCF itself, the output format of band dispersion (`Band.dispersion`), and the cube file format are not covered.

## Notation [Design]

The notation of [Notation, units, and terminology](../../../conventions/notation.md) and [k points and Bloch sums](kpoints-bloch.md) is used. Additions are as follows.

| Symbol | Meaning |
|---|---|
| \\( N\_k \\) | number of k points (equal to the number of cells of the Born–von Kármán period) |
| \\( c\_{I\alpha,\mu}(\mathbf k) \\) | LCAO coefficients; \\( C(\mathbf k)^\dagger S(\mathbf k) C(\mathbf k) = I \\) |
| \\( \varepsilon\_{\mu\mathbf k} \\) | Kohn–Sham eigenvalues; \\( \mu\_{\mathrm{F}} \\) is the chemical potential ("Fermi level" in the output) |
| \\( \psi\_{\mu\mathbf k} \\) | \\( \psi\_{\mu\mathbf k}(\mathbf r) = \sum\_{\mathbf L} e^{i\mathbf k\cdot\mathbf L}\sum\_{I\alpha} c\_{I\alpha,\mu}(\mathbf k)\\, \phi\_{I\alpha}(\mathbf r - \mathbf R\_I - \mathbf L) \\) |
| \\( w\_{n\mathbf L} \\) | Wannier function (the \\( n \\)-th in cell \\( \mathbf L \\)); \\( N\_W \\) is the number per cell |
| \\( U(\mathbf k) \\) | transformation from the band space to Wannier functions (\\( N\_{\mathrm{band}}\times N\_W \\), \\( U^\dagger U = I \\)) |

\\( \psi\_{\mu\mathbf k} \\) is normalized per cell: \\( \langle\psi\_{\mu\mathbf k}\vert\psi\_{\nu\mathbf k'}\rangle = N\_k\delta\_{\mathbf k\mathbf k'}\delta\_{\mu\nu} \\).
The Bloch-sum phase is attached **only to lattice vectors** (atomic positions are not included; [k points and Bloch sums](kpoints-bloch.md)).
With this convention one can choose \\( c(\mathbf k + \mathbf G) = c(\mathbf k) \\), so shifting a k point by a reciprocal lattice vector does not require relabeling the coefficients.

Wannier functions and the Hamiltonian matrix elements (hoppings) in terms of them are, commonly for MLWF and CWF,

$$
w\_{n\mathbf L}(\mathbf r) = \frac{1}{N\_k}\sum\_{\mathbf k} e^{-i\mathbf k\cdot\mathbf L}\sum\_{\mu}\psi\_{\mu\mathbf k}(\mathbf r)\\, U\_{\mu n}(\mathbf k),
\qquad
t\_{mn}(\mathbf L) \equiv \langle w\_{m\mathbf 0}\vert\hat H\vert w\_{n\mathbf L}\rangle
{}= \frac{1}{N\_k}\sum\_{\mathbf k} e^{-i\mathbf k\cdot\mathbf L}\left[U(\mathbf k)^\dagger E(\mathbf k) U(\mathbf k)\right]\_{mn}
$$

(\\( E(\mathbf k) = \mathrm{diag}(\varepsilon\_{\mu\mathbf k}) \\)). Wannier interpolation is \\( H^{W}(\mathbf k) = \sum\_{\mathbf L} e^{i\mathbf k\cdot\mathbf L}\\, t(\mathbf L) \\),
with the same sign as the SCF Bloch sum.

---

## 1. Maximally localized Wannier functions (MLWF)

### 1.1 Spread functional [Public]

The spread of Marzari–Vanderbilt (PRB **56**, 12847)

$$
\Omega = \sum\_{n=1}^{N\_W}\left[\langle r^2\rangle\_n - \lvert\langle\mathbf r\rangle\_n\rvert^2\right]
$$

is evaluated by finite differences on the k-point grid. At each k point, the difference vectors \\( \mathbf b \\) to neighboring k points and weights \\( w\_b \\) are chosen to satisfy

$$
\sum\_{\mathbf b} w\_b\\, b\_\alpha b\_\beta = \delta\_{\alpha\beta}
$$

and from the overlap matrices

$$
M\_{mn}(\mathbf k,\mathbf b) = \langle u\_{m\mathbf k}\vert u\_{n,\mathbf k+\mathbf b}\rangle
{}= \frac{1}{N\_k}\langle\psi\_{m\mathbf k}\vert e^{-i\mathbf b\cdot\mathbf r}\vert\psi\_{n,\mathbf k+\mathbf b}\rangle
$$

(\\( u = e^{-i\mathbf k\cdot\mathbf r}\psi \\), represented after multiplication by \\( U(\mathbf k) \\)),

$$
\bar{\mathbf r}\_n = -\frac{1}{N\_k}\sum\_{\mathbf k,\mathbf b} w\_b\\, \mathbf b\\, \mathrm{Im}\ln M\_{nn}(\mathbf k,\mathbf b),
$$

$$
\Omega\_{\mathrm I} = \frac{1}{N\_k}\sum\_{\mathbf k,\mathbf b} w\_b\left(N\_W - \sum\_{mn}\lvert M\_{mn}\rvert^2\right),\quad
\Omega\_{\mathrm{OD}} = \frac{1}{N\_k}\sum\_{\mathbf k,\mathbf b} w\_b\sum\_{m\ne n}\lvert M\_{mn}\rvert^2,\quad
\Omega\_{\mathrm D} = \frac{1}{N\_k}\sum\_{\mathbf k,\mathbf b} w\_b\sum\_n\left(-\mathrm{Im}\ln M\_{nn} - \mathbf b\cdot\bar{\mathbf r}\_n\right)^2
$$

and \\( \Omega = \Omega\_{\mathrm I} + \Omega\_{\mathrm D} + \Omega\_{\mathrm{OD}} \\) (corresponding to Eqs. 31–36 of the same paper).
\\( \Omega\_{\mathrm I} \\) does not depend on the choice of \\( U(\mathbf k) \\) (the gauge); it is determined by the band subspace alone.

The gauge-dependent part is minimized by updating \\( U(\mathbf k) \leftarrow U(\mathbf k)\exp\Delta W(\mathbf k) \\) with an anti-Hermitian matrix \\( \Delta W(\mathbf k) \\). The gradient is

$$
G(\mathbf k) = 4\sum\_{\mathbf b} w\_b\left(\mathcal A[R] - \mathcal S[T]\right),\qquad
R\_{mn} = M\_{mn}M\_{nn}^{\*},\quad T\_{mn} = \frac{M\_{mn}}{M\_{nn}}\left(\mathrm{Im}\ln M\_{nn} + \mathbf b\cdot\bar{\mathbf r}\_n\right)
$$

(\\( \mathcal A[X] = (X - X^\dagger)/2 \\), \\( \mathcal S[X] = (X + X^\dagger)/2i \\)), and for steepest descent
\\( \Delta W = \dfrac{\epsilon}{4 w\_{\mathrm{tot}}}\\, G \\), \\( w\_{\mathrm{tot}} = \sum\_{\mathbf b} w\_b \\) (Sec. IV of the same paper).

### 1.2 Disentanglement [Public]

Following Souza–Marzari–Vanderbilt (PRB **65**, 035109), at each k point an \\( N\_W \\)-dimensional subspace is selected from the bands within the outer energy window,
and \\( \Omega\_{\mathrm I} \\) is minimized. States within the inner window are always included in the subspace (frozen).
At each iteration, the \\( N\_W \\) vectors corresponding to the largest eigenvalues of the matrix \\( Z(\mathbf k) \\), built from projections onto the subspaces of neighboring k points, are selected,
and oscillations are suppressed by linear mixing with \\( Z \\) of the previous iteration (mixing coefficient `Wannier.Dis.Mixing.Para`).
The convergence criterion applies `Wannier.Dis.Conv.Criterion` to the difference of \\( \Omega\_{\mathrm I} \\) between consecutive iterations (Å²) (manual 47.2.3).
Afterwards, the gauge-dependent part is minimized as in Sec. 1.1.

### 1.3 Matrix elements in LCAO [Public] [Convention]

[Public] Weng–Ozaki–Terakura (PRB **79**, 235118) write \\( M \\) and the initial projection \\( A \\) needed to construct MLWFs in terms of LCAO coefficients and integrals between atomic orbitals.
In the convention of this note,

$$
M\_{mn}(\mathbf k,\mathbf b) = \sum\_{I\alpha,J\beta}\sum\_{\mathbf L} c\_{I\alpha,m}^{\*}(\mathbf k)\\, c\_{J\beta,n}(\mathbf k + \mathbf b)\\,
e^{i(\mathbf k+\mathbf b)\cdot\mathbf L}\\, e^{-i\mathbf b\cdot\mathbf R\_I}\\,
\left\langle\phi\_{I\alpha}(\mathbf r - \mathbf R\_I)\middle\vert e^{-i\mathbf b\cdot(\mathbf r - \mathbf R\_I)}\middle\vert\phi\_{J\beta}(\mathbf r - \mathbf R\_J - \mathbf L)\right\rangle
$$

(even if \\( \mathbf k + \mathbf b \\) falls outside the k-point grid, under the convention above \\( c(\mathbf k+\mathbf b) = c(\mathbf k + \mathbf b - \mathbf G) \\) can be used as is).

[Convention] The existing implementation evaluates it as follows.

- The atom-pair integral \\( \langle\phi\_{I\alpha}\vert e^{-i\mathbf b\cdot(\mathbf r - \mathbf R\_I)}\vert\phi\_{J\beta}\rangle \\) is a numerical integral on the SCF real-space grid
  (the sum over grid points of the [real-space grid](grid-poisson-xc.md) multiplied by the volume element).
- \\( M \\) is the average of the evaluation with the phase reference on the bra-side atom and the evaluation of the same integral reinterpreted as the ket-side pair (for continuous integrals the two are equal; they differ only by the grid error).
- In noncollinear calculations, the sum over the two spinor components of the formula above is taken.

[Public] The initial projection is the overlap between Bloch states and localized functions (projectors) \\( g\_n \\) placed at positions \\( \mathbf p\_n \\):

$$
A\_{\mu n}(\mathbf k) = \frac{1}{N\_k}\langle\psi\_{\mu\mathbf k}\vert g\_n\rangle
{}= \sum\_{I\alpha}\sum\_{\mathbf L} e^{-i\mathbf k\cdot\mathbf L}\\, c\_{I\alpha,\mu}^{\*}(\mathbf k)\\, \langle\phi\_{I\alpha}(\mathbf r - \mathbf R\_I - \mathbf L)\vert g\_n(\mathbf r - \mathbf p\_n)\rangle .
$$

The projector \\( g\_n \\) is a PAO (or a hybrid of PAOs), and the overlap can be computed in the same way as the [two-center integrals](two-center-overlap-kinetic.md).
For isolated bands, the initial \\( U \\) is the polar decomposition of \\( A \\) (if \\( A = W\Sigma V^\dagger \\), then \\( U = WV^\dagger \\)); for entangled bands, the initial subspace is selected with the same projection.

### 1.4 Specification of projectors [Public] [Convention]

[Public] Atomic species for projectors are added to `Definition.of.Atomic.Species` in the same format as ordinary species (the pseudopotential field is required only by the format and is not used),
and each line of the `<Wannier.Initial.Projectors` block contains "`species name-orbital name`, the 3 coordinates of the center, the 3 components of the local z axis, the 3 components of the local x axis" (manual 47.2.1).
Since only one radial function per angular momentum is used, a specification such as `s1p1d1f1` is recommended.

The orbital names and numbers of projectors are as follows. Hybrids are linear combinations of PAOs in the local coordinate system (in the order of the [real spherical harmonics](pao-basis.md)).

| Orbital name | Number | Content |
|---|---|---|
| `s`, `px`, `py`, `pz`, `dz2`, `dx2-y2`, `dxy`, `dxz`, `dyz`, `fz3`, `fxz2`, `fyz2`, `fzx2`, `fxyz`, `fx3-3xy2`, `f3yx2-y3` | 1 | that single PAO |
| `p`, `d`, `f` | 3, 5, 7 | all components of that \\( l \\) |
| `sp`, `sp2`, `sp3` | 2, 3, 4 | standard hybrids of s and p (for \\( sp^3 \\), each coefficient has magnitude \\( 1/2 \\)) |
| `sp3dz2`, `sp3deg` | 5, 6 | hybrids including s, p and d (\\( d\_{z^2} \\), \\( d\_{x^2-y^2} \\)) |

Note: the table in manual 47.2.1 gives the coefficient of the \\( sp^3 \\) hybrid as \\( 1/\sqrt 2 \\), but since it is a sum of 4 components, normalization requires \\( 1/2 \\).
[Convention] The existing implementation also uses \\( 1/2 \\). Also, in the same table the descriptions of `dxz` and `dyz` read "\\( d\_{xy} \\)", which is a typo; they are \\( d\_{xz} \\) and \\( d\_{yz} \\), respectively.

[Convention]

- For the local coordinate system, Euler angles are obtained from the given z and x axes, and the projectors are rotated by the rotation matrices of the real spherical harmonics for \\( l = 1, 2, 3 \\) (s is rotation invariant).
- The coordinate unit `Wannier.Initial.Projectors.Unit` is `Ang` (default), `AU`, or `FRAC`. With `FRAC`, each component is wrapped into \\( [0, 1] \\) before conversion to Cartesian coordinates.
- If the total number of projectors does not match `Wannier.Func.Num`, it is an error when the initial guess is used (`Wannier.Initial.Guess on`).

### 1.5 Band selection [Convention]

- The window edges are given in eV relative to the chemical potential. At each k point and spin, the bands in the outer window \\( [\mu\_{\mathrm F} + E\_{\mathrm{ob}},\ \mu\_{\mathrm F} + E\_{\mathrm{ot}}) \\) (lower end included, upper end excluded) are selected. The inner window is treated the same way.
- The number of bands may differ between k points; the number of bands in the file output, \\( N\_{\mathrm{win}} \\), is its maximum.
- It is an error if at some k point the number of bands in the outer window is less than \\( N\_W \\), or the number of bands in the inner window exceeds \\( N\_W \\).
- If the difference between the upper and lower ends of the inner window is less than \\( 10^{-6} \\) eV, the inner window is not used (no freezing). Upper end \\( < \\) lower end is an error, and an inner window extending outside the outer window is also an error.

### 1.6 Shells of b vectors [Public] [Convention]

[Public] Shells (sets of \\( \mathbf b \\) of equal length) are searched in order of increasing distance on the k-point grid, and weights satisfying the 6 equations \\( \sum w\_b b\_\alpha b\_\beta = \delta\_{\alpha\beta} \\) are found
(the method of Eq. 26 of Mostofi et al., the Wannier90 paper arXiv:0708.0650). The upper limit on the shells searched is `Wannier.MaxShells` (manual 47.2.1).

[Convention] Shells are added one at a time; if a singular value below \\( 10^{-5} \\) appears in the 6 × (number of shells) coefficient matrix (linear dependence on existing shells), that shell is discarded and the search proceeds to the next;
the search ends when the least-squares solution satisfies the 6 equations within \\( 10^{-5} \\) in each component.

### 1.7 Minimization procedure [Public] [Convention]

[Public] `Wannier.Minimizing.Scheme` is 0 = steepest descent (SD), 1 = conjugate gradient (CG, with line search by the secant method), 2 = SD followed by CG.
The SD step is `Wannier.Minimizing.StepLength` (\\( \epsilon \\) in Sec. 1.1), and the initial step and maximum number of steps of the CG secant method are
`Wannier.Minimizing.Secant.StepLength` and `Wannier.Minimizing.Secant.Steps`.
The convergence criterion applies `Wannier.Minimizing.Conv.Criterion` to the norm of the gradient. In scheme 2, SD and CG each run up to `Wannier.Minimizing.Max.Steps` times (manual 47.2.3).

[Convention] The branch of \\( \mathrm{Im}\ln M\_{nn} \\) is chosen to be continuous with the previous iteration (if \\( \lvert M\_{nn}\rvert < 10^{-8} \\), only the branch correction is applied).

### 1.8 Keywords and defaults [Public] [Convention]

The defaults are values confirmed in the existing implementation. Those differing from the manual are annotated.

| Keyword | Default | Meaning |
|---|---|---|
| `Wannier.Func.Calc` | off | MLWF calculation |
| `Wannier.Func.Num` | 0 (effectively required) | \\( N\_W \\) |
| `Wannier.Outer.Window.Bottom` / `.Top` | −10.0 / 10.0 eV | outer window; the manual says "no default" |
| `Wannier.Inner.Window.Bottom` / `.Top` | 0.0 / 0.0 eV | inner window (not used if the values are equal) |
| `Wannier.Initial.Guess` | on | initial guess by projection |
| `Wannier.Initial.Projectors.Unit` | `Ang` | `Ang`, `AU`, `FRAC` |
| `Wannier.Kgrid` | 4 4 4 | k-point grid for MLWF; the manual says "no default" |
| `Wannier.MaxShells` | 12 | upper limit of the shell search |
| `Wannier.Dis.SCF.Max.Steps` | 200 | maximum iterations of disentanglement |
| `Wannier.Dis.Conv.Criterion` | \\( 10^{-8} \\) | ditto, convergence criterion (Å²) |
| `Wannier.Dis.Mixing.Para` | 0.5 | ditto, mixing coefficient |
| `Wannier.Minimizing.Scheme` | 0 | 0 = SD, 1 = CG, 2 = SD->CG |
| `Wannier.Minimizing.StepLength` | 2.0 | SD step |
| `Wannier.Minimizing.Secant.Steps` | 5 | maximum number of secant steps |
| `Wannier.Minimizing.Secant.StepLength` | 2.0 | initial secant step |
| `Wannier.Minimizing.Conv.Criterion` | \\( 10^{-8} \\) | convergence criterion on the gradient norm |
| `Wannier.Minimizing.Max.Steps` | 200 | maximum iterations |
| `Wannier.Readin.Overlap.Matrix` | off | restart by rereading `.mmn`; the manual says "default on", but the implementation uses off |
| `Wannier.Interpolated.Bands` | off | interpolated bands (requires `Band.dispersion on`; error otherwise) |
| `Wannier.Function.Plot` | off | cube output of MLWFs |
| `Wannier.Function.Plot.SuperCells` | 0 0 0 | range of cells to plot; error if \\( 2l\_i + 1 \\) exceeds the corresponding component of `Wannier.Kgrid` |
| `Wannier90.fileout` | off | output of files for Wannier90 |

[Convention] On restart, it is checked that the outer window and the k-point grid are the same as at saving time. The inner window, initial guess, and convergence conditions may be changed.

### 1.9 Output [Public] [Convention]

[Public] (manual 47.2.2, 47.2.5, 47.2.7)

| File | Content |
|---|---|
| `System.Name.mmn` | \\( M\_{mn}(\mathbf k,\mathbf b) \\). Wannier90 format. Line 2 has \\( N\_{\mathrm{win}} \\), number of k points, number of \\( \mathbf b \\), number of spins. Blocks are ordered spin -> k point -> \\( \mathbf b \\); each block starts with the k-point index, the index of \\( \mathbf k+\mathbf b \\), and 3 lattice integers |
| `System.Name.amn` | \\( A\_{\mu n}(\mathbf k) \\). Wannier90 format (band index, Wannier index, k-point index, real part, imaginary part) |
| `System.Name.eigen` | chemical potential, number of bands, eigenvalues and eigenvectors (for restart) |
| `System.Name.eig` | eigenvalues for Wannier90 (eV relative to the chemical potential) |
| `System.Name.HWR` | \\( t\_{mn}(\mathbf L) \\). \\( \mathbf L \\) in the Wigner–Seitz supercell and degeneracies, lattice vectors (bohr), chemical potential and values in Hartree |
| `System.Name.Wannier_Band` | interpolated bands (same format as `.Band`) |
| `System.Name.mlwf{s}_{n}_{r,i}.cube` | real and imaginary parts of MLWFs; \\( s \\) is the spin, \\( n \\) the MLWF index |
| `System.Name.win` | Wannier90 control file (when `Wannier90.fileout on`) |

The standard output shows disentanglement (lines ending in `DISE`, \\( \Omega\_{\mathrm I} \\) and its change), minimization (`CONV`, gradient norm and \\( \Delta\Omega \\)),
the breakdown of the spread (`SPRD`, \\( \Omega\_{\mathrm I} \\), \\( \Omega\_{\mathrm D} \\), \\( \Omega\_{\mathrm{OD}} \\), total), and the centers (`CENT`, Å and Å²).

Note: manual 47.2.2 writes the cube file name as `.mlwfl_…`, but [Convention] the actual name is `.mlwf{s}_{n}_r.cube`.
Likewise, the extension of the interpolated bands is `.Wannier_Band`.

---

## 2. Closest Wannier functions (CWF)

### 2.1 Method [Public]

Technical note CWF (Ozaki PRB **110**, 125115, Eqs. 2–18) is rewritten in the notation of this note.
From guiding functions \\( Q\_{\mathbf 0 p} \\) (\\( p = 1, \dots, N\_{\mathrm{CWF}} \\), localized in the home cell) and a window function \\( w(\varepsilon) \\), the projection matrix

$$
a\_{\mu p}(\mathbf k) = w(\varepsilon\_{\mu\mathbf k})\\, \frac{1}{N\_k}\langle\psi\_{\mu\mathbf k}\vert Q\_{\mathbf 0 p}\rangle
$$

(\\( N\_{\mathrm{band}}\times N\_{\mathrm{CWF}} \\), \\( N\_{\mathrm{CWF}} \le N\_{\mathrm{band}} \\)) is built.
The sum of squared distances between the window-weighted projection functions \\( L\_{\mathbf R p} \\) and orthonormal Wannier functions \\( W\_{\mathbf R p} \\) (coefficients \\( B(\mathbf k) \\), \\( B^\dagger B = I \\)),

$$
F[B] = \frac{1}{N\_k}\sum\_{\mathbf k}\left\lVert A(\mathbf k) - B(\mathbf k)\right\rVert\_{\mathrm F}^2
$$

is minimized by \\( B = U \\) from the polar decomposition \\( A = UP \\) of \\( A(\mathbf k) \\) (from the singular value decomposition \\( A = W\Sigma V^\dagger \\), \\( U = WV^\dagger \\), \\( P = V\Sigma V^\dagger \\)), with

$$
F\_{\min} = \frac{1}{N\_k}\sum\_{\mathbf k}\sum\_p\left(\sigma\_{p\mathbf k} - 1\right)^2
$$

No iteration is needed. The properties are as follows.

- If all singular values of \\( A \\) are positive, \\( U \\) is unique and independent of the phases of the Bloch states and of the phase and degeneracy freedoms of the singular vectors (no gauge problem).
- If \\( A^\dagger A \\) is nonsingular, \\( U = A(A^\dagger A)^{-1/2} \\), i.e., the same as windowed Löwdin orthogonalization. The polar decomposition does not form \\( (A^\dagger A)^{-1/2} \\), so it is stable even when ill-conditioned.
- If the window is narrow and fewer than \\( N\_{\mathrm{CWF}} \\) states fall inside it at some k point, a singular value becomes 0, \\( U \\) is no longer unique, and symmetry breaks. This is prevented by adding a small constant \\( \delta \\) to the window function.
- The hoppings and Wannier interpolation are the formulas of the "Notation" section with \\( U \\) inserted.
- Since \\( N\_{\mathrm{CWF}} \\) functions are built from \\( N\_{\mathrm{band}} \\) states, disentanglement is always done automatically by the window function.

### 2.2 Window functions [Public] [Convention]

[Public] The two types of manual 47.1.1. Let \\( x = \varepsilon - \mu\_{\mathrm F} \\) (relative to the chemical potential).

**Polynomial window** (`Poly`): give \\( \epsilon\_0 < \epsilon\_1 \le \epsilon\_2 < \epsilon\_3 \\) (eV) in `CWF.disentangling.Erange`, and

$$
w\_P(x) = \begin{cases}
0 & x \le \epsilon\_0 \\\\
P\\!\left(\dfrac{\epsilon\_1 - x}{\epsilon\_1 - \epsilon\_0}\right) & \epsilon\_0 < x \le \epsilon\_1 \\\\
1 & \epsilon\_1 < x \le \epsilon\_2 \\\\
P\\!\left(\dfrac{x - \epsilon\_2}{\epsilon\_3 - \epsilon\_2}\right) & \epsilon\_2 < x \le \epsilon\_3 \\\\
0 & \epsilon\_3 < x
\end{cases},\qquad P(y) = 1 - 6y^2 + 8y^3 - 3y^4 .
$$

Since \\( P(0) = 1 \\), \\( P'(0) = 0 \\), \\( P(1) = P'(1) = P''(1) = 0 \\), the window is \\( C^1 \\) at the inner edges and \\( C^2 \\) at the outer edges.
Since it is exactly 0 in the outer intervals, those states can be excluded from the calculation. [Convention] \\( \delta \\) is not added to the polynomial window.

**Fermi window** (`Fermi`): give \\( \epsilon\_0, \epsilon\_1 \\) in `CWF.disentangling.Erange` and \\( k\_BT\_0, k\_BT\_1 \\) (eV) in `CWF.disentangling.smearing.kBT`, and

$$
w\_F(x) = \frac{1}{\left(1 + e^{(\epsilon\_0 - x)/k\_BT\_0}\right)\left(1 + e^{(x - \epsilon\_1)/k\_BT\_1}\right)} + \delta,
\qquad \delta = \texttt{CWF.disentangling.smearing.bound} .
$$

Note: Eq. 4 of technical note CWF has the form of the sum of two Fermi functions minus 1, \\( f\_0 + f\_1 - 1 \\), which differs from the product \\( f\_0 f\_1 \\) above
(they are nearly the same if the window is sufficiently wider than \\( k\_BT \\)). [Convention] The existing implementation uses the product form, as in the manual.

### 2.3 Guiding functions [Public] [Convention]

Selected from three types with `CWF.Guiding.Orbitals` (manual 47.1.1–47.1.4, technical note Sec. III B).

- **AO** (atomic orbitals): the PAOs themselves. In `<CWF.Guiding.AOs`, write one line per species with names "`1s`, `2s`, ..., `1p`, ..." and 0/1 in the order of the basis \\( (l, \zeta) \\) (the order of the [PAO basis](pao-basis.md)).
  All \\( (l,\zeta) \\) of the basis must be given explicitly as 0 or 1.
  The number of CWFs per atom is \\( \sum\_{(l,\zeta)\ \mathrm{selected}}(2l+1) \\).
- **HO** (hybrid orbitals): in the space of the selected AOs, diagonalize the on-site density-matrix block of that atom \\( \rho\_{\mathbf 0 I\alpha,\mathbf 0 I\beta} \\), and
  use, in order of decreasing occupation, as many eigenvectors as there are AOs. For the same selected AOs, the AO and HO CWFs span the same subspace.
  [Convention] The density matrix is the per-spin SCF density matrix (with Fermi-distribution occupations); no window function is used. Since the PAOs are orthonormal on site ([PAO basis](pao-basis.md)),
  orthonormality of the eigenvectors under the Euclidean inner product directly gives orthonormality of the functions. In spin-polarized calculations they are built per spin.
- **MO** (embedded molecular orbitals): all atoms are divided into groups \\( g \\) with `<CWF.MO.Grouped.Atoms`, and the right singular vectors of the group block
  \\( \Lambda\_g = \sum\_{\mathbf R g'}\rho\_{\mathbf 0 g,\mathbf R g'}\\, S\_{\mathbf R g',\mathbf 0 g} \\) (the diagonal group block of \\( \rho S \\)) are used.
  As in Eqs. 31–33 of the technical note, \\( \Lambda\_g^\dagger\Lambda\_g \\) is diagonalized to obtain real singular vectors \\( y\_{g\nu} \\),
  which are normalized with the on-site overlap within the group \\( S\_{\mathbf 0 g,\mathbf 0 g} \\), and the expectation values \\( \langle\bar y\vert\hat H\vert\bar y\rangle \\) are computed.
  The number used per group is specified with `<CWF.MO.Selection`.

[Convention] Details of MO guiding functions:

- \\( \rho \\) is the density matrix with the window function \\( w(\varepsilon) \\) used in place of the occupations (as described in manual 47.1.4). In spin-polarized calculations the average over both spins (for both \\( \rho \\) and \\( H \\)) is used.
- The selected ones are the specified number from the top **in order of decreasing singular value**. The technical note says they are sorted and selected by expectation value, but the existing implementation uses the singular-value order.
- `System.Name.CWF_MOs_info` lists, per group, the index, expectation value (eV relative to the chemical potential), and singular value (displayed as "population") in order of decreasing singular value.
- Group numbers must be consecutive starting from 1 (the number of groups is the number of indices used, and the indices are used directly as subscripts).
  The lines of `CWF.MO.Selection` are also ordered from group 1.

[Design] Since the guiding function \\( Q \\) is a linear combination of PAOs, \\( \langle\psi\_{\mu\mathbf k}\vert Q\rangle/N\_k = \left[C(\mathbf k)^\dagger S(\mathbf k)\\, q\right]\_{\mu} \\) can be used (\\( q \\) is the PAO coefficients of \\( Q \\) placed in the global ordering).

### 2.4 Band range and k points [Convention]

- For the polynomial window, the range of band indices falling in the interval where the window is nonzero is determined from **the eigenvalues at the first k point of the SCF k-point grid (spin 1)**, and the same index range is used at all k points.
  It is an error if the number of bands in the range is less than \\( N\_{\mathrm{CWF}} \\).
- For the Fermi window, all bands (\\( N\_{\mathrm{band}} \\) = number of basis functions) are used.
- For periodic systems (`scf.EigenvalueSolver Band`), the calculation is done on the SCF k-point grid (`scf.Kgrid`); of k points paired by time reversal, only one is computed, with weight 2.
  Since the hoppings are built for \\( \mathbf L = \sum\_i l\_i\mathbf a\_i \\), \\( l\_i = -(N\_i - 1)/2, \dots, (N\_i - 1)/2 \\) (without a Wigner–Seitz supercell or degeneracies),
  **the number of k points in each direction must be odd** (manual 47.1.7). In collinear calculations only the real part of the hoppings is stored.
- For isolated systems (`Cluster`), only the Γ point is used, and only \\( \mathbf L = \mathbf 0 \\).
- It is an error if \\( 2l\_i + 1 \\) of `CWF.Plot.SuperCells` exceeds the number of k points.

### 2.5 Output of the distance function and effective charges [Public] [Convention]

[Convention] The output distance function (`DM_func`) is \\( \sum\_{\mathbf k} \omega\_{\mathbf k}\sum\_p(\sigma\_{p\mathbf k} - 1)^2 / \sum\_{\mathbf k}\omega\_{\mathbf k} \\) (\\( \omega\_{\mathbf k} \\) are the k-point weights),
doubled in non-spin-polarized calculations. The accompanying "per (CWF*spin)" is that value divided by \\( 2N\_{\mathrm{CWF}} \\).
Manual 47.1.1 states that a value above about 0.7 is an indication that the symmetry of the guiding functions does not match the Bloch states [Public].

[Public] From the CWF occupations (technical note Eq. 34)

$$
n\_p = \frac{1}{N\_k}\sum\_{\mathbf k\mu} f(\varepsilon\_{\mu\mathbf k})\\, \lvert u\_{\mu p}(\mathbf k)\rvert^2
$$

the effective charge of an atom \\( \kappa\_I = Z\_I^{\mathrm{val}} - \sum\_{p\in I} n\_p \\) is obtained (\\( n\_p \\) doubled in non-spin-polarized calculations).
The `.out` file contains the occupations per CWF and per atom (per spin, sum, difference) and the overall sum (manual 47.1.8).
\\( f \\) is the Fermi distribution at the same electronic temperature as the SCF.

### 2.6 Keywords and defaults [Public] [Convention]

| Keyword | Default | Meaning |
|---|---|---|
| `CWF.Calc` | off | CWF calculation (when `scf.EigenvalueSolver` is `Cluster` or `Band`) |
| `CWF.Guiding.Orbitals` | `HO` | `AO`, `HO`, `MO` |
| `CWF.window.function` | `Poly` | `Poly`, `Fermi` |
| `CWF.disentangling.Erange` | Poly: −1000 −100 0 3 / Fermi: −1000 0 (eV) | window edges; the same keyword takes a different number of values |
| `CWF.disentangling.smearing.kBT` | 0.3 0.3 (eV) | \\( k\_BT\_0, k\_BT\_1 \\) of the Fermi window |
| `CWF.disentangling.smearing.bound` | \\( 10^{-10} \\) | \\( \delta \\) of the Fermi window; the manual says \\( 10^{-12} \\) (see the note below) |
| `<CWF.Guiding.AOs` | — | selection for AO and HO |
| `<CWF.MO.Grouped.Atoms` / `<CWF.MO.Selection` | — | MO groups and numbers |
| `CWF.Dis.vs.H` | off | text output of distances and hoppings |
| `CWF.fileout.flag` | off | cube output of CWFs |
| `CWF.Plot.SuperCells` | 0 0 0 | range of cells to plot |
| `CWF.fileout.Num` / `<CWF.fileout.cube` | 0 / — | indices of atoms (AO, HO) or groups (MO) to plot |
| `CWF.Wannier90.fileout` | off | `.mmn`, `.amn`, `.eig`, `.win` for Wannier90 |

Note: [Convention] if `Wannier.Func.Calc on` is specified at the same time, the default of `CWF.disentangling.smearing.bound` is reread as \\( 10^{-12} \\) when the MLWF keywords are read.
In the usual case where only CWF is used, the default is \\( 10^{-10} \\).

### 2.7 Output [Public]

(manual 47.1.1, 47.1.4, 47.1.6, 47.1.9, 47.1.10)

| File | Content |
|---|---|
| `System.Name.CWF.Band` | Wannier-interpolated bands (along the `Band.dispersion` path) |
| `System.Name.CWF.Dis_vs_H` | spin, atom GA, atom GB, \\( l\_1, l\_2, l\_3 \\), CWF indices \\( i, j \\), distance between centers (Å), hopping (eV). No spin column in noncollinear calculations |
| `System.Name.CWF.Hop` | binary of hoppings (for auxiliary programs) |
| `System.Name.CWF_MOs_info` | list of MO guiding functions (Sec. 2.3) |
| `System.Name.cwf_{s}_{atom}_{p}_r.cube` | collinear CWFs (always real) |
| `System.Name.cwf_{atom}_{p}.nccube` | noncollinear CWFs (in the order \\( \alpha \\) real part, \\( \alpha \\) imaginary part, \\( \beta \\) real part, \\( \beta \\) imaginary part) |
| `System.Name_{s}.amn` etc. | for Wannier90 (with spin suffix; no suffix for noncollinear) |

---

## 3. Band unfolding

### 3.1 Formulation [Public]

Lee–Yamada-Takamura–Ozaki (JPCM **25**, 345501, Eqs. 1–26 of arXiv:1211.4690) is rewritten in the notation of this note.

Eigenstates of a supercell (lattice vectors \\( \mathbf L \\), atoms \\( I \\)) are represented in the Brillouin zone of a smaller reference cell (lattice vectors \\( \mathbf r \\), volume ratio \\( n = V\_{\mathrm{super}}/V\_{\mathrm{ref}} \\)).
As preparation,

- each atom \\( I \\) of the supercell is assigned an "atom type" \\( g(I) \\) within the reference cell (the group number of `Unfolding.Map`);
- the reference-cell lattice vector containing position \\( \mathbf x \\) is \\( \mathbf r(\mathbf x) = \sum\_i\lfloor\mathbf b\_i^{\mathrm{ref}}\cdot(\mathbf x - \mathbf o)/2\pi\rfloor\\, \mathbf a\_i^{\mathrm{ref}} \\) (\\( \mathbf o \\) is the origin of the reference cell).

When the supercell eigenvalue problem is solved at a reference-cell k point \\( \mathbf k \\) (the supercell Bloch phase \\( e^{i\mathbf k\cdot\mathbf L} \\) is the same as at \\( \mathbf k \\) folded into the supercell BZ),
the unfolding weight of state \\( \mu \\) at \\( \mathbf k \\) is

$$
W\_\mu(\mathbf k) = \sum\_{I\alpha}\sum\_{J\beta,\mathbf L}\ \sum\_{M:\\, g(M) = g(J)}
c^{\*}\_{I\alpha,\mu}(\mathbf k)\\, c\_{M\beta,\mu}(\mathbf k)\\,
e^{i\mathbf k\cdot\left[\mathbf r(\mathbf R\_J + \mathbf L) - \mathbf r(\mathbf R\_M)\right]}\\, S\_{I\alpha,J\beta}(\mathbf L)
$$

and the spectral weight is \\( W\_\mu(\mathbf k)/n \\). The sum runs over all orbitals \\( (J\beta, \mathbf L) \\) overlapping with the home-cell orbital \\( I\alpha \\),
each reinterpreted as "the translation of atom \\( M \\) of the same type \\( g(J) \\) to the reference cell \\( \mathbf r(\mathbf R\_J + \mathbf L) \\)"
(\\( \beta \\) is the same orbital index for \\( J \\) and \\( M \\); atoms in the same group must have the same species and the same basis).
Owing to the two approximations of the original paper (substituting the actual supercell positions for the positions of the reference-cell orbitals), \\( S \\) is just the overlap of the supercell calculation.

Properties (checked with `tools/check_wannier_unfolding.py`):

- Sum rule: the sum over the \\( n \\) points obtained by shifting \\( \mathbf k \\) by supercell reciprocal lattice vectors, \\( \sum\_{\mathbf G} W\_\mu(\mathbf K + \mathbf G)/n = 1 \\) (from \\( \mathbf c^\dagger S(\mathbf K)\mathbf c = 1 \\)). It holds exactly even with defects.
- For a perfect supercell (no defects), \\( W\_\mu/n \\) is 0 or 1, and the energies of states with weight 1 coincide with the reference-cell bands. In this case the result does not depend on the choice of origin \\( \mathbf o \\).
- Orbital decomposition: \\( W\_\mu = \sum\_{M\beta} W\_{\mu,M\beta} \\) (the partial sum of the formula above with \\( M\beta \\) fixed) shows which orbitals form a band (Eq. 26 of the original paper).

Note: Eq. 25 of the original paper (and Eq. 16 of technical note unfolding_otw) attaches the complex conjugate to the coefficient \\( C\_M \\) on the reinterpreted side.
Used literally with the Bloch convention of this note (\\( C^\dagger S C = 1 \\), phase on the ket side), the sum rule above fails for systems with defects (a deviation of about 1.5 in the check).
The conjugate is attached to the coefficient on the bra side (the home-cell orbital \\( I\alpha \\)). [Convention] The existing implementation does the same.

### 3.2 Origin of the reference cell and mapping [Public] [Convention]

[Public] The reference cell is given with `<Unfolding.ReferenceVectors` (unit `Atoms.UnitVectors.Unit`; if omitted, the supercell itself),
and the group number of each atom with `<Unfolding.Map` (any positive integers; the values and order of the numbers are arbitrary). The origin can be given explicitly with `<Unfolding.ReferenceOrigin` (manual 55.2, 55.3).

[Convention]

- If `Unfolding.Map` is omitted, every atom forms its own group. Group numbers are renumbered contiguously in order of appearance.
- Default origin: placed at the center of the supercell minus the center of the reference cell, and each Cartesian component is further shifted by −0.1 bohr.
  If the fractional coordinates of some atom (or its periodic image) relative to the reference cell are within \\( \tfrac12(0.1,0.1,0.1)\cdot\mathbf b\_i^{\mathrm{ref}}/2\pi \\) of integers in all 3 components, the origin is shifted by −0.1 bohr once more.
  When the origin is given explicitly, no shift is applied.
- It is an error if different periodic images of the same atom fall into the same reference cell (e.g., when the reference cell is larger than the supercell in some direction).
- Note: manual 55.3 explains that the origin is determined automatically so that "atoms of the same group fall, as far as possible, one into each reference cell" and "the number of reference cells containing atoms is minimal", but
  the implementation in the version checked uses the fixed rule above. For surfaces and strongly distorted systems, giving the origin explicitly is recommended.

### 3.3 k points, energies, and output [Public] [Convention]

[Public] (manual 55.1)

- In `<Unfolding.kpoint`, `Unfolding.Nkpoint` points are listed in fractional coordinates with respect to the reciprocal lattice vectors of the **reference cell** (name, 3 coordinates). The polyline connecting consecutive points is the path.
- `Unfolding.desired_totalnkpt` is a guideline for the number of points along the whole path.
- The energy range is `Unfolding.LowerBound`, `Unfolding.UpperBound` (eV relative to the chemical potential).

[Convention]

- With total path length \\( D \\) (bohr⁻¹) and the guideline number of points \\( N\_d \\): if \\( N\_d \\) is at most the number of segments, each segment gets 2 points (both ends);
  otherwise segment \\( i \\) gets \\( \max(1, \lfloor D\_i/\Delta\rfloor) + 1 \\) points, \\( \Delta = D/(N\_d - 1) \\), equally spaced including both ends.
  Points at segment junctions are computed twice and have the same horizontal-axis distance.
- As an alternative format, `UnfoldBand.Nkpath` and `<UnfoldBand.kpath` (each line: 3 start coordinates, 3 end coordinates, 2 names) are also accepted (not documented in the manual).
- Only states with \\( \mu\_{\mathrm F} + E\_{\mathrm{lower}} \le \varepsilon \le \mu\_{\mathrm F} + E\_{\mathrm{upper}} \\) (both ends included) are output.
- The eigenvalue problem is solved with the rule that replaces negative eigenvalues of \\( S \\) by \\( 10^{-14} \\) (rule A′ of [generalized eigenvalue problem](generalized-eigen.md)).
- \\( \lvert\mathrm{Re}\\, W\_\mu\rvert/n \\) is written as the total weight and \\( \mathrm{Re}\\, W\_{\mu,M\beta}/n \\) (replaced by 0 if negative) as the orbital decomposition.
- Disabled when `scf.EigenvalueSolver` is other than `Cluster` or `Band`. For calculations with orbital optimization, \\( H \\) and \\( S \\) in the contracted basis are used.

Output files [Public]:

| File | Content |
|---|---|
| `System.Name.unfold_totup` (`_totdn`; `_tot` for noncollinear) | cumulative k distance (bohr⁻¹), energy (eV, relative to the chemical potential), total weight |
| `System.Name.unfold_orbup` (`_orbdn`; `_orb` for noncollinear) | the same 2 columns followed by per-orbital weights (ordering given in `.out`) |
| `System.Name.unfold_plotexample` | gnuplot example |

[Public] The intensity map is made by broadening the weight \\( w \\) with \\( L(k,E) = w/\left[(k/\Delta\_k)^2 + (E/\Delta\_E)^2 + 1\right] \\) (manual 55.4, auxiliary program).
The absolute values have no physical meaning. It can be used with noncollinear, spin-orbit, and DFT+U calculations without additional keywords (manual 55.5).

### 3.4 Unfolding of optical transition weights [Public]

Technical note unfolding_otw (Lee et al. PRB **102**, 075143) describes a method that splits the Kubo–Greenwood optical conductivity into contributions per k point and state (partial optical conductivity),
and transfers them to the reference-cell BZ with the weights of Sec. 3.1. In the notation of this note, for the partial optical conductivity of the supercell

$$
\sigma\_{ab}(\mathbf K, J, \omega) = \frac{-i}{\Omega}\sum\_{J'}\frac{f\_{\mathbf KJ} - f\_{\mathbf KJ'}}{\varepsilon\_{\mathbf KJ} - \varepsilon\_{\mathbf KJ'}}\\,
\frac{\langle\mathbf KJ\vert\hat P\_a\vert\mathbf KJ'\rangle\langle\mathbf KJ'\vert\hat P\_b\vert\mathbf KJ\rangle}{\varepsilon\_{\mathbf KJ} - \varepsilon\_{\mathbf KJ'} + \omega + i\eta}
$$

(\\( \Omega \\) is the supercell volume; for degeneracies the fraction is replaced by the energy derivative of the occupation), the partial optical conductivity in the reference cell is obtained
by transferring the spectral-function tensor with the same weight \\( W\_J(\mathbf k)/n \\) as in Sec. 3.1 and replacing the volume by that of the reference cell \\( \Omega\_{\mathrm{ref}} = \Omega/n \\):

$$
\sigma\_{ab}(\mathbf k, j, \omega) = \frac{\Omega}{n\\,\Omega\_{\mathrm{ref}}}\sum\_{\mathbf G}\delta\_{\mathbf k - \mathbf G,\mathbf K}\\, W\_J(\mathbf k)\\, \sigma\_{ab}(\mathbf K, J, \omega)
{}= \sum\_{\mathbf G}\delta\_{\mathbf k - \mathbf G,\mathbf K}\\, W\_J(\mathbf k)\\, \sigma\_{ab}(\mathbf K, J, \omega)
$$

(Eqs. 15–17 of the original paper rewritten with the normalization of \\( W \\) in this note; the factor \\( (L/l)^2 \\) in the original depends on its definitions of \\( L \\) and \\( l \\)).
The result integrated over a frequency interval is overlaid on the band diagram.
The momentum matrix elements can be written with LCAO coefficients and \\( \langle\phi\vert\nabla\vert\phi\rangle \\) between atomic orbitals, and decomposed into sets of four orbitals (Eqs. 19–24 of the original paper).

[Convention] The version checked has no input keyword that outputs unfolded optical transition weights (optical conductivity and unfolding can each be computed independently).

---

## 4. Natural atomic orbital (NAO) analysis and O(N) orbital analysis

### 4.1 Method [Public]

NAOs and natural populations (NP), the first stage of the natural bond orbital method of Weinhold et al. (Reed–Weinstock–Weinhold, JCP **83**, 735 (1985)), are computed (manual 36).
From the density matrix \\( P \\) (per spin) and the overlap \\( S \\), they are built in the following order.

1. For each atom \\( A \\) and angular momentum \\( l \\), with \\( P\_{Al} \\), \\( S\_{Al} \\) averaged over the \\( (2l+1) \\) components (\\( N\_l\times N\_l \\) in \\( \zeta \\)),
   solve \\( P\_{Al}\mathbf n = w\\, S\_{Al}\mathbf n \\) to obtain the pre-NAOs and occupations \\( w \\).
2. Those with occupation at or above the threshold form the natural minimal basis (NMB), the rest the natural Rydberg basis (NRB).
3. The NMB is orthogonalized between atoms by occupation-weighted symmetric orthogonalization (OWSO); the NRB is Schmidt-orthogonalized against the NMB, and then the NRB functions are orthogonalized among themselves.
4. In the orthogonalized basis, the atomic-block density matrix is again diagonalized per \\( (A, l) \\) to obtain the NAOs. The NP of an atom is its trace.

[Public] Ohwaki–Otani–Ozaki (JCP **140**, 244105) used the truncated-cluster idea of O(N) methods to obtain NAOs and NPs
from only \\( P \\), \\( S \\) of the cluster around specified central atoms, enabling local analysis of large systems.

### 4.2 Input and conventions [Public] [Convention]

| Keyword | Default | Meaning |
|---|---|---|
| `NBO.switch` | `off` | `on1`: use the density matrix of the whole system in a Cluster calculation. `on2`: central atoms only in a Krylov (O(N)) calculation |
| `NAO.threshold` | 0.85 | occupation threshold separating NMB and NRB (electrons per spin orbital) |
| `NBO.Num.CenterAtoms` | 1 | number of central atoms for `on2` (error if it exceeds the total number of atoms) |
| `<NBO.CenterAtoms` | — | indices of the central atoms (effective only for `on2`) |
| `NAO.only` | on | compute only NAOs and NPs |
| `NHO.fileout`, `NBO.fileout` | off | cube output of natural hybrid orbitals and natural bond orbitals |
| `NBO.SmallCell.Switch`, `<NBO.SmallCell` | off | restrict the output range in fractional coordinates |

[Convention]

- `on1` works only with `scf.EigenvalueSolver Cluster`, `on2` only with `Krylov`. For other combinations, a message is printed and nothing is done.
  `on3` and `on4` (for periodic systems) are accepted, but their processing is disabled in the version checked.
- In spin-polarized (collinear) calculations, the threshold is `NAO.threshold`/2 and is compared with the per-spin occupations (in non-spin-polarized calculations, `NAO.threshold` as is).
- The run stops if an eigenvalue below \\( 10^{-14} \\) appears in the atomic block \\( S\_{Al} \\).
- The cluster of `on2` is the same truncated cluster as in the Krylov method (the central atom, atoms with overlapping basis, and atoms further outside) (see [Open](#open)).
- The results are written to standard output as the NP per atom and its total, and the occupations, energies and PAO expansion coefficients of the NAOs per atom.

---

## 5. Orbital optimization of the basis (generation of contracted PAOs)

### 5.1 Method [Public]

This is the method of technical note 5 (Ozaki PRB **67**, 155108, Eqs. 1–12) and manual 26. Basis functions are expressed as linear combinations of primitive PAOs (radial functions of the same atom and same \\( l \\)),
and their coefficients (contraction coefficients) \\( a \\) are optimized variationally. In the notation of this note,

$$
\phi\_{I,l\zeta m}(\mathbf r) = \sum\_{p=0}^{N^{\mathrm{prim}}\_l - 1} a\_{I,l\zeta,p}\\, \chi\_{I,lpm}(\mathbf r),\qquad
\chi\_{I,lpm}(\mathbf r) = R^{(I)}\_{lp}(r)\\, S\_{lm}(\hat{\mathbf r})
$$

Since \\( a \\) does not depend on \\( m \\) (restricted optimization), the total energy is invariant under rotation of the system.
After solving the SCF with LCAO coefficients \\( c \\), by the same reasoning as the force theorem,

$$
\frac{\partial E\_{\mathrm{tot}}}{\partial a\_{I\alpha,p}} = 2\sum\_{\sigma}\sum\_{J\beta}
\left(\rho^{\sigma}\_{I\alpha,J\beta}\\, \langle\chi\_{Ip}\vert\hat H^{\sigma}\vert\phi\_{J\beta}\rangle - E^{\sigma}\_{I\alpha,J\beta}\\, \langle\chi\_{Ip}\vert\phi\_{J\beta}\rangle\right)
$$

is obtained (\\( \rho \\) is the density matrix and \\( E \\) the energy density matrix, both in the contracted basis, summed over all real-space atom pairs and periodic images; in non-spin-polarized calculations the spin sum is taken as twice).
This gradient is determined solely by \\( \rho \\), \\( E \\) in the contracted basis and \\( H \\), \\( S \\) in the primitive basis; no additional derivatives are needed.
In restricted optimization, the sum of the gradients over the \\( (2l+1) \\) values of \\( m \\) is taken as the gradient for all \\( m \\) of that \\( (l,\zeta) \\).

(checked) For a finite cluster with a fixed Hamiltonian, this gradient was confirmed to agree with the difference quotient of the band energy to about \\( 10^{-8} \\) (`tools/check_wannier_unfolding.py`).

### 5.2 Specification [Public] [Convention]

[Public] The basis is written as "number of primitives `>` number of contracted functions", e.g. `H5.0-s4>1`, `C5.0-s4>1p4>1` ([PAO basis](pao-basis.md)).
The primitive functions are the first \\( N^{\mathrm{prim}}\_l \\) functions of the PAO file. `orbitalOpt.Method` selects per atom (`Atoms`) or per species (`Species`, coefficients shared among atoms with the same species name).

[Convention]

- `orbitalOpt.Method` is `Off` (default), `Atoms`, `Species`, `Atoms2`, or `Species2`. `Atoms2` and `Species2` optimize only the contraction coefficients of the **last** \\( \zeta \\) for each \\( l \\), fixing the rest at their initial values (not documented in the manual).
- For `Species`, the gradient of a species is the sum of the gradients of all atoms of that species.
- The contraction coefficients are normalized at every update so that \\( \sum\_{pq}a\_p a\_q\\, S^{(I)}\_{pq} = 1 \\) (on-site overlap of the primitive functions). Since the PAOs are orthonormal, this is effectively \\( \sum\_p a\_p^2 = 1 \\).
- One orbital-optimization step is done every `orbitalOpt.scf.maxIter` SCF iterations, or earlier when the SCF converges (energy change below `scf.criterion` and residual norm below \\( 10^{-7} \\)).
- The "gradient norm" of the convergence criterion is the sum of squares of the gradient components (no square root; counted per \\( m \\) component). When it falls below `orbitalOpt.criterion` or the number of steps exceeds `orbitalOpt.Opt.maxIter`, one final run of `orbitalOpt.scf.maxIter` SCF iterations is done and the procedure ends.
- In MD and geometry optimization, the optimization is redone every `orbitalOpt.per.MDIter` steps (at steps 1, 1+P, 1+2P, ...).
- Incompatible features: NEGF (input error), `RMM-DIISH` and `RMM-DIISV` ([density mixing](density-mixing.md)), stress ([Total energy and forces](energy-force.md)).

### 5.3 Update rules [Public] [Convention]

[Public] Steps before `orbitalOpt.StartPulay` use steepest descent, and later steps a quasi-Newton method (`EF`: eigenvector following, `DIIS`), with `orbitalOpt.HistoryPulay` history entries (manual 26).

[Convention] With the number of history entries \\( d = \min(\texttt{HistoryPulay},\ \mathrm{step} - \texttt{StartPulay} + 1) \\) (1 if less than 1),

| Mode | \\( d = 1 \\) (steepest descent) | \\( d > 1 \\) |
|---|---|---|
| `Atoms` | \\( a \leftarrow a - \lambda g \\), \\( \lambda = \min(4\\, s/\lVert g\rVert^2,\ 0.05) \\) | DIIS extrapolation \\( \bar a - \gamma\bar g \\). \\( \gamma \\) is 0.1, 0.2, 0.4, 0.5 when the gradient norm exceeds \\( 10^{-2}, 10^{-3}, 10^{-4}, 10^{-5} \\), and 0.7 below that (DIIS even if `EF` is specified) |
| `Species` | \\( a \leftarrow a - (s/\max\_i\lvert g\_i\rvert)\\, g \\) | first build \\( \bar a \\), \\( \bar g \\) by DIIS; for `EF`, compute a quasi-Newton step with the BFGS-updated Hessian (eigenvalues bounded below by 0.1), shrink it to 0.001 if its largest component exceeds 0.001, and subtract it from \\( \bar a \\). For `DIIS`, \\( \bar a - \gamma\bar g \\) (\\( \gamma \\) is 0.05, 0.10, 0.15, 0.20, 0.40) |

Here \\( s \\) is `orbitalOpt.SD.step` and \\( \lVert g\rVert^2 \\) is the gradient norm above.

### 5.4 Initial coefficients [Public] [Convention]

[Public] Eqs. 11–12 of technical note 5: for each atom, the Kohn–Sham equation of a cluster including neighboring atoms is solved non-self-consistently, and the initial coefficients are built with Fermi-distribution weighting.

[Convention] `orbitalOpt.InitCoes` is `Symmetrical` (default), `Free`, or `Simple`.

- `Symmetrical`, `Free` (the same processing in the version checked): the cluster of the atom and its neighbors is diagonalized in the primitive basis, and with Fermi-distribution weighting at a local chemical potential determined from the number of electrons in the cluster and
  an electronic temperature of 1500 K, the on-site density matrix of that atom per \\( l \\) (summed over \\( m \\)) is built.
  The eigenvectors of its diagonalization, in order of decreasing occupation, become the contraction coefficients for \\( \zeta = 0, 1, \dots \\), with signs chosen so that the component of largest magnitude is positive
  (to lift degeneracies, \\( -10^{-9}(N\_l - p)^2 \\) is added to the diagonal).
- `Simple`: \\( a\_{l\zeta,p} = \delta\_{\zeta p} \\) (the \\( \zeta \\)-th primitive function itself).

### 5.5 Keywords and defaults [Public] [Convention]

| Keyword | Default | Meaning |
|---|---|---|
| `orbitalOpt.Method` | `Off` | `Off`, `Atoms`, `Species`, `Atoms2`, `Species2` |
| `orbitalOpt.Opt.Method` | `DIIS` | `DIIS`, `EF` |
| `orbitalOpt.InitCoes` | `Symmetrical` | initial coefficients |
| `orbitalOpt.SD.step` | 0.001 | steepest-descent step \\( s \\) |
| `orbitalOpt.HistoryPulay` | 20 | number of history entries; the manual says 15 |
| `orbitalOpt.StartPulay` | 10 | step at which quasi-Newton starts; the manual says 1 |
| `orbitalOpt.scf.maxIter` | 40 | SCF iterations per step |
| `orbitalOpt.Opt.maxIter` | 100 | maximum number of steps |
| `orbitalOpt.per.MDIter` | 1000000 | interval for redoing the optimization in MD |
| `orbitalOpt.criterion` | \\( 10^{-4} \\) | threshold on the gradient norm (displayed as (Hartree/bohr)²) |
| `CntOrb.fileout` | off | output of contracted PAOs (only for `Atoms` or `Species`) |
| `Num.CntOrb.Atoms` / `<Atoms.Cont.Orbitals` | 1 / — | number and indices of atoms to output |

### 5.6 Output of contracted PAO files [Public] [Convention]

[Public] With `CntOrb.fileout on`, `{species name}_{atom index}.pao` is written for each specified atom, and can be used directly as a basis (manual 26).

[Convention] The file structure is as follows.

1. Description lines (original PAO file name, basis specification, input file name).
2. `<Contraction.coefficients` ... `Contraction.coefficients>`: each line `Atom= L= Mul= p= coefficient` (only the \\( m = 0 \\) component; \\( p \\) beyond the number of primitives is 0).
3. A copy of the original PAO file up to just before `<pseudo.atomic.orbitals.L=0`.
4. The radial-function tables for the \\( l \\) used in the basis: the first \\( N^{\mathrm{cont}}\_l \\) columns are the contracted \\( \sum\_p a\_{l\zeta,p}R\_{lp}(r) \\),
   and the remaining columns (\\( \zeta \ge N^{\mathrm{cont}}\_l \\)) are the original primitive functions \\( R\_{l\zeta} \\) themselves (not orthogonalized against the contracted functions).
5. Tables for \\( l \\) not used in the basis (\\( l \\) beyond \\( l\_{\max} \\)) are kept as is.

Since the reader ([PAO basis](pao-basis.md)) uses the first \\( N\_l \\) functions of the table, specifying \\( N\_l \\) at most the number of contracted functions gives the optimized basis.

---

## 6. Edge cases [Convention] [Design]

- **Not enough states in the window**: in MLWF, a k point with fewer than \\( N\_W \\) bands in the outer window, and in CWF, a polynomial-window band range smaller than \\( N\_{\mathrm{CWF}} \\), are both errors.
  In CWF, if the Fermi window is narrow and \\( \delta \\) is small, singular values become of order \\( \delta \\); the CWFs keep their symmetry but the interpolation quality degrades (the Si example in technical note Sec. IV A).
- **Guiding functions not matching the band symmetry**: \\( F \\) becomes large and the interpolated bands are distorted. With two temperatures \\( k\_BT\_0 \ne k\_BT\_1 \\), bands to be discarded can be treated sharply and bands to be included smoothly (technical note Fig. 6).
- **Selecting AOs with nearly zero occupation in HO**: the eigenvalues of the on-site density matrix become degenerate near 0 and may mix with directions of unselected AOs. [Design] A reimplementation diagonalizes the matrix restricted to the subspace of the selected AOs.
- **Even number of k points for CWF**: the range of \\( \mathbf L \\) for the hoppings is not symmetric, so the existing implementation accepts only odd numbers.
- **Mixing different species in a group in unfolding**: the number and ordering of orbitals do not match, and the formula is undefined. [Design] Make this an error in input validation. Put impurities in separate groups (manual 55.2).
- **Atoms on the boundary of the reference cell**: a slight difference in the origin changes \\( \mathbf r(\cdot) \\), and for defective systems the weights change. For perfect crystals the result does not depend on the origin.
- **Negative orbital weights in a nonorthogonal basis**: the orbitally decomposed weights can be negative (the total weight and the sum rule are correct). The existing implementation replaces them by 0 only at output.
- **Number of contracted functions equal to the number of primitives in orbital optimization**: `s1>1` is the same as `s1`. The degrees of freedom to optimize are \\( N^{\mathrm{prim}} - 1 \\) excluding normalization.

## 7. Design: requirements for a modular reimplementation [Design]

- **Common input**: per-k-point \\( \varepsilon\_{\mu\mathbf k} \\), \\( C(\mathbf k) \\), \\( S(\mathbf k) \\), the chemical potential, real-space \\( S(\mathbf L) \\), \\( H(\mathbf L) \\), and the density matrix.
  It should be recomputable from the [sparse block format](../design/sparse-block-format.md) and the k-point results, independently of the internal state of the SCF.
- **Wannier core**: split into (a) building the projection matrix \\( A(\mathbf k) \\) (with interchangeable guiding functions, projectors, and window functions), (b) polar decomposition (CWF),
  (c) \\( M(\mathbf k,\mathbf b) \\) and MV/SMV minimization (MLWF), and (d) building hoppings, interpolation, and occupations from \\( U(\mathbf k) \\).
  (d) is shared by MLWF and CWF and takes the set of \\( \mathbf L \\) (a symmetric box, or a Wigner–Seitz supercell with degeneracies) as an argument.
- **Window functions** are kept as pure functions, with the polynomial window and the Fermi window (product form) as defaults, and the \\( f\_0 + f\_1 - 1 \\) form of the technical note also selectable.
- **Integrals for \\( M \\)**: the existing implementation uses sums over the real-space grid, but a reimplementation may evaluate them as two-center integrals (the Fourier method of [two-center integrals](two-center-overlap-kinetic.md) with the phase \\( e^{-i\mathbf b\cdot\mathbf r} \\) added). The difference is absorbed by the fixture tolerances.
- **File formats**: `.mmn`, `.amn`, `.eig` follow the Wannier90 formats, and a test checks that Wannier90 produces the same results.
- **Unfolding** is split into the preprocessing "reference cell, groups, origin -> per-atom \\( \mathbf r(\cdot) \\) and group" and the computation of "weights per k".
  The weights are O(N) by scanning the sparse-block \\( S(\mathbf L) \\) over atom pairs and can be parallelized independently per state. The origin defaults to the rule of the existing implementation; explicit specification is recommended.
- **NAO** is pure linear algebra "extract atomic blocks -> pre-NAO -> OWSO -> NAO", with the input \\( P \\), \\( S \\) given either from the whole system or from a cluster (the O(N) version differs only in the part that builds the cluster).
- **Orbital optimization** is an optimization outside the SCF loop, split into the gradient computation (from \\( \rho \\), \\( E \\) and the primitive-basis \\( H \\), \\( S \\)) and the update rules (SD, DIIS, BFGS).
  Contraction is represented as a linear combination of the radial-function tables of the [PAO basis](pao-basis.md), and the contracted tables are passed directly to the two-center integral and grid computations.
- **Handling of defaults**: defaults that differ between the manual and the implementation (`Wannier.Kgrid`, the windows, `Wannier.Readin.Overlap.Matrix`, `CWF.disentangling.smearing.bound`, `orbitalOpt.HistoryPulay`, `orbitalOpt.StartPulay`) follow the implementation values in compatibility mode.

## 8. Tests [Design]

| Kind | Content | Tolerance (guideline) |
|---|---|---|
| analytic | polynomial window: \\( P(0)=1 \\), \\( P'(0)=0 \\), \\( P(1)=P'(1)=P''(1)=0 \\), continuity at interval edges | \\( 10^{-12} \\) |
| analytic | polar decomposition: \\( F[U] = N\_k^{-1}\sum(\sigma-1)^2 \\), and \\( F[B] \ge F[U] \\) for any partial unitary \\( B \\) | \\( 10^{-12} \\) |
| analytic | if \\( A \\) has full rank, \\( U = A(A^\dagger A)^{-1/2} \\) | \\( 10^{-8} \\) |
| invariant | CWF hoppings are invariant under random phases applied to the Bloch states | \\( 10^{-12} \\) |
| invariant | CWF hoppings are real for collinear, time-reversal-symmetric systems | \\( 10^{-10} \\) |
| analytic | all orbitals as guiding functions, \\( w \equiv 1 \\): the interpolated bands agree with the original bands at grid k points | \\( 10^{-10} \\) |
| invariant | MLWF: \\( \Omega\_{\mathrm I} \\) is independent of the gauge (random unitaries) | \\( 10^{-12} \\) |
| analytic | MLWF: for a complete set (point-like orthonormal orbitals), \\( \Omega\_{\mathrm I}=0 \\), minimization gives \\( \Omega\to 0 \\), and the centers coincide with the orbital positions | \\( \Omega<10^{-5} \\), centers \\( 10^{-3} \\) |
| analytic | the b-vector weights satisfy \\( \sum w\_b b\_\alpha b\_\beta = \delta\_{\alpha\beta} \\) (simple cubic, fcc, hexagonal) | \\( 10^{-10} \\) |
| invariant | unfolding sum rule \\( \sum\_{\mathbf G}W/n = 1 \\) (with defects and displacements) | \\( 10^{-10} \\) |
| analytic | perfect supercell: weights 0/1, energies of weight-1 states agree with the reference-cell bands, independent of the origin | \\( 10^{-10} \\) |
| analytic | reference cell = supercell (`Unfolding.ReferenceVectors` omitted): the total weight of every state is 1 | \\( 10^{-10} \\) |
| analytic | the orbital-optimization gradient agrees with the difference quotient of the total energy (or of the band energy with fixed \\( H \\)) | relative \\( 10^{-6} \\) |
| analytic | normalization of the contraction coefficients \\( a^{\mathsf T}S^{(I)}a = 1 \\); the contracted columns of the output `.pao` are \\( \sum\_p a\_pR\_p \\) | \\( 10^{-12} \\) |
| cross-check | Si (diamond) MLWF (manual example, 8 functions): \\( \Omega\_{\mathrm I} \\), total spread, centers | spread 1%, centers \\( 10^{-3} \\) Å |
| cross-check | CWF of GaAs (conditions of manual 47.1.1): distance function, hoppings, CWF occupations | hoppings \\( 10^{-4} \\) eV, occupations \\( 10^{-5} \\) |
| cross-check | unfolding weights of a 2×2 supercell of 2D SiC (perfect and with a Si vacancy) (manual 55.2) | \\( 10^{-6} \\) |
| cross-check | NP of ethylene carbonate (manual 36 example) | \\( 10^{-4} \\) electrons |
| cross-check | orbital optimization of methane (manual 26 example): the gradient-norm history of the first few steps, final total energy | gradient norm relative \\( 10^{-3} \\), energy \\( 10^{-6} \\) Hartree |
| mutual | the spreads and interpolated bands obtained by feeding the MLWF/CWF `.mmn`, `.amn`, `.eig` to Wannier90 agree with the results of the main program | spread 1% |

Check: `uv run --with numpy --with scipy python tools/check_wannier_unfolding.py` verifies, on a 1D Gaussian-orbital model (nonorthogonal),
the unfolding sum rule, 0/1 weights for perfect crystals, and origin independence; that the sum rule fails with the conjugate placed as in the original paper;
minimality, gauge invariance, agreement with Löwdin, and interpolation of the CWF polar decomposition; smoothness of the polynomial window; invariance and minimization of \\( \Omega\_{\mathrm I} \\) for MLWF;
and the orbital-optimization gradient (all pass).

## [Open]

- Whether the MLWF calculation can be used with `scf.EigenvalueSolver Cluster` (isolated systems), and if so, how b vectors are handled.
- Whether the details of MLWF disentanglement (construction of the \\( Z \\) matrix, treatment of frozen states, choice of the initial subspace) agree with Wannier90. Confirmation with fixtures is needed.
- How the Wigner–Seitz supercell of `.HWR` is constructed (whether the range and degeneracy determination are the same as in Wannier90).
- Details of the secant method in MLWF CG and of the SD->CG switch in scheme 2 (step acceptance conditions).
- Whether to inherit, as the reimplementation default, the existing rule of determining the band-index range for the CWF polynomial window from the eigenvalues at a single k point (using all bands inside the window at each k point would be more natural).
- How HO and MO guiding functions (which may be complex) are built in noncollinear CWF calculations, and how spin indices are handled.
- Whether the normalization (k-point weights) in the output of CWF effective charges is the same for Cluster and Band.
- NAO: the intent of the rule halving the threshold for spin polarization (since the maximum per-spin occupation is 1, 0.425 may be too low).
- The exact extent of the cluster in the O(N) version of NAO (`on2`) (its relation to the truncated cluster defined in the Krylov-method chapter), and the processing of NHO and NBO (`NAO.only off`). This chapter covers only NAO and NP.
- Orbital optimization: the meaning of `orbitalOpt.Force.Skip`; the intent of how initial coefficients are built for `Atoms2` and `Species2` (placing the eigenvector of largest occupation into the last \\( \zeta \\));
  the unit displayed for `orbitalOpt.criterion` (displayed as (Hartree/bohr)², but actually the sum of squares of the gradient with respect to the contraction coefficients).
- For the automatic determination of the unfolding origin, which to adopt as the reimplementation default: the explanation in manual 55.3 (determination by search) or the implementation in the version checked (fixed rule).
