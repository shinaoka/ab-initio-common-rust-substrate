# Energy, forces and stress for noncollinear and spin-orbit calculations

> Status: first version (awaiting review). This extends what [Total energy and forces](energy-force.md) wrote only for collinear spin (real matrices)
> to noncollinear spin (`scf.SpinPolarization NC`) and spin-orbit coupling (`scf.SpinOrbit.Coupling on`).
> The equations re-derive, in the notation of this note, technical note 2 (Ver. 1.0, hereafter TN2), technical note 1 (Ver. 1.3, hereafter TN1), Kübler et al. (J. Phys. F **18**, 469 (1988)) and
> manual 4.0 (Secs. 37, 38, 39, 41, 43).
> Rules not written in the public documents are supplemented as [Convention] from the behavior of the existing implementation.
> Algebraic identities are checked numerically with `tools/check_noncollinear_energy_force.py`.

## Scope of this chapter [Design]

- Covered: the spin 2×2 density matrix and energy-weighted density matrix, evaluation of each energy component, forces, availability of stress,
  contribution of the noncollinear exchange-correlation potential to forces, energies of spin-orientation constraints and Zeeman terms,
  input of the initial spin configuration, and output of local spin and orbital moments.
- Parts identical to the collinear case (screened interatomic terms, atomic-density correction, grid treatment, sign and units) are not repeated;
  see [Total energy and forces](energy-force.md).
- The equations for the noncollinear grid density and XC potential themselves are handled in Secs. 2.4 and 4.6 of [Real-space grid, Poisson, exchange-correlation](grid-poisson-xc.md),
  the spin-orbit matrix elements in Sec. 5 of [Nonlocal pseudopotential](nonlocal-pseudopotential.md),
  the ordering of spinor matrices and Bloch sums in [k-points and Bloch sums](kpoints-bloch.md), and the storage of the density matrix and initial density in
  [Density matrix and charge mixing](density-mixing.md).
- Details of DFT+U (`Uhub`), spin-orientation constraints, Zeeman terms, orbital moments and initial magnetic configuration are handled in [DFT+U, constrained DFT and magnetic features](dftu-magnetism.md).
  Here only what is needed from the viewpoint of energy and force contractions is summarized, and the symbols of the equations are aligned.

## Symbols [Design]

In addition to the symbols of [Total energy and forces](energy-force.md), the following are used.

| Symbol | Meaning |
|---|---|
| \\( \sigma, \sigma' \in \\{\uparrow,\downarrow\\} \\) | Spin with the global \\( z \\) axis as quantization axis |
| \\( c\_{\sigma I\alpha,\nu}(\mathbf k) \\) | Spinor eigenvector (length \\( 2N \\); first half \\( \uparrow \\), second half \\( \downarrow \\)) |
| \\( \rho^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) \\) | Spin 2×2 real-space density matrix (complex) |
| \\( \varepsilon^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) \\) | Energy-weighted density matrix of the same form |
| \\( X^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) \\) | Matrix element \\( \langle I\alpha\sigma,\mathbf 0\rvert\hat X\lvert J\beta\sigma',\mathbf L\rangle \\) of a spinor operator \\( \hat X \\) |
| \\( \underline n(\mathbf r) \\) | 2×2 spin density on the grid (components \\( n\_{\sigma\sigma'} \\)) |
| \\( \mathbf m(\mathbf r) \\), \\( \hat{\mathbf m} \\) | Spin density vector and its direction \\( (\sin\theta\cos\varphi, \sin\theta\sin\varphi, \cos\theta) \\) |
| \\( n'\_\uparrow, n'\_\downarrow \\) | Densities in the local spin frame, \\( \frac12(n \pm \lVert\mathbf m\rVert) \\) |
| \\( \underline V(\mathbf r) \\) | 2×2 effective potential on the grid |
| \\( \boldsymbol\sigma = (\sigma\_x,\sigma\_y,\sigma\_z) \\), \\( \sigma\_0 \\) | Pauli matrices and the 2×2 identity |
| \\( \langle \rho, X\rangle \\) | Contraction defined in Sec. 1.4 |

"Direction", as in manual Sec. 37, refers to the direction of the spin density \\( \mathbf m \\), not of the electron's magnetic moment
(the physical magnetic moment points the opposite way).

## 1. Spin 2×2 density matrix

### 1.1 Definition [Public] [Convention]

[Public] TN2 Eq. 1–3: the density operator is built from the occupation sum of spinors \\( \psi\_\nu = (\psi\_{\nu\uparrow}, \psi\_{\nu\downarrow}) \\).
With the phase convention of this note ([sparse block format](../design/sparse-block-format.md), [k-points and Bloch sums](kpoints-bloch.md)), the real-space density matrix is defined as

$$
\rho^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k}\sum\_\nu f\_\nu(\mathbf k)\\, e^{+i\mathbf k\cdot\mathbf L}\\,
\overline{c\_{\sigma I\alpha,\nu}(\mathbf k)}\\, c\_{\sigma' J\beta,\nu}(\mathbf k)
$$

$$
\varepsilon^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k}\sum\_\nu f\_\nu(\mathbf k)\\, e\_\nu(\mathbf k)\\, e^{+i\mathbf k\cdot\mathbf L}\\,
\overline{c\_{\sigma I\alpha,\nu}(\mathbf k)}\\, c\_{\sigma' J\beta,\nu}(\mathbf k)
$$

\\( f\_\nu \le 1 \\) (one electron per spinor), and the factor corresponding to the collinear \\( g\_s \\) is 1.

[Convention] The density matrix of the existing implementation is exactly this expression (order \\( \bar c\\, c \\), phase \\( e^{+i\mathbf k\cdot\mathbf L} \\), real and imaginary parts stored separately).
From Hermiticity

$$
\rho^{\sigma'\sigma}\_{J\beta,I\alpha}(-\mathbf L) = \overline{\rho^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L)}
$$

\\( \rho^{\downarrow\uparrow} \\) is determined by \\( \rho^{\uparrow\downarrow} \\), and the stored components are, as in "Noncollinear case" of
[Density matrix and charge mixing](density-mixing.md), \\( \mathrm{Re}\\,\rho^{\uparrow\uparrow} \\), \\( \mathrm{Re}\\,\rho^{\downarrow\downarrow} \\),
\\( \mathrm{Re}\\,\rho^{\uparrow\downarrow} \\), \\( \mathrm{Im}\\,\rho^{\uparrow\downarrow} \\), plus \\( \mathrm{Im}\\,\rho^{\uparrow\uparrow} \\), \\( \mathrm{Im}\\,\rho^{\downarrow\downarrow} \\) when time-reversal-breaking terms are present.

### 1.2 Time reversal and imaginary parts [Public] (derived)

Since the basis functions are real, the time-reversal operator acts as \\( \Theta = (i\sigma\_y)K \\) (\\( K \\) is complex conjugation) without changing orbital indices.
When the Hamiltonian is time-reversal symmetric and Kramers pairs are occupied together (nonmagnetic), \\( \Theta\rho\Theta^{-1} = \rho \\) and Hermiticity give

$$
\rho^{\uparrow\uparrow}\_{ab}(\mathbf L) = \overline{\rho^{\downarrow\downarrow}\_{ab}(\mathbf L)},\qquad
\rho^{\uparrow\downarrow}\_{ab}(\mathbf L) = -\rho^{\uparrow\downarrow}\_{ba}(-\mathbf L)
$$

(\\( a = I\alpha \\), \\( b = J\beta \\)). That is, even in the nonmagnetic case, spin-orbit coupling gives rise to
\\( \mathrm{Im}\\,\rho^{\uparrow\uparrow} = -\mathrm{Im}\\,\rho^{\downarrow\downarrow} \ne 0 \\) and to a \\( \rho^{\uparrow\downarrow} \\) antisymmetric in the orbitals.
Without spin-orbit coupling, and if all spin-dependent terms are products of orbitally real symmetric matrices and spin matrices,
\\( \rho^{\sigma\sigma'}\_{ab}(\mathbf L) = \rho^{\sigma\sigma'}\_{ba}(-\mathbf L) \\) (symmetric in the orbitals).
Both are checked in checks 5 and 6 of `tools/check_noncollinear_energy_force.py`.

### 1.3 Correspondence with the spin density on the grid [Public] [Convention]

[Public] TN2 Eq. 3, 33–35: the 2×2 density on the grid is \\( n\_{\sigma\sigma'}(\mathbf r) = \sum\_\nu f\_\nu\\, \psi\_{\nu\sigma}(\mathbf r)\\,\overline{\psi\_{\nu\sigma'}(\mathbf r)} \\),
\\( \underline n = \frac12(n\\,\sigma\_0 + \mathbf m\cdot\boldsymbol\sigma) \\) (same as Sec. 2.4 of [grid](grid-poisson-xc.md); \\( m\_y = -2\\,\mathrm{Im}\\,n\_{\uparrow\downarrow} \\)).

Written in terms of \\( \rho \\) of Sec. 1.1, the spin indices are swapped.

$$
n\_{\sigma\sigma'}(\mathbf r) = \sum\_{\mathbf R}\sum\_{\mathbf L}\sum\_{ab} \rho^{\sigma'\sigma}\_{ab}(\mathbf L)\\,
\phi\_a(\mathbf r - \mathbf R\_I - \mathbf R)\\,\phi\_b(\mathbf r - \mathbf R\_J - \mathbf L - \mathbf R)
$$

(\\( \mathbf R \\) are periodic images; the sum is over the atom-pair range). The diagonal components are determined only by the real part of \\( \rho^{\sigma\sigma} \\) (the imaginary part is antisymmetric under \\( a\leftrightarrow b \\) and cancels), and
the off-diagonal component is \\( n\_{\uparrow\downarrow} = \overline{\sum \rho^{\uparrow\downarrow}\phi\phi} \\).

[Convention] The existing implementation builds \\( \sum\rho^{\uparrow\downarrow}\phi\phi \\) and then flips the sign of the imaginary part to obtain \\( n\_{\uparrow\downarrow} \\). This is the same as the equation above.

### 1.4 Contraction rule [Public] (derived) [Convention]

For the matrix elements \\( X^{\sigma\sigma'}\_{ab}(\mathbf L) \\) of a spinor operator \\( \hat X \\) (with Bloch sum \\( X(\mathbf k) = \sum\_{\mathbf L} e^{+i\mathbf k\cdot\mathbf L}X(\mathbf L) \\)),

$$
\sum\_{\mathbf k} w\_{\mathbf k}\sum\_\nu f\_\nu\\, c\_\nu^\dagger X(\mathbf k)\\, c\_\nu
{}= \sum\_{\sigma\sigma'}\sum\_{\mathbf L}\sum\_{ab} \rho^{\sigma\sigma'}\_{ab}(\mathbf L)\\, X^{\sigma\sigma'}\_{ab}(\mathbf L)
{}\equiv \langle\rho, X\rangle
$$

holds (note that this is a sum of products without complex conjugation). If \\( X \\) is Hermitian, the value is real, and combining \\( \downarrow\uparrow \\) into \\( \uparrow\downarrow \\) by Hermiticity gives

$$
\langle\rho, X\rangle = \sum\_{\mathbf L, ab}\Bigl[
\sum\_{\sigma}\bigl(\mathrm{Re}\\,\rho^{\sigma\sigma}\\,\mathrm{Re}\\,X^{\sigma\sigma} - \mathrm{Im}\\,\rho^{\sigma\sigma}\\,\mathrm{Im}\\,X^{\sigma\sigma}\bigr)
{}+ 2\bigl(\mathrm{Re}\\,\rho^{\uparrow\downarrow}\\,\mathrm{Re}\\,X^{\uparrow\downarrow} - \mathrm{Im}\\,\rho^{\uparrow\downarrow}\\,\mathrm{Im}\\,X^{\uparrow\downarrow}\bigr)\Bigr]
$$

(the sum is over all stored atom pairs and lattice vectors). For spin-independent real symmetric \\( X = O\otimes\sigma\_0 \\) (\\( T \\), \\( S \\), \\( V^{\mathrm{na}} \\), scalar \\( V^{\mathrm{NL}} \\)), this simplifies to

$$
\langle\rho, O\otimes\sigma\_0\rangle = \sum\_{\mathbf L,ab}\mathrm{Re}\bigl(\rho^{\uparrow\uparrow}\_{ab} + \rho^{\downarrow\downarrow}\_{ab}\bigr)(\mathbf L)\\, O\_{ab}(\mathbf L)
$$

(the contribution of \\( \mathrm{Im}\\,\rho^{\sigma\sigma} \\) vanishes by Hermiticity).

[Convention] All energy and force contractions of the existing implementation take one of these two forms (`Re·Re − Im·Im`, with a factor 2 on the spin-off-diagonal part).
Numerical check: check 1 (agreement with the band energy, agreement of the decomposed form).

## 2. Energy components

### 2.1 List [Public] [Convention]

Output labels and their order are the same as in the collinear case ("Output energy components" in [Total energy and forces](energy-force.md)). The noncollinear definitions are as follows.

| Label | Noncollinear definition | Remarks |
|---|---|---|
| `Uele` | \\( \sum\_{\mathbf k}w\_{\mathbf k}\sum\_\nu f\_\nu e\_\nu \\) | Sum of spinor eigenvalues (factor 1) |
| `Ukin` | \\( \langle\rho, T\otimes\sigma\_0\rangle \\) | Only \\( \mathrm{Re}(\rho^{\uparrow\uparrow}+\rho^{\downarrow\downarrow}) \\) |
| `Una` | \\( \langle\rho, V^{\mathrm{na}}\otimes\sigma\_0\rangle \\), or \\( \Delta V\sum\_p n\\, V^{\mathrm{na}} \\) | \\( n = n\_{\uparrow\uparrow}+n\_{\downarrow\downarrow} \\) |
| `Unl` | \\( \langle\rho, V^{\mathrm{NL}}\rangle \\) | With spin-orbit, all four blocks (Sec. 2.2) |
| `UH0`, `Ucore` | Same as collinear | Functions of atomic positions only |
| `UH1` | \\( \frac12\Delta V\sum\_p (n - n^{a})\\,\delta V\_H \\) | \\( n \\) is the total density |
| `Uxc0`, `Uxc1` | Majority and minority spin parts in the local spin frame (Sec. 2.3) | Not the global-frame \\( \uparrow\downarrow \\) |
| `Uhub` | DFT+U | [DFT+U, constrained DFT and magnetic features](dftu-magnetism.md) |
| `Ucs` | Spin-orientation constraint (Sec. 2.4) | `scf.Constraint.NC.Spin` |
| `Uzs`, `Uzo` | Zeeman terms of spin and orbital moments (Sec. 2.4) | `scf.NC.Zeeman.Spin`, `scf.NC.Zeeman.Orbital` |
| `Uef`, `UvdW`, `Uch` | Same as collinear (electric field, dispersion, core hole) | The `Uch` matrix uses the general contraction of Sec. 1.4 |

`Utot` is the same sum as in the collinear case.

### 2.2 Nonlocal term and spin-orbit coupling [Public] [Convention]

[Public] The matrix elements \\( V^{\mathrm{NL},\sigma\sigma'}\_{ab}(\mathbf L) \\) of the j-dependent pseudopotential (TN2 Eq. 40–58) are built with the equations of
Sec. 5.3 of [Nonlocal pseudopotential](nonlocal-pseudopotential.md). The imaginary part of the diagonal spin blocks is the \\( L\_z \\) term, and
the \\( \uparrow\downarrow \\) block is the \\( L\_- \\) term.

[Convention]

- `Unl` is the general form \\( \langle\rho, V^{\mathrm{NL}}\rangle \\) of Sec. 1.4 and includes the spin-orbit part (the imaginary part of the diagonal spin blocks and the \\( \uparrow\downarrow \\) block).
- Only the nonlocal part of the pseudopotential enters `Unl`. The potentials of DFT+U, constraints and Zeeman terms are
  added to the same imaginary-part matrices in the Hamiltonian, but in the energy they are evaluated separately as `Uhub`, `Ucs`, `Uzs`, `Uzo`.
- In noncollinear calculations without spin-orbit coupling (`scf.SpinOrbit.Coupling off`, or no j-dependent files),
  \\( V^{\mathrm{NL}} \\) is scalar (j-averaged, Sec. 4 of the nonlocal chapter), and the energy is \\( \langle\rho, V^{\mathrm{NL}}\otimes\sigma\_0\rangle \\).

### 2.3 Exchange-correlation: local spin frame [Public] [Convention]

[Public] Kübler et al. (1988) and TN2 Eq. 4–5, 17–25: the exchange-correlation energy is obtained by applying the collinear functional to the densities \\( n'\_\uparrow, n'\_\downarrow \\)
in the local frame that diagonalizes the 2×2 density at each point. Including the partial core correction,

$$
E\_{\mathrm{xc}} = \int d^3r\\; \bigl(n + n^{\mathrm{pcc}}\bigr)\\,
\epsilon\_{\mathrm{xc}}\bigl(n'\_\uparrow + n^{\mathrm{pcc}}\_\uparrow,\\; n'\_\downarrow + n^{\mathrm{pcc}}\_\downarrow;\\; \nabla n'\_\uparrow, \nabla n'\_\downarrow\bigr)
$$

(the treatment of GGA gradients is in Sec. 4.6 of [grid](grid-poisson-xc.md): \\( n'\_\sigma \\) are regarded as scalar fields on the grid).
Since \\( n'\_\sigma \\) are eigenvalues of \\( \underline n \\), they are invariant under global spin rotations.

[Convention]

- `Uxc0` and `Uxc1` are the collinear expressions ("Terms evaluated on the grid" in [Total energy and forces](energy-force.md)) with \\( n\_\uparrow, n\_\downarrow \\)
  replaced by the local-frame \\( n'\_\uparrow, n'\_\downarrow \\) (\\( n'\_\uparrow \ge n'\_\downarrow \\)):
  \\( U\_{\mathrm{xc}0} = \Delta V\sum\_p (n'\_\uparrow + n^{\mathrm{pcc}}\_\uparrow)\\,\epsilon\_{\mathrm{xc}} \\) (+ half of the atomic-density correction), and `Uxc1` likewise with \\( n'\_\downarrow \\).
- \\( n^{\mathrm{pcc}}\_\uparrow = n^{\mathrm{pcc}}\_\downarrow = \frac12 n^{\mathrm{pcc}} \\). For species with \\( +1 \\) (\\( -1 \\)) in `<scf.pcc.opencore`,
  the whole is added to the majority (minority) spin of the local frame.
- The atomic-density correction \\( \Delta E^{(0)} \\) is spin-independent, so as in the collinear case it is added half to each.
- Therefore `Uxc0` and `Uxc1` can be individually compared with a collinear calculation only when \\( \mathbf m \\) points along \\( +z \\) at all points (\\( n\_{\uparrow\uparrow} \ge n\_{\downarrow\downarrow} \\)).
  The sum `Uxc0 + Uxc1` can always be compared.
- Points with \\( \lVert\mathbf m\rVert < 10^{-14} \\) are handled as in Sec. 2.4 of [grid](grid-poisson-xc.md).

### 2.4 Constraints and Zeeman terms [Public] [Convention]

[Public] Manual Sec. 41: `scf.Constraint.NC.Spin` (`on`, `on2`, `off`; default `off`) adds a penalty term that raises the energy when the spin direction deviates
from the initial direction. `on` constrains only the direction; `on2` constrains both direction and magnitude.
The strength is `scf.Constraint.NC.Spin.v` (eV, default 0), and the target atoms are given by column 12 of `Atoms.SpeciesAndCoordinates` (1 means target).
Manual Sec. 43 and TN2 Eq. 139–150: `scf.NC.Zeeman.Spin` and `scf.NC.Mag.Field.Spin` (Tesla) add a Zeeman term for the spin moment, and
`scf.NC.Zeeman.Orbital` and `scf.NC.Mag.Field.Orbital` add a Zeeman term for the orbital moment.
The field direction is given per atom by columns 8 and 9 (spin) and columns 10 and 11 (orbital), and acts only on atoms whose column 12 is 1.

[Convention] Details are in Secs. 5–7 of [DFT+U, constrained DFT and magnetic features](dftu-magnetism.md); they are summarized here. These terms are defined with the per-atom spin 2×2 occupation \\( \underline N\_I \\) (same orientation convention as the grid density; \\( T^{(I)} \\) of that chapter).

$$
A^{\sigma\sigma'}\_I = \sum\_{i\in I}\sum\_{J,\mathbf L}\sum\_{j} \rho^{\sigma\sigma'}\_{Ii,Jj}(\mathbf L)\\, S\_{Ii,Jj}(\mathbf L),\qquad
\underline N\_I = \tfrac12\bigl(\overline{\underline A\_I} + \underline A\_I^{\mathsf T}\bigr)
$$

This is the DFT+U 'dual' occupation matrix traced over orbitals, a Hermitian 2×2 matrix
(definition of the 'dual' occupation: Han, Ozaki, Yu, PRB **74**, 045110 (2006)). Let the eigenvalues of \\( \underline N\_I \\) be \\( N^\uparrow\_I \ge N^\downarrow\_I \\), and
\\( \mathbf s\_I = (\mathrm{Re}\\,N^{\uparrow\downarrow}\_I,\ -\mathrm{Im}\\,N^{\uparrow\downarrow}\_I,\ \frac12(N^{\uparrow\uparrow}\_I - N^{\downarrow\downarrow}\_I)) \\) (spin with \\( \hbar = 1 \\)).
\\( \hat{\mathbf n}^0\_I \\) is the unit vector of the input angles (after applying the reversal rule of Sec. 5).

- Constraint (`on`):

$$
U\_{\mathrm{cs}} = v\sum\_{I}\bigl\lVert \underline N\_I - \underline N^0\_I\bigr\rVert\_F^2,\qquad
\underline N^0\_I = \tfrac12\Bigl[(N^\uparrow\_I + N^\downarrow\_I)\\,\sigma\_0 + (N^\uparrow\_I - N^\downarrow\_I)\\,\hat{\mathbf n}^0\_I\cdot\boldsymbol\sigma\Bigr]
$$

  (\\( \lVert\cdot\rVert\_F \\) is the Frobenius norm, and \\( v \\) is `scf.Constraint.NC.Spin.v` converted to Hartree).
  The target \\( \underline N^0\_I \\) keeps the current magnitude and sets only the direction to the input direction.
  The potential also includes in the derivative the dependence of \\( \underline N^0\_I \\) on \\( \underline N\_I \\) (through the eigenvalues).
- Constraint (`on2`): \\( N^\uparrow\_I - N^\downarrow\_I \\) is replaced by the input \\( \lvert N^{\uparrow}\_{I,\mathrm{in}} - N^{\downarrow}\_{I,\mathrm{in}}\rvert \\) (columns 6 and 7) (the total \\( N^\uparrow\_I + N^\downarrow\_I \\) is the current value).
- Spin Zeeman term: \\( U\_{\mathrm{zs}} = -B\_s\sum\_I \hat{\mathbf n}^0\_I\cdot\mathbf s\_I \\). \\( B\_s \\) is `scf.NC.Mag.Field.Spin` converted with
  \\( 1\ \mathrm{T} = 1/(2.35051742\times10^{5}) \\) a.u. The sign is such that the energy decreases when the spin aligns with \\( \hat{\mathbf n}^0\_I \\).
- Orbital Zeeman term: \\( U\_{\mathrm{zo}} = -\frac12 B\_o\sum\_I \hat{\mathbf n}^{0}\_{o,I}\cdot\mathbf M\_I \\). \\( \mathbf M\_I \\) is the orbital-moment vector of Sec. 6.2 (\\( \mathbf L\_I \\) of that chapter),
  \\( B\_o \\) is `scf.NC.Mag.Field.Orbital` converted to a.u. in the same way, and \\( \hat{\mathbf n}^{0}\_{o,I} \\) is the direction from columns 10 and 11.
- Input conditions: all of these are available only in noncollinear calculations and require `scf.Hubbard.Occupation` to be `dual` (default).
  The constraint and the spin Zeeman term cannot be used together. The orbital Zeeman term requires spin-orbit coupling.
  Inputs that do not satisfy these conditions stop with an error.

[Public] The spin component \\( s\_y = \frac{i}{2}(N\_{\alpha\beta} - N\_{\beta\alpha}) = -\mathrm{Im}\\,N\_{\alpha\beta} \\) of TN2 Eq. 148–150 agrees with \\( \mathbf s\_I \\) above.

### 2.5 Band-energy identity [Public] (derived)

From Sec. 1.4, for the same eigenvectors, even before convergence,

$$
U\_{\mathrm{ele}} = \bigl\langle \rho,\ T\otimes\sigma\_0 + V^{\mathrm{na}}\otimes\sigma\_0 + V^{\mathrm{NL}} + H^{\mathrm{grid}} + V^{\mathrm{U}} \bigr\rangle,\qquad
H^{\mathrm{grid},\sigma\sigma'}\_{ab}(\mathbf L) = \Delta V\sum\_p \phi\_a\\, V\_{\sigma\sigma'}(\mathbf r\_p)\\,\phi\_b
$$

holds (\\( V^{\mathrm{U}} \\) is the matrix of DFT+U, constraints and Zeeman terms; \\( \underline V = (\delta V\_H\\,[+V^{\mathrm{na}}])\sigma\_0 + \underline V\_{\mathrm{xc}} \\)).
By the index swap of Sec. 1.3, the grid term is \\( \langle\rho, H^{\mathrm{grid}}\rangle = \Delta V\sum\_p \mathrm{Tr}[\underline V\\,\underline n] \\) (check 3).
The double-counting form of TN2 Eq. 6, \\( E = E\_{\mathrm{band}} - \frac12\int n V\_H - \int\mathrm{Tr}(\underline V\_{\mathrm{xc}}\underline n) + E\_{\mathrm{xc}} \\), is a restatement of this identity.

## 3. Forces [Public] [Convention]

### 3.1 Overall form [Public] (derived)

\\( \mathbf F\_K = -\partial U\_{\mathrm{tot}}/\partial\mathbf R\_K \\) (sign, units and output are the same as collinear). For the spinor generalized eigenvalue problem
\\( H^{(2)}c = e\\,S^{(2)}c \\), \\( S^{(2)} = \sigma\_0\otimes S \\), at fixed occupations

$$
\partial\sum\_\nu f\_\nu e\_\nu = \langle\rho, \partial H\rangle - \langle\varepsilon, \partial S^{(2)}\rangle,\qquad
\langle\varepsilon, \partial S^{(2)}\rangle = \sum\_{\mathbf L,ab}\mathrm{Re}\bigl(\varepsilon^{\uparrow\uparrow}\_{ab} + \varepsilon^{\downarrow\downarrow}\_{ab}\bigr)(\mathbf L)\\,\partial S\_{ab}(\mathbf L)
$$

(the spinor version of TN1 Eq. 71–75; check 2). Hence only \\( \mathrm{Re}(\varepsilon^{\uparrow\uparrow}+\varepsilon^{\downarrow\downarrow}) \\) is needed for the Pulay term.
[Convention] The existing implementation also uses only this sum for the Pulay term.

The decomposition of forces is the same as in the collinear case, with the "sum over spin \\( g\_s\sum\_\sigma \\)" of each term replaced by the contraction of Sec. 1.4.

| Term | Contribution to \\( \partial U\_{\mathrm{tot}}/\partial\mathbf R\_K \\) (noncollinear form) |
|---|---|
| Pulay | \\( -\langle\varepsilon, \partial S\otimes\sigma\_0\rangle \\) |
| Kinetic, projector-expansion \\( V^{\mathrm{na}} \\) | \\( \langle\rho, \partial O\otimes\sigma\_0\rangle \\), only \\( \mathrm{Re}(\rho^{\uparrow\uparrow}+\rho^{\downarrow\downarrow}) \\) |
| Nonlocal (including spin-orbit) | \\( \langle\rho, \partial V^{\mathrm{NL}}\rangle \\) (Sec. 3.2) |
| Grid: movement of basis functions | \\( \langle\rho, \partial H^{\mathrm{grid}}\rangle \\) (Sec. 3.3) |
| Grid: \\( n^{a} \\), \\( n^{\mathrm{pcc}} \\), grid-mode \\( V^{\mathrm{na}} \\), \\( \Delta E^{(0)} \\) | Same as collinear (\\( n^{\mathrm{pcc}} \\): Sec. 3.4) |
| `Ucore + UH0` | Same as collinear |
| `Ucs`, `Uzs`, `Uzo`, `Uhub` | Sec. 3.5 |

### 3.2 Nonlocal and spin-orbit forces [Public] (derived) [Convention]

The derivatives of the nonlocal matrix elements \\( \partial V^{\mathrm{NL},\sigma\sigma'}\_{ab}/\partial\mathbf R\_X \\) are built as complex numbers for each of the four spin blocks with the equations of Sec. 6 of [Nonlocal pseudopotential](nonlocal-pseudopotential.md)
(including the derivative of the damping function; \\( P \\) is the 2×2×orbital matrix of the j projection). The force is the contraction of these with \\( \rho \\)

$$
\frac{\partial U\_{\mathrm{nl}}}{\partial\mathbf R\_X} = \sum\_{\mathbf L,ab}\Bigl[\sum\_\sigma\bigl(\mathrm{Re}\\,\rho^{\sigma\sigma}\\,\mathrm{Re}\\,\partial V^{\sigma\sigma} - \mathrm{Im}\\,\rho^{\sigma\sigma}\\,\mathrm{Im}\\,\partial V^{\sigma\sigma}\bigr)
{}+ 2\bigl(\mathrm{Re}\\,\rho^{\uparrow\downarrow}\\,\mathrm{Re}\\,\partial V^{\uparrow\downarrow} - \mathrm{Im}\\,\rho^{\uparrow\downarrow}\\,\mathrm{Im}\\,\partial V^{\uparrow\downarrow}\bigr)\Bigr]
$$

The spin-orbit part of the force enters through \\( \mathrm{Im}\\,\rho^{\sigma\sigma} \\) and \\( \rho^{\uparrow\downarrow} \\).

[Convention] The existing implementation uses this general form when any of spin-orbit coupling, DFT+U, constraints or Zeeman terms is present,
and otherwise uses the same form as collinear with only \\( \mathrm{Re}(\rho^{\uparrow\uparrow}+\rho^{\downarrow\downarrow}) \\) (in the latter case it is mathematically identical to the general form).

### 3.3 Grid terms and the noncollinear XC potential [Public] (derived) [Convention]

The grid term due to movement of basis functions is the derivative of \\( \Delta V\sum\_p\mathrm{Tr}[\underline V\\,\underline n] \\) at fixed \\( \rho \\) and \\( \underline V \\).
The row and column contributions are complex conjugates of each other by Hermiticity, so

$$
\frac{\partial E}{\partial\mathbf R\_K}\Big|\_{\mathrm{basis}}
{}= -2\\,\mathrm{Re}\sum\_{\alpha\in K}\sum\_{J,\mathbf L,\beta}\sum\_{\sigma\sigma'} \rho^{\sigma\sigma'}\_{K\alpha,J\beta}(\mathbf L)\\;
\Delta V\sum\_p \nabla\phi\_{K\alpha}(\mathbf r\_p - \mathbf R\_K)\\, V\_{\sigma\sigma'}(\mathbf r\_p)\\, \phi\_{J\beta}(\mathbf r\_p - \mathbf R\_J - \mathbf L)
$$

(corresponding to the factor 2 of collinear TN1 Eq. 67, 79). Here \\( V\_{\sigma\sigma'} = \delta\_{\sigma\sigma'}(\delta V\_H\\,[+V^{\mathrm{na}}]\\,[+V\_{\mathrm{ef}}]) + V\_{\mathrm{xc},\sigma\sigma'} \\) and
\\( \underline V\_{\mathrm{xc}} = \bar v\\,\sigma\_0 + \Delta v\\,\hat{\mathbf m}\cdot\boldsymbol\sigma \\) (Sec. 4.6 of [grid](grid-poisson-xc.md)).

**Why the angles of the local frame need not be differentiated.** \\( E\_{\mathrm{xc}} \\) depends on \\( \underline n \\) only through the eigenvalues \\( n'\_\pm \\).
The derivative of the eigenvalues of a 2×2 Hermitian matrix is \\( \partial n'\_\pm/\partial n\_{\sigma'\sigma} = (P\_\pm)\_{\sigma\sigma'} \\) with \\( P\_\pm = \frac12(\sigma\_0 \pm \hat{\mathbf m}\cdot\boldsymbol\sigma) \\) (eigenprojectors), so

$$
\frac{\partial E\_{\mathrm{xc}}}{\partial n\_{\sigma'\sigma}(\mathbf r)} = v'\_\uparrow (P\_+)\_{\sigma\sigma'} + v'\_\downarrow (P\_-)\_{\sigma\sigma'} = V\_{\mathrm{xc},\sigma\sigma'}(\mathbf r)
$$

and the potential of TN2 Eq. 18–19, "computed in the local frame and transformed back to the global frame", is itself the exact derivative of \\( E\_{\mathrm{xc}} \\) (check 4).
Changes of the angles \\( (\theta,\varphi) \\) are contained in \\( P\_\pm \\), and no separate term appears in the force. The same holds for GGA as long as \\( v'\_\sigma \\) is the (discretized) functional derivative with respect to \\( n'\_\sigma \\).
However, at points with \\( \lVert\mathbf m\rVert \to 0 \\), \\( n'\_\pm \\) is not a smooth function of \\( \underline n \\) (\\( n'\_+ - n'\_- = \lVert\mathbf m\rVert \\)), and the GGA \\( \nabla n'\_\sigma \\) has a kink there, so
agreement with finite differences may deteriorate in their vicinity (an [Open] item).

[Convention] The existing implementation prepares \\( (V\_{\uparrow\uparrow}, V\_{\downarrow\downarrow}, \mathrm{Re}\\,V\_{\downarrow\uparrow}, \mathrm{Im}\\,V\_{\downarrow\uparrow}) \\) at each grid point and
contracts them with the changes of the density components due to movement of basis functions in the form \\( \sum\partial n\_{\uparrow\uparrow}V\_{\uparrow\uparrow} + \partial n\_{\downarrow\downarrow}V\_{\downarrow\downarrow} + 2\\,\mathrm{Re}(\partial n\_{\uparrow\downarrow}V\_{\downarrow\uparrow}) \\). This is the same as the equation above.

### 3.4 Partial-core-correction forces [Public] (derived) [Convention]

Since \\( n^{\mathrm{pcc}} \\) enters both spins of the local frame half each, the exact derivative is

$$
\frac{\partial E\_{\mathrm{xc}}}{\partial\mathbf R\_K}\Big|\_{\mathrm{pcc}} = \Delta V\sum\_p \tfrac12\bigl(v'\_\uparrow + v'\_\downarrow\bigr)\\,\frac{\partial n^{\mathrm{pcc}}\_K}{\partial\mathbf R\_K}
{}= \Delta V\sum\_p \tfrac12\mathrm{Tr}\\,\underline V\_{\mathrm{xc}}\\,\frac{\partial n^{\mathrm{pcc}}\_K}{\partial\mathbf R\_K}
$$

For open-core species (`<scf.pcc.opencore` is \\( \pm1 \\)), the exact form multiplies by the potential \\( v'\_\uparrow \\) (\\( v'\_\downarrow \\)) of the majority (minority) spin of the local frame.

[Convention] For ordinary species the existing implementation uses \\( \mathrm{Tr}\\,\underline V\_{\mathrm{xc}} = V\_{\mathrm{xc},\uparrow\uparrow} + V\_{\mathrm{xc},\downarrow\downarrow} \\) (same as the exact expression).
For open-core species, even in the noncollinear case, it multiplies twice the global-frame \\( V\_{\mathrm{xc},\uparrow\uparrow} \\) (\\( +1 \\)) or \\( V\_{\mathrm{xc},\downarrow\downarrow} \\) (\\( -1 \\)) by \\( \frac12\partial n^{\mathrm{pcc}}\_K \\).
This agrees with the exact derivative only at points where \\( \hat{\mathbf m} \\) points along \\( +z \\).

### 3.5 Forces of DFT+U, constraints and Zeeman terms [Convention]

These energies are functions of the occupations \\( \underline N\_I \\) (and their orbital components) built from \\( \rho \\) and \\( S \\).
The part due to changes of \\( \rho \\) combines into the Pulay term through the potential put into the Hamiltonian; the remainder is, from the derivative of \\( S \\) at fixed \\( \rho \\),
\\( \sum\_I\mathrm{Tr}[(\partial U/\partial\underline N\_I)\\,\partial\underline N\_I/\partial\mathbf R] \\).
The existing implementation evaluates this with the same expression as the DFT+U ('dual') force, using an effective potential matrix to which the constraint and Zeeman potentials are added.
The expression is in Sec. 11 of [DFT+U, constrained DFT and magnetic features](dftu-magnetism.md) (\\( \tfrac12(V\partial S + \partial S V) \\) contracted with \\( \rho \\) by the contraction of Sec. 1.4).

### 3.6 Force-check results of the existing implementation [Convention]

The result files of the distributed force check (`openmx -forcetest`, manual Sec. 72) record, for noncollinear inputs, the difference between the analytic forces and
central differences with \\( \pm 3\times10^{-4} \\) bohr. As a guideline, they are of the following order (all terms, maximum over components of atom 1).

| Input | Content | Analytic vs. numerical difference (Hartree/bohr) |
|---|---|---|
| `CH_NC.dat` | Noncollinear, no spin-orbit (molecule) | \\( 4\times10^{-4} \\) |
| `Cr2.dat` | Noncollinear, spin-orbit, DFT+U, constraint (Cluster) | \\( 1\times10^{-4} \\) |
| `GaAs_NC_SO.dat` | Noncollinear, spin-orbit (Band) | \\( 1\times10^{-3} \\) |
| `GaAs_NC_SO+U.dat` | Noncollinear, spin-orbit, DFT+U | \\( 3\times10^{-4} \\) |

These values include errors from grid discretization (egg-box effect) and SCF convergence. They are used as a reference for the upper bound of the reimplementation's tolerances.

## 4. Stress [Convention] [Design]

[Convention] The existing implementation cannot combine noncollinear calculations (`scf.SpinPolarization NC`) with stress calculation. When stress is requested
(`scf.stress.tensor on`, optimizations that move the cell, NPT molecular dynamics), it stops with an error at input checking.
This is independent of whether spin-orbit coupling is present. Other combinations that cannot be used with stress are orbital optimization and the DFT+U 'full' projector
("Output and gradient with respect to lattice vectors" in [Total energy and forces](energy-force.md)).

[Design] The first version of the reimplementation likewise disallows combining noncollinear and stress, and rejects the same inputs with an error. The outlook for future support is recorded here.

- The matrix-element terms (kinetic, nonlocal, projector-expansion \\( V^{\mathrm{na}} \\)) and the Pulay term only require replacing
  \\( g\_s\sum\_\sigma\rho^\sigma \\) in the collinear table (TN23 Eq. 16–31, 62–67) by the contraction of Sec. 1.4. The spin-orbit terms can be handled in the same way with
  \\( \partial/\partial\tau\_\gamma\\,\tau\_\eta \\) for each two-center factor.
- The exchange-correlation volume term and the LDA term have the same form as collinear when written in \\( n'\_\sigma \\). The additional GGA part can also be written with finite-difference gradients of \\( n'\_\sigma \\), but
  the treatment of points with \\( \lVert\mathbf m\rVert\to 0 \\) must be decided.
- The stress of the DFT+U, constraint and Zeeman terms has not been examined.

## 5. Initial spin configuration [Public] [Convention]

[Public] Manual Sec. 37: in noncollinear calculations, each line of `Atoms.SpeciesAndCoordinates` has 13 columns.

| Column | Content |
|---|---|
| 1–2 | Serial number, species |
| 3–5 | Coordinates |
| 6–7 | Initial numbers of up and down electrons \\( N^{\uparrow}\_{I,\mathrm{in}}, N^{\downarrow}\_{I,\mathrm{in}} \\) |
| 8–9 | Spin direction \\( (\theta, \varphi) \\) (degrees). Direction of the initial spin density, and target direction of the constraint and spin Zeeman term |
| 10–11 | Orbital-moment direction (degrees). Target direction of the orbital Zeeman term |
| 12 | Whether to apply the constraint and Zeeman terms (1/0) |
| 13 | Whether to promote DFT+U orbital polarization (`on`/`off`) |

[Convention]

- Angles are read in degrees and converted to radians. When the absolute value of \\( \theta \\) (both spin and orbital) is smaller than \\( 10^{-10} \\) degrees,
  a uniform random number of width \\( 10^{-5} \\) degrees (centered at 0) is added before conversion. To the spin \\( \theta, \varphi \\), a uniform random number of width \\( 10^{-10} \\) radians is further added after conversion.
  This slightly displaces exactly symmetric configurations; the effect on results is negligible.
- When \\( N^{\uparrow}\_{I,\mathrm{in}} < N^{\downarrow}\_{I,\mathrm{in}} \\), the spin direction is replaced by the opposite direction \\( -\hat{\mathbf n} \\), and
  \\( N^{\uparrow}\_{I,\mathrm{in}} \\) and \\( N^{\downarrow}\_{I,\mathrm{in}} \\) are swapped. The orbital-moment direction is reversed at the same time.
  The new angles are the polar coordinates of \\( -\hat{\mathbf n} \\) (\\( \theta\in[0,\pi] \\), \\( \varphi\in[-\pi/2, 3\pi/2] \\), Sec. 6.3).
  Example: an atom with columns 6 and 7 equal to 6.0, 8.0 and angles \\( (42^\circ, 79^\circ) \\) is treated as \\( N^\uparrow = 8 \\), \\( N^\downarrow = 6 \\), direction \\( (138^\circ, 259^\circ) \\).
  As a result, "up" is thereafter always the majority spin along the input direction. The target directions of the constraint and Zeeman terms are also the reversed directions.
- The initial density is built from the total density, given by the atomic density of atom \\( I \\) scaled by \\( (N^\uparrow\_{I,\mathrm{in}} + N^\downarrow\_{I,\mathrm{in}})/Z\_I \\), and
  \\( \mathbf m \\) with magnitude scaled by \\( (N^\uparrow\_{I,\mathrm{in}} - N^\downarrow\_{I,\mathrm{in}})/Z\_I \\) along the input direction
  (Sec. 6 of [Density matrix and charge mixing](density-mixing.md)).
- For the treatment of lines missing column 8 onward and for the orbital-polarization flag, see Sec. 8 of [DFT+U, constrained DFT and magnetic features](dftu-magnetism.md).

## 6. Output of local spin and orbital moments [Public] [Convention]

### 6.1 Per-atom spin moment [Convention]

The per-atom 2×2 Mulliken-type occupation \\( \underline M\_I \\) is built from \\( \underline A\_I \\) of Sec. 2.4 as follows.

$$
M^{\sigma\sigma}\_I = \mathrm{Re}\\,A^{\sigma\sigma}\_I,\qquad M^{\uparrow\downarrow}\_I = \overline{A^{\uparrow\downarrow}\_I},\qquad M^{\downarrow\uparrow}\_I = A^{\uparrow\downarrow}\_I
$$

(the rows of \\( \underline A\_I \\) are the orbitals of atom \\( I \\), the columns all orbitals of the partner; for the per-orbital decomposition the sum over \\( i\in I \\) is not taken). \\( \underline M\_I \\) is
converted to a vector by the same rule as Sec. 1.3, \\( \mathbf m = (2\\,\mathrm{Re}\\,M^{\uparrow\downarrow}, -2\\,\mathrm{Im}\\,M^{\uparrow\downarrow}, M^{\uparrow\uparrow}-M^{\downarrow\downarrow}) \\), and
\\( N^\uparrow\_I, N^\downarrow\_I = \frac12(\mathrm{Tr}\\,\underline M\_I \pm \lVert\mathbf m\_I\rVert) \\) and the direction \\( (\theta\_I, \varphi\_I) \\) of \\( \mathbf m\_I \\) are obtained.

- This \\( \underline M\_I \\) differs from \\( \underline N\_I \\) (the Hermitian part) of Sec. 2.4 in that its off-diagonal components use only the rows of atom \\( I \\).
  The difference between the two comes from the part of \\( \rho^{\uparrow\downarrow} \\) antisymmetric in the orbitals. By Sec. 1.2, the difference is 0 in systems without spin-orbit coupling.
  In systems with spin-orbit coupling, a difference can arise at low-symmetry atoms, and an apparent per-atom moment may appear even in time-reversal-symmetric systems
  (confirmed in a simple model: check 6; when atomic sites have high symmetry, as in GaAs, it is 0 by symmetry). The sum over all atoms agrees between the two.
- The output moment is \\( N^\uparrow\_I - N^\downarrow\_I \\) (\\( \mu\_B \\), assuming \\( g = 2 \\)) and is always nonnegative. Its direction is expressed by angles.
- The spin moment of the whole system is the magnitude and direction of \\( \sum\_I (N^\uparrow\_I - N^\downarrow\_I)\\,\hat{\mathbf m}\_I \\).
- The obtained \\( N^\uparrow\_I, N^\downarrow\_I, \theta\_I, \varphi\_I \\) are also used, e.g., as initial values for the next step of MD or geometry optimization
  (the target directions of constraints and Zeeman terms remain the input values).

### 6.2 Orbital moment [Public] [Convention]

[Public] Manual Sec. 39, TN2 Eq. 151, 160–162: in noncollinear calculations, the per-atom orbital moment is obtained as the expectation value of the orbital angular momentum in the localized basis.
In the real basis the angular-momentum matrix elements are purely imaginary, so it is 0 without spin-orbit coupling.

[Convention] The existing implementation builds \\( \mathbf M\_I \\) from \\( \langle\phi\_{Im}\rvert\hat{\mathbf L}\lvert\phi\_{In}\rangle \\) on atom \\( I \\) (a one-center matrix) and
the 'dual'-type product (average of \\( \rho S \\) and \\( S\rho \\)) of \\( \mathrm{Im}\\,\rho^{\uparrow\uparrow} \\), \\( \mathrm{Im}\\,\rho^{\downarrow\downarrow} \\), and outputs its magnitude and angles.
The expression and signs are in Sec. 7 of [DFT+U, constrained DFT and magnetic features](dftu-magnetism.md). The angle convention is the same as Sec. 6.3 below.

### 6.3 Output format [Public] [Convention]

[Public] Manual Secs. 37 and 39: the Mulliken section of `System.Name.out` contains the total spin moment and its angles (degrees), a per-atom table of
`Up`, `Down`, `Sum`, `Diff`, `theta`, `phi`, and a sum line. The subsequent per-orbital decomposition also carries angles.
The orbital-moment section contains the total orbital moment and its angles, per-atom values and angles, and the per-orbital decomposition.

[Convention]

- `Up`, `Down` are \\( N^\uparrow\_I, N^\downarrow\_I \\); `Sum` is their sum, `Diff` their difference; `theta`, `phi` are the angles (degrees) of \\( \hat{\mathbf m}\_I \\).
  "Total spin moment (muB)" on the first line is the total spin moment of Sec. 6.1.
- Angles are computed with \\( \theta\in[0,\pi] \\), and \\( \varphi \\) as \\( \arcsin(y/\sqrt{x^2+y^2}) \\) for \\( x \ge 0 \\) and \\( \pi - \arcsin(\cdot) \\) for \\( x < 0 \\),
  i.e. in the range \\( [-\pi/2, 3\pi/2] \\). If the magnitude is below \\( 10^{-14} \\), \\( \theta = \pi/2 \\), \\( \varphi = 0 \\); if the \\( xy \\) component is below \\( 10^{-14} \\), \\( \varphi = 0 \\).
  The value converted to degrees is written as the (signed) remainder modulo 360, so \\( \varphi \\) takes values in the range \\( -90 \\) to \\( 270 \\).
- The standard output of each SCF iteration shows, per atom (the first 20 atoms with the default output level), `MulP` up, down, sum, difference and \\( (\theta,\varphi) \\),
  the orbital-moment magnitude `Ml` and angles, and the magnitude `Ml+s` and angles of the combined spin and orbital vector \\( (N^\uparrow\_I - N^\downarrow\_I)\hat{\mathbf m}\_I + \mathbf M\_I \\).

## 7. Design [Design]

- **Type of the density matrix.** The reimplementation holds \\( \rho^{\sigma\sigma'}(\mathbf L) \\) and \\( \varepsilon^{\sigma\sigma'}(\mathbf L) \\) as complex 2×2 spin blocks per atom pair
  (the noncollinear form of the [sparse block format](../design/sparse-block-format.md)), also storing \\( \downarrow\uparrow \\).
  The 4–6 real arrays of the existing implementation are converted only at fixture input/output.
- **A single contraction function.** All of energy, forces and Pulay are written with \\( \langle\rho, X\rangle \\) of Sec. 1.4 (\\( X \\) is a 2×2 block or \\( O\otimes\sigma\_0 \\)),
  and the collinear case passes through the same function as the special case \\( \rho = \rho^\sigma\oplus\rho^{\sigma'} \\) with only diagonal spin.
  For spin-independent \\( X \\), the shortcut \\( \mathrm{Re}(\rho^{\uparrow\uparrow}+\rho^{\downarrow\downarrow}) \\) may be used, but the result must agree with the general form.
- **Orientation of the grid density.** The grid holds \\( \underline n \\) (convention of Sec. 1.3, \\( n\_{\uparrow\downarrow} = (m\_x - i m\_y)/2 \\)), and
  the indices are swapped when assembling from \\( \rho \\). Transformation to the local frame is done only inside the energy and XC evaluation; the storage format on the grid is unchanged.
- **Splitting of `Uxc0` and `Uxc1`.** For comparison with fixtures, they are split, as in the existing implementation, into majority and minority spins of the local frame.
- **Open-core forces.** The default is the exact derivative (\\( v'\_\uparrow \\) or \\( v'\_\downarrow \\) of Sec. 3.4), and a compatibility mode for comparison with the existing implementation
  (global-frame diagonal components) is provided.
- **Mulliken spin moment.** The default is the same definition as the existing implementation (Sec. 6.1), with the definition by the Hermitian part available as an option.
- **Stress.** Combining noncollinear and stress is an error in the first version (Sec. 4).

## 8. Tests [Design]

### 8.1 Algebraic checks (`tools/check_noncollinear_energy_force.py`)

| No. | Content | Expected result |
|---|---|---|
| 1 | Spinor band energy \\( = \langle\rho, H\rangle \\), real/imaginary decomposed form, Hermiticity of \\( \rho \\) | \\( 10^{-10} \\) |
| 2 | Pulay identity (spinor, \\( S^{(2)} = \sigma\_0\otimes S \\)), and sufficiency of \\( \mathrm{Re}(\varepsilon^{\uparrow\uparrow}+\varepsilon^{\downarrow\downarrow}) \\) alone | \\( 10^{-7} \\) (finite difference), \\( 10^{-12} \\) (difference between the two forms) |
| 3 | \\( \Delta V\sum\mathrm{Tr}[\underline V\\,\underline n] = \langle\rho, H^{\mathrm{grid}}\rangle \\), Hermiticity of \\( \underline n \\) (index swap of Sec. 1.3) | \\( 10^{-10} \\) |
| 4 | \\( \underline V\_{\mathrm{xc}} \\) computed in the local frame is the derivative of \\( E\_{\mathrm{xc}} \\) (model functional, finite differences with Hermitian perturbations) | \\( 10^{-7} \\) |
| 5 | Global spin rotation without spin-orbit: eigenvalues and the \\( \mathrm{Re}(\rho^{\uparrow\uparrow}+\rho^{\downarrow\downarrow}) \\) contraction are invariant, atomic moments rotate as vectors | \\( 10^{-10} \\) |
| 6 | Time-reversal symmetric (with spin-orbit): \\( \rho^{\uparrow\downarrow} \\) is antisymmetric in the orbitals, the Hermitian-part atomic moment is 0 | \\( 10^{-8} \\) |
| 7 | Agreement of the two ways of writing the constraint target \\( \underline N^0 \\), and recovery of the direction | Rounding error |
| 8 | With \\( H^{\uparrow\downarrow} = 0 \\) and real \\( H^{\sigma\sigma} \\), the spinor result agrees with two collinear ones | \\( 10^{-10} \\) |

### 8.2 Invariants of the full calculation

| Type | Content | Tolerance guideline |
|---|---|---|
| Reduction to collinear | `NC` without spin-orbit and all atoms at \\( \theta = 0 \\) vs. `On` with the same input: `Utot` and all components (including `Uxc0`, `Uxc1`, for systems with \\( n\_{\uparrow\uparrow}\ge n\_{\downarrow\downarrow} \\) at all points), forces, moments | \\( 10^{-8} \\) Hartree, \\( 10^{-6} \\) Hartree/bohr (within SCF convergence) |
| Reduction to nonmagnetic | `NC` without spin-orbit and \\( N^\uparrow = N^\downarrow \\) vs. `Off` | Same as above |
| Global spin rotation | Without spin-orbit, rotate the directions of all atoms by the same rotation: `Utot`, each component and forces invariant, moment directions rotate | Same as above; angles \\( 10^{-4} \\) degrees |
| Time reversal | With spin-orbit, reverse the directions of all atoms (\\( \theta\to\pi-\theta \\), \\( \varphi\to\varphi+\pi \\)): `Utot` and forces invariant, spin and orbital moments reversed | Same as above |
| Non-invariance under rotation | With spin-orbit, rotating the directions changes the energy (magnetic anisotropy; sign and magnitude checked against fixtures) | — |
| Identity | \\( U\_{\mathrm{ele}} = \langle\rho, H\rangle \\) (Sec. 2.5) | \\( 10^{-10} \\) |
| Sum of forces | The matrix-element terms (including spin-orbit) and the Pulay term sum to 0 | \\( 10^{-10} \\) |
| Constraint | Setting `scf.Constraint.NC.Spin.v` to 0 agrees with no constraint; `Ucs` is 0 for a solution aligned with the target directions | \\( 10^{-10} \\) |

### 8.3 Finite differences

- Central difference of `Utot` moving one atom by \\( \pm 3\times10^{-4} \\) bohr vs. analytic forces (all terms).
  Targets: small molecules with spin-orbit (e.g. Pt₂, Cr₂ corresponding to `Cr2.dat`) and a noncollinear crystal (corresponding to `GaAs_NC_SO.dat`).
  The tolerance is bounded above by the level of the existing implementation (Sec. 3.6), starting from \\( 10^{-3} \\) Hartree/bohr.
- Per-term checks: matrix-element terms at fixed \\( \rho \\), \\( \varepsilon \\) (`Ukin`, `Unl` including the spin-orbit part) and the Pulay term to \\( 10^{-7} \\).
- Grid terms are examined both for configurations without points where \\( \lVert\mathbf m\rVert \\) becomes 0 (e.g. ferromagnetically aligned configurations) and with such points,
  and the effect of the GGA kink is recorded.

### 8.4 Comparison with fixtures

The following are written out from the existing implementation and compared.

- Energies: `Uele`, `Ukin`, `UH0`, `UH1`, `Una`, `Unl`, `Uxc0`, `Uxc1`, `Ucore`, `Uhub`, `Ucs`, `Uzs`, `Uzo`, `Utot`.
- Intermediate quantities: converged \\( \rho^{\sigma\sigma'}(\mathbf L) \\) (6 components), \\( \mathrm{Re}\\,\varepsilon^{\sigma\sigma} \\), \\( V^{\mathrm{NL},\sigma\sigma'}(\mathbf L) \\) including spin-orbit,
  and \\( \underline n \\) (global frame) and \\( \underline V\_{\mathrm{xc}} \\) on the grid. From these, energy components and forces can be recomputed without solving the eigenvalue problem.
- Output: \\( N^\uparrow, N^\downarrow, \theta, \varphi \\) of the Mulliken section, total spin moment, orbital moments, `<coordinates.forces`.
- Targets: the distributed input examples `Cr2.dat` (Cluster, spin-orbit, DFT+U, constraint) and `GaAs.dat` (Band, spin-orbit),
  and the force-check inputs `CH_NC.dat`, `GaAs_NC_SO.dat`, `GaAs_NC_SO+U.dat`.

| Quantity | Tolerance guideline (same grid, same basis) |
|---|---|
| Each energy component, `Utot` | \\( 10^{-6} \\) Hartree |
| \\( N^\uparrow, N^\downarrow \\), moments | \\( 10^{-5} \\) |
| Angles | \\( 10^{-3} \\) degrees (only atoms with moments of \\( 0.1\\,\mu\_B \\) or more) |
| Forces | \\( 10^{-5} \\) Hartree/bohr |

## [Open]

- Whether the force of `Uzo` includes the dependence of the orbital moment through \\( S \\) (the dual-type product) (to be reconciled with the treatment in Sec. 11 of [DFT+U, constrained DFT and magnetic features](dftu-magnetism.md)).
- Whether to make the existing implementation's behavior of using global-frame components the reimplementation default for partial-core-correction forces of open-core species (Sec. 3.4).
- How large an apparent moment the row-only definition of the Mulliken spin moment (Sec. 6.1) actually produces in low-symmetry systems with spin-orbit coupling (verification on real systems).
- Whether the GGA \\( \underline V\_{\mathrm{xc}} \\) is a discretely exact derivative (shared with [Total energy and forces](energy-force.md)), and
  the treatment of the kink of \\( n'\_\sigma \\) at points with \\( \lVert\mathbf m\rVert\to 0 \\).
- How to allocate the GGA per-spin energy density in splitting `Uxc0` and `Uxc1` (shared with the collinear chapter).
- Reproducibility of the random numbers added to angles (confirming, when generating fixtures, that their effect is sufficiently smaller than the tolerances).
- Noncollinear energies (including the Harris type) and forces with solvers such as `DC`, `Krylov`, `NEGF`, `Cluster-LNO`.
- Noncollinear stress (absent in the existing implementation; whether to implement the outlook of Sec. 4).
- (Resolved) The contraction in Sec. 6 of [Nonlocal pseudopotential](nonlocal-pseudopotential.md) was unified with \\( \sum\rho^{\sigma\sigma'}\_{ab}(\mathbf L)V^{\sigma\sigma'}\_{ab}(\mathbf L) \\) (without complex conjugation) of Sec. 1.4 of this chapter.
