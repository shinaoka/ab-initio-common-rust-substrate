# Total energy and forces

> Status: draft (awaiting review). The equations re-derive, in the notation of this note, the methods of technical note 1 (Ver. 1.3, hereafter TN1) and
> the technical note on stress (Shiihara–Ozaki, hereafter TN23).
> The parts checked numerically are in `tools/check_energy_force.py`.

## Symbols [Design]

In addition to the symbols of [Notation](../../../conventions/notation.md) and [sparse block format](../design/sparse-block-format.md), the following are used.

| Symbol | Meaning |
|---|---|
| \\( \sigma \\) | Spin (\\( \uparrow, \downarrow \\) in the collinear case) |
| \\( \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L) \\) | Real-space density matrix (per spin) |
| \\( \varepsilon^{\sigma}\_{I\alpha,J\beta}(\mathbf L) \\) | Energy-weighted density matrix (defined below) |
| \\( O\_{I\alpha,J\beta}(\mathbf L) \\) | Matrix element for atom pair \\( (I, J, \mathbf L) \\): \\( S \\) (overlap), \\( T \\) (kinetic), \\( V^{\mathrm{na}} \\), \\( V^{\mathrm{NL}} \\), etc. |
| \\( \boldsymbol\tau\_{IJ\mathbf L} = \mathbf R\_J + \mathbf L - \mathbf R\_I \\) | Interatomic vector |
| \\( Z\_I \\) | Number of valence electrons of atom \\( I \\) (core charge of the pseudopotential) |
| \\( n\_\sigma(\mathbf r) \\), \\( n = \sum\_\sigma n\_\sigma \\) | Valence density |
| \\( n^{a}\_I(r) \\), \\( n^{a} = \sum\_{I,\mathbf L} n^{a}\_I \\) | Valence density of the isolated atom (density in the PAO file) and its superposition |
| \\( \delta n = n - n^{a} \\) | Difference density |
| \\( V^{a}\_{H,I}(r) \\) | Hartree potential produced by \\( n^{a}\_I \\); \\( Z\_I/r \\) far away |
| \\( \delta V\_H \\) | Hartree potential produced by \\( \delta n \\) |
| \\( V^{\mathrm{na}}\_I = V^{\mathrm{loc}}\_I + V^{a}\_{H,I} \\) | Neutral-atom potential ([neutral atom](neutral-atom.md)) |
| \\( n^{\mathrm{pcc}} \\) | Partial-core-correction density (sum over atoms) |
| \\( \mathbf r\_p \\), \\( \Delta V \\) | Real-space grid points and the volume per grid point |
| \\( \Omega \\) | Unit-cell volume |
| \\( r\_{c,I} \\) | PAO cutoff radius of atom \\( I \\) |

The sum over atom pairs \\( \sum\_{I,J,\mathbf L} \\) runs over \\( I \\) in the unit cell and
all \\( (J, \mathbf L) \\) satisfying \\( \lVert\boldsymbol\tau\_{IJ\mathbf L}\rVert \le r\_{c,I} + r\_{c,J} \\) (including \\( J = I \\), \\( \mathbf L = 0 \\)).
Below, the sum \\( \sum\_{I,J,\mathbf L}\sum\_{\alpha\beta} \\) is sometimes abbreviated as \\( \sum\_{\mathrm{pairs}} \\).

[Convention] In O(N) solvers (`DC`, `Krylov`, `DC-LNO`) the density matrix is not symmetric, \\( \rho\_{I\alpha,J\beta}(\mathbf L) \ne \rho\_{J\beta,I\alpha}(-\mathbf L) \\), so
the pair sums for energy and forces are taken over **ordered pairs**, counting \\( (I,J,\mathbf L) \\) and \\( (J,I,-\mathbf L) \\) as separate pairs
(one must not reconstruct one of them by Hermiticity and double it). See "Asymmetry of the density matrix" in [O(N) methods](on-methods.md) for details.

### Density matrix and energy-weighted density matrix [Public]

From eigenvectors \\( c^{\sigma}\_{n}(\mathbf k) \\), eigenvalues \\( e^{\sigma}\_{n}(\mathbf k) \\), occupations \\( f^{\sigma}\_{n}(\mathbf k) \\), and
\\( \mathbf k \\)-point weights \\( w\_{\mathbf k} \\) (summing to 1), define

$$
\rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k} \sum\_n f^{\sigma}\_n(\mathbf k)\\, \mathrm{Re}\left[ e^{i\theta(\mathbf k,\mathbf L)}\\, c^{\sigma\ast}\_{I\alpha,n}(\mathbf k)\\, c^{\sigma}\_{J\beta,n}(\mathbf k) \right]
$$

$$
\varepsilon^{\sigma}\_{I\alpha,J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k} \sum\_n f^{\sigma}\_n(\mathbf k)\\, e^{\sigma}\_n(\mathbf k)\\, \mathrm{Re}\left[ e^{i\theta(\mathbf k,\mathbf L)}\\, c^{\sigma\ast}\_{I\alpha,n}(\mathbf k)\\, c^{\sigma}\_{J\beta,n}(\mathbf k) \right]
$$

(TN1 Eq. 4, 74. The phase \\( \theta \\) follows the Bloch-sum convention of the [sparse block format](../design/sparse-block-format.md).)
In spin-unpolarized calculations, \\( \rho^{\uparrow} = \rho^{\downarrow} \\) is stored for one spin only and doubled in the energy [Convention].
The equations below assume the collinear case (including spin-unpolarized), a real basis and real \\( \rho \\). The noncollinear case (spin 2×2 density matrix, spin-orbit coupling) is treated in [Energy, forces and stress for noncollinear and spin-orbit calculations](noncollinear-energy-force.md).

## Decomposition of the total energy [Public]

The naive decomposition of TN1 Eq. 7 is

$$
E = E\_{\mathrm{kin}} + E^{\mathrm{loc}}\_{\mathrm{ec}} + E^{\mathrm{NL}}\_{\mathrm{ec}} + E\_{\mathrm{ee}} + E\_{\mathrm{xc}} + E\_{\mathrm{cc}}
$$

Here \\( E^{\mathrm{loc}}\_{\mathrm{ec}} = \int n \sum\_I V^{\mathrm{loc}}\_I \\),
\\( E\_{\mathrm{ee}} = \frac12\iint n(\mathbf r)n(\mathbf r')/\lvert\mathbf r-\mathbf r'\rvert \\), and
\\( E\_{\mathrm{cc}} = \frac12\sum\_{I\ne J} Z\_I Z\_J/\lvert\mathbf R\_I-\mathbf R\_J\rvert \\) (including periodic images).
To make the cancellation among long-range terms explicit, substitute \\( n = n^{a} + \delta n \\) and rearrange.
From

$$
E\_{\mathrm{ee}} = \tfrac12\int n^{a} V^{a}\_H + \int \delta n\\, V^{a}\_H + \tfrac12\int \delta n\\, \delta V\_H,\qquad
\int n\\, V^{\mathrm{na}} = E^{\mathrm{loc}}\_{\mathrm{ec}} + \int n^{a} V^{a}\_H + \int \delta n\\, V^{a}\_H
$$

(\\( V^{\mathrm{na}} = \sum\_I V^{\mathrm{na}}\_I \\), \\( V^{a}\_H = \sum\_I V^{a}\_{H,I} \\)), one obtains

$$
E^{\mathrm{loc}}\_{\mathrm{ec}} + E\_{\mathrm{ee}} + E\_{\mathrm{cc}}
{}= \underbrace{\int n\\, V^{\mathrm{na}}}\_{E\_{\mathrm{na}}}
{}+ \underbrace{\tfrac12\int \delta n\\, \delta V\_H}\_{E\_{\delta\mathrm{ee}}}
{}+ \underbrace{E\_{\mathrm{cc}} - \tfrac12\int n^{a} V^{a}\_H}\_{E\_{\mathrm{scc}}}
$$

(TN1 Eq. 19–22). Therefore

$$
E = E\_{\mathrm{kin}} + E\_{\mathrm{na}} + E^{\mathrm{NL}}\_{\mathrm{ec}} + E\_{\delta\mathrm{ee}} + E\_{\mathrm{xc}} + E\_{\mathrm{scc}}
$$

(TN1 Eq. 23, TN23 Eq. 1).

- Each pair of \\( E\_{\mathrm{scc}} = \frac12\sum\_{I,J}\left[ \frac{Z\_I Z\_J}{\lvert\cdot\rvert}(1-\delta\_{IJ}^{\mathbf L=0}) - \int n^{a}\_I V^{a}\_{H,J} \right] \\)
  is exactly 0 at distances where the supports of \\( n^{a}\_I \\) and \\( n^{a}\_J \\) do not overlap
  (because the potential outside a spherically symmetric charge \\( Z\_I \\) equals that of a point charge). That is, it is a screened short-range pair potential.
  The \\( I = J \\), \\( \mathbf L = 0 \\) term of the sum is the atomic self-energy \\( -\frac12\int n^{a}\_I V^{a}\_{H,I} \\), a constant independent of the configuration.
- \\( V^{\mathrm{na}}\_I \\) in \\( E\_{\mathrm{na}} \\) is 0 outside the neutral atom, so the matrix elements involve only a finite number of atom pairs.
- The long-range electrostatic interaction is carried only by \\( E\_{\delta\mathrm{ee}} \\), which is solved by FFT ([grid](grid-poisson-xc.md)).
- Exchange-correlation, including the partial core correction, is
  \\( E\_{\mathrm{xc}} = \int (n + n^{\mathrm{pcc}})\\, \epsilon\_{\mathrm{xc}}(n\_\uparrow + \tfrac12 n^{\mathrm{pcc}}, n\_\downarrow + \tfrac12 n^{\mathrm{pcc}}) \\) (TN1 Eq. 15; in GGA the gradients also enter as arguments).

## Output energy components [Convention]

In the "Total energy (Hartree)" section of the output file (`System.Name.out`), values are printed with the following labels.
A period follows each label (`Utot.` etc.). The unit is Hartree.

| Label | Definition (notation of this note) | Evaluation method |
|---|---|---|
| `Uele` | \\( \sum\_\sigma\sum\_{\mathbf k} w\_{\mathbf k}\sum\_n f^{\sigma}\_n e^{\sigma}\_n \\) (band energy; not included in `Utot`) | Sum of eigenvalues |
| `Ukin` | \\( \sum\_\sigma \sum\_{\mathrm{pairs}} \rho^{\sigma} T \\) | Matrix trace (two-center integrals) |
| `UH0` | \\( -\frac12\sum\_{I,J,\mathbf L} s\_{IJ} \int n^{a}\_I V^{a}\_{H,J} \\) (including the self term \\( J=I, \mathbf L=0 \\); \\( s\_{IJ} \\) below) | Numerical integration per atom pair |
| `UH1` | \\( E\_{\delta\mathrm{ee}} = \frac12\Delta V\sum\_p \delta n(\mathbf r\_p)\\, \delta V\_H(\mathbf r\_p) \\) | Sum on the grid |
| `Una` | \\( E\_{\mathrm{na}} \\) | Matrix trace (projector expansion) or sum on the grid |
| `Unl` | \\( \sum\_\sigma \sum\_{\mathrm{pairs}} \rho^{\sigma} V^{\mathrm{NL}} \\) | Matrix trace (products of two-center integrals) |
| `Uxc0`, `Uxc1` | Per-spin parts of \\( E\_{\mathrm{xc}} \\) (below) | Sum on the grid + atomic-density correction |
| `Ucore` | \\( \frac12\sum\_{I,(J,\mathbf L)\ne(I,0)} Z\_I Z\_J / \lVert\boldsymbol\tau\_{IJ\mathbf L}\rVert \\) | Sum over atom pairs |
| `Uhub` | DFT+U energy | [DFT+U, constrained DFT and magnetic features](dftu-magnetism.md) (noncollinear: [noncollinear](noncollinear-energy-force.md)) |
| `Ucs` | Spin-orientation constraint energy | [DFT+U, constrained DFT and magnetic features](dftu-magnetism.md) Sec. 5, [noncollinear](noncollinear-energy-force.md) Sec. 2.4 |
| `Uzs`, `Uzo` | Zeeman terms of spin and orbital magnetic moments | [DFT+U, constrained DFT and magnetic features](dftu-magnetism.md) Sec. 6, [noncollinear](noncollinear-energy-force.md) Sec. 2.4 |
| `Uef` | External electric-field energy \\( \Delta V\sum\_p n(\mathbf r\_p) V\_{\mathrm{ef}}(\mathbf r\_p) \\) (electrons only; same definition as Sec. 3.3 of [ESM, external field, core levels](esm-field-corehole.md)) | Sum on the grid |
| `UvdW` | DFT-D2 or D3 dispersion energy (two-body terms only; the D3 three-body term is not implemented; [Exact exchange, dispersion correction](exx-vdw.md)) | Sum over atom pairs |
| `Uch` | Penalty term for creating a core hole | [ESM, external field, core levels](esm-field-corehole.md) Sec. 5.2 |
| `Utot` | Sum below | |
| `UpV` | Product of applied pressure and volume (`MD.applied.pressure`) | |
| `Enpy` | `Utot` + `UpV` (enthalpy) | |

The total energy is

$$
U\_{\mathrm{tot}} = U\_{\mathrm{kin}} + U\_{\mathrm{H0}} + U\_{\mathrm{H1}} + U\_{\mathrm{na}} + U\_{\mathrm{nl}} + U\_{\mathrm{xc0}} + U\_{\mathrm{xc1}} + U\_{\mathrm{core}} + U\_{\mathrm{hub}} + U\_{\mathrm{cs}} + U\_{\mathrm{zs}} + U\_{\mathrm{zo}} + U\_{\mathrm{ef}} + U\_{\mathrm{vdW}} + U\_{\mathrm{ch}}
$$

The correspondence with the symbols of TN1 is \\( E\_{\mathrm{scc}} \leftrightarrow U\_{\mathrm{core}} + U\_{\mathrm{H0}} \\),
\\( E\_{\delta\mathrm{ee}} \leftrightarrow U\_{\mathrm{H1}} \\), \\( E\_{\mathrm{xc}} \leftrightarrow U\_{\mathrm{xc0}} + U\_{\mathrm{xc1}} \\).

- Even with an electronic temperature, `Utot` does not include the electronic entropy term \\( -TS \\).
- Without spin polarization, `Uxc0` and `Uxc1` are equal, each being half of \\( E\_{\mathrm{xc}} \\).
  For collinear spin, \\( U\_{\mathrm{xc}\sigma} = \Delta V\sum\_p (n\_\sigma + \tfrac12 n^{\mathrm{pcc}})\\, \epsilon\_{\mathrm{xc}}(\mathbf r\_p) \\)
  (\\( \epsilon\_{\mathrm{xc}} \\) is the energy per electron determined by the total density) plus half of the correction below each.
- When `scf.EigenvalueSolver` is `Cluster-LNO`, `Utot` is constructed by a different formula (the Harris-type form described later),
  and the output labels are a different sequence including `UH_dc`, `Uxc0_dc`, `Uxc1_dc`.

### Terms evaluated as matrix traces [Convention]

`Ukin`, `Una` (with projector expansion), and `Unl` are sums of products of the density matrix and matrix elements over the stored range of atom pairs.

$$
U\_X = g\_s \sum\_{\sigma} \sum\_{I,J,\mathbf L}\sum\_{\alpha\beta} \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\, X\_{I\alpha,J\beta}(\mathbf L)
$$

\\( g\_s \\) is 2 without spin polarization and 1 otherwise.
\\( T \\) is constructed in [Two-center integrals](two-center-overlap-kinetic.md), \\( V^{\mathrm{NL}} \\) in [Nonlocal pseudopotential](nonlocal-pseudopotential.md),
and \\( V^{\mathrm{na}} \\) in the [neutral atom](neutral-atom.md) chapter.
When `scf.ProExpn.VNA` (default on) is off, `Una` is evaluated as the grid sum \\( \Delta V\sum\_p n(\mathbf r\_p) V^{\mathrm{na}}(\mathbf r\_p) \\).

### Terms evaluated on the grid [Convention]

From grid values \\( n\_\sigma(\mathbf r\_p) \\), \\( n^{a}(\mathbf r\_p) \\), \\( n^{\mathrm{pcc}}(\mathbf r\_p) \\), \\( \delta V\_H(\mathbf r\_p) \\),
\\( \epsilon\_{\mathrm{xc}}(\mathbf r\_p) \\), integrals are evaluated by the trapezoidal rule (simple sum over grid points multiplied by \\( \Delta V \\)).

- `UH1`: \\( \frac12\Delta V\sum\_p (n - n^{a})\\, \delta V\_H \\).
- `Uxc`: \\( \Delta V\sum\_p \sum\_\sigma (n\_\sigma + \frac12 n^{\mathrm{pcc}})\\, \epsilon\_{\mathrm{xc}} \\).
  \\( \epsilon\_{\mathrm{xc}} \\) on the grid points is obtained with the same functional as the Hamiltonian for `Uele`, on the same grid, with the same finite-difference gradients (for GGA).

### Atomic-density correction to exchange-correlation [Convention]

When `scf.Exc0.correction.flag` (default on) is on, the following \\( \Delta E^{(0)} \\) is added to `Uxc0 + Uxc1` (half each).

$$
\Delta E^{(0)} = -\Delta V\sum\_p \left(n^{a} + n^{\mathrm{pcc}}\right)\epsilon^{\mathrm{CA}}\\!\left(n^{a} + n^{\mathrm{pcc}}\right)\Big\rvert\_{\mathbf r\_p}
{}+ \sum\_{I}\int\_{\lVert\mathbf x\rVert < r\_{c,I}} d^3x\\; n^{a}\_I(\mathbf x)\\, \epsilon^{\mathrm{CA}}\\!\left(n^{a}(\mathbf R\_I + \mathbf x)\right)
$$

- \\( \epsilon^{\mathrm{CA}} \\) is always the spin-unpolarized Ceperley–Alder LDA, regardless of the functional chosen for the calculation
  (details of the parametrization are to be fixed in the [grid](grid-poisson-xc.md) chapter).
- The first term is a sum on the coarse real-space grid, and the second is a fine atom-centered quadrature. Within a sphere centered on atom \\( I \\),
  the radial direction uses 150-point Gauss–Legendre on \\( [0, r\_{c,I}] \\) and the angular direction uses a Lebedev rule
  (number of points `scf.Num.Leb.Grid`, default 590).
  \\( n^{a}(\mathbf R\_I + \mathbf x) \\) is the sum of \\( n^{a}\_J \\) of neighboring atoms (atoms within the atom-pair range), each \\( n^{a}\_J \\) being 0 outside its own cutoff.
- The intent is to re-evaluate the exchange-correlation energy of the superposition of atomic densities with an accurate atom-centered quadrature rather than on the grid
  (reducing errors due to grid coarseness, in particular errors that vary with the position of atoms relative to the grid).
- The density of the second term does not include \\( n^{\mathrm{pcc}} \\), while that of the first term does. For atoms with partial core correction,
  the two terms are not two evaluations of the same quantity, and a difference \\( \int[(n^{a}+n^{\mathrm{pcc}})\epsilon(n^{a}+n^{\mathrm{pcc}}) - n^{a}\epsilon(n^{a})] \\) remains.
  This is recorded as behavior of the existing implementation. Whether the reimplementation should also match this definition is [Open].

### Screened interatomic terms: Ucore and UH0 [Public] [Convention]

[Public] TN1 states that since the second term of \\( E\_{\mathrm{scc}} \\) is a function of the interatomic distance only, it is tabulated numerically and interpolated (after TN1 Eq. 22).

[Convention] The existing implementation does not use a distance table; instead it evaluates the following numerical integral for each atom pair every time.

$$
J\_{IJ}(d) = \int d^3 r\\; n^{a}\_I(\mathbf r)\\, V^{a}\_{H,J}(\mathbf r - d\\,\hat{\mathbf z})
$$

- The integral is done in cylindrical coordinates with the pair axis as \\( z \\). The support \\( [-r\_{c,I}, r\_{c,I}] \\) of \\( n^{a}\_I \\) is divided equally in the \\( z \\) direction
  (target spacing \\( \pi/\sqrt{9600} \approx 0.032 \\) bohr; the spacing is adjusted so that there is an odd number of points including both ends),
  the cylindrical radial direction \\( [0, r\_{c,I}] \\) uses 150-point Gauss–Legendre, and the azimuth is analytic, \\( 2\pi \\).
  \\( n^{a}\_I \\) is the atomic density from the PAO file, and \\( V^{a}\_{H,J} \\) is the atomic Hartree potential from the pseudopotential file
  (\\( Z\_J/r \\) outside the range), both interpolated in the radial direction.
- \\( U\_{\mathrm{H0}} \\) takes, for each pair, the average of \\( J\_{IJ} \\) and \\( J\_{JI} \\), multiplies by \\( -\frac12 \\), and sums (symmetrization).
- Pairs of distinct atoms (\\( J \ne I \\) or \\( \mathbf L \ne 0 \\)) are multiplied by a constant per species pair

$$
s\_{IJ} = \frac{Z\_I Z\_J / (r\_{c,I} + r\_{c,J})}{J\_{IJ}(r\_{c,I} + r\_{c,J})}
$$

  (1 for the self term). As a result, at the edge of the atom-pair range \\( d = r\_{c,I} + r\_{c,J} \\),
  \\( Z\_I Z\_J/d - s\_{IJ} J\_{IJ}(d) = 0 \\) holds numerically as well, and the energy is continuous when pairs enter or leave the range.
  \\( s\_{IJ} \\) is computed once for the initial structure and kept fixed, and is not differentiated in forces.
- \\( U\_{\mathrm{core}} \\) is the sum of point charges over the same range of atom pairs. Pairs outside the range enter neither \\( U\_{\mathrm{core}} \\) nor \\( U\_{\mathrm{H0}} \\), but
  since the screened pair sum is 0 outside the range, the total is unaffected by the truncation.
  On the other hand, \\( U\_{\mathrm{core}} \\) and \\( U\_{\mathrm{H0}} \\) taken separately are both large numbers that depend on the truncation. Comparisons are made on their sum.

[Design] In the reimplementation, \\( J\_{IJ}(d) \\) may be precomputed per species pair as a table in \\( d \\) and interpolated by splines
(the TN1 method). In that case too, the endpoint normalization \\( s\_{IJ} \\) and derivatives consistent with the table interpolation are used.
The cylindrical-grid quadrature has been confirmed, for a Gaussian density and its Hartree potential \\( Z\\,\mathrm{erf}(\sqrt{a}r)/r \\),
to agree with the analytic solution \\( Z\_I Z\_J\\,\mathrm{erf}(\sqrt{\mu}d)/d \\) (\\( \mu = a\_I a\_J/(a\_I+a\_J) \\)) to a relative \\( 10^{-10} \\).

### Other terms

- [Public] `UvdW`: Grimme's DFT-D2 and D3 (manual Sec. 54). [Convention] Only two-body terms; the D3 three-body term (Axilrod–Teller–Muto type) is not implemented. Since it is a function of atomic positions only, forces and stress are added by analytic derivatives (equations in "Forces and stress" of [Exact exchange, dispersion correction](exx-vdw.md)).
- [Convention] The definitions and forces of `Uhub` (DFT+U occupation definitions: onsite, full, dual), `Ucs`, `Uzs`, `Uzo` are given in
  [DFT+U, constrained DFT and magnetic features](dftu-magnetism.md) (for the noncollinear case, [Energy, forces and stress for noncollinear and spin-orbit calculations](noncollinear-energy-force.md)),
  and `Uch` is defined in Sec. 5 of [ESM, external field, core levels](esm-field-corehole.md).
- [Convention] `Uef` is the grid sum of the product of the electron density and the field potential; it does not include the interaction between nuclei (core charges) and the field.

## Consistency identities

### Band energy and double counting [Convention]

The Hamiltonian matrix elements are

$$
H^{\sigma} = T + V^{\mathrm{na}} + V^{\mathrm{NL}} + \langle \phi \rvert \delta V\_H + v\_{\mathrm{xc},\sigma} \lvert \phi \rangle\_{\mathrm{grid}} + \cdots
$$

(the last term is the matrix element built from a sum on the grid), so for \\( \rho \\) built from the same eigenvectors,
\\( U\_{\mathrm{ele}} = g\_s\sum\_\sigma\sum\_{\mathrm{pairs}} \rho^{\sigma} H^{\sigma} \\) holds identically.
Using this, when SCF has converged,

$$
U\_{\mathrm{tot}} = U\_{\mathrm{ele}}
{}- \Delta V\sum\_p n\\,\delta V\_H + U\_{\mathrm{H1}}
{}- \Delta V\sum\_p \sum\_\sigma n\_\sigma v\_{\mathrm{xc},\sigma} + U\_{\mathrm{xc0}} + U\_{\mathrm{xc1}}
{}+ U\_{\mathrm{core}} + U\_{\mathrm{H0}} + \cdots
$$

(Harris type, double-counting form) holds. For `Cluster-LNO`, `Utot` is built in this form, and
\\( -\Delta V\sum n\\,\delta V\_H \\) is output as `UH_dc` and \\( -\Delta V\sum n\_\sigma v\_{\mathrm{xc},\sigma} \\) as `Uxc0_dc`, `Uxc1_dc`.
The double-counting terms use only the valence density \\( n\_\sigma \\) (because \\( n^{\mathrm{pcc}} \\) does not contribute to the Hamiltonian eigenvalues).
The double-counting Hartree term uses \\( n \\) rather than \\( \delta n \\) because \\( \delta V\_H \\) in the Hamiltonian acts on all valence electrons.

During SCF (when the density used to build the Hamiltonian differs from the density built from its eigenvectors), the two expressions do not agree.
Tests compare them after sufficient convergence.

### Other identities [Design]

- A single isolated atom (a box large enough that periodic images do not enter the atom-pair range): \\( U\_{\mathrm{core}} = 0 \\),
  \\( U\_{\mathrm{H0}} = -\frac12\int n^{a}\_I V^{a}\_{H,I} \\).
- At the edge of the atom-pair range: the screened pair energy is 0 by the definition of \\( s\_{IJ} \\).
- \\( \Delta V\sum\_p \delta n(\mathbf r\_p) \approx 0 \\) (since the system is neutral; within grid errors).

## Forces [Public]

### Sign and units [Convention]

- The force is \\( \mathbf F\_K = -\partial U\_{\mathrm{tot}}/\partial\mathbf R\_K \\) (Hartree/bohr).
  The `<coordinates.forces` block of the output file lists the atom index, element symbol, coordinates (Å), and the three force components (Hartree/bohr).
  The MD history file (`.md`) also outputs forces with the same sign.
- In geometry optimization and `MD.Type NOMD`, no correction making the sum of forces zero is applied.
  In molecular dynamics (`NVE`, `NVT_VS`, `NVT_VS2`, `NVT_VS4`, `NVT_NH`, `NVT_Langevin`), before time evolution,
  the mean force is subtracted component-wise over unfixed atoms so that the sum becomes 0.
  [Convention] The output forces are the values before correction ([Geometry optimization and molecular dynamics](geometry-md.md) Sec. 10). The `.md` file contains forces before correction, and
  the `<coordinates.forces` block of `.out` is written only for `NOMD` and geometry optimization, where no correction is applied, so corrected forces never appear in the output.
- Forces are evaluated with the final SCF density matrix as analytic derivatives of the total energy.

### Overall form [Public]

The derivative with respect to \\( \mathbf R\_K \\) splits as follows, following TN1 Eq. 59–81.
Contributions through changes of the density matrix all combine into derivatives of the overlap matrix (the Pulay term) by using the derivative of the orthonormality condition of the eigenvectors \\( c^\dagger S c = 1 \\),
\\( \partial c^\dagger S c + c^\dagger S \partial c = -c^\dagger (\partial S) c \\) (TN1 Eq. 73).

$$
\frac{\partial U\_{\mathrm{tot}}}{\partial \mathbf R\_K}
{}= \underbrace{- g\_s\sum\_\sigma\sum\_{\mathrm{pairs}} \varepsilon^{\sigma}\\, \frac{\partial S}{\partial\mathbf R\_K}}\_{\text{Pulay}}
{}+ \underbrace{g\_s\sum\_\sigma\sum\_{\mathrm{pairs}} \rho^{\sigma} \frac{\partial (T + V^{\mathrm{na}} + V^{\mathrm{NL}})}{\partial\mathbf R\_K}}\_{\text{matrix elements}}
{}+ \underbrace{\frac{\partial U\_{\mathrm{grid}}}{\partial\mathbf R\_K}\Big\rvert\_{\rho}}\_{\text{grid}}
{}+ \frac{\partial (U\_{\mathrm{core}} + U\_{\mathrm{H0}})}{\partial\mathbf R\_K}
{}+ \cdots
$$

The sum of the first two terms, the derivative of the matrix elements at fixed \\( \rho \\) and the Pulay term,
corresponds to the formula for the derivative of \\( \sum f e \\), \\( \partial \sum\_n f\_n e\_n = \mathrm{Tr}[\rho\\, \partial H] - \mathrm{Tr}[\varepsilon\\, \partial S] \\) (at fixed \\( f \\)).
This formula has been confirmed numerically for the generalized eigenvalue problem in `tools/check_energy_force.py`.

### Derivatives of matrix elements [Public]

- Two-center integrals \\( O(\boldsymbol\tau) \\) (\\( S \\), \\( T \\)) are functions of \\( \boldsymbol\tau = \mathbf R\_J + \mathbf L - \mathbf R\_I \\) only, so
  \\( \partial O/\partial\mathbf R\_I = -\nabla\_{\boldsymbol\tau} O \\), \\( \partial O/\partial\mathbf R\_J = +\nabla\_{\boldsymbol\tau} O \\), and 0 for all other atoms.
  \\( \nabla\_{\boldsymbol\tau} O \\) uses the analytic derivatives of Sec. 6 of [Two-center integrals](two-center-overlap-kinetic.md).
  On-site pairs (\\( I = J \\) with \\( \boldsymbol\tau = 0 \\)) do not depend on the configuration and do not contribute.
- The nonlocal term \\( V^{\mathrm{NL}}\_{I\alpha,J\beta} = \sum\_K\sum\_{ij}\langle\phi\_{I\alpha}\vert\beta^K\_i\rangle D^K\_{ij}\langle\beta^K\_j\vert\phi\_{J\beta}\rangle \\)
  depends on the positions of three atoms \\( I, J, K \\). Each factor is a two-center integral, so

$$
\frac{\partial V^{\mathrm{NL}}\_{I\alpha,J\beta}}{\partial \mathbf R\_K}
{}= \sum\_{ij}\left[\frac{\partial\langle\phi\_{I\alpha}\vert\beta^K\_i\rangle}{\partial\mathbf R\_K} D^K\_{ij}\langle\beta^K\_j\vert\phi\_{J\beta}\rangle
{}+ \langle\phi\_{I\alpha}\vert\beta^K\_i\rangle D^K\_{ij}\frac{\partial\langle\beta^K\_j\vert\phi\_{J\beta}\rangle}{\partial\mathbf R\_K}\right]
$$

  i.e. the derivative of the two-center integral is used factor by factor. When \\( K \\) coincides with \\( I \\) or \\( J \\), the same rule is applied factor by factor.
  Translational invariance gives \\( \partial/\partial\mathbf R\_I + \partial/\partial\mathbf R\_J + \partial/\partial\mathbf R\_K = 0 \\), so
  computing two of them yields the third (TN1 Eq. 62).
- When the neutral-atom term is built by projector expansion (TN1 Eq. 32–34, [neutral atom](neutral-atom.md)),
  \\( \langle\phi\vert V^{\mathrm{na}}\_K \bar R\_{l\zeta} Y\_{lm}\rangle \\) is differentiated in the same way as a two-center integral.
  The derivatives for combinations computed directly without the expansion (one-center and two-center cases) are [Open].

### Derivatives of grid terms [Public]

The grid points \\( \mathbf r\_p \\) stay fixed and do not move with the atoms; only the functions attached to atoms (basis functions, \\( n^{a}\_K \\), \\( n^{\mathrm{pcc}}\_K \\),
and \\( V^{\mathrm{na}}\_K \\) in grid mode) move, and the discretized energy is differentiated exactly (TN1 Eq. 63–69, 79–80).
Therefore the analytic forces agree with finite differences of the discretized total energy for any grid spacing.

The derivative of a function \\( g(\mathbf r\_p - \mathbf R\_K) \\) is \\( \partial g/\partial\mathbf R\_K = -\nabla g(\mathbf r\_p - \mathbf R\_K) \\).

1. **Movement of basis functions** (TN1 Eq. 67, 79): the change in density at fixed \\( \rho \\)

$$
\frac{\partial n\_\sigma(\mathbf r\_p)}{\partial\mathbf R\_K}\Big\rvert\_{\rho}
{}= -2 \sum\_{\alpha}\sum\_{J,\mathbf L,\beta} \rho^{\sigma}\_{K\alpha,J\beta}(\mathbf L)\\, \nabla\phi\_{K\alpha}(\mathbf r\_p - \mathbf R\_K)\\, \phi\_{J\beta}(\mathbf r\_p - \mathbf R\_J - \mathbf L)
$$

   (the factor 2 combines, by the symmetry of \\( \rho \\), the cases where \\( K \\) appears on the left or the right)
   is multiplied by the effective potential \\( V^{\sigma}\_{\mathrm{eff}} = \delta V\_H + v\_{\mathrm{xc},\sigma} \\) (\\( + V^{\mathrm{na}} \\) when `Una` is evaluated on the grid,
   \\( + V\_{\mathrm{ef}} \\) when an electric field is present) and summed with \\( g\_s \Delta V\sum\_p\sum\_\sigma \\).
   The indirect change through \\( \delta V\_H \\) combines into this single term because \\( E\_{\delta\mathrm{ee}} \\) is a quadratic form in \\( \delta n \\) (TN1 Eq. 64–65).
2. **Movement of atomic densities** (TN1 Eq. 66, 68, 80): \\( E\_{\delta\mathrm{ee}} \\) also depends on \\( n^{a} \\) through \\( \delta n = n - n^{a} \\), giving
   \\( -\Delta V\sum\_p \delta V\_H(\mathbf r\_p)\\, \partial n^{a}\_K(\mathbf r\_p - \mathbf R\_K)/\partial\mathbf R\_K \\).
3. **Movement of the partial core correction** (TN1 Eq. 69, 80): \\( \frac12\Delta V\sum\_p\sum\_\sigma v\_{\mathrm{xc},\sigma}(\mathbf r\_p)\\, \partial n^{\mathrm{pcc}}\_K/\partial\mathbf R\_K \\).
4. **Movement of the neutral-atom potential in grid mode**: \\( \Delta V\sum\_p n(\mathbf r\_p)\\, \partial V^{\mathrm{na}}\_K(\mathbf r\_p - \mathbf R\_K)/\partial\mathbf R\_K \\).
5. **Derivative of the atomic-density correction \\( \Delta E^{(0)} \\)** [Convention]: the first term gives
   \\( -\Delta V\sum\_p v^{\mathrm{CA}}(n^{a}+n^{\mathrm{pcc}})\\, \partial(n^{a}\_K + n^{\mathrm{pcc}}\_K)/\partial\mathbf R\_K \\)
   (\\( v^{\mathrm{CA}} = d(n\epsilon^{\mathrm{CA}})/dn \\)). In the second term the atom-centered quadrature points move with the central atom, so
   it depends only on the relative position of the center \\( I \\) and the neighbor \\( J \\);

$$
\frac{\partial}{\partial \mathbf R\_J}\int n^{a}\_I(\mathbf x)\\, \epsilon^{\mathrm{CA}}(n^{a}(\mathbf R\_I + \mathbf x))\\, d^3x
{}= -\int n^{a}\_I(\mathbf x)\\, \frac{d\epsilon^{\mathrm{CA}}}{dn}\\, \nabla n^{a}\_J(\mathbf R\_I + \mathbf x - \mathbf R\_J)\\, d^3x
$$

   is added to \\( J \\), and its negative to \\( I \\) (evaluated at the same quadrature points).

For GGA, unless \\( v\_{\mathrm{xc}} \\) is the partial derivative with respect to \\( n\_p \\) of the discretized energy \\( \Delta V\sum\_p f(n\_p, \nabla\_h n\_p) \\) (\\( \nabla\_h \\) is the finite-difference gradient described in the stress section below)
(including the transpose of the finite-difference gradient), expression 1 above does not agree with finite differences.
Whether \\( v\_{\mathrm{xc}} \\) of the existing implementation is discretely exact in this sense is [Open].

### Derivatives of screened interatomic terms [Convention]

- \\( U\_{\mathrm{core}} \\): \\( \partial/\partial\mathbf R\_I \\) is \\( -\sum\_{(J,\mathbf L)} Z\_I Z\_J\\, (\mathbf R\_I - \mathbf R\_J - \mathbf L)/d^3 \\) (\\( d = \lVert\boldsymbol\tau\rVert \\)).
- \\( U\_{\mathrm{H0}} \\): the derivative of \\( J\_{IJ}(d) \\) with respect to \\( d \\) is evaluated on the same cylindrical grid from the interpolation of the radial derivative of \\( V^{a}\_{H,J} \\),
  and multiplied by the unit vector along the pair axis. The coefficient \\( s\_{IJ} \\) is treated as a constant.
  The derivative of \\( J\_{IJ} \\) is obtained with \\( n^{a}\_I \\) fixed and \\( V^{a}\_{H,J} \\) moved; \\( J\_{JI} \\) is handled likewise, and the two are symmetrized.

### Electronic temperature and forces [Design]

Even when the occupations \\( f \\) follow the Fermi–Dirac distribution, the analytic forces above are derivatives of the free energy including the entropy term,
\\( A = U\_{\mathrm{tot}} - T\_e S \\) (\\( S = -g\_s\sum w\_{\mathbf k}\sum\_n [f\ln f + (1-f)\ln(1-f)] \\))
(from stationarity with respect to \\( f \\) and conservation of the electron number). Since they are not derivatives of `Utot` itself,
for metals or high electronic temperatures finite differences of `Utot` and the analytic forces disagree at \\( O(T\_e) \\) (same point as the remark in TN1).
Force tests use insulators or compare with finite differences of \\( A \\).

## Stress [Public]

### Definition [Public]

Following TN23, under a strain \\( \boldsymbol\epsilon \\) (a symmetric 3×3 matrix), position vectors are transformed as \\( \mathbf r \to (1+\boldsymbol\epsilon)\mathbf r \\),
lattice vectors and atomic positions are moved by the same transformation, and the shapes of basis functions, projectors, atomic densities, etc. are unchanged (only their centers move) (TN23 Eq. 5–9).
The first-order coefficient of the energy

$$
\Sigma\_{\gamma\eta} = \frac{\partial U\_{\mathrm{tot}}}{\partial\epsilon\_{\gamma\eta}}\Big\rvert\_{\boldsymbol\epsilon = 0}
$$

is computed (TN23 Eq. 10–11). Reciprocal lattice vectors change by \\( (1+\boldsymbol\epsilon)^{-T} \\), and \\( \mathbf k\cdot\mathbf L \\), \\( \mathbf G\cdot\mathbf r \\) are invariant.
Integrals over the unit cell are multiplied by \\( \det(1+\boldsymbol\epsilon) \\), so grid sums produce a "volume term" \\( \delta\_{\gamma\eta}\times \\) (that term).

The strain derivative of a function of atomic positions \\( g(\boldsymbol\tau) \\) is, by the chain rule,

$$
\frac{\partial g}{\partial\epsilon\_{\gamma\eta}} = \frac{\partial g}{\partial\tau\_\gamma}\\,\tau\_\eta
$$

(TN23 Eq. 14–15). Since \\( \boldsymbol\epsilon \\) is symmetric, \\( \Sigma \\) is symmetrized at the end [Convention].

### Contributions of each term [Public]

| Term | \\( \partial/\partial\epsilon\_{\gamma\eta} \\) | Source |
|---|---|---|
| Matrix elements (\\( T \\), \\( V^{\mathrm{NL}} \\), projector-expansion \\( V^{\mathrm{na}} \\)) | \\( g\_s\sum\_\sigma\sum\_{\mathrm{pairs}}\rho^{\sigma}\sum\_{\text{factors}} \frac{\partial O}{\partial\tau\_\gamma}\tau\_\eta \\). Three-center terms use, for each two-center factor, that factor's interatomic vector | TN23 Eq. 16–31 |
| Pulay (overlap) | \\( -g\_s\sum\_\sigma\sum\_{\mathrm{pairs}}\varepsilon^{\sigma}\frac{\partial S}{\partial\tau\_\gamma}\tau\_\eta \\) | TN23 Eq. 62–67 |
| \\( E\_{\delta\mathrm{ee}} \\) | \\( \delta\_{\gamma\eta}E\_{\delta\mathrm{ee}} + \int \frac{\partial\delta n}{\partial\epsilon\_{\gamma\eta}}\delta V\_H + 4\pi\Omega\sum\_{\mathbf G\ne 0}\lvert\delta n(\mathbf G)\rvert^2\frac{G\_\gamma G\_\eta}{G^4} \\) | TN23 Eq. 32–45 |
| \\( E\_{\mathrm{xc}} \\) (LDA) | \\( \delta\_{\gamma\eta}E\_{\mathrm{xc}} + \sum\_\sigma\int v\_{\mathrm{xc},\sigma}\left(\frac{\partial n\_\sigma}{\partial\epsilon\_{\gamma\eta}} + \frac12\frac{\partial n^{\mathrm{pcc}}}{\partial\epsilon\_{\gamma\eta}}\right) \\) | TN23 Eq. 46–48 |
| \\( E\_{\mathrm{xc}} \\) (additional GGA part) | \\( -\Delta V\sum\_\sigma\sum\_p \frac{\partial f}{\partial\lVert\nabla n\_\sigma\rVert}\frac{\partial\_\gamma n\_\sigma\\, \partial\_\eta n\_\sigma}{\lVert\nabla n\_\sigma\rVert}\Big\rvert\_{\mathbf r\_p} \\) | TN23 Appendix A |
| \\( E\_{\mathrm{na}} \\) in grid mode | \\( \delta\_{\gamma\eta}E\_{\mathrm{na}} + \int V^{\mathrm{na}}\frac{\partial n}{\partial\epsilon\_{\gamma\eta}} + \int n \frac{\partial V^{\mathrm{na}}}{\partial\epsilon\_{\gamma\eta}} \\) | Same idea as TN23 Sec. 4.2 |
| \\( U\_{\mathrm{core}} \\) | \\( -\frac12\sum\_{I,(J,\mathbf L)} Z\_I Z\_J\frac{\tau\_\gamma\tau\_\eta}{d^3} \\) | TN23 Eq. 57–59 |
| \\( U\_{\mathrm{H0}} \\) | \\( -\frac12\sum s\_{IJ}\frac{\partial J\_{IJ}}{\partial\tau\_\gamma}\tau\_\eta \\) (symmetrized) | TN23 Eq. 60–61 |
| `UpV` | \\( P\Omega\\,\delta\_{\gamma\eta} \\) (\\( P \\) is `MD.applied.pressure`) | [Convention] |

- The strain derivative of a function on the grid is taken by pulling the function back to the original grid, \\( \partial g((1+\boldsymbol\epsilon)(\mathbf r - \mathbf t))/\partial\epsilon\_{\gamma\eta} = \partial\_\gamma g\\,(\mathbf r - \mathbf t)\_\eta \\)
  (TN23 Eq. 36, 38–40). \\( \partial n/\partial\epsilon \\) also includes the change of \\( \rho \\), but as for forces this combines into the Pulay term.
- The last term of \\( E\_{\delta\mathrm{ee}} \\) comes from the strain derivative of \\( 1/\lvert\mathbf G\rvert^2 \\), \\( 2G\_\gamma G\_\eta/G^4 \\) (TN23 Eq. 43).
  \\( \delta n(\mathbf G) = \Omega^{-1}\int\_\Omega \delta n\\, e^{-i\mathbf G\cdot\mathbf r} \\).
- GGA: the existing implementation obtains the gradient at grid points from central differences along the lattice-vector directions, \\( \Delta\_i n = n(\mathbf r\_p + \mathbf h\_i) - n(\mathbf r\_p - \mathbf h\_i) \\)
  (\\( \mathbf h\_i = \mathbf a\_i/N\_i \\)), by solving the 3×3 linear system \\( \mathbf h\_i\cdot\nabla n = \frac12\Delta\_i n \\) (TN23 Appendix A, Eq. 73–77).
  With this finite-difference gradient, applying strain with grid values fixed changes the gradient to \\( (1+\boldsymbol\epsilon)^{-T}\nabla n \\), which gives the additional GGA part in the table.
  The part due to changes of grid values combines into \\( v\_{\mathrm{xc}} \\) including the transpose of the finite-difference gradient (TN23 Eq. 91–93).
- The exchange-correlation volume term does not include the atom-centered quadrature part of the atomic-density correction (the atom-centered part depends only on the relative positions of pairs, so
  like other pair terms it is treated with \\( \partial/\partial\tau\_\gamma\\,\tau\_\eta \\)) [Convention].

> **Notes on TN23**:
> - The volume term in Eq. 34 and Eq. 44 is printed as \\( \delta\_{\gamma\eta}\int\delta n\\,\delta V\_H \\), but
>   the change of the integration domain of \\( E\_{\delta\mathrm{ee}} = \frac12\int\delta n\\,\delta V\_H \\) gives \\( \delta\_{\gamma\eta}E\_{\delta\mathrm{ee}} \\) (with the factor \\( \frac12 \\)).
>   For a periodic array of hard Gaussian charges, the finite difference of the strained energy agrees with the sum of the three terms to \\( 10^{-10} \\) only
>   with the factor \\( \frac12 \\) (`tools/check_energy_force.py`). The existing implementation also adds `UH1` once to the diagonal [Convention].
> - The numerator \\( Z\_I - Z\_J \\) in Eq. 57–59 is considered a typo for \\( Z\_I Z\_J \\). Also, since \\( \partial(1/\lvert\mathbf t\rvert)/\partial\epsilon\_{\gamma\eta} = -t\_\gamma t\_\eta/\lvert\mathbf t\rvert^3 \\),
>   the right-hand side of Eq. 59 needs a minus sign. The table above includes the minus sign (the existing implementation also has the minus sign [Convention]).

### Output and gradient with respect to lattice vectors [Convention]

- Computed with `scf.stress.tensor on` (default off), or in optimizations that move the cell and in NPT MD.
- Standard output shows the 3×3 matrix \\( \Sigma\_{\gamma\eta}/\Omega \\) (after symmetrization) as "Stress tensor (Hartree/bohr^3)".
  The sign of \\( \Sigma \\) is such that the diagonal components are positive when the cell wants to shrink (stretched state).
  The pressure is used as \\( P = -\frac13\mathrm{tr}(\Sigma)/\Omega \\) (positive in a compressed state).
- Cannot be combined with orbital optimization, noncollinear spin, or the DFT+U 'full' projector (stops at input checking).
- [Design] For the matrix \\( A \\) whose rows are the lattice vectors (\\( A\_{ij} \\) is the \\( j \\) component of \\( \mathbf a\_i \\)), strain gives \\( A \to A(1+\boldsymbol\epsilon)^T \\), so
  the gradient with respect to the lattice vectors is

$$
\frac{\partial U\_{\mathrm{tot}}}{\partial A\_{ij}} = \sum\_k (A^{-1})\_{ki}\\, \Sigma\_{kj}
$$

  This corresponds to the quantity that manual Sec. 19.2 calls the "gradient with respect to lattice vectors".

[Open]

- In what format stress and the gradient with respect to lattice vectors are written to the output file `System.Name.out`.
- Details of the stress for combinations of \\( V^{\mathrm{na}} \\) without projector expansion, and for collinear DFT+U (dual) ([DFT+U, constrained DFT and magnetic features](dftu-magnetism.md) Sec. 11).

[Convention] Resolved items:

- The applied pressure `UpV` enters the stress only for cell optimization (`OptC*`, `RFC*`), where `MD.applied.pressure.flag` is multiplied and added to the Cartesian diagonal components. NPT does not use `MD.applied.pressure` but `MD.CellPressureControl` ([Geometry optimization and molecular dynamics](geometry-md.md) Secs. 5.1 and 6.7).
- Equations for the vdW stress are in "Forces and stress" of [Exact exchange, dispersion correction](exx-vdw.md) (the two-body term and the term through the D3 coordination number).
- Noncollinear calculations cannot be combined with stress ([noncollinear](noncollinear-energy-force.md) Sec. 4).

## Tests [Design]

### Input conditions

- Make the SCF convergence criterion sufficiently strict (`scf.criterion` around \\( 10^{-10} \\)), because forces and stress have non-variational errors.
- Finite-difference tests are done with insulators (e.g. H₂O molecule, Si crystal) at low electronic temperature, or compared with the free energy \\( A \\).
- Tests that enable only one term (setting the contributions of other terms to 0) are also prepared. This follows the same idea as the force check of the existing implementation
  (`openmx -forcetest`, flags 0–8 select kinetic, nonlocal, neutral atom, difference Hartree, exchange-correlation, electric field, DFT+U; manual Sec. 72).

### Invariants and finite differences

| Type | Content | Tolerance guideline |
|---|---|---|
| Finite difference | Central difference \\( -[U(+)-U(-)]/(2h) \\) of `Utot` moving one atom by \\( \pm 3\times 10^{-4} \\) bohr in each direction vs. analytic force | \\( 10^{-4} \\) Hartree/bohr (all terms); \\( 10^{-5} \\) per term |
| Finite difference | Two- and three-center terms only (`Ukin`, `Unl`, projector-expansion `Una`) and the Pulay term: derivatives at fixed \\( \rho \\), \\( \varepsilon \\) | \\( 10^{-7} \\) |
| Finite difference | `Ucore + UH0` only | \\( 10^{-6} \\) (order of the cylindrical-grid quadrature error) |
| Finite difference | Central difference of `Utot` with strain \\( \pm 10^{-4} \\) vs. \\( \Sigma\_{\gamma\eta} \\) (6 components); off-diagonal with \\( \epsilon\_{\gamma\eta} = \epsilon\_{\eta\gamma} = h/2 \\) | \\( 10^{-4} \\) Hartree (per unit cell) |
| Finite difference | Derivative of \\( E(V) \\) under uniform volume change vs. \\( -P \\) | Same order as above |
| Invariant | Sum of forces: the two-center, three-center, Pulay and `Ucore + UH0` contributions each sum to 0 | \\( 10^{-10} \\) |
| Invariant | Sum of forces: sum of grid contributions (error because the grid does not move with the atoms) | Confirm it decreases with the grid cutoff; about \\( 10^{-3} \\) Hartree/bohr at 150 Ry as a guideline |
| Invariant | Translate all atoms by an integer multiple of the grid spacing: energy and forces invariant | \\( 10^{-10} \\) |
| Invariant | Translate all atoms by an arbitrary amount: the energy change (egg-box effect) decreases with the grid cutoff | Confirm the trend |
| Invariant | Rotation of an isolated molecule: energy invariant within the egg-box effect, forces transform by the same rotation | Within the egg-box effect |
| Invariant | Stress tensor is symmetric | Rounding error |
| Identity | \\( U\_{\mathrm{ele}} = g\_s\sum\rho H \\) | \\( 10^{-10} \\) |
| Identity | Harris-type `Utot` after convergence vs. ordinary `Utot` | \\( 10^{-8} \\) |
| Identity | Isolated atom: `Ucore` = 0, `UH0` = self term | \\( 10^{-10} \\) |
| Identity | At the edge of the atom-pair range, the screened pair energy is 0 | \\( 10^{-10} \\) |
| Analytic | \\( J\_{IJ}(d) \\) for Gaussian densities vs. \\( Z\_I Z\_J\\,\mathrm{erf}(\sqrt\mu d)/d \\) | \\( 10^{-8} \\) (relative) |
| Analytic | \\( \partial\sum fe = \mathrm{Tr}[\rho\partial H] - \mathrm{Tr}[\varepsilon\partial S] \\) for the generalized eigenvalue problem | \\( 10^{-8} \\) |
| Analytic | Hartree stress of periodic hard Gaussian charges | \\( 10^{-8} \\) |

The last three analytic tests are implemented in `tools/check_energy_force.py`.

### Comparison with fixtures

The following quantities are written out from the existing implementation and compared.

- Energies: `Uele`, `Ukin`, `UH0`, `UH1`, `Una`, `Unl`, `Uxc0`, `Uxc1`, `Ucore`, `UvdW`, `Utot` (values in the output file).
  Both on and off for `scf.Exc0.correction.flag`, and both on and off for `scf.ProExpn.VNA`.
- Forces: values of `<coordinates.forces` and per-term contributions (analytic forces for each flag of `-forcetest`).
- Stress: the stress tensor in standard output.
- Intermediate quantities: converged \\( \rho \\), \\( \varepsilon \\) (blocks per atom pair), and \\( \delta V\_H \\), \\( v\_{\mathrm{xc},\sigma} \\), \\( n^{a} \\) on the grid.
  With these, the energy components and forces can be recomputed and compared without solving the eigenvalue problem.
- Targets: H₂O molecule (`Cluster`), Si crystal (`Band`, including stress), bcc Fe (collinear spin, `Band`).

| Quantity | Tolerance guideline (same grid, same basis) |
|---|---|
| Each energy component | \\( 10^{-6} \\) Hartree |
| `Utot` | \\( 10^{-6} \\) Hartree |
| Forces | \\( 10^{-5} \\) Hartree/bohr |
| Stress \\( \Sigma/\Omega \\) | \\( 10^{-7} \\) Hartree/bohr³ |

If the number of grid points or the choice of origin differs from the existing implementation, energy differences widen to the order of the egg-box effect, so
the tolerances will be re-determined once the grid rules of the [grid](grid-poisson-xc.md) chapter are fixed.

## [Open]

- Resolved: energy components and forces for noncollinear spin and spin-orbit coupling are treated in [noncollinear](noncollinear-energy-force.md); `Uhub`, `Ucs`, `Uzs`, `Uzo` in [DFT+U](dftu-magnetism.md); and `Uch` in [ESM, core levels](esm-field-corehole.md).
- Forces of the partial core correction for pseudopotentials with "open core" treatment (where \\( v\_{\mathrm{xc}} \\) of only one spin may be used). For the noncollinear case, Sec. 3.4 of [noncollinear](noncollinear-energy-force.md).
- Whether the GGA \\( v\_{\mathrm{xc}} \\) is a discretely exact derivative including the transpose of the finite-difference gradient.
- Whether the reimplementation should also follow the choice of not including the partial-core-correction density in the atom-centered quadrature of the atomic-density correction.
- Derivatives of matrix elements for neutral-atom combinations that do not use the projector expansion.
- Format of stress and the lattice-vector gradient in the output file.
