# Real-space grid, Poisson, and exchange-correlation

> Status: prototype (awaiting review). The formulas rewrite technical note 1 (Ver. 1.3, Eq. 15, 21, 47–58, 63–80), technical note 2,
> and Sections 11 and 14 of the manual in the notation of this note; rules not described in the public documents are supplemented as [Convention] from the behavior of the existing implementation.
> Numerical checks are done in `tools/check_grid_poisson_xc.py`.

This chapter covers the construction of the real-space grid that divides the unit cell uniformly, the representation of the density on it,
the Poisson equation for the deviation \\( \delta n \\) from the atomic density, the exchange-correlation (XC) potential,
and the matrix elements from the potential on the grid. For the neutral-atom potential \\( V\_{\mathrm{na}} \\), see
[Neutral-atom potential](neutral-atom.md); for the energy decomposition and forces overall, see [Total energy and forces](energy-force.md).

## 1. Real-space grid

### 1.1 Grid points [Public]

Following technical note 1 (Eq. 48–50), each lattice vector \\( \mathbf a\_i \\) (\\( i = 1,2,3 \\)) is divided into \\( N\_i \\) equal parts.

$$
\mathbf h\_i = \frac{\mathbf a\_i}{N\_i},\qquad
\Delta V = \lvert \mathbf h\_1\cdot(\mathbf h\_2\times\mathbf h\_3)\rvert = \frac{\Omega}{N\_1N\_2N\_3},\qquad
\mathbf r\_{\mathbf n} = \mathbf O + \sum\_i n\_i \mathbf h\_i,\quad n\_i = 0,\dots,N\_i-1
$$

\\( \Omega \\) is the unit-cell volume and \\( \mathbf O \\) is the grid origin.
The reciprocal vectors of \\( \mathbf h\_i \\) are \\( \mathbf g\_i = N\_i \mathbf b\_i \\) (\\( \mathbf b\_i \\) are the reciprocal lattice vectors of the unit cell,
\\( \mathbf a\_i\cdot\mathbf b\_j = 2\pi\delta\_{ij} \\)).
Functions on the grid are treated as cyclic with period \\( N\_i \\); the value at a periodic-image position \\( \mathbf r\_{\mathbf n} + \mathbf L \\) equals the value at \\( \mathbf r\_{\mathbf n} \\).

### 1.2 Determination of the division numbers \\( N\_i \\) [Convention]

From the value \\( E \\) of `scf.energycutoff` (Ry, default 220), the target division number for each direction,

$$
N\_i^{\ast} = \frac{2\sqrt{E}}{\lVert\mathbf b\_i\rVert}
$$

is computed, and \\( N\_i \\) is the integer of the form \\( 2^p3^q5^r7^s \\) (\\( p,q,r,s \ge 0 \\)) that minimizes \\( \lvert \ln N - \ln N\_i^{\ast}\rvert \\).
Since the nearest number is chosen on a logarithmic scale, \\( N\_i < N\_i^{\ast} \\) can occur (the actual cutoff becomes smaller than the specified value).

The "actual cutoff" for each direction is

$$
E\_i = \left(\frac{\lVert\mathbf g\_i\rVert}{2}\right)^2 = \left(\frac{N\_i\lVert\mathbf b\_i\rVert}{2}\right)^2 \quad(\mathrm{Ry})
$$

where \\( \lVert\mathbf g\_i\rVert/2 \\) is the maximum wavenumber the grid can represent in that direction (the Nyquist wavenumber).
For a rectangular-parallelepiped unit cell, \\( N\_i^{\ast} = \sqrt{E}\\,\lVert\mathbf a\_i\rVert/\pi \\), and the grid spacing is about \\( \pi/\sqrt{E} \\) bohr
(0.257 bohr at 150 Ry). Examples (output of `tools/check_grid_poisson_xc.py`):

| Unit cell (bohr) | \\( E \\) (Ry) | \\( (N\_1,N\_2,N\_3) \\) | \\( E\_i \\) (Ry) |
|---|---|---|---|
| cubic 10 | 150 | (40, 40, 40) | 157.91 |
| cubic 10 | 220 | (48, 48, 48) | 227.40 |
| rectangular 5 × 7.3 × 21 | 200 | (24, 32, 96) | 227.40, 189.65, 206.25 |
| fcc (\\( a = 5.4 \\)) | 220 | (15, 15, 15) | 228.46 |

> **Note on technical note 1 Eq. 47**: That equation writes the per-direction cutoff as \\( \tfrac12\lVert\mathbf g\_i\rVert^2 \\), but
> the value adopted and output by the existing implementation is the \\( E\_i = (\lVert\mathbf g\_i\rVert/2)^2 \\) Ry above (\\( \tfrac18\lVert\mathbf g\_i\rVert^2 \\) in Hartree), a factor of 4 different.
> In the reimplementation, the determination of \\( N\_i \\) is defined by the rule above, and Eq. 47 is referenced only for its intent of "making the cutoffs of the three directions equal".

[Public] Section 14.1 of the manual states that 150 to 200 Ry is adequate in many cases, and that 300 Ry or more may be needed when the potential surface is flat
or for small energy differences between spin configurations.

### 1.3 Fixing the division numbers (`scf.Ngrid`)

[Public] Section 14.2 of the manual: giving `scf.Ngrid n1 n2 n3` fixes \\( N\_i \\) to those values, and `scf.energycutoff` is not used.
This is for avoiding discontinuities due to changes of \\( N\_i \\) when drawing an energy curve while varying the lattice constant.

[Convention]

- Each value of `scf.Ngrid` must be of the form \\( 2^p3^q5^r7^s \\); otherwise it is an input error. When all three are 0 (default), nothing is fixed.
- When fixed, the requested cutoff for reporting is replaced by \\( (E\_1+E\_2+E\_3)/3 \\).
- In lattice optimization (cell optimization), the \\( N\_i \\) determined first are kept fixed thereafter.

### 1.4 Grid origin [Convention]

By default, the origin is placed so that the centroid of the grid points coincides with the arithmetic mean of the atomic positions.

$$
\mathbf O = \frac{1}{N\_{\mathrm{at}}}\sum\_I \mathbf R\_I - \sum\_i \frac{N\_i-1}{2}\\,\mathbf h\_i
$$

- The origin is determined only once at the start of the calculation (the first step of MD or geometry optimization), and is not changed thereafter even when atoms move.
- Giving `scf.fixed.grid x y z` (bohr) fixes the origin to that value ([Public] Section 14.3 of the manual:
  to keep the relative position of atoms and grid in a series of calculations and make grid-induced errors consistent; the value is the origin printed to standard output, used as is).
- [Convention] With ESM (`ESM.switch`), the determination of the origin changes for the normal-direction component (the arithmetic mean of the atomic positions is truncated to an integer multiple of the grid spacing before being put into the formula above).
  The rule is in Section 2.2 of [ESM, external field, and core levels](esm-field-corehole.md).

### 1.5 Box for isolated systems (clusters) [Convention]

When the input has no `<Atoms.UnitVectors` block, a rectangular box is generated automatically for an isolated system
(this is an error for band calculations, fractional coordinates, and ESM).

1. For each atom \\( I \\), using the cutoff radius \\( r\_{c,I} \\) written in the PAO file of its species (`radial.cutoff.pao`),
   compute \\( \min\_I (R\_{I\alpha} - r\_{c,I}) \\) and \\( \max\_I (R\_{I\alpha} + r\_{c,I}) \\) for each Cartesian component \\( \alpha \\).
2. Set the edge length to 1.10 times that difference, and take the lattice vectors parallel to the Cartesian axes.
3. With grid spacing \\( h = \pi/\sqrt{E + 0.001} \\) (\\( E \\) in Ry), find the smallest integer \\( j \\) of the form \\( 2^p3^q5^r7^s \\) that is at least \\( \lfloor L/h\rfloor \\) and satisfies \\( jh \ge L \\),
   and enlarge the edge length to \\( jh \\).

A box built in this way gives exactly \\( N\_i = j \\) under the rule of Section 1.2 (confirmed by the script).
The box is treated as is with periodic boundary conditions. For special treatment of the Poisson equation, see Section 3.3.

## 2. Density on the grid

### 2.1 Electron density [Public]

Technical note 1 (Eq. 6, 15): from the basis functions \\( \phi\_{I\alpha}(\mathbf r) = R\_{l\zeta}(r) S\_{lm}(\hat{\mathbf r}) \\) ([PAO basis](pao-basis.md)) and
the density matrix \\( \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L) \\),

$$
n\_\sigma(\mathbf r) = \sum\_{I\alpha}\sum\_{J\beta,\mathbf L} \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\,
\phi\_{I\alpha}(\mathbf r - \mathbf R\_I)\\, \phi\_{J\beta}(\mathbf r - \mathbf R\_J - \mathbf L)
$$

is evaluated directly at the grid points \\( \mathbf r\_{\mathbf n} \\). \\( n = n\_\uparrow + n\_\downarrow \\).

[Convention]

- For each periodic image of each atom, a list of grid points satisfying \\( \lVert \mathbf r\_{\mathbf n} - \mathbf R\_I - \mathbf L\rVert \le r\_{c,I} \\) is built
  (\\( r\_{c,I} \\) is the PAO cutoff radius of the basis specification). The contribution of an atom pair \\( (I, J+\mathbf L) \\) is computed only on the grid points in the intersection of these two spheres.
- In the spin-unpolarized case, two components are held with \\( n\_\uparrow = n\_\downarrow = n/2 \\).
- The density on the grid is not renormalized to the number of electrons (the difference between \\( \Delta V\sum n \\) and the number of electrons is left as is as a discretization error).

### 2.2 Superposition of atomic densities and \\( \delta n \\) [Public]

Technical note 1 (Eq. 6):

$$
n^{(a)}(\mathbf r) = \sum\_{I,\mathbf L} n^{(a)}\_I(\lVert\mathbf r - \mathbf R\_I - \mathbf L\rVert),\qquad
\delta n(\mathbf r) = n(\mathbf r) - n^{(a)}(\mathbf r)
$$

[Convention]

- \\( n^{(a)}\_I \\) is the atomic valence density (radial function) stored in the PAO file. On reading, \\( 4\pi\int r^2 n^{(a)}\_I\\,dr \\) is computed as a sum on the PAO logarithmic mesh,
  and the density is scaled by a constant to match the number of valence electrons (the valence charge \\( Z\_I \\) of the pseudopotential).
- On the grid, values are placed only at the same "grid points with \\( r \le r\_{c,I} \\)" as in Section 1.
- The assignment to spins when \\( n^{(a)} \\) is used as the initial density is treated in [Density matrix and charge mixing](density-mixing.md).

### 2.3 Partial core correction (PCC) density [Public] [Convention]

[Public] Technical note 1 (Eq. 15): for XC, the density with half of \\( n\_{\mathrm{pcc}} \\) added to each spin is used.

$$
E\_{\mathrm{xc}} = \int d^3r\\, \bigl(n(\mathbf r) + n\_{\mathrm{pcc}}(\mathbf r)\bigr)\\,
\varepsilon\_{\mathrm{xc}}\bigl(n\_\uparrow + \tfrac12 n\_{\mathrm{pcc}},\\, n\_\downarrow + \tfrac12 n\_{\mathrm{pcc}}\bigr)
$$

[Convention]

- \\( n\_{\mathrm{pcc}} = \sum\_{I,\mathbf L} n^{\mathrm{pcc}}\_I \\), where \\( n^{\mathrm{pcc}}\_I \\) is the partial core density in the pseudopotential (VPS) file.
  On the grid, as for the atomic density, it is placed only at grid points with \\( r \le r\_{c,I} \\).
- `scf.partialCoreCorrection` must always be on; off is an error (for pseudopotentials without PCC, \\( n^{\mathrm{pcc}}\_I = 0 \\)).
- In spin-polarized calculations, \\( +1 \\), \\( -1 \\), or \\( 0 \\) can be specified per species in the `<scf.pcc.opencore` block;
  for species with \\( +1 \\) (\\( -1 \\)), all of \\( n^{\mathrm{pcc}}\_I \\) is added to the up (down) spin instead of half to each (treatment of open-shell cores).

### 2.4 Noncollinear spin [Public] [Convention]

[Public] Technical note 2 (Eq. 3, 22–25, 34–35): each grid point holds a 2×2 density matrix

$$
\underline n(\mathbf r) = \begin{pmatrix} n\_{11} & n\_{12}\\\\ n\_{21} & n\_{22}\end{pmatrix}
{}= \tfrac12\left(n\\,\sigma\_0 + \mathbf m\cdot\boldsymbol\sigma\right),\qquad n\_{21} = \overline{n\_{12}}
$$

which is diagonalized locally to obtain the densities \\( n'\_{\uparrow}, n'\_{\downarrow} \\) in the local spin frame and the magnetization direction \\( (\theta,\varphi) \\).
In the notation of this note,

$$
\mathbf m = (2\\,\mathrm{Re}\\,n\_{12},\\; -2\\,\mathrm{Im}\\,n\_{12},\\; n\_{11}-n\_{22}),\qquad
n'\_{\uparrow,\downarrow} = \tfrac12\left(n \pm \lVert\mathbf m\rVert\right),
$$

and \\( \theta, \varphi \\) are the polar and azimuthal angles of \\( \hat{\mathbf m} \\) (\\( \hat{\mathbf m} = (\sin\theta\cos\varphi, \sin\theta\sin\varphi, \cos\theta) \\)).

[Convention] When \\( \lVert\mathbf m\rVert < 10^{-14} \\), \\( \theta = \pi/2,\ \varphi = 0 \\) is used (in this case \\( n'\_\uparrow = n'\_\downarrow \\), so the result is unaffected).
The angles are obtained in the ranges \\( \theta \in [0,\pi] \\), \\( \varphi \in [-\pi/2, 3\pi/2) \\).

## 3. Hartree potential (Poisson equation by FFT)

### 3.1 Periodic systems [Public]

Technical note 1 (Eq. 52–54): from the discrete Fourier coefficients of \\( \delta n \\),

$$
\delta\tilde n(\mathbf G) = \frac{1}{N\_1N\_2N\_3}\sum\_{\mathbf n} \delta n(\mathbf r\_{\mathbf n})\\, e^{-i\mathbf G\cdot\mathbf r\_{\mathbf n}},\qquad
\mathbf G = \sum\_i k\_i \mathbf b\_i
$$

the potential is

$$
\delta\tilde V\_H(\mathbf G) = \frac{4\pi}{\lVert\mathbf G\rVert^2}\\,\delta\tilde n(\mathbf G),\qquad
\delta V\_H(\mathbf r\_{\mathbf n}) = \sum\_{\mathbf G} \delta\tilde V\_H(\mathbf G)\\, e^{i\mathbf G\cdot\mathbf r\_{\mathbf n}}
$$

(Hartree atomic units; \\( \delta V\_H \\) is the potential energy felt by one electron, with the convention that the electron density is positive).
\\( \delta\tilde n \\) is a discrete Fourier coefficient, related to the continuous Fourier transform \\( \tilde f \\) of this note ([Notation](../../../conventions/notation.md)) by
\\( \delta\tilde n(\mathbf G) \approx \tilde{\delta n}(\mathbf G)/\Omega \\).
The phase \\( e^{-i\mathbf G\cdot\mathbf O} \\) due to the origin \\( \mathbf O \\) cancels between the forward and inverse transforms, so \\( \mathbf r\_{\mathbf n} \\) in the exponent may be replaced by \\( \sum\_i n\_i\mathbf h\_i \\).

[Convention]

- The integer \\( k\_i \\) is obtained from the index \\( 0 \le \kappa\_i < N\_i \\) as \\( k\_i = \kappa\_i \\) if \\( \kappa\_i < \lfloor N\_i/2\rfloor \\),
  and \\( k\_i = \kappa\_i - N\_i \\) otherwise. Only the real part of the inverse transform is used.
  - For even \\( N\_i \\), this is the same as the usual symmetric choice (placing the Nyquist component on the \\( -N\_i/2 \\) side), and the result is real.
  - For odd \\( N\_i \\), \\( \kappa\_i = (N\_i-1)/2 \\) maps to \\( k\_i = -(N\_i+1)/2 \\), which differs from the symmetric choice
    (\\( k\_i \in [-(N\_i-1)/2, (N\_i-1)/2] \\)) only in the highest-wavenumber component. Because a pair of conjugate components is assigned different \\( \lVert\mathbf G\rVert \\),
    the inverse transform has an imaginary part, and taking the real part amounts to averaging the kernel over that pair.
    The difference appears only in components at the resolution limit and vanishes as the cutoff is raised (example: edge 10 bohr, \\( N = 15 \\), \\( e^{-3r^2} \\) gives
    a difference of \\( 5\times10^{-5} \\) Ha in \\( E\_{\delta\mathrm H} \\); 0 for \\( N=16 \\)).
  - [Design] The reimplementation adopts this mapping as is (for comparison with fixtures), and provides the symmetric choice as an option.
- The \\( \mathbf G = 0 \\) component is set to 0. That is, the unit-cell average of \\( \delta V\_H \\) is 0.
  When the system is charged (`scf.system.charge`), this amounts to placing a uniform compensating charge.
- The forward transform uses \\( e^{-i\cdots} \\) and the inverse \\( e^{+i\cdots} \\); the normalization \\( 1/(N\_1N\_2N\_3) \\) is on the forward side.

### 3.2 Hartree energy [Public]

Technical note 1 (Eq. 21, 51):

$$
E\_{\delta\mathrm{H}} = \frac12\int\_\Omega \delta n\\, \delta V\_H\\, d^3r \approx \frac{\Delta V}{2}\sum\_{\mathbf n} \delta n(\mathbf r\_{\mathbf n})\\, \delta V\_H(\mathbf r\_{\mathbf n})
{}= \frac{\Omega}{2}\sum\_{\mathbf G\neq 0}\frac{4\pi}{\lVert\mathbf G\rVert^2}\lvert\delta\tilde n(\mathbf G)\rvert^2
$$

The last equality is the discrete Parseval identity and is exact on the grid.

### 3.3 Isolated systems: spherical Coulomb cutoff [Convention]

With `scf.coulomb.cutoff on` (default off), the kernel is replaced by

$$
\delta\tilde V\_H(\mathbf G) = \frac{4\pi}{G^2}\left(1 - \cos(G R\_c)\right)\delta\tilde n(\mathbf G)\quad(\mathbf G\ne 0),\qquad
\delta\tilde V\_H(0) = 2\pi R\_c^2\\,\delta\tilde n(0)
$$

\\( R\_c \\) is half the length of the shortest lattice vector of the unit cell. This is the Fourier transform of the Coulomb kernel truncated at distance \\( R\_c \\),
\\( \Theta(R\_c - r)/r \\); if the charge is confined within a range of radius about \\( R\_c/2 \\),
it gives the isolated-system potential in that range without the influence of periodic images (the \\( \mathbf G=0 \\) value is the \\( G\to0 \\) limit).
By default (off), the periodic formula of Section 3.1 is used as is, even for the automatically generated box of an isolated system (Section 1.5).
The rationale is that \\( \delta n \\) is a difference from neutral atoms, so contributions of higher order than dipole dominate, and the influence of periodic images is suppressed by the vacuum margin of the box.

[Design] The Poisson equation for surface systems by the ESM method (`ESM.switch`) and external electric fields are treated in [ESM, external field, and core levels](esm-field-corehole.md), and the electrode boundary conditions of NEGF in Section 7 of [Electron transport by NEGF](negf.md); they are outside the scope of this chapter.

## 4. Exchange-correlation

### 4.1 List of functionals [Public] [Convention]

[Public] Section 11 of the manual: selected with `scf.XcType`.

| `scf.XcType` | Functional | Spin |
|---|---|---|
| `LDA` | Ceperley–Alder (Perdew–Zunger interpolation, PRB **23**, 5048) | `scf.SpinPolarization off` only |
| `LSDA-CA` | Spin-polarized version of the above | off, on, NC |
| `LSDA-PW` | Perdew–Wang 92 correlation (PRB **45**, 13244) + Slater exchange | off, on, NC |
| `GGA-PBE` | Perdew–Burke–Ernzerhof (PRL **77**, 3865) | off, on, NC |

[Convention] Combining `LDA` with spin polarization (on, NC) is an input error. The input also accepts the values `HSE` and `NLX`, but
[Convention] the grid XC evaluator computes no exchange-correlation for them (the input is accepted but effectively unreachable), so the calculation is not meaningful.
The reimplementation stops as "not supported" ("Status of the existing implementation" in [Exact exchange and dispersion correction](exx-vdw.md)). For some double-counting rules of LDA+U, there is a variant that uses the average density of both spins in the XC evaluation (treating spin polarization as 0);
this is treated in the DFT+U chapter.

### 4.2 Formulas of each functional [Public]

The formulas of the original papers above are used. Notation in this note: \\( r\_s = (3/(4\pi n))^{1/3} \\), \\( \zeta = (n\_\uparrow - n\_\downarrow)/n \\),

$$
f(\zeta) = \frac{(1+\zeta)^{4/3} + (1-\zeta)^{4/3} - 2}{2^{4/3}-2},\qquad
\varepsilon\_x^{\mathrm{unif}}(n) = -\frac34\left(\frac{3}{\pi}\right)^{1/3} n^{1/3}
$$

- Exchange (common to all functionals) uses spin scaling \\( E\_x[n\_\uparrow,n\_\downarrow] = \tfrac12\left(E\_x[2n\_\uparrow] + E\_x[2n\_\downarrow]\right) \\);
  in LDA the per-particle value is \\( \varepsilon\_x = \sum\_\sigma n\_\sigma\\,\varepsilon\_x^{\mathrm{unif}}(2n\_\sigma)/n \\).
  For `LDA` and `LSDA-CA`, writing the same quantity as \\( \varepsilon\_x^{\mathrm{P}} + f(\zeta)(\varepsilon\_x^{\mathrm F} - \varepsilon\_x^{\mathrm P}) \\)
  (\\( \varepsilon\_x^{\mathrm F} = 2^{1/3}\varepsilon\_x^{\mathrm P} \\)) is equivalent.
- PZ81 correlation: for the unpolarized (P) and fully polarized (F) cases, \\( \gamma/(1+\beta\_1\sqrt{r\_s}+\beta\_2 r\_s) \\) for \\( r\_s \ge 1 \\),
  and \\( A\ln r\_s + B + C r\_s\ln r\_s + D r\_s \\) for \\( r\_s < 1 \\); the two are interpolated with \\( f(\zeta) \\).
  Coefficients: P: \\( (\gamma,\beta\_1,\beta\_2,A,B,C,D) = (-0.1423, 1.0529, 0.3334, 0.0311, -0.048, 0.0020, -0.0116) \\),
  F: \\( (-0.0843, 1.3981, 0.2611, 0.01555, -0.0269, 0.0007, -0.0048) \\).
- PW92 correlation: the three coefficient sets of Table I of the original paper (\\( \varepsilon\_c(r\_s,0) \\), \\( \varepsilon\_c(r\_s,1) \\), \\( -\alpha\_c(r\_s) \\)) and
  \\( \varepsilon\_c = \varepsilon\_c(r\_s,0) + \alpha\_c f(\zeta)(1-\zeta^4)/f''(0) + [\varepsilon\_c(r\_s,1)-\varepsilon\_c(r\_s,0)] f(\zeta)\zeta^4 \\).
  [Convention] \\( f''(0) = 1.70992093416137 \\) (exact value \\( 8/(9(2^{4/3}-2)) \\)) is used (libxc's `LDA_C_PW` uses 1.709921).
- PBE: exchange enhancement factor \\( F\_x(s) = 1+\kappa-\kappa/(1+\mu s^2/\kappa) \\), \\( \kappa = 0.804 \\), \\( \mu = \beta\pi^2/3 \\),
  \\( s = \lvert\nabla n\rvert/(2k\_F n) \\) (applied to \\( 2n\_\sigma \\) via spin scaling).
  Correlation is the PW92 \\( \varepsilon\_c \\) plus the gradient correction \\( H(r\_s,\zeta,t) \\), with \\( \beta = 0.06672455060314922 \\),
  \\( \gamma = (1-\ln 2)/\pi^2 \\), \\( t = \lvert\nabla n\rvert/(2\phi k\_s n) \\), \\( \phi = \tfrac12[(1+\zeta)^{2/3}+(1-\zeta)^{2/3}] \\).

It is confirmed by the script that these formulas agree with libxc 7.1.2 to a relative error of \\( 10^{-15} \\) (LDA) and \\( 3\times10^{-7} \\) (PBE; due to the difference in the number of digits of the PW92 coefficients on the libxc side).

### 4.3 Treatment of small densities [Convention]

As thresholds that can affect physical quantities, the existing implementation uses the following (densities in bohr⁻³, values after adding PCC).

| Functional | Rule |
|---|---|
| `LDA` | \\( \varepsilon\_{xc} = v\_{xc} = 0 \\) for \\( n \le 10^{-15} \\) |
| `LSDA-CA` | 0 for \\( n < 10^{-15} \\). If \\( \lvert\zeta\rvert > 1 \\), clamp to \\( \zeta = \pm(1-10^{-15}) \\) |
| `LSDA-PW` | 0 for \\( n < 10^{-14} \\). In correlation, clamp \\( \zeta > 0.99 \\) to 0.99 (replaced by \\( -0.99 \\) only when \\( \zeta < -1 \\)) |
| `GGA-PBE` | At points where the valence density \\( n\_\uparrow + n\_\downarrow < 10^{-14} \\), \\( \varepsilon\_{xc} = v\_{xc} = 0 \\) and the gradient is also 0. Elsewhere, each spin density is raised to at least \\( 0.5\times10^{-14} \\), \\( \lvert\nabla n\rvert \\) and \\( 2\lvert\nabla n\_\sigma\rvert \\) to at least \\( 10^{-14} \\), and \\( \zeta \\) is clamped as in `LSDA-PW` |

Only for `GGA-PBE` is the test made on the valence density before adding PCC. [Design] The reimplementation uses the same thresholds when comparing with fixtures;
otherwise it suffices to confirm by tests that the difference in thresholds affects only points with densities of order \\( 10^{-14} \\) bohr⁻³.

### 4.4 XC energy and potential on the grid [Public]

Technical note 1 (Eq. 55, 58): at each grid point, form \\( \tilde n\_\sigma = n\_\sigma + \tfrac12 n\_{\mathrm{pcc}} \\) and set

$$
E\_{\mathrm{xc}}^{\mathrm{grid}} = \Delta V\sum\_{\mathbf n}\left(n + n\_{\mathrm{pcc}}\right)\varepsilon\_{\mathrm{xc}}(\tilde n\_\uparrow,\tilde n\_\downarrow,\nabla\tilde n\_\uparrow,\nabla\tilde n\_\downarrow),\qquad
v\_{\mathrm{xc},\sigma} = \frac{\partial e\_{\mathrm{xc}}}{\partial \tilde n\_\sigma} - \nabla\cdot\frac{\partial e\_{\mathrm{xc}}}{\partial \nabla\tilde n\_\sigma}
$$

(\\( e\_{\mathrm{xc}} = (n+n\_{\mathrm{pcc}})\varepsilon\_{\mathrm{xc}} \\) is the energy per volume). In LDA there is no second term.

### 4.5 GGA gradients: 27-point polynomial fit [Public] [Convention]

[Public] Section 11 of the manual: GGA density gradients are evaluated in real space by fitting a polynomial through a grid point and its 26 surrounding points (27 points in all).

[Convention] Specifically:

1. Take the 27 points at displacements \\( \mathbf d\_{jkl} = j\mathbf h\_1 + k\mathbf h\_2 + l\mathbf h\_3 \\) (\\( j,k,l\in\\{-1,0,1\\} \\)) from the central grid point.
2. Make a linear combination of the 27 monomials \\( x^p y^q z^s \\) (\\( p,q,s\in\\{0,1,2\\} \\)) in the Cartesian coordinates \\( \mathbf d = (x,y,z) \\) match the values at the 27 points
   exactly (a 27×27 linear system; the coefficient matrix depends only on the grid, so its inverse is built only once).
3. The gradient is given by the first-order coefficients of \\( x, y, z \\).

This method is exact for polynomials of total degree at most 2, and on orthogonal grids coincides with the central difference \\( (f\_{+} - f\_{-})/(2h) \\) in each direction (confirmed by the script).
The density whose gradient is taken is \\( \tilde n\_\sigma \\), including PCC.

The second term of the potential, \\( \nabla\cdot(\partial e\_{\mathrm{xc}}/\partial\nabla\tilde n\_\sigma) \\), is also computed by applying the same fit to the three components of
\\( \partial e\_{\mathrm{xc}}/\partial\nabla\tilde n\_\sigma \\) obtained at each grid point.
Hence the potential depends on values up to ±2 points from the center in each direction.

[Design] FFT-based gradients are also possible, but for comparison with the existing implementation, the 27-point method above is the default.

### 4.6 Noncollinear XC [Public] [Convention]

[Public] Technical note 2 (Eq. 18–19): \\( v\_{\mathrm{xc},\uparrow}, v\_{\mathrm{xc},\downarrow} \\) are computed from the local-frame densities \\( n'\_\uparrow, n'\_\downarrow \\)
and transformed back to the global frame.
With \\( \bar v = \tfrac12(v\_{\uparrow}+v\_{\downarrow}) \\), \\( \Delta v = \tfrac12(v\_\uparrow - v\_\downarrow) \\),

$$
\underline V\_{\mathrm{xc}} = \bar v\\,\sigma\_0 + \Delta v\\,\hat{\mathbf m}\cdot\boldsymbol\sigma
{}= \begin{pmatrix} \bar v + \Delta v\cos\theta & \Delta v\sin\theta\\, e^{-i\varphi}\\\\ \Delta v\sin\theta\\, e^{i\varphi} & \bar v - \Delta v\cos\theta\end{pmatrix}
$$

It is confirmed by the script that \\( \mathrm{Tr}(\underline V\_{\mathrm{xc}}\\,\delta\underline n) = v\_\uparrow\delta n'\_\uparrow + v\_\downarrow\delta n'\_\downarrow \\) (to first order) holds.

[Convention]

- In GGA, the local-frame densities \\( n'\_\uparrow(\mathbf r), n'\_\downarrow(\mathbf r) \\) are regarded as scalar fields on the grid and their gradients are taken as in Section 4.5
  (an approximation that does not include the spatial variation of the magnetization direction in the gradient).
- PCC is added half to each of \\( n'\_\uparrow, n'\_\downarrow \\) (for open-shell-core species, to one side by the rule of Section 2.3).
- The Hartree potential is added only to the diagonal components (11 and 22).

### 4.7 Accurate evaluation of the XC energy of the atomic density (\\( E\_{\mathrm{xc}}^{(0)} \\) correction) [Convention]

The grid sum \\( E\_{\mathrm{xc}}^{\mathrm{grid}} \\) depends on the relative position of atoms and grid (egg-box effect).
By default (`scf.Exc0.correction.flag on`), the existing implementation subtracts from the grid the CA-LDA energy of the superposition of atomic densities \\( \hat n = n^{(a)} + n\_{\mathrm{pcc}} \\),
and adds it back with a fine atom-centered spherical integration.

$$
E\_{\mathrm{xc}} = \Delta V\sum\_{\mathbf n}\Bigl[(n+n\_{\mathrm{pcc}})\varepsilon\_{\mathrm{xc}}(\tilde n\_\uparrow,\tilde n\_\downarrow,\dots)
{}- \hat n\\,\varepsilon^{\mathrm{CA}}\_{xc}(\hat n)\Bigr] + E\_{\mathrm{xc}}^{(0)}
$$

$$
E\_{\mathrm{xc}}^{(0)} = \sum\_I \int\_{\lVert\mathbf r-\mathbf R\_I\rVert < r\_{c,I}} \hat n\_I(\mathbf r - \mathbf R\_I)\\, \varepsilon^{\mathrm{CA}}\_{xc}\bigl(\hat n(\mathbf r)\bigr)\\, d^3r,\qquad
\hat n\_I = n^{(a)}\_I + n^{\mathrm{pcc}}\_I
$$

- \\( \varepsilon^{\mathrm{CA}}\_{xc} \\) is the unpolarized CA-LDA (Section 4.2, threshold \\( 10^{-15} \\)) regardless of `scf.XcType`.
  Since the support of each \\( \hat n\_I \\) lies within its \\( r\_{c,I} \\), \\( \sum\_I \hat n\_I = \hat n \\) makes the second term a partition of \\( \int\hat n\\,\varepsilon^{\mathrm{CA}}(\hat n) \\).
- The integral is performed, centered at atom \\( I \\), with the product of a 150-point Gauss–Legendre rule on the radial interval \\( [0, r\_{c,I}] \\) and an angular Lebedev rule (`scf.Num.Leb.Grid`, default 590 points).
  \\( \hat n(\mathbf r) \\) is the sum of \\( \hat n\_J \\) over atoms \\( J \\) (including periodic images) within interatomic distance \\( r\_{c,I}+r\_{c,J} \\) of \\( I \\).
- Since \\( \hat n \\) does not change during SCF, the potential (Hamiltonian) does not change. Only the energy and forces change.
- When off, \\( E\_{\mathrm{xc}} = E\_{\mathrm{xc}}^{\mathrm{grid}} \\).

[Design] The reimplementation defaults to on, and compares both on and off with fixtures.

## 5. Matrix elements of the grid potential

### 5.1 Effective potential [Public] [Convention]

[Public] Technical note 1 (Eq. 30–31): the local potential handled on the grid is

$$
V\_\sigma^{\mathrm{grid}}(\mathbf r) = \delta V\_H(\mathbf r) + v\_{\mathrm{xc},\sigma}(\mathbf r)\\;\bigl[+\\, V\_{\mathrm{na}}(\mathbf r)\bigr]
$$

[Convention] \\( V\_{\mathrm{na}} \\) is put on the grid only when `scf.ProExpn.VNA off` (with the default on, it is computed separately by projector expansion).
When an external electric field is used, the potential of that uniform field is also added (treated in Section 3 of [ESM, external field, and core levels](esm-field-corehole.md)). In the noncollinear case, \\( \delta V\_H \\) (and \\( V\_{\mathrm{na}} \\)) are added to the diagonal 11 and 22 components,
and the 12 component is \\( \underline V\_{\mathrm{xc}} \\) only.

### 5.2 Matrix elements [Public]

Technical note 1 (Eq. 57–58):

$$
H^{\mathrm{grid},\sigma}\_{I\alpha,J\beta}(\mathbf L) = \Delta V \sum\_{\mathbf n}
\phi\_{I\alpha}(\mathbf r\_{\mathbf n} - \mathbf R\_I)\\, V^{\mathrm{grid}}\_\sigma(\mathbf r\_{\mathbf n})\\, \phi\_{J\beta}(\mathbf r\_{\mathbf n} - \mathbf R\_J - \mathbf L)
$$

(\\( \mathbf r\_{\mathbf n} \\) are the periodically extended grid points). In the noncollinear case, \\( \sigma \\) is replaced by the 2×2 component \\( (s,s') \\).
These matrix elements are added to the other terms of \\( H(\mathbf L) \\) in the [sparse block format](../design/sparse-block-format.md) (kinetic energy, nonlocal, and projector-expanded \\( V\_{\mathrm{na}} \\)).

[Convention]

- The sum is taken only over the "grid points in the intersection of the two cutoff spheres" of Section 2.1. Atom pairs are those with \\( \lVert\mathbf R\_J+\mathbf L-\mathbf R\_I\rVert \le r\_{c,I}+r\_{c,J} \\).
- The values of the basis functions on grid points are obtained by interpolating the radial-function tables and are stored in single-precision floating point. The relative error of the matrix elements is of order \\( 10^{-7} \\).
  [Design] The reimplementation defaults to double precision, and the tolerance for comparison with fixtures allows for this difference.
- There are no devices such as a double grid that uses coarse and fine grids separately. There is a single grid, used for density, Poisson, XC, and matrix elements alike
  (the exception is the atom-centered spherical integration of Section 4.7).

## 6. Contributions to forces and stress (formulas only) [Public] [Convention]

[Public] Technical note 1 (Eq. 63–69, 79–80). The grid-related force terms with respect to the position \\( \boldsymbol\tau\_k \\) of atom \\( k \\) are

$$
\frac{\partial E}{\partial\boldsymbol\tau\_k}\Big|\_{\mathrm{grid}}
{}= \sum\_\sigma\sum\_{\alpha\in k}\sum\_{J\beta,\mathbf L} 2\rho^\sigma\_{k\alpha,J\beta}(\mathbf L)\\,\Delta V\sum\_{\mathbf n}
\frac{\partial\phi\_{k\alpha}(\mathbf r\_{\mathbf n}-\boldsymbol\tau\_k)}{\partial\boldsymbol\tau\_k}\\, V^{\mathrm{grid}}\_\sigma(\mathbf r\_{\mathbf n})\\, \phi\_{J\beta}(\mathbf r\_{\mathbf n}-\mathbf R\_J-\mathbf L)
{}- \Delta V\sum\_{\mathbf n}\delta V\_H\\,\frac{\partial n^{(a)}\_k}{\partial\boldsymbol\tau\_k}
{}+ \Delta V\sum\_{\mathbf n}\sum\_\sigma v\_{\mathrm{xc},\sigma}\\,\frac{\partial n^{\mathrm{pcc}}\_{k,\sigma}}{\partial\boldsymbol\tau\_k}
$$

(\\( n^{\mathrm{pcc}}\_{k,\sigma} \\) is the part assigned to spin \\( \sigma \\) in Section 2.3; the factor 2 comes from the symmetry of the density matrix).
[Convention] The PCC force in noncollinear calculations uses the components of \\( v\_{\mathrm{xc}} \\) in the global frame rather than the local spin frame, so for open-shell-core species it differs from the exact derivative (Section 3.4 of [Noncollinear](noncollinear-energy-force.md)).
The basis-function derivative terms are accompanied separately by a term from the change of \\( \rho \\) (the energy-weighted density matrix term). This is treated in [Total energy and forces](energy-force.md).

[Convention] When the \\( E\_{\mathrm{xc}}^{(0)} \\) correction (Section 4.7) is on, the following is further added:

$$
{}-\Delta V\sum\_{\mathbf n} v^{\mathrm{CA}}\_{xc}(\hat n)\\,\frac{\partial\hat n\_k}{\partial\boldsymbol\tau\_k} + \frac{\partial E^{(0)}\_{\mathrm{xc}}}{\partial\boldsymbol\tau\_k},\qquad
\frac{\partial E^{(0)}\_{\mathrm{xc}}}{\partial\boldsymbol\tau\_k}
{}= \sum\_{I}\int \hat n\_I(\mathbf r-\mathbf R\_I)\\,\frac{d\varepsilon^{\mathrm{CA}}\_{xc}}{dn}(\hat n)\left[\frac{\partial\hat n\_k(\mathbf r-\boldsymbol\tau\_k)}{\partial\boldsymbol\tau\_k}(1-\delta\_{Ik}) - \delta\_{Ik}\sum\_{J\neq k}\frac{\partial\hat n\_J(\mathbf r-\mathbf R\_J)}{\partial\mathbf R\_J}\right] d^3r
$$

(\\( v^{\mathrm{CA}} = d(n\varepsilon^{\mathrm{CA}})/dn \\). The \\( I = k \\) term of the second formula incorporates, via the change of variables \\( \mathbf u = \mathbf r - \boldsymbol\tau\_k \\),
the fact that the integration region moves together with atom \\( k \\)).

Stress: \\( E\_{\delta\mathrm H} \\) is differentiated with respect to lattice strain via the \\( \mathbf G \\)-space formula of Section 3.2, and \\( E\_{\mathrm{xc}} \\) via the displacement of grid points, the change of \\( \Delta V \\), and the GGA gradient terms.
Details are deferred to technical note 23 and [Total energy and forces](energy-force.md) ([Open]).

## 7. Tests [Design]

`tools/check_grid_poisson_xc.py` contains reference implementations of the analytic items below.

| Kind | Content | Tolerance (guideline) |
|---|---|---|
| Analytic | Division-number rule: the 4 examples in the table of Section 1.2, and the isolated-system box of Section 1.5 (example: extent 12, 12, 14.3 bohr, 150 Ry → edges 13.8515, 13.8515, 16.1601 bohr, \\( N = 54, 54, 63 \\)) | exact (integer match) |
| Analytic | Periodic Poisson: in a distorted unit cell, \\( \delta n = \cos(\mathbf G\_0\cdot\mathbf r) + c \\) → \\( \delta V\_H = 4\pi\cos(\mathbf G\_0\cdot\mathbf r)/G\_0^2 \\) (the constant \\( c \\) vanishes) | \\( 10^{-12} \\) |
| Analytic | Spherical Coulomb cutoff: cube of edge 24 bohr, 96³ points, unit-charge Gaussian \\( (\alpha/\pi)^{3/2}e^{-\alpha r^2} \\) (\\( \alpha = 1.5 \\)); \\( \delta V\_H = \mathrm{erf}(\sqrt\alpha r)/r \\) for \\( r < L/2 - 4 \\) from the center, self-energy \\( \sqrt{\alpha/(2\pi)} \\) | \\( 10^{-8} \\) |
| Analytic | Solving the same Gaussian without cutoff (Section 3.1, mean subtracted) gives a self-energy of 0.4296 (isolated value 0.4886): confirms the size of the periodic-image effect | \\( 10^{-4} \\) |
| Analytic | 27-point gradient: exact for gradients of polynomials of total degree at most 2 on a distorted grid; matches central differences on an orthogonal grid | \\( 10^{-10} \\) |
| Analytic | XC: compare \\( \varepsilon\_{xc} \\) and \\( v\_{xc} \\) with the reference values below (libxc 7.1.2, unpolarized, Hartree) | relative \\( 10^{-6} \\) |
| Analytic | Noncollinear: the formula for the 12 component of \\( \underline V\_{\mathrm{xc}} \\), \\( \mathrm{Tr}(\underline V\delta\underline n) = \sum v\_\sigma\delta n'\_\sigma \\) | \\( 10^{-12} \\) |
| Invariant | Charge: \\( \Delta V\sum\_{\mathbf n} n \\) converges to the number of electrons as \\( E \\) is raised. \\( \Delta V\sum\delta n \\) is close to the system charge with reversed sign | about \\( 10^{-4} \\) at 200 Ry (to be decided) |
| Invariant | Translation of the grid origin: shifting all atoms and `scf.fixed.grid` by the same vector leaves the energy invariant (since the grid points move by the same amount as the atoms). The same holds for the default origin (Section 1.4), since the origin moves with the centroid of the atoms | rounding error |
| Analytic | Wavenumber mapping for odd \\( N\_i \\) (Section 3.1): matches the symmetric choice for \\( N=16 \\); differs only in the highest-wavenumber component for \\( N=15 \\) | \\( 10^{-12} \\) |
| Invariant | The energy variation (egg box) when the origin is shifted by fractions of \\( \mathbf h\_i \\) decreases with \\( E \\), and is smaller with the \\( E\_{\mathrm{xc}}^{(0)} \\) correction on than off | trend check |
| Invariant | Forces: agreement between forces including the grid terms and finite differences of the total energy with `scf.fixed.grid` fixed | about \\( 10^{-4} \\) Ha/bohr |
| Comparison | Compare with fixtures of the existing implementation's output (division numbers, actual cutoff, grid origin, \\( E\_{\delta\mathrm H} \\), \\( E\_{\mathrm{xc}} \\), grid matrix elements). Examples: isolated H₂O (LDA, GGA-PBE), Si crystal (GGA-PBE), bcc Fe (LSDA-CA, spin-polarized), noncollinear Fe chain | energy \\( 10^{-6} \\) Ha, matrix elements \\( 10^{-6} \\) Ha (allowing for the single-precision difference; to be decided) |

XC reference values (computed with libxc 7.1.2; `LDA,PZ` is Slater exchange + PZ81, `LDA,PW` is Slater exchange + PW92):

| Functional | Input | \\( \varepsilon\_{xc} \\) | \\( v\_{xc} \\) |
|---|---|---|---|
| PZ | \\( n = 10^{-4} \\) | −0.0495735123 | −0.0644958518 |
| PZ | \\( n = 10^{-2} \\) | −0.1970983191 | −0.2564000609 |
| PZ | \\( n = 1 \\) | −0.8091965677 | −1.0635669021 |
| PZ | \\( n = 10 \\) (\\( r\_s < 1 \\) branch) | −1.6819531872 | −2.2215516600 |
| PW | \\( n = 10^{-2} \\) | −0.1968153660 | −0.2560329456 |
| PW | \\( n = 1 \\) | −0.8097590800 | −1.0642022422 |
| PZ (polarized) | \\( n\_\uparrow = 0.08,\ n\_\downarrow = 0.02 \\) | −0.4167263619 | ↑ −0.5749254057, ↓ −0.4326289466 |
| PW (polarized) | \\( n\_\uparrow = 0.08,\ n\_\downarrow = 0.02 \\) | −0.4175794798 | ↑ −0.5772617587, ↓ −0.4297856767 |
| PBE | \\( n = 0.1,\ \lvert\nabla n\rvert = 0.1 \\) | −0.3969182816 | \\( \partial e/\partial n \\) = −0.5149085316, \\( \partial e/\partial(\lvert\nabla n\rvert^2) \\) = −0.0156917755 |
| PBE | \\( n = 1,\ \lvert\nabla n\rvert = 0.5 \\) | −0.8097745831 | −1.0641456355, −0.0001245773 |

The PBE \\( v \\) follows the libxc convention (\\( \partial e/\partial n \\) and \\( \partial e/\partial\sigma \\), \\( \sigma = \lvert\nabla n\rvert^2 \\)).
The difference between PZ and PW at \\( n = 1 \\) (\\( 6\times10^{-4} \\) Ha) is large enough to detect a mix-up of functionals.

## [Open]

- The exact formulas of each grid-related stress term (correspondence with technical note 23).
- Concrete tolerance values for the comparison fixtures (in particular the matrix-element differences due to single-precision basis-function values).
- The correction for charged isolated systems when `scf.coulomb.cutoff` is not used (whether a Makov–Payne-type correction is applied). Our current understanding is "no correction", but this needs to be confirmed in the energy output.
