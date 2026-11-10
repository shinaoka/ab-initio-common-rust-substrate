# ESM method, external electric field, charge doping, virtual atoms, and core levels

> Status: first draft (awaiting review). The ESM method is based on Otani–Sugino (PRB **73**, 115407 (2006)) and Ohwaki et al. (JCP **136**, 134101 (2012)),
> core levels on technical note 24 (Ozaki–Lee, PRL **118**, 026401 (2017)), and the input on manual Secs. 30–32, 49, and 60;
> everything is rederived in the notation of this note. Rules not described in the public documents are supplemented as [Convention] from the behavior of the existing implementation.
> Numerical checks are done with `tools/check_esm_field_corehole.py`.

## Scope of this chapter

This chapter treats five features added externally to the standard periodic-boundary calculation ([Grid, Poisson, XC](grid-poisson-xc.md), [Total energy and forces](energy-force.md)).

| Feature | Main input | Quantities changed |
|---|---|---|
| Charge doping | `scf.system.charge` | Electron number, treatment of the \\( \mathbf G=0 \\) component of the Hartree potential |
| ESM method | `ESM.switch` etc. | Boundary condition of the Poisson equation in one direction (\\( \delta V\_H \\)) |
| Uniform electric field | `scf.Electric.Field` | Adds a sawtooth potential to the local potential on the grid |
| Virtual atom | PAO and VPS files with fractional `AtomSpecies` | Species data only (the flow of the calculation is unchanged) |
| Core hole (XPS) | `scf.core.hole`, `<core.hole.state` | Penalty term in the Hamiltonian, splitting of the Hartree potential |

Electrode boundary conditions of NEGF are out of scope.

## Notation [Design]

In addition to the symbols of [Notation](../../../conventions/notation.md), [Grid](grid-poisson-xc.md), and [Total energy and forces](energy-force.md), the following are used.

| Symbol | Meaning |
|---|---|
| \\( Q \\) | `scf.system.charge` (excess charge of the system with the charge of one electron taken as \\( -1 \\); positive means hole doping) |
| \\( \delta n = n - n^{a} \\) | Difference between the valence electron density and the superposition of atomic densities ([Grid](grid-poisson-xc.md) Sec. 2.2) |
| \\( \delta V\_H \\) | Hartree potential created by \\( \delta n \\), signed as the potential energy of one electron (\\( \nabla^2\delta V\_H = -4\pi\delta n \\)) |
| \\( \mathbf a\_1 \\), \\( a = \lVert\mathbf a\_1\rVert \\) | Lattice vector taken along the surface normal in ESM, and its length |
| \\( z \\) | Coordinate along the ESM normal. The cell center is 0, \\( -z\_0 \le z \le z\_0 \\), \\( z\_0 = a/2 \\) |
| \\( \mathbf G\_\parallel \\), \\( G\_\parallel = \lVert\mathbf G\_\parallel\rVert \\) | In-plane reciprocal lattice vector |
| \\( g\_m = 2\pi m/a \\) | Wavenumber along the normal (\\( m \\) is an integer) |
| \\( \Delta \\) | Magnitude of the core-hole penalty (Hartree) |
| \\( V\_{\mathrm{ef}} \\) | Sawtooth potential of the uniform electric field (potential energy of an electron) |

---

## 1. Charge doping

### 1.1 Electron number [Public] [Convention]

[Public] Manual Sec. 31: an excess charge is given by `scf.system.charge` (default 0; positive for hole doping, negative for electron doping; fractional values allowed).

[Convention] The number of electrons to occupy is \\( N\_{\mathrm e} = \sum\_I Z\_I - Q \\) (\\( Z\_I \\) is the valence charge of the pseudopotential).
For the determination of the chemical potential see [Density matrix and charge mixing](density-mixing.md) Sec. 1, and for the initial density Sec. 6 of the same chapter
(the electron number of the initial density does not reflect \\( Q \\)).

### 1.2 Compensating charge [Public] [Convention]

[Public] Manual Sec. 31: since the Poisson equation is solved by FFT, the excess charge is neutralized by a uniform background charge of opposite sign.
Therefore total energies of different charge states include the electrostatic energy with the background charge and cannot be compared as they are
(the correct treatment of isolated systems is in manual Sec. 61).

[Convention] In the notation of this note:

- \\( \int\_\Omega\delta n\\,d^3r = N\_{\mathrm e} - \sum\_I Z\_I = -Q \\). The periodic Poisson equation ([Grid](grid-poisson-xc.md) Sec. 3.1) sets
  the \\( \mathbf G = 0 \\) component to 0, so it effectively solves for the density \\( \delta n \\) minus a uniform \\( -Q/\Omega \\).
  `UH1` \\( = \tfrac12\Delta V\sum\_{\mathbf n}\delta n\\,\delta V\_H \\) includes the interaction with this background charge.
- There are no correction terms such as Makov–Payne.
- With `scf.coulomb.cutoff on`, the spherical Coulomb cutoff ([Grid](grid-poisson-xc.md) Sec. 3.3) retains the \\( \mathbf G=0 \\) component as
  \\( 2\pi R\_c^2\\,\delta\tilde n(0) \\), so no background charge is placed. Charged isolated systems are treated this way.
- In the ESM method (Sec. 2), `on2`, `on3`, and `on4` can treat charged systems (the excess charge is canceled by the charge induced on the ideal metal).
- The input keyword `Atoms.NetCharge` is read but not used in the calculation.

---

## 2. ESM method (effective screening medium)

### 2.1 Idea [Public]

Otani–Sugino (2006): outside both ends \\( z = \pm z\_0 \\) of the cell containing a slab along the normal direction (\\( z \\)), either vacuum with permittivity 1 or
an ideal metal with permittivity \\( \infty \\) is placed, and the Poisson equation is solved with periodic boundary conditions in-plane and these boundary conditions along the normal.
For \\( \delta n(\mathbf G\_\parallel, z) \\), the in-plane 2D Fourier transform,

$$
\left(\frac{d^2}{dz^2} - G\_\parallel^2\right)\delta V\_H(\mathbf G\_\parallel, z) = -4\pi\\,\delta n(\mathbf G\_\parallel, z)
$$

is solved with the Green's function \\( g(G\_\parallel; z, z') \\) satisfying the conditions of the media as

$$
\delta V\_H(\mathbf G\_\parallel, z) = \int\_{-z\_0}^{z\_0} g(G\_\parallel; z, z')\\, \delta n(\mathbf G\_\parallel, z')\\, dz'
$$

[Public] Ohwaki et al. (2012): in OpenMX, the pair of nuclei and atomic density \\( n^{a} \\) (the neutral-atom potential \\( V^{\mathrm{na}} \\)) is
short-ranged and the long-range electrostatic interaction is carried by \\( \delta n \\) alone ([Total energy and forces](energy-force.md)), so
the ESM Green's function is applied only to \\( \delta n \\).

### 2.2 Input and cell [Public] [Convention]

[Public] Manual Sec. 49.1:

| Keyword | Value | Default | Meaning |
|---|---|---|---|
| `ESM.switch` | `off`, `on1`, `on2`, `on3`, `on4` | `off` | Medium at lower end / cell / upper end. `on1` = vacuum / vacuum, `on2` = metal / metal, `on3` = vacuum / metal, `on4` = `on2` with a potential difference |
| `ESM.direction` | `x`, `y`, `z` | `x` | Cartesian axis of the normal direction |
| `ESM.potential.diff` | Real (eV) | 0 | Potential difference between the two metals in `on4` |
| `ESM.buffer.range` | Real (Å) | 10.0 | Distance kept between atoms and the medium (see Sec. 2.7) |
| `ESM.wall.switch`, `.position`, `.height` | Sec. 2.6 | | Wall for MD |

- `on1` is restricted to electrically neutral systems (`scf.system.charge` 0). `on2`, `on3`, and `on4` can treat charged systems.

[Convention]

- The normal direction is taken along the **first lattice vector** \\( \mathbf a\_1 \\). \\( \mathbf a\_1 \\) must be parallel to the axis of `ESM.direction` with positive orientation,
  and \\( \mathbf a\_2, \mathbf a\_3 \\) must have zero component along that axis (violations are input errors).
  Input without `<Atoms.UnitVectors` (automatic box for isolated systems) is also an error.
- For atoms given in fractional coordinates, the coordinate along \\( \mathbf a\_1 \\) must be in \\( [0,1] \\) (it is not wrapped as in other directions; it is an error).
- For `on1`, `on2`, `on3`, a nonzero `ESM.potential.diff` is an input error. Charge neutrality for `on1` is not checked.
- **Position of the media**: the lower medium is the plane perpendicular to \\( \mathbf a\_1 \\) through the grid origin, and the upper one is that plane shifted by \\( \mathbf a\_1 \\).
  With ESM, the determination of the grid origin ([Grid](grid-poisson-xc.md) Sec. 1.4) changes as follows. For the component along the normal axis,
  the arithmetic mean \\( c \\) of the atomic positions is replaced by the integer multiple of the grid spacing \\( h = a/N\_1 \\), \\( \mathrm{trunc}(c/h)\\,h \\) (truncation toward 0), and
  then the origin is determined by the formula of Sec. 1.4. The other two components are as in Sec. 1.4.
  Hence the planes of the media are not exactly at \\( x = 0, a \\) but are shifted by up to about a grid spacing depending on the center of mass of the atoms.
  If `scf.fixed.grid` is given, that origin is used.
- The coordinate \\( z \\) is \\( z = x\_\perp - z\_0 \\), where \\( x\_\perp \\) is the "distance from the grid origin" along the normal axis.
  Grid points \\( n\_1 = 0, \dots, N\_1-1 \\) are at \\( x\_\perp = n\_1 h \\); the point \\( z = -z\_0 \\) (lower plane) is a grid point, and the upper plane is a grid point of the periodic image.

[Design] It is recommended to use the usual in-plane k-points and a single k-point along the normal (the existing implementation does not check this).

### 2.3 Green's functions [Public]

The three boundary conditions of Otani–Sugino (2006) are written with the signs of this note (\\( \delta V\_H \\) is the potential energy of an electron; \\( \delta n \\) is the electron density, positive).
The ideal-metal planes are placed at \\( z\_1 = z\_0 \\) (the cell edges).

**`on1` (vacuum | vacuum)**: zero field at \\( z\to\pm\infty \\).

$$
g(G\_\parallel; z,z') = \frac{2\pi}{G\_\parallel} e^{-G\_\parallel\lvert z-z'\rvert},\qquad
g(0; z,z') = -2\pi\lvert z-z'\rvert
$$

The arbitrary constant at \\( G\_\parallel = 0 \\) is fixed by the form above (additive constant 0).

**`on2` (metal | metal, grounded)**: \\( \delta V\_H = 0 \\) at \\( z = \pm z\_1 \\).

$$
g(G\_\parallel; z,z') = \frac{4\pi\\,\sinh\bigl(G\_\parallel(z\_< + z\_1)\bigr)\sinh\bigl(G\_\parallel(z\_1 - z\_>)\bigr)}{G\_\parallel\sinh(2G\_\parallel z\_1)},\qquad
g(0; z,z') = \frac{4\pi\\,(z\_< + z\_1)(z\_1 - z\_>)}{2z\_1}
$$

(\\( z\_< = \min(z,z') \\), \\( z\_> = \max(z,z') \\)).

**`on3` (vacuum | metal)**: vacuum at the lower end (\\( z\to-\infty \\)), grounded metal at the upper end \\( z = z\_1 \\).

$$
g(G\_\parallel; z,z') = \frac{2\pi}{G\_\parallel}\left(e^{-G\_\parallel\lvert z-z'\rvert} - e^{G\_\parallel(z+z'-2z\_1)}\right),\qquad
g(0; z,z') = 4\pi\\,(z\_1 - z\_>)
$$

For `on3` at \\( G\_\parallel = 0 \\), the field on the lower side is 0 and all field lines from the total charge \\( q = \int\delta n \\) terminate on the upper metal
(\\( d\delta V\_H/dz = -4\pi q \\), \\( z \to z\_1 \\)).

**`on4` (metal | metal + potential difference)**: to the solution of `on2`, the charge-independent linear function

$$
V\_{\mathrm{bias}}(z) = \frac{\Delta\Phi}{2}\left(1 - \frac{z}{z\_1}\right),\qquad \Delta\Phi = \texttt{ESM.potential.diff}\ (\mathrm{eV} \to \mathrm{Hartree})
$$

is added to the \\( \mathbf G\_\parallel = 0 \\) component [Convention]. That is, the potential energy of an electron is \\( \Delta\Phi \\) at the lower metal and 0 at the upper one.
Electrons are pulled toward the upper side (the side of lower electron energy). The field magnitude is \\( \Delta\Phi/a \\).

The script confirms that the three Green's functions are symmetric (\\( g(z,z') = g(z',z) \\)) and satisfy their respective boundary conditions.

### 2.4 Numerical method [Convention]

The choices that affect the results of the existing implementation are as follows (the flow of computation is not described).

- **Interpolation**: the phase of \\( \delta n \\) along the normal is shifted to the coordinate system \\( z = x\_\perp - z\_0 \\).
   That is, a trigonometric interpolation on the interval \\( [-z\_0, z\_0) \\), \\( \delta n(\mathbf G\_\parallel, z) = \sum\_m n\_m(\mathbf G\_\parallel)\\, e^{i g\_m z} \\), is built
   (the range of \\( m \\) follows the same mapping as [Grid](grid-poisson-xc.md) Sec. 3.1).
- **Integration**: for each \\( \mathbf G\_\parallel \\), the Green's-function integral is performed **analytically over the interval \\( [-z\_0, z\_0] \\)**.
   Using \\( e^{\pm i g\_m z\_0} = (-1)^m \\), for example for `on1`,

$$
\delta V\_H(\mathbf G\_\parallel, z) = \sum\_m n\_m\left[\frac{4\pi e^{ig\_m z}}{G^2} - \frac{2\pi(-1)^m}{G\_\parallel G^2}
\left((G\_\parallel + ig\_m)\\,e^{G\_\parallel(z - z\_0)} + (G\_\parallel - ig\_m)\\,e^{-G\_\parallel(z + z\_0)}\right)\right],\quad G^2 = G\_\parallel^2 + g\_m^2
$$

$$
\delta V\_H(0, z) = -2\pi n\_0\left(z^2 + z\_0^2\right) + \sum\_{m\ne0} n\_m\left[\frac{4\pi e^{ig\_m z}}{g\_m^2} - \frac{4\pi i(-1)^m z}{g\_m} - \frac{4\pi(-1)^m}{g\_m^2}\right]
$$

   The first term is the periodic solution, and the rest are homogeneous solutions arising from the interval ends (\\( e^{\pm G\_\parallel z} \\), or \\( 1, z \\)).
   Likewise for `on2` and `on3`, the Green's function is split into "free-space term + sum of products \\( e^{\pm G\_\parallel z}e^{\pm G\_\parallel z'} \\)" (polynomials in \\( z, z' \\) for \\( G\_\parallel = 0 \\)) and
   evaluated with elementary integrals such as \\( \int\_{-z\_0}^{z\_0} e^{s G\_\parallel z'} e^{ig\_m z'}dz' = 2(-1)^m\sinh(sG\_\parallel z\_0)/(sG\_\parallel + ig\_m) \\).
- **Edge smoothing**: for the \\( \mathbf G\_\parallel = 0 \\) component only, the values at the 9 points straddling the cell edge \\( z = \pm z\_0 \\) (the edge grid point and 4 points on each side) are
   replaced by a cubic polynomial (Hermite interpolation) matching the values and slopes at the two outermost points (\\( \pm 4h \\) from the edge).
   For the slopes, the derivative of the analytic terms excluding the periodic term (the quadratic in \\( n\_0 \\) and the linear terms) is used.
   This removes the discontinuity that the trigonometric interpolation creates at the edge where the cell touches its periodic image along the normal.
   Since the region near the edge is vacuum (Sec. 2.7), the effect on the energy should be small ([Open]: the magnitude of the effect has not been measured).

[Design] When there is no density at the edges (a buffer region as in Sec. 2.7 exists), the analytic interval integration of the "Integration" item above agrees with the method
"add a homogeneous solution \\( A e^{G\_\parallel z} + B e^{-G\_\parallel z} \\) (\\( A + Bz \\) for \\( G\_\parallel=0 \\)) to the periodic solution \\( V\_{\mathrm{per}} \\)
and determine \\( A, B \\) from the boundary conditions at the edges". The boundary conditions are

| Boundary | \\( G\_\parallel > 0 \\) | \\( G\_\parallel = 0 \\) |
|---|---|---|
| Vacuum (upper end) | \\( V'(z\_0) = -G\_\parallel V(z\_0) \\) | \\( V'(z\_0) = -2\pi q \\) |
| Vacuum (lower end) | \\( V'(-z\_0) = G\_\parallel V(-z\_0) \\) | `on1`: \\( V'(-z\_0) = 2\pi q \\), `on3`: \\( V'(-z\_0) = 0 \\) |
| Metal | \\( V(\pm z\_1) = 0 \\) | \\( V(\pm z\_1) = 0 \\) |

and for `on1` at \\( G\_\parallel = 0 \\) the constant is fixed by \\( V(z\_0) + V(-z\_0) = -4\pi z\_0 q \\) (the convention \\( g = -2\pi\lvert z-z'\rvert \\)).
For the \\( n\_0 \\) term the particular solution \\( -2\pi n\_0 z^2 \\) is used. A reimplementation may use either method; the script confirms that both agree with direct numerical integration.

### 2.5 Energy and forces [Convention]

- `UH1` is \\( \tfrac12\Delta V\sum\_{\mathbf n}\delta n\\,\delta V\_H \\) as in the periodic case, with the ESM solution for \\( \delta V\_H \\).
  Other energy terms are unchanged. The interaction with the charge induced on the ideal metal is included in this \\( \tfrac12\int\delta n\\,\delta V\_H \\).
- For `on4`, \\( V\_{\mathrm{bias}} \\) also enters `UH1` as part of \\( \delta V\_H \\), so the energy contains \\( \tfrac12\int\delta n\\,V\_{\mathrm{bias}} \\).
  On the other hand, \\( V\_{\mathrm{bias}} \\) enters the Hamiltonian and forces in full, so `Utot` for `on4` is not variational with respect to the density.
- Forces are obtained simply by evaluating the grid terms of [Grid](grid-poisson-xc.md) Sec. 6 (including \\( \delta V\_H \\)) with the ESM \\( \delta V\_H \\); there is no ESM-specific term.
  For `on1`–`on3` the Green's functions are symmetric, so energy and forces are consistent.
- Stress is not supported for ESM (it is computed but does not include the ESM contribution).
- Separately, with `MD.Artificial_Force on`, `MD.Artificial_Grad` (Hartree/bohr, default 0) is added to the normal component of the energy gradient of atom 1
  (effective only when ESM is on; for testing).

### 2.6 Wall for MD (`ESM.wall.*`) [Public] [Convention]

[Public] Manual Sec. 49.1: in MD of solid-liquid interfaces, a wall potential is placed to keep atoms from approaching the upper medium.
`ESM.wall.position` (Å) places the wall origin, given as its inward offset from the plane bounding the cell on the upper side, and `ESM.wall.height` is the potential at the point 1.0 Å inside from the wall origin (eV, default 100.0).

[Convention]

- Effective only with `ESM.wall.switch on` (default off). With the normal coordinate of atom \\( I \\) as \\( x\_I \\), the upper plane as \\( x\_b \\) (Sec. 2.2), and
  \\( x\_w = x\_b - \texttt{ESM.wall.position} \\), atoms with \\( d = x\_I - x\_w > 0 \\) receive the force

$$
E\_{\mathrm{wall}}(d) = H\left(\frac{d}{1.89\\,\mathrm{bohr}}\right)^3,\qquad F\_I = -\frac{dE\_{\mathrm{wall}}}{dd} = -\frac{3Hd^2}{(1.89)^3}
$$

  along the normal (\\( H \\) is `ESM.wall.height` in Hartree; 1.89 bohr is used as an approximation of 1 Å).
- **\\( E\_{\mathrm{wall}} \\) is not included in `Utot`** (only the force is added).
- The default of `ESM.wall.position` is 1.0 Å (the note in the manual's input example says 10.0). Negative values are input errors.

### 2.7 `ESM.buffer.range` [Public] [Convention]

[Public] Manual Sec. 49.1: the width (Å, default 10.0) kept between atoms and the medium so that wave functions do not overlap the medium.

[Convention] In this version, `ESM.buffer.range` is read only with `ESM.wall.switch on` and **is not used in the calculation**
(the check whether atoms enter the buffer region is disabled).

[Design] A reimplementation checks that the cutoff sphere of each atom (PAO cutoff radius) does not come within `ESM.buffer.range` of the medium planes,
and issues a warning on violation.

---

## 3. Uniform electric field (`scf.Electric.Field`)

### 3.1 Overview [Public]

Manual Sec. 30: a uniform electric field is represented by a sawtooth potential and can be used in SCF and geometry optimization. Three components are given to `scf.Electric.Field`
in GV/m (\\( 10^9 \\) V/m) (default 0 0 0). In periodic systems, if the discontinuity plane of the sawtooth lies inside the material it can be numerically unstable, whereas
for molecules the discontinuity plane lies in vacuum and usually causes no problem.

### 3.2 Potential [Convention]

- The unit conversion is \\( 1\ \mathrm{GV/m} = 0.1944688\times10^{-2} \\) Hartree/(\\( e \\) bohr) (relative difference \\( 10^{-6} \\) from the CODATA value \\( 1.944690\times10^{-3} \\)).
- The three input components \\( (E\_1, E\_2, E\_3) \\) are treated as components along the lattice vectors \\( \mathbf a\_1, \mathbf a\_2, \mathbf a\_3 \\), respectively
  (the field direction shown in standard output is the unit vector of \\( \sum\_i E\_i\mathbf a\_i/\lVert\mathbf a\_i\rVert \\)).
- The value at grid point \\( \mathbf n = (n\_1,n\_2,n\_3) \\) (\\( 0 \le n\_i < N\_i \\)) is

$$
V\_{\mathrm{ef}}(\mathbf n) = \sum\_{i=1}^{3} E\_i\left(n\_i\lVert\mathbf h\_i\rVert + O\_i - c\_i\right)
$$

  where \\( \mathbf h\_i \\) is the grid step ([Grid](grid-poisson-xc.md) Sec. 1.1), \\( O\_i \\) is the \\( i \\)-th Cartesian component of the grid origin, and
  \\( c\_i \\) is the \\( i \\)-th Cartesian component of the arithmetic mean of the atomic positions.
  For a rectangular cell with lattice vectors parallel to the Cartesian axes, this is \\( V\_{\mathrm{ef}}(\mathbf r) = \mathbf E\cdot(\mathbf r - \mathbf c) \\).
  For a general cell, the expression mixes the distance \\( n\_i\lVert\mathbf h\_i\rVert \\) along \\( \mathbf a\_i \\) with the Cartesian origin components and is not a uniform field.
- \\( V\_{\mathrm{ef}} \\) is the potential energy of an electron and is added to the local potential on the grid
  (in the noncollinear case to the diagonal 11 and 22 components). Electrons are pushed opposite to \\( \mathbf E \\), so \\( \mathbf E \\) is in the usual physical field direction (the direction in which positive charges are pushed).
  [Public] The manual's description "defined with respect to the field acting on electrons" is read in this sense.
- The discontinuity plane of the sawtooth is between \\( n\_i = N\_i - 1 \\) and \\( n\_i = 0 \\) of the next periodic image in each direction, i.e. the cell face through the grid origin. The jump is \\( E\_i a\_i \\), with no smoothing.
- \\( V\_{\mathrm{ef}} \\) is built at the start of each SCF (the start of each MD or geometry-optimization step) with the center of mass of the atoms at that time, and is fixed during the SCF. The grid origin remains fixed from the start of the calculation.
- Simultaneous use with ESM is not prohibited ([Open]: the meaning of the combination has not been verified).

### 3.3 Energy [Convention]

The output label `Uef` is

$$
U\_{\mathrm{ef}} = \Delta V\sum\_{\mathbf n} n(\mathbf r\_{\mathbf n})\\, V\_{\mathrm{ef}}(\mathbf n)
$$

and is included in `Utot` (table in [Total energy and forces](energy-force.md)). It is a term of the electron density only, and **does not include the interaction of the nuclei (valence charges \\( Z\_I \\)) with the field,
\\( -\sum\_I Z\_I V\_{\mathrm{ef}}(\mathbf R\_I) \\)**. Since the center of mass \\( \mathbf c \\) is subtracted, \\( U\_{\mathrm{ef}} \\) contains a constant of the form \\( -N\_{\mathrm e}\\,\mathbf E\cdot\mathbf c \\).

### 3.4 Forces [Convention]

The force consists only of the basis-function-derivative terms of [Grid](grid-poisson-xc.md) Sec. 6 evaluated with the local grid potential including \\( V\_{\mathrm{ef}} \\).
The force on the nuclei \\( Z\_I\mathbf E \\) is not added, and the dependence of the center of mass \\( \mathbf c \\) on the atomic positions is not differentiated.
As a result,

- even for a neutral molecule in a field, the total force does not vanish but is approximately \\( -N\_{\mathrm e}\mathbf E \\) (the force on the electrons);
- on the other hand, since \\( V\_{\mathrm{ef}} \\) moves with the center of mass, \\( U\_{\mathrm{ef}} \\) does not change under rigid translation, and finite differences of `Utot` with respect to atomic positions do not agree with the forces.

[Design] A reimplementation keeps the existing definition as a "compatibility mode" for comparison with fixtures, but by default adds the nuclear term
\\( -\sum\_I Z\_I V\_{\mathrm{ef}}(\mathbf R\_I) \\) to the energy, adds \\( Z\_I\mathbf E \\) to the forces, and replaces the constant shift by the center of mass with a fixed reference point.
Then, for neutral systems, the total energy is translationally invariant and forces and energy are consistent (provided no atoms are near the discontinuity plane).

---

## 4. Virtual atoms (fractional nuclear charge)

### 4.1 Construction [Public]

Manual Sec. 32: a pseudopotential with fractional nuclear charge is generated with ADPACK. In the ADPACK input, `AtomSpecies` (nuclear charge),
`total.electron`, `valence.electron`, and the occupation numbers are given as fractions (e.g. nuclear charge 6.2 between C and N, 4.2 valence electrons).
The basis functions must also be generated for the virtual atom with the same fractional nuclear charge, because the neutral-atom potential is built from the atomic density in the PAO file.
As an example, there is the input `DIA8-VA.dat`, in which one of the 8 C atoms of diamond is replaced by a virtual atom with nuclear charge 6.2.

### 4.2 Treatment in OpenMX [Convention]

There is no special processing for virtual atoms; the values in the PAO and VPS files are used as is. The places where fractional values pass through are as follows.

- `AtomSpecies` of the PAO file and the VPS file must agree (compared as reals; mismatch is an input error).
  The exception is when `AtomSpecies` of the VPS is 0 (empty atom).
- The valence charge of a species is \\( Z\_I = \texttt{valence.electron} + \texttt{AtomSpecies} - \texttt{total.electron} \\) (values from the VPS file) and can be fractional.
  This \\( Z\_I \\) is used for the electron number (Sec. 1.1), the point charges of `Ucore`, and the normalization of the atomic density ([Grid](grid-poisson-xc.md) Sec. 2.2:
  scaled so that \\( 4\pi\int r^2 n^{a}\_I dr = Z\_I \\)).
- Where an integer atomic number is needed (labels etc.), 1 is used if \\( 0 < \texttt{AtomSpecies} < 1 \\), and otherwise the value with the fractional part truncated.

### 4.3 Overriding the electron number per species [Convention]

Although not in manual Sec. 32, giving the `<Number.of.Electrons.Species` block (each line a species name and a real number)
changes the valence charge of a species without creating new files.

- \\( Z\_I \\) of that species is replaced by the given value \\( Z'\_I \\).
- The local potential \\( V^{\mathrm{loc}}\_I \\) of the VPS is multiplied by \\( Z'\_I/Z\_I \\).
- The atomic density is renormalized to \\( Z'\_I \\), from which the atomic Hartree potential and \\( V^{\mathrm{na}}\_I \\) are rebuilt
  (the usual procedure of [Neutral atom](neutral-atom.md)).
- The nonlocal terms, partial core density, and PAO basis are unchanged.

[Open] The status of this feature in the public documents (description in other sections of the manual) has not been confirmed.

---

## 5. Absolute binding energies of core levels (XPS)

### 5.1 Definition of the binding energy [Public]

Technical note 24 (Eqs. 1–4) and manual Sec. 60.1: from the total energy \\( E\_i^{(0)}(N) \\) of the initial state (the \\( N \\)-electron ground state) and
the total energy \\( E\_f^{(0)}(N-1) \\) of the final state with one hole created in a core level (\\( N-1 \\) electrons),

$$
E\_b^{\mathrm{(gas)}} = E\_f^{(0)}(N-1) - E\_i^{(0)}(N),\qquad
E\_b^{\mathrm{(bulk)}} = E\_f^{(0)}(N-1) - E\_i^{(0)}(N) + \mu\_0,\qquad
E\_b^{\mathrm{(metal)}} = E\_f^{(0)}(N) - E\_i^{(0)}(N)
$$

\\( \mu\_0 \\) is the chemical potential of the initial-state calculation (`Chemical Potential (Hartree)` in the output).

- The solid-state formula is derived by setting the vacuum level of the spectrometer to the chemical potential of the sample plus the work function of the spectrometer (the chemical potentials align by electrical contact),
  and noting that an arbitrary shift \\( \Delta\mu \\) of the reference of the chemical potential cancels between \\( E\_i \\) and \\( E\_f \\) as \\( N\Delta\mu - (N-1)\Delta\mu \\).
  Hence neither the vacuum level nor the work function is needed, and the reference is \\( \mu\_0 \\) of the initial state.
- The metal formula follows from Janak's theorem \\( E\_f(N-1) - E\_f(N) = -\mu\_0 \\); it uses the neutralized final state (the hole charge screened by conduction electrons),
  and converges much faster with respect to cell size. It cannot be used for gapped systems (technical note 24: off by 1.2 eV for NH₃).
- \\( E\_f \\) is the DFT energy plus the penalty term (Sec. 5.2) (`Utot` in the output, including `Uch`).

### 5.2 Penalty functional [Public] [Convention]

[Public] Technical note 24 (Eqs. 5–10):

$$
E\_f = E\_{\mathrm{DFT}} + E\_{\mathrm{pen}},\qquad
E\_{\mathrm{pen}} = \frac{1}{V\_B}\int\_B d^3k\sum\_\mu f\_\mu(\mathbf k)\\,\langle\psi\_\mu(\mathbf k)\rvert\hat P\lvert\psi\_\mu(\mathbf k)\rangle,\qquad
\hat P = \lvert R\\,\Phi\_{JM}\rangle\\,\Delta\\,\langle R\\,\Phi\_{JM}\rvert
$$

\\( \Phi\_{JM} \\) is the angular part (two-component spinor) of the Dirac equation in a spherically symmetric potential; for \\( J = l+\tfrac12 \\), \\( M = m+\tfrac12 \\),

$$
\Phi\_{JM} = \sqrt{\frac{l+m+1}{2l+1}}\\,Y\_{lm}\\,\alpha + \sqrt{\frac{l-m}{2l+1}}\\,Y\_{l,m+1}\\,\beta,
$$

and for \\( J = l-\tfrac12 \\), \\( M = m-\tfrac12 \\),

$$
\Phi\_{JM} = \sqrt{\frac{l-m+1}{2l+1}}\\,Y\_{l,m-1}\\,\alpha - \sqrt{\frac{l+m}{2l+1}}\\,Y\_{lm}\\,\beta
$$

(\\( \alpha, \beta \\) are up and down spin). Penalizing the occupation of the target state with a large \\( \Delta \\) creates a hole within the SCF,
and screening by core and valence electrons, spin-orbit coupling, and exchange are included self-consistently. \\( \hat P \\) is added to the Kohn–Sham equation.
Technical note 24 states that \\( \Delta = 100 \\) Ry was used.

[Convention] The concrete form in the existing implementation is as follows.

**Projection functions**: for the atom \\( C \\) where the hole is created (the first value of `<core.hole.state`, the serial number of the atom), at the specified \\( l \\) (the second value `s`, `p`, `d`, `f`),
the **basis functions of the first radial function (\\( \zeta = 0 \\))** \\( \chi\_{\mu} = R\_{l0}(r)S\_{l\mu} \\) (\\( \mu = 0,\dots,2l \\), ordering of the [PAO basis](pao-basis.md)) are used.
Since a pseudopotential including the core as valence and a matching basis are used (e.g. `C_PBE19_1s` and `C7.0_1s-…`, and for the final state `C7.0_1s_CH-…` optimized for the hole),
the first radial function of the \\( l \\) channel is the target core orbital.

**Penalty magnitude**:

$$
\Delta = \texttt{scf.base.penalty.value} - \varepsilon\_{\mathrm{core}}
$$

`scf.base.penalty.value` is in Hartree (default 20.0), and \\( \varepsilon\_{\mathrm{core}} \\) is the \\( (l, \mu=0) \\) value from the eigenvalue table of the PAO file of the species of atom \\( C \\)
that begins with `Lmax=` (each line `l mu` \\( l \\) \\( \mu \\) eigenvalue, in Hartree)
(e.g. 1s of `C6.0_1s.pao` is \\( -10.022 \\) Ha, so \\( \Delta \approx 30.0 \\) Ha). If the table has no value, \\( \varepsilon\_{\mathrm{core}} = 0 \\).
This differs from the "100 Ry" of the technical note.

**Specification of the state (third value \\( j \\))**: \\( 0 \le j \le 2(2l+1) \\).

| Case | \\( j = 1, \dots, 2l+1 \\) | \\( j = 2l+2, \dots, 2(2l+1) \\) | \\( j = 0 \\) |
|---|---|---|---|
| Collinear (`scf.SpinPolarization on`), or noncollinear with \\( l \\)-dependent pseudopotential | Up-spin real orbital \\( \mu = j-1 \\) | Down-spin real orbital \\( \mu = j-2l-2 \\) | All orbitals of the \\( l \\) shell, both spins |
| Noncollinear with \\( j \\)-dependent pseudopotential (`scf.SpinOrbit.Coupling on`) | \\( J = l+\tfrac12 \\), in the order \\( M = J, J-1, \dots, -J \\) (\\( 2l+2 \\) states) | Followed by \\( J = l-\tfrac12 \\), \\( M = J, \dots, -J \\) | All \\( (J,M) \\) (whole shell) |

The boundary between the first and second columns of the table moves to between \\( 2l+2 \\) and \\( 2l+3 \\) only in the \\( j \\)-dependent case (since there are \\( 2l+2 \\) states with \\( J = l+\tfrac12 \\)).
For `s` there is only \\( J = \tfrac12 \\), and \\( j = 1, 2 \\) are \\( M = \pm\tfrac12 \\).

> **Note on manual Sec. 60.2**: the text gives the orbital numbers as "1 to \\( 4l+1 \\)", but both the table in that section and the existing implementation use
> 1 to \\( 2(2l+1) \\) (\\( j = 0 \\) is also accepted).

**Matrix form**: the target state is represented by a unit vector \\( \mathbf c \\) including spin (\\( 2(2l+1) \\)-dimensional in real orbital \\( \mu \\) and spin), and
\\( D = \Delta\\,\mathbf c\mathbf c^\dagger \\). For the collinear specification \\( \mathbf c \\) is the unit vector \\( e\_{(\mu,\sigma)} \\);
for the \\( j \\)-dependent specification it is \\( \Phi\_{JM} \\) above transformed to real spherical harmonics
(\\( Y\_{lm} = \sum\_\mu \overline{U^{(l)}\_{\mu m}}\\,S\_{l\mu} \\), \\( U^{(l)} \\) of the [PAO basis](pao-basis.md)); for \\( j = 0 \\), \\( D = \Delta\\,I \\).
The matrix elements between basis functions are

$$
H^{\mathrm{CH},\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) = f(d\_{IJ})\sum\_{\mathbf L'}\sum\_{\mu\nu}
S\_{I\alpha,\\,(C+\mathbf L')\mu}\\; D^{\sigma\sigma'}\_{\mu\nu}\\; S\_{(C+\mathbf L')\nu,\\,J\beta+\mathbf L}
$$

\\( S \\) is the overlap integral ([Two-center integrals](two-center-overlap-kinetic.md), including the on-site block with \\( I = C \\)),
\\( f \\) is the same damping function as for the nonlocal pseudopotential ([Nonlocal pseudopotential](nonlocal-pseudopotential.md) Sec. 3.4, \\( r\_{\mathrm{cut}} = r\_{c,I}+r\_{c,J} \\)),
and the range of atom pairs is also the same. The sum over \\( \mathbf L' \\) runs over **all periodic images** of atom \\( C \\) (there is one hole in each periodic image of the supercell).
This matrix is added to the Hamiltonian (in the collinear case to the up and down spin blocks; in the noncollinear case to all 2×2 components), and with spin-orbit coupling the imaginary part is added as well.

**Energy and forces**: `Uch` \\( = \sum\_{\sigma\sigma'}\sum\_{\mathrm{pairs}}\rho^{\sigma'\sigma}H^{\mathrm{CH},\sigma\sigma'} \\) (the same trace rule as `Unl` in [Total energy and forces](energy-force.md)) and
is included in `Utot`. Forces are added from the derivatives of \\( S \\) and of the damping function, in the same form as the nonlocal term.

**Prerequisites**: `scf.core.hole on` requires spin polarization (`on` or `nc`); `off` is an input error.
For shells other than `s`, manual Sec. 60.2 requires `nc` (with spin-orbit coupling).

### 5.3 Periodic images and splitting of the Hartree potential [Public] [Convention]

[Public] Technical note 24 and manual Sec. 60.3: the final state of a gapped system is treated as a charged state (not neutralized, because the escape time of the photoelectron is short), and
to remove the electrostatic interaction between holes in periodic images, the final-state density is split as \\( n\_f = n\_i + \Delta n \\) (\\( n\_i \\) is the periodic density of the initial state,
\\( \int\Delta n = -1 \\)), with

$$
V\_H = V\_H^{(\mathrm P)}[n\_i] + V\_H^{(\mathrm{NP})}[\Delta n],\qquad
\tilde V\_H^{(\mathrm{NP})}(\mathbf G) = \frac{4\pi}{G^2}\left(1 - \cos(GR\_c)\right)\Delta\tilde n(\mathbf G)
$$

\\( V\_H^{(\mathrm P)} \\) is the usual periodic solution and \\( V\_H^{(\mathrm{NP})} \\) is the solution with the spherical Coulomb cutoff (\\( R\_c \\) is half the length of the shortest lattice vector).
If \\( \Delta n \\) fits in a sphere of radius \\( R \\), the hole is electrostatically isolated from its periodic images when \\( R\_c \ge 2R \\) and \\( 4R < L \\).

[Convention]

- This splitting becomes effective when all three of `scf.coulomb.cutoff on`, `scf.restart on`, and **`scf.restart.filename` given explicitly** hold
  (regardless of the value of `scf.core.hole`). Then the density of the restart file (the density of the initial-state SCF) is kept as \\( n\_i \\), and at the start of the SCF
  \\( V\_H^{(\mathrm P)} \\) is solved once with periodic boundaries from \\( \delta n\_i = n\_i - n^{a} \\) (when `scf.core.hole on`).
  In each subsequent SCF iteration the cutoff solution is applied to \\( n - n\_i \\), giving \\( \delta V\_H = V\_H^{(\mathrm P)} + V\_H^{(\mathrm{NP})}[n - n\_i] \\).
  The density of the restart file also serves as the initial density of the SCF.
- `UH1` is \\( \tfrac12\Delta V\sum\_{\mathbf n}\delta n\\,\delta V\_H \\) (\\( \delta n = n - n^{a} \\), \\( \delta V\_H \\) the sum above).
  Since the two kernels differ, the cross terms of this do not agree with the symmetric splitting of the electrostatic energy (\\( \tfrac12 n\_i V^{(\mathrm P)}[n\_i] + \Delta n\\,V^{(\mathrm P)}[n\_i] + \tfrac12\Delta n\\,V^{(\mathrm{NP})}[\Delta n] \\)). This is recorded as the definition of the existing implementation.
- For gas-phase molecules (manual Sec. 60.2), no restart file is used, and both initial and final states apply the cutoff to the whole \\( \delta n \\) with `scf.coulomb.cutoff on` ([Grid](grid-poisson-xc.md) Sec. 3.3).
- If the three conditions for splitting hold but `scf.core.hole off`, \\( V\_H^{(\mathrm P)} \\) is not computed and stays 0 (the periodic part of the Hartree potential is missing). The same happens if reading the restart file fails.
- **Combination with ESM** (manual Sec. 60.4): if `ESM.switch` is on, `scf.core.hole on`, and the restart file was read,
  \\( V\_H^{(\mathrm P)} \\) is solved with ESM (Sec. 2), after which ESM is disabled for the rest of that calculation (including the wall force).
  Hence the final state must be combined with `scf.coulomb.cutoff on` (if off, the whole \\( \delta n \\) is subsequently solved with periodic boundaries).
  Since ESM is used only for the neutral initial-state density, `on1` is used.

### 5.4 ΔSCF procedure [Public] [Convention]

The procedure of manual Secs. 60.2–60.4 is organized from the standpoint of a reimplementation. The initial- and final-state calculations use the same atomic configuration, cell, grid (`scf.Ngrid` or the same cutoff),
k-points, electronic temperature, and spin treatment. In the final state, only the species of the hole atom is switched to the basis for the hole (`…_CH`).

| System | Initial state | Final state | Binding energy |
|---|---|---|---|
| Gas-phase molecule | `scf.coulomb.cutoff on`, spin polarization on/nc | Same as left + `scf.system.charge 1`, `scf.core.hole on`, `<core.hole.state` | \\( E\_b^{\mathrm{(gas)}} \\) |
| Gapped solid | Periodic boundaries, write restart file | `scf.restart on`, `scf.restart.filename` (initial state), `scf.coulomb.cutoff on`, `scf.core.hole on`, `scf.system.charge 1` | \\( E\_b^{\mathrm{(bulk)}} \\) (\\( \mu\_0 \\) of the initial state) |
| Metal | Same as above | Same as above with `scf.system.charge 0` | \\( E\_b^{\mathrm{(metal)}} \\) |
| Surface / 2D system | `ESM.switch on1`, write restart file | Final state of the gapped solid + `ESM.switch on1` | \\( E\_b^{\mathrm{(bulk)}} \\) |

Numerical examples from the public documents (Hartree→eV conversion checked by the script): C-1s of C₂H₂ is \\( -66.084858926233 - (-76.787732114928) \\) Ha \\( = 291.24 \\) eV (experiment 291.14 eV),
C-1s of TiC (216 atoms, metal) is 281.55 eV, and 2p₃/₂ and 2p₁/₂ of Si (\\( \mu\_0 = -0.201641583 \\) Ha) are 99.41 and 100.13 eV.
As a guide to convergence, a hole-hole distance of 15–27 Å (about 500 atoms) is needed for gapped systems and about 64 atoms for metals (technical note 24).
With ESM, N-1s of 2D NBe₂ converges at a vacuum length of 10 Å (Table 14 of manual Sec. 60.4).

### 5.5 Input keywords [Public] [Convention]

| Keyword | Default | Meaning |
|---|---|---|
| `scf.core.hole` | off | Add the penalty term |
| `<core.hole.state` … `core.hole.state>` | — | One line: atom number, `s`/`p`/`d`/`f`, state number \\( j \\) (Sec. 5.2) |
| `scf.base.penalty.value` | 20.0 (Hartree) | Reference value of \\( \Delta \\) [Convention] |
| `scf.coulomb.cutoff` | off | Spherical Coulomb cutoff |
| `scf.restart`, `scf.restart.filename` | off, `System.Name` | Reading of the initial-state density (condition for the splitting of Sec. 5.3) |
| `scf.system.charge` | 0 | Charge of the hole in the final state (1 for gapped systems and gases, 0 for metals) |

---

## 6. Requirements for implementations [Design]

- **Pluggable Poisson solver**: abstract the Poisson equation as a function "\\( \delta n \\) (on the grid) → \\( \delta V\_H \\) (on the grid)",
  so that periodic boundaries, spherical cutoff, ESM (`on1`–`on4`), and difference mode (Sec. 5.3) can be selected as implementations of the same type.
  The ESM implementation takes the 1D Green's function (Sec. 2.3) as an argument, and edge smoothing (the "Edge smoothing" item of Sec. 2.4) is separated as an option for comparison with the existing implementation.
- **Registration of external potentials**: external terms such as the uniform field, the potential difference of ESM `on4`, and the wall are registered as modules returning three things: "contribution to the grid potential", "contribution to the energy", and "forces on atoms",
  with switching between compatibility mode and a physically consistent mode (Sec. 3.4, `on4` in Sec. 2.5, the wall energy in Sec. 2.6).
- **Penalty term**: the core-hole \\( H^{\mathrm{CH}} \\) has the same three-center structure as the nonlocal pseudopotential ("projector = function on some atom"), so
  implement it within the same framework of atom-pair range, damping function, and force calculation, only supplying "the basis functions of the specified atom" as projectors.
  \\( \mathbf c \\) (the target state) treats the collinear specification and the \\( (J,M) \\) specification as the same \\( 2(2l+1) \\)-dimensional vector.
- **Virtual atoms**: keep species data (\\( Z\_I \\), \\( V^{\mathrm{loc}} \\), \\( n^{a} \\)) as reals, and do not create processing that depends on integer atomic numbers.
- **Restart files**: for the initial-state density used in difference mode, check that the grid (number of divisions, origin) matches.

## 7. Tests [Design]

`tools/check_esm_field_corehole.py` contains reference implementations of the analytic items below.

| Kind | Content | Tolerance (guide) |
|---|---|---|
| Analytic | Field unit conversion: 1 GV/m \\( = 1.944690\times10^{-3} \\) Ha/(\\( e \\) bohr). Difference from the coefficient of the existing implementation | Relative \\( 2\times10^{-6} \\) |
| Analytic | ESM Green's functions (`on1`–`on3`, \\( G\_\parallel = 0, 0.25, 1.3 \\) bohr⁻¹): with trigonometric interpolation on a 24 bohr interval with 192 points, the "periodic solution + homogeneous solution" of Sec. 2.4 agrees with direct numerical integration. The closed form for `on1` (Sec. 2.4) also agrees | \\( 10^{-6} \\) |
| Analytic | Symmetry of the Green's functions, zero at metal planes, field \\( -4\pi q \\) on the metal side for `on3` at \\( G\_\parallel=0 \\), \\( V\_{\mathrm{bias}}(-z\_1) = \Delta\Phi \\) and \\( V\_{\mathrm{bias}}(z\_1) = 0 \\) for `on4` | \\( 10^{-12} \\) (field by finite difference \\( 10^{-3} \\)) |
| Analytic | Sawtooth: jump \\( E a \\) at the cell face; \\( U\_{\mathrm{ef}} = \mathbf E\cdot(\mathbf x\_0 - \mathbf c) \\) for a small lump of unit charge | \\( 10^{-10} \\) |
| Analytic | Core-hole projection: \\( \sum\_{J,M}\mathbf c\mathbf c^\dagger = I \\) (\\( l = 0\dots3 \\)). Examples of diagonal components: p \\( (J,M) = (\tfrac32,\tfrac32) \\) has \\( \tfrac12 \\) each on \\( p\_x\uparrow, p\_y\uparrow \\). \\( (\tfrac12,-\tfrac12) \\) has \\( \tfrac13 \\) each on \\( p\_x\uparrow, p\_y\uparrow, p\_z\downarrow \\). d \\( (\tfrac52,\tfrac12) \\) has \\( \tfrac35 \\) on \\( d\_{3z^2-r^2}\uparrow \\) and \\( \tfrac15 \\) each on \\( d\_{xz}\downarrow, d\_{yz}\downarrow \\). \\( (\tfrac52,\tfrac32) \\) has \\( \tfrac25 \\) on \\( d\_{xz}\uparrow, d\_{yz}\uparrow \\) and \\( \tfrac1{10} \\) on \\( d\_{x^2-y^2}\downarrow, d\_{xy}\downarrow \\) | \\( 10^{-10} \\) |
| Analytic | Arithmetic of binding energies (numbers from the public documents in Sec. 5.4) | 0.005 eV |
| Invariant | Charge doping: \\( \Delta V\sum\delta n = -Q \\) (within grid discretization error); with `scf.coulomb.cutoff`, `Utot` of a charged isolated system converges with cell size, while with periodic boundaries it varies as \\( L^{-1} \\) | Trend check |
| Invariant | ESM `on1`: for a neutral symmetric slab, the \\( \mathbf G\_\parallel=0 \\) potential is equal at both ends. For a slab with a dipole, the difference between the ends is \\( 4\pi \\) × (dipole per unit area) | \\( 10^{-6} \\) Ha |
| Invariant | ESM `on1`–`on3`: agreement between forces and finite differences of the total energy (atoms outside the buffer region) | About \\( 10^{-4} \\) Ha/bohr |
| Invariant | Electric field (physical mode of Sec. 3.4): total force 0 for neutral molecules, translational invariance of the total energy, agreement of forces with finite differences. In compatibility mode, total force \\( \approx -N\_{\mathrm e}\mathbf E \\) | \\( 10^{-4} \\) |
| Invariant | Core hole: in the final state the occupation of the target state \\( \mathrm{Tr}(\rho\\,\cdot) \\) is close to 0, and `Uch` \\( \ge 0 \\) is small. Doubling \\( \Delta \\) changes \\( E\_b \\) by at most 0.01 eV | 0.01 eV |
| Comparison | Build fixtures from the existing implementation's examples: `Al-Si111_ESM.dat` (`on3`, \\( Q = 0.005 \\)) and its `on1`, `on2`, `on4` variants, `Nitro_Benzene.dat` (field 10 GV/m), `DIA8-VA.dat` (virtual atom), `C2H2.dat`/`C2H2-CH.dat`, `N1Be2_IS.dat`/`N1Be2_FS.dat` (ESM + core hole) | Energy \\( 10^{-6} \\) Ha, binding energy 0.01 eV (to be decided) |

## [Open]

- The magnitude of the effect of ESM edge smoothing (the "Edge smoothing" item of Sec. 2.4) on energies and forces, and whether to make it the default in a reimplementation.
- The physical meaning of the energy definition for `on4` (`UH1` contains \\( \tfrac12\int\delta n\\,V\_{\mathrm{bias}} \\)), and the treatment of the correct free energy (the work term charge × potential) under constant electrode potential.
- The treatment of the uniform field in non-orthogonal cells (the formula of Sec. 3.2 is not a uniform field except for orthogonal cells). Confirmation of the intent of the existing implementation (treating the input components as along the lattice vectors).
- The meaning of combining the uniform field with ESM, and the field with charged systems.
- Whether `<Number.of.Electrons.Species` is described in the public documents, and its intended use.
- For the \\( j \\)-dependent core-hole specification, that the phases of the transformation to real spherical harmonics agree exactly with the existing implementation (diagonal components and \\( \sum\mathbf c\mathbf c^\dagger = I \\) are confirmed; phases of off-diagonal components are to be confirmed with fixtures).
- The size of the error in binding energies caused by `UH1` in the difference mode of Sec. 5.3 (asymmetric cross terms).
- The treatment of stress in ESM calculations (not supported by the existing implementation).
