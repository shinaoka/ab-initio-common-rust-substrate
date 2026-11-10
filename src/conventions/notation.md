# Notation, units and terminology

These conventions are shared by all method families. A method chapter may add symbols of its own,
but must not redefine the ones below.

## Units [Convention]

- Internal units are Hartree atomic units: lengths in bohr, energies in Hartree, wave numbers in bohr⁻¹.
- Input files may use Å, eV or Ry; they are converted to internal units when read.
- When a cutoff energy is given in Ry, the corresponding wave number is \\( k = \sqrt{E\_{\mathrm{Ry}}} \\) bohr⁻¹
  (\\( E\_{\mathrm{Ha}} = k^2/2 \\), \\( E\_{\mathrm{Ry}} = k^2 \\)).

## Symbols [Design]

| Symbol | Meaning |
|---|---|
| \\( I, J \\) | Atoms (atom index within the unit cell) |
| \\( \mathbf{R}\_I \\) | Position of atom \\( I \\) |
| \\( \mathbf{L} \\) | Lattice vector (periodic image) |
| \\( \alpha = (l, \zeta, m) \\) | Orbital on an atom. \\( l \\) is the angular momentum, \\( \zeta \\) the multiplicity (0-based), \\( m \\) the magnetic quantum number |
| \\( R\_{l\zeta}(r) \\) | Radial function |
| \\( Y\_{lm} \\) | Complex spherical harmonics (Condon–Shortley phase) |
| \\( S\_{lm} \\) | Real spherical harmonics (ordering: see [PAO basis](../methods/lcao-openmx/modules/pao-basis.md)) |
| \\( r\_c \\) | Cutoff radius of a PAO |
| \\( \boldsymbol\tau \\) | Interatomic vector in a two-center integral |

## Fourier transform convention [Design]

This book uses

$$
\tilde f(\mathbf k) = \int d^3 r\\, e^{-i\mathbf k\cdot\mathbf r} f(\mathbf r),\qquad
f(\mathbf r) = \frac{1}{(2\pi)^3}\int d^3 k\\, e^{i\mathbf k\cdot\mathbf r} \tilde f(\mathbf k)
$$

OpenMX technical note 1 attaches \\( (2\pi)^{-3} \\) to the forward transform, so take care when
comparing coefficients. The radial spherical Bessel transform is

$$
\tilde R\_l(k) = \int\_0^\infty dr\\, r^2 j\_l(kr) R\_l(r),\qquad
R\_l(r) = \frac{2}{\pi}\int\_0^\infty dk\\, k^2 j\_l(kr) \tilde R\_l(k)
$$
