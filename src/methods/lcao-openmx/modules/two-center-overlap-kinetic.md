# Two-center integrals: overlap and kinetic energy

> Status: prototype (awaiting review). The formulas re-derive the method of technical note 1 (Eqs. 36–46) using this note's
> [Fourier transform convention](../../../conventions/notation.md).

## Definition

Let \\( \phi\_a \\) be a basis function of atom \\( I \\), \\( \phi\_b \\) a basis function of atom \\( J \\),
and the interatomic vector \\( \boldsymbol\tau = \mathbf R\_J + \mathbf L - \mathbf R\_I \\).

$$
S\_{ab}(\boldsymbol\tau) = \int d^3 r\\, \phi\_a^{\ast}(\mathbf r)\\, \phi\_b(\mathbf r - \boldsymbol\tau),\qquad
T\_{ab}(\boldsymbol\tau) = \int d^3 r\\, \phi\_a^{\ast}(\mathbf r) \left(-\tfrac12 \nabla^2\right) \phi\_b(\mathbf r - \boldsymbol\tau)
$$

## Method [Public]

As in technical note 1, the computation is done in wave-number space.

### 1. Transformation of radial functions

Writing \\( \phi\_a(\mathbf r) = R\_a(r) Y\_{lm}(\hat{\mathbf r}) \\) with complex spherical harmonics, the plane-wave expansion

$$
e^{-i\mathbf k\cdot\mathbf r} = 4\pi \sum\_{LM} (-i)^L j\_L(kr)\\, Y\_{LM}(\hat{\mathbf k})\\, Y^{\ast}\_{LM}(\hat{\mathbf r})
$$

and the orthonormality of spherical harmonics give

$$
\tilde\phi\_a(\mathbf k) = 4\pi (-i)^{l}\\, Y\_{lm}(\hat{\mathbf k})\\, \tilde R\_a(k),\qquad
\tilde R\_a(k) = \int\_0^\infty dr\\, r^2 j\_l(kr) R\_a(r)
$$

### 2. Overlap

Using Parseval's identity \\( \int f^{\ast} g = (2\pi)^{-3}\int \tilde f^{\ast}\tilde g \\) and the translation
\\( \widetilde{\phi\_b(\cdot - \boldsymbol\tau)}(\mathbf k) = e^{-i\mathbf k\cdot\boldsymbol\tau}\tilde\phi\_b(\mathbf k) \\),
and expanding \\( e^{-i\mathbf k\cdot\boldsymbol\tau} \\) again in plane waves, one obtains

$$
S\_{ab}(\boldsymbol\tau) = 8 \sum\_{LM} i^{\\,l - l' - L}\\; G(lm;\\, l'm',\\, LM)\\; Y^{\ast}\_{LM}(\hat{\boldsymbol\tau})\\; I^{(0)}\_{L}(\tau)
$$

$$
I^{(n)}\_{L}(\tau) = \int\_0^\infty dk\\, k^{2+n}\\, j\_L(k\tau)\\, \tilde R\_a(k)\\, \tilde R\_b(k)
$$

Here \\( l' m' \\) are the quantum numbers of \\( \phi\_b \\), \\( \tau = \lVert\boldsymbol\tau\rVert \\), and

$$
G(lm;\\, l'm',\\, LM) = \int d\Omega\\; Y^{\ast}\_{lm}\\, Y\_{l'm'}\\, Y\_{LM}
$$

is the Gaunt coefficient. The factor 8 is \\( (4\pi)^3/(2\pi)^3 \\).

- The Gaunt coefficient is nonzero only when \\( l + l' + L \\) is even and the triangle condition \\( |l - l'| \le L \le l + l' \\) holds,
  so \\( i^{\\,l-l'-L} \\) is \\( \pm 1 \\). The sum is finite.
- Terms other than \\( M = m - m' \\) vanish.

### 3. Kinetic energy

Since \\( -\tfrac12\nabla^2 \\) corresponds to multiplication by \\( k^2/2 \\) in wave-number space,

$$
T\_{ab}(\boldsymbol\tau) = 4 \sum\_{LM} i^{\\,l - l' - L}\\; G(lm;\\, l'm',\\, LM)\\; Y^{\ast}\_{LM}(\hat{\boldsymbol\tau})\\; I^{(2)}\_{L}(\tau)
$$

> **Check**: the overlap formula (several pairs with \\( l, l' \le 2 \\)) and the kinetic-energy formula (s–s)
> were confirmed with `tools/check_two_center.py` to agree to 6 or more digits with direct integration on a 3D grid and with analytic solutions,
> for Gaussian radial functions \\( r^l e^{-\alpha r^2} \\).

> **Note on Eq. 46 of technical note 1**: in that equation, the argument of the spherical Bessel function is printed as \\( kr \\), but
> from the correspondence with the overlap formula (Eq. 44), \\( k\lvert\boldsymbol\tau\rvert \\) is believed to be correct.
> The derivation above is consistent with that.

### 4. Transformation to real spherical harmonics

In the real basis \\( \phi^{R}\_{\mu} = \sum\_m U^{(l)}\_{\mu m}\phi\_{m} \\) ([PAO basis](pao-basis.md)),

$$
S^{R}\_{\mu\nu} = \sum\_{m m'} \overline{U^{(l)}\_{\mu m}}\\; U^{(l')}\_{\nu m'}\\; S\_{m m'}
$$

The result is real (the imaginary part is rounding error only). One may also compute directly in real arithmetic using real Gaunt coefficients and real spherical harmonics
\\( S\_{LM}(\hat{\boldsymbol\tau}) \\).

[Convention] The existing implementation computes in the complex basis, then applies the above transformation and keeps only the real part.

## Range of atom pairs [Convention]

The atom pairs for which the overlap and kinetic energy are computed are all \\( (I, J, \mathbf L) \\) satisfying

$$
\lVert \boldsymbol\tau \rVert \le r\_{c,I} + r\_{c,J}
$$

- \\( r\_c \\) is the value of `radial.cutoff.pao` in the PAO file of each species ([PAO basis](pao-basis.md)).
- The inequality includes equality. Pairs at a distance of exactly \\( r\_{c,I} + r\_{c,J} \\) are included, but their values are 0.
- The atom itself (\\( J = I \\), \\( \mathbf L = 0 \\)) is included as a pair.
- Periodic images (\\( J = I \\), \\( \mathbf L \ne 0 \\)) are also treated as separate pairs if they satisfy the condition.
- For pairs outside the range, \\( S = T = 0 \\) and they are not stored in the matrix.

[Design] The draft of the [sparse block format](../design/sparse-block-format.md) uses the strict inequality \\( < \\).
Since the values of pairs on the boundary are 0, physical quantities do not change whichever is adopted.
Matching \\( \le \\) is necessary only when comparing even the ordering and number of atom pairs with the existing implementation.

## Numerical parameters of the wave-number integration [Convention]

### Input keywords [Public]

Manual Section 8.2 has the following three keywords.

| Keyword | Default | Meaning (summary of the manual's description) |
|---|---|---|
| `1DFFT.EnergyCutoff` | 3600 (Ry) | Upper limit of the wave number used in the Fourier transform of radial functions, given as an energy |
| `1DFFT.NumGridK` | 900 | Number of points of the wave-number grid on which the transformed radial functions are stored |
| `1DFFT.NumGridR` | 900 | Number of points of the real-space radial grid used for the numerical integration of the transform |

The upper limit of the wave number is \\( k\_{\max} = \sqrt{E\_{\mathrm{cut}}[\mathrm{Ry}]} \\) ([notation](../../../conventions/notation.md)). With the default it is 60 bohr⁻¹.

### Usage in the existing implementation

Below, \\( N\_K \\) denotes `1DFFT.NumGridK` and \\( N\_R \\) denotes `1DFFT.NumGridR`.

**(a) Spherical Bessel transform of radial functions** \\( \tilde R(k) = \int dr\\, r^2 j\_l(kr) R(r) \\)

| Item | Value |
|---|---|
| Integration range | \\( r\_{\min} \le r \le r\_c + 0.5 \\) bohr (\\( r\_{\min} \\) is the first grid point of the PAO file) |
| Quadrature | \\( N\_R \\) equally spaced intervals in \\( r \\), trapezoidal rule. \\( R(r) \\) at each point is interpolated from the table |
| Prefactor | None (as defined in [notation](../../../conventions/notation.md); no \\( \sqrt{2/\pi} \\) etc.) |
| Stored wave numbers | \\( k\_j = j\\, k\_{\max}/N\_K \\), \\( j = 0, 1, \dots, N\_K - 1 \\) |
| Outside storage | \\( \tilde R = 0 \\) for \\( k > k\_{N\_K-1} \\). Interpolated in between |
| Targets | Only the radial functions used by the basis specification (the first \\( N\_l \\) of each \\( l \\)) |

**(b) Wave-number integration for overlap and kinetic energy** \\( I^{(n)}\_L(\tau) = \int dk\\, k^{2+n} j\_L(k\tau)\tilde R\_a\tilde R\_b \\)

| Item | Value |
|---|---|
| Integration range | \\( k\_{\min} \le k \le k\_{\max} \\), \\( k\_{\min} = 10^{-6} \\) bohr⁻¹ |
| Quadrature | \\( N\_R \\) equally spaced intervals in \\( k \\) (\\( N\_R + 1 \\) points), trapezoidal rule |
| \\( \tilde R(k) \\) | Interpolation of the table from (a) |
| Spherical Bessel functions | \\( L = 0, \dots, 2\max(l\_{\max,I}, l\_{\max,J}) \\) evaluated directly at each point |

- `1DFFT.NumGridR` is used both for the number of \\( r \\) intervals in (a) and for the number of \\( k \\) intervals in (b).
  With the default, the step in (b) is \\( 60/900 \approx 0.067 \\) bohr⁻¹.
- Gauss–Legendre quadrature is not used. The "128-point Gauss–Legendre" in technical note 15 concerns other integrals (projector expansion, etc.).
- Since \\( \tilde R \\) near \\( k\_{\max} \\) is sufficiently small, the effect of the storage range of (a) falling one step short of \\( k\_{\max} \\) is negligible.

[Design] The reimplementation may freely choose quadrature and grids as long as the same accuracy is achieved
(e.g. the Toyoda–Ozaki fast spherical Bessel transform (CPC **181**, 277), Gaussian quadrature rules).
However, the upper limit \\( k\_{\max} \\) is determined from `1DFFT.EnergyCutoff` by the formula above, and the default of 3600 Ry is kept.
`1DFFT.NumGridK` and `1DFFT.NumGridR` are accepted as input and interpreted as specifying the resolution of the integration.
\\( \tilde R(k) \\) is tabulated only once per species and not recomputed per atom pair.

**Accuracy estimate**: for actual PAOs (the first s function of C6.0, H5.0, O6.0), comparing the s–s overlap on the default grid above with
direct integration in real space gave differences of at most \\( 4\times 10^{-7} \\) for \\( 0 \le \tau \le 12 \\) bohr.
The residual of the wave-number integration at \\( \tau \ge r\_{c,a}+r\_{c,b} \\) is about \\( 10^{-8} \\) (`tools/check_pao_file.py`).

## Special cases

### \\( \tau = 0 \\) [Convention]

Since \\( j\_L(0) = \delta\_{L0} \\), only the \\( L = 0 \\) term remains, and

$$
S\_{ab}(0) = \delta\_{ll'}\delta\_{mm'}\\,\frac{2}{\pi}\int\_0^\infty dk\\, k^2 \tilde R\_a \tilde R\_b = \delta\_{ll'}\delta\_{mm'}\int\_0^\infty dr\\, r^2 R\_a R\_b,
\qquad
T\_{ab}(0) = \delta\_{ll'}\delta\_{mm'}\\,\frac{1}{\pi}\int\_0^\infty dk\\, k^4 \tilde R\_a \tilde R\_b
$$

The direction \\( \hat{\boldsymbol\tau} \\) is undefined, but the \\( L = 0 \\) term does not depend on direction.

- The existing implementation replaces \\( \tau < 10^{-10} \\) bohr by \\( \tau = 10^{-10} \\) bohr and computes with the general formula.
  The \\( L \ge 1 \\) terms then effectively vanish, in agreement with the above formula. The direction may be arbitrary.
- For the pair of an atom with itself (\\( J = I \\), \\( \mathbf L = 0 \\)), the derivative with respect to the interatomic vector (below) is set to 0.
- [Design] The reimplementation treats \\( \tau < 10^{-10} \\) bohr as \\( \tau = 0 \\) and uses the closed formula above.
  The overlap between different \\( \zeta \\) of the same \\( l \\) on the same atom is close to 0 owing to the orthogonality of the table ([PAO basis](pao-basis.md)).

### \\( \tau \ge r\_{c,a} + r\_{c,b} \\)

Exactly 0. Numerical integration in wave-number space leaves a small residual (above), so atom pairs in this range are not computed.

### Small random numbers added to the overlap matrix [Convention]

The existing implementation adds to each element of the overlap matrix a uniform random number of width \\( 10^{-13} \\) (\\( [-5\times10^{-14}, 5\times10^{-14}] \\)).
This is presumably intended to break numerical degeneracies.

[Design] The reimplementation does not add them. Comparison tolerances are taken larger than \\( 10^{-12} \\).

## Derivatives with respect to the interatomic vector (for force calculations)

[Public] The gradient of the product of \\( f\_L(\tau) = I^{(n)}\_L(\tau) \\) and \\( Y^{\ast}\_{LM}(\hat{\boldsymbol\tau}) \\) is

$$
\nabla\_{\boldsymbol\tau}\left[f\_L(\tau)\\,Y^{\ast}\_{LM}(\hat{\boldsymbol\tau})\right]
{}= f\_L'(\tau)\\,Y^{\ast}\_{LM}\\,\hat{\boldsymbol\tau}
{}+ \frac{f\_L(\tau)}{\tau}\left(\partial\_\theta Y^{\ast}\_{LM}\\,\mathbf e\_\theta + \frac{1}{\sin\theta}\partial\_\varphi Y^{\ast}\_{LM}\\,\mathbf e\_\varphi\right)
$$

with \\( f\_L'(\tau) = \int dk\\, k^{3+n} j\_L'(k\tau)\tilde R\_a\tilde R\_b \\). Technical note 1 states that
these derivatives are evaluated analytically in force calculations.

[Convention]

- The existing implementation computes and stores, together with \\( S \\) and \\( T \\), the three components of the derivative with respect to the position of the central atom,
  \\( \partial S\_{ab}/\partial \mathbf R\_I = -\nabla\_{\boldsymbol\tau} S\_{ab} \\) (likewise for \\( T \\)), on the same wave-number grid.
- At the poles (about \\( \lvert\sin\theta\rvert < 10^{-13} \\); pairs with \\( \boldsymbol\tau \\) parallel to the \\( z \\) axis),
  the azimuth is set to 0 and the \\( \frac{1}{\sin\theta}\partial\_\varphi \\) term is dropped.
  Since \\( Y\_{LM} \propto \sin^{\lvert M\rvert}\theta \\), the limiting value of this term is nonzero only for \\( \lvert M\rvert = 1 \\),
  and that contribution (the \\( y \\) component at azimuth 0) is missing. For example, placing a pair of an s orbital and a \\( p\_y \\) orbital on the \\( z \\) axis,
  the correct \\( \partial S/\partial \tau\_y = f(\tau)/\tau \\) (\\( f \\) is the radial part) is nonzero, but this treatment gives 0.
- [Design] The reimplementation either includes the limiting value correctly, or writes the derivatives in Cartesian components to avoid the singularity.
  Therefore, when comparing derivatives with the existing implementation, either avoid structures containing atom pairs parallel to the \\( z \\) axis, or allow for this difference.
  [Open] How much this difference shows up in the forces of the existing implementation (it is often canceled by symmetry) is checked by finite differences
  in the [energy and forces](energy-force.md) chapter.
- Matrix elements of the orbital angular momentum operator are also computed in the same routine, but they are outside the scope of this note.

## Tests

| Kind | Content | Approximate tolerance |
|---|---|---|
| Analytic | Two normalized s-type Gaussians \\( e^{-\alpha r^2} \\), \\( e^{-\beta r^2} \\): \\( S = \left(\frac{2\sqrt{\alpha\beta}}{\alpha+\beta}\right)^{3/2} e^{-\mu\tau^2} \\), \\( \mu = \frac{\alpha\beta}{\alpha+\beta} \\) | \\( 10^{-10} \\) |
| Analytic | Kinetic energy of the same pair: \\( T = \mu\\,(3 - 2\mu\tau^2)\\, S \\) | \\( 10^{-10} \\) |
| Analytic | Overlap of p- and d-type Gaussians (closed form by the Gaussian product theorem) | \\( 10^{-10} \\) |
| Analytic | \\( S\_{ab} = \delta\_{ll'}\delta\_{mm'}\int r^2 R\_a R\_b\\,dr \\) at \\( \tau = 0 \\) | \\( 10^{-12} \\) |
| Analytic | Pairs with \\( \tau < 10^{-10} \\) bohr agree with the \\( \tau = 0 \\) formula and are independent of direction | \\( 10^{-12} \\) |
| Invariant | \\( S\_{ab}(\boldsymbol\tau) = \overline{S\_{ba}(-\boldsymbol\tau)} \\), likewise for \\( T \\) | Rounding error |
| Invariant | Under rotation of the whole structure, the matrices transform by the rotation matrices of real spherical harmonics | \\( 10^{-10} \\) |
| Invariant | Pairs with \\( \tau \ge r\_{c,a}+r\_{c,b} \\) are not in the list. Pairs with exactly \\( \tau = r\_{c,a}+r\_{c,b} \\) are in the list in the existing implementation | Exact |
| Invariant | Agreement of analytic derivatives with finite differences (including \\( \boldsymbol\tau \\) near the poles) | About \\( 10^{-6} \\) |
| Convergence | For actual PAOs, difference between the default wave-number grid and a finer grid (or direct integration in real space) | \\( 10^{-6} \\) |
| Comparison | For actual PAOs (e.g. C6.0, O6.0, Fe6.0S), compare with the output of the existing implementation (overlap matrix in scfout, etc.) at several \\( \boldsymbol\tau \\) | \\( 10^{-6} \\) (allowing for differences in quadrature; the random contribution \\( 10^{-13} \\) is negligible) |

> **Check**: comparison with analytic solutions is done with `tools/check_two_center.py` (above), and the accuracy of the default grid for actual PAOs with
> `tools/check_pao_file.py` (given the path of a PAO file as an argument).
