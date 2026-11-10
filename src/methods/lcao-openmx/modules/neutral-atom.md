# Neutral-atom potential and projector expansion

> Status: prototype (awaiting review). The formulas rederive, in the notation of this note, the methods of technical note 15 (corresponding to Ozaki–Kino, PRB **72**, 045121 (2005); hereafter TN15) and
> technical note 1 (Ver. 1.3; hereafter TN1).
> The parts checked numerically are in `tools/check_neutral_atom.py`.

## Notation [Design]

In addition to the symbols of [Notation](../../../conventions/notation.md) and [Total energy and forces](energy-force.md), the following are used.

| Symbol | Meaning |
|---|---|
| \\( Z\_I \\) | Valence charge of the species of atom \\( I \\) (the core charge of the pseudopotential; \\( Z\_v \\) of Section 2.2 of [Nonlocal pseudopotential](nonlocal-pseudopotential.md)) |
| \\( V^{\mathrm{loc}}\_I(r) \\) | Local pseudopotential. \\( -Z\_I/r \\) at large distance |
| \\( n^{a}\_I(r) \\) | Valence density of the confined isolated atom. \\( 4\pi\int\_0^\infty r^2 n^{a}\_I\\, dr = Z\_I \\) |
| \\( V^{a}\_{H,I}(r) \\) | Hartree potential of \\( n^{a}\_I \\) |
| \\( V^{\mathrm{na}}\_I(r) \\) | Neutral-atom potential |
| \\( r\_{c,I} \\) | PAO cutoff radius (`radial.cutoff.pao` in the PAO file) |
| \\( \boldsymbol\tau\_{IK} = \mathbf R\_K - \mathbf R\_I \\) | Interatomic vector (lattice vectors of periodic images are included) |
| \\( \bar R\_{K,l\zeta}(r) \\) | Radial function of the projector expansion of atom \\( K \\) (Section 3) |
| \\( c\_{K,l\zeta} \\) | Its \\( V^{\mathrm{na}}\_K \\)-weighted norm (Section 3) |
| \\( \beta^{K}\_{l\zeta m}(\mathbf r) = V^{\mathrm{na}}\_K(r)\\, \bar R\_{K,l\zeta}(r)\\, S\_{lm}(\hat{\mathbf r}) \\) | Projector (with atom \\( K \\) as origin) |
| \\( L\_{\max} \\), \\( N\_{\mathrm{rad}} \\) | Maximum angular momentum of the projector expansion, and the number of radial functions per \\( l \\) |

## 1. Neutral-atom potential

### 1.1 Definition and short range [Public]

Following TN1 Eq. 18 and TN15 (after Eq. 9), the neutral-atom potential of atom \\( I \\) is defined as

$$
V^{\mathrm{na}}\_I(r) = V^{\mathrm{loc}}\_I(r) + V^{a}\_{H,I}(r),\qquad
V^{a}\_{H,I}(r) = \int d^3 r'\\, \frac{n^{a}\_I(r')}{\lvert\mathbf r - \mathbf r'\rvert}
$$

\\( n^{a}\_I \\) is the atomic valence density obtained under the same confinement potential as used to generate the PAOs, and
vanishes for \\( r > r\_{c,I} \\). Since the Hartree potential of a spherically symmetric charge distribution is

$$
V^{a}\_{H,I}(r) = \frac{4\pi}{r}\int\_0^r dr'\\, r'^2 n^{a}\_I(r') + 4\pi\int\_r^\infty dr'\\, r' n^{a}\_I(r')
$$

we have \\( V^{a}\_{H,I}(r) = Z\_I/r \\) for \\( r \ge r\_{c,I} \\) (Gauss's law).
On the other hand, the local pseudopotential coincides with \\( -Z\_I/r \\) outside the generation radius of the pseudopotential.
Hence \\( V^{\mathrm{na}}\_I \\) vanishes for \\( r \ge r\_{c,I} \\), and is a short-range potential with cutoff radius \\( r\_{c,I} \\).
Physically, this corresponds to the fact that a neutral atom (ionic core plus confined valence electron cloud) produces zero potential outside.

By this property, the matrix elements of \\( \sum\_I V^{\mathrm{na}}\_I \\) are determined by a finite number of atom combinations only,
and the long-range electrostatic interaction remains only in the Hartree term of the difference density \\( \delta n = n - \sum\_I n^{a}\_I \\) (Section 6).

### 1.2 Input data [Convention]

- \\( n^{a}\_I \\): the `<valence.charge.density` block of the PAO file. Each line holds \\( \ln r \\), \\( r \\), \\( n^{a}(r) \\) (bohr⁻³)
  (the grid is the same as for the radial-function tables; [PAO basis](pao-basis.md)).
- The density read in is scaled by a constant so that \\( 4\pi\int r^2 n^{a}\\, dr = Z\_I \\).
  Some distributed files have tables whose integral does not match \\( Z\_I \\) (example: about 4.5 for C6.0, with \\( Z = 4 \\)).
  Therefore the renormalization cannot be omitted.
- \\( V^{\mathrm{loc}}\_I \\): the third column of the `<Pseudo.Potentials` block of the VPS file (Section 2.4 of [Nonlocal pseudopotential](nonlocal-pseudopotential.md)).
- \\( r\_{c,I} \\): `radial.cutoff.pao` in the PAO file. Not the number in the basis name (e.g. `C6.0`).
- For empty atoms (species with `AtomSpecies` equal to 0), \\( V^{\mathrm{na}} = 0 \\) and \\( n^{a} = 0 \\).

### 1.3 Construction of the radial table [Convention]

\\( V^{\mathrm{na}}\_I \\) is held as a table on the radial grid \\( \\{r\_i\\} \\) of the VPS file. The procedure is as follows.

1. From the renormalized \\( n^{a}\_I \\), compute \\( V^{a}\_{H,I}(r\_i) \\) by numerical integration of the formula in Section 1.1
   (\\( Z\_I/r \\) outside the density table).
2. **Asymptotic correction**: at the grid point \\( r\_{\ast} \\) closest to \\( r\_{c,I} \\), compute

$$
s\_I = -\frac{V^{\mathrm{loc}}\_I(r\_{\ast})}{V^{a}\_{H,I}(r\_{\ast})}
$$

   and multiply all of \\( V^{a}\_{H,I} \\) by \\( s\_I \\) (no correction when \\( Z\_I \\) is 0).
   From the argument of Section 1.1, \\( s\_I \\) is 1 within numerical error (\\( \lvert s\_I - 1\rvert \lesssim 10^{-4} \\) for the distributed files of C, O, H, Fe).
   With this correction, \\( V^{\mathrm{loc}}\_I + V^{a}\_{H,I} \\) becomes 0 numerically as well at \\( r\_{c,I} \\).
3. **Smooth truncation**:

$$
V^{\mathrm{na}}\_I(r) =
\begin{cases}
\dfrac{V^{\mathrm{loc}}\_I(r) + s\_I V^{a}\_{H,I}(r)}{1 + e^{20\\,(r - r\_{c,I})}} & r \le r\_{c,I} \\\\
0 & r > r\_{c,I}
\end{cases}
$$

   (\\( r \\) in bohr). The Fermi-type factor is 1/2 at \\( r\_{c,I} \\), but since the bracket is nearly 0 owing to the correction in step 2,
   the discontinuity at \\( r\_{c,I} \\) is negligible (for the distributed files, \\( \lvert V^{\mathrm{loc}} + sV^{a}\_H\rvert \\) within 1 bohr inside \\( r\_{c,I} \\) is of order \\( 10^{-4} \\) Hartree or less).
   Values for \\( r > r\_{c,I} \\) are treated as exactly 0, independent of interpolation.

\\( V^{a}\_{H,I} \\) itself (after correction) is also used in the screened interatomic term of [Total energy and forces](energy-force.md).

[Design] The interpolation scheme is not specified. For \\( r < r\_1 \\) (inside the first grid point), \\( V^{\mathrm{na}} \\) is extended smoothly in the form \\( a r^2 + b \\)
(since it is regular at the origin). The existing implementation also uses an extension of this form [Convention].

### 1.4 Wavenumber-space tables [Convention]

The wavenumber-space integrals of Sections 4 and 3.5 use the spherical Bessel transforms of the radial functions of \\( V^{\mathrm{na}}\_I \\) and \\( \beta^{K} \\),

$$
\tilde V^{\mathrm{na}}\_I(k) = \int\_0^{\infty} dr\\, r^2 j\_0(kr)\\, V^{\mathrm{na}}\_I(r)
$$

(the convention of [Notation](../../../conventions/notation.md)). The existing implementation truncates the real-space integral at \\( r\_{c,I} + 0.5 \\) bohr, and
tabulates on 300 Gauss–Legendre points in the wavenumber interval \\( 10^{-6} \le k \le k\_{\max} \\) (\\( k\_{\max} = \sqrt{E} \\), where \\( E \\) is `1DFFT.EnergyCutoff`; 60 bohr⁻¹ at the default 3600 Ry).
The wavenumber integrals are also performed on the same 300 points.
This grid is separate from the uniform grid for overlap and kinetic energy in [Two-center integrals](two-center-overlap-kinetic.md).

## 2. Classification of matrix elements [Public] [Convention]

The neutral-atom term of the Hamiltonian is

$$
V^{\mathrm{na}}\_{I\alpha,J\beta}(\mathbf L) = \sum\_{K} \langle \phi\_{I\alpha} \vert V^{\mathrm{na}}\_K \vert \phi\_{J\beta}(\cdot - \boldsymbol\tau\_{IJ\mathbf L}) \rangle
$$

where the sum runs over all atoms \\( K \\) (including periodic images) whose \\( V^{\mathrm{na}}\_K \\) support overlaps both basis functions.
The computational method is split according to the relation of \\( K \\) to \\( I \\) and \\( J \\), as follows.

| Case | Example | Method | Source |
|---|---|---|---|
| Basis pair on the same atom (\\( J = I \\), \\( \mathbf L = 0 \\)) | \\( K = I \\) (one-center), \\( K \ne I \\) (two-center) | Directly in wavenumber space (Section 4) | TN15, [Convention] |
| Basis pair on different atoms, with \\( K = I \\) or \\( K = J \\) | two-center | Projector expansion (Section 3). Exact, since the basis functions lie in the space spanned by the radial functions of the expansion | TN15 Eq. 16, [Convention] |
| Basis pair on different atoms, with \\( K \ne I, J \\) | three-center | Projector expansion (approximate) | TN15 Eq. 13 |
| With `scf.ProExpn.VNA off` | all | Sum on the real-space grid (Section 5) | manual, [Convention] |

[Public] TN15 states that, for the two-center case with \\( K = I \\) or \\( K = J \\), the property that the expansion reproduces \\( V^{\mathrm{na}}\_K\vert\phi\rangle \\) exactly (TN15 Eq. 16) is used,
and that the case \\( I = J \ne K \\) is computed separately by treating the product of \\( V^{\mathrm{na}}\_K \\) and the radial functions in wavenumber space.

[Convention] The concrete split in the existing implementation (`scf.ProExpn.VNA on`, default) is as follows.

- **Blocks on the same atom** (\\( J = I \\), \\( \mathbf L = 0 \\)): for all \\( K \\), including \\( K = I \\), the sum of the direct computations of Section 4 is used.
  Values from the projector expansion are not used.
- **Other blocks** (\\( J \ne I \\) or \\( \mathbf L \ne 0 \\); including pairs of periodic images of the same atom):
  the projector-expansion sum over all \\( K \\), including \\( K = I \\) and \\( K = J \\), is taken and multiplied by the smooth truncation factor of Section 3.6.
  The terms with \\( K = I, J \\) are not computed by a separate formula (since the expansion is exact, the result is the same).

## 3. Projector expansion

### 3.1 Form of the expansion [Public]

Following TN15 Eq. 13–15 (TN1 Eq. 32–34), \\( V^{\mathrm{na}}\_K \\) of atom \\( K \\) is approximated by the separable operator

$$
\hat V^{\mathrm{na}}\_K = \sum\_{l=0}^{L\_{\max}} \sum\_{m=-l}^{l} \sum\_{\zeta=1}^{N\_{\mathrm{rad}}}
\frac{\lvert \beta^{K}\_{l\zeta m} \rangle \langle \beta^{K}\_{l\zeta m} \rvert}{c\_{K,l\zeta}},\qquad
\beta^{K}\_{l\zeta m}(\mathbf r) = V^{\mathrm{na}}\_K(r)\\, \bar R\_{K,l\zeta}(r)\\, S\_{lm}(\hat{\mathbf r})
$$

The radial functions \\( \bar R\_{l\zeta} \\) (\\( K \\) omitted below) are obtained from the seed functions \\( \chi\_{l\zeta} \\) (Section 3.2) by Gram–Schmidt orthogonalization
with respect to the \\( V^{\mathrm{na}}\_K \\)-weighted inner product

$$
(f, g)\_V = \int\_0^\infty dr\\, r^2 f(r)\\, V^{\mathrm{na}}\_K(r)\\, g(r)
$$

that is:

$$
\bar R\_{l\zeta} = \chi\_{l\zeta} - \sum\_{\eta < \zeta} \frac{(\bar R\_{l\eta}, \chi\_{l\zeta})\_V}{c\_{l\eta}}\\, \bar R\_{l\eta},\qquad
c\_{l\zeta} = (\bar R\_{l\zeta}, \bar R\_{l\zeta})\_V .
$$

**Property (TN15 Eq. 16)**: for a function \\( f(\mathbf r) = R(r) S\_{lm}(\hat{\mathbf r}) \\) (\\( l \le L\_{\max} \\)) whose radial function lies in \\( \mathrm{span}\\{\bar R\_{l\zeta}\\}\_\zeta \\),
\\( \hat V^{\mathrm{na}}\_K f = V^{\mathrm{na}}\_K f \\) holds exactly.
Indeed, writing \\( R = \sum\_\zeta a\_\zeta \bar R\_{l\zeta} \\), the orthogonality \\( (\bar R\_{l\zeta'}, \bar R\_{l\zeta})\_V = \delta\_{\zeta\zeta'} c\_{l\zeta} \\) and
the orthogonality of spherical harmonics give \\( \langle \beta^{K}\_{l'\zeta' m'} \vert f \rangle = \delta\_{ll'}\delta\_{mm'} a\_{\zeta'} c\_{l\zeta'} \\), so

$$
\hat V^{\mathrm{na}}\_K f = \sum\_{\zeta'} V^{\mathrm{na}}\_K \bar R\_{l\zeta'} S\_{lm}\\, a\_{\zeta'} = V^{\mathrm{na}}\_K f
$$

Therefore, if the basis functions on atom \\( K \\) lie in the expansion space, then for any \\( \phi \\),
\\( \langle \phi \vert \hat V^{\mathrm{na}}\_K \vert \phi\_{K\alpha} \rangle = \langle \phi \vert V^{\mathrm{na}}\_K \vert \phi\_{K\alpha} \rangle \\).
Three-center terms converge to the exact values as \\( L\_{\max} \\) and \\( N\_{\mathrm{rad}} \\) are increased (TN15 Fig. 2).

- \\( V^{\mathrm{na}} \\) does not necessarily have a definite sign (it is mostly negative for many atoms), so \\( (\cdot,\cdot)\_V \\) is indefinite as an inner product.
  \\( c\_{l\zeta} \\) can be negative or close to 0. Note that this expansion does not need to preserve the sign of \\( V^{\mathrm{na}}\_K \\).
- Because it is separable, the three-center matrix elements become sums of products of two-center integrals \\( \langle \phi \vert \beta \rangle \\),
  and the same assembly as for the [Nonlocal pseudopotential](nonlocal-pseudopotential.md) can be used (the projector energy \\( D \\) corresponds to the diagonal \\( 1/c \\)).
- \\( \hat V^{\mathrm{na}}\_K \\) is rotationally invariant (since the sum over \\( m \\) is taken for each \\( l \\)). Replacing \\( S\_{lm} \\) by complex \\( Y\_{lm} \\) gives the same operator.

### 3.2 Seed functions [Public] [Convention]

[Public] TN15 states that the pseudo-atomic orbitals of the ground and excited states, obtained under the same confinement potential as \\( n^{a} \\), are used as seed functions.

[Convention] The seed functions of the existing implementation are as follows. Let the radial functions stored in the PAO file be \\( R^{\mathrm{PAO}}\_{l\zeta} \\)
(\\( l = 0,\dots,L\_p \\), \\( \zeta = 0,\dots,M\_p - 1 \\); \\( L\_p \\) is `PAO.Lmax` and \\( M\_p \\) is `PAO.Mul`),
let \\( u(r) = V^{\mathrm{na}}\_K(r)/10 + 10^{-13} \\), and count \\( \zeta \\) from 0.

| Condition | Seed function \\( \chi\_{l\zeta}(r) \\) |
|---|---|
| \\( l \le L\_p \\), \\( \zeta < M\_p \\) | \\( R^{\mathrm{PAO}}\_{l\zeta}(r) \\) |
| \\( l \le L\_p \\), \\( \zeta \ge M\_p \\) | \\( R^{\mathrm{PAO}}\_{l0}(r)\\, u(r)^{\zeta} \\) |
| \\( l > L\_p \\), \\( \zeta < M\_p \\) | \\( R^{\mathrm{PAO}}\_{L\_p\zeta}(r)\\, r^{\\,l - L\_p} \\) |
| \\( l > L\_p \\), \\( \zeta \ge M\_p \\) | \\( R^{\mathrm{PAO}}\_{L\_p, M\_p-1}(r)\\, u(r)^{\\,\zeta - M\_p + 1} \\) |

- The PAOs used are not the basis selected for the calculation (e.g. `s2p2d1`) but all of those stored in the file, from the first one.
  The distributed PAO files have \\( M\_p = 15 \\) (only one has 20 and only one has 5), so with the default \\( N\_{\mathrm{rad}} = 12 \\),
  all seeds for \\( l \le L\_p \\) are the PAOs themselves.
- Since the basis functions are the leading radial functions of the PAO file (or linear combinations of them), if the number of selected radial functions for each \\( l \le \min(L\_p, L\_{\max}) \\)
  is at most \\( N\_{\mathrm{rad}} \\), the two-center terms are exact by the property of Section 3.1.
- The seeds for \\( l > L\_p \\) are the \\( l = L\_p \\) PAOs multiplied by \\( r^{l-L\_p} \\) to match the behavior \\( r^l \\) near the origin.
  The default \\( L\_{\max} \\) (Section 3.4) is usually larger than \\( L\_p \\), so these rows are actually used in default calculations.
- The radial functions are interpolated onto the VPS radial grid (following the rules for out-of-range values in [PAO basis](pao-basis.md)).

### 3.3 Numerical treatment of the orthogonalization [Convention] [Design]

[Convention] The key points of the procedure in the existing implementation are as follows.

- Before orthogonalization, the seed functions are normalized with the ordinary norm \\( \int r^2 \chi^2 dr \\).
- The orthogonalization is exactly the formula of Section 3.1 (classical Gram–Schmidt, using the original seeds \\( \chi\_{l\zeta} \\) in the coefficients).
- Terms with \\( \lvert c\_{l\zeta}\rvert < 10^{-15} \\) are discarded by setting \\( 1/c\_{l\zeta} = 0 \\). For empty atoms, all terms are discarded.
- Radial integrals are performed with 300 Gauss–Legendre points over the interval from the first to the last point of the VPS radial grid.
- The choice of normalization (a constant multiple of \\( \bar R \\)) cancels in \\( \lvert\beta\rangle\langle\beta\rvert/c \\), so it does not affect the expansion.

[Design] In the reimplementation, the orthogonalization may be done by modified Gram–Schmidt or by diagonalizing the Gram matrix of \\( (\cdot,\cdot)\_V \\)
(in exact arithmetic these give the same expansion). However, when the seed functions become nearly linearly dependent (\\( c \\) very small),
\\( 1/c \\) becomes large and rounding errors grow, so a safeguard such as the threshold above must always be included.

### 3.4 Truncation values [Public] [Convention]

[Public] From convergence results for the H₂O molecule, TN15 recommends \\( L\_{\max} \\) = "maximum angular momentum of the basis + 4" and \\( N\_{\mathrm{rad}} = 4 \\).
For the water molecule, convergence is reached at \\( L\_{\max} \approx 6 \\), \\( N\_{\mathrm{rad}} \approx 4 \\), and both are stated to require at least 2.

[Convention] The defaults of the existing implementation are larger than these and can be changed by input keywords.

| Keyword | Default | Meaning |
|---|---|---|
| `scf.ProExpn.VNA` | on | Whether to use the projector expansion. If off, computed on the real-space grid (described in Section 8.2 of the manual) |
| `scf.BufferL.VNA` | 6 | \\( L\_{\max} = \max\_{\text{all species}}(\text{maximum } l \text{ of the basis}) + \\) this value |
| `scf.RadialF.VNA` | 12 | \\( N\_{\mathrm{rad}} \\) |

- \\( L\_{\max} \\) and \\( N\_{\mathrm{rad}} \\) are common to all species (not per species; determined by the maximum basis \\( l \\) over all species appearing in the calculation).
  Example: if the basis goes up to `s2p2d1`, \\( L\_{\max} = 8 \\).
- The latter two keywords are not found in the keyword list of the manual.

[Design] The reimplementation matches the defaults of the existing implementation (6 and 12). Comparisons with fixtures are also done at the defaults.

### 3.5 Overlap with projectors [Public]

\\( \langle \phi\_{I\alpha} \vert \beta^{K}\_{l\zeta m} \rangle \\) is the overlap formula of [Two-center integrals](two-center-overlap-kinetic.md)
with the radial function of \\( \phi\_b \\) replaced by \\( V^{\mathrm{na}}\_K \bar R\_{K,l\zeta} \\) and the interatomic vector replaced by \\( \boldsymbol\tau\_{IK} \\) (the method of TN1).
The spherical Bessel transform of \\( V^{\mathrm{na}}\_K \bar R\_{K,l\zeta} \\) is tabulated on the same grid as in Section 1.4.
It is exactly 0 for \\( \lVert\boldsymbol\tau\_{IK}\rVert \ge r\_{c,I} + r\_{c,K} \\).

### 3.6 Three-center sum and smooth truncation [Convention]

[Convention] Blocks between different atoms are

$$
V^{\mathrm{na}}\_{I\alpha,J\beta}(\mathbf L) = d\big(\lVert\boldsymbol\tau\_{IJ\mathbf L}\rVert;\\, r\_{c,I} + r\_{c,J}\big)
\sum\_{K} \sum\_{l\zeta m} \frac{\langle \phi\_{I\alpha} \vert \beta^{K}\_{l\zeta m} \rangle \langle \beta^{K}\_{l\zeta m} \vert \phi\_{J\beta}(\cdot-\boldsymbol\tau\_{IJ\mathbf L}) \rangle}{c\_{K,l\zeta}}
$$

\\( K \\) (including periodic images) runs over atoms within atom-pair range of both \\( I \\) and \\( J \\) ([sparse block format](../design/sparse-block-format.md);
distance at most the sum of \\( r\_c \\)). Outside this range one of the overlaps is 0, so the truncation of the sum is exact.

The truncation factor \\( d \\) is a \\( C^2 \\) step function with a buffer width of 1 bohr. With \\( x = R - R\_c \\),

$$
d(R; R\_c) =
\begin{cases}
1 & R \le R\_c - 1 \\\\
{}-x^3\\,(10 + 15x + 6x^2) & R\_c - 1 < R \le R\_c \\\\
0 & R > R\_c
\end{cases}
$$

(\\( R \\) in bohr). It is the quintic polynomial with value 1 and vanishing first and second derivatives at \\( x = -1 \\), and vanishing value, first, and second derivatives at \\( x = 0 \\).

- **Reason**: the exact matrix element \\( \langle\phi\_I\vert V^{\mathrm{na}}\_K\vert\phi\_J\rangle \\) is 0 for \\( \lVert\boldsymbol\tau\_{IJ\mathbf L}\rVert \ge r\_{c,I} + r\_{c,J} \\),
  where the supports of \\( \phi\_I \\) and \\( \phi\_J \\) do not overlap, but the separable approximation loses this property,
  and the value drops discontinuously to 0 at the edge of the atom-pair range. The factor \\( d \\) makes the energy and forces continuous as pairs enter or leave the range.
- It is not applied to blocks on the same atom (the direct computation of Section 2) (also, \\( d(0) = 1 \\)).
- Since the factor multiplies the matrix elements themselves, it appears consistently in the Hamiltonian, `Una`, and the forces. For the forces, its derivative is also needed (Section 7).

## 4. Basis pairs on the same atom and \\( V^{\mathrm{na}} \\) of another atom [Public] [Convention]

The formula for the case that TN15 treats in wavenumber space (\\( I = J \\), \\( \mathbf L = 0 \\), arbitrary \\( K \\)) is derived in the conventions of this note.

Take atom \\( I \\) as the origin, \\( \boldsymbol\tau = \boldsymbol\tau\_{IK} \\),
\\( \phi\_a = R\_a(r) Y\_{lm} \\), \\( \phi\_b = R\_b(r) Y\_{l'm'} \\). Expanding the product of the angular parts as

$$
Y^{\ast}\_{lm} Y\_{l'm'} = \sum\_{LM} G(lm;\\, l'm',\\, LM)\\, Y^{\ast}\_{LM}
$$

(\\( G \\) is the Gaunt coefficient of [Two-center integrals](two-center-overlap-kinetic.md); the sum is only over \\( M = m - m' \\)),
writing \\( V^{\mathrm{na}}\_K(\lvert\mathbf r - \boldsymbol\tau\rvert) \\) as a Fourier integral and using the plane-wave expansion twice, one obtains

$$
\langle \phi\_a \vert V^{\mathrm{na}}\_K(\cdot - \boldsymbol\tau) \vert \phi\_b \rangle
{}= 8 \sum\_{LM} G(lm;\\, l'm',\\, LM)\\, Y^{\ast}\_{LM}(\hat{\boldsymbol\tau})
\int\_0^\infty dk\\, k^2 j\_L(k\tau)\\, \tilde P^{ab}\_L(k)\\, \tilde V^{\mathrm{na}}\_K(k)
$$

$$
\tilde P^{ab}\_L(k) = \int\_0^\infty dr\\, r^2 j\_L(kr)\\, R\_a(r) R\_b(r)
$$

Key points of the derivation: \\( \int d^3r\\, f(r) Y^{\ast}\_{LM}(\hat{\mathbf r}) e^{i\mathbf k\cdot\mathbf r} = 4\pi i^L Y^{\ast}\_{LM}(\hat{\mathbf k}) \tilde f\_L(k) \\)
cancels against the coefficient \\( (-i)^L \\) of the expansion of \\( e^{-i\mathbf k\cdot\boldsymbol\tau} \\), so no phase remains.
The coefficient 8 is \\( (4\pi)^3/(2\pi)^3 \\) (including the \\( 4\pi \\) of the three-dimensional Fourier transform \\( 4\pi\tilde V^{\mathrm{na}}\_K \\)).

- \\( \tau = 0 \\) (\\( K = I \\), one-center): since \\( j\_L(0) = \delta\_{L0} \\),
  \\( \langle\phi\_a\vert V^{\mathrm{na}}\_I\vert\phi\_b\rangle = \delta\_{ll'}\delta\_{mm'}\int\_0^\infty dr\\, r^2 R\_a R\_b V^{\mathrm{na}}\_I \\).
  The direction \\( \hat{\boldsymbol\tau} \\) is undefined, so this case is treated specially.
- In the real spherical harmonics basis, the result is made real by the same transformation as in Section 4 of [Two-center integrals](two-center-overlap-kinetic.md).
- It is exactly 0 for \\( \tau \ge r\_{c,I} + r\_{c,K} \\).
- [Convention] The existing implementation tabulates \\( \tilde P^{ab}\_L \\) on the same Gauss–Legendre wavenumber grid as \\( \tilde V^{\mathrm{na}}\_K \\), and performs the \\( k \\) integral on the same grid.

> **Check**: For Gaussian radial functions and an analytic model \\( V^{\mathrm{na}} \\) (see the test section below),
> the formula above agrees with direct integration on a three-dimensional grid to 6 or more digits for 4 combinations of \\( (l, l') \\) (item (2) of `tools/check_neutral_atom.py`).

## 5. Evaluation on the grid (`scf.ProExpn.VNA off`) [Public] [Convention]

[Public] According to the manual (Section 8.2), setting `scf.ProExpn.VNA` off computes the matrix elements of \\( V^{\mathrm{na}} \\) on the real-space grid.

[Convention] In this case, \\( \sum\_K V^{\mathrm{na}}\_K(\lvert\mathbf r\_p - \mathbf R\_K\rvert) \\) is built on the grid points,
added to \\( \delta V\_H + v\_{\mathrm{xc},\sigma} \\) as an effective potential, and turned into matrix elements by the method of Section 5 of [Grid](grid-poisson-xc.md).
`Una` is \\( \Delta V\sum\_p n(\mathbf r\_p) V^{\mathrm{na}}(\mathbf r\_p) \\) ([Total energy and forces](energy-force.md)).
With the default on, \\( V^{\mathrm{na}} \\) is not put on the grid.

TN15 (Fig. 3, Table I) shows that with grid evaluation the energy oscillates depending on the position of the atoms relative to the grid (egg-box effect),
and that the projector expansion reduces this substantially.

## 6. Role in the total-energy decomposition [Public]

The derivation of the total-energy decomposition is in [Total energy and forces](energy-force.md). Here only the points relating to the neutral-atom potential are summarized (TN1 Eq. 18–23).

- \\( E\_{\mathrm{na}} = \int n \sum\_K V^{\mathrm{na}}\_K \\) is, with the projector expansion, the sum of products of the density matrix and the matrix elements of Sections 2–3 (`Una`).
- The long-range part of the electrostatic energy remains only in
  \\( E\_{\delta\mathrm{ee}} = \frac12\int\delta n\\,\delta V\_H \\) of the difference density \\( \delta n = n - \sum\_I n^{a}\_I \\) and its Hartree potential \\( \delta V\_H \\), which is solved by FFT ([Grid](grid-poisson-xc.md)).
  \\( n^{a}\_I \\) is the same density renormalized in Section 1.2.
- The repulsion between ionic cores and the Hartree terms between atomic densities combine into the screened short-range pair potential

$$
E\_{\mathrm{scc}} = \frac12\sum\_{I,J}\left[\frac{Z\_I Z\_J}{\lVert\mathbf R\_I - \mathbf R\_J\rVert}(1 - \delta\_{IJ}) - \int d^3r\\, n^{a}\_I\\, V^{a}\_{H,J}\right]
$$

  (TN1 Eq. 22; \\( J \\) includes periodic images, and \\( I = J \\) means the same image).
  Each pair term vanishes at distances where the supports of the atomic densities do not overlap. The evaluation in the existing implementation (`Ucore` and `UH0`, including normalization at the endpoint) is in
  the "Screened interatomic term" section of [Total energy and forces](energy-force.md).

[Design] The Coulomb integral of two spherically symmetric charges can also be computed as the one-dimensional integral

$$
\int d^3r\\, n^{a}\_I(\mathbf r)\\, V^{a}\_{H,J}(\mathbf r - \mathbf d) = \frac{2}{\pi}\int\_0^\infty dk\\, j\_0(kd)\\, N\_I(k) N\_J(k),\qquad
N\_I(k) = 4\pi\int\_0^\infty dr\\, r^2 j\_0(kr)\\, n^{a}\_I(r)
$$

(since \\( N\_I(0) = Z\_I \\), it tends to \\( Z\_I Z\_J/d \\) at large distance). It can be used as an alternative when building a distance table.
Agreement with the analytic solution to \\( 10^{-14} \\) for Gaussian charges is confirmed (item (1) of `tools/check_neutral_atom.py`).

## 7. Forces [Public] [Convention]

The overall form (derivatives of matrix elements at fixed density matrix, and the Pulay term with the energy-weighted density matrix) is in
the "Forces" section of [Total energy and forces](energy-force.md) (TN1 Eq. 59–61, 75). Here the derivatives of the neutral-atom matrix elements are given.
\\( g\_s \\) and the density-matrix symbols follow that chapter.

### 7.1 Projector-expansion blocks

Since \\( A^{K}\_{I\alpha,p} = \langle \phi\_{I\alpha} \vert \beta^{K}\_p \rangle \\) (\\( p = (l\zeta m) \\)) is a function of \\( \boldsymbol\tau\_{IK} \\) only,

$$
\frac{\partial A^{K}\_{I\alpha,p}}{\partial \mathbf R\_K} = \nabla\_{\boldsymbol\tau} A^{K}\_{I\alpha,p}(\boldsymbol\tau\_{IK}),\qquad
\frac{\partial A^{K}\_{I\alpha,p}}{\partial \mathbf R\_I} = -\nabla\_{\boldsymbol\tau} A^{K}\_{I\alpha,p}(\boldsymbol\tau\_{IK})
$$

(\\( \nabla\_{\boldsymbol\tau} \\) is the analytic derivative of Section 6 of [Two-center integrals](two-center-overlap-kinetic.md); when \\( K = I \\) the two cancel to 0).
Writing the matrix element of Section 3.6 as \\( V = d\\, M \\), \\( M = \sum\_K\sum\_p A^{K}\_{I\alpha,p} c\_p^{-1} A^{K}\_{J\beta,p} \\), the derivative with respect to atom \\( X \\) is

$$
\frac{\partial V^{\mathrm{na}}\_{I\alpha,J\beta}}{\partial\mathbf R\_X}
{}= d \sum\_K\sum\_p \frac{1}{c\_p}\left[\frac{\partial A^{K}\_{I\alpha,p}}{\partial\mathbf R\_X} A^{K}\_{J\beta,p} + A^{K}\_{I\alpha,p}\frac{\partial A^{K}\_{J\beta,p}}{\partial\mathbf R\_X}\right]
{}+ M\\, d'(\tau\_{IJ})\\, \frac{\partial \tau\_{IJ}}{\partial \mathbf R\_X}
$$

Here \\( \tau\_{IJ} = \lVert\boldsymbol\tau\_{IJ\mathbf L}\rVert \\), \\( \partial\tau\_{IJ}/\partial\mathbf R\_J = \hat{\boldsymbol\tau}\_{IJ\mathbf L} = -\partial\tau\_{IJ}/\partial\mathbf R\_I \\), and
\\( d'(R) = -30x^2(1+x)^2 \\) (for \\( R\_c - 1 < R \le R\_c \\); 0 otherwise).
The contribution to the energy is the sum of this multiplied by \\( g\_s\sum\_\sigma\rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L) \\).
[Convention] The existing implementation includes the term with the derivative of the truncation factor.

### 7.2 Blocks on the same atom

\\( W^{IK}\_{\alpha\beta}(\boldsymbol\tau\_{IK}) = \langle \phi\_{I\alpha} \vert V^{\mathrm{na}}\_K \vert \phi\_{I\beta} \rangle \\) (the formula of Section 4) is a function of \\( \boldsymbol\tau\_{IK} \\) only.

$$
\frac{\partial W^{IK}}{\partial \mathbf R\_K} = \nabla\_{\boldsymbol\tau} W^{IK},\qquad
\frac{\partial W^{IK}}{\partial \mathbf R\_I} = -\nabla\_{\boldsymbol\tau} W^{IK}\quad(K \ne I)
$$

and the \\( K = I \\) term is a configuration-independent constant. For \\( \nabla\_{\boldsymbol\tau} \\), the gradient formula of Section 6 of [Two-center integrals](two-center-overlap-kinetic.md) is applied as is
to each term \\( f\_L(\tau) Y^{\ast}\_{LM}(\hat{\boldsymbol\tau}) \\) of the Section 4 formula (\\( f\_L(\tau) = \int dk\\, k^2 j\_L(k\tau)\tilde P^{ab}\_L\tilde V^{\mathrm{na}}\_K \\)),
using \\( f\_L'(\tau) = \int dk\\, k^3 j\_L'(k\tau)\tilde P^{ab}\_L\tilde V^{\mathrm{na}}\_K \\).
[Convention] The existing implementation also uses this analytic derivative, and drops the \\( \varphi \\)-direction term when \\( \sin\theta \\) is small.

With this, the "derivatives of matrix elements for combinations not using the projector expansion" listed under [Open] in [Total energy and forces](energy-force.md)
are determined by the formulas above for the default `scf.ProExpn.VNA on`.

### 7.3 Others

- The force for grid evaluation (Section 5) is the grid term (fourth term) of [Total energy and forces](energy-force.md).
- The force of \\( E\_{\mathrm{scc}} \\) is the "derivative of the screened interatomic term" in [Total energy and forces](energy-force.md).
- By translational invariance, for each block the sum of the derivatives with respect to the atoms involved is 0 (sum over \\( X = I, J, K \\) in Section 7.1, over \\( X = I, K \\) in Section 7.2).

## 8. Summary of numerical parameters [Convention]

| Item | Value |
|---|---|
| Truncation of \\( V^{\mathrm{na}} \\) | 0 at \\( r\_c \\) (`radial.cutoff.pao`). Steepness of the Fermi-type factor 20 bohr⁻¹ |
| Asymptotic correction | \\( V^{\mathrm{loc}} + sV^{a}\_H = 0 \\) at the VPS grid point closest to \\( r\_c \\) |
| Projector expansion | \\( L\_{\max} = \max(l\_{\mathrm{basis}}) + 6 \\), \\( N\_{\mathrm{rad}} = 12 \\) (`scf.BufferL.VNA`, `scf.RadialF.VNA`) |
| Radial integrals for orthogonalization | Full range of the VPS grid, 300 Gauss–Legendre points |
| Discarding small \\( c \\) | \\( \lvert c\rvert < 10^{-15} \\) |
| Spherical Bessel transform | Real space \\( r \le r\_c + 0.5 \\) bohr; wavenumbers \\( [10^{-6}, \sqrt{E\_{\mathrm{1DFFT}}}] \\) with 300 Gauss–Legendre points |
| Truncation factor | \\( C^2 \\) quintic polynomial, buffer 1 bohr, \\( R\_c = r\_{c,I} + r\_{c,J} \\) |

## Tests

### Analytic model [Design]

From a Gaussian atomic density and an error-function local potential,

$$
n^{a}(r) = Z\left(\frac{\alpha}{\pi}\right)^{3/2} e^{-\alpha r^2},\qquad
V^{\mathrm{loc}}(r) = -Z\frac{\mathrm{erf}(\sqrt{\beta}\\,r)}{r}\quad(\beta > \alpha)
$$

one obtains \\( V^{\mathrm{na}}(r) = Z\\,[\mathrm{erf}(\sqrt\alpha\\, r) - \mathrm{erf}(\sqrt\beta\\, r)]/r \\).
It does not have finite support, but since it decays like a Gaussian, it can be used to test everything except the truncation of Section 1.3.

| Kind | Content | Tolerance (guideline) |
|---|---|---|
| Analytic | \\( V^{a}\_H \\) from the radial integral of Section 1.1 vs. the model \\( Z\\,\mathrm{erf}(\sqrt\alpha r)/r \\) | \\( 10^{-10} \\) |
| Analytic | Coulomb integral of spherical charges in Section 6 vs. \\( Z\_1 Z\_2\\,\mathrm{erf}(\sqrt{\mu}d)/d \\) (\\( \mu = a\_1 a\_2/(a\_1+a\_2) \\)) | \\( 10^{-12} \\) |
| Analytic | Formula of Section 4 vs. direct integration on a three-dimensional grid (Gaussian \\( R\_a, R\_b \\), \\( l, l' \le 2 \\)) | \\( 10^{-6} \\) (determined by grid accuracy) |
| Analytic | \\( \tau \to 0 \\) limit of the Section 4 formula vs. \\( \int r^2 R\_a R\_b V^{\mathrm{na}} \\) | \\( 10^{-10} \\) |
| Analytic | Truncation factor: value, first and second derivatives at \\( x = -1, 0 \\) | rounding error |

### Projector expansion [Design]

| Kind | Content | Tolerance (guideline) |
|---|---|---|
| Invariant | \\( (\bar R\_{l\zeta}, \bar R\_{l\zeta'})\_V = \delta\_{\zeta\zeta'} c\_{l\zeta} \\) | \\( 10^{-12} \\) (relative) |
| Invariant | \\( \hat V^{\mathrm{na}}\chi = V^{\mathrm{na}}\chi \\) for seed functions \\( \chi \\): radial inner products \\( \sum\_\zeta (\bar R\_\zeta, \chi)\_V (\bar R\_\zeta, g)\_V / c\_\zeta = (\chi, g)\_V \\) for arbitrary \\( g \\) | \\( 10^{-10} \\) (relative) |
| Invariant | Two-center (\\( K = I \\)) expansion value vs. direct three-dimensional integration | grid accuracy; 6 digits in item (3) of `tools/check_neutral_atom.py` |
| Convergence | Three-center values approach direct integration as \\( L\_{\max} \\), \\( N\_{\mathrm{rad}} \\) increase | trend check (example below) |
| Invariant | Matrix-element symmetry \\( V^{\mathrm{na}}\_{I\alpha,J\beta}(\mathbf L) = V^{\mathrm{na}}\_{J\beta,I\alpha}(-\mathbf L) \\) | rounding error |
| Invariant | Under rotation of the structure, the matrix transforms by the rotation matrices of real spherical harmonics | \\( 10^{-10} \\) |
| Invariant | At the edge of the atom-pair range (\\( \tau = r\_{c,I} + r\_{c,J} \\)), blocks between different atoms are 0, and smooth within 1 bohr inside | exact |

Convergence example (item (3) of `tools/check_neutral_atom.py`; s- and p-type Gaussian functions, \\( V^{\mathrm{na}} \\) on a third center):
with \\( N\_{\mathrm{rad}} = 4 \\) and \\( L\_{\max} = 0, 2, 4, 6 \\), the relative errors with respect to direct integration are approximately
\\( 6\times10^{-2} \\), \\( 6\times10^{-3} \\), \\( 2\times10^{-4} \\), \\( 1\times10^{-5} \\).
The effect of increasing \\( N\_{\mathrm{rad}} \\) from 2 to 6 is of order \\( 10^{-3} \\); that the effect of \\( L\_{\max} \\) is larger agrees with the trend of TN15 Fig. 2.

### Tests on actual data [Design]

| Kind | Content | Tolerance (guideline) |
|---|---|---|
| Invariant | \\( 4\pi\int r^2 n^{a} = Z \\) after renormalization | \\( 10^{-10} \\) |
| Invariant | \\( s \\) of Section 1.3 is close to 1 (distributed files for C, O, H, Fe, etc.) | \\( \lvert s-1\rvert < 10^{-3} \\) |
| Invariant | \\( \lvert V^{\mathrm{loc}} + sV^{a}\_H\rvert \\) is small within 1 bohr inside \\( r\_c \\) | \\( 10^{-3} \\) Hartree |
| Finite difference | Central differences of `Una` at fixed \\( \rho \\) under atomic displacements vs. the analytic forces of Section 7 (including the derivative of the truncation factor) | \\( 10^{-7} \\) Hartree/bohr |
| Finite difference | The above agreement holds also for displacements crossing the edge of the atom-pair range | same as above |
| Invariant | The force contributions of `Una` sum to 0 | \\( 10^{-10} \\) |

### Comparison with fixtures [Design]

- Dump and compare \\( V^{\mathrm{na}}\_I(r) \\), \\( \bar R\_{l\zeta} \\), \\( c\_{l\zeta} \\) (per species) of the existing implementation.
  Since the normalization of \\( \bar R \\) may differ between implementations, compare the combination \\( \beta\beta/c \\) (e.g. the radial kernel \\( \sum\_\zeta V\bar R\_\zeta(r)V\bar R\_\zeta(r')/c\_\zeta \\)).
- Matrix elements: compare \\( V^{\mathrm{na}} \\) for same-atom blocks and different-atom blocks in the H₂O molecule (`H6.0-s1p1`, `O6.0-s1p1`, etc.) and Si crystal.
- Energy: `Una` (both `scf.ProExpn.VNA` on and off).
- Tolerances: \\( 10^{-8} \\) Hartree for \\( V^{\mathrm{na}}(r) \\), \\( 10^{-7} \\) Hartree for matrix elements, and \\( 10^{-6} \\) Hartree for `Una` as guidelines
  (with the same \\( L\_{\max} \\), \\( N\_{\mathrm{rad}} \\), and the same wavenumber grid).

## [Open]

- Details of the numerical integration of \\( V^{a}\_H \\) in Section 1.3 and the quadrature of the spherical Bessel transforms in Sections 1.4 and 3.5 (the existing implementation uses a dedicated quadrature for each).
  We expect them to change results only within tolerance, but this is not settled until confirmed by comparison with fixtures.
- Whether \\( \tilde P^{ab}\_L \\) of Section 4 is built for each combination of selected basis functions or for all radial functions of the PAO file (the result is the same; only efficiency is affected).
- Whether, in calculations using orbital optimization (contraction), the condition for the contracted basis functions to lie in the expansion space (the number of primitive functions used in the contraction is at most \\( N\_{\mathrm{rad}} \\)) is always satisfied.
- Treatment with noncollinear spin and spin-orbit coupling (since \\( V^{\mathrm{na}} \\) is spin-independent, we expect that the same term is simply added to the spin-diagonal components, but this is unconfirmed).
- Contribution to the stress (the stress section of [Total energy and forces](energy-force.md) states that projector-expansion terms are handled per two-center factor, but the treatment of the truncation factor \\( d \\) and of the direct-computation terms of Section 4 is unconfirmed).
