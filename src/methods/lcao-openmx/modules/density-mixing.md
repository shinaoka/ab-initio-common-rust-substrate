# Density Matrix and Charge Mixing

> Status: draft (awaiting review). The formulas for occupations, the density matrix, mixing, convergence criteria and extrapolation are technical notes 1, 8 and 12 and
> manual 4.0 (Sec. 8.2, Ch. 16, Ch. 17) rewritten in the notation of this note.
> Rules not written in the public documents (defaults, normalization of the residual norm, weight update rules, combination of convergence criteria, etc.)
> were confirmed as behavior of the existing implementation and are marked [Convention].
> Numerical properties are checked with `tools/check_density_mixing.py`.

## Scope of this chapter

This chapter covers the part of the SCF loop from after the eigenvalue problem is solved until the Hamiltonian of the next iteration is built.

1. Determine the occupations and the chemical potential from the eigenvalues.
2. Build the real-space density matrix \\( \rho(\mathbf L) \\) and the energy-weighted density matrix \\( E(\mathbf L) \\) from the eigenvectors and occupations.
3. Build the electron density \\( n(\mathbf r) \\) on the grid from \\( \rho \\).
4. Use the difference between input and output (the residual) to build the input of the next iteration (density matrix, density, potential, or Hamiltonian) (mixing).
5. Check convergence.

The eigenvalue problem itself is in [Generalized eigenvalue problem and regularization of S](generalized-eigen.md), k points and Bloch sums are in [k points and Bloch sums](kpoints-bloch.md),
and grid integration and the Poisson equation are in [Real-space grid, Poisson, exchange-correlation](grid-poisson-xc.md).

## Overview of the SCF loop [Public] [Convention]

[Public] According to technical note 12 and manual Sec. 16.1, the SCF repeats "build the Hamiltonian from the input density (or density matrix) -> solve the eigenvalue problem ->
build the output density matrix -> mix input and output to build the next input", and the mixed quantity and method are selected with `scf.Mixing.Type`.

[Convention] The order of one iteration in the existing implementation is as follows.

1. Build the Hartree and exchange-correlation potentials from the input density \\( n^{\mathrm{in}} \\) and assemble the Hamiltonian \\( H(\mathbf L) \\).
   In method RMM-DIISV the potential is mixed at this stage; in method RMM-DIISH the Hamiltonian is mixed right after this.
2. Solve the eigenvalue problem at each k point and spin.
3. Determine the chemical potential, and build the occupations, the band energy \\( U\_{\mathrm{ele}} \\), the output density matrix \\( \rho^{\mathrm{out}} \\), and the energy-weighted density matrix.
4. Check convergence (see below). If converged, exit the loop.
5. Depending on the method, mix the density matrix or the density and build the input density \\( n^{\mathrm{in}} \\) of the next iteration on the grid.

Since the convergence check comes before mixing, the residual norm used in the check is the value from the most recent mixing (the residual of the previous iteration).

## 1. Occupations and chemical potential

### Number of electrons [Public] [Convention]

[Public] According to the manual (Sec. 8.2 `scf.system.charge`, Ch. 31), `scf.system.charge` is the extra charge of the system; a positive value means hole doping, a negative value electron doping, and the default is 0.

[Convention] The number of electrons to be occupied is

$$
N\_{\mathrm{e}} = \sum\_{I} Z\_I^{\mathrm{val}} - Q,\qquad Q = \texttt{scf.system.charge}
$$

where \\( Z\_I^{\mathrm{val}} \\) is the number of valence electrons (core charge) of the pseudopotential of atom \\( I \\).

### Fermi–Dirac distribution [Public] [Convention]

[Public] According to manual Sec. 8.2, `scf.ElectronicTemperature` is the electronic temperature (K), with default 300 K.

[Convention] The occupations follow the Fermi–Dirac distribution

$$
f(\varepsilon) = \frac{1}{1 + e^{x}},\qquad x = \frac{\varepsilon - \mu}{k\_{\mathrm B} T}
$$

No other smearing functions (Gaussian, Methfessel–Paxton, etc.) are available.

- The temperature is converted to Hartree units with \\( k\_{\mathrm B} = 8.617251324\times 10^{-5} \\) eV/K and 1 Hartree = 27.2113845 eV.
- To avoid exponent overflow, \\( x \\) is clipped to \\( [-60, 60] \\) before evaluation (for Band).
- States whose occupation is smaller than about \\( 10^{-13} \\) may be omitted from the density-matrix sum (the effect on the result is at the level of rounding error).

### Chemical potential [Convention]

\\( \mu \\) is determined from the electron-number condition

$$
N(\mu) = g\_s \sum\_{\mathbf k} w\_{\mathbf k} \sum\_{\sigma}\sum\_{\nu=1}^{N\_{\mathrm{st}}} f\left(\varepsilon\_{\nu\sigma}(\mathbf k)\right) = N\_{\mathrm e}
$$

where \\( w\_{\mathbf k} \\) are the k-point weights normalized to \\( \sum\_{\mathbf k} w\_{\mathbf k} = 1 \\) ([k points and Bloch sums](kpoints-bloch.md)),
and \\( g\_s \\) is 2 for non-spin-polarized and 1 for collinear spin-polarized and noncollinear calculations.

- **The solver is bisection.** The initial interval is \\( [-20, 20] \\) Hartree (\\( [-30, 30] \\) for Cluster),
  the stopping condition is \\( \lvert N(\mu) - N\_{\mathrm e}\rvert < 10^{-12} \\), and the maximum number of iterations is 1000.
  Since \\( N(\mu) \\) is monotonically increasing, bisection always converges.
- [Convention] This interval and tolerance are the values for Cluster and Band. The O(N) solvers (`DC`, `Krylov`, `DC-LNO`) use their own initial interval and tolerance per solver ("Number of electrons and chemical potential" in [O(N) methods](on-methods.md)).
  In the central region of NEGF, no search under the electron-number constraint is done; the number of electrons is determined by the chemical potentials of the leads (only the Band stage follows a different rule; [Electron transport with NEGF](negf.md)).
  Also, \\( \rho \\) of the O(N) solvers does not satisfy the Hermiticity below (\\( \rho\_{I\alpha,J\beta}(\mathbf L) \ne \rho\_{J\beta,I\alpha}(-\mathbf L) \\); [O(N) methods](on-methods.md)).
- The existing implementation first finds \\( \mu \\) at 1200 K, uses it as the first trial point, and then performs bisection again from the interval \\( [-20, 20] \\) at the specified temperature.
  Since the final result is determined solely by the electron-number condition at the specified temperature, a single bisection (or a combination with Newton's method) suffices in a reimplementation.
- In collinear spin-polarized calculations, a **common** \\( \mu \\) is used for up and down spins. The number of electrons per spin is not fixed; only the total number of electrons is matched.
- \\( N\_{\mathrm{st}} \\) is the number of eigenstates computed; the existing implementation computes \\( N\_{\mathrm e}/2 \\) plus extra states (bounded by the number of basis functions).
  If the occupation of the \\( N\_{\mathrm{st}} \\)-th state is not negligible at high temperature, the number of states must be increased ([Design]: warn when the occupation exceeds a threshold).

> **Note**: For an insulator with \\( k\_{\mathrm B}T \\) sufficiently smaller than the gap, \\( N(\mu) \\) is nearly flat inside the gap.
> Since the stopping condition is on the electron number, the resulting \\( \mu \\) stops somewhere in the gap, not necessarily at its center
> (`tools/check_density_mixing.py`: for a two-level system at 300 K with a gap of 0.4 Hartree, it stopped at a point 0.1 Hartree away from the center).
> Occupations and density are not affected, but when \\( \mu \\) is used in output or post-processing, this indeterminacy must be included in the tolerance.

### Fixed spin moment [Convention]

No input keyword was found in the existing implementation that fixes the numbers of up- and down-spin electrons separately in collinear calculations (i.e., uses a separate chemical potential per spin).
To control the spin configuration, the spin-direction constraint of noncollinear calculations (`scf.Constraint.NC.Spin`) or the Zeeman term (`scf.NC.Zeeman.Spin`) is used.
These add terms to the Hamiltonian; the occupations are determined in the same way as above.

### Band energy [Convention]

$$
U\_{\mathrm{ele}} = g\_s \sum\_{\mathbf k} w\_{\mathbf k} \sum\_{\sigma\nu} f\left(\varepsilon\_{\nu\sigma}(\mathbf k)\right) \varepsilon\_{\nu\sigma}(\mathbf k)
$$

is computed at every iteration. It is used in the convergence check and in adjusting the mixing weight.

### Automatic control of the electronic temperature [Convention]

When `scf.Mixing.Control.Temp` (logical, default off) is set to on, the electronic temperature is temporarily raised when the SCF oscillates.

- Let \\( \lVert R\rVert \\) be the residual norm (NormRD, see below). If the average of \\( \lVert R\rVert^2 \\) over the last 5 iterations is smaller than the current \\( \lVert R\rVert^2 \\) (i.e., the residual is increasing),
  \\( \lVert R\rVert > 0.01 \\), and more than 6 iterations have been done, the temperature is set to 10 times the specified value (only once per SCF).
- When \\( \lVert R\rVert \le 0.01 \\) afterwards, the temperature is reset to the specified value.
- Each time the temperature is changed, the history of Pulay-type mixing is reset and the procedure is restarted from `scf.Mixing.StartPulay`.

## 2. Real-space density matrix

### Definition and phase [Public] [Convention]

[Public] Consistent with Eq. 4 of technical note 1 and the convention of [k points and Bloch sums](kpoints-bloch.md),
\\( O(\mathbf k) = \sum\_{\mathbf L} e^{+i\mathbf k\cdot\mathbf L} O(\mathbf L) \\), the density matrix is built as

$$
\rho^{\sigma}\_{I\alpha, J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k} \sum\_{\nu} f\_{\nu\sigma}(\mathbf k)\\,
\mathrm{Re}\left[ e^{+i\mathbf k\cdot\mathbf L}\\, \overline{c\_{I\alpha,\nu\sigma}(\mathbf k)}\\, c\_{J\beta,\nu\sigma}(\mathbf k) \right]
$$

where \\( c \\) are eigenvectors normalized as \\( c^\dagger S(\mathbf k) c = 1 \\), and \\( \mathbf L \\) is the periodic image of the column-side (\\( J \\)) atom.

[Convention]

- As in this formula, the existing implementation accumulates only the real part multiplied by \\( \cos(\mathbf k\cdot\mathbf L) \\) and \\( \sin(\mathbf k\cdot\mathbf L) \\),
  and divides by the sum of k-point weights at the end (integer weights are used, so normalization is done all at once at the end).
- Atomic positions do not enter the phase.
- \\( \rho \\) is built only on the same set of atom pairs \\( (I, J, \mathbf L) \\) as the [sparse block format](../design/sparse-block-format.md).
  Both \\( (I, J, \mathbf L) \\) and \\( (J, I, -\mathbf L) \\) are stored.
- In non-spin-polarized calculations, only the density matrix for one spin (occupations \\( f \le 1 \\)) is stored. The total number of electrons and the total density are twice that.
- For Cluster (Γ point only), \\( w = 1 \\) and the phase is 1.

[Design] With time-reversal symmetry (\\( O(\mathbf L) \\) real in collinear calculations), the contributions of \\( \pm\mathbf k \\) are complex conjugates of each other, so
using \\( P(\mathbf k) = \sum\_\nu f\_\nu c\_\nu c\_\nu^\dagger \\) one can also write

$$
\rho(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k}\\, e^{-i\mathbf k\cdot\mathbf L}\\, P(\mathbf k)
$$

This is exactly the inverse transform of the Bloch sum \\( P(\mathbf k) = \sum\_{\mathbf L} e^{+i\mathbf k\cdot\mathbf L}\rho(\mathbf L) \\) (when the k-point sum is exact).
Hence products of real-space matrices correspond to products in k space, and the relations among \\( \rho \\), \\( S \\), \\( H \\) (the sum rules and idempotency below) can be tested by real-space convolutions.

### Noncollinear case [Convention]

In noncollinear calculations, from the spinor eigenvectors (first half: up-spin coefficients, second half: down-spin coefficients),

$$
\rho^{\sigma\sigma'}\_{I\alpha, J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k} \sum\_{\nu} f\_{\nu}(\mathbf k)\\,
e^{+i\mathbf k\cdot\mathbf L}\\, \overline{c\_{\sigma I\alpha,\nu}(\mathbf k)}\\, c\_{\sigma' J\beta,\nu}(\mathbf k)
$$

is built. It is complex, and the existing implementation stores the following components.

- The four components \\( \mathrm{Re}\\,\rho^{\uparrow\uparrow} \\), \\( \mathrm{Re}\\,\rho^{\downarrow\downarrow} \\), \\( \mathrm{Re}\\,\rho^{\uparrow\downarrow} \\), \\( \mathrm{Im}\\,\rho^{\uparrow\downarrow} \\) (always).
- \\( \mathrm{Im}\\,\rho^{\uparrow\uparrow} \\), \\( \mathrm{Im}\\,\rho^{\downarrow\downarrow} \\) (when any of spin-orbit coupling, LDA+U, spin constraint, or Zeeman term is present).

\\( \rho^{\downarrow\uparrow} \\) is obtained from the Hermiticity \\( \rho^{\sigma'\sigma}\_{J\beta,I\alpha}(-\mathbf L) = \overline{\rho^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L)} \\).
This is the same formula that [k points and Bloch sums](kpoints-bloch.md) lists as the "natural form", and the existing implementation also uses this form.
The imaginary parts of the spin-diagonal components are nonzero only when time reversal is broken, and do not contribute to the density on the grid (the product of real basis functions is symmetric under \\( a\leftrightarrow b \\)), but
they are used for orbital magnetic moments and the LDA+U occupation matrix.

### Energy-weighted density matrix [Convention]

For the force calculation (the term with the derivative of the overlap matrix), the matrix of the same form multiplied by the eigenvalues,

$$
E^{\sigma}\_{I\alpha, J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k} \sum\_{\nu} f\_{\nu\sigma}(\mathbf k)\\, \varepsilon\_{\nu\sigma}(\mathbf k)\\,
\mathrm{Re}\left[ e^{+i\mathbf k\cdot\mathbf L}\\, \overline{c\_{I\alpha,\nu\sigma}(\mathbf k)}\\, c\_{J\beta,\nu\sigma}(\mathbf k) \right]
$$

is built at the same time. Its use is described in [Total energy and forces](energy-force.md). \\( E \\) is not mixed (it is always the output of the last eigenvalue problem).

### Sum rules [Design]

The following follow immediately from the definitions (abbreviating \\( \sum\_{\mathbf L}\sum\_{I\alpha J\beta} \\) as \\( \sum \\)). They are used as test invariants.

$$
g\_s \sum\_{\sigma} \sum \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\, S\_{I\alpha,J\beta}(\mathbf L) = N\_{\mathrm e},\qquad
g\_s \sum\_{\sigma} \sum \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\, H^{\sigma}\_{I\alpha,J\beta}(\mathbf L) = U\_{\mathrm{ele}},\qquad
g\_s \sum\_{\sigma} \sum E^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\, S\_{I\alpha,J\beta}(\mathbf L) = U\_{\mathrm{ele}}
$$

The first follows from \\( c^\dagger S(\mathbf k) c = 1 \\), the second from \\( c^\dagger H(\mathbf k) c = \varepsilon \\).
Furthermore, for an insulator at \\( T = 0 \\) (occupations 0 or 1), \\( P S P = P \\) in k space gives in real space

$$
\sum\_{\mathbf L\_1 + \mathbf L\_2 + \mathbf L\_3 = \mathbf L} \rho(\mathbf L\_1)\\, S(\mathbf L\_2)\\, \rho(\mathbf L\_3) = \rho(\mathbf L)
$$

If the k points form an \\( N\_1\times N\_2\times N\_3 \\) Γ-centered mesh, this holds exactly as a cyclic convolution with \\( \mathbf L \\) identified with period \\( N\_i \\)
(note that for a shifted mesh, \\( e^{i\mathbf k\cdot N\_i\mathbf a\_i} \ne 1 \\), so the convolution is not cyclic).
However, since \\( \rho \\) is built only within the range of the sparse block format (\\( \lVert\boldsymbol\tau\rVert < r\_{c,I} + r\_{c,J} \\)),
for \\( \rho \\) in an actual calculation the equality is violated by the amount of the dropped out-of-range components.

## 3. Electron density on the grid [Public] [Convention]

[Public] As in technical note 1 (around Eq. 4),

$$
n^{\sigma}(\mathbf r) = \sum\_{I}\sum\_{(J,\mathbf L)}\sum\_{\alpha\beta} \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\,
\phi\_{I\alpha}(\mathbf r - \mathbf R\_I)\\, \phi\_{J\beta}(\mathbf r - \mathbf R\_J - \mathbf L)
$$

is evaluated at each point of the real-space grid (the basis functions are real). For the grid, the treatment of periodic images, and how the basis-function values on the grid are stored, see
[Real-space grid, Poisson, exchange-correlation](grid-poisson-xc.md).

[Convention]

- In non-spin-polarized calculations, \\( n^{\uparrow} = n^{\downarrow} \\), each built from the density matrix for one spin (the total density is the sum).
- In noncollinear calculations, the four quantities \\( n^{\uparrow\uparrow} \\), \\( n^{\downarrow\downarrow} \\), \\( \mathrm{Re}\\,n^{\uparrow\downarrow} \\), \\( \mathrm{Im}\\,n^{\uparrow\downarrow} \\) are stored on the grid.
  In the exchange-correlation calculation, the 2×2 density matrix is diagonalized at each point to rotate to the local spin axis (details in the grid chapter).
- The superposition of atomic densities \\( n\_{\mathrm{atom}}(\mathbf r) \\) (the initial density and the reference for the Hartree potential) is also stored separately, in the form split in half per spin.

## 4. Charge mixing

### Notation [Design]

The input of iteration \\( m \\) is written \\( x^{\mathrm{in}}\_m \\), the output obtained from it \\( x^{\mathrm{out}}\_m \\), and the residual \\( R\_m = x^{\mathrm{out}}\_m - x^{\mathrm{in}}\_m \\).
Depending on the method, \\( x \\) is the density matrix \\( \rho(\mathbf L) \\), the Fourier components \\( n(\mathbf G) \\) of the density on the grid, the Fourier components \\( V(\mathbf G) \\) of the potential,
or the Hamiltonian \\( H(\mathbf L) \\). The history is ordered newest first as \\( m = n, n-1, \dots \\).

The Fourier components of grid quantities are defined by the unnormalized discrete Fourier transform

$$
\tilde n(\mathbf G) = \sum\_{\mathbf p} n(\mathbf r\_{\mathbf p})\\, e^{-i\mathbf G\cdot\mathbf r\_{\mathbf p}},\qquad
n(\mathbf r\_{\mathbf p}) = \frac{1}{N\_{\mathrm g}} \sum\_{\mathbf G} \tilde n(\mathbf G)\\, e^{+i\mathbf G\cdot\mathbf r\_{\mathbf p}}
$$

(\\( N\_{\mathrm g} = N\_1 N\_2 N\_3 \\) is the total number of grid points, and \\( \mathbf G \\) is a reciprocal lattice vector represented by a triple of integers from about \\( -N\_i/2 \\) to about \\( N\_i/2 \\)).
[Public] Eq. 7 of technical note 12 attaches \\( 1/N\_{\mathrm g} \\) to the forward transform, but since mixing is linear the result does not change.
Only the value of the residual norm depends on the normalization, so the [Convention] normalizations below refer to this definition.

### List of methods [Public] [Convention]

[Public] Based on the description of `scf.Mixing.Type` in manual Sec. 16.1 and Sec. 8.2. The [Convention] entries in the remarks column are facts confirmed in the existing implementation.

| Keyword value | Mixed quantity | Space | Remarks |
|---|---|---|---|
| `Simple` | density matrix | real space (matrix elements) | [Convention] default |
| `RMM-DIIS` | density matrix | real space | [Convention] the Hamiltonian is mixed in the first few iterations |
| `GR-Pulay` | density matrix | real space | [Convention] in the Pulay stage, mixing and no mixing alternate |
| `Kerker` | per-spin density \\( \tilde n^{\sigma}(\mathbf G) \\) | reciprocal space | |
| `RMM-DIISK` | per-spin density \\( \tilde n^{\sigma}(\mathbf G) \\) | reciprocal space | recommended by the manual as the first choice |
| `RMM-DIISV` | density-dependent local potential \\( \tilde V^{\sigma}(\mathbf G) \\) | reciprocal space | |
| `RMM-DIISH` | Hamiltonian matrix elements | real space | for LDA+U and constrained calculations |
| `RMM-ADIIS` | density matrix | real space | [Convention] present in the implementation but not documented in the manual ([Open]) |

[Convention]

- `RMM-DIISH` and `RMM-DIISV` cannot be combined with orbital optimization (`orbitalOpt.*`). With `RMM-DIIS`, `scf.Mixing.StartPulay` must be 4 or larger.
- The reciprocal-space methods mix the total density itself (per spin). Technical note 12 writes them in terms of the difference \\( \delta n \\) from the sum of atomic densities, but
  since the atomic densities do not change during the SCF and the mixing coefficients always sum to 1, mixing either gives the same result.
- The noncollinear reciprocal-space methods mix the three components \\( \tilde n^{\uparrow\uparrow} \\), \\( \tilde n^{\downarrow\downarrow} \\), \\( \tilde n^{\uparrow\downarrow} \\) (complex) with the same coefficients.
- \\( V^{\sigma} \\) of `RMM-DIISV` is the sum of the Hartree potential (the part due to the difference from the sum of atomic densities) and the exchange-correlation potential
  (density-independent parts such as the neutral-atom potential are not included).
- In methods that mix the density matrix, \\( n^{\mathrm{in}} \\) on the grid is built from the mixed \\( \rho^{\mathrm{in}} \\). In reciprocal-space methods,
  \\( n^{\mathrm{out}} \\) is built from the output \\( \rho^{\mathrm{out}} \\), Fourier transformed, mixed, and inverse transformed to give \\( n^{\mathrm{in}} \\) (\\( \rho \\) itself is not mixed).

### Simple mixing [Public] [Convention]

[Public] Eq. 1 of technical note 12 is simple mixing, which builds the next input as a linear combination of input and output.

[Convention] The weight \\( w \\) (`scf.*.Mixing.Weight`) of the existing implementation is **the fraction of the output**:

$$
x^{\mathrm{in}}\_{n+1} = w\\, x^{\mathrm{out}}\_n + (1 - w)\\, x^{\mathrm{in}}\_n = x^{\mathrm{in}}\_n + w R\_n
$$

Since \\( \alpha \\) in Eq. 1 of technical note 12 is the fraction of the input, read it as \\( \alpha = 1 - w \\).
Simple mixing of the density matrix uses the same \\( w \\) for all stored components (spin components, noncollinear real and imaginary parts).

- In the first iteration of the first MD step, \\( w = 1 \\) (the output is used as the next input as is).
- From the second MD or geometry-optimization step on, the initial weight of each SCF starts from `scf.Max.Mixing.Weight`.

### Residual norm [Public] [Convention]

[Public] According to manual Sec. 16.1, `NormRD` in the standard output of each iteration is the magnitude of the residual.

[Convention] The existing implementation stores the squared norm \\( \mathcal N \\) of the residual, and the displayed NormRD is its square root \\( \sqrt{\mathcal N} \\).
The definition of \\( \mathcal N \\) differs by method.

| Method / stage | \\( \mathcal N \\) |
|---|---|
| density matrix (`Simple`, `RMM-DIIS`, `GR-Pulay`) | \\( \sum\_{\text{components}}\sum\_{(I,J,\mathbf L)}\sum\_{\alpha\beta} \left(R\_{I\alpha,J\beta}(\mathbf L)\right)^2 \\) (no normalization) |
| `Kerker`, and RMM-DIISK before Pulay starts | \\( \frac{1}{N\_{\mathrm g}} \sum\_{\sigma}\sum\_{\mathbf G} K\_{0.2}(\mathbf G)\\, \lvert\tilde R^{\sigma}(\mathbf G)\rvert^2 \\) |
| Pulay stage of RMM-DIISK | \\( \frac{1}{N\_{\mathrm g} N\_s} \sum\_{\sigma}\sum\_{\mathbf G} K\_{0.1}(\mathbf G)\\, \lvert\tilde R^{\sigma}(\mathbf G)\rvert^2 \\) |
| RMM-DIISV | \\( \frac{1}{10^{5} N\_{\mathrm g} N\_s} \sum\_{\sigma}\sum\_{\mathbf G} K\_{0.1}(\mathbf G)\\, \lvert\tilde R\_V^{\sigma}(\mathbf G)\rvert^2 \\) |
| simple mixing of the Hamiltonian | \\( \frac{1}{N\_{\mathrm{atom}}} \sum (R\_H)^2 \\) |
| RMM-DIISH | \\( \frac{1}{N\_{\mathrm{atom}}} \sum m\_{I\alpha} (R\_H)^2 \\) (\\( m \\) is the metric below) |

Here \\( N\_s \\) is the number of spin components mixed in reciprocal space (1, 2, 3), \\( N\_{\mathrm{atom}} \\) is the number of atoms in the unit cell,
and \\( K\_{\gamma} \\) is the Kerker-type metric (next section). The value of the metric for the \\( G = 0 \\) component is 1 for Kerker and RMM-DIISV, and 0.5 in the Pulay stage of RMM-DIISK.
In the first iteration of each SCF, \\( \mathcal N \\) is taken to be 1 (the initial value is 100).

Since \\( \mathcal N \\) has different units and normalizations depending on the method, NormRD values cannot be compared across methods.

### Automatic weight adjustment [Convention]

The weight \\( w \\) of simple mixing (density matrix, Kerker, Hamiltonian, and iterations of RMM-DIISK and RMM-DIISV without Pulay) is
updated at every iteration according to the change of the residual norm. With the current and previous squared norms \\( \mathcal N\_0 \\), \\( \mathcal N\_1 \\),

$$
w \leftarrow
\begin{cases}
w\\, \dfrac{\mathcal N\_1}{\max(\mathcal N\_1 - \mathcal N\_0,\\, 10^{-9})} & (\mathcal N\_0 < \mathcal N\_1\text{: residual decreased})\\\\[2mm]
w\\, \dfrac{\mathcal N\_1}{\mathcal N\_1 + \mathcal N\_0} & (\mathcal N\_0 > \mathcal N\_1\text{: residual increased})\\\\[2mm]
w & (\mathcal N\_0 = \mathcal N\_1)
\end{cases}
$$

and the result is clipped to \\( [w\_{\min}, w\_{\max}] \\) (`scf.Min.Mixing.Weight`, `scf.Max.Mixing.Weight`).

- The slower the residual decreases (\\( \mathcal N\_0 \approx \mathcal N\_1 \\)), the larger the weight; if the residual increases, the weight is at least halved.
- The initial value is `scf.Init.Mixing.Weight`.
- There is also a mechanism that doubles the upper bound (up to 0.7) only when the residual keeps decreasing with the weight stuck at the upper bound for a very large number (1000) of iterations, but in practice it never acts.

Only the stage of RMM-DIISV before Pulay starts (Kerker mixing of the potential) follows a different rule, based on the direction of change of the last three band energies \\( U\_0, U\_1, U\_2 \\) and the increase or decrease of the residual:

| Direction of change of \\( U \\) | Residual | New \\( w \\) |
|---|---|---|
| same (\\( \mathrm{sgn}(U\_0 - U\_1) = \mathrm{sgn}(U\_1 - U\_2) \\)) | decreased | \\( 1.6\\, w \\) |
| same | increased | \\( w / 1.6 \\) |
| opposite | decreased | \\( 1.2\\, w \\) |
| opposite | increased | \\( w / 2 \\) |

(clipped to the same range).

### Pulay-type coefficients (RMM-DIIS) [Public] [Convention]

[Public] Eqs. 2–6 of technical note 12 (Pulay, CPL **73**, 393 (1980)): the norm of the linear combination of the last \\( p \\) residuals
\\( \bar R = \sum\_m \alpha\_m R\_m \\) is minimized under the condition \\( \sum\_m \alpha\_m = 1 \\).
With a Lagrange multiplier \\( \lambda \\),

$$
\begin{pmatrix} A & \mathbf 1 \\\\ \mathbf 1^{\mathsf T} & 0 \end{pmatrix}
\begin{pmatrix} \boldsymbol\alpha \\\\ \lambda/2 \end{pmatrix} = \begin{pmatrix} \mathbf 0 \\\\ 1 \end{pmatrix},\qquad
A\_{mm'} = \langle R\_m \vert R\_{m'} \rangle
$$

and in closed form \\( \boldsymbol\alpha = A^{-1}\mathbf 1 / (\mathbf 1^{\mathsf T}A^{-1}\mathbf 1) \\).

[Convention] With the resulting \\( \alpha \\), the existing implementation adds a residual correction to the linear combination of inputs:

$$
x^{\mathrm{in}}\_{n+1} = \sum\_{m} \alpha\_m x^{\mathrm{in}}\_m + \beta\\, \bar R
$$

Eq. 6 of technical note 12 has the form \\( \beta = 0 \\), but the existing implementation selects \\( \beta \\) from the following table according to the current squared norm \\( \mathcal N \\).

| Method | \\( \beta \\) as a function of \\( \mathcal N \\) |
|---|---|
| RMM-DIIS (density matrix), RMM-DIISH | \\( \mathcal N \ge 10^{-1} \\): 0.5, \\( \ge 10^{-2} \\): 0.6, \\( \ge 10^{-3} \\): 0.7, \\( \ge 10^{-4} \\): 0.8, otherwise: 1.0 |
| RMM-DIISK, RMM-DIISV | \\( \mathcal N \ge 10^{-2} \\): 0.5, \\( \ge 10^{-4} \\): 0.6, \\( \ge 10^{-6} \\): 0.7, \\( \ge 10^{-8} \\): 0.8, \\( \ge 10^{-10} \\): 0.9, otherwise: 1.0 |
| GR-Pulay | \\( \mathcal N \ge 10^{-2} \\): 0.2, \\( \ge 10^{-4} \\): 0.3, \\( \ge 10^{-6} \\): 0.4, \\( \ge 10^{-8} \\): 0.5, otherwise: 1.0 |

- The inner product \\( \langle R\_m\vert R\_{m'}\rangle \\) is the sum of products over all components for the density matrix, and the sum with the metric below for reciprocal space and the Hamiltonian.
- The linear system is small (dimension = number of history entries + 1), so it is solved directly. If NaN or infinity appears in the solution,
  it is replaced by \\( (\alpha\_n, \alpha\_{n-1}) = (0.1, 0.9) \\) using only the latest two (\\( (0.05, 0.95) \\) for RMM-DIISH).
- The history length is bounded by `scf.Mixing.History` (default 5). On entering the Pulay stage it starts with 2 entries and grows by one per iteration.

[Design] \\( A \\) approaches singularity as the residuals approach linear dependence. A reimplementation may add stabilization such as normalizing by the diagonal elements before solving,
or discarding old history when the smallest eigenvalue falls below a threshold (the same purpose as the replacement on NaN in the existing implementation).

### Stages of RMM-DIIS (density matrix) [Convention]

Let \\( P \\) be `scf.Mixing.StartPulay` (default 6) and \\( n \\) the SCF iteration number.

1. \\( n \le \lfloor P/2\rfloor \\): the density matrix is not mixed (the output is kept as is); instead, simple mixing of the Hamiltonian is done before solving the eigenvalue problem:
   \\( H^{\mathrm{in}}\_{n} = w\\, H[n^{\mathrm{in}}\_n] + (1-w)\\, H^{\mathrm{in}}\_{n-1} \\). The weight is adjusted automatically (the residual is \\( H[n^{\mathrm{in}}\_n] - H^{\mathrm{in}}\_{n-1} \\)).
   This stage is skipped when starting from a restart file and when orbital optimization is used.
2. \\( \lfloor P/2\rfloor < n < P \\): simple mixing of the density matrix (weight adjusted automatically).
3. \\( n = P \\): simple mixing with a fixed (unadjusted) weight. The residual history is rearranged for Pulay.
4. \\( n > P \\): Pulay mixing (starting with 2 history entries).

GR-Pulay likewise does simple mixing for \\( n < P \\) and simple mixing with a fixed weight at \\( n = P \\), and enters the alternating procedure above from \\( n > P \\).

### GR-Pulay [Convention]

This method corresponds to the "guaranteed residual reduction" type Pulay method of Bowler and Gillan (Chem. Phys. Lett. **325**, 473 (2000)).
Before Pulay starts, it is simple mixing of the density matrix (there is no Hamiltonian-mixing stage as in RMM-DIIS; the first iteration has \\( w = 1 \\)). In the Pulay stage, iterations are grouped in pairs and the following are done alternately.

- Odd: store the current output \\( \rho^{\mathrm{out}} \\), put its residual into the history, and use the output as the next input **without mixing**.
- Even: since the input is the previous output \\( \tilde\rho = \rho^{\mathrm{out}}\_{\text{prev}} \\), form the residual \\( R\_0 = \rho^{\mathrm{out}}[\tilde\rho] - \tilde\rho \\),
  find the Pulay coefficients from \\( \\{(\tilde\rho, R\_0)\\} \\) and the history \\( \\{(x^{\mathrm{in}}\_m, R\_m)\\} \\), and set \\( x^{\mathrm{in}} = \alpha\_0 \tilde\rho + \sum\_m \alpha\_m x^{\mathrm{in}}\_m + \beta\bar R \\).
  The coefficients are obtained in closed form after dividing \\( A \\) by \\( A\_{00} \\).

GR-Pulay in the existing implementation does not mix the noncollinear imaginary components (the imaginary parts of the spin-diagonal components when spin-orbit coupling etc. is present).

### Kerker mixing [Public] [Convention]

[Public] Eqs. 7–10 of technical note 12 (Kerker, PRB **23**, 3082 (1981)): simple mixing in reciprocal space with a component-dependent weight;
the weight \\( g(\mathbf G) = \lvert\mathbf G\rvert^2/(\lvert\mathbf G\rvert^2 + q\_0^2) \\) suppresses the update of long-wavelength components (the main cause of charge sloshing).
\\( q\_0 = \gamma\lvert\mathbf G\_{\min}\rvert \\), where \\( \gamma \\) is `scf.Kerker.factor`.

[Convention] The weight in the existing implementation has a form that does not vanish at long wavelength. With the length of the shortest reciprocal lattice vector \\( b\_{\min} = \min\_i\lVert\mathbf b\_i\rVert \\)
(\\( \mathbf b\_i\cdot\mathbf a\_j = 2\pi\delta\_{ij} \\)) and \\( q\_0 = \gamma\\, b\_{\min} \\), define the metric

$$
K\_{\eta}(\mathbf G) = \frac{G^2 + q\_0^2}{G^2 + (\eta q\_0)^2}\quad (\mathbf G \ne 0)
$$

Then Kerker mixing is

$$
\tilde n^{\mathrm{in}}\_{n+1}(\mathbf G) = \tilde n^{\mathrm{in}}\_n(\mathbf G) + \frac{w}{\sqrt{K\_{0.2}(\mathbf G)}}\\, \tilde R\_n(\mathbf G),\qquad K\_{0.2}(0) = 1
$$

The mixing fraction is \\( w \\) at large \\( G \\), \\( 0.2\\, w \\) as \\( G \to 0 \\), and \\( w \\) at \\( G = 0 \\).

- The weight \\( w \\) is adjusted automatically by the rule above.
- In the first iteration of the first MD step, the existing implementation **keeps the input as is** (output fraction 0; \\( K \\) is 1 in this iteration).
  As a result, the second iteration has the same Hamiltonian as the first. [Design] A reimplementation may make this iteration simple mixing with \\( w = 1 \\)
  (adopting the output as is) or Kerker mixing. The SCF history merely shifts by one iteration.
- Since the weight only drops to \\( 0.2w \\) at long wavelength, for large metallic unit cells it does not suppress sloshing as much as textbook Kerker mixing (weight 0 as \\( G \to 0 \\))
  (`tools/check_density_mixing.py`: for a Thomas–Fermi-type response \\( \chi \propto q\_{\mathrm{TF}}^2/G^2 \\),
  the amplification factor of the longest-wavelength component is \\( 1 - 0.2w(1 + q\_{\mathrm{TF}}^2/G\_{\min}^2) \\), which diverges for \\( G\_{\min} \ll q\_{\mathrm{TF}} \\) unless \\( w \\) is made sufficiently small).
  In practice the automatic weight adjustment copes by lowering \\( w \\).

### RMM-DIISK [Public] [Convention]

[Public] Eqs. 11–12 of technical note 12: the Pulay method is applied to the reciprocal-space density, with the Kerker metric \\( 1/g(\mathbf G) \\) in the inner product.
\\( 1/g = (G^2 + q\_0^2)/G^2 \\) weights long-wavelength residuals heavily.

[Convention]

- The metric is \\( K\_{0.1}(\mathbf G) \\), regularized so as not to diverge as \\( G \to 0 \\), and 0.5 at \\( G = 0 \\):
  \\( \langle R\_m\vert R\_{m'}\rangle = \sum\_{\sigma}\sum\_{\mathbf G} K\_{0.1}(\mathbf G)\\, \mathrm{Re}\left[\overline{\tilde R\_m^{\sigma}(\mathbf G)}\\, \tilde R\_{m'}^{\sigma}(\mathbf G)\right] \\).
- Before Pulay starts (\\( n < P \\)), Kerker mixing as above is used. The residuals during this period also enter the history, so the history is available from the start of the Pulay stage.
- The correction term in the Pulay stage is divided back by the metric:
  \\( \tilde n^{\mathrm{in}} = \sum\_m \alpha\_m \tilde n^{\mathrm{in}}\_m + \beta\\, \bar R' \\), where
  \\( \bar R' = \sum\_m\alpha\_m \tilde R\_m \\) if \\( \mathcal N \le 1 \\), and \\( \bar R' = \sum\_m\alpha\_m \tilde R\_m / \sqrt{K\_{0.1}} \\) if \\( \mathcal N > 1 \\) (further suppressing long wavelengths).
- With `scf.Mixing.EveryPulay` denoted \\( E \\) (default 1), in the Pulay stage Pulay mixing is done only once every \\( E \\) iterations, and in between
  Kerker-type mixing \\( \tilde n^{\mathrm{in}}\_{n+1} = \tilde n^{\mathrm{in}}\_n + w\\, \tilde R\_n / K\_{0.1}(\mathbf G) \\) (\\( 2w \\) at \\( G = 0 \\)) is done.
  Each time Pulay mixing is done, \\( w \\) is reset to \\( \max(0.1\\, w\_{\max},\\, w\_{\min}) \\) for the next Kerker-type mixing.
- The history count has a quirk: in the first iteration of the Pulay stage, \\( \min(P-1, h) \\) entries (when \\( P \le 6 \\)) or 2 entries (when \\( P > 6 \\)) are used,
  and thereafter one entry is added per iteration up to `scf.Mixing.History` (\\( h \\)).
  (For \\( P \le 6 \\) the count decreases at the second iteration, so in a reimplementation it suffices to start from 2 entries and increase regardless of \\( P \\).)

[Public] Manual Sec. 16.1 states that `scf.Mixing.EveryPulay` affects RMM-DIISK, RMM-DIISV and RMM-DIISH.
[Convention] In the existing implementation it affects only RMM-DIISK and RMM-DIISV (RMM-DIISH does Pulay mixing at every iteration in the Pulay stage).

### RMM-DIISV [Convention]

The same procedure as RMM-DIISK is applied to the potential \\( \tilde V^{\sigma}(\mathbf G) \\) instead of the density. It is mixed right before the Hamiltonian is assembled.

- The metric is \\( K\_{0.1} \\), but 1 at \\( G = 0 \\). The normalization of the squared norm is as in the table above (divided by \\( 10^5 \\)).
- Kerker-type mixing before Pulay starts is \\( \tilde V^{\mathrm{in}}\_{n+1} = \tilde V^{\mathrm{in}}\_n + (w/\sqrt{K\_{0.1}})\\, \tilde R\_V \\), and the weight is adjusted by the separate rule (above) using the trend of the band energy.
- In the first iteration the new potential is used as is.
- The \\( \beta \\) table of the Pulay stage, `scf.Mixing.EveryPulay`, and the dividing back of the correction term are the same as for RMM-DIISK.

### RMM-DIISH [Convention]

The Hamiltonian matrix elements \\( H(\mathbf L) \\) are mixed right before solving the eigenvalue problem.
Since LDA+U and constraint terms enter the Hamiltonian through the density matrix (occupation matrix), this is more stable than mixing only the density (recommended in manual Sec. 16.1).

- The input is the Hamiltonian \\( H^{\mathrm{in}}\_{n-1} \\) mixed in the previous iteration, the output is \\( H[n^{\mathrm{in}}\_n] \\) assembled from the current input density, and the residual is their difference.
- Before Pulay starts (\\( n < P \\)), simple mixing of the Hamiltonian (weight adjusted automatically). The first iteration has \\( w = 1 \\).
- The inner product in the Pulay stage is the sum multiplied by the metric per orbital \\( \alpha \\) of atom \\( I \\) (row index)

$$
m\_{I\alpha} = \frac{5}{d^2 + 5},\qquad d = \left\lvert H^{\uparrow\uparrow}\_{I\alpha,I\alpha}(\mathbf 0) - \mu \right\rvert\ \ (\text{Hartree})
$$

  (\\( H \\) is the most recent input, \\( \mu \\) the chemical potential). Orbitals with energies farther from the Fermi level are weighted less.
- In noncollinear calculations, \\( H^{\uparrow\uparrow} \\), \\( H^{\downarrow\downarrow} \\), \\( \mathrm{Re}\\,H^{\uparrow\downarrow} \\), \\( \mathrm{Im}\\,H^{\uparrow\downarrow} \\), and
  all imaginary-part matrices due to spin-orbit coupling etc. are mixed with the same coefficients.

### Automatic determination of the Kerker factor [Public] [Convention]

[Public] Manual Sec. 16.2: if `scf.Kerker.factor` is not specified, it is determined automatically from the mean and anisotropy of the lengths of the reciprocal lattice vectors (Eqs. 21–23 of the same section).
The manual's formula appears to form the mean and differences with squared lengths of the reciprocal lattice vectors and multiply by a factor of 4, but the fraction structure as displayed is hard to read, and the details do not match the implementation below.

[Convention] The existing implementation (Kerker, RMM-DIISK, RMM-DIISV, when not specified) uses the lengths of the reciprocal lattice vectors \\( b\_i = \lVert\mathbf b\_i\rVert \\) to set

$$
\bar b = \frac{b\_1 + b\_2 + b\_3}{3},\qquad
\Delta b = \frac{\lvert b\_2 - b\_1\rvert + \lvert b\_3 - b\_2\rvert + \lvert b\_1 - b\_3\rvert}{3},\qquad
\gamma = \frac{0.5}{b\_{\min}}\left(\frac{2\Delta b}{\bar b} + 1\right)
$$

Hence \\( q\_0 = \gamma b\_{\min} = 0.5\\,(2\Delta b/\bar b + 1) \\) bohr⁻¹, independent of the unit-cell size;
it is 0.5 bohr⁻¹ for an isotropic unit cell, approaches 1.5 bohr⁻¹ for a unit cell extremely long in one direction (an elongated unit cell, with only one short reciprocal lattice vector),
and is at most 2.5 bohr⁻¹ for a unit cell extremely short in one direction (with only one long reciprocal lattice vector)
(example: 1.32 bohr⁻¹ for a 5×5×40 bohr rectangular box).
When specified, \\( q\_0 = \gamma\\, b\_{\min} \\), which depends on the unit-cell size through \\( b\_{\min} \\).

### Parameters and defaults [Public] [Convention]

[Public] defaults are those documented in manual Sec. 8.2; [Convention] are values confirmed in the existing implementation (the two agree).

| Keyword | Default | Meaning |
|---|---|---|
| `scf.Mixing.Type` | `Simple` [Convention] | method |
| `scf.Init.Mixing.Weight` | 0.3 [Public] | initial weight of simple mixing (fraction of the output) |
| `scf.Min.Mixing.Weight` | 0.001 [Public] | lower bound of the weight |
| `scf.Max.Mixing.Weight` | 0.4 [Public] | upper bound of the weight; initial value from the second MD/optimization step on |
| `scf.Kerker.factor` | automatic [Public] | \\( \gamma \\); automatic determination above |
| `scf.Mixing.History` | 5 [Convention] | upper bound of the Pulay-type history; the manual recommends around 30 |
| `scf.Mixing.StartPulay` | 6 [Public] | iteration at which Pulay mixing starts |
| `scf.Mixing.EveryPulay` | 1 [Public] | interval of Pulay steps in RMM-DIISK and RMM-DIISV |
| `scf.Mixing.Control.Temp` | off [Convention] | automatic control of the electronic temperature |
| `scf.criterion` | \\( 10^{-6} \\) Hartree [Public] | threshold of the convergence criterion |
| `scf.maxIter` | 40 [Public] | maximum number of iterations |
| `scf.ElectronicTemperature` | 300 K [Public] | electronic temperature |
| `scf.system.charge` | 0 [Public] | extra charge |
| `scf.ExtCharge.History` | 3 [Convention] | number of history entries used for density extrapolation |

### Adjustment during the SCF [Public] [Convention]

[Public] Manual Sec. 16.3: if, during the calculation, a file named `System.Name` with `_SCF_keywords` appended is placed in the working directory,
the file is read at the end of each SCF iteration and `scf.maxIter`, `scf.Min.Mixing.Weight`, `scf.Max.Mixing.Weight`, `scf.Kerker.factor`,
and `scf.Mixing.StartPulay` are updated. Setting `scf.maxIter` to a negative value stops the calculation.

[Convention] The existing implementation also reads `scf.criterion` from this file. The file remains after being read and is reread at every iteration.

[Convention] From the second MD or geometry-optimization step on, if the run could be restarted from the density of the previous step,
then when \\( \mathcal N < 0.01 \\) is reached in an iteration before `scf.Mixing.StartPulay`, Pulay mixing starts from the next iteration (only during the SCF of that step).

## 5. Convergence criterion [Public] [Convention]

[Public] Manual Sec. 8.2 and Sec. 16.1: `scf.criterion` (in Hartree, default \\( 10^{-6} \\)) is the threshold for the change of the band energy between iterations,
\\( \mathrm{d}U\_{\mathrm{ele}} \\). After `scf.maxIter` (default 40) iterations, the SCF is terminated even if not converged.

[Convention] The criterion of the existing implementation looks not only at the change of the band energy but also at the residual norm. With threshold \\( c \\),
convergence is declared when all of

$$
n \ge 2,\qquad \left\lvert U\_{\mathrm{ele},n} - U\_{\mathrm{ele},n-1}\right\rvert < c,\qquad \mathcal N < 10\\, c
$$

are satisfied. \\( \mathcal N \\) is the most recent squared norm (the square of NormRD, with the per-method normalization in the table above);
note that it compares a quantity of a different dimension with \\( c \\), which has the dimension of energy.

- When orbital optimization is used, the residual-norm condition is not used.
- A terminated run is treated as unconverged, and in the next MD step the density of the previous step is not reused if \\( \mathcal N > 0.3 \\).

[Design] A reimplementation uses the same criterion (both energy and residual norm) as the default, with the residual-norm normalization per method matching the table above.
If the criterion is changed, tests should confirm that the tolerance on the total energy (about \\( 10^{-6} \\) Hartree) is met.

## 6. Initial density and restart [Public] [Convention]

[Convention] In the first SCF without a restart file, the input density is the superposition of atomic densities.

- The density of atom \\( I \\) is the valence density of the pseudo-atom (spherically symmetric, integrating to \\( Z\_I^{\mathrm{val}} \\)) scaled by \\( (N\_I^{\uparrow} + N\_I^{\downarrow})/Z\_I^{\mathrm{val}} \\) using
  the up and down electron numbers \\( N\_I^{\uparrow} \\), \\( N\_I^{\downarrow} \\) given in the atomic-coordinate field of the input,
  with \\( N\_I^{\uparrow}/(N\_I^{\uparrow}+N\_I^{\downarrow}) \\) assigned to up spin and the rest to down spin (half each in the unpolarized case).
  Hence the number of electrons of the initial density is the input \\( \sum\_I (N\_I^{\uparrow}+N\_I^{\downarrow}) \\), and `scf.system.charge` is not used at this stage
  (the correct number of electrons is obtained from the first eigenvalue problem).
- The sum of atomic densities \\( n\_{\mathrm{atom}} \\) used as the reference for the Hartree potential is the sum of neutral-atom densities without this scaling.
- In noncollinear calculations, the 2×2 density \\( n^{\sigma\sigma'} \\) is built along the spin direction (angles) of each atom given in the input:
  from the total density \\( n \\) and magnetization \\( \mathbf m \\) of each atom, \\( n^{\uparrow\uparrow} = (n + m\_z)/2 \\), \\( n^{\downarrow\downarrow} = (n - m\_z)/2 \\),
  \\( \mathrm{Re}\\,n^{\uparrow\downarrow} = m\_x/2 \\), \\( \mathrm{Im}\\,n^{\uparrow\downarrow} = -m\_y/2 \\).
- Since the first iteration builds the Hamiltonian from this density, no density-matrix input is needed.

[Public] Manual Sec. 17.1: with `scf.restart on`, the SCF is restarted from the results of a previous calculation saved in the `System.Name_rst` directory.

[Convention] The restart uses the difference of the grid density per spin component from the atomic density, \\( \Delta n^{\sigma} = n^{\sigma} - n\_{\mathrm{atom}}/2 \\), and the Hamiltonian matrix elements.

- The new input density is \\( n^{\sigma} = n\_{\mathrm{atom}}/2 + \Delta n^{\sigma} \\) (\\( n\_{\mathrm{atom}} \\) rebuilt at the current atomic positions).
  If the number of grid points differs from that at saving time, no restart is done.
- Even when restarted, the saved Hamiltonian is normally not used; the Hamiltonian is rebuilt from the loaded density
  (when spin-orbit coupling, LDA+U, constraints, or a Zeeman term are present, the saved Hamiltonian is used as the Hamiltonian of the first iteration).
- In MD and geometry optimization, from the second step on the run automatically restarts from the result of the previous step (disabled by `scf.restart alloff`).

## 7. Density extrapolation in MD and geometry optimization [Public] [Convention]

[Public] Technical note 8: regarding the atomic coordinates as a single vector \\( \mathbf x \\) (of dimension \\( 3N\_{\mathrm{atom}} \\)),
find the coefficients that best approximate the new coordinates by a linear combination \\( \sum\_m \alpha\_m \mathbf x\_m \\) of the coordinates of the past \\( p \\) steps (Eqs. 1–3; least squares under the condition \\( \sum\_m\alpha\_m = 1 \\)),
and extrapolate the density with the same coefficients (Eq. 4):

$$
n^{\mathrm{in}}\_{n+1} = n\_{\mathrm{atom}}[\mathbf x\_{n+1}] + \sum\_{m=n-p+1}^{n} \alpha\_m\\, \Delta n^{\mathrm{out}}\_m,\qquad
\Delta n^{\mathrm{out}}\_m = n^{\mathrm{out}}\_m - n\_{\mathrm{atom}}[\mathbf x\_m]
$$

The technical note states that \\( p = 3 \\) is optimal in many cases.
The coefficients are obtained from the stationarity condition of \\( \lVert \sum\_m\alpha\_m\mathbf x\_m - \mathbf x\_{n+1}\rVert^2 - \lambda(1 - \sum\_m\alpha\_m) \\),

$$
\begin{pmatrix} X & \mathbf 1 \\\\ \mathbf 1^{\mathsf T} & 0 \end{pmatrix}
\begin{pmatrix} \boldsymbol\alpha \\\\ \lambda/2 \end{pmatrix} = \begin{pmatrix} \mathbf b \\\\ 1 \end{pmatrix},\qquad
X\_{mm'} = \mathbf x\_m\cdot\mathbf x\_{m'},\quad b\_m = \mathbf x\_m\cdot\mathbf x\_{n+1}
$$

(the display of the right-hand side of Eq. 3 in the technical note is garbled, but from the stationarity condition the right-hand side is \\( \mathbf b \\)).

[Convention] The existing implementation determines the coefficients as follows; unlike the technical note, it is an **unconstrained** least-squares fit.

- The number of history entries is `scf.ExtCharge.History` (\\( p \\), default 3). Extrapolation is done only when the MD step number exceeds \\( p \\) and \\( p > 1 \\).
  Until then, \\( (\alpha\_n, \alpha\_{n-1}, \dots) = (1, 0, \dots) \\) (the difference density of the previous step is used as is).
- The normal equation \\( X\boldsymbol\alpha = \mathbf b \\) is solved via the eigendecomposition of \\( X \\) (a pseudo-inverse discarding directions with eigenvalue magnitude below \\( 10^{-15} \\)).
- If the magnitude of any resulting coefficient exceeds 10, or \\( \lvert\sum\_m\alpha\_m - 1\rvert > 10^{-3} \\), extrapolation is abandoned and \\( (1, 0, \dots) \\) is restored.
  (Since the coordinates are large in magnitude and the displacements small, the sum is usually close to 1 even without the constraint.)
- Extrapolation is applied to the density on the grid (per spin component); the Hamiltonian and density matrix are not extrapolated.
- When an atom leaves the unit cell and its fractional coordinates are wrapped back into \\( [0, 1) \\), the coordinate history is shifted by the same amount so that the extrapolation is continuous.

[Design] A reimplementation uses the constrained formula of the technical note, with the same safeguards as the existing implementation (an upper bound on the coefficient magnitudes and \\( (1, 0, \dots) \\) on failure).
For a trajectory with uniform acceleration (coordinates quadratic in time), both methods give \\( \boldsymbol\alpha = (3, -3, 1) \\), and the new coordinates are reproduced exactly.

## Requirements for implementations [Design]

- Occupations and chemical potential are pure functions whose only inputs are the array of eigenvalues (k point, spin, state), weights, number of electrons, and temperature.
- Density-matrix construction receives the structure of the sparse block format (the set of atom pairs) and fills in only the values. \\( \rho \\) and \\( E \\) are built at the same time.
- Mixing abstracts "the vector space of the mixed quantity" (inner product, metric, linear combination) so that each method (Simple / Pulay / Kerker) can be applied
  uniformly to the density matrix, the reciprocal-space density, the potential, and the Hamiltonian.
  The history length, the \\( \beta \\) table, the weight-adjustment rule, and the Pulay starting iteration are held as configuration values.
- The default rules can be switched between reproducing the SCF history of the existing implementation (NormRD at each iteration) and prioritizing convergence speed.

## Tests

`tools/check_density_mixing.py` checks the analytic, invariant and linear-model items below.

| Kind | Content | Tolerance (guideline) |
|---|---|---|
| analytic | two-level system (\\( \varepsilon\_1 = -0.3 \\), \\( \varepsilon\_2 = 0.1 \\) Hartree, unpolarized): compare \\( \mu \\) for \\( N\_{\mathrm e} = 3 \\) with the value obtained by solving \\( 2(f\_1 + f\_2) = 3 \\) with a different root finder | \\( \mu \\): \\( 10^{-9} \\) Hartree |
| analytic | same system with \\( N\_{\mathrm e} = 2 \\): the electron number matches and \\( \mu \\) lies in the gap (the deviation from the center depends on temperature and stopping condition) | electron number \\( 10^{-11} \\) |
| analytic | bisection stopping condition \\( \lvert N(\mu) - N\_{\mathrm e}\rvert < 10^{-12} \\) | electron number \\( 10^{-11} \\) |
| invariant | electron-number sum rule \\( g\_s\sum\rho S = N\_{\mathrm e} \\) (1D two-orbital chain with overlap, metallic temperature and low temperature) | \\( 10^{-12} \\) |
| invariant | band-energy sum rule \\( g\_s\sum\rho H = g\_s\sum E S = U\_{\mathrm{ele}} \\) | \\( 10^{-12} \\) |
| invariant | Hermiticity \\( \rho\_{J\beta,I\alpha}(-\mathbf L) = \rho\_{I\alpha,J\beta}(\mathbf L) \\) (Cluster and Band; does not hold for O(N) solvers; in noncollinear calculations with complex conjugation and spin exchange) | rounding error |
| invariant | idempotency \\( \rho S \rho = \rho \\) for an insulator at low temperature (cyclic convolution on a Γ-centered mesh) | \\( 10^{-10} \\) |
| invariant | in a noncollinear system with \\( \rho^{\uparrow\downarrow} \\) set to 0 (two copies of a collinear Hamiltonian), the density matrix agrees with the collinear calculation | \\( 10^{-12} \\) |
| invariant | merging \\( \pm\mathbf k \\) into one point with weight 2 by time reversal does not change the density matrix | \\( 10^{-12} \\) |
| linear model | \\( x^{\mathrm{out}} - x^\ast = J(x^{\mathrm{in}} - x^\ast) \\) (with specified eigenvalues \\( \lambda \\) of \\( J \\)): simple mixing converges for \\( w < 2/(1 - \lambda\_{\min}) \\) and diverges beyond | convergence/divergence decision |
| linear model | in the same model, the Pulay method with a sufficiently long history converges within dimension + a few iterations | number of iterations |
| linear model | Pulay coefficients: the solution of the bordered linear system agrees with the closed form \\( A^{-1}\mathbf 1/(\mathbf 1^{\mathsf T}A^{-1}\mathbf 1) \\), and \\( \sum\alpha = 1 \\) | \\( 10^{-8} \\) |
| linear model | Thomas–Fermi-type response \\( \tilde R \propto -(q\_{\mathrm{TF}}^2/G^2) \\): textbook Kerker mixing (\\( q\_0 = q\_{\mathrm{TF}} \\), \\( w = 1 \\)) converges in one iteration, and simple mixing (\\( w = 0.3 \\)) diverges | number of iterations |
| analytic | extrapolation coefficients: when the coordinates move quadratically in time, both constrained and unconstrained give \\( (3, -3, 1) \\) | \\( 10^{-8} \\) |
| analytic | automatic Kerker factor: \\( q\_0 = 0.5 \\) bohr⁻¹ for a cubic lattice, and the value of the formula above for tetragonal and orthorhombic lattices | \\( 10^{-12} \\) |
| cross-check | fixtures from the existing implementation: for small systems (H₂O molecule, Si crystal, Fe crystal collinear and noncollinear), \\( \mu \\), \\( \rho(\mathbf L) \\), \\( E(\mathbf L) \\) given the eigenvalues and eigenvectors | \\( \mu \\): \\( 10^{-8} \\) Hartree (excluding the indeterminacy within the gap), \\( \rho \\), \\( E \\): \\( 10^{-10} \\) |
| cross-check | with the same fixtures, the result of one mixing step and the residual norm given pairs of input/output density matrices (or densities) and the history | \\( 10^{-10} \\) (relative) |
| regression | converge the SCF of small systems with each method; the total energy agrees regardless of method, and the number of iterations is comparable to the existing implementation | total energy \\( 10^{-6} \\) Hartree; number of iterations as a guideline |

## [Open]

- The explicit formulas of `RMM-ADIIS` (DIIS using the augmented Roothaan–Hall energy function) and its treatment in the manual.
- How the chemical potential is determined and the density matrix built in the O(N) methods (DC, Krylov) and NEGF is covered in [O(N) methods](on-methods.md) and [Electron transport with NEGF](negf.md) (this chapter targets Cluster and Band).
- For the noncollinear density on the grid, the sign correspondence between \\( \mathrm{Im}\\,n^{\uparrow\downarrow} \\) and the magnetization \\( m\_y \\) (in the initial density, \\( \mathrm{Im}\\,n^{\uparrow\downarrow} = -m\_y/2 \\)).
  Whether this agrees with the derivation from the definition of the density matrix needs to be checked together with the grid chapter.
- The formats of the restart file and the extrapolation history file (may be newly defined in this note; candidates for moving to [Design]).
- Treatment of inputs that force the occupation of specific states to 0, such as `empty.occupation.*` and `empty.states.*` (used e.g. for core-hole calculations).
