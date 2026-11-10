# DFT+U, constrained DFT, and magnetic features

> Status: first draft (awaiting review). The equations rewrite technical note 6 (LDA+U, Ver. 1.0; hereafter TN6), technical note 2 (noncollinear; hereafter TN2),
> Han–Ozaki–Yu (PRB **74**, 045110 (2006)), Ryee–Han (Sci. Rep. **8**, 9559 (2018); JPCM **30**, 275802 (2018)),
> and manual Secs. 37–43 in the notation of this note; conventions not found in the public documents were confirmed against the existing implementation.
> The consistency of the equations (energy and potential, forces, invariants of the Coulomb tensor, closed form of the constraint term) is
> checked numerically with `tools/check_dftu_magnetism.py`.

## Scope of this chapter

This chapter collects the features that add terms to the Hamiltonian through the occupation matrices of localized orbitals on atoms.

- DFT+U: the simplified form (Dudarev type, `scf.DFTU.Type 1`) and the general form (Liechtenstein type, `scf.DFTU.Type 2`),
  three definitions of the occupation matrix (onsite, full, dual), four double-counting corrections, the energy `Uhub`, the potential, and forces.
- Enhancement of orbital polarization (the orbital-polarization flag in `Atoms.SpeciesAndCoordinates`).
- Constraint on the direction (and magnitude) of noncollinear spins (`scf.Constraint.NC.Spin`, energy `Ucs`).
- Zeeman terms for spin and orbital magnetic moments (`scf.NC.Zeeman.Spin`, `scf.NC.Zeeman.Orbital`; energies `Uzs`, `Uzo`).
- Evaluation of orbital magnetic moments, and the procedure for the magnetic anisotropy energy (MAE) by the second-variation method.
- Input of the initial magnetic configuration (spin, angle, and flag columns of the coordinate lines).

The noncollinear density and XC on the grid are in [Grid, Poisson, and exchange-correlation](grid-poisson-xc.md), the density matrix and restart in
[Density matrix and charge mixing](density-mixing.md), and the sum of the total energy and the force framework in [Total energy and forces](energy-force.md).
The matrix elements of spin-orbit coupling themselves are not treated here.

## Notation [Design]

In addition to the symbols of [Notation](../../../conventions/notation.md), [Sparse block format](../design/sparse-block-format.md), and [Total energy and forces](energy-force.md), the following are used.

| Symbol | Meaning |
|---|---|
| \\( \rho^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) \\) | Real-space density matrix. \\( \sum\_{\mathbf k} w\_{\mathbf k}\sum\_\nu f\_\nu e^{i\mathbf k\cdot\mathbf L}\\, \overline{c\_{\sigma I\alpha,\nu}}\\, c\_{\sigma' J\beta,\nu} \\) ([Density matrix](density-mixing.md)) |
| \\( s = (I, l, \zeta) \\) | Shell: the set of \\( 2l+1 \\) orbitals of atom \\( I \\) with angular momentum \\( l \\) and multiplicity \\( \zeta \\) |
| \\( a, b, c, d \\) | Index of real spherical harmonics within a shell (ordering of the [PAO basis](pao-basis.md), 0 to \\( 2l \\)) |
| \\( n^{\sigma\sigma'}\_{s,ab} \\) | Occupation matrix of shell \\( s \\) (defined below) |
| \\( N^{\sigma\sigma'}\_s = \sum\_a n^{\sigma\sigma'}\_{s,aa} \\), \\( N\_s = N^{\uparrow\uparrow}\_s + N^{\downarrow\downarrow}\_s \\) | Occupation number of the shell (per spin component, total) |
| \\( v^{\sigma\sigma'}\_{s,ab} \\) | Potential matrix of the shell (\\( = \partial E / \partial n^{\sigma\sigma'}\_{s,ab} \\)) |
| \\( U\_s, J\_s \\) | U and J of the shell (input in eV, Hartree internally) |
| \\( U\_{abcd} = \langle ab \vert V \vert cd\rangle \\) | Coulomb tensor in real spherical harmonics (\\( a, c \\) are electron 1, \\( b, d \\) are electron 2) |
| \\( \boldsymbol\sigma = (\sigma\_x, \sigma\_y, \sigma\_z) \\) | Pauli matrices |
| \\( T^{(I)} \\) | 2×2 spin density matrix of atom \\( I \\) (defined below) |
| \\( \hat{\mathbf e}(\theta,\varphi) = (\sin\theta\cos\varphi, \sin\theta\sin\varphi, \cos\theta) \\) | Direction vector |

**Orientation of matrix elements and derivatives** [Design]. Since the density matrix is defined as \\( \rho\_{ab} = \sum f\\, \overline{c\_a}\\, c\_b \\), the band energy is
\\( \sum\_{ab} H\_{ab}\\, \rho\_{ab} \\) (sum of elementwise products without transposition). In this chapter, for an energy \\( E(n) \\) the potential is defined as
\\( v\_{ab} = \partial E/\partial n\_{ab} \\) (partial derivative treating each component of \\( n \\) as independent), and \\( v \\) is added to the Hamiltonian as is (without transposition).
In collinear calculations \\( n \\) is real symmetric, so transposition cannot be distinguished, but in noncollinear calculations the sign is reversed unless this orientation is respected
(e.g. \\( E = \tfrac12 U\\,[\mathrm{Tr}\\,n - \mathrm{Tr}\\,n^2] \\) gives \\( v = U(\tfrac12 I - n^{\mathsf T}) \\), the same as the equation after TN2 Eq. 115).

**Real-space matrix product** [Design]. The product of matrices \\( A, B \\) in atom-pair blocks is written as
\\( (AB)\_{I\alpha,J\beta}(\mathbf L) = \sum\_{K,\mathbf L'}\sum\_\gamma A\_{I\alpha,K\gamma}(\mathbf L')\\, B\_{K\gamma,J\beta}(\mathbf L-\mathbf L') \\).
This corresponds to the k-space product \\( A(\mathbf k)B(\mathbf k) \\) (section on the sum rule in [Density matrix](density-mixing.md)). In the noncollinear case each spin block is multiplied by
the spin-independent \\( S \\). \\( (\cdot)\_{II} \\) denotes the diagonal block of atom \\( I \\) at \\( \mathbf L = 0 \\), and \\( (\cdot)\_s \\) further denotes the part belonging to shell \\( s \\).

## 1. Occupation matrix [Public] [Convention]

[Public] TN6 Eqs. 5–16, TN2 Eqs. 103–114: the occupation matrix is the expectation value of a projection operator,

$$
n^{\sigma\sigma'}\_{s,ab} = \sum\_\nu f\_\nu \langle \psi\_\nu \vert \hat n^{\sigma\sigma'}\_{s,ab} \vert \psi\_\nu\rangle
$$

and there are three kinds depending on the choice of operator. Using the dual orbitals \\( \vert\tilde\phi\rangle = \sum S^{-1}\vert\phi\rangle \\),
onsite corresponds to \\( \vert\tilde\phi\_a\rangle\langle\tilde\phi\_b\vert \\), full to \\( \vert\phi\_a\rangle\langle\phi\_b\vert \\),
and dual to \\( \tfrac12(\vert\tilde\phi\_a\rangle\langle\phi\_b\vert + \vert\phi\_a\rangle\langle\tilde\phi\_b\vert) \\).
In terms of the density matrix and the overlap matrix (with the index orientation of this note):

| `scf.Hubbard.Occupation` | Occupation matrix | Sum rule \\( \sum\_I \mathrm{Tr}\\, n\_{II} = N\_{\mathrm e} \\) |
|---|---|---|
| `onsite` | \\( n\_{II} = \rho\_{II} \\) (the diagonal block at \\( \mathbf L = 0 \\) itself) | Does not hold |
| `full` | \\( n\_{II} = (S\rho S)\_{II} \\) | Does not hold |
| `dual` (default) | \\( n\_{II} = \tfrac12(\rho S + S\rho)\_{II} \\) | Holds |

The occupation matrix \\( n\_s \\) of shell \\( s \\) is its \\( (2l+1)\times(2l+1) \\) part (in the noncollinear case \\( 2(2l+1) \\)-dimensional, including the 2×2 spin structure).
In the collinear case it is defined per spin (\\( \sigma = \sigma' \\)).

[Convention] Index orientation: the dual of the existing implementation is \\( n^{\sigma\sigma'}\_{ab} = \tfrac12[(\rho^{\sigma\sigma'} S)\_{ab} + (S\rho^{\sigma\sigma'})\_{ab}] \\), the same as TN2 Eq. 114.
TN6 Eq. 16 is written in the transposed form, but TN6 treats only the collinear (real symmetric) case, so there is no difference.
In the noncollinear case, the dual and full \\( n \\) are Hermitian as \\( 2(2l+1) \\)-dimensional matrices (\\( n^{\downarrow\uparrow} = (n^{\uparrow\downarrow})^\dagger \\)).
For noncollinear onsite, \\( n^{\uparrow\downarrow}\_{ab} = \rho^{\uparrow\downarrow}\_{Ia,Ib}(0) \\) and \\( n^{\downarrow\uparrow}\_{ab} = \overline{\rho^{\uparrow\downarrow}\_{Ib,Ia}(0)} \\).

[Convention] Others:

- In non-spin-polarized calculations \\( \rho \\) is for one spin ([Density matrix](density-mixing.md)), and so is \\( n \\). The energy is doubled (see below).
- When orbital optimization (`orbitalOpt.Method`) is used, U acts on the uncontracted (primitive) orbitals. The occupation matrix is built after transforming back to the density matrix of the primitive orbitals with the contraction coefficients (TN6 Sec. 6).
- `full` cannot be combined with stress calculation (`scf.stress.tensor`); the input check stops. Noncollinear calculations also cannot be combined with stress ([Total energy and forces](energy-force.md)).

### Position within the SCF [Convention]

- The occupation matrix and the potential \\( v \\) are built from the (input-side) density matrix **after** charge mixing. This \\( v \\) enters the Hamiltonian of the next iteration,
  and `Uhub`, `Ucs`, `Uzs` output in that iteration are evaluated with the same \\( n \\).
- After convergence, the occupation matrix is rebuilt from the final density matrix and written to the file `System.Name.DM_onsite`
  (occupation matrices per atom and spin, their eigenvalues and eigenvectors; the heading is "Occupation Number in LDA+U and Constraint DFT").
- If any of LDA+U, the constraint, or a Zeeman term is present, the first iteration on restart uses the saved Hamiltonian ([Density matrix](density-mixing.md) Sec. 6).

## 2. Simplified form (`scf.DFTU.Type 1`) [Public] [Convention]

[Public] TN6 Eqs. 2–4, TN2 Eq. 101 (rotationally invariant form of Dudarev et al.):

$$
E\_U = \frac12 \sum\_s U\_s \sum\_\sigma \left[ \mathrm{Tr}\\, n^{\sigma\sigma}\_s - \mathrm{Tr}\\, (n^{\sigma\sigma}\_s n^{\sigma\sigma}\_s) \right]
\quad\text{(collinear)},\qquad
E\_U = \frac12 \sum\_s U\_s \left[ \mathrm{Tr}\\, n\_s - \mathrm{Tr}\\, (n\_s n\_s) \right]
\quad\text{(noncollinear, }2(2l+1)\text{-dimensional)}
$$

$$
v^{\sigma\sigma'}\_{s,ab} = U\_s \left( \tfrac12 \delta\_{\sigma\sigma'}\delta\_{ab} - n^{\sigma'\sigma}\_{s,ba} \right)
$$

- The sum is only within each shell; there are no terms between different shells (different \\( l \\) or \\( \zeta \\)). The input \\( U \\) is used in the sense of the Dudarev form
  \\( U\_{\mathrm{eff}} = U - J \\) (J is not read).
- [Convention] For non-spin-polarized calculations, `Uhub` is twice the one-spin expression (the potential remains the one-spin expression).
  [Public] Manual Sec. 40.1 requires `on` or `nc` for `scf.SpinPolarization` with DFT+U, but the existing implementation runs non-polarized calculations with the simplified form without stopping.
- For idempotent \\( n\_s \\) (eigenvalues 0 or 1), \\( E\_U = 0 \\).

## 3. General form (`scf.DFTU.Type 2`) [Public] [Convention]

### 3.1 Energy and potential

[Public] The form of Liechtenstein et al., with the extension by Ryee–Han (Sci. Rep. 2018, JPCM 2018) to double counting and the noncollinear case.
For each shell, the Hartree–Fock-type interaction

$$
E\_{\mathrm{HF}}[n] = \frac12 \sum\_{abcd} U\_{abcd}\\, \mathcal N\_{ac}\\, \mathcal N\_{bd}
{}- \frac12 \sum\_{abcd} U\_{abdc} \sum\_{\sigma\sigma'} n^{\sigma\sigma'}\_{ac}\\, n^{\sigma'\sigma}\_{bd},\qquad
\mathcal N = n^{\uparrow\uparrow} + n^{\downarrow\downarrow}
$$

is formed (taking the real part of the whole), and the double-counting correction \\( E\_{\mathrm{dc}} \\) is subtracted. In the collinear case, setting \\( n^{\uparrow\downarrow} = 0 \\) gives

$$
E\_{\mathrm{HF}} = \frac12\sum U\_{abcd}\left(n^{\uparrow}\_{ac} n^{\downarrow}\_{bd} + n^{\downarrow}\_{ac} n^{\uparrow}\_{bd}\right)
{}+ \frac12\sum \left(U\_{abcd} - U\_{abdc}\right)\sum\_\sigma n^{\sigma}\_{ac} n^{\sigma}\_{bd}
$$

which, [Convention], is the same as the collinear expression of the existing implementation. The potential is

$$
v^{\sigma\sigma'}\_{ac} = \delta\_{\sigma\sigma'} \sum\_{bd} U\_{abcd}\\, \mathcal N\_{bd} - \sum\_{bd} U\_{abdc}\\, n^{\sigma'\sigma}\_{bd} - v^{\mathrm{dc},\sigma\sigma'}\_{ac}
$$

(the derivative of \\( E \\) using \\( U\_{abcd} = U\_{badc} \\)).

### 3.2 Double-counting correction (`scf.dc.Type`) [Public] [Convention]

| `scf.dc.Type` | \\( E\_{\mathrm{dc}} \\) | Potential correction |
|---|---|---|
| `sFLL` (default) | \\( \tfrac12 U N(N-1) - \tfrac12 J \left[ \sum\_{\sigma\sigma'} N^{\sigma\sigma'} N^{\sigma'\sigma} - N \right] \\) | \\( v^{\mathrm{dc},\sigma\sigma'}\_{ac} = \delta\_{ac}\left[ U(N-\tfrac12)\delta\_{\sigma\sigma'} - J(N^{\sigma'\sigma} - \tfrac12\delta\_{\sigma\sigma'}) \right] \\) |
| `cFLL` | \\( \tfrac12 U N(N-1) - \tfrac14 J N(N-2) \\) | \\( v^{\mathrm{dc},\sigma\sigma'}\_{ac} = \delta\_{ac}\delta\_{\sigma\sigma'}\left[ U(N-\tfrac12) - J(\tfrac12 N - \tfrac12) \right] \\) |
| `sAMF` | None (\\( E = E\_{\mathrm{HF}}[\tilde n] \\)) | None (\\( v = v\_{\mathrm{HF}}[\tilde n] \\)) |
| `cAMF` | None (\\( E = E\_{\mathrm{HF}}[\tilde n] \\)) | None (\\( v = v\_{\mathrm{HF}}[\tilde n] \\)) |

Here \\( N, N^{\sigma\sigma'} \\) are shell values. AMF (around mean field) uses \\( \tilde n \\), the occupation matrix with its average subtracted.

- `sAMF`: \\( \tilde n^{\sigma\sigma'} = n^{\sigma\sigma'} - \dfrac{N^{\sigma\sigma'}}{2l+1} I \\) (in the noncollinear case the same form is subtracted from the spin-off-diagonal components as well).
- `cAMF`: \\( \tilde n^{\sigma\sigma} = n^{\sigma\sigma} - \dfrac{N}{2(2l+1)} I \\); spin-off-diagonal components are left unchanged.

In the collinear case \\( \sum\_{\sigma\sigma'} N^{\sigma\sigma'}N^{\sigma'\sigma} = (N^\uparrow)^2 + (N^\downarrow)^2 \\), so sFLL becomes the usual
\\( \tfrac12 UN(N-1) - \tfrac12 J\sum\_\sigma N^\sigma(N^\sigma-1) \\).
In AMF, \\( v = v\_{\mathrm{HF}}[\tilde n] \\) is the exact derivative of \\( E\_{\mathrm{HF}}[\tilde n(n)] \\) because, for a rotationally invariant tensor, the component of
\\( \mathrm{Tr}\\, v\_{\mathrm{HF}}[\tilde n] \\) corresponding to the subtracted average vanishes (confirmed by check).

- (checked) sFLL with \\( J = 0 \\) agrees with the simplified form with \\( U = F^0 \\) in both energy and potential. sFLL and sAMF are invariant under global spin rotation.
- [Public] Manual 40.1.2: the leading `s` or `c` distinguishes whether the LDA/GGA part is treated by a spin-density functional (s) or a charge-density functional (c).
  [Convention] For `cFLL` and `cAMF`, from the **second SCF iteration onward**, exchange-correlation on the grid is evaluated with the density obtained by averaging the two spin densities (\\( n\_\uparrow = n\_\downarrow = n/2 \\))
  (the first iteration uses spin-polarized XC). The simplified form (Type 1) is implicitly equivalent to sFLL.
- [Public] The manual states that the enhancement of orbital polarization (Sec. 4) is not recommended in combination with `cFLL` or `cAMF`.

### 3.3 Coulomb tensor [Public] [Convention]

[Public] The matrix elements in complex spherical harmonics \\( Y\_{lm} \\) (Condon–Shortley) are, in terms of the Slater integrals \\( F^k \\) and Gaunt coefficients,

$$
\langle m\_1 m\_2 \vert V \vert m\_3 m\_4\rangle = \delta\_{m\_1+m\_2,\\, m\_3+m\_4} \sum\_{k=0,2,\dots,2l} c^k(m\_1, m\_3)\\, c^k(m\_4, m\_2)\\, F^k,
$$

$$
c^k(m, m') = \sqrt{\tfrac{4\pi}{2k+1}} \int d\Omega\\, Y^{\ast}\_{lm} Y\_{k,m-m'} Y\_{lm'}
{}= (-1)^{m} (2l+1) \begin{pmatrix} l & k & l\\ 0&0&0\end{pmatrix}\begin{pmatrix} l & k & l\\ -m & m-m' & m'\end{pmatrix}
$$

The transformation to real spherical harmonics \\( S\_a = \sum\_m T\_{am} Y\_{lm} \\) (\\( U^{(l)} \\) of the [PAO basis](pao-basis.md)) is
\\( U\_{abcd} = \sum \overline{T\_{a m\_1}}\\,\overline{T\_{b m\_2}}\\, T\_{c m\_3} T\_{d m\_4} \langle m\_1 m\_2 \vert V \vert m\_3 m\_4\rangle \\), and the result is real.
The averages are \\( U = (2l+1)^{-2}\sum\_{ab} U\_{abab} = F^0 \\) and \\( U - J = [2l(2l+1)]^{-1}\sum\_{ab}(U\_{abab} - U\_{abba}) \\).

[Convention] \\( F^k \\) are determined from the input \\( U, J \\) as follows (`scf.Yukawa off`, default).

| \\( l \\) | \\( F^0 \\) | \\( F^2 \\) | \\( F^4 \\) | \\( F^6 \\) |
|---|---|---|---|---|
| 0 | \\( U \\) | — | — | — |
| 1 | \\( U \\) | \\( 5J \\) | — | — |
| 2 | \\( U \\) | \\( 14J/(1+r) \\) | \\( rF^2 \\) | — |
| 3 | \\( U \\) | \\( 6435J/(286 + 195\cdot 0.668 + 250\cdot 0.494) \\) | \\( 0.668\\, F^2 \\) | \\( 0.494\\, F^2 \\) |

\\( r \\) is `scf.Slater.Ratio` (default 0.625). The f-orbital ratios 0.668 and 0.494 are fixed and cannot be changed from the input.
All satisfy \\( J = (F^2+F^4)/14 \\) (d), \\( F^2/5 \\) (p), \\( (286F^2 + 195F^4 + 250F^6)/6435 \\) (f) (checked).

[Convention] The existing implementation's transformation to real spherical harmonics agrees with the ordering and signs (all coefficients positive) of the [PAO basis](pao-basis.md).
Since \\( U\_{abcd} \\) depends on the signs of the individual real orbitals, a reimplementation must also use the same \\( T \\) as the PAO basis.

### 3.4 Estimation of J by Yukawa-type screening (`scf.Yukawa on`) [Public] [Convention]

[Public] Manual 40.2.2: only the input \\( U \\) is used; the screening length \\( \lambda \\) of the screened Coulomb potential \\( e^{-\lambda r}/r \\) is determined
so that \\( F^0 = U \\), and \\( F^2, F^4, F^6 \\) and \\( J \\) are built with the same \\( \lambda \\). `Hund.J.values` and `scf.Slater.Ratio` are ignored.

$$
F^k(\lambda) = (2k+1)\lambda \int\_0^\infty\\!\\!\int\_0^\infty dr\\, dr'\\, r^2 r'^2 R^2(r) R^2(r')\\, i\_k(\lambda r\_<)\\, \tilde k\_k(\lambda r\_>)
$$

\\( R \\) is the PAO radial function of the shell, \\( i\_k \\) is the modified spherical Bessel function (\\( i\_0(x) = \sinh x/x \\)), and \\( \tilde k\_k \\) is the function of the second kind satisfying the same recurrence as \\( \tilde k\_0(x) = e^{-x}/x \\)
(\\( e^{-\lambda R}/R = \lambda\sum\_k (2k+1) i\_k(\lambda r\_<)\tilde k\_k(\lambda r\_>) P\_k(\cos\gamma) \\); for \\( \lambda\to 0 \\) the usual \\( r\_<^k/r\_>^{k+1} \\). Checked).

[Convention] Numerical choices (those affecting the result):

- \\( \lambda \\) is found by bisection on the interval \\( [0.1, \lambda\_{\max}] \\) bohr⁻¹. \\( \lambda\_{\max} \\) is 10 if \\( U \ge 1.5 \\) Hartree, 20 if \\( U \ge 0.5 \\) Hartree, and 40 otherwise.
  It stops when \\( \lvert F^0 - U\rvert \le 0.0005 \\) eV; if not converged in 200 iterations, a warning is issued and the calculation proceeds.
- The radial integrals use Simpson's rule on \\( r \in [0, r\_c + 1] \\) bohr (\\( r\_c \\) is the PAO cutoff radius), with an outer step of 0.01 bohr and a step of 0.001 bohr for the inner cumulative integral.
- \\( J \\) is obtained from the relations in the table above (p: \\( F^2/5 \\), d: \\( (F^2+F^4)/14 \\), f: \\( (286F^2+195F^4+250F^6)/6435 \\)), and \\( \lambda \\), \\( U \\), \\( J \\), \\( F^k \\), \\( F^4/F^2 \\) are written to the output.

### 3.5 Applicability of the general form [Convention]

- The tensor is built only for shells with \\( U > 10^{-5} \\) Hartree or \\( J > 10^{-5} \\) Hartree; nothing is added for other shells.
- **The general form of the existing implementation handles only shells with \\( \zeta = 0, 1 \\) correctly.** If nonzero U or J is given for a shell with \\( \zeta \ge 2 \\), the computation of the shell's starting position is incorrect and the result is undefined.
- The general form in non-spin-polarized calculations is not defined in the existing implementation (the manual also requires spin polarization).
- The general form is effective only with `scf.Hubbard.U on`; the default of `scf.DFTU.Type` is 1.

## 4. Enhancement of orbital polarization [Public] [Convention]

[Public] Manual 40.1.3, TN6 Sec. 5: in the first few SCF iterations, the occupation matrix is replaced with an orbitally polarized form
so that the calculation more easily reaches an orbitally polarized solution. The input is the flag in the last column of `Atoms.SpeciesAndCoordinates` (`on`, `off`).

[Convention] Rules of the existing implementation:

- Conditions: `scf.Hubbard.U on`, SCF iteration number 1–8, the density matrix has not been read from a restart file (so it does not act from the second step of MD or geometry optimization onward),
  and the flag of the atom is not `off`.
- Target: among the shells of the atom with \\( l \ge 2 \\), **only those with \\( \zeta = 0 \\)**. The matrix is the corresponding part of the onsite matrix \\( \rho\_{II}(0) \\),
  regardless of the choice of `scf.Hubbard.Occupation`.
- Procedure: the matrix (for each case below) is diagonalized, and only if the sum of eigenvalues \\( N \\) exceeds 0.3, occupation 1 is given to the \\( \lfloor N\rfloor \\) eigenvectors with the largest eigenvalues,
  the fractional part \\( N - \lfloor N\rfloor \\) to the next one, and 0 to the rest; the matrix is rebuilt and replaces the corresponding block of the occupation matrix of that iteration (the one used for the potential and energy).
  The electron number is conserved.
  - Non-spin-polarized: the \\( (2l+1) \\)-dimensional one-spin matrix. With the flag `ex`, the occupations of the level given the fractional part and of the lowest of the levels with occupation 1 are exchanged (first excited configuration).
  - Collinear spin-polarized: the \\( 2(2l+1) \\)-dimensional matrix formed by placing the blocks of both spins side by side (spin off-diagonal 0) is diagonalized, and only the spin-diagonal blocks are put back.
    **The fractional part is given not to the "next level" but to the level one further down** (behavior of the existing implementation; one level with occupation 0 lies in between).
  - Noncollinear: \\( 2(2l+1) \\)-dimensional complex matrix. The fractional part goes to the next level. The rebuilt matrix \\( P \\) is stored transposed with respect to the orbital indices of each spin block
    (\\( n^{\sigma\sigma'}\_{ba} = P^{\sigma\sigma'}\_{ab} \\)).
  - The flag `ex` is not supported for collinear spin-polarized and noncollinear cases; all occupations of the corresponding block become 0.
- Before diagonalization, the existing implementation adds random perturbations of about \\( 10^{-13} \\) to \\( 10^{-12} \\) to the matrix elements. Their effect on the result is negligible, so they are unnecessary in a reimplementation [Design].

## 5. Constraint on the direction of noncollinear spins (`scf.Constraint.NC.Spin`) [Public] [Convention]

[Public] Manual Sec. 41: a penalty functional is added that drives the spin direction of each atom toward its initial direction (and, with `on2`, also the magnitude toward its initial value).
The strength is `scf.Constraint.NC.Spin.V` (eV, default 0), and atoms to be constrained have column 12 of the coordinate line set to 1.

[Convention] The 2×2 spin density matrix of atom \\( I \\) is built from the diagonal sum of the dual occupation matrices of all orbitals of that atom,

$$
T^{(I)}\_{\sigma\sigma'} = \sum\_{\alpha \in I} n^{\sigma'\sigma}\_{I,\alpha\alpha}
$$

(the charge-density convention, the same orientation as \\( n\_{\sigma\sigma'} \\) in TN2 Eq. 3, corresponding to the transpose of \\( n \\)). The total electron number and magnetization are

$$
N\_I = \mathrm{Tr}\\, T,\qquad m\_{I,c} = \mathrm{Tr}\\,(\sigma\_c T) \quad\Leftrightarrow\quad
\mathbf m = (2\\,\mathrm{Re}\\, T\_{\uparrow\downarrow},\ -2\\,\mathrm{Im}\\, T\_{\uparrow\downarrow},\ T\_{\uparrow\uparrow} - T\_{\downarrow\downarrow})
$$

which is the same convention as Sec. 2.4 of [Grid](grid-poisson-xc.md). The target matrix, using the initial direction \\( \hat{\mathbf e}\_0 = \hat{\mathbf e}(\theta\_0, \varphi\_0) \\) (Sec. 8 below), is

$$
T\_0 = \tfrac12\left( N\_I\\, \sigma\_0 + M\\, \hat{\mathbf e}\_0\cdot\boldsymbol\sigma \right),\qquad
M = \begin{cases} \lVert \mathbf m\_I \rVert & (\texttt{on}) \\ M\_0 = \lvert N^{\uparrow}\_{I,0} - N^{\downarrow}\_{I,0}\rvert & (\texttt{on2}) \end{cases}
$$

(\\( N^{\uparrow,\downarrow}\_{I,0} \\) are the input initial occupations). With the rotation matrix \\( U(\theta\_0,\varphi\_0) \\) of TN2 Eq. 16 it can also be written as
\\( T\_0 = U^\dagger \mathrm{diag}(\lambda\_+, \lambda\_-) U \\) (\\( \lambda\_\pm = (N\_I \pm M)/2 \\)) (checked). The penalty uses the Frobenius norm,

$$
E\_{\mathrm{cs}} = V \sum\_{I:\\,\text{constrained}} \lVert T^{(I)} - T\_0^{(I)} \rVert\_F^2
{}= \frac{V}{2}\sum\_I \lVert \mathbf m\_I - M\hat{\mathbf e}\_0 \rVert^2
$$

so that \\( E\_{\mathrm{cs}} = V\lVert\mathbf m\rVert^2 (1 - \cos\gamma) \\) for `on` (\\( \gamma \\) is the angle between \\( \mathbf m \\) and \\( \hat{\mathbf e}\_0 \\)), and \\( \tfrac V2 \lVert \mathbf m - M\_0\hat{\mathbf e}\_0\rVert^2 \\) for `on2`.
The output `Ucs` is this value. The potential, using the effective field \\( \mathbf h = \partial E\_{\mathrm{cs}}/\partial \mathbf m \\), is

$$
v^{\sigma\sigma'}\_{I,\alpha\beta} \mathrel{+}= \delta\_{\alpha\beta}\\, (\mathbf h\cdot\boldsymbol\sigma)\_{\sigma\sigma'}\quad(\alpha \in I),\qquad
\mathbf h = \begin{cases} V\left( 2\mathbf m - \hat{\mathbf m}\\,(\mathbf m\cdot\hat{\mathbf e}\_0) - \lVert\mathbf m\rVert\\,\hat{\mathbf e}\_0 \right) & (\texttt{on}) \\ V(\mathbf m - M\_0\hat{\mathbf e}\_0) & (\texttt{on2}) \end{cases}
$$

which is the exact derivative of \\( E\_{\mathrm{cs}} \\) with respect to \\( n \\) (checked). \\( v \\) is added as a dual Hamiltonian (Sec. 6).

- [Convention] The existing implementation uses a potential that includes, via the chain rule, the dependence of the target \\( T\_0 \\) on \\( T \\) (through the eigenvalues for `on`, through \\( N\_I \\) for `on2`).
  Whether that derivative agrees with the exact gradient above within rounding error has not been checked ([Open]).
- [Public] Manual Sec. 41: the constraint is a numerical constraint, not a physical external field, and can be combined with DFT+U. To avoid local minima, the recommended usage is to rotate from the easy axis in steps of about 30°, restarting from the previous result.

## 6. Zeeman terms (`scf.NC.Zeeman.Spin`, `scf.NC.Zeeman.Orbital`) [Public] [Convention]

[Public] TN2 Eqs. 139–150, manual Sec. 43: the interaction of a uniform magnetic field with the spin and orbital magnetic moments, \\( E\_z = \sum\_I \mathbf B^s\_I\cdot\mathbf s\_I + \tfrac12 \sum\_I \mathbf B^o\_I\cdot\mathbf l\_I \\).
The field magnitudes are `scf.NC.Mag.Field.Spin` and `scf.NC.Mag.Field.Orbital` (Tesla, default 0), and the directions are given by the angles in columns 8 and 9 (spin) and 10 and 11 (orbital) of the coordinate lines.

[Convention]

- Units: \\( B[\text{a.u.}] = B[\mathrm T] / 2.35051742\times 10^{5} \\) (relative difference \\( 6\times10^{-8} \\) from the CODATA 2018 value \\( 2.35051757\times10^5 \\) T).
- They act **only on atoms whose column 12 is 1** (both spin and orbital). The direction is the initial direction of that atom (after the correction in Sec. 8).
- Spin: with \\( \mathbf s\_I = \tfrac12 \mathbf m\_I \\) (from \\( T^{(I)} \\) of Sec. 5) and \\( \mathbf b\_I = -B\_s\\, \hat{\mathbf e}\_{s,I} \\),

$$
U\_{\mathrm{zs}} = \sum\_I \mathbf b\_I\cdot \mathbf s\_I = -\tfrac12 B\_s \sum\_I \hat{\mathbf e}\_{s,I}\cdot\mathbf m\_I,\qquad
v^{\sigma\sigma'}\_{I,\alpha\beta} \mathrel{+}= \delta\_{\alpha\beta}\\, \tfrac12 (\mathbf b\_I\cdot\boldsymbol\sigma)\_{\sigma\sigma'}
$$

  A positive \\( B\_s \\) acts to align the spin moment displayed by OpenMX (\\( \mathbf m \\)) with \\( \hat{\mathbf e}\_s \\).
  [Public] As stated in manual Sec. 37, the moment displayed by OpenMX follows the convention of opposite sign to the physical magnetic moment of the electron.
- Orbital: with \\( \mathbf a\_I = -\tfrac12 B\_o\\, \hat{\mathbf e}\_{o,I} \\) and \\( \mathbf L\_I \\) the orbital moment of Sec. 7,

$$
U\_{\mathrm{zo}} = \sum\_I \mathbf a\_I\cdot\mathbf L\_I,\qquad
v^{\sigma\sigma}\_{I,\alpha\beta} \mathrel{+}= \mathbf a\_I\cdot\langle \phi\_{I\alpha}\vert \hat{\mathbf L}\vert \phi\_{I\beta}\rangle \quad(\text{both spins, spin-diagonal only})
$$

  (\\( \langle\phi\_{I\alpha}\vert\hat L\_c\vert\phi\_{I\beta}\rangle \\) is purely imaginary; Sec. 7). This \\( v \\) is also added as a dual Hamiltonian, so
  its contribution to \\( \sum H\rho \\) agrees with the dual \\( \mathbf a\cdot\mathbf L \\). However, \\( \mathbf L\_I \\) in \\( U\_{\mathrm{zo}} \\) is evaluated from the output density matrix of the preceding diagonalization (before mixing).
- Input checks: all require noncollinear (`scf.SpinPolarization NC`) and `scf.Hubbard.Occupation dual` (the same applies to the constraint).
  The orbital Zeeman term requires spin-orbit coupling (`scf.SpinOrbit.Coupling on`). The spin Zeeman term and `scf.Constraint.NC.Spin` cannot be used simultaneously.
  If the conditions are not met, the program stops.
- `Uzs` and `Uzo` are included in `Utot` ([Total energy and forces](energy-force.md)). Subtract them when comparing with the energy without the external field.

## 7. Orbital magnetic moment [Public] [Convention]

[Public] TN2 Eqs. 151–162, manual Sec. 39: in noncollinear calculations, the expectation value of the orbital angular momentum of each atom is evaluated in the localized basis. Without spin-orbit coupling it is 0.

[Convention] The matrix elements of the orbital angular momentum on atom \\( I \\) are taken as

$$
\langle \phi\_{I\alpha} \vert \hat L\_c \vert \phi\_{I\beta}\rangle = i\\, \Lambda^{c}\_{\alpha\beta},\qquad
\Lambda^c\_{\alpha\beta} = \Big(\int r^2 R\_{\alpha} R\_{\beta}\\, dr\Big)\\, \mathrm{Im}\\,\langle S\_{l a}\vert \hat L\_c \vert S\_{l b}\rangle\\,\delta\_{l\_\alpha l\_\beta}
$$

(\\( \hbar = 1 \\); the real \\( \Lambda \\) is antisymmetric), and from the dual projection of the imaginary part of the spin-diagonal components of the density matrix,

$$
\tilde n^{\sigma}\_I = \tfrac12\left( \mathrm{Im}\\,\rho^{\sigma\sigma}\\, S + S\\, \mathrm{Im}\\,\rho^{\sigma\sigma}\right)\_{II}
$$

it is defined as

$$
L\_{I,c} = -\sum\_{\sigma}\sum\_{\alpha\beta\in I} \Lambda^c\_{\alpha\beta}\\, \tilde n^{\sigma}\_{I,\alpha\beta}
{}= \mathrm{Re}\sum\_\sigma\sum\_{\alpha\beta}\langle\phi\_{I\alpha}\vert\hat L\_c\vert\phi\_{I\beta}\rangle\\, n^{\sigma\sigma}\_{I,\alpha\beta}
$$

(\\( n \\) is dual. \\( \mathrm{Re}\\,n \\) is symmetric and cancels against the antisymmetric \\( \Lambda \\)).
For example, for p orbitals (in the order \\( x, y, z \\)), \\( \langle x\vert \hat L\_z\vert y\rangle = -i \\) and \\( \langle y\vert\hat L\_z\vert x\rangle = i \\) (the same as TN2 Eq. 162).

- \\( \mathbf L\_I \\) is displayed as is, in units of μB (the physical orbital magnetic moment is \\( -\mu\_B \mathbf L \\); the display follows the opposite-sign convention).
- It is computed after the diagonalization in each SCF iteration from the output (pre-mixing) density matrix. The total \\( \sum\_I \mathbf L\_I \\) and its angles are written to standard output, and
  the magnitude and angles per atom and the decomposition by orbital and multiplicity are written to the file `System.Name.OM` (also included in the "Orbital moments" section of `System.Name.out`).
- The magnitude and angles \\( (\theta, \varphi) \\) per atom are the polar coordinates of \\( \mathbf L\_I \\).

## 8. Input of the initial magnetic configuration [Public] [Convention]

[Public] Manual Sec. 37: the columns of each line of `Atoms.SpeciesAndCoordinates` are as follows.

| Column | Collinear (`on`, `off`) | Noncollinear (`NC`) |
|---|---|---|
| 1 | Serial number | Serial number |
| 2 | Species | Species |
| 3–5 | Coordinates | Coordinates |
| 6 | Initial up-spin occupation \\( N^\uparrow\_0 \\) | Same as left |
| 7 | Initial down-spin occupation \\( N^\downarrow\_0 \\) | Same as left |
| 8 | Orbital-polarization flag | Spin angle \\( \theta\_s \\) (degrees) |
| 9 | — | Spin angle \\( \varphi\_s \\) (degrees) |
| 10 | — | Orbital (moment and field) angle \\( \theta\_o \\) (degrees) |
| 11 | — | Orbital angle \\( \varphi\_o \\) (degrees) |
| 12 | — | Flag activating the constraint and Zeeman terms (integer, 1 enables) |
| 13 | — | Orbital-polarization flag |

[Convention]

- Initial density: the atomic density is scaled by \\( N^\uparrow\_0 + N^\downarrow\_0 \\) and distributed to up and down spins ([Density matrix](density-mixing.md) Sec. 6).
  In the noncollinear case, a 2×2 density with magnetization along \\( \hat{\mathbf e}(\theta\_s, \varphi\_s) \\) is built.
- **Reassignment to the majority-spin direction**: in the noncollinear case, for lines with \\( N^\uparrow\_0 < N^\downarrow\_0 \\), \\( N^\uparrow\_0 \\) and \\( N^\downarrow\_0 \\) are swapped,
  and the spin and orbital directions are each reversed (\\( \hat{\mathbf e} \to -\hat{\mathbf e} \\)). Hence the internal "initial direction" is always the direction of the majority-spin magnetization, and
  \\( \hat{\mathbf e}\_0 \\) of the constraint (Sec. 5) and the Zeeman terms (Sec. 6) is also this reassigned direction.
- Small angle fluctuations: when \\( \lvert\theta\rvert < 10^{-10} \\) degrees, a random number smaller than \\( 10^{-5} \\) degrees is added to \\( \theta\_s \\) and \\( \theta\_o \\).
  In addition, random numbers smaller than \\( 10^{-10} \\) rad are added to \\( \theta\_s, \varphi\_s \\) (to avoid the singularity of the azimuthal angle). The effect on the result is negligible.
- The orbital-polarization flag is `off`, `on`, or `ex` (case-insensitive), and is read only with `scf.Hubbard.U on`. Any other string stops the program.
  [Public] The manual states that the default when omitted is `off`, but the existing implementation does not initialize the value for lines lacking the column (the value from the previous line remains). It is safe to specify it explicitly for all atoms.
- For lines with missing columns, the existing implementation does not update the values of the missing columns (the initial values of the angle and flag arrays remain). [Design] A reimplementation uses angle 0, flag 0, and orbital polarization `off` as defaults.

## 9. DFT+U input [Public] [Convention]

| Keyword | Value | Default | Meaning |
|---|---|---|---|
| `scf.Hubbard.U` | on/off | off | Use DFT+U |
| `scf.DFTU.Type` | 1/2 | 1 | Simplified form / general form |
| `scf.Hubbard.Occupation` | onsite/full/dual | dual | Definition of the occupation matrix (also used for the constraint and Zeeman terms) |
| `<Hubbard.U.values` | Block | — | U per shell (eV) |
| `<Hund.J.values` | Block | — | J per shell (eV; Type 2 with `scf.Yukawa off`) |
| `scf.dc.Type` | sFLL/sAMF/cFLL/cAMF | sFLL | Double-counting correction (Type 2) |
| `scf.Slater.Ratio` | Real | 0.625 | \\( F^4/F^2 \\) for d (Type 2) |
| `scf.Yukawa` | on/off | off | Estimate J and \\( F^k \\) from the screening length (Type 2) |
| `scf.Constraint.NC.Spin` | off/on/on2 | off | Constraint on spin direction (and magnitude) |
| `scf.Constraint.NC.Spin.V` | Real (eV) | 0 | Strength of the constraint |
| `scf.NC.Zeeman.Spin`, `scf.NC.Zeeman.Orbital` | on/off | off | Zeeman terms |
| `scf.NC.Mag.Field.Spin`, `scf.NC.Mag.Field.Orbital` | Real (T) | 0 | Magnetic field magnitude |

[Convention]

- Option strings are case-insensitive.
- Each line of `<Hubbard.U.values` lists, after the species name, "label value" pairs for the \\( (l, \zeta) \\) of the basis of that species,
  in the order \\( l = 0, 1, \dots \\) and, within each \\( l \\), \\( \zeta = 0, 1, \dots \\), **as many as in the basis specification (e.g. `s2p2d2f1`)**
  (e.g. `Fe 1s 0.0 2s 0.0 1p 0.0 2p 0.0 1d 5.0 2d 0.0 1f 0.0`). The labels (`1d` etc.) are skipped and not checked.
  If the order and number do not match, values are assigned to other shells. As many lines as there are species are read, in any species order. `<Hund.J.values` has the same format.
- Conversion from eV to Hartree divides by 27.2113845.
- With `scf.DFTU.Type` 1, `scf.dc.Type` and J are not read.

## 10. Incorporation into the Hamiltonian [Public] [Convention]

[Public] TN6 Eqs. 21–28, TN2 Eqs. 122–124: the shell potentials \\( v\_s \\) are collected into a per-atom block-diagonal matrix \\( V \\) (the part of \\( V\_{II} \\) for shell \\( s \\) is \\( v\_s \\);
zero between shells and between atoms), and depending on the kind of projection

| Occupation matrix | Contribution to the Hamiltonian \\( H^U \\) | Range of atom pairs |
|---|---|---|
| onsite | \\( H^U\_{II}(0) = V\_{II} \\) | The atom itself only |
| full | \\( H^U = S V S \\) | Two-step reach of \\( S \\) (only components contained in the existing set of atom pairs) |
| dual | \\( H^U = \tfrac12 (V S + S V) \\), i.e. \\( H^U\_{I\alpha,J\beta}(\mathbf L) = \tfrac12[ (V\_{II} S\_{IJ}(\mathbf L))\_{\alpha\beta} + (S\_{IJ}(\mathbf L) V\_{JJ})\_{\alpha\beta} ] \\) | Same as \\( S \\) |

is added. In the noncollinear case the same expression is applied to each spin block \\( V^{\sigma\sigma'} \\) of \\( V \\), and the imaginary part is added to the imaginary part of the Hamiltonian (the same place as the spin-orbit coupling term).
The \\( v \\) of the constraint and Zeeman terms (Secs. 5 and 6) is added to the same \\( V \\) and enters the Hamiltonian through the projection of `scf.Hubbard.Occupation` (restricted to dual by the input check).

- By definition, \\( \sum\_{ab} H^U\_{ab}\rho\_{ab} = \sum\_s\sum\_{ab} v\_{s,ab}\\, n\_{s,ab} \\) holds (checked; for collinear non-polarized, both sides are for one spin).
- [Convention] For full, \\( S V S \\) is built only for components whose \\( (I, J, \mathbf L) \\) is contained in the set of atom pairs of the sparse block format (out-of-range components are discarded).

## 11. Forces and stress [Public] [Convention]

[Public] TN6 Eqs. 29–30: \\( \partial E\_U/\partial\tau = \sum v\\, \partial n/\partial \tau \\). The change of \\( n \\) through \\( \rho \\) is already included
through \\( H^U \\) in the Hamiltonian and the \\( \partial S \\) term of the energy-weighted density matrix (the Pulay identity in [Total energy and forces](energy-force.md)),
so what remains is the change through \\( S \\) with \\( \rho \\) fixed. With \\( V \\) fixed,

$$
\frac{\partial E\_U}{\partial \mathbf R\_K}\bigg|\_{\rho} = g\_s \sum\_{\sigma\sigma'} \sum\_{\mathrm{pairs}} \rho^{\sigma\sigma'}\_{ab}\\, \frac{\partial H^{U,\sigma\sigma'}\_{ab}}{\partial \mathbf R\_K}\bigg|\_{V},\qquad
\frac{\partial H^U}{\partial \mathbf R\_K} = \begin{cases} 0 & \text{onsite}\\ \partial S\\, V S + S V \partial S & \text{full}\\ \tfrac12(V\\,\partial S + \partial S\\, V) & \text{dual}\end{cases}
$$

(\\( g\_s \\) is 2 for non-polarized and 1 otherwise; the real part is taken). Only atom-pair blocks containing atom \\( K \\) have nonzero \\( \partial S/\partial\mathbf R\_K \\),
and the derivative of the overlap matrix ([Two-center integrals](two-center-overlap-kinetic.md)) is used directly. The force is \\( \mathbf F\_K = -\partial E/\partial\mathbf R\_K \\).
This is confirmed by finite differences (`tools/check_dftu_magnetism.py`, dual and full, collinear and noncollinear).

[Convention]

- The constraint and Zeeman terms also enter the same \\( V \\), so they are included in the force by the same formula. The \\( \Lambda \\) of the orbital Zeeman term is only in the atom's own block and does not depend on positions.
- The existing implementation adds the above terms for both full and dual, for both collinear and noncollinear (onsite gives 0).
- Stress: full cannot be combined with stress (stops at input). Stress calculation itself is not supported for noncollinear. The contribution to the stress for collinear dual is [Open].

## 12. Magnetic anisotropy energy (second-variation method) [Public] [Convention]

[Public] Manual Sec. 42: for a ferromagnet, (1) converge a collinear calculation without spin-orbit coupling, (2) convert its restart file
to a noncollinear calculation with `scf.restart c2n` and diagonalize once including spin-orbit coupling with the spin direction `scf.Restart.Spin.Angle.Theta`, `scf.Restart.Spin.Angle.Phi` (degrees), and
(3) take the difference of `Uele` (band energy) between directions as the MAE (force theorem; the double-counting terms depending only on the input density do not depend on the direction, so `Uele` rather than `Utot` is compared).

[Convention]

- The defaults of `scf.Restart.Spin.Angle.Theta` and `...Phi` are 0.
- In the conversion, the collinear grid densities \\( n\_\uparrow, n\_\downarrow \\) are rotated to the same direction \\( \hat{\mathbf e}(\theta,\varphi) \\) at all grid points:
  \\( \underline n = U^\dagger \mathrm{diag}(n\_\uparrow, n\_\downarrow) U = \tfrac12[(n\_\uparrow+n\_\downarrow)\sigma\_0 + (n\_\uparrow - n\_\downarrow)\hat{\mathbf e}\cdot\boldsymbol\sigma] \\)
  (\\( U \\) is the matrix of TN2 Eq. 16, charge-density convention).
- Details of the conversion and restart (treatment of the density matrix and Hamiltonian, relation to `scf.maxIter`) are [Open].

## Edge cases [Convention] [Design]

- A shell with all U equal to 0 adds nothing in the simplified form (\\( v = 0 \\)) and is not handled in the general form if below the threshold \\( 10^{-5} \\) Hartree.
- Shells with \\( \zeta \ge 2 \\) in the general form, and the non-polarized general form: undefined in the existing implementation (Sec. 3.5). [Design] A reimplementation raises an error at the input stage.
- The assignment of the fractional part in the collinear spin-polarized enhancement of orbital polarization, and the lack of support for `ex` (Sec. 4), are behaviors of the existing implementation. [Design] They are reproducible only for comparison; the default is "fractional part to the next level; `ex` is an error except for non-polarized".
- For the constraint `on`, \\( \hat{\mathbf m} \\) is undefined for atoms with \\( \lVert\mathbf m\rVert \to 0 \\). Since \\( E\_{\mathrm{cs}} \\) and \\( \mathbf h \\) go continuously to 0 in \\( \lVert\mathbf m\rVert \\), set \\( \mathbf h = 0 \\) for \\( \lVert\mathbf m\rVert < 10^{-12} \\) [Design].
- Since the occupation matrix is built from the density matrix after mixing, its eigenvalues may fall outside \\( [0,1] \\) during the SCF. The energy and potential expressions are used as is (no clipping).
- For onsite and full, \\( \sum \mathrm{Tr}\\,n \neq N\_{\mathrm e} \\), and shell occupations may exceed 1.

## Requirements for implementations [Design]

- **Occupation-matrix projection**: make `Projector` (onsite, full, dual) an enumeration with three functions: (a) building \\( n \\) of each shell (or each atom) from \\( \rho \\) and \\( S \\),
  (b) building \\( H^U \\) from \\( V \\) and \\( S \\), and (c) building the force contribution from \\( V \\), \\( \rho \\), and \\( \partial S \\).
  Tests guarantee that the three satisfy the same convention (\\( \sum H^U\rho = \sum v n \\)).
- **Interaction functional**: an interface that takes the shell \\( n\_s \\) and returns \\( (E\_s, v\_s) \\). Implementations are `Dudarev { u }` and
  `General { tensor, u, j, dc: {SFLL, SAMF, CFLL, CAMF} }`. Collinear is handled by the same code as a special case of noncollinear (spin off-diagonal 0), and for non-polarized the caller multiplies by \\( g\_s = 2 \\).
- **Coulomb tensor**: separate a pure function building \\( U\_{abcd} \\) from \\( (l, F^k) \\), a function building \\( F^k \\) from \\( (U, J, r) \\), and the Yukawa estimation. The transformation matrix to real spherical harmonics is shared with the PAO basis.
- **Per-atom spin terms**: implement the constraint (`on`, `on2`) and the Zeeman terms as terms that take the atom's \\( T^{(I)} \\) (and \\( \mathbf L\_I \\)) and return \\( (E, \mathbf h) \\);
  the conversion to \\( v \\) (placing \\( \mathbf h\cdot\boldsymbol\sigma \\) on the diagonal of all orbitals of the atom) is shared. The constraint uses the exact gradient.
- **Connection to the SCF**: make the order explicit: \\( \rho \\) after mixing → \\( n \\) → \\( (E, V) \\) → \\( H \\) of the next iteration, and guarantee by types that the energy and Hamiltonian are built from the same \\( n \\).
  Return `Uhub`, `Ucs`, `Uzs`, `Uzo` separately.
- **Input checks**: check all conditions of Secs. 8 and 9 (whether noncollinear, dual, and spin-orbit coupling are required; exclusivity of the constraint and spin Zeeman; full and stress; number and order of `Hubbard.U.values`) at read time,
  and also check the labels (`1d` etc.) against the basis (stricter than the existing implementation).
- **Output**: allow output equivalent to `System.Name.DM_onsite` (occupation matrices per atom and spin with their eigenvalues and eigenvectors) and the decomposition of orbital moments.
- **Switches for comparison with the existing implementation**: make the position of the fractional part in the enhancement of orbital polarization, the transposition in the noncollinear case, and the handling of `ex` reproducible only in comparison mode.

## Tests [Design]

### Analytic checks and invariants (confirmed with `tools/check_dftu_magnetism.py`)

| Content | Tolerance |
|---|---|
| Real spherical harmonics (p, d) agree with the PAO basis table | \\( 10^{-14} \\) |
| \\( U\_{abcd} \\) is real; symmetries \\( U\_{abcd} = U\_{badc} = U\_{cdab} = U\_{cbad} \\) | \\( 10^{-14} \\) |
| \\( (2l+1)^{-2}\sum U\_{abab} = F^0 \\); \\( J \\) obtained from the averages equals the input \\( J \\) (l = 1, 2, 3) | \\( 10^{-14} \\) |
| \\( U\_{abcd} \\) is invariant under an arbitrary rotation \\( D \\) | \\( 10^{-11} \\) |
| \\( v = \partial E/\partial n \\) for the simplified and general forms (four double countings), collinear and noncollinear (central difference, relative) | \\( 10^{-8} \\) |
| sFLL with \\( J = 0 \\) = simplified form with \\( U = F^0 \\) (energy and potential) | \\( 10^{-12} \\) |
| \\( E\_U = 0 \\) of the simplified form for idempotent \\( n \\); simplified form, sFLL, sAMF invariant under global spin rotation | \\( 10^{-12} \\) |
| Dual sum rule \\( \sum\mathrm{Tr}\\,n = N\_{\mathrm e} \\), \\( \sum H^U\rho = \sum v n \\) (dual and full, collinear and noncollinear) | \\( 10^{-12} \\) |
| Formula for \\( \partial E\_U/\partial\mathbf R \\) at fixed \\( \rho \\) vs. central difference (dual and full, collinear and noncollinear) | \\( 10^{-8} \\) |
| Constraint: closed form of \\( V\lVert T-T\_0\rVert\_F^2 \\), \\( v = \mathbf h\cdot\boldsymbol\sigma \\) vs. central difference (on, on2) | \\( 10^{-8} \\) |
| Spin Zeeman: \\( v = \tfrac12\mathbf b\cdot\boldsymbol\sigma \\) vs. central difference | \\( 10^{-8} \\) |
| Occupation rule for orbital polarization conserves the electron number | \\( 10^{-12} \\) |
| Yukawa expansion and the \\( \lambda\to 0 \\) limit | \\( 10^{-12} \\), \\( 10^{-3} \\) |

### Additional unit tests [Design]

- For a small two-atom system with given \\( \rho \\) and \\( S \\), compare onsite, full, and dual \\( n \\) with hand calculations (\\( 10^{-14} \\)).
- Orbital Zeeman term: for an atom with only a p shell, \\( \Lambda^z \\) reproduces \\( \langle x\vert L\_z\vert y\rangle = -i \\), and the occupation of the \\( L\_z \\) eigenstate \\( (x + iy)/\sqrt2 \\) determines the sign of \\( L\_z = \pm 1 \\).
- Constraint `on`: \\( E\_{\mathrm{cs}} = 0 \\) and \\( \mathbf h = 0 \\) when \\( \mathbf m \parallel \hat{\mathbf e}\_0 \\); \\( E = 2V\lVert\mathbf m\rVert^2 \\) when antiparallel.
- Yukawa estimation: for a known radial function (e.g. hydrogen-like 3d), the bisection for \\( \lambda \\) satisfies \\( F^0 = U \\) within 0.0005 eV.

### Comparison with the existing implementation (fixtures)

| System | Settings | Quantities compared | Tolerance |
|---|---|---|---|
| NiO (AFM, collinear, example of manual Sec. 40) | Type 1, dual, U on Ni 1d | `Uhub`, `Utot`, Ni spin moment, eigenvalues of `DM_onsite` | Total energy \\( 10^{-6} \\) Hartree, moment \\( 10^{-4} \\) μB, occupation \\( 10^{-5} \\) |
| Same | onsite, full | `Uhub`, forces (also compared with finite differences) | Same, forces \\( 10^{-5} \\) Hartree/bohr |
| Same | Type 2, sFLL, sAMF, cFLL, cAMF, `scf.Yukawa` on/off | `Uhub`, output \\( \lambda \\), \\( J \\), \\( F^k \\) | \\( 10^{-6} \\) Hartree, \\( \lambda \\) to \\( 10^{-4} \\) |
| FeO (noncollinear + spin-orbit, equivalent to `work/FeO_NC.dat`) | Type 1, dual | `Uhub`, spin and orbital moments and their angles | \\( 10^{-6} \\) Hartree, \\( 10^{-4} \\) μB, angles \\( 10^{-3} \\) degrees |
| Same | Constraint on, on2 (V = 0.5 eV), Zeeman spin and orbital (100 T) | `Ucs`, `Uzs`, `Uzo`, angles | \\( 10^{-6} \\) Hartree |
| FePt (manual Sec. 42) | c2n, θ = 0° and 90° | Difference of `Uele` (MAE) | \\( 10^{-6} \\) Hartree |

- Comparison is done in two stages: a one-shot comparison that builds the Hamiltonian once from the same density matrix (eliminating differences in the SCF path), and a comparison after SCF convergence.
- Fixtures using the enhancement of orbital polarization are compared in comparison mode (reproducing the behavior of Sec. 4).

## [Open]

- Whether the derivative of the target \\( T\_0 \\) used by the existing implementation in the constraint potential agrees with the exact gradient (if not, `Ucs` and the angles after SCF convergence will differ).
- How the XC averaging of `cFLL` and `cAMF` works with GGA (confirmed for LDA), and whether the averaged density is also used for the energy evaluation of `Uxc`.
- The contribution to the stress for collinear dual and onsite.
- The value of U when the `<Hubbard.U.values` block is absent (should be 0 but unconfirmed).
- How quantities other than the grid density (density matrix, Hamiltonian, occupation matrix) are handled on restart for the second-variation method (`c2n`).
- Whether the transposition of the orbital indices in the noncollinear enhancement of orbital polarization is intended (reproduced in comparison mode, but whether to make it the default behavior requires discussion).
- Whether the rules written here (especially Secs. 3.5 and 4) are the same in upstream OpenMX 3.9/4.0 and in the GPU version consulted.
