# Geometry optimization, molecular dynamics, and NEB

> Status: first draft (awaiting review). For every option of `MD.Type` (geometry optimization, lattice optimization, molecular dynamics, NEB),
> this chapter rewrites the equations of the public documents (manual sections 17–20 and 51, technical note 7 "Geometry Optimization", the NEB technical note, standard papers)
> in the notation of this note, and adds conventions and defaults absent from the public documents, as confirmed from the existing implementation.
> The parts checked numerically are in `tools/check_geometry_md.py`.

## Scope of this chapter

- This chapter covers the part that, after the SCF has converged and the forces (and stress) have been obtained, advances the atomic coordinates, lattice vectors, and velocities to the next step.
- The definitions and signs of forces and stress are in [Total energy and forces](energy-force.md), the initial density (extrapolation) for the SCF of the next step is in
  sections 6 and 7 of [Density matrix and charge mixing](density-mixing.md), and the grid origin is in section 1.4 of [Real-space grid](grid-poisson-xc.md).
  This chapter only gives entry points to them.
- For each value of `MD.Type`, it specifies what the variables are, by which equations they are updated, and what is output.
  The goal of a reimplementation is "to converge to the same structure for the same input" and "to generate trajectories of the same statistical ensemble";
  step-by-step agreement with the existing implementation is required only where explicitly stated (compatibility mode for comparison).

## Notation [Design]

In addition to the symbols of [Notation](../../../conventions/notation.md), the following are used. Internal units are atomic units (bohr, Hartree, electron mass \\( m\_e \\),
time unit \\( \hbar/E\_h \\)).

| Symbol | Meaning |
|---|---|
| \\( N \\) | Number of atoms (in the unit cell) |
| \\( \mathbf R\_I \\), \\( \mathbf x \in \mathbb R^{3N} \\) | Position of atom \\( I \\), and the vector concatenating the coordinates of all atoms |
| \\( U \\) | Total energy `Utot` ([Total energy and forces](energy-force.md)) |
| \\( \mathbf g = \partial U/\partial\mathbf x \\), \\( \mathbf F = -\mathbf g \\) | Gradient and force (Hartree/bohr) |
| \\( m\_I \\) | Mass of atom \\( I \\) (internally in units of \\( m\_e \\)) |
| \\( \mathbf v\_I \\), \\( K = \sum\_I \tfrac12 m\_I \lvert\mathbf v\_I\rvert^2 \\) | Velocity and nuclear kinetic energy `Ukc` |
| \\( \Delta t \\), \\( n \\) | Time step and step number (1-based; cumulative across restarts) |
| \\( A \\) | 3×3 matrix whose **rows** are the lattice vectors (\\( A\_{ij} \\) is the \\( j \\) component of \\( \mathbf a\_i \\)) |
| \\( \Omega = \det A \\) | Unit-cell volume |
| \\( \boldsymbol\sigma\_i \\), \\( S \\) | \\( \boldsymbol\sigma\_1 = \mathbf a\_2\times\mathbf a\_3 \\) etc. (cyclic). The matrix having these as rows is \\( S = \Omega A^{-\mathsf T} \\) |
| \\( \mathbf s\_I \\) | Fractional coordinates \\( \mathbf R\_I = \mathbf O + A^{\mathsf T}\mathbf s\_I \\) (\\( \mathbf O \\) is the grid origin) |
| \\( \Sigma\_{\gamma\eta} \\) | Stress of [Total energy and forces](energy-force.md) (\\( \Sigma/\Omega \\) is the output stress tensor; diagonal positive under tension) |
| \\( k\_B \\), \\( T \\), \\( T\_0 \\) | Boltzmann constant, kinetic temperature, target temperature |
| \\( \zeta \\), \\( \eta \\), \\( Q \\) | Nosé–Hoover friction coefficient, its time integral (thermostat coordinate), thermostat mass |
| \\( W \\), \\( P \\) | Barostat mass and external pressure |

## 1. Overall flow [Public] [Convention]

[Public] Manual sections 18 and 20: geometry optimization and MD select `MD.Type` in the same input file and, at each step, repeat
SCF → forces → coordinate update.

[Convention] The processing of one step is in the following order.

1. Converge the SCF at the current coordinates (and lattice), and obtain `Utot`, forces, and, if needed, the stress.
2. Append one block with the current coordinates and forces (and the velocities held at that point) to the `.md` file (section 10).
3. In molecular dynamics, apply the correction that makes the sum of forces zero (below).
4. Update coordinates, lattice, and velocities according to `MD.Type`, and append one line to the `.ene` file.
5. Only in molecular dynamics (`NVE`, `NVT_*`, `NPT_*`), move atoms that have left the unit cell back inside (section 7). Not done in geometry optimization.
6. Write the restart input file `System.Name.dat#` from the current state (section 10).

- The step number is the cumulative \\( n = n\_{\mathrm{local}} + \\) `MD.Current.Iter` (`MD.Current.Iter` defaults to 0 and is written by the restart file).
  The iteration ends when the convergence criterion is met or \\( n \\) reaches `MD.maxIter` (default 1).
  `MD.maxIter` is the upper bound on the **cumulative** number of steps including steps before the restart.
- At the step with \\( n = \\) `MD.maxIter`, the coordinates are not updated even if not converged (the run ends with the structure and forces of the last SCF still corresponding).
- In molecular dynamics (`NVE`, `NVT_VS`, `NVT_VS2`, `NVT_VS4`, `NVT_NH`, `NVT_Langevin`), before the update,
  the **simple average** of the gradient (not mass-weighted) is subtracted per Cartesian component, over the unfixed components:
  \\( g\_{I\alpha} \leftarrow g\_{I\alpha} - \frac{1}{N\_\alpha}\sum\_{J\\,\text{free}} g\_{J\alpha} \\) (\\( N\_\alpha \\) is the number of atoms whose \\( \alpha \\) component is free).
  NPT and geometry optimization do not apply this correction (same content as "Signs and units" in [Total energy and forces](energy-force.md)).
- From the second step on, the SCF restarts from the density of the previous step and extrapolates the density (sections 6 and 7 of [Density matrix and charge mixing](density-mixing.md)).
  The grid origin \\( \mathbf O \\), once determined at the first step, is not changed ([Real-space grid](grid-poisson-xc.md) section 1.4).

## 2. Units and constants [Convention]

| Quantity | Input unit | Conversion to internal units |
|---|---|---|
| Time step `MD.TimeStep` | fs (default 0.5; negative is an input error) | \\( 1\ \mathrm{fs} = 41.3411 \\) |
| Mass (4th column of `Definition.of.Atomic.Species`) | unified atomic mass unit | \\( 1\ \mathrm u = 1.660538921/0.000910938291\ m\_e \approx 1822.888\ m\_e \\) |
| Velocity (`MD.Init.Velocity`, `.md` output) | m/s | \\( 1\ \mathrm{m/s} = 0.4571028\times 10^{-6} \\) |
| Pressure (`MD.applied.pressure`, `MD.CellPressureControl`) | GPa | \\( 1\ \mathrm{GPa} = 0.3398827\times 10^{-4}\ \mathrm{Hartree/bohr^3} \\) |
| Temperature | K | \\( k\_B = 8.617251324\times10^{-5}\ \mathrm{eV/K} \\), \\( 1\ \mathrm{Hartree} = 27.2113845\ \mathrm{eV} \\) |
| Length | Å or bohr | \\( 1\ \mathrm{bohr} = 0.529177249\ \text{Å} \\) |

- The GPa conversion factor differs from the CODATA value (\\( 3.39893\times10^{-5} \\)) by a relative \\( 3\times10^{-5} \\). For comparison, use the factor of the existing implementation.
- The default mass is the atomic mass of the most abundant isotope of each element (e.g. H 1.00782503207, C 12.00, O 15.99491461956, Si 27.9769265325, Fe 55.9349375).
  [Public] Manual section 20.9: it can be overridden by the 4th column of `Definition.of.Atomic.Species`.
- [Convention] There are two definitions of temperature, used in different places.
  - \\( T\_A = K/(\tfrac32 k\_B N) \\): divides by the total number of atoms \\( N \\). Used for the `.ene` temperature, thermostat control, and NPT.
  - \\( T\_B = 2K/(k\_B N\_{\mathrm{free}}) \\): \\( N\_{\mathrm{free}} \\) is the number of unfixed components.
    Used for the temperature in the `.md` header line of `NVE`, `NVT_VS`, `NVT_NH`, `NVT_VS2`, `NVT_Langevin`.
  - In both cases \\( K \\) is summed only over unfixed components (all components for the `NVT_VS` family; section 6.4).

## 3. List of `MD.Type` [Public] [Convention]

[Public] Manual sections 18–20 and 51. The default is `NOMD` [Convention].

| Value | Kind | Variables | Method | Section |
|---|---|---|---|---|
| `NOMD` | single SCF | — | Coordinates are not moved. Writes one line to `.ene` and coordinates and forces to `.out`, then ends | — |
| `Opt` | geometry optimization | atomic coordinates | steepest descent (variable prefactor) | 4.2 |
| `DIIS` | ditto | ditto | GDIIS + identity Hessian | 4.4 |
| `BFGS` | ditto | ditto | GDIIS + BFGS inverse Hessian | 4.5 |
| `RF` | ditto | ditto | GDIIS + BFGS Hessian + rational function (RF) step | 4.6 |
| `EF` | ditto | ditto | GDIIS + BFGS Hessian + eigenvalue lower bound | 4.7 |
| `OptC1`–`OptC4` | lattice optimization | lattice (fractional coordinates fixed) | steepest descent | 5.2 |
| `OptC5`–`OptC7` | ditto | lattice and atomic coordinates | steepest descent | 5.2 |
| `RFC5`–`RFC7` | ditto | lattice and fractional coordinates | GDIIS + BFGS + RF | 5.3 |
| `EvsLC` | energy vs. lattice constant | lattice | scale the lattice by a fixed increment | 5.5 |
| `DF` | for the delta factor | volume | compute at 7 volumes | 5.5 |
| `NVE` | molecular dynamics | coordinates, velocities | velocity Verlet | 6.2 |
| `NVT_VS`, `NVT_VS2`, `NVT_VS4` | ditto | ditto | velocity scaling (global, per atom, per atom group) | 6.4 |
| `NVT_NH` | ditto | ditto | Nosé–Hoover | 6.5 |
| `NVT_Langevin` | ditto | ditto | Langevin thermostat | 6.6 |
| `NPT_VS_PR`, `NPT_VS_WV`, `NPT_NH_PR`, `NPT_NH_WV` | ditto | coordinates, velocities, lattice | velocity scaling / Nosé–Hoover with Parrinello–Rahman / Wentzcovitch | 6.7 |
| `NEB` | reaction path | coordinates of all images | nudged elastic band + DIIS/BFGS | 8 |

[Convention] Other handling:

- `RFC1`–`RFC4` are accepted, but the program prints "not supported" and exits.
- `Constraint_DIIS` and `Opt_LBFGS` are accepted but perform no update at all (the SCF is repeated on the same structure up to `MD.maxIter`).
  A reimplementation may treat them as input errors [Design].
- The program exits if `MD.Opt.DIIS.History` exceeds 19.

## 4. Geometry optimization (atomic coordinates)

### 4.1 Convergence criterion and maximum displacement [Public] [Convention]

[Public] Manual section 18.2: `MD.Opt.criterion` (default 0.0003 Hartree/bohr) is compared with the maximum force acting on the atoms.

[Convention] How the "maximum force" is taken differs between methods. Fixed components (section 4.10) are excluded.

| Method | Quantity compared |
|---|---|
| `Opt` (steepest descent) | maximum of per-atom vector magnitudes \\( \max\_I \lVert\mathbf g\_I\rVert \\) |
| `DIIS`, `BFGS`, `RF`, `EF` | maximum of per-component absolute values \\( \max\_{I\alpha}\lvert g\_{I\alpha}\rvert \\) |

- The test is made with the gradient after the SCF of that step; if satisfied, the run ends without updating the coordinates at that step.
- Maximum displacement per step: in the quasi-Newton family, if the per-component maximum of the update \\( \max\_{I\alpha}\lvert\Delta x\_{I\alpha}\rvert \\)
  exceeds \\( d\_{\max} = 0.5 \\) bohr, the whole update vector is scaled by \\( d\_{\max}/\max\lvert\Delta x\rvert \\) (direction preserved).
  In steepest descent, it is scaled in the same way if \\( s\max\_I\lVert\mathbf g\_I\rVert \\) (below) exceeds 0.5 bohr.
- [Convention] If the GDIIS linear equation cannot be solved (singular coefficient matrix), the existing implementation exits **as if converged**.
  A reimplementation discards the history and falls back to steepest descent, or reports an error [Design].

### 4.2 Steepest descent (`Opt`) [Public] [Convention]

[Public] Manual section 18.1: steepest descent whose prefactor changes according to the convergence behavior.

[Convention] The update is \\( \mathbf x\_{n+1} = \mathbf x\_n - c\\, s\_n\\, \mathbf g\_n \\), where \\( c \\) is the reduction factor for the maximum displacement of section 4.1. The prefactor \\( s \\) (bohr²/Hartree) is as follows.

- Reference value \\( s\_0 = (2\ \mathrm{fs})^2/(1\ \mathrm u) = (82.6822)^2/1822.888 \approx 3.750 \\). The range is restricted to \\( [0.04 s\_0,\ 15 s\_0] \\).
- At the first step (and at the first step after switching to SD, section 4.8): \\( s = 1/(1.5\times 0.529177249) \approx 1.260 \\)
  (the value at which the atom with the largest force moves by about 1/1.5 of \\( \lVert\mathbf g\rVert \\) Å).
- Thereafter, with \\( U\_n \\) the `Utot` of each step,
  - if \\( U\_n > U\_{n-1} \\), then \\( s \leftarrow s/2 \\);
  - otherwise, if \\( U\_{n-1} < U\_{n-2} \\) and \\( U\_n < U\_{n-1} \\), and \\( n \bmod 4 = 1 \\), then \\( s \leftarrow 2.5 s \\).
  - The energy history used for the comparison is recorded from the second step (the value of the first step does not enter the comparison).

### 4.3 GDIIS (common part) [Public]

[Public] Technical note 7, Eqs. 6–10 (the GDIIS of Csaszar–Pulay, called RMM-DIIS in the technical note).
From the coordinates \\( \mathbf x^{(m)} \\) and gradients \\( \mathbf g^{(m)} \\) of the last \\( p \\) steps, find the coefficients
that minimize \\( \lVert\sum\_m c\_m\mathbf g^{(m)}\rVert^2 \\) subject to \\( \sum\_m c\_m = 1 \\):

$$
\begin{pmatrix} B & \mathbf 1 \\\\ \mathbf 1^{\mathsf T} & 0 \end{pmatrix}
\begin{pmatrix} \mathbf c \\\\ \lambda/2 \end{pmatrix}
{}= \begin{pmatrix} \mathbf 0 \\\\ 1 \end{pmatrix},\qquad
B\_{mm'} = \mathbf g^{(m)}\cdot\mathbf g^{(m')}
$$

$$
\tilde{\mathbf x} = \sum\_m c\_m \mathbf x^{(m)},\qquad
\tilde{\mathbf g} = \sum\_m c\_m \mathbf g^{(m)},\qquad
\mathbf x\_{n+1} = \tilde{\mathbf x} + \Delta\mathbf x(\tilde{\mathbf g})
$$

[Convention]

- \\( p = \min(k, \\) `MD.Opt.DIIS.History`\\( ) \\), where \\( k \\) is the number of steps since entering the quasi-Newton stage (section 4.8). The default of `MD.Opt.DIIS.History` is 2.
  However, the coordinates and gradients of the steepest-descent stage are also in the history.
- \\( B \\) is divided by the maximum absolute value of its elements before solving (\\( \mathbf c \\) is unchanged). It is solved with a direct method for symmetric indefinite matrices.
- Gradients of fixed components are stored in the history as 0.
- The closed form is \\( \mathbf c = B^{-1}\mathbf 1/(\mathbf 1^{\mathsf T}B^{-1}\mathbf 1) \\) (checked in `tools/check_geometry_md.py`).

### 4.4 `DIIS` [Public] [Convention]

[Public] Technical note 7, Eq. 11: the Hessian is approximated by the identity, \\( \mathbf x\_{n+1} = \tilde{\mathbf x} - \alpha\tilde{\mathbf g} \\).

[Convention] \\( \alpha \\) (bohr²/Hartree) is determined by \\( r = \lVert\tilde{\mathbf g}\rVert/N \\) (the 3N-dimensional Euclidean norm divided by the number of atoms).

| \\( r \\) | \\( > 10^{-2} \\) | \\( (10^{-3}, 10^{-2}] \\) | \\( (10^{-4}, 10^{-3}] \\) | \\( (10^{-5}, 10^{-4}] \\) | \\( (10^{-6}, 10^{-5}] \\) | \\( \le 10^{-6} \\) |
|---|---|---|---|---|---|---|
| \\( \alpha \\) | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 |

### 4.5 `BFGS` [Public] [Convention]

[Public] Technical note 7, Eqs. 12–14: the inverse Hessian \\( H^{-1} \\) is updated by BFGS, and \\( \mathbf x\_{n+1} = \tilde{\mathbf x} - H^{-1}\tilde{\mathbf g} \\).
With \\( \Delta\mathbf x = \mathbf x^{(n)} - \mathbf x^{(n-1)} \\), \\( \Delta\mathbf g = \mathbf g^{(n)} - \mathbf g^{(n-1)} \\) (differences of the raw history, not of the GDIIS averages),

$$
H^{-1}\_{\mathrm{new}} = H^{-1} + \frac{\Delta\mathbf x\cdot\Delta\mathbf g + \Delta\mathbf g\cdot H^{-1}\Delta\mathbf g}{(\Delta\mathbf x\cdot\Delta\mathbf g)^2}\\, \Delta\mathbf x\\,\Delta\mathbf x^{\mathsf T}
{}- \frac{H^{-1}\Delta\mathbf g\\,\Delta\mathbf x^{\mathsf T} + \Delta\mathbf x\\,\Delta\mathbf g^{\mathsf T}H^{-1}}{\Delta\mathbf x\cdot\Delta\mathbf g}
$$

[Convention] At the first step of the quasi-Newton stage, \\( H^{-1} = I \\) (in units of bohr²/Hartree), and it is updated at every subsequent step.
The update is not stopped even when \\( \Delta\mathbf x\cdot\Delta\mathbf g \le 0 \\). `MD.Opt.Init.Hessian` is not used (manual section 18.3 also states that the model Hessian cannot be used with BFGS or DIIS).

### 4.6 `RF` [Public] [Convention]

[Public] Technical note 7, Eqs. 15–19 (the rational function method of Banerjee et al. and Baker): the Hessian \\( H \\) is updated directly by BFGS,

$$
H\_{\mathrm{new}} = H + \frac{\Delta\mathbf g\\,\Delta\mathbf g^{\mathsf T}}{\Delta\mathbf x\cdot\Delta\mathbf g} - \frac{H\Delta\mathbf x\\,(H\Delta\mathbf x)^{\mathsf T}}{\Delta\mathbf x\cdot H\Delta\mathbf x}
$$

and the update is obtained from the lowest eigenvalue \\( \lambda \\) of the augmented Hessian and its eigenvector:

$$
\begin{pmatrix} H & \tilde{\mathbf g} \\\\ \tilde{\mathbf g}^{\mathsf T} & 0\end{pmatrix}
\begin{pmatrix} \Delta\mathbf x \\\\ 1\end{pmatrix} = \lambda\begin{pmatrix} \Delta\mathbf x \\\\ 1\end{pmatrix},\qquad
\mathbf x\_{n+1} = \tilde{\mathbf x} + \Delta\mathbf x
$$

An equivalent form is \\( \Delta\mathbf x = -(H - \lambda I)^{-1}\tilde{\mathbf g} \\), \\( \lambda = \tilde{\mathbf g}\cdot\Delta\mathbf x < \lambda\_{\min}(H) \\) (`tools/check_geometry_md.py`).
The eigenvector is normalized by dividing by its last component.

[Convention]

- The initial Hessian (section 4.9) is rebuilt at steps 1, 31, 61, … of the quasi-Newton stage (\\( k \bmod 30 = 1 \\)).
- The BFGS denominators \\( \Delta\mathbf x\cdot\Delta\mathbf g \\) and \\( \Delta\mathbf x\cdot H\Delta\mathbf x \\), if smaller than \\( 10^{-10} \\), are replaced by \\( 10^{-10} \\) and **the update continues**.
  At steps where the curvature condition is violated (\\( \Delta\mathbf x\cdot\Delta\mathbf g < 0 \\)), a large positive term \\( 10^{10}\\,\Delta\mathbf g\Delta\mathbf g^{\mathsf T} \\) enters the Hessian.
  This differs from the description in manual section 18.2, "no update if the inner product is negative". A reimplementation makes skipping the update the default [Design].

### 4.7 `EF` [Public] [Convention]

[Public] Technical note 7, Eqs. 20–21 (Baker's eigenvector-following method used for minimization): the BFGS-updated \\( H \\) is diagonalized,
its eigenvalues are made positive, \\( H'^{-1} \\) is built, and \\( \mathbf x\_{n+1} = \tilde{\mathbf x} - f\\, H'^{-1}\tilde{\mathbf g} \\).

[Convention]

- The initial Hessian (section 4.9) is built only at the first step of the quasi-Newton stage.
- The BFGS update is performed only when \\( \Delta\mathbf x\cdot\Delta\mathbf g > 0 \\).
- Eigenvalue lower bound: eigenvalues smaller than 0.10 if \\( N \le 4 \\), otherwise 0.005 Hartree/bohr², are replaced by that value
  (technical note 7 gives 0.02 as the value of Ver. 3.3, but these are the current values).
- Prefactor \\( f \\): \\( f = 2.0 \\) at the first step of the quasi-Newton stage, and thereafter \\( f \leftarrow 0.95 f \\) each time `Utot` rises above that of the previous step.

### 4.8 Switching from steepest descent [Public] [Convention]

[Public] Manual section 18.2: steepest descent is used up to step `MD.Opt.StartDIIS`, and from there the method switches to quasi-Newton.
To avoid convergence to a saddle point, increase this value so that the switch happens after steepest descent has passed the inflection point.

[Convention]

| Keyword | Default (existing implementation) | Stated in manual | Meaning |
|---|---|---|---|
| `MD.Opt.DIIS.History` | 2 | 3 | upper bound on the GDIIS history length \\( p \\) |
| `MD.Opt.StartDIIS` | 1 | 5 | quasi-Newton from this step |
| `MD.Opt.EveryDIIS` | 200000 | 200 | when one stage has continued for this many steps, switch to the other |
| `MD.Opt.criterion` | 0.0003 | 0.0003 | Hartree/bohr |
| `MD.maxIter` | 1 | 1 | upper bound on the cumulative number of steps |

- While \\( n < \\) `MD.Opt.StartDIIS`, steepest descent (section 4.2). At \\( n = \\) `MD.Opt.StartDIIS` the quasi-Newton stage is entered,
  and the in-stage step number \\( k \\) is counted from 1 (the GDIIS history length and the Hessian rebuild are determined by this \\( k \\)).
- When the stage is switched, the steepest-descent prefactor returns to the initial value of section 4.2.
- The defaults stated in the manual (3rd column of the table) do not agree with the existing implementation. A reimplementation uses the values of the existing implementation as defaults [Design].

### 4.9 Initial Hessian [Public] [Convention]

[Public] Manual section 18.3: `MD.Opt.Init.Hessian` is `Schlegel` (default) or `iden`.
`Schlegel` uses **only the bond stretches** of Schlegel's model (bond stretch, bend, torsion).
Sources of the model: H. B. Schlegel, Theor. Chim. Acta **66**, 333 (1984); J. M. Wittbrodt and H. B. Schlegel, J. Mol. Struct. (Theochem) **398–399**, 55 (1997).

[Convention]

- `iden`: \\( H = I \\) (Hartree/bohr²).
- `Schlegel`: pair contributions are added to \\( H = 0.1\\, I \\). For atom \\( I \\) and its neighbors \\( (J, \mathbf L) \\)
  (\\( \lVert\boldsymbol\tau\_{IJ\mathbf L}\rVert \le r\_{c,I} + r\_{c,J} \\), the set of atom pairs of the [sparse block format](../design/sparse-block-format.md); \\( (J,\mathbf L) \ne (I, 0) \\)), with
  \\( r = \lVert\mathbf R\_I - \mathbf R\_J - \mathbf L\rVert \\), \\( \mathbf u = (\mathbf R\_I - \mathbf R\_J - \mathbf L)/r \\),

$$
k\_{IJ} = \frac{1.734}{(r - B\_{p(I)p(J)})^3}\ \ (\text{Hartree/bohr}^2),\qquad
H\_{I\alpha,I\beta} \mathrel{+}= k\_{IJ} u\_\alpha u\_\beta,\qquad
H\_{I\alpha,J\beta} \mathrel{+}= -k\_{IJ} u\_\alpha u\_\beta
$$

  (added for all ordered pairs \\( (I,J) \\), so the matrix is symmetric). \\( p(I) \\) is the period of the element: 1 for H, 2 for Z ≤ 10, 3 for Z ≤ 18, 4 for Z ≤ 36, 5 for Z ≤ 54, 6 for Z ≤ 86, 7 for Z ≤ 103,
  and 0 for the noble gases (He, Ne, Ar, Kr, Xe, Rn). \\( B \\) (bohr) is a symmetric table.

| \\( p \\) | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 |
| 1 | | −0.2573 | 0.3401 | 0.6937 | 0.7126 | 0.8335 | 0.9491 | 1.0 |
| 2 | | | 0.9652 | 1.2843 | 1.4625 | 1.6549 | 1.7190 | 2.0 |
| 3 | | | | 1.6925 | 1.8238 | 2.1164 | 2.3185 | 2.5 |
| 4 | | | | | 2.0203 | 2.2137 | 2.5206 | 2.7 |
| 5 | | | | | | 2.3718 | 2.5110 | 2.7 |
| 6 | | | | | | | 2.5200 | 2.7 |
| 7 | | | | | | | | 2.9 |

- The pair contributions satisfy \\( H\mathbf t = 0 \\) for a uniform translation \\( \mathbf t \\); only the diagonal term 0.1 gives curvature along the translation directions (`tools/check_geometry_md.py`).
  The contribution of a pair with a periodic image of the same atom (\\( J = I \\), \\( \mathbf L \ne 0 \\)) cancels to 0 in the diagonal block.
- The rows and columns of fixed components (section 4.10) are set to 0 and the diagonal to 1.
- Because the range of atom pairs is determined by the sum of PAO cutoff radii, distant pairs (\\( r \\) of 10 bohr or more) also have small contributions.
  For pairs with \\( r \le B \\), \\( k \\) becomes negative or diverges (does not occur at ordinary bond lengths). A reimplementation puts a lower bound on \\( r - B \\) [Design].
- There is a third option `FF` (builds the Hessian by fitting a Lennard-Jones-type force field to the forces), but it is not described in the public documents and is not treated in this chapter ([Open]).

### 4.10 Constraints `MD.Fixed.XYZ` [Public] [Convention]

[Public] Manual section 18.4: each line of the block is "atom number, x, y, z flags", where 1 is fixed and 0 is free. The default is all free.
Usable with all geometry optimizations and MD.

[Convention]

- As many lines as atoms are read (the number in the first column is skipped; lines are assigned to atoms in order). Too few lines is an input error.
- Geometry optimization: the gradient of fixed components is treated as 0 and the coordinates are not moved. They are excluded from the maximum-force and maximum-displacement tests.
- MD: fixed components do not have their positions updated, and their velocities are kept at 0 (except for the exception of section 6.1).
- In `RFC5`–`RFC7`, the flags act as constraints on the **components of the fractional coordinates** (along \\( \mathbf a\_1,\mathbf a\_2,\mathbf a\_3 \\)) (section 5.3).

### 4.11 Restart [Public] [Convention]

[Public] Manual sections 17.3 and 18.5: at each step `System.Name.dat#` (an input file containing the last structure) is written,
and the approximate Hessian is also saved at every step and taken over automatically when restarting from `System.Name.dat#`. The number of iterations does not change when the run is split.

[Convention] For geometry optimization (`Opt`, `DIIS`, `BFGS`, `RF`, `EF`), `geoopt.restart on` is added to `.dat#`.
When `geoopt.restart on`, the optimization state saved in `System.Name_rst/` (which stage, step count within the stage, prefactor, energy history, coordinate and gradient history, Hessian) is read at the first step.
The storage format may be chosen freely in this note [Design].

## 5. Lattice optimization

### 5.1 Variables and gradients [Public] [Convention]

[Public] Manual section 19: lattice optimization is available only for collinear spin (including +U) and cannot be combined with DFT-D2/D3.
Convergence is tested by comparing the maximum of the energy gradient with respect to atomic coordinates or lattice vectors with `MD.Opt.criterion`.
Manual section 19.4: when `MD.applied.pressure` (GPa, default 0) is given, the enthalpy \\( H = U + P\Omega \\) is minimized.

[Convention]

- The gradient with respect to the lattice vectors is the expression \\( \partial U/\partial A\_{ij} = \sum\_k (A^{-1})\_{ki}\Sigma\_{kj} \\) of [Total energy and forces](energy-force.md).
- In lattice optimization the stress calculation is enabled automatically, and the applied pressure enters the stress:
  \\( \Sigma\_{\gamma\gamma} \mathrel{+}= P\Omega f\_\gamma \\) (\\( f\_\gamma \\) is the \\( \gamma \\)-th entry of `MD.applied.pressure.flag`, default `1 1 1`).
  \\( U\_{pV} = P\Omega \\) (independent of the flags) is added to the energy, and `Enpy` = `Utot` + `UpV` is output (this resolves the [Open] item on `UpV` in [Total energy and forces](energy-force.md): it enters the stress only in lattice optimization).
- [Public] Manual section 19.4 explains that `MD.applied.pressure.flag` applies pressure per axis for an orthorhombic cell and, otherwise, in the direction perpendicular to the other two lattice vectors.
  Since the existing implementation adds to the Cartesian diagonal components \\( \Sigma\_{\gamma\gamma} \\) as above, this explanation is exact only when the lattice vectors lie along the coordinate axes.
- When the lattice is moved, atoms are moved by \\( \mathbf R\_I = \mathbf O + A^{\mathsf T}\mathbf s\_I \\), preserving the fractional coordinates \\( \mathbf s\_I \\) (measured from the origin \\( \mathbf O \\)).
- At each step of lattice optimization, the set of atom pairs (neighbor list) is rebuilt (same handling as ordinary geometry optimization and MD).

### 5.2 `OptC1`–`OptC7` (steepest descent) [Public] [Convention]

[Public] Table in manual section 19.1: `OptC1` has no constraint on the lattice and fixed internal coordinates; `OptC2` also fixes the angles and changes only the lengths;
`OptC3` further imposes \\( \lvert\mathbf a\_1\rvert = \lvert\mathbf a\_2\rvert = \lvert\mathbf a\_3\rvert \\); `OptC4` imposes \\( \lvert\mathbf a\_1\rvert = \lvert\mathbf a\_2\rvert \ne \lvert\mathbf a\_3\rvert \\);
`OptC5` optimizes the lattice and internal coordinates simultaneously; `OptC6` fixes \\( \mathbf a\_3 \\); `OptC7` fixes \\( \mathbf a\_2, \mathbf a\_3 \\).

[Convention] With \\( D\_{ij} = \partial H/\partial A\_{ij} \\) (0 for fixed components) and \\( s \\) the prefactor:

| Method | Convergence quantity \\( G\_{\max} \\) | Update |
|---|---|---|
| `OptC1` | \\( \max\_{ij}\lvert D\_{ij}\rvert \\) | \\( A\_{ij} \leftarrow A\_{ij} - c s D\_{ij} \\) (free components) |
| `OptC2` | \\( \max\_i \lvert d\_i\rvert \\), \\( d\_i = (\mathbf a\_i\cdot\mathbf D\_i)/\lVert\mathbf a\_i\rVert \\) | \\( \mathbf a\_i \leftarrow (1 - c s d\_i)\mathbf a\_i \\) (only \\( \mathbf a\_i \\) with all 3 components free) |
| `OptC3` | same as above (fixed flags are not consulted) | same as above, then the 3 lengths are set to their average |
| `OptC4` | same as above | same as above, then the lengths of \\( \mathbf a\_1, \mathbf a\_2 \\) are set to their average |
| `OptC5`–`OptC7` | \\( \max(\max\_{ij}\lvert D\_{ij}\rvert,\ \max\_{I\alpha}\lvert g\_{I\alpha}\rvert) \\) | Lattice as in `OptC1`. Atoms are moved with the new lattice preserving \\( \mathbf s\_I \\), then \\( \mathbf R\_I \leftarrow \mathbf R\_I - 0.2\\, c s\\, \mathbf g\_I \\) |

- `OptC6` and `OptC7` are `OptC5` with all components of \\( \mathbf a\_3 \\) (and \\( \mathbf a\_2 \\)) additionally fixed.
- The prefactor \\( s \\) has the same initial value and range as in section 4.2, but \\( H = U + P\Omega \\) is used as the energy and the adjustment rules differ:
  if \\( H\_n > H\_{n-1} \\) and \\( n \bmod 3 = 1 \\), then \\( s \leftarrow s/5 \\); if \\( H\_{n-1} < H\_{n-2} \\) and \\( H\_n < H\_{n-1} \\) and \\( n \bmod 4 = 1 \\), then \\( s \leftarrow 2.5 s \\).
- \\( c \\) is the reduction factor \\( 0.5/(sG\_{\max}) \\) applied when \\( s G\_{\max} \\) exceeds 0.5 (same idea as section 4.1; lattice components are also compared in bohr).
- The `OptC2` update moves \\( \mathbf a\_i \\) by \\( \mathbf a\_i \\) itself multiplied by \\( c s d\_i \\), not by the unit vector along \\( \mathbf a\_i \\), so
  the change in length is \\( c s d\_i\lVert\mathbf a\_i\rVert \\) (dimensionally inconsistent, but a convention of the existing implementation).
- In `OptC3` and `OptC4`, the length-equalizing operation is performed even if the initial structure does not satisfy the condition (there is no input check).

### 5.3 `RFC5`–`RFC7` [Public] [Convention]

[Public] Manual section 19.1: the lattice and internal coordinates are optimized simultaneously by a hybrid of RF and DIIS (Hessian updated by BFGS).
For arbitrary constraints, the manual recommends combining `RFC5` with `MD.Fixed.Cell.Vectors`.

[Convention]

- The variables are the fractional coordinates \\( \mathbf s\_I \\) (\\( 3N \\) of them) and the Cartesian components of the lattice vectors that are moved (3 vectors for `RFC5`, \\( \mathbf a\_1,\mathbf a\_2 \\) for `RFC6`, \\( \mathbf a\_1 \\) for `RFC7`).
  The gradients are \\( \partial H/\partial s\_{Ii} = \sum\_\alpha g\_{I\alpha} A\_{i\alpha} \\) and \\( D\_{ij} \\).
- The update applies the GDIIS of section 4.3 and the RF step of section 4.6 to this set of variables (the treatment of the BFGS denominators is also the same as in section 4.6).
  The maximum-displacement test compares the maximum change of the updated Cartesian coordinates and lattice components (bohr) with 0.5.
- The convergence quantity is \\( \max(\max\_{ij}\lvert D\_{ij}\rvert,\ \max\_{I\alpha}\lvert g\_{I\alpha}\rvert) \\) (tested with the Cartesian gradient).
- Before entering the quasi-Newton stage (section 4.8), `OptC5`-type steepest descent is used.
- Initial Hessian (built only at the first step of the quasi-Newton stage):
  - `iden`: identity matrix (dimension \\( 3N + 9 \\)).
  - `Schlegel`: start with diagonal 0.1 for the atomic part and diagonal 0.05 for the lattice part, add the pair contributions of section 4.9, then transform the atomic part
    to fractional coordinates (\\( H^{(s)}\_{Ii,Jj} = \sum\_{\alpha\beta} A\_{i\alpha}H\_{I\alpha,J\beta}A\_{j\beta} \\)),
    and build the lattice-lattice and atom-lattice parts from the second derivative of the same bond-stretch model with respect to strain.
    The details of these cross terms affect only the convergence path, so this chapter does not specify them ([Open]).
- For fixed fractional-coordinate components, the Hessian rows and columns are set to 0 and the diagonal to 200.
- Lattice components fixed by `MD.Fixed.Cell.Vectors` have their gradients set to 0.

### 5.4 `MD.Fixed.Cell.Vectors` [Public] [Convention]

[Public] Manual section 19.3: a 3×3 array of flags; row \\( i \\) is \\( \mathbf a\_i \\), columns are the x, y, z components. 1 is fixed, 0 is free.
Usable with `OptC1`, `OptC2`, `OptC5`, `RFC5`.

[Convention] If this block is present, it takes precedence over the constraints determined by `OptC6`, `OptC7`, `RFC6`, `RFC7`.
In `OptC2`–`OptC4`, since the update is per vector, only vectors with all 3 components free move (the length-equalizing operation of `OptC3` and `OptC4` does not consult the flags).

### 5.5 `EvsLC` and `DF` [Convention]

- `EvsLC`: after each step, using the initial lattice \\( A^{(0)} \\), set \\( \mathbf a\_i \leftarrow \mathbf a\_i \pm (\delta/100)\\,\mathbf a\_i^{(0)} \\).
  \\( \delta \\) is `MD.EvsLC.Step` (%, default 0.4); the sign is given by `MD.EvsLC.flag` (3 integers: 1 to expand, −1 to shrink, 0 to fix; default `1 1 1`).
  Fractional coordinates are preserved. The lattice vectors (Å) and `Utot` of each step are appended to `System.Name.EvsLC`.
- `DF` (for the delta factor): `MD.maxIter` must be 7, otherwise it is an input error. Step 1 uses the input lattice; steps 2–7 use lattices
  with the volume scaled by 0.94, 0.96, 0.98, 1.02, 1.04, 1.06 (the lattice scaled isotropically by the 1/3 power of the volume ratio),
  and the volume per atom (Å³) and `Utot` per atom (eV) of each step are appended to `System.Name.DF`.

## 6. Molecular dynamics

### 6.1 Initial velocities [Public] [Convention]

[Public] Manual section 20.8: each line of the `MD.Init.Velocity` block is "atom number, \\( v\_x, v\_y, v\_z \\) (m/s)"; it can be combined with `MD.Fixed.XYZ`.

[Convention]

- If `MD.Init.Velocity` is present, it is taken as the current velocity (fixed components set to 0).
  `MD.Init.Velocity.Prev` (same format) is the velocity one step earlier (half step), written by the restart file.
- If neither is present and `MD.Type` has a thermostat (`NVT_*`, `NPT_*`), initial velocities are generated randomly:
  1. For each atom, \\( \lVert\mathbf v\_I\rVert = \sqrt{3k\_BT\_1/m\_I} \\) (\\( T\_1 \\) is the temperature in the first line of `MD.TempControl`), with direction given by a normalized uniform random vector in \\( [-\tfrac12,\tfrac12]^3 \\).
  2. Set the total momentum to 0: \\( \mathbf v\_I \leftarrow \mathbf v\_I - \bar{\mathbf P}/m\_I \\), \\( \bar{\mathbf P} = \frac1N\sum\_J m\_J\mathbf v\_J \\).
  3. Scale everything by a constant so that \\( T\_A \\) (section 2; \\( K \\) excluding fixed components) equals exactly \\( T\_1 \\).
  - The random sequence is the C library pseudo-random generator (seed 1) and depends on the environment. For comparison, use inputs that specify `MD.Init.Velocity` explicitly [Design].
  - In this case, the velocities of fixed components are not set to 0 (they do not enter \\( K \\), but they enter the temperature in the `NVT_VS` family; section 11).
- Otherwise (`NVE`, etc.), the initial velocity is 0.

### 6.2 `NVE` (velocity Verlet) [Public] [Convention]

[Public] Manual section 20.1: microcanonical (NVE) MD. Conservation is checked with columns 14–17 of `.ene` (`Ukc`, `Utot`, `Utot+Ukc`, Fermi level).

[Convention] Velocity Verlet. With \\( \mathbf a\_n = -\mathbf g\_n/m \\) (after the correction of section 1), at step \\( n \\):

$$
\mathbf v\_n = \mathbf v\_{n-1/2} + \tfrac12\mathbf a\_n\Delta t\ \ (n\ge 2),\qquad
\mathbf v\_{n+1/2} = \mathbf v\_n + \tfrac12\mathbf a\_n\Delta t,\qquad
\mathbf x\_{n+1} = \mathbf x\_n + \mathbf v\_{n+1/2}\Delta t
$$

- At \\( n = 1 \\), the input velocity is taken as \\( \mathbf v\_1 \\). \\( K \\) and the temperature \\( T\_A \\) are computed from \\( \mathbf v\_n \\) and written to `.ene`.
- The velocity held at the end of the step is \\( \mathbf v\_{n+1/2} \\) (this half-step velocity appears in the `.md` of the next step).
- Fixed components have neither position nor velocity updated.
- `tools/check_geometry_md.py` checks that for a harmonic oscillator the energy error is bounded and \\( O(\Delta t^2) \\).

### 6.3 Temperature control table `MD.TempControl` [Public] [Convention]

[Public] Manual sections 20.2 and 20.3. The first number of the block is the number of lines \\( L \\); the format of the following lines depends on the method.

- Velocity scaling (`NVT_VS`, `NVT_VS2`, `NVT_VS4`, `NPT_VS_PR`, `NPT_VS_WV`): each line is "\\( N\_i \\) (last step of the interval), \\( I\_i \\) (scaling interval), \\( T\_i \\) (K), \\( \alpha\_i \\)".
- Nosé–Hoover and Langevin (`NVT_NH`, `NVT_Langevin`, `NPT_NH_PR`, `NPT_NH_WV`): each line is "\\( S\_i \\) (step), \\( T\_i \\) (K)".
- In all cases, the target temperature between lines is linearly interpolated in the step number.

[Convention] How the interpolation intervals are taken differs between methods. \\( (N\_0, T\_0) = (0, T\_1) \\) (or \\( (S\_0, T\_0) = (0, T\_1) \\)) is prepended.

| Method | Interval containing step \\( n \\) | Outside the intervals |
|---|---|---|
| Velocity scaling (NVT, NPT) | \\( N\_{i-1} < n \le N\_i \\) | no scaling (for the `NVT_VS` family only the integration of section 6.4; what the manual calls the transition to NVE) |
| `NVT_NH` | \\( S\_i \le n < S\_{i+1} \\) (\\( i = 1,\dots,L-1 \\)) | the thermostat is switched off and the system advances by velocity Verlet (section 6.2) |
| `NVT_Langevin` | same as above | the target temperature stays at the last set value; the thermostat keeps acting |
| `NPT_NH_*` | \\( S\_{i-1} < n \le S\_i \\) | the target temperature stays at the last value; the thermostat keeps acting |

### 6.4 Velocity scaling `NVT_VS`, `NVT_VS2`, `NVT_VS4` [Public] [Convention]

[Public] Manual section 20.2: in interval \\( i \\), the velocities are multiplied by \\( s \\) every \\( I\_i \\) steps,
so that the temperature after scaling is \\( T\_0 + (T - T\_0)\alpha\_i \\) (\\( \alpha = 0 \\) matches the target temperature; \\( \alpha = 1 \\) does nothing):

$$
s = \sqrt{\frac{T\_0 + (T - T\_0)\alpha\_i}{T}}
$$

(The equation in the HTML manual has a broken square root and is hard to read, but the factor that brings the temperature to this value includes the square root. `tools/check_geometry_md.py`.)
`NVT_VS2` scales each atom by its per-atom temperature, and `NVT_VS4` (manual section 20.4) scales each group by its per-group temperature.
Atom groups are specified by `MD.num.AtomGroup` (number of groups) and `MD.AtomGroup` (each line "atom number, group number"; group numbers consecutive from 1),
and all groups follow the same `MD.TempControl`.

[Convention]

- Form of the integration (common to the 3 methods): with \\( \mathbf u \\) the held velocity,
  - \\( n = 1 \\): \\( \mathbf x\_2 = \mathbf x\_1 + \mathbf u\_1\Delta t + \tfrac12\mathbf a\_1\Delta t^2 \\), \\( \mathbf u \leftarrow \mathbf u\_1 + \mathbf a\_1\Delta t \\). No scaling.
  - \\( n \ge 2 \\): \\( \mathbf u \leftarrow \mathbf u + \mathbf a\_n\Delta t \\) → subtract the velocity average (below) → compute the temperature → scaling test and \\( \mathbf u \leftarrow s\mathbf u \\) →
    \\( \mathbf x\_{n+1} = \mathbf x\_n + \mathbf u\Delta t + \tfrac12\mathbf a\_n\Delta t^2 \\).
  - This is **not** velocity Verlet. Without scaling, \\( \mathbf x\_{n+1} - 2\mathbf x\_n + \mathbf x\_{n-1} = \tfrac12(3\mathbf a\_n - \mathbf a\_{n-1})\Delta t^2 \\),
    and for a harmonic oscillator \\( \omega \\) the characteristic roots have absolute value \\( \sqrt{1 - (\omega\Delta t)^2/2} \\), i.e. energy is lost at a rate of \\( (\omega\Delta t)^2/2 \\) per step
    (0.45% per step for an oscillation with period 33 fs and \\( \Delta t = 0.5 \\) fs). Velocity scaling compensates for this dissipation. Checked in `tools/check_geometry_md.py`.
- Velocity average: for the free components, the simple average of the velocity (no mass weighting) is subtracted per Cartesian component. This does not coincide with zero total momentum (when masses differ).
- Temperature: `NVT_VS` uses \\( T\_A \\), but with \\( K \\) summed over **all components** (including fixed components).
  `NVT_VS2` uses the per-atom \\( T\_I = K\_I/(\tfrac32k\_B) \\) (computed from the velocity before subtracting the average), and `NVT_VS4` the per-group \\( T\_g = K\_g/(\tfrac32 k\_B N\_g) \\) (after subtracting the average).
- Scaling is performed only at steps with \\( (n - N\_{i-1}) \bmod I\_i = 0 \\) in interval \\( i \\). The target temperature is the interpolated value at that step.
- At the final step (\\( n = \\) `MD.maxIter`), neither scaling nor the position update is performed.

### 6.5 Nosé–Hoover `NVT_NH` [Public] [Convention]

[Public] Manual section 20.3: thermostat mass `NH.Mass.HeatBath` (default 20.0, in units of unified atomic mass unit × bohr²). The Hamiltonian of the extended system is output to `.ene` as the conserved quantity.
The equations of motion are of the form of Nosé (Mol. Phys. **52**, 255 (1984)) and Hoover (Phys. Rev. A **31**, 1695 (1985)),

$$
\dot{\mathbf R}\_I = \mathbf v\_I,\qquad
\dot{\mathbf v}\_I = \frac{\mathbf F\_I}{m\_I} - \zeta\mathbf v\_I,\qquad
\dot\zeta = \frac{2K - gk\_BT\_0}{Q},\qquad
\dot\eta = \zeta
$$

and \\( H' = U + K + \tfrac12 Q\zeta^2 + gk\_BT\_0\eta \\) is conserved (for constant \\( T\_0 \\)).

[Convention]

- \\( g = 3N \\) (counted with the total number of atoms even if there are fixed components); \\( Q \\) is `NH.Mass.HeatBath` converted to \\( m\_e \\)·bohr².
- Discretization (velocity-Verlet type). When the forces of step \\( n \\) are available:
  1. **Correction** (\\( n \ge 2 \\)): from the predicted value \\( \zeta' \\) of the previous step and the half-step velocity \\( \mathbf v\_{n-1/2} \\), the simultaneous equations
     \\( \mathbf v\_n = \mathbf v\_{n-1/2} + \tfrac12\Delta t\\,(\mathbf a\_n - \zeta\_n\mathbf v\_n) \\),
     \\( \zeta\_n = \zeta' + \Delta t\\,(K(\mathbf v\_n) - \tfrac12 gk\_BT\_0)/Q \\)
     are solved by Newton's method (at most 21 iterations; the details of the convergence test matter only for bitwise agreement and are not specified).
  2. Subtract the velocity average (same as section 6.4), compute \\( K \\), \\( T\_A \\), \\( H' \\), and write them to `.ene`.
  3. **Prediction** (except at the final step): \\( \mathbf v\_{n+1/2} = \mathbf v\_n + \tfrac12\Delta t\\,(\mathbf a\_n - \zeta\_n\mathbf v\_n) \\), \\( \mathbf x\_{n+1} = \mathbf x\_n + \mathbf v\_{n+1/2}\Delta t \\),
     \\( \zeta' = \zeta\_n + \Delta t\\,(K(\mathbf v\_n) - \tfrac12 gk\_BT\_0)/Q \\), \\( \eta\_{n+1} = \eta\_n + \zeta'\Delta t \\).
  - Combining 1 and 3 gives \\( \zeta\_{n+1} = \zeta\_n + \tfrac{\Delta t}{2}[\dot\zeta\_n + \dot\zeta\_{n+1}] \\) (trapezoidal rule), consistent with \\( \dot\zeta = (2K - gk\_BT\_0)/Q \\) of the continuous equations above.
    `tools/check_geometry_md.py` checks that this discretization conserves \\( H' \\) to \\( O(\Delta t^2) \\).
- At \\( n = 1 \\), start from \\( \zeta = \eta = 0 \\) (the input `NH.czeta` etc. are used only on restart, i.e. for cumulative \\( n \ge 2 \\)).
- Outside the intervals (section 6.3), \\( \zeta \\) is not used and the system advances by velocity Verlet. The \\( H' \\) in `.ene` is computed with the \\( \zeta, \eta \\) at that time.
- Columns 20, 21, 22 of `.ene` are \\( \zeta \\), the predicted value \\( \zeta' \\), and \\( H' \\).

### 6.6 Langevin `NVT_Langevin` [Convention]

There is no description in the public documents ([Convention] only). The equation of motion is the standard Langevin equation

$$
m\_I\dot{\mathbf v}\_I = \mathbf F\_I - m\_I\gamma\mathbf v\_I + \mathbf R\_I(t),\qquad
\langle R\_{I\alpha}(t)R\_{J\beta}(t')\rangle = 2m\_I\gamma k\_BT\_0\\,\delta\_{IJ}\delta\_{\alpha\beta}\delta(t-t')
$$

\\( \gamma \\) is `Langevin.Friction.Factor` (inverse of the atomic time unit, default 0.001, i.e. about 41 ps⁻¹).

- At each step, a Gaussian random number with standard deviation \\( \sqrt{2m\_I\gamma k\_BT\_0/\Delta t} \\) (one branch of Box–Muller) is added to the force on the free components. Below, \\( \tilde{\mathbf a} \\) is the acceleration including this term.
- Discretization:
  - \\( n = 1 \\): \\( \mathbf x\_2 = \mathbf x\_1 + \mathbf u\Delta t + \tfrac12\Delta t^2(\tilde{\mathbf a}\_1 - \gamma\mathbf u) \\), \\( \mathbf u \leftarrow \mathbf u + \tilde{\mathbf a}\_1\Delta t \\).
  - \\( n \ge 2 \\): \\( \mathbf u \leftarrow \mathbf u + \tfrac12(\tilde{\mathbf a}\_n + \tilde{\mathbf a}\_{n-1})\Delta t - \gamma(\mathbf x\_n - \mathbf x\_{n-1}) \\),
    \\( \mathbf x\_{n+1} = \mathbf x\_n + (\mathbf u\Delta t + \tfrac12\tilde{\mathbf a}\_n\Delta t^2)/(1 + \tfrac12\gamma\Delta t) \\).
  - The position is updated even at the final step (unlike the other methods). The velocity average is not subtracted.
- \\( K \\) and \\( T\_A \\) are computed with the updated \\( \mathbf u \\).
- The random numbers continue the same C library pseudo-random sequence as in section 6.1. Trajectories cannot be compared, so statistics (mean temperature, etc.) are compared [Design].

### 6.7 NPT `NPT_VS_PR`, `NPT_VS_WV`, `NPT_NH_PR`, `NPT_NH_WV` [Public] [Convention]

[Public] Manual sections 20.6 and 20.7: the thermostat is velocity scaling (VS) or Nosé–Hoover (NH);
the barostat is Parrinello–Rahman (PR, J. Appl. Phys. **52**, 7182 (1981)) or Wentzcovitch (WV, Phys. Rev. B **44**, 2358 (1991)).
Keywords:

| Keyword | Default | Meaning |
|---|---|---|
| `MD.TempControl` | — | section 6.3 |
| `MD.CellPressureControl` | — | after the number of lines, lines "step, pressure (GPa)". Linearly interpolated in between. Positive pressure compresses |
| `NPT.Mass.Barostat` | automatic | barostat mass \\( W \\) (unified atomic mass unit) |
| `NPT.LatticeRestriction` | `none` | `none`, `cubic`, `ortho` |
| `MD.TempTolerance` | 100 K | no scaling if \\( \lvert T - T\_0\rvert \\) is below this (VS) |
| `MD.applied.pressure.flag` | `1 1 1` | components to which pressure is applied |
| `NH.Mass.HeatBath` | 20 | NH thermostat mass (same units as section 6.5) |

Constraints (`MD.Fixed.XYZ`) cannot be combined with NPT (manual section 20.6).

[Convention] The variables are the fractional coordinates \\( \mathbf s\_I \\), their time derivatives \\( \dot{\mathbf s}\_I = A^{-\mathsf T}\mathbf v\_I \\), the lattice \\( A \\), and its velocity \\( \dot A \\)
(the lattice velocity is stored in the `Atoms.Unitvectors.Velocity` block of the restart file and is initially 0).
With \\( G = AA^{\mathsf T} \\) (\\( G\_{ij} = \mathbf a\_i\cdot\mathbf a\_j \\)), the atomic equation of motion is

$$
\ddot{\mathbf s}\_I = -\frac{A^{-\mathsf T}\mathbf g\_I}{m\_I} - G^{-1}\dot G\\,\dot{\mathbf s}\_I\ \ [-\ \zeta\dot{\mathbf s}\_I\ \text{(NH)}]
$$

The lattice equation of motion (PR) uses the internal pressure tensor

$$
\Pi\_{\alpha\beta} = \frac1\Omega\sum\_I m\_I v\_{I\alpha}v\_{I\beta} - \frac{\Sigma\_{\alpha\beta}}{\Omega},\qquad \mathbf v\_I = A^{\mathsf T}\dot{\mathbf s}\_I
$$

and reads

$$
W\ddot A\_{i\alpha} = \sum\_\beta\left(\Pi\_{\alpha\beta} - P f\_\alpha\delta\_{\alpha\beta}\right)\sigma\_{i\beta}
$$

(\\( \sigma\_{i\beta} \\) are the elements of \\( S \\); \\( \partial\Omega/\partial A\_{i\beta} = \sigma\_{i\beta} \\), `tools/check_geometry_md.py`).
For `cubic`, the off-diagonal elements of the matrix in parentheses on the right-hand side are set to 0 and the diagonal to its average; for `ortho`, the off-diagonal elements are set to 0.

- Time integration (one step; semi-implicit Euler):
  1. Compute \\( \ddot{\mathbf s} \\) from the current quantities, and \\( \dot{\mathbf s} \leftarrow \dot{\mathbf s} + \ddot{\mathbf s}\Delta t \\).
  2. Compute \\( K \\) and \\( T\_A \\) with the new \\( \dot{\mathbf s} \\) and the current \\( A \\).
     VS: within an interval, if \\( n = 1 \\), or if \\( (n - N\_{i-1}) \bmod I\_i = 0 \\) and \\( \lvert T\_0 - T\rvert > \\) `MD.TempTolerance`, multiply \\( \dot{\mathbf s} \\) by \\( s \\) (factor 1 if \\( T < 10^{-8}T\_0 \\)).
     NH: \\( \zeta \leftarrow \zeta + Nk\_B(T\_A - T\_0)\Delta t/Q \\), \\( \eta \leftarrow \eta + \zeta\Delta t \\) (starting from \\( \zeta = \eta = 0 \\) at \\( n = 1 \\)).
  3. \\( \mathbf s \leftarrow \mathbf s + \dot{\mathbf s}\Delta t \\).
  4. Compute \\( \ddot A \\) with the new \\( \dot{\mathbf s} \\), the current \\( A \\), and the stress; \\( \dot A \leftarrow \dot A + \ddot A\Delta t \\), \\( A \leftarrow A + \dot A\Delta t \\).
  5. \\( \mathbf R\_I = \mathbf O + A^{\mathsf T}\mathbf s\_I \\), \\( \mathbf v\_I = A^{\mathsf T}\dot{\mathbf s}\_I \\).
- External pressure \\( P \\): \\( (0, P\_1) \\) is prepended to the lines \\( (S\_i, P\_i) \\) of `MD.CellPressureControl`, and the value is interpolated for \\( S\_{i-1} < n \le S\_i \\).
  The value beyond the step of the last line is undefined ([Convention]; a reimplementation keeps the last value [Design]). `MD.applied.pressure` is not used in NPT.
- The NH thermostat is \\( \dot\zeta = Nk\_B(T\_A - T\_0)/Q = (2K - 3Nk\_BT\_0)/(3Q) \\). This corresponds to the form of section 6.5 with mass \\( 3Q \\),
  and the conserved quantity contains \\( \tfrac12(3Q)\zeta^2 \\). The value output by the existing implementation uses \\( \tfrac12 Q\zeta^2 + 3Nk\_BT\_0\eta \\), so it is not strictly conserved (item 5 of `tools/check_geometry_md.py`).
- Default \\( W \\): for PR, \\( W = \frac{3}{4\pi^2}\sum\_I m\_I \\); for WV, \\( W = \frac{3}{4\pi^2}\sum\_I m\_I\times\frac{0.01}{\Omega^{2/3}} \\) (\\( \Omega \\) is the initial volume, bohr³).
  [Public] Eq. 27 in manual section 20.6 agrees with the PR expression. The expression for the WV default in section 20.7 is garbled and cannot be read.
- The WV lattice equation of motion is the Euler–Lagrange equation of the Lagrangian whose lattice kinetic energy is \\( \tfrac12 W\sum\_{ij}(\dot{\mathbf a}\_i\cdot\dot{\mathbf a}\_j)F\_{ij} \\), \\( F = SS^{\mathsf T} \\) (current lattice),

$$
W\sum\_j\left(F\_{ij}\ddot{\mathbf a}\_j + \dot F\_{ij}\dot{\mathbf a}\_j\right) = \mathbf f\_i + \tfrac12 W\sum\_{jk}(\dot{\mathbf a}\_j\cdot\dot{\mathbf a}\_k)\\,\frac{\partial F\_{jk}}{\partial\mathbf a\_i},\qquad
f\_{i\alpha} = \sum\_\beta\left(\Pi\_{\alpha\beta} - Pf\_\alpha\delta\_{\alpha\beta}\right)\sigma\_{i\beta}
$$

  solved for \\( \ddot{\mathbf a}\_i \\) by multiplying by \\( F^{-1} = G/\Omega^2 \\) (\\( \mathbf f\_i \\) is the same as the PR right-hand side;
  \\( \partial\boldsymbol\sigma\_j/\partial a\_{i\alpha} \\) is \\( (\sigma\_{i\alpha}\boldsymbol\sigma\_j - \sigma\_{j\alpha}\boldsymbol\sigma\_i)/\Omega \\)).
  The lattice restrictions (`cubic`, `ortho`) are applied to the matrix in parentheses in \\( \mathbf f \\). [Convention] The existing implementation computes \\( \dot F\_{ij} \\) as
  \\( 2\dot{\boldsymbol\sigma}\_i\cdot\boldsymbol\sigma\_j \\) (not symmetrized).
- Output: one line per step in `System.Name.npt.ene`: "step, time \\( n\Delta t \\) (fs), volume (bohr³), \\( T\_A \\) (K), internal pressure \\( -\mathrm{tr}\\,\Sigma/(3\Omega) \\) (GPa, factor 29.42×10³),
  `Utot`, `Ukc`, `Ubaro`, `Wbaro` = \\( P\Omega \\), their sum (for NH the thermostat energy is also added)" (Hartree).
  For PR, `Ubaro` is \\( \tfrac12 W\sum\_{i}\lVert\dot{\mathbf a}\_i\rVert^2 \\). For WV, `Ubaro` is a quantity built from \\( F^{(0)} \\) of the initial lattice and does not agree with the kinetic energy above ([Open]).
  Two header lines (method name, time step, \\( W \\) (unified atomic mass unit)) are placed at the top of the file.
- In NPT, `.ene` is not written. The velocity columns of `.md` are 0 (section 10).

### 6.8 MD restart [Convention]

The following are written to `System.Name.dat#`; using it as input as is continues the calculation.

- `MD.Current.Iter` (that step number), `<Atoms.SpeciesAndCoordinates` (current coordinates, in the same units as the input), `scf.restart on`, `scf.fixed.grid` (origin).
- `<Atoms.UnitVectors` (current lattice); for NPT, `<Atoms.Unitvectors.Velocity` (\\( \dot A \\)), and for WV, `<NPT.WV.F0`.
- For MD, `<MD.Init.Velocity` (held velocity, m/s). For `NVE`, `NVT_NH`, `NVT_Langevin`, `NPT_NH_*`, additionally `<MD.Init.Velocity.Prev` (half-step velocity).
  For NH, `NH.R` (\\( \eta \\)), `NH.nzeta` (\\( \zeta' \\)), `NH.czeta` (\\( \zeta \\)).

## 7. Wrapping atoms under periodic boundary conditions [Convention]

- After the update of each step of molecular dynamics (`NVE`, `NVT_*`, `NPT_*`), atoms whose fractional coordinates measured from the origin \\( \mathbf O \\),
  \\( s\_i = (\mathbf R - \mathbf O)\cdot\mathbf b\_i/(2\pi) \\) (\\( \mathbf b\_i \\) are the reciprocal lattice vectors), are out of range are moved back:
  if \\( s\_i > 1 \\), \\( s\_i \leftarrow s\_i - \lfloor s\_i\rfloor \\); if \\( s\_i < 0 \\), \\( s\_i \leftarrow s\_i + \lfloor -s\_i\rfloor + 1 \\).
  Exactly \\( s\_i = 1 \\) is left as is (exactly \\( s\_i = -1 \\) becomes 1).
- For atoms moved back, the coordinate history used for density extrapolation and the GDIIS coordinate history are shifted by the same lattice vector, so that the extrapolation and history are continuous
  ([Density matrix and charge mixing](density-mixing.md) section 7). Whether atoms were moved back is passed to the SCF of the next step and affects how it restarts from the quantities of the previous step (details [Open]).
- Geometry optimization and lattice optimization do not wrap. NEB does not wrap either (section 8).
- For isolated systems (`Cluster`, etc.), the same processing is applied with respect to the box. When a molecule straddles the box boundary it gets split, so a reimplementation is recommended not to wrap for isolated systems [Design].
- When fractional coordinates (`FRAC`) are used in the input, they are also brought into \\( [0, 1) \\) when read (same rule).

## 8. NEB

### 8.1 Input [Public] [Convention]

[Public] Manual section 51: after optimizing the reactant and product structures separately, `MD.Type NEB` finds the minimum energy path connecting them.
The same unit cell, cutoff, basis, pseudopotentials, and electronic temperature are used. For periodic systems, give coordinates in `Ang` or `AU` and do not use `FRAC` (images become discontinuous).

| Keyword | Default (existing implementation) | Stated in manual | Meaning |
|---|---|---|---|
| `MD.NEB.Number.Images` | 10 | 10 | number of images \\( P \\) excluding the end points |
| `MD.NEB.Spring.Const` | 0.4 | 0.1 | spring constant \\( k \\) (Hartree/bohr²) |
| `MD.NEB.Parallel.Number` | 0 (automatic) | — | number of groups for image parallelization (section 8.4) |
| `MD.Opt.DIIS.History` | 2 [Open] | 7 | DIIS history |
| `MD.Opt.StartDIIS` | 1 | 5 | DIIS + BFGS from this step |
| `MD.Opt.criterion` | 0.0003 | 0.0003 | Hartree/bohr |
| `MD.maxIter` | 1 | 1 | upper bound on optimization iterations |
| `MD.Fixed.XYZ` | none | — | components to fix |

[Convention]

- The reactant is given by `<Atoms.SpeciesAndCoordinates` and the product by `<NEB.Atoms.SpeciesAndCoordinates` (same format).
- The initial path is linear interpolation in Cartesian coordinates, \\( \mathbf x\_p = \mathbf x\_0 + \frac{p}{P+1}(\mathbf x\_{P+1} - \mathbf x\_0) \\) (\\( p = 1,\dots,P \\)).
  End points given in `FRAC` are brought into \\( [0,1) \\) when read, so a motion crossing a periodic boundary becomes discontinuous (the reason for the manual's caution).
  The minimum-image convention is not used.
- [Public] Manual section 51.5: with `scf.restart on`, the coordinates of each image can be given by `<NEB1.Atoms.SpeciesAndCoordinates` … `<NEBP.Atoms.SpeciesAndCoordinates`.
- Each image (and the two end points, only the first time) is run as an independent ordinary SCF calculation, producing `System.Name.dat_p` (input) and `System.Name_p.out` (\\( p = 0,\dots,P+1 \\)).

### 8.2 Tangent and forces [Public] [Convention]

[Public] The NEB technical note (T. Ozaki, 2011) and manual section 51.1: the improved tangent (Eqs. 8–11), perpendicular force (Eq. 4), and parallel spring force (Eq. 12)
of G. Henkelman and H. Jónsson, J. Chem. Phys. **113**, 9978 (2000) are used. There is no climbing image.

With \\( V\_p \\) the energy of image \\( p \\), \\( \boldsymbol\tau^+ = \mathbf x\_{p+1} - \mathbf x\_p \\), \\( \boldsymbol\tau^- = \mathbf x\_p - \mathbf x\_{p-1} \\),

$$
\boldsymbol\tau\_p = \begin{cases}
\boldsymbol\tau^+ & V\_{p-1} < V\_p < V\_{p+1} \\\\
\boldsymbol\tau^- & V\_{p+1} < V\_p < V\_{p-1} \\\\
\Delta V^{\max}\boldsymbol\tau^+ + \Delta V^{\min}\boldsymbol\tau^- & \text{extremum, } V\_{p-1} \le V\_{p+1} \\\\
\Delta V^{\min}\boldsymbol\tau^+ + \Delta V^{\max}\boldsymbol\tau^- & \text{extremum, } V\_{p-1} > V\_{p+1}
\end{cases}
$$

(\\( \Delta V^{\max/\min} = \max/\min(\lvert V\_{p+1} - V\_p\rvert, \lvert V\_{p-1} - V\_p\rvert) \\)) is normalized to give \\( \hat{\boldsymbol\tau}\_p \\), and
the gradient used for the optimization is

$$
\mathbf G\_p = \mathbf g\_p - (\mathbf g\_p\cdot\hat{\boldsymbol\tau}\_p)\hat{\boldsymbol\tau}\_p - k\left(\lVert\boldsymbol\tau^+\rVert - \lVert\boldsymbol\tau^-\rVert\right)\hat{\boldsymbol\tau}\_p
$$

(\\( \mathbf F^{\mathrm{NEB}}\_p = -\mathbf G\_p \\)).

[Convention]

- An "extremum" is the case \\( V\_{p+1} \ge V\_p \le V\_{p-1} \\) or \\( V\_{p+1} < V\_p > V\_{p-1} \\).
  Combinations of equalities such as \\( V\_p = V\_{p+1} > V\_{p-1} \\) match none of the cases; the tangent becomes 0 and the normalization diverges (section 11).
- \\( \mathbf G\_p \\) of fixed components is set to 0. The gradient is the SCF gradient of each image as is, without the zero-sum correction.

### 8.3 Optimization [Public] [Convention]

[Public] Manual section 51.3 and the technical note: a hybrid of DIIS and BFGS on the \\( 3NP \\)-dimensional vector collecting the coordinates of all images.

[Convention]

- Convergence test: \\( \max\_{p,I}\lVert\mathbf G\_{p,I}\rVert < \\) `MD.Opt.criterion` (per-atom vector magnitude).
- \\( n < \\) `MD.Opt.StartDIIS`: steepest descent \\( \mathbf x \leftarrow \mathbf x - s\mathbf G \\). The prefactor is restricted to \\( [0.005s\_0, 10s\_0] \\) with reference \\( s\_0 = (4\ \mathrm{fs})^2/1836.1526 \\),
  and initially \\( s = 1/(5\times0.529177249) \\). If \\( \lVert\mathbf G\rVert \\) (Euclidean norm of the whole vector) increased from the previous step and \\( n \bmod 4 = 1 \\), then \\( s/4 \\);
  if it decreased twice in a row and \\( n \bmod 4 = 1 \\), then \\( 1.2s \\). There is no maximum-displacement limit.
- \\( n \ge \\) `MD.Opt.StartDIIS`: the GDIIS of section 4.3 (with \\( \mathbf G \\) instead of \\( \mathbf g \\), all images together) builds \\( \tilde{\mathbf x} \\), \\( \tilde{\mathbf G} \\), and with
  the BFGS inverse Hessian of section 4.5 (\\( I \\) at step `MD.Opt.StartDIIS`, no curvature-condition check),
  \\( \mathbf x \leftarrow \tilde{\mathbf x} - f\\, H^{-1}\tilde{\mathbf G} \\). \\( f \\) is determined by the maximum force \\( F\_{\max} \\):

| \\( F\_{\max} \\) | \\( >0.7 \\) | \\( >0.3 \\) | \\( >0.1 \\) | \\( >0.05 \\) | \\( >0.02 \\) | \\( >0.01 \\) | \\( >0.005 \\) | \\( >0.002 \\) | \\( >0.001 \\) | otherwise |
|---|---|---|---|---|---|---|---|---|---|---|
| \\( f \\) | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 |

- If the per-component maximum of the update exceeds 0.2 bohr, the whole update is scaled down (unlike the 0.5 bohr of geometry optimization).
- The end images are not moved. `tools/check_geometry_md.py` checks, on a 2D model potential, that this method converges, the image spacing becomes uniform, and the highest image approaches the saddle-point energy.

### 8.4 Parallelization, output, restart [Public] [Convention]

[Public] Manual sections 51.3, 51.4, 51.7:

- `MD.NEB.Parallel.Number` is the number of MPI process groups to which images are assigned. If not specified, parallelization is first over images, and when there are more processes than images, each image is also parallelized internally.
  Setting it to 1 disables image parallelization, and all processes compute images one at a time (minimum memory). [Convention] If the number of groups exceeds the number of processes or the number of images, it reverts to automatic (0).
- `System.Name.neb.opt`: optimization history (iteration, SD prefactor, maximum force, maximum displacement, \\( \lVert\mathbf G\rVert \\), sum of image energies).
  [Convention] The header of the maximum-displacement column says Å, but the values are in bohr.
- `System.Name.neb.ene`: per image, "index (0 and \\( P+1 \\) are the end points), total energy, distance to the neighboring image, cumulative distance from the end".
  [Convention] Columns 2 and 3 are distances (bohr) according to the header, but the values are the **squared distance** (bohr²) and its cumulative sum. In addition, a 5th column gives the difference in x coordinate between atom 2 and atom 1 (Å).
- `System.Name.neb.xyz`: coordinates of each image (XYZ format; the header line has the image index and energy). [Convention] `System.Name.neb.md` is also written, in the same format as `.md` (section 10), with the force columns \\( -\mathbf G\_p \\) and velocity columns 0.
- `System.Name.dat#`: restart input containing the coordinates of all images (`NEBp.Atoms.SpeciesAndCoordinates`) and `scf.restart on`.
  [Convention] The energies of the end points are saved in `System.Name_0_rst/`, and the end-point SCFs are skipped on restart.

## 9. Relation to the SCF [Convention]

- Carried over between steps: the density of the previous step (restart), the history for density extrapolation ([Density matrix and charge mixing](density-mixing.md) sections 6 and 7),
  and the grid origin ([Real-space grid](grid-poisson-xc.md) section 1.4). The mixing weight is reset to `scf.Max.Mixing.Weight` at the start of each step (section 4 of the same chapter).
- At each step, the set of atom pairs ([sparse block format](../design/sparse-block-format.md)) and
  matrix elements such as two-center integrals are rebuilt to follow the atomic motion. In lattice optimization and `DF`, the number of real-space grid divisions is fixed to that of the first step,
  and only the grid spacing changes with the lattice vectors (the effective cutoff changes). The treatment in NPT is [Open].
- Even at steps where the SCF did not converge, the update continues with the forces at that point. `scf_conv=` in the `.md` header line becomes 0.
- Forces are analytic derivatives at the final density matrix of the SCF ([Total energy and forces](energy-force.md)). The accuracy of conservation in MD depends on the SCF convergence criterion, so in MD `scf.criterion` should be made sufficiently small.

## 10. Output files [Public] [Convention]

### `.md`

[Public] Manual section 18.6: for each step, a block consisting of a "number of atoms" line, a header line, and per-atom lines (14 columns).
After a free-form string, the header line carries `Cell_Vectors=` (the 9 x, y, z values of \\( \mathbf a\_1,\mathbf a\_2,\mathbf a\_3 \\), Å) and `scf_conv=` (1 if the SCF converged).
Each atom line is "element symbol, x, y, z (Å), \\( F\_x, F\_y, F\_z \\) (Hartree/bohr), \\( v\_x, v\_y, v\_z \\) (m/s), net charge (electron charge taken as negative), magnetic moment (\\( \mu\_B \\)), spin-direction angles \\( \theta, \phi \\) (degrees)".

[Convention]

- Free-form part of the header line: for MD (`NVE`, `NVT_VS`, `NVT_NH`, `NVT_VS2`, `NVT_Langevin`, `NPT_*`),
  `time= ` (\\( (n-1)\Delta t \\), fs), `Energy= ` (`Utot`), `Temperature= `, `(Given Temp.= ` \\( T\_0 \\)`)`.
  For `NVT_VS4`, per-group `Temperature[g]= ` instead of the temperature. Otherwise (geometry optimization, etc.), only `time=` and `Energy=`.
- `.md` is written right after the SCF, **before** the coordinate update, so the coordinates and forces correspond to the same structure. The velocity columns are the velocity held at that point
  (the half-step \\( \mathbf v\_{n-1/2} \\) for `NVE`; \\( \mathbf u \\) after the update of the previous step for the `NVT_VS` family).
  For geometry optimization and `NPT_*`, the velocity columns are 0. The temperature is \\( T\_B \\) of section 2 (except for `NPT_*` and `NVT_VS4`).
- The force columns are the values **before** the MD zero-sum correction (the correction is applied after writing `.md`).
- The net charge and magnetic moment are computed from the per-atom electron numbers (up and down), but at which stage these values are taken (initial configuration, or charge analysis after SCF) is [Open]. They are excluded from comparison [Design].
- `System.Name.md2` contains only the structure of the last step (overwritten).

### `.ene`

[Convention] When `MD.Type` is `NOMD`, geometry optimization, lattice optimization, `EvsLC`, `DF`, `NVE`, or `NVT_*`, one line is appended per step (34 columns, space-separated).
The columns required for compatibility with a reimplementation are as follows. The other columns are internal quantities of the existing implementation and are not compared.

| Column | Content |
|---|---|
| 1 | step \\( n \\) |
| 2 | time \\( (n-1)\Delta t \\) (fs; written using `MD.TimeStep` even for geometry optimization) |
| 14 | `Ukc` (Hartree) |
| 15 | `Utot` (Hartree) |
| 16 | `Utot` + `Ukc` |
| 17 | chemical potential (Hartree). Column 32 is the same |
| 18 | target temperature \\( T\_0 \\) (K) |
| 19 | temperature \\( T\_A \\) (K) |
| 20, 21, 22 | \\( \zeta \\), \\( \zeta' \\), Nosé–Hoover conserved quantity (`NVT_NH`) |
| 23–25, 26 | x, y, z of \\( \mathbf a\_1 \\) and z of \\( \mathbf a\_3 \\) (Å) |
| 30 | maximum force (definition of the method of that step, section 4.1) |
| 31 | twice the total spin magnetic moment |

[Public] Manual section 20.1 explains columns 1, 2, and 14–17.

### Others

- `System.Name.out`: contains, for each step of geometry optimization, "iteration, prefactor, maximum force, maximum displacement (Å), `Utot`" (for lattice optimization, `Enpy` and volume (Å³) are added; manual section 19.4)
  as the table "History of geometry optimization" / "History of cell optimization",
  and the coordinates and forces of the last structure as the `<coordinates.forces` block. [Convention] The prefactor column in quasi-Newton rows remains the steepest-descent prefactor (meaningless).
- `System.Name.npt.ene`: section 6.7. `System.Name.neb.*`: section 8.4. `System.Name.EvsLC`, `System.Name.DF`: section 5.5.
- `System.Name.dat#`: sections 4.11, 6.8, 8.4.

## 11. Edge cases [Convention]

- With `MD.maxIter` = 1 (default), even if geometry optimization or MD is specified, the run ends after a single SCF (coordinates are not updated).
- If geometry optimization is already converged at the first step, the run ends there.
- If the GDIIS coefficient matrix is singular (e.g. the same gradient twice in a row), the run ends as "converged" (section 4.1).
- Violation of the curvature condition: `RF` replaces the denominator by \\( 10^{-10} \\) and continues updating; `BFGS` and NEB update unconditionally (diverging if the denominator is near 0); `EF` skips the update.
- If the `Schlegel` Hessian has an atom pair with \\( r \le B \\), negative curvature or divergence occurs.
- If MD with a thermostat starts from random initial velocities and there are fixed components, the velocities of fixed components are not 0 and keep entering the `NVT_VS` temperature.
- In `NVT_VS2`, for an atom with all components fixed, \\( T\_I = 0 \\) and the scaling factor is undefined (there are no components to multiply, so the result is not affected).
- In `NVT_NH`, if the first step of `MD.TempControl` is greater than 1, the target temperature used in the prediction at \\( n = 1 \\) is unset (same as 0).
- In NEB, if neighboring images have equal energies, the tangent may be undefined (section 8.2).
- In NPT steps beyond the range of `MD.CellPressureControl`, the external pressure is undefined (section 6.7).
- An atom exactly at fractional coordinate 1 is not wrapped (section 7).

## Requirements for implementations [Design]

- **Structure**: separate the "calculator that returns forces and stress" from the "integrator/optimizer that advances the state". Make the integrator/optimizer a pure function that receives only
  \\( (\mathbf x, A, \mathbf g, \Sigma, U) \\) and its own state (velocities, thermostat, history, Hessian),
  with no dependence on the SCF, grid, or file I/O. This enables unit tests with analytic model potentials.
- **State saving**: serialize the optimizer/integrator state as a single structure and write it to the restart file (there is no need to match the format of the existing implementation's `_rst` directory).
  Test that a restarted calculation agrees step by step with an uninterrupted calculation.
- **Default method and compatibility mode**: on the following points, allow switching between "recommended", which behaves differently from the existing implementation, and "compatible", which matches it.

| Item | Recommended (default) | Compatible |
|---|---|---|
| Integration with velocity scaling | velocity Verlet + scaling | form of section 6.4 |
| Removal of velocity average | zero total momentum (mass-weighted) | simple average |
| Curvature violation in RF's BFGS | skip the update | denominator set to \\( 10^{-10} \\) |
| Singular GDIIS | discard history and continue | end as converged |
| NPT NH conserved quantity | quantity corresponding to the actual equations of motion | uses \\( \tfrac12 Q\zeta^2 \\) |
| WV \\( \dot F \\) | \\( \dot{\boldsymbol\sigma}\_i\cdot\boldsymbol\sigma\_j + \boldsymbol\sigma\_i\cdot\dot{\boldsymbol\sigma}\_j \\) | \\( 2\dot{\boldsymbol\sigma}\_i\cdot\boldsymbol\sigma\_j \\) |
| Wrapping in isolated systems | no | yes |
| Distances in `.neb.ene` | distance | squared distance |

- **Temperature definition**: make \\( T = 2K/(k\_B N\_{\mathrm{dof}}) \\) a single function, with \\( N\_{\mathrm{dof}} \\) selectable from "\\( 3N \\)", "number of free components", and "number of free components − 3 (center-of-mass constraint)". Which one each output uses is managed by a table (section 2).
- **Random numbers**: for Langevin and initial velocities, use an environment-independent generator whose seed can be specified in the input.
- **Units**: define the input/output conversion factors of section 2 in one place.
- **Constraints**: distinguish fixed flags by type into two kinds: Cartesian components (atomic-coordinate optimization, MD) and fractional components (`RFC*`).
- **NEB**: make the SCF of each image callable independently (parallelization is the responsibility of the outer layer), and save and reuse the SCF results of the end points. For periodic systems, provide as an option a preprocessing step that matches end-point coordinates by the minimum image.

## Tests [Design]

`tools/check_geometry_md.py` checks those items in the "Analytic" and "Invariant" rows below that can be written with models.

| Kind | Content | Typical tolerance |
|---|---|---|
| Analytic | unit conversions (fs, u, m/s, GPa) and values from independent constants | relative \\( 10^{-5} \\) (\\( 10^{-4} \\) for GPa) |
| Analytic | velocity Verlet for a harmonic oscillator: energy error bounded, and 1/4 when \\( \Delta t \\) is halved | ratio \\( 4 \pm 0.2 \\) |
| Analytic | characteristic roots \\( \lvert z\rvert = \sqrt{1-(\omega\Delta t)^2/2} \\) of the integration of section 6.4, and the energy decay rate in simulation | relative 5% on the decay rate |
| Analytic | Nosé–Hoover (discretization of section 6.5): fluctuation of the extended-system energy is \\( O(\Delta t^2) \\) | ratio \\( 4 \pm 1 \\), fluctuation at most \\( 0.05\\,k\_BT \\) (\\( \Delta t = 0.5 \\) fs) |
| Analytic | NPT_NH-type thermostat: the conserved quantity contains \\( \tfrac32 Q\zeta^2 \\) (with a high-accuracy ODE solution) | \\( 10^{-6}\\,k\_BT \\) |
| Analytic | temperature after velocity scaling is \\( T\_0 + (T-T\_0)\alpha \\) | \\( 10^{-9} \\) K |
| Analytic | GDIIS coefficients: bordered equation vs. closed form | \\( 10^{-10} \\) |
| Analytic | BFGS secant conditions \\( H\_{\mathrm{new}}\Delta\mathbf x = \Delta\mathbf g \\), \\( H^{-1}\_{\mathrm{new}}\Delta\mathbf g = \Delta\mathbf x \\), and the two being mutual inverses | \\( 10^{-9} \\) |
| Analytic | RF: the step from the eigenvector of the augmented Hessian is \\( -(H-\lambda)^{-1}\mathbf g \\), \\( \lambda = \mathbf g\cdot\Delta\mathbf x < \lambda\_{\min}(H) \\) | \\( 10^{-10} \\) |
| Analytic | on a quadratic function (9 dimensions), `Opt`, `DIIS`, `BFGS`, `RF`, `EF` all converge to `MD.Opt.criterion` | converges |
| Analytic | \\( S = \Omega A^{-\mathsf T} \\), \\( \partial\Omega/\partial A = S \\) | \\( 10^{-7} \\) |
| Analytic | NEB: on a 2D model, converges, the relative standard deviation of image spacings is small, and the highest point is close to the saddle-point energy | 2%, 0.05 |
| Invariant | initial velocities have zero total momentum and temperature \\( T\_1 \\) | \\( 10^{-12} \\) |
| Invariant | the model Hessian is symmetric, and its pair part is 0 on translations | \\( 10^{-12} \\) |
| Invariant | wrapping: \\( s = 1.25 \to 0.25 \\), \\( -0.25 \to 0.75 \\), \\( 1.0 \\) kept | exact |
| Invariant | in geometry optimization, coordinates of fixed components are bitwise unchanged | exact |
| Invariant | restart: coordinates and velocities of a calculation stopped at \\( n \\) steps and restarted from `.dat#` vs. a continuous calculation (same method, same random seed) | \\( 10^{-10} \\) bohr (when the SCF initial density is the same) |
| Invariant | fluctuation of `Utot` + `Ukc` in an NVE real DFT calculation (H₂O molecule, \\( \Delta t = 0.5 \\) fs, `scf.criterion` \\( 10^{-10} \\)) | about \\( 10^{-4} \\) Hartree over 100 steps (trend check) |
| Comparison | existing-implementation fixture: given the same initial structure and the same sequence of forces (recorded forces), the one-step update of `Opt`, `DIIS`, `BFGS`, `RF`, `EF` | \\( 10^{-8} \\) bohr |
| Comparison | likewise, the one-step coordinates, velocities, and \\( \zeta \\) of `NVE`, `NVT_VS` (compatibility mode), `NVT_NH` | \\( 10^{-10} \\) (relative) |
| Comparison | one-step lattice and coordinates of `OptC1`, `OptC5`, `RFC5` (given stress and forces) | \\( 10^{-8} \\) bohr |
| Comparison | `Schlegel` initial Hessian (H₂O, Si crystal) | \\( 10^{-10} \\) |
| Regression | `RF` optimization of the H₂O molecule: converged structure (bond lengths, angle) and number of iterations | bond length \\( 10^{-3} \\) Å, angle 0.1°, iteration count indicative |
| Regression | `RFC5` for Si crystal: converged lattice constant | \\( 10^{-3} \\) Å |
| Regression | NEB: barrier height for the manual's example (rotation of C₂H₄, etc.) | 10 meV |

Since the existing implementation has no path for "supplying forces externally", the comparison fixtures are built by running a small system including SCF for 2–3 steps,
recording the coordinates, forces, and stress of each step (`.md`, `.out`), and reproducing the one-step update from them.

## [Open]

- The formula of `MD.Opt.Init.Hessian FF` (initial Hessian by force-field fitting). It is not described in the public documents and no usage example has been found, so the proposal is to exclude it from the reimplementation.
- The exact formula for the lattice-lattice and atom-lattice parts of the `Schlegel` initial Hessian for `RFC5`–`RFC7` (affects only the convergence path).
- The effective default of `MD.Opt.DIIS.History` in NEB. The NEB input reader uses 3 as the default, but it appears to be overwritten with 2 by the input processing of each image's SCF; this needs to be confirmed by running it.
- The definition of the WV `Ubaro` (column 8 of `.npt.ene`) and the purpose of `NPT.WV.F0` (not used in the equations of motion).
- The meaning of columns 3–13, 27–29, and 33–34 of `.ene` (internal quantities; whether to make them compatibility targets).
- Behavior when `MD.Fixed.XYZ` is given with NPT (the manual says they cannot be combined, but there is no input check).
- Differences from upstream OpenMX (CPU version): what was consulted is a derivative of the GPU version, and it is unconfirmed whether the form of the `NVT_VS`-family integration (section 6.4) and the defaults are the same upstream.
- The unit of `Langevin.Friction.Factor` (judged to be the inverse of the atomic time unit, but without support from the public documents).
