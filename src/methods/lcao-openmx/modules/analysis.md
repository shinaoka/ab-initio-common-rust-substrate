# Analysis features: bands, density of states, charge analysis, exchange interactions, Berry phase, optical response

> Status: first draft (awaiting review). The equations and methods of the manual 4.0 (Sections 22, 23, 25, 33–35, 44–45, 50, 53, 56, 58, 59, 62, 70–71),
> the technical notes (TN9 "Charge Analysis", TN11 "Electric Polarization by Berry Phase", "Calculation of Optical Matrix Elements",
> "Dielectric Function and Conductivity", "Analysis of spin texture in the k-space"), and the papers
> (Han–Ozaki–Yu, PRB **70**, 184421 (2004); Terasawa et al., JPSJ **88**, 114706 (2019); Sawahata et al., JJAP **57**, 030309 (2018);
> Kotaka et al., JJAP **52**, 035204 (2013)) have been rewritten in the notation of this note.
> File layouts, defaults, and numerical choices were confirmed as behavior of the existing implementation and are marked [Convention].
> The Mulliken sum rule, the tetrahedron method, the Berry phase, the Chern number, and the relations among the exchange-coupling formulas are checked numerically with `tools/check_analysis.py`.

## Scope of this chapter

This chapter covers the analyses (post-processing) performed after the SCF has finished. In order of priority:

1. The format of the `.scfout` file, which is the input to post-processing, and the OpenMX API (Section 1).
2. Charge analysis: Mulliken, Voronoi, and ESP charges (Section 4).
3. Density of states (DOS) and partial density of states (PDOS), and the post-processing program DosMain (Section 3).
4. Exchange coupling constants \\( J\_{ij} \\) (jx, Section 6).
5. Electric polarization by the Berry phase (polB, Section 7).

In addition, the band structure (Section 2), orbital-decomposed total energy (Section 5), Chern number, Berry curvature, and \\( \mathbb Z\_2 \\) invariant (Section 8),
spin texture in k-space (Section 9), optical conductivity and dielectric function (Section 10), and output of quantities on the grid (cube files, molecular orbitals, STM, work function, Section 11) are treated briefly.

The definitions of the Mulliken spin and orbital moments in noncollinear calculations are in Section 6 of [Noncollinear DFT energy and forces](noncollinear-energy-force.md) and are not repeated here.
NEGF outputs (transmission, eigenchannels, etc.) are in [NEGF](negf.md). Natural orbital analysis (NPA, NBO, LNAO), CWF/MLWF, band unfolding,
and COHP are out of scope (see the [Open] section).

## Notation [Design]

The notation of [k-points and Bloch sums](kpoints-bloch.md) and [Density matrix and charge mixing](density-mixing.md) is followed.

| Symbol | Meaning |
|---|---|
| \\( O\_{I\alpha,J\beta}(\mathbf L) \\) | Real-space matrix ([sparse block format](../design/sparse-block-format.md)). \\( \mathbf L \\) is the periodic image of the column-side atom \\( J \\) |
| \\( O(\mathbf k) = \sum\_{\mathbf L} e^{i\mathbf k\cdot\mathbf L} O(\mathbf L) \\) | Bloch sum (phase on lattice vectors only) |
| \\( c\_{a\mu}(\mathbf k) \\) | Eigenvectors of the generalized eigenvalue problem (\\( c^\dagger S(\mathbf k) c = 1 \\)). \\( a = (I\alpha) \\) is the global orbital index |
| \\( \rho^{\sigma}(\mathbf L) \\) | Real-space density matrix (for one spin in the non-spin-polarized case) |
| \\( \mu \\), \\( T \\) | Chemical potential (Hartree), electronic temperature (K) |
| \\( f(\varepsilon) \\) | Fermi distribution \\( 1/(1+e^{(\varepsilon-\mu)/k\_{\mathrm B}T}) \\) |
| \\( g\_s \\) | Spin degeneracy (2 for non-spin-polarized, 1 otherwise) |
| \\( \Omega \\) | Unit-cell volume (bohr³) |
| \\( \mathbf a\_i \\), \\( \mathbf b\_i \\) | Lattice vectors and reciprocal-lattice vectors (\\( \mathbf a\_i\cdot\mathbf b\_j = 2\pi\delta\_{ij} \\)) |
| \\( D^{d}\_{I\alpha,J\beta}(\mathbf L) \\) | Position-operator matrix \\( \langle\phi\_{I\alpha}\rvert (r\_d - R\_{I,d}) \lvert\phi\_{J\beta}(\mathbf L)\rangle \\) (\\( d = x, y, z \\), measured from the center of atom \\( I \\)) |
| \\( \nabla^{d}\_{I\alpha,J\beta}(\mathbf L) \\) | Derivative-operator matrix \\( \langle\phi\_{I\alpha}\rvert \partial\_d \lvert\phi\_{J\beta}(\mathbf L)\rangle \\) |

\\( \lvert\phi\_{J\beta}(\mathbf L)\rangle \\) is shorthand for \\( \phi\_{J\beta}(\mathbf r - \mathbf R\_J - \mathbf L) \\).

## 1. The `.scfout` file and the OpenMX API

### 1.1 Role [Public]

Manual Section 70.1: setting the input keyword `HS.fileout` (`on`/`off`, default `off`) to `on` writes
`System.Name.scfout` at the end of the calculation. It contains the Hamiltonian, overlap matrix, density matrix, and matrix elements of the position and momentum operators,
and the post-processing programs (jx, polB, kSpin, Z2FH, calB, and the reading example) take only this file as input.
The manual does not give the file layout and only refers to the source of the reading example. The layout below consists of facts confirmed from the existing implementation.

### 1.2 Binary layout [Convention]

- Unformatted binary; integers are 4 bytes, reals are 8 bytes (double precision), in the byte order of the writing machine. There are no record markers.
- A "block matrix" below consists, for atoms \\( I = 1, \dots, N\_{\mathrm{atom}} \\) in order and for each element \\( h = 0, \dots, n\_I \\) of their neighbor list,
  of the \\( n\_{\mathrm{orb}}(I)\times n\_{\mathrm{orb}}(J\_h) \\) blocks laid out in **row-major** order (rows = orbitals of atom \\( I \\), columns = orbitals of the partner \\( J\_h \\); the column index varies fastest).
  The total number of elements of one block matrix is \\( \sum\_I \sum\_{h=0}^{n\_I} n\_{\mathrm{orb}}(I)\\, n\_{\mathrm{orb}}(J\_h) \\).

| Order | Contents | Type × count |
|---|---|---|
| 1 | Header: \\( N\_{\mathrm{atom}} \\), \\( s + 4v \\), three NEGF atom counts, \\( N\_{\mathrm{cell}} \\) | int × 6 |
| 2 | Maximum order of the position operator \\( p\_{\max} \\) (1 in OpenMX 4.0) | int × 1 |
| 3 | Cartesian translation vectors of periodic images (bohr), \\( (\cdot, x, y, z) \\) for each \\( r \\), \\( r = 0, \dots, N\_{\mathrm{cell}} \\) | double × \\( 4(N\_{\mathrm{cell}}+1) \\) |
| 4 | Integer coordinates of periodic images, \\( (\cdot, n\_1, n\_2, n\_3) \\) for each \\( r \\) | int × \\( 4(N\_{\mathrm{cell}}+1) \\) |
| 5 | Number of orbitals of each atom \\( n\_{\mathrm{orb}}(I) \\) | int × \\( N\_{\mathrm{atom}} \\) |
| 6 | Number of neighbor atoms of each atom \\( n\_I \\) (excluding itself) | int × \\( N\_{\mathrm{atom}} \\) |
| 7 | Indices of neighbor atoms (1-based): \\( n\_I + 1 \\) per atom | int × \\( \sum\_I (n\_I+1) \\) |
| 8 | Periodic-image indices \\( r \\) of neighbor atoms: \\( n\_I + 1 \\) per atom | int × \\( \sum\_I (n\_I+1) \\) |
| 9 | Lattice vectors \\( \mathbf a\_1, \mathbf a\_2, \mathbf a\_3 \\) (bohr), each as \\( (\cdot, x, y, z) \\) | double × 12 |
| 10 | Reciprocal-lattice vectors \\( \mathbf b\_1, \mathbf b\_2, \mathbf b\_3 \\) (including \\( 2\pi \\), bohr⁻¹), each as \\( (\cdot, x, y, z) \\) | double × 12 |
| 11 | Atomic positions (Cartesian, bohr), \\( (\cdot, x, y, z) \\) per atom | double × \\( 4N\_{\mathrm{atom}} \\) |
| 12 | Hamiltonian \\( H \\): one block matrix per spin component \\( \sigma = 0, \dots, s \\) (Hartree) | double |
| 13 | (only when \\( s = 3 \\)) Additional imaginary parts of the Hamiltonian: block matrices for 3 components | double |
| 14 | Overlap matrix \\( S \\): one block matrix | double |
| 15 | Position operator \\( D^{d} \\): block matrices in the order \\( d = x, y, z \\), and for each \\( d \\) in the order of orders \\( 1, \dots, p\_{\max} \\) | double |
| 16 | Derivative operator \\( \nabla^{d} \\): block matrices in the order \\( d = x, y, z \\) | double |
| 17 | Density matrix \\( \rho \\): one block matrix per spin component \\( \sigma = 0, \dots, s \\) | double |
| 18 | Additional imaginary parts of the density matrix: **always** block matrices for 2 components (written even in collinear calculations) | double |
| 19 | Solver index (table below) | int × 1 |
| 20 | \\( \mu \\) (Hartree), \\( T \\) (K), 3 components of the core-charge dipole (Debye), 3 components of the uniform-background-charge dipole (Debye), number of valence electrons \\( N\_{\mathrm e} \\), total spin \\( S \\) | double × 10 |
| 21 | Number of lines of the input file \\( n\_{\mathrm{line}} \\) | int × 1 |
| 22 | Each line of the input file: \\( n\_{\mathrm{line}} \\) fixed-length records of 256 bytes | char × \\( 256\\,n\_{\mathrm{line}} \\) |

Meaning of each item:

- **Second header element**: the spin type \\( s \\) (0: non-spin-polarized, 1: collinear spin-polarized, 3: noncollinear) and the format version \\( v \\) combined as \\( s + 4v \\).
  OpenMX 3.9 and 4.0 use \\( v = 3 \\), so the value is one of 12, 13, 15. A reader decomposes it as \\( s = (\cdot) \bmod 4 \\), \\( v = \lfloor(\cdot)/4\rfloor \\), and
  refuses to read if \\( v \ne 3 \\) (\\( v = 0 \\) is 3.8 and earlier; 1 and 2 are development-version formats).
- **Byte-order detection**: when the second element exceeds \\( 4\times 3 + 3 = 15 \\), the existing reader regards the whole file as byte-swapped and reads each value with its bytes reversed.
  If the value is still out of range after reversal, it is an error.
- **NEGF atom counts**: the numbers of atoms in the central, left, and right regions; 0 for non-NEGF calculations. The reader does not use them.
- **Periodic images**: \\( N\_{\mathrm{cell}} = (2N\_c+1)^3 - 1 \\), where \\( N\_c \\) is the maximum integer coordinate in each direction sufficient to include all neighbor atoms.
  Index \\( r = 0 \\) is the origin cell, and \\( r \ge 1 \\) enumerate \\( [-N\_c, N\_c]^3 \\) minus the origin in lexicographic order with \\( n\_1 \\) outermost and \\( n\_3 \\) innermost.
  A reader should not depend on the index order but use the integer coordinates of item 4.
- **4-element vectors**: items 3, 4, 9, 10, 11 do not use the first element (its contents are undefined); elements 1–3 are \\( x, y, z \\) (or \\( n\_1, n\_2, n\_3 \\)).
- **Neighbor list**: element 0 of atom \\( I \\) is the atom itself (atom index \\( I \\), periodic image \\( r = 0 \\)). Elements 1 onward are all \\( (J, \mathbf L) \\) satisfying
  \\( \lVert \mathbf R\_J + \mathbf L - \mathbf R\_I\rVert \le r\_{c,I} + r\_{c,J} \\) (equality included), in unspecified order.
  If \\( (J, \mathbf L) \\) is a neighbor of \\( I \\), then \\( (I, -\mathbf L) \\) is a neighbor of \\( J \\). This is the same set of atom pairs as in [sparse block format](../design/sparse-block-format.md).
- **\\( H \\) of item 12**: in collinear calculations, real matrices for \\( \sigma = 0 \\) (up) and \\( \sigma = 1 \\) (down, only for \\( s = 1 \\)).
  In noncollinear calculations (\\( s = 3 \\)), the 4 components \\( (H\_0, H\_1, H\_2, H\_3) \\), together with the 3 components \\( (X^{\uparrow\uparrow}, X^{\downarrow\downarrow}, X^{\uparrow\downarrow}) \\) of item 13
  (imaginary parts due to spin-orbit coupling, DFT+U, spin constraints, Zeeman terms, etc.), give

  $$
  H^{\uparrow\uparrow} = H\_0 + iX^{\uparrow\uparrow},\qquad H^{\downarrow\downarrow} = H\_1 + iX^{\downarrow\downarrow},\qquad
  H^{\uparrow\downarrow} = H\_2 + i\left(H\_3 + X^{\uparrow\downarrow}\right),\qquad H^{\downarrow\uparrow}\_{J\beta,I\alpha}(-\mathbf L) = \overline{H^{\uparrow\downarrow}\_{I\alpha,J\beta}(\mathbf L)}
  $$

  The assembly of the \\( 2N\times 2N \\) matrix in k-space is as in [k-points and Bloch sums](kpoints-bloch.md).
- **Timing of items 12 and 17**: they are written at the end of the calculation (including geometry optimization and MD). \\( H \\) is the last diagonalized Hamiltonian (built from the input density of the last SCF iteration),
  and \\( \rho \\) is the output density matrix obtained by diagonalizing it. If the SCF has converged, the two are consistent to within the convergence criterion.
- **Orbital contraction**: when orbital optimization (contraction) is used, \\( n\_{\mathrm{orb}} \\), \\( H \\), and \\( S \\) refer to the contracted basis.
- **\\( D \\) of item 15**: \\( \langle\phi\_{I\alpha}\rvert (r\_d - R\_{I,d})^{p} \lvert\phi\_{J\beta}(\mathbf L)\rangle \\) (\\( p = 1, \dots, p\_{\max} \\)), with the origin at the center of the **row-side atom** \\( I \\).
  It is evaluated not analytically but as a sum over the same real-space grid as the Poisson equation (multiplied by the grid-point volume), so its accuracy depends on the grid fineness (`scf.energycutoff`).
- **\\( \nabla \\) of item 16**: \\( \langle\phi\_{I\alpha}\rvert \partial\_d \lvert\phi\_{J\beta}(\mathbf L)\rangle \\) (the momentum-operator matrix is \\( -i \\) times this). This is also a grid sum.
- **\\( \rho \\) of items 17 and 18**: in collinear calculations, \\( \rho^\sigma(\mathbf L) \\) of [Density matrix and charge mixing](density-mixing.md) (for one spin in the non-spin-polarized case), and item 18 is 0.
  In noncollinear calculations, item 17 is \\( (\mathrm{Re}\\,\rho^{\uparrow\uparrow}, \mathrm{Re}\\,\rho^{\downarrow\downarrow}, \mathrm{Re}\\,\rho^{\uparrow\downarrow}, \mathrm{Im}\\,\rho^{\uparrow\downarrow}) \\),
  and item 18 is \\( (\mathrm{Im}\\,\rho^{\uparrow\uparrow}, \mathrm{Im}\\,\rho^{\downarrow\downarrow}) \\).
- **Solver index of item 19**: corresponds to the value of `scf.EigenvalueSolver`.

  | Value | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|
  | Solver | Recursion | Cluster | Band | NEGF | DC | GDC | Cluster-DIIS | Krylov | Cluster2 | EGAC | DC-LNO | Cluster-LNO |

  For `Band` with `scf.Kgrid` equal to `1 1 1`, the calculation goes through the Cluster path (Γ point with periodic boundaries) ([k-points and Bloch sums](kpoints-bloch.md)), but 3 (Band) is written to the file.
- **Item 20**: \\( T \\) is the value of `scf.ElectronicTemperature` (K). The dipoles are in Debye (1 \\( e \\)·bohr = 2.54174776 Debye);
  the core-charge term is \\( \sum\_I Z\_I^{\mathrm{val}} \mathbf R\_I \\), and the background-charge term is
  \\( -\int\_\Omega (Q/\Omega)\\,\mathbf r\\,d^3r \\) for the uniform background \\( Q/\Omega \\) of charge doping (`scf.system.charge` \\( = Q \\)), evaluated as a grid sum. The total spin \\( S \\) is, in the collinear case, Mulliken's \\( \frac12\sum\_I (M^\uparrow\_I - M^\downarrow\_I) \\) (Section 4.1).
- **Item 22**: each line is a string obtained by reading the input file line by line up to 255 characters (including the newline character; contents after the terminating NUL are undefined).
  Lines longer than 255 characters are split into multiple records. A reader extracts up to the NUL and concatenates to reconstruct the input file.
  If no input file name is given, items 21 and 22 are absent.

### 1.3 Correspondence with the sparse block format [Design]

The contents of `.scfout` map directly onto the [sparse block format](../design/sparse-block-format.md) of this note.

- Atom pair \\( (I, J\_h, \mathbf n\_h) \\): \\( J\_h \\) from item 7, and \\( \mathbf n\_h \\) from \\( r \\) of item 8 converted to integer coordinates via item 4.
- Block values: each block of items 12–18 (row-major). If the memory layout of the sparse block format is column-major, transpose on reading.
- Noncollinear: combine into complex spin blocks using the formula above.

A new implementation provides a **reader and writer** of `.scfout` as a converter to and from the sparse block format (all analyses below take the sparse block format as input).
For compatibility with existing post-processing programs, the writer reproduces the version-3 layout exactly.

### 1.4 OpenMX API [Public] [Convention]

[Public] Manual Section 71: a mechanism for calling OpenMX as a library from external programs (C, Fortran, Python);
static libraries (C, Fortran) and a Python extension module are provided. The typical use is to split an MPI communicator and
compute multiple input files in parallel.

[Convention] There is a single public entry point, taking three arguments (MPI communicator, input file name, number of OpenMP threads) and returning an integer (0 for normal termination).
The Fortran version receives the communicator as a Fortran handle, and the Python version is a single function taking the same three arguments.

- When called, it performs the calculation described in the input file (SCF, geometry optimization, MD, NEB, etc.) on the given communicator exactly as in a normal run,
  and all output is written to the usual files (`.out`, `.scfout`, cube files, etc.).
- There is no function that returns in-memory matrices or physical quantities to the caller. Quantities used for post-processing are received via files such as `.scfout`.

[Design] The API of a new implementation has functions that return \\( H \\), \\( S \\), \\( \rho \\) in sparse block format and the structure without going through files.
An entry point identical to the existing API ("pass an input file and run") is provided as a thin compatibility layer.

## 2. Band structure

### 2.1 Input [Public] [Convention]

[Public] Manual Section 22: specified by `Band.dispersion` (`on`/`off`, default `off`), the number of paths `Band.Nkpath`, and
the path block `<Band.kpath` (each line: number of points \\( N \\), start point \\( (q\_1, q\_2, q\_3) \\), end point \\( (q\_1, q\_2, q\_3) \\), labels of the start and end points).
`<Band.KPath.UnitCell` can give a different unit cell (e.g., the conventional cell of an fcc or bcc lattice) in which the k-point coordinates are interpreted. If omitted, `Atoms.UnitVectors` is used.

[Convention]

- Coordinates are fractional with respect to the reciprocal-lattice vectors. When `<Band.KPath.UnitCell` is given, coordinates written in terms of the reciprocal-lattice vectors of that cell are
  converted to fractional coordinates of the cell used in the calculation (by converting \\( \mathbf k \\) to Cartesian and then \\( q\_i = \mathbf k\cdot\mathbf a\_i/2\pi \\)).
- The points of each path are \\( \mathbf q\_j = \mathbf q\_{\mathrm s} + (\mathbf q\_{\mathrm e} - \mathbf q\_{\mathrm s})\\,(j-1)/(N-1) \\), \\( j = 1, \dots, N \\),
  **including both ends**, equally spaced in fractional coordinates. A shared endpoint of adjacent paths appears twice. \\( N = 1 \\) is not allowed (division by zero).
- At each point the [generalized eigenvalue problem](generalized-eigen.md) is solved and all eigenvalues are output (post-processing rule A′).

### 2.2 Output `System.Name.Band` [Convention]

Text, written in the following order.

1. \\( n\\;\\; s\\;\\; \mu \\): \\( n \\) is the total number of orbitals (written as \\( 2N \\) in the noncollinear case, with \\( s \\) written as 0), \\( \mu \\) in Hartree.
2. The 9 components of the reciprocal-lattice vectors (Cartesian components of \\( \mathbf b\_1, \mathbf b\_2, \mathbf b\_3 \\), including \\( 2\pi \\), bohr⁻¹).
3. The number of paths, followed for each path by \\( N \\), start point, end point, and two labels.
4. In the order path, point, spin: a line \\( n\\;\\; q\_1\\;\\; q\_2\\;\\; q\_3 \\) and a line with \\( n \\) eigenvalues (Hartree, \\( \mu \\) not subtracted).

[Public] The conversion program bandgnu13 produces gnuplot data with the distance along the path as the horizontal axis and \\( \varepsilon - \mu \\) (eV) as the vertical axis (manual Section 22).

## 3. Density of states and partial density of states

### 3.1 Calculation stages [Public] [Convention]

[Public] Manual Section 23.1: with `Dos.fileout on`, eigenvalue and eigenvector information for the DOS is written after the SCF to `System.Name.Dos.val` (text) and
`System.Name.Dos.vec` (binary), and the DOS and PDOS are obtained with a separate program, DosMain.
`Dos.Erange` is the energy range (eV, with \\( \mu \\) as origin), and `Dos.Kgrid` is the number of k-point divisions. The DOS can also be output with O(N) methods (Gaussian broadening only).
Section 23.2: `DosGauss.fileout on` computes the DOS on the fly with Gaussian broadening without writing eigenvectors; it cannot be used together with `Dos.fileout`.

[Convention]

| Keyword | Default | Remarks |
|---|---|---|
| `Dos.fileout` | `off` | |
| `Dos.Erange` | `-20.0 20.0` (eV) | Relative to \\( \mu \\) |
| `Dos.Kgrid` | value of `scf.Kgrid` | |
| `DosGauss.fileout` | `off` | Usable only with Band (including periodic Γ-point calculations) |
| `DosGauss.Num.Mesh` | 200 | Number of energy divisions |
| `DosGauss.Width` | 0.2 (eV) | Width \\( a \\) of the Gaussian |

DOS calculation in the Band case:

- The k-points form a Γ-centered grid \\( q\_{i,l} = l/N\_i \\) (\\( l = 0, \dots, N\_i - 1 \\)), of this form regardless of `scf.Generation.Kpoint`
  (in directions with 2 or more divisions, the same small shift \\( \delta \\) as in [k-points and Bloch sums](kpoints-bloch.md) is added). Time-reversal reduction is **not** performed.
  The whole grid is needed because the tetrahedron method refers to periodically neighboring points.
- The eigenvalue problem is solved at each k-point and spin (regularization by rule A′), and the range of state indices \\( [l\_{\min}, l\_{\max}] \\) falling in the energy window is determined in common for all k-points.
  At each k-point the range runs from the state just before the first state exceeding \\( \mu + E\_{\min} \\) to the first state exceeding \\( \mu + E\_{\max} \\); after taking the union over all k-points,
  a margin of \\( \max(\lfloor n/20\rfloor, 10) \\) states is added on each side (clipped to \\( 1 \le l \le n \\)).
- The weight of state \\( \mu \\) (k-point \\( \mathbf k \\), spin \\( \sigma \\)) on each orbital is the Mulliken gross-population weight

  $$
  w\_{a\mu}(\mathbf k) = \mathrm{Re}\left[\overline{c\_{a\mu}(\mathbf k)}\\,\bigl(S(\mathbf k)\\,c\_{\mu}(\mathbf k)\bigr)\_a\right],\qquad \sum\_a w\_{a\mu}(\mathbf k) = 1
  $$

  (the same partitioning as the Mulliken charges of Section 4.1).

### 3.2 Format of `.Dos.val` and `.Dos.vec` [Convention]

`.Dos.val` is a keyword-format text file with the following items.

| Keyword | Contents |
|---|---|
| `mode` | 1: Cluster/Band (per eigenstate), 5: DC, 6: NEGF, 7: `DosGauss.fileout` (per energy grid point) |
| `NonCol` | 0: collinear, 1: noncollinear |
| `N` | Total number of orbitals \\( n \\) |
| `Nspin` | \\( s \\) (written as 1 in the noncollinear case) |
| `Erange` | `Dos.Erange` (Hartree) |
| `irange` | \\( l\_{\min}\\; l\_{\max} \\) (1-based) |
| `Kgrid` | \\( N\_1\\; N\_2\\; N\_3 \\) (`1 1 1` for Cluster) |
| `atomnum` | Number of atoms |
| `<WhatSpecies` … `WhatSpecies>` | Species index of each atom (0-based) |
| `SpeciesNum` | Number of species |
| `<Spe_Total_CNO` … `Spe_Total_CNO>` | Number of orbitals of each species |
| `MaxL` | Maximum angular momentum |
| `<Spe_Num_CBasis` … `Spe_Num_CBasis>` | Number of radial functions (multiplicity) for each species and each \\( l \\) |
| `ChemP` | \\( \mu \\) (Hartree) |
| `<SpinAngle` … `SpinAngle>` | (noncollinear only) Direction \\( (\theta, \varphi) \\) of the Mulliken spin of each atom (radians) |
| `<Eigenvalues` … `Eigenvalues>` | One line per k-point and spin: \\( l\_1\\; l\_2\\; l\_3 \\) (integer indices of the k-point) followed by \\( \varepsilon\_{l\_{\min}}, \dots, \varepsilon\_{l\_{\max}} \\) (Hartree) |

`.Dos.vec` is a sequence of records in the order k-point, spin, state \\( l\_{\min} \dots l\_{\max} \\); each record is int × 3 (integer indices of the k-point) followed by
4-byte real × \\( n \\) (weights \\( w\_{a\mu} \\), in orbital order). The order of k-points is determined by the work distribution, so a reader identifies k-points by their integer indices.
In the noncollinear case each record has 4-byte real × \\( 2n \\): for each atom, \\( n\_{\mathrm{orb}}(I) \\) "majority-spin" weights and
\\( n\_{\mathrm{orb}}(I) \\) "minority-spin" weights, with the direction of that atom's Mulliken spin as the quantization axis.

### 3.3 Computation in DosMain [Public] [Convention]

[Public] Manual Section 23.1: DosMain reads `.Dos.val` and `.Dos.vec` and interactively asks for the method (tetrahedron method or Gaussian broadening),
total DOS or PDOS, the target atoms of the PDOS, and the Gaussian width \\( a \\) (eV, function \\( \exp(-(E/a)^2) \\)).
In the output, column 1 is the energy (eV, with \\( \mu \\) as origin), column 2 the DOS (eV⁻¹), and column 3 the integrated value. With spin polarization, the up and down DOS and integrated values are listed.

[Convention]

- **Energy grid**: \\( E\_j = E\_{\min} + (E\_{\max} - E\_{\min})\\,j/(N\_E - 1) \\), \\( j = 0, \dots, N\_E - 1 \\) (relative to \\( \mu \\)).
  For the tetrahedron method \\( N\_E = \lfloor (E\_{\max} - E\_{\min})/1\\,\mathrm{meV}\rfloor \\); for the Gaussian method \\( N\_E = \lfloor 5(E\_{\max} - E\_{\min})/a\rfloor \\).
- **Gaussian**:

  $$
  D\_\sigma(E) = \frac{1}{N\_k}\sum\_{\mathbf k}\sum\_{\nu} W\_{\nu\sigma}(\mathbf k)\\, \frac{1}{a\sqrt\pi}\exp\left[-\left(\frac{E + \mu - \varepsilon\_{\nu\sigma}(\mathbf k)}{a}\right)^2\right]
  $$

  where \\( W = \sum\_a w\_a \\) (= 1) for the total DOS, and the weights of the target orbitals for the PDOS. The Gaussian is truncated at about \\( \pm 3a \\) from its center (\\( 3a/\Delta E + 3 \\) grid points).
  The unit is eV⁻¹ per spin per unit cell (\\( \nu \\) is the state index, \\( \mu \\) the chemical potential).
- **Tetrahedron method**: each small cube of the k-point grid (8 vertices, numbered \\( 4i + 2j + k \\) with \\( i, j, k \in \\{0, 1\\} \\)) is divided into 6 tetrahedra sharing the diagonal 0–7,
  \\( \\{0,1,2,5\\}, \\{1,2,3,5\\}, \\{2,3,5,7\\}, \\{0,2,4,5\\}, \\{2,4,5,6\\}, \\{2,5,6,7\\} \\) (wrapped periodically).
  The vertex energies of each tetrahedron are sorted as \\( e\_1 \le e\_2 \le e\_3 \le e\_4 \\), and the linearly interpolated DOS (Blöchl et al., PRB **49**, 16223 (1994), Appendix C, normalized to 1 per tetrahedron)

  $$
  g(E) = \begin{cases}
  \dfrac{3(E - e\_1)^2}{e\_{21}e\_{31}e\_{41}} & e\_1 < E < e\_2 \\\\[2mm]
  \dfrac{1}{e\_{31}e\_{41}}\left[3e\_{21} + 6(E - e\_2) - \dfrac{3(e\_{31} + e\_{42})(E - e\_2)^2}{e\_{32}e\_{42}}\right] & e\_2 < E < e\_3 \\\\[2mm]
  \dfrac{3(e\_4 - E)^2}{e\_{41}e\_{42}e\_{43}} & e\_3 < E < e\_4
  \end{cases}
  $$

  (\\( e\_{ij} = e\_i - e\_j \\), 0 outside the range) is summed and divided by \\( 6N\_k \\). Blöchl's correction term is not used.
  For the total DOS the weight of each state is taken as constant, equal to its value at the reference vertex of the small cube; for the PDOS, the formula with vertex weights linearly interpolated within the tetrahedron (the derivative of Appendix B of the same paper) is used.
- **Integrated value**: the DOS integrated from \\( E\_0 \\) to \\( E\_j \\) by Simpson's rule.
- **Output**: in the non-spin-polarized case, \\( 2D \\) and its integral are output. In the spin-polarized case the order is \\( D\_\uparrow \\), \\( -D\_\downarrow \\) (**with a negative sign**), \\( \int D\_\uparrow \\), \\( \int D\_\downarrow \\).
  File names are `System.Name.DOS.Tetrahedron`, `System.Name.DOS.Gaussian`, and for the PDOS `System.Name.PDOS.Tetrahedron.atom<I>.<orbital name>` etc.
  The PDOS is split into files per atom, per orbital resolved by \\( l \\), multiplicity, and \\( m \\) (`s1`, `p1` (\\( p\_x \\)), etc.), and the sum for that atom.

### 3.4 `DosGauss.fileout` [Public] [Convention]

[Public] Manual Section 23.2: the DOS is computed directly with Gaussian broadening without writing the eigenvector file. The tetrahedron method cannot be used.

[Convention] `.Dos.val` has `mode 7`, `Kgrid` is `1 1 1`, `<Eigenvalues` lists the energy grid points, and `.Dos.vec` contains, for each grid point,
the per-orbital weights accumulated with Gaussian broadening. DosMain reads them with the same procedure.

## 4. Charge analysis

### 4.1 Mulliken charges [Public] [Convention]

[Public] TN9 Eq. 1–6 (notation rewritten): the Mulliken gross population of orbital \\( a = (I\alpha) \\) is

$$
M^{\sigma}\_{I\alpha} = \sum\_{\mathbf L}\sum\_{J\beta} \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\, S\_{I\alpha,J\beta}(\mathbf L)
$$

and satisfies the sum rule \\( g\_s\sum\_\sigma\sum\_{I\alpha} M^\sigma\_{I\alpha} = N\_{\mathrm e} \\) (the same as the sum rule of [Density matrix and charge mixing](density-mixing.md)).
The atomic population is \\( M^\sigma\_I = \sum\_{\alpha\in I} M^\sigma\_{I\alpha} \\). Written in k-space, it equals the sum of the weights \\( w\_{a\mu} \\) of Section 3.1 weighted by the occupation numbers.

[Public] Manual Section 35.1: the results are written by default to `System.Name.out`, with `Up spin`, `Down spin`, `Sum`, `Diff` per atom and
a per-orbital decomposition (including sums over \\( m \\), and over \\( m \\) and multiplicity).

[Convention]

- No input keyword is needed; it is also computed in every SCF iteration (it is also used for atomic spins in Section 6 of [Noncollinear](noncollinear-energy-force.md), DFT+U initialization, etc.).
- The "charge" in the output is the electron population (a positive number). The net charge of an atom is \\( Z\_I^{\mathrm{val}} - (M^\uparrow\_I + M^\downarrow\_I) \\), which is written in the atom lines of cube files (Section 11.1).
- In the non-spin-polarized case, the same \\( M\_I \\) (for one spin) is printed for `Up spin` and `Down spin`.
- In the noncollinear case, \\( N^\uparrow, N^\downarrow, \theta, \varphi \\) are constructed from a 2×2 Mulliken-type matrix (Section 6.1 of [Noncollinear DFT energy and forces](noncollinear-energy-force.md)).
- The total spin (used for item 20 of `.scfout` etc.) is, in the collinear case, \\( S = \frac12\sum\_I (M^\uparrow\_I - M^\downarrow\_I) \\).
- When orbital contraction is used, \\( \rho \\) and \\( S \\) in the contracted basis are used.

### 4.2 Voronoi charges (fuzzy cells) [Public] [Convention]

[Public] TN9 Eq. 7–16, manual Section 35.2: space is partitioned smoothly into atoms by Becke's fuzzy cells (Becke, JCP **88**, 2547 (1988); Becke–Dickson, JCP **89**, 2993 (1988)),

$$
N^\sigma\_I = \int d^3r\\, w\_I(\mathbf r)\\, n^\sigma(\mathbf r),\qquad
w\_I(\mathbf r) = \frac{P\_I(\mathbf r)}{\sum\_J P\_J(\mathbf r)},\qquad
P\_I(\mathbf r) = \prod\_{J\ne I} s\left(\mu\_{IJ}(\mathbf r)\right)
$$

$$
\mu\_{IJ}(\mathbf r) = \frac{\lvert\mathbf r - \mathbf R\_I\rvert - \lvert\mathbf r - \mathbf R\_J\rvert}{\lvert\mathbf R\_I - \mathbf R\_J\rvert},\qquad
s(\mu) = \tfrac12\left[1 - p^{(k)}(\mu)\right],\qquad p(x) = \tfrac32 x - \tfrac12 x^3
$$

(\\( p^{(k)} \\) is \\( p \\) composed \\( k \\) times). No correction for atomic size is applied, so the partition is determined by the structure alone.
The keyword is `Voronoi.charge` (`on`/`off`, default `off`); `Up spin`, `Down spin`, `Sum`, `Diff` and the fuzzy-cell volume (Å³) are output per atom.

[Convention]

- **The number of compositions is 6**. TN9 states \\( k = 3 \\) (Becke's recommended value), but the existing implementation composes \\( p \\) 6 times in total.
  The more compositions, the sharper the boundaries, approaching Voronoi cells. \\( \sum\_I w\_I = 1 \\) holds for either number (`tools/check_analysis.py`).
- Ranges of \\( P\_I \\) and the denominator sum: for grid points \\( \mathbf r \\) belonging to atom \\( I \\), the denominator sum runs over atoms \\( J \\) in the neighbor list of atom \\( I \\) (Section 1.2, including itself),
  and the product in each \\( P\_J \\) runs over atoms in the neighbor list of \\( J \\) (excluding itself). These are not sums/products over all atoms.
  Since the neighbor list is determined by the PAO cutoff radii, strictly the partition depends not on "the structure alone" but also on the cutoffs. If the denominator is below \\( 10^{-14} \\), \\( w\_I = 0 \\).
- The integral is a sum over the real-space grid points contained in the cutoff spheres of the basis functions of atom \\( I \\) (multiplied by the grid-point volume). \\( n^\sigma \\) is the valence electron density (excluding core charge).
  In the noncollinear case, the total density and the magnetization vector are integrated with \\( w\_I \\) and converted to \\( N^\uparrow, N^\downarrow \\) and a direction.
- The volume is \\( \int w\_I\\,d^3r \\) (converted to Å³).

### 4.3 Fitting to the electrostatic potential (ESP charges) [Public] [Convention]

[Public] TN9 Eq. 17–21, manual Section 35.3: the electrostatic potential of point charges \\( Q\_I \\) placed on each atom,
\\( V^{\mathrm{ESP}}(\mathbf r) = \sum\_I Q\_I/\lvert\mathbf r - \mathbf R\_I\rvert \\), is fitted by least squares to the DFT potential.
The sample points are real-space grid points between two shells formed by the unions of spheres centered on the atoms with radii \\( \lambda\_1 \\) and \\( \lambda\_2 \\) times the van der Waals radii.
A total-charge constraint (or total-charge and dipole constraints) is imposed with Lagrange multipliers (Singh–Kollman; Besler–Merz–Kollman).
Run the SCF with `level.of.fileout` at 1 or higher to produce `System.Name.out` and `System.Name.vhart.cube`, and then run the program esp as
`esp System.Name -c <0|1> -s λ1 λ2`. It is intended for small molecules and is not suited to large molecules or crystals.

[Convention]

- **Fitted quantity**: the **sign-reversed** value of `.vhart.cube`, \\( \delta V\_{\mathrm H}(\mathbf r) \\) (the Hartree potential of the difference \\( \delta n \\) between the electron density and the superposition of atomic densities, Section 11.1),
  i.e., \\( -\delta V\_{\mathrm H} \\). This is the electrostatic potential of the charge \\( -\delta n \\); since the electrostatic potential of the superposition of neutral atoms (core charges and atomic densities) vanishes outside the atoms,
  \\( Q\_I \\) is the net charge of the atom (negative if electrons are gained). TN9 says "the DFT Hartree potential", but precisely it is this difference quantity.
- **Equations**: with \\( b\_I(\mathbf r) = \sum\_{\mathbf n} 1/\lvert\mathbf r - \mathbf R\_I - \mathbf L\_{\mathbf n}\rvert \\) (the sum over periodic images is over a finite range determined by the grid of the cube file),
  \\( A\_{IJ} = \sum\_p b\_I(\mathbf r\_p)\\,b\_J(\mathbf r\_p) \\), and \\( B\_I = -\sum\_p \delta V\_{\mathrm H}(\mathbf r\_p)\\,b\_I(\mathbf r\_p) \\),
  - `-c 0`: constrain \\( \sum\_I Q\_I = 0 \\) (the total charge \\( Q\_{\mathrm{tot}} \\) is **always 0**; `scf.system.charge` is not used).
  - `-c 1`: additionally constrain \\( \sum\_I Q\_I \mathbf R\_I = \mathbf d \\) (\\( \mathbf d \\) is the total dipole written in `.out`, converted to atomic units).

  The augmented linear system is solved. The default is `-c 0`.
- **Default scale factors**: if `-s` is omitted, \\( (\lambda\_1, \lambda\_2) = (1.4, 5.0) \\) (the manual's example uses `1.4 2.0`). It is an error unless \\( \lambda\_2 > \lambda\_1 \ge 0 \\).
- **van der Waals radii**: Bondi's values (Å; e.g., H 1.20, C 1.70, N 1.55, O 1.52, F 1.47, Si 2.10, P 1.80, S 1.80, Cl 1.75), looked up by atomic number.
  If an element without a registered value is present (e.g., most transition metals), the program terminates with an error.
- Atomic positions and the lattice are read from the cube-file header. The dipole (Debye) is \\( \sum\_I Q\_I\mathbf R\_I \\).
- The "RMS" in the output is a diagnostic value whose definition differs from the standard root mean square (the square root of the sum of squared differences over sample points, atoms, and periodic images, divided by the number of sample points and the grid-point volume); do not use it for comparison.
- Another scheme (`-l`: imposing only the total-charge constraint on local charges given in a file) also remains, but it is not described in the manual and is out of scope.

## 5. Orbital-decomposed total energy [Public] [Convention]

[Public] Manual Section 25: with `Energy.Decomposition on` (default `off`), the total energy decomposed into atoms and orbitals is output.

$$
U\_{\mathrm{tot}} = \sum\_\sigma\left(U\_{\mathrm{kin}}^\sigma + U\_{v}^\sigma + U\_{\mathrm{con}}^\sigma\right) + U\_{\mathrm{core}} + U\_{\mathrm{H0}} + U\_{\mathrm{vdW}}
$$

where \\( U\_{\mathrm{con}} \\) is a term from a formally added constant potential. The feature is described as experimental.

[Convention] For each orbital \\( a = (I\alpha) \\) and spin \\( \sigma \\),

$$
U\_{\mathrm{kin},a}^\sigma = \sum\_{\mathbf L, b}\rho^\sigma\_{ab}(\mathbf L)\\,T\_{ab}(\mathbf L),\qquad
U\_{v,a}^\sigma = \sum\_{\mathbf L, b}\rho^\sigma\_{ab}(\mathbf L)\left[H^\sigma\_{ab}(\mathbf L) - T\_{ab}(\mathbf L)\right],\qquad
U\_{\mathrm{con},a}^\sigma = c\sum\_{\mathbf L, b}\rho^\sigma\_{ab}(\mathbf L)\\,S\_{ab}(\mathbf L)
$$

(\\( T \\) is the kinetic-energy matrix; both spins are the same in the non-spin-polarized case). The constant \\( c \\) is the difference between the total energy minus the core–core, \\( E\_{\mathrm{H0}} \\) (neutral-atom), and
van der Waals terms, and the band energy, divided by the number of electrons; it is chosen so that the sum reproduces the total energy
(for the definitions of each term see [Total energy and forces](energy-force.md)). The core–core, \\( E\_{\mathrm{H0}} \\), and vdW terms are split only per atom and placed on the first orbital of that atom.
In the noncollinear case, products of the imaginary parts and spin off-diagonal components of \\( \rho \\) and \\( H \\) are added to \\( U\_v \\).
The output is `System.Name.DecE`, with a per-atom table (\\( U\_{\mathrm{tot}} \\), per spin, \\( U\_{\mathrm{kin}} \\), \\( U\_v \\), \\( U\_{\mathrm{con}} \\), \\( U\_{\mathrm{core}} + U\_{\mathrm{H0}} \\), \\( U\_{\mathrm{vdW}} \\), Hartree) and a per-orbital table.
It cannot be used together with orbital contraction.

## 6. Exchange coupling constants (jx)

### 6.1 Formula [Public]

This is the formula based on the magnetic force theorem of Liechtenstein et al. (JMMM **67**, 65 (1987)), extended by Han–Ozaki–Yu (PRB **70**, 184421 (2004)) to nonorthogonal localized bases.
With the exchange-splitting matrix of atom \\( i \\),

$$
\Delta\_i = H^{\uparrow}\_{ii}(\mathbf 0) - H^{\downarrow}\_{ii}(\mathbf 0)\qquad(n\_{\mathrm{orb}}(i)\times n\_{\mathrm{orb}}(i)\text{ one-center block})
$$

and using the parts of the eigenvectors \\( c^\sigma\_{\mu}(\mathbf k) \\) on atoms \\( i \\) and \\( j \\),

$$
J\_{ij} = \frac14\\,\frac{1}{N\_k}\sum\_{\mathbf k}\sum\_{\mu\nu}
\frac{f(\varepsilon^\uparrow\_\mu) - f(\varepsilon^\downarrow\_\nu)}{\varepsilon^\downarrow\_\nu - \varepsilon^\uparrow\_\mu}\\,
\left[c^{\uparrow\dagger}\_{\mu}\Delta\_j c^{\downarrow}\_{\nu}\right]\left[c^{\downarrow\dagger}\_{\nu}\Delta\_i c^{\uparrow}\_{\mu}\right]
$$

(corresponding to manual Section 45.1 Eq. 42, 47; the argument \\( \mathbf k \\) is omitted). This is the sum over periodic images \\( \sum\_{\mathbf L} J\_{i\mathbf 0, j\mathbf L} \\).
The individual \\( J\_{i\mathbf 0,j\mathbf L} \\) are obtained from the Green's function

$$
G^\sigma\_{ij}(z;\mathbf L) = \frac{1}{N\_k}\sum\_{\mathbf k} e^{-i\mathbf k\cdot\mathbf L}\sum\_\mu \frac{c^\sigma\_{i\mu}(\mathbf k)\\,c^{\sigma\dagger}\_{j\mu}(\mathbf k)}{z - \varepsilon^\sigma\_\mu(\mathbf k)}
\qquad\left(= \left[(zS - H^\sigma)^{-1}\right]\_{i\mathbf 0, j\mathbf L}\right)
$$

and the poles of the continued-fraction expansion of the Fermi function (Ozaki, PRB **75**, 035123 (2007); Terasawa et al., JPSJ **88**, 114706 (2019)) as

$$
J\_{i\mathbf 0, j\mathbf L} = \frac{k\_{\mathrm B}T}{2}\sum\_{p=1}^{N\_P} R\_p\\, \mathrm{Re}\\,\mathrm{Tr}\left[\Delta\_i\\, G^\downarrow\_{ij}(z\_p;\mathbf L)\\, \Delta\_j\\, G^\uparrow\_{ji}(z\_p;-\mathbf L)\right],\qquad
z\_p = \mu + i\\,k\_{\mathrm B}T\\,\zeta\_p
$$

(manual Eq. 44–46). The poles \\( \zeta\_p \\) and residues \\( R\_p \\) are obtained from the generalized eigenvalue problem of order \\( M = 2N\_P \\), \\( A\mathbf x = \lambda B\mathbf x \\)
(\\( A \\) tridiagonal with zero diagonal and off-diagonal elements \\( -\tfrac12 \\), \\( B = \mathrm{diag}(1, 3, 5, \dots, 2M-1) \\), \\( \mathbf x^{\mathsf T}B\mathbf x = 1 \\)),
from the negative eigenvalues \\( \lambda\_p \\) (the \\( N\_P \\) smallest) and the first components \\( x\_{p,1} \\) of the eigenvectors, as

$$
\zeta\_p = -1/\lambda\_p,\qquad R\_p = -\tfrac14\\, x\_{p,1}^2\\, \zeta\_p^2
$$

which gives \\( 1/(1+e^{x}) = \tfrac12 + \sum\_p R\_p\left[(x - i\zeta\_p)^{-1} + (x + i\zeta\_p)^{-1}\right] \\)
(`tools/check_analysis.py`: error \\( 10^{-15} \\) for \\( \lvert x\rvert \le 30 \\) with \\( N\_P = 60 \\); the low poles coincide with the Matsubara poles \\( \zeta = (2p-1)\pi \\), \\( R\_p = -1 \\)).
\\( \sum\_{\mathbf L} J\_{i\mathbf 0,j\mathbf L} \\) agrees with the periodic-sum formula above (to \\( 10^{-17} \\) with the same script).

### 6.2 Sign and normalization conventions [Design] (checked)

Compared with the magnetic force theorem (the second-order change of the grand potential \\( \Omega \\) at fixed \\( \mu \\) with respect to rotations of the spin directions),
\\( J \\) of the \\( \frac14 \\) formula above (Band mode) corresponds to

$$
{}-\frac{\partial^2\Omega}{\partial\theta\_i\\,\partial\theta\_j}\bigg\rvert\_{\theta=0} = 2J\_{ij}
\quad\Longleftrightarrow\quad
E = -\sum\_{i\ne j} J\_{ij}\\,\mathbf e\_i\cdot\mathbf e\_j = -2\sum\_{i<j}J\_{ij}\\,\mathbf e\_i\cdot\mathbf e\_j
$$

(\\( \mathbf e\_i \\) is the unit vector of the spin direction; confirmed to an accuracy of \\( 10^{-8} \\) on the 3-site nonorthogonal cluster of `tools/check_analysis.py`).
\\( J > 0 \\) is ferromagnetic. The spin magnitude is absorbed into \\( J \\).

### 6.3 Program and configuration file [Public] [Convention]

[Public] Manual Sections 45.2–45.6: run the SCF with `HS.fileout on`, then run `jx System.Name.scfout jx.config`.
Only collinear calculations (Cluster or Band) are supported; noncollinear calculations are out of scope.

[Convention] The configuration file is in keyword format, and the following items are read.

| Keyword | Default | Contents |
|---|---|---|
| `Flag.PeriodicSum` | `off` | For Band: `on` gives the periodic sum \\( \sum\_{\mathbf L}J\_{i\mathbf 0,j\mathbf L} \\), `off` the individual \\( J\_{i\mathbf 0,j\mathbf L} \\) |
| `Num.Poles` | 60 | Number of poles \\( N\_P \\) (used only with `PeriodicSum off`) |
| `Num.Kgrid` | `1 1 1` | Number of k-point divisions (0 or less is an error) |
| `Num.ij.pairs` | 0 | Number of pairs to compute |
| `Bunch.ij.pairs` | `Num.ij.pairs` | Number of pairs computed together at once (does not affect results) |
| `<ijpairs.cellid` … | — | Each line `i j n1 n2 n3` (atom indices 1-based, \\( \mathbf L = \sum n\_l\mathbf a\_l \\)) |
| `<ijpairs.nocellid` … | — | Each line `i j` |
| `Flag.LNO` | `off` | Project \\( \Delta \\) onto localized natural orbitals (LNO) (Band only) |
| `LNO.occ.cutoff` | 0.0 | Occupation threshold for LNO |

- For Band with `PeriodicSum off`, only `<ijpairs.cellid` is read. Otherwise `<ijpairs.nocellid` takes precedence; if absent, `<ijpairs.cellid` is used with its lattice vectors ignored.
- **Solver dispatch**: if the solver index in `.scfout` is Cluster (2) or Cluster-DIIS (7), Cluster mode is used (Γ point, single MPI process only); if Band (3), Band mode;
  otherwise nothing is computed.
- **k-points**: \\( q\_{i,l} = -\tfrac12 + (2l+1)/(2N\_i) \\) (Monkhorst–Pack type; Γ is not included if \\( N\_i \\) is even), all with weight \\( 1/N\_k \\), no reduction.
- **Eigenvalue problem**: orthogonalization by eigendecomposition of \\( S(\mathbf k) \\), without regularization of small eigenvalues of \\( S \\) (different from all rules of [generalized eigenvalue problem](generalized-eigen.md)).
- **Temperature and chemical potential**: \\( \mu \\) and \\( T \\) from `.scfout`, and \\( k\_{\mathrm B} = 3.166813628\times10^{-6} \\) Hartree/K.
- **Numerics of individual \\( J \\)**: \\( \eta = 10^{-10} \\) Hartree is added to the imaginary part of \\( z\_p \\) (negligible effect).
- **Cluster mode**: the same form of sum is taken at the Γ point only, but **the value is \\( \tfrac12 \\) times that of the Band-mode formula**
  (because the exchange splitting is taken as \\( \tfrac12\Delta \\) and the prefactor as \\( \tfrac12 \\), giving \\( \tfrac18 \\) overall). Solving the same Γ-point problem in Band mode (`Num.Kgrid 1 1 1`, periodic sum) gives twice the value.
  Which normalization is intended cannot be determined from the public documents ([Open]). A new implementation uses the Band-mode normalization (Section 6.2) by default and offers the Cluster normalization as a compatibility option.
- **LNO**: for each atom and spin, the one-center block of \\( (\rho S)\_{ii} \\) (summed over partner atoms) is diagonalized, and with the projector \\( P^\sigma\_i \\) onto eigenvectors whose occupation exceeds the threshold,
  \\( \Delta\_i \to P^\uparrow\_i\Delta\_iP^\downarrow\_i \\) (on the \\( j \\) side, \\( P^\downarrow\_j\Delta\_jP^\uparrow\_j \\)). The eigenvectors are right eigenvectors of a nonsymmetric matrix.
- **Output**: standard output only; for each pair, \\( J \\) is printed in meV and mRy (\\( 1\\,\mathrm{Ha} = 2000\\,\mathrm{mRy} \\)). For Band with `PeriodicSum off`, the integer coordinates of the lattice vector are also printed.
- The periodic sum for `i = j` includes the \\( \mathbf L = \mathbf 0 \\) (self) term. When used as a physical coupling, exclude it from the sum of the individual \\( J \\).

## 7. Electric polarization by the Berry phase (polB)

### 7.1 Formula [Public]

The Berry-phase formula of King-Smith–Vanderbilt (PRB **47**, 1651 (1993)) and Resta (RMP **66**, 899 (1994)) is written in LCAO following TN11 Eq. 1–11.
For the direction of the reciprocal-lattice vector \\( \mathbf b\_1 \\), with a string of k-points along \\( \mathbf b\_1 \\), \\( \mathbf k\_j = \mathbf k\_\perp + (j/N\_1)\mathbf b\_1 \\), \\( j = 0, \dots, N\_1 \\), and
\\( \mathbf b = \mathbf b\_1/N\_1 \\), the overlap matrix between occupied states

$$
M\_{\mu\nu}(\mathbf k, \mathbf k + \mathbf b) = \langle\psi\_{\mu\mathbf k}\rvert e^{-i\mathbf b\cdot\mathbf r}\lvert\psi\_{\nu,\mathbf k+\mathbf b}\rangle\quad(\text{per unit cell})
$$

is formed, and from the argument of the product of its determinants the Berry phase of each string,

$$
\varphi(\mathbf k\_\perp) = -\arg\prod\_{j=0}^{N\_1-1}\det M(\mathbf k\_j, \mathbf k\_{j+1})
$$

is obtained; from its string average \\( \bar\varphi\_1 \\), the electronic contribution to the dipole (per unit cell, \\( e > 0 \\), electron charge \\( -e \\))

$$
\mathbf d\_{\mathrm{el}} = -\frac{g\_s\\,e}{2\pi}\sum\_{i=1}^{3}\bar\varphi\_i\\,\mathbf a\_i,\qquad \mathbf P = \frac{\mathbf d\_{\mathrm{el}} + \mathbf d\_{\mathrm{ion}}}{\Omega}
$$

is obtained. In the collinear spin-polarized case, the \\( \varphi \\) of each spin are added and \\( g\_s = 1 \\).

In LCAO, since the Bloch sum has phases on lattice vectors only, one can choose \\( c(\mathbf k + \mathbf G) = c(\mathbf k) \\), and the coefficients at \\( \mathbf k\_0 \\) can be used as is for the end of the string \\( \mathbf k\_{N\_1} = \mathbf k\_0 + \mathbf b\_1 \\)
(periodic gauge). The matrix elements are

$$
M\_{\mu\nu}(\mathbf k, \mathbf k + \mathbf b) = \sum\_{ab}\overline{c\_{a\mu}(\mathbf k)}\\,T\_{ab}\\,c\_{b\nu}(\mathbf k + \mathbf b),\qquad
T\_{I\alpha,J\beta} = \sum\_{\mathbf L} e^{i(\mathbf k+\mathbf b)\cdot\mathbf L}\\,\langle\phi\_{I\alpha}\rvert e^{-i\mathbf b\cdot\mathbf r}\lvert\phi\_{J\beta}(\mathbf L)\rangle
$$

and TN11 Eq. 10–11 expands \\( e^{-i\mathbf b\cdot\mathbf r} \\) to first order about the center of the row-side atom:

$$
\langle\phi\_{I\alpha}\rvert e^{-i\mathbf b\cdot\mathbf r}\lvert\phi\_{J\beta}(\mathbf L)\rangle \approx X\_{I\alpha,J\beta}(\mathbf L) \equiv
e^{-i\mathbf b\cdot\mathbf R\_I}\left[S\_{I\alpha,J\beta}(\mathbf L) - i\\,\mathbf b\cdot\mathbf D\_{I\alpha,J\beta}(\mathbf L)\right]
$$

(\\( \mathbf D \\) is the position-operator matrix of item 15 of `.scfout`). The error is \\( O(\lvert\mathbf b\rvert^2) \\) and vanishes as the string division is refined.
`tools/check_analysis.py` confirms, on a Rice–Mele chain with point-like orbitals, that this construction (phases on lattice vectors only and periodic gauge) agrees with Resta's supercell formula.

### 7.2 Choices of the existing implementation [Convention]

- **Input**: run `polB System.Name.scfout` and enter, in order, the numbers of k-point divisions \\( N\_1\\; N\_2\\; N\_3 \\) (each at least 1) and the directions to compute (0 or 1 for each of the 3 reciprocal-lattice directions) (manual Section 44).
  There are no restrictions on the solver or the spin type (0, 1, 3).
- **Strings**: for direction \\( i \\), strings are formed for all combinations of fractional coordinates \\( l/N \\) (\\( l = 0, \dots, N - 1 \\), Γ-centered) in the other two directions, and the points on a string are also \\( j/N\_i \\).
  The number of strings is the product of the division numbers of the other two directions.
- **Symmetrization of \\( T \\)**: from \\( X \\) above, \\( T \\) is built in two ways and averaged (the sum is over all neighbor pairs \\( (I, J, \mathbf L) \\)).

  $$
  T = \tfrac12\left(T^{(1)} + T^{(2)}\right),\qquad
  T^{(1)}\_{I\alpha,J\beta} = \sum\_{\mathbf L} e^{i(\mathbf k+\mathbf b)\cdot\mathbf L}X\_{I\alpha,J\beta}(\mathbf L),\qquad
  T^{(2)}\_{J\beta,I\alpha} = \sum\_{\mathbf L} e^{-i\mathbf k\cdot\mathbf L}X\_{I\alpha,J\beta}(\mathbf L)
  $$

  \\( T^{(2)} \\) approximates \\( \langle\phi\_{J\beta}\rvert e^{-i\mathbf b\cdot\mathbf r}\lvert\phi\_{I\alpha}(-\mathbf L)\rangle \\) by an expansion about the center of \\( I \\) instead of \\( J \\)
  (since \\( S \\) and \\( \mathbf D \\) are real, it takes a form without complex conjugation); averaging the two cancels the bias from the choice of the center of the first-order expansion.
- **Number of occupied states**: \\( \lfloor N\_{\mathrm e}/2\rfloor \\) for non-spin-polarized, \\( \lfloor N\_{\mathrm e}\rfloor \\) for noncollinear, and for collinear spin-polarized, per spin, the number of states with \\( \varepsilon \le \mu \\) at that k-point.
  In the non-spin-polarized and noncollinear cases, if there is a k-point where this number differs from the number of states with \\( \varepsilon \le \mu \\), a "metal" warning is issued (the calculation continues).
- **Eigenvalue problem**: the generalized eigenvalue problem is solved as is (no regularization). In the noncollinear case, \\( 2N \\)-dimensional spinor matrices are used, and the same \\( T \\) is used for both spins.
- **String average**: for \\( \psi = \arg\prod\det M \in (-\pi, \pi] \\) of each string, a reference \\( \psi\_0 = \arg\sum\_{\text{strings}} e^{i\psi} \\) is computed,
  each \\( \psi \\) is shifted by an integer multiple of \\( 2\pi \\) into \\( (\psi\_0 - \pi, \psi\_0 + \pi] \\), and then averaged (\\( \bar\varphi = -\overline{\psi} \\)).
- **Output**: the electronic, core (ionic), and background-charge contributions and the total dipole (Debye) and polarization (μC/cm²) are printed to standard output as 3 Cartesian components.
  The core and background-charge contributions use the values of item 20 of `.scfout` as is (the core contribution is not folded using the \\( 2\pi \\) ambiguity).
  The conversions are 1 \\( e \\)·bohr = 2.54174776 Debye and \\( P[\mu\mathrm C/\mathrm{cm}^2] = 5721.52891433\\,d[e\\,\mathrm{bohr}]/\Omega[\mathrm{bohr}^3] \\).
  \\( e\mathbf a\_i/\Omega \\) is also printed as the polarization quantum, but in the non-spin-polarized case the actual ambiguity is \\( g\_s e\mathbf a\_i/\Omega \\) (twice that).
- Since the polarization has an ambiguity of integer multiples of the polarization quantum, only differences are meaningful (e.g., the Born effective charge \\( Z^\*\_{\alpha\beta} = (\Omega/e)\\,\Delta P\_\alpha/\Delta u\_\beta \\)) (manual Section 44).

### 7.3 Dipole during SCF [Convention]

The dipole of an isolated system, \\( \sum\_I Z\_I^{\mathrm{val}}\mathbf R\_I - \int n(\mathbf r)\\,\mathbf r\\,d^3r \\) (including the background-charge term), is computed after the SCF as a grid sum and printed in Debye to `.out`.
The existing implementation also has a polarization path within the SCF that uses an analytic expression for the Berry connection, but its correspondence to input keywords is not described in the manual, so it is marked [Open].

## 8. Chern number, Berry curvature, \\( \mathbb Z\_2 \\) invariant

### 8.1 Method of Fukui–Hatsugai–Suzuki [Public]

Fukui–Hatsugai–Suzuki (JPSJ **74**, 1674 (2005)) and manual Section 58: a two-dimensional k-plane is divided into a mesh, and from the link variables between adjacent k-points,
\\( U(\mathbf k, \mathbf k') = \det M(\mathbf k, \mathbf k')/\lvert\det M\rvert \\) (\\( M \\) is the same overlap matrix as in Section 7.1, over the selected bands),
the lattice field strength of each small cell (plaquette),

$$
F = \arg\left[U(\mathbf k\_1, \mathbf k\_2)\\,U(\mathbf k\_2, \mathbf k\_3)\\,U(\mathbf k\_3, \mathbf k\_4)\\,U(\mathbf k\_4, \mathbf k\_1)\right]\in(-\pi, \pi]
$$

is computed, giving the Chern number \\( C = \frac{1}{2\pi}\sum F \\). If the mesh is fine enough, \\( C \\) is an integer; the Berry curvature is \\( F \\) divided by the plaquette area.
The \\( \mathbb Z\_2 \\) invariant (Fukui–Hatsugai, JPSJ **76**, 053702 (2007); Sawahata et al., JJAP **57**, 030309 (2018); manual Section 59) is obtained,
on half of the Brillouin zone of a time-reversal-invariant plane, after fixing the gauge of the Kramers pairs at the time-reversal-invariant momenta (TRIM), as
\\( D = \frac{1}{2\pi}\left[\sum\_{\text{half}} F - \oint A\right] \bmod 2 \\).

### 8.2 Existing implementation (calB, Z2FH) [Convention]

- Both accept only noncollinear (\\( s = 3 \\)) `.scfout` files. The overlap matrices are built with the same first-order expansion and symmetrization as in Section 7.2.
- **calB**: interactively reads (a) whether to compute the Berry curvature of a single band (diagonal element of \\( M \\)) or the sum over bands up to a specified number (determinant), (b) the number of bands (number of valence electrons if 0 or less),
  (c) the normal direction of the plane (0/1 for each of the 3 reciprocal-lattice directions), and (d) the two in-plane division numbers. The plane passes through Γ and divides \\( [0,1)\times[0,1) \\) equally.
  The output `BerryC<number of bands>.dat` contains a header with the division numbers, the number of bands, and the Chern number, followed by the fractional coordinates of the center of each plaquette and \\( F \\) divided by the area (Å²).
  As the area, the **Cartesian component** along the normal direction (not the magnitude) of the cross product of the two in-plane reciprocal-lattice vectors is used. Also, the coordinate in column 2 is incorrect when the two division numbers differ.
  The Chern number itself is not affected by these.
- **Z2FH**: interactively reads (a) the division number \\( m \\) along the half direction, (b) whether to compute the \\( k = 0 \\) plane, and (c) the starting plane (6 planes: \\( 0 \\) and \\( \pi \\) in each direction).
  The two in-plane directions are divided as \\( [-\tfrac12, \tfrac12] \\) (\\( 2m \\) divisions) and \\( [0, \tfrac12] \\) (\\( m \\) divisions); at the 4 TRIM points
  one member of a Kramers pair is generated from the other by \\( \Theta(\psi\_\uparrow, \psi\_\downarrow) = (\overline{\psi\_\downarrow}, -\overline{\psi\_\uparrow}) \\), and the states on the boundary of the half are also generated with \\( \Theta \\).
  The sum of \\( n = \frac{1}{2\pi}\left[\sum\_{\text{edges}}\arg U - \arg\prod U\right] \\) over plaquettes is taken mod 2, and \\( (\nu\_0; \nu\_1\nu\_2\nu\_3) \\) is written to `Z2.dat` from the values of the 6 planes.
  For links along the second in-plane direction (\\( m \\) divisions, width \\( \tfrac12 \\)), the \\( \mathbf b \\) used in the expansion of \\( e^{-i\mathbf b\cdot\mathbf r} \\) is \\( 1/m \\) rather than the actual step \\( 1/(2m) \\).
  Since the Bloch phases use the correct k, the effect is only on the magnitude of the position-operator correction. The conditions under which the integer result changes are [Open]; a new implementation uses the correct \\( \mathbf b \\).

## 9. Spin texture in k-space (kSpin) [Public] [Convention]

[Public] Manual Section 56, the technical note "Analysis of spin texture in the k-space", Kotaka et al. (JJAP **52**, 035204 (2013)):
from the `.scfout` of a noncollinear Band calculation, the expectation values of the spin, decomposed per atom and orbital, are computed for each state at specified k-points.
There are modes for tracing constant-energy lines (FermiLoop), a mesh on a plane (GridCalc), bands along a path (BandDispersion), and specified k-points (MulPOnly),
and a separate program that aggregates the results over atoms (MulPCalc).

[Convention]

- Input is in keyword format (`Filename.scfout`, `Filename.outdata`, `Calc.Type`, `Energy.Range` (eV, relative to \\( \mu \\), default \\( -0.5\dots0.5 \\)),
  `Band.Nkpath` and `<Band.kpath`, `Search.kCentral`, `Calc.Type.3mesh`, `k-plane.1stStep`, `k-plane.2ndStep`, `kRange.3mesh`, `k-plane.height`, `kRange.height`,
  `MulP.Vec.Scale`, `Spin.Degenerate`, `Calc.Bandbyband`, etc.). The solver is restricted to Band and the spin to noncollinear.
- The spin density matrix of state \\( \mu \\) per orbital \\( a \\) is built in Mulliken gross-population form,

  $$
  n^{\sigma\sigma'}\_{a\mu} = \overline{c\_{\sigma a,\mu}}\left(S(\mathbf k)\\,c\_{\sigma'\mu}\right)\_a
  $$

  summed over atoms (or orbitals) to give \\( \mathbf m = (2\\,\mathrm{Re}\\,n^{\uparrow\downarrow}, 2\\,\mathrm{Im}\\,n^{\uparrow\downarrow}, n^{\uparrow\uparrow} - n^{\downarrow\downarrow}) \\)
  (in units of the expectation value of the Pauli matrices, not multiplied by \\( \hbar/2 \\)), and output after multiplication by `MulP.Vec.Scale`. The sum over all orbitals equals \\( \langle\psi\rvert\boldsymbol\sigma\lvert\psi\rangle \\).
- Energies are output as \\( \varepsilon - \mu \\) in eV, and k-points in Cartesian coordinates (bohr⁻¹; Å⁻¹ in MulPCalc).
- FermiLoop searches a coarse mesh for plaquettes where the sign changes, finds points with \\( \varepsilon(\mathbf k) = E \\) (\\( E \\) is the midpoint of `Energy.Range`) using a fine mesh and Brent's method (enabled by default, 5 iterations), and connects them into lines.

## 10. Optical conductivity and dielectric function (CDDF) [Public] [Convention]

[Public] Manual Section 62, the technical notes "Calculation of Optical Matrix Elements" and "Dielectric Function and Conductivity":
after the SCF, the optical conductivity tensor \\( \sigma\_{ab}(\omega) \\) and the dielectric tensor \\( \varepsilon\_{ab}(\omega) \\) are computed with a Kubo–Greenwood-type formula.
The momentum matrix elements are built either from the LCAO \\( \langle\phi\rvert\nabla\lvert\phi'\rangle \\) (approach 0) or from the \\( [H, \mathbf r] \\) form (approach 1, including the nonlocal-pseudopotential contribution, also supporting noncollinear).

[Convention] Keywords and defaults:

| Keyword | Default | Contents |
|---|---|---|
| `CDDF.start` | `off` | |
| `CDDF.FWHM` | 0.2 eV | Broadening \\( \eta \\) (**the FWHM value itself** is used as \\( i\eta \\) in the Lorentzian denominator) |
| `CDDF.maximum_energy` | 10.0 eV | Upper limit of output frequencies |
| `CDDF.minimum_energy` | 0.0 eV | Affects only the step \\( (E\_{\max} - E\_{\min})/N \\) (the start of the grid does not move) |
| `CDDF.additional_maximum_energy` | 0.0 eV | Extends the range of included bands |
| `CDDF.frequency.grid.total_number` | 10000 | \\( N \\) |
| `CDDF.approach` | 0 | 0: \\( \langle\phi\rvert\nabla\lvert\phi'\rangle \\) (collinear only), 1: the \\( [H,\mathbf r] \\) form (collinear and noncollinear) |
| `CDDF.shift_ChemP` | 0.0 eV | Shift \\( \mu \\) |
| `CDDF.Kgrid` | value of `scf.Kgrid` | Must be `1 1 1` if and only when `scf.Kgrid` is `1 1 1` |

- Frequencies are \\( \omega\_j = j\\,\Delta \\) (\\( j \ge 1 \\)); the value at \\( \omega = 0 \\) is the linear extrapolation \\( 2y\_1 - y\_2 \\) from \\( j = 1, 2 \\).
- k-points are Monkhorst–Pack type \\( -\tfrac12 + (2l+1)/(2N) \\) (0 if \\( N = 1 \\)), and \\( \pm\mathbf k \\) are merged into one point with weight 2.
- The momentum matrix elements of approach 0 are \\( p^{a}\_{mn}(\mathbf k) \propto c^\dagger\_m\\,\nabla^a(\mathbf k)\\,c\_n \\) from \\( \nabla^a(\mathbf k) \\), the Bloch sum taken with phases on lattice vectors only.
- The interband terms are, with \\( \Delta\varepsilon = \varepsilon\_m - \varepsilon\_n \\) and \\( \Delta f = f\_n - f\_m \\) (occupation doubled in the non-spin-polarized case),
  \\( \sigma\_{ab} \propto \sum \frac{\Delta f}{\Delta\varepsilon}\\, p^a\_{mn}\overline{p^b\_{mn}}\\,\frac{i}{\Delta\varepsilon + \omega + i\eta} \\) and
  \\( \varepsilon\_{ab} - \delta\_{ab} \propto \sum \frac{-i}{\Delta\varepsilon}\times(\text{the same term}) \\). The prefactor is \\( e^2/(\hbar a\_0) = 4599847.99 \\) S/m divided by \\( \Omega\sum w\_{\mathbf k} \\) (for \\( \sigma \\)),
  further divided by \\( \varepsilon\_0 E\_{\mathrm h}/\hbar = 366044.288 \\) S/m (for \\( \varepsilon \\)).
  The intraband (Drude) term enters in the noncollinear path in a form using \\( -f'(\varepsilon) \\), but not in the collinear path, which removes pairs with \\( \Delta f = 0 \\) at the outset.
- Output: `System.Name.cd_re`, `.cd_im` (\\( \sigma \\), S/m), `.df_re`, `.df_im` (\\( \varepsilon \\); **not including \\( \delta\_{ab} \\)**).
  Column 1 is the energy (eV), columns 2–10 are \\( xx, xy, xz, yx, yy, yz, zx, zy, zz \\), column 11 is 1/3 of the trace; in the collinear spin-polarized case, 9 components per spin follow.
  In the files of derived quantities (refractive index, extinction coefficient, absorption coefficient, reflectivity, transmissivity), \\( \delta\_{ab} \\) is added.
- The phase of \\( \sigma \\) (orientation of \\( i \\)) and the index order of \\( p^a\overline{p^b} \\) affect the sign of the off-diagonal components. Comparison with the formulas of the public documents is [Open].

## 11. Output of quantities on the grid

### 11.1 Cube files [Public] [Convention]

[Public] Manual Sections 9 and 21: if `level.of.fileout` (0–3, default 1) is 1 or higher, the electron density and potentials are output in Gaussian cube format.

[Convention]

| File | Contents | Condition |
|---|---|---|
| `.tden.cube` | Total electron density \\( n^\uparrow + n^\downarrow \\) (valence, bohr⁻³) | level ≥ 1 |
| `.den0.cube`, `.den1.cube` | Up and down densities | Spin-polarized, noncollinear |
| `.sden.cube` | \\( n^\uparrow - n^\downarrow \\) | Same as above |
| `.dden.cube` | \\( n - n\_{\mathrm{atom}} \\) (difference from the superposition of atomic densities) | level ≥ 1 |
| `.v0.cube`, `.v1.cube` | Effective potential per spin \\( \delta V\_{\mathrm H} + V\_{\mathrm{xc},\sigma} \\) (+ external electric field) (Hartree). Includes \\( V\_{\mathrm{NA}} \\) only when `scf.ProExpn.VNA` is `off` (not included with the default `on`) | level ≥ 1 |
| `.vhart.cube` | \\( \delta V\_{\mathrm H} \\): the Hartree potential of \\( \delta n = n - n\_{\mathrm{atom}} \\) (\\( \mathbf G = 0 \\) component set to 0) (Hartree) | level ≥ 1 |
| `.vna.cube`, `.vxc0.cube`, `.vxc1.cube`, `.grid` | Neutral-atom potential, exchange-correlation potential, coordinates of grid points (Å) | level ≥ 2 |
| `.ncsden.xsf`, `.nc.xsf`, `.nco.xsf` | Noncollinear spin-density vector field, per-atom spin and orbital moments (XCrySDen) | Noncollinear |

- **Header**: lines 1–2 are comments (for molecular orbitals, the eigenvalue and \\( \mu \\)), line 3 is the number of atoms and the grid origin, lines 4–6 are the number of grid points and the grid step vector in each direction (bohr; the number of points is positive),
  followed per atom by the atomic number, the net charge (\\( Z^{\mathrm{val}} \\) minus the Mulliken electron count), and the position (bohr).
- **Data**: the first direction is outermost and the third innermost (fastest varying), 6 values per line, with a line break after each column along the third direction.
- \\( V\_{\mathrm{NA}} \\) in `.v0` is short-ranged and therefore 0 in vacuum. The reference of the potential (the arbitrary constant) is determined by the treatment of \\( \mathbf G = 0 \\) of \\( \delta V\_{\mathrm H} \\);
  when ESM ([ESM, electric field, core hole](esm-field-corehole.md)) is used, the vacuum value is meaningful.

### 11.2 Molecular orbitals and LCAO coefficients [Public] [Convention]

[Public] Manual Sections 33 and 34: for Cluster, `level.of.fileout 2` prints the LCAO coefficients and eigenvalues (Hartree) to `.out`.
With `MO.fileout on` and `num.HOMOs`, `num.LUMOs`, and for periodic systems `MO.Nkpoint` and `<MO.kpoint` (fractional coordinates), molecular orbitals near the HOMO and LUMO are output to cube files.

[Convention] The defaults of `num.HOMOs`, `num.LUMOs`, and `MO.Nkpoint` are all 1. Effective only for Cluster, Band, and Cluster-DIIS. The HOMO is the last state with occupation exceeding \\( \tfrac12 \\).
The real and imaginary parts of \\( \psi(\mathbf r) = \sum\_a c\_a\phi\_a(\mathbf r) \\) (\\( c^\dagger Sc = 1 \\), no additional normalization, unit bohr^{-3/2}) are written to
`System.Name.homo<k>_<spin>_<index>_r.cube` and `_i.cube` (without the k index for Cluster; the imaginary part only with spin-orbit coupling or noncollinear), and likewise for the LUMO.
For periodic systems, at each grid point the periodic images \\( \mathbf L \\) of the basis functions are summed with Bloch phases.

### 11.3 Partial charge density for STM images [Public] [Convention]

[Public] Manual Section 53: in the Tersoff–Hamann approximation, the charge density of the states within an energy window measured from \\( \mu \\) is output to `System.Name.pden.cube`.
`partial.charge on`, `partial.charge.energy.window` (eV; positive for the conduction-band side, negative for the valence-band side), and `Dos.fileout` (or `DosGauss.fileout`) are required.

[Convention] The default window is 0.0 eV. Usable only with Cluster and Band. A density matrix is built with the per-state weight \\( \lvert f(\varepsilon) - f(\varepsilon - \Delta E)\rvert \\) (\\( \Delta E \\) is the window, \\( f \\) the Fermi distribution at the calculation temperature)
and converted to a density on the grid. It is doubled in the non-spin-polarized case; in the spin-polarized case `.pden0.cube` and `.pden1.cube` are also output.

### 11.4 Work function [Public] [Convention]

[Public] Manual Section 50: for a slab model, \\( W = \phi\_\infty - \mu \\) from the vacuum value \\( \phi\_\infty \\) of the in-plane averaged electrostatic potential and \\( \mu \\)
(based on \\( E[N-1] - E[N] = -\mu \\); an approximation for gapped systems). For polarized surfaces, ESM is required.

[Convention] The in-plane average is computed with the auxiliary program gcube2oned: given a cube file (usually `.v0.cube`) and an axis index, it outputs
the grid-point index, the position along the axis (Å), and the in-plane average (in the cube's units).
When axis 2 or 3 is chosen, column 1 of the output (the index) does not advance correctly (the position in column 2 is correct).

## Requirements for implementations [Design]

- **Input in sparse block format**: all analyses are stateless functions taking as input the structure (atomic positions, lattice, neighbor list) and \\( H^\sigma(\mathbf L) \\), \\( S(\mathbf L) \\),
  \\( \rho^\sigma(\mathbf L) \\), and if needed \\( \mathbf D(\mathbf L) \\), \\( \nabla(\mathbf L) \\) in [sparse block format](../design/sparse-block-format.md). `.scfout` is confined to one place as a converter (reader and writer) to and from these.
- **Common components**: the Bloch sum, k-point generation (Γ-centered, Monkhorst–Pack, paths), and the generalized eigenvalue problem (regularization rule selected by argument; the existing implementation's rules are A′ for DOS and bands, none for jx and polB)
  use those of [k-points and Bloch sums](kpoints-bloch.md) and [generalized eigenvalue problem](generalized-eigen.md).
- **Per-analysis functions** (with input/output units stated explicitly, atomic units internally):
  - Mulliken: \\( (\rho, S) \to M\_{I\alpha}^\sigma \\). Make it possible to also return the per-state weights \\( w\_{a\mu} \\).
  - Voronoi: integrate with \\( w\_I \\) from the density grid and atomic positions. Make the number of compositions of \\( p \\) (default 6, TN9's 3 as an option) and the neighbor range (neighbor list/all atoms) arguments.
  - ESP: make the sampling shells, the constraint (charge only/charge and dipole), and the total charge \\( Q\_{\mathrm{tot}} \\) (default 0, `scf.system.charge` as an option) arguments.
  - DOS: tetrahedron-method and Gaussian DOS/PDOS from arrays of eigenvalues and weights. The energy grid is given as an argument.
  - jx: \\( J \\) from \\( (H^\uparrow, H^\downarrow, S) \\) and a list of pairs. Make the normalization (Section 6.2, with Cluster's \\( \tfrac12 \\) as a compatibility option), periodic sum/individual, and the number of poles arguments.
  - Berry phase: \\( \varphi \\) from \\( (H, S, \mathbf D) \\) and a string specification. Make the link evaluation selectable (first-order expansion + symmetrization / higher order). The Chern number and \\( \mathbb Z\_2 \\) use the same link function.
  - Optical response: \\( \sigma\_{ab}(\omega) \\), \\( \varepsilon\_{ab}(\omega) \\) from \\( (H, S, \nabla) \\).
- **Handling quirks of the existing implementation**: quirks needed for comparison (the \\( \tfrac12 \\) of Cluster jx, the 6 compositions of Voronoi, the default ESP scale factor \\( \lambda\_2 = 5.0 \\), the \\( \mathbf b \\) of Z2FH, the area of calB, the negative sign of the down-spin DOS) are
  reproduced via compatibility options, while the defaults are the physically correct (or public-document-conforming) choices. Which was used is recorded in the output.
- **Output formats**: the `.Band`, `.Dos.val`, `.Dos.vec`, and cube formats are reproduced for compatibility. When adding new formats (e.g., HDF5 or JSON), units and conventions are written in the file.

## Edge cases [Convention] [Design]

- **Berry phase in metals**: if the number of occupied states changes between k-points, the polarization cannot be defined. The existing implementation warns and continues, but a new implementation treats this as an error by default.
- **Degenerate eigenvalues**: the jx periodic-sum and Cluster formulas become 0/0 when \\( \varepsilon^\downarrow\_\nu = \varepsilon^\uparrow\_\mu \\) (e.g., when an atom without exchange splitting has degenerate levels).
  Pairs with small \\( \lvert\Delta\varepsilon\rvert \\) are replaced by the limiting value \\( -f'(\varepsilon) \\) [Design].
- **Non-spin-polarized or noncollinear input to jx**: the existing implementation reads nonexistent components without checking. A new implementation treats this as an error.
- **\\( T = 0 \\)**: the jx and occupation formulas assume \\( T > 0 \\).
- **Tetrahedron method for DOS**: the formula is undefined for tetrahedra with equal vertex energies (\\( e\_{ij} = 0 \\)) and at points where \\( E \\) exactly equals a vertex energy.
  A new implementation handles degenerate tetrahedra with the limiting formula.
- **ESP**: elements without a registered van der Waals radius, no sample points within the shells, periodic systems (large effect of truncating the sum over periodic images).
- **Voronoi**: an isolated atom whose neighbor list contains only itself (\\( w\_I = 1 \\)).
- **Version and byte order of `.scfout`**: versions other than 3, byte-swapped files (a reader can read them by reversing), files without the input-file section.
- **Band paths**: a path with 1 point.

## Tests [Design]

`tools/check_analysis.py` confirms the following on small models (all pass).

| Kind | Contents | Tolerance |
|---|---|---|
| Invariant | Mulliken sum rule \\( \sum\_I M\_I = N\_{\mathrm e} \\) on a nonorthogonal 1D chain | \\( 10^{-12} \\) |
| Invariant | The per-state Mulliken weights sum to 1, and their occupation-weighted sum equals the gross population of \\( \rho S \\) | \\( 10^{-12} \\) |
| Invariant | Fuzzy-cell weights sum to 1 (3 and 6 compositions) | \\( 10^{-14} \\) |
| Analytic | Simple-cubic s band (12³ grid): the integral of the tetrahedron DOS is 1, and the integral up to \\( E\_F \\) agrees with counting on a fine grid | \\( 2\times10^{-3} \\), \\( 5\times10^{-3} \\) |
| Invariant | The integral of the Gaussian DOS is 1 | \\( 10^{-6} \\) |
| Analytic | Rice–Mele chain: the Berry phase with lattice-vector-only phases and \\( e^{-i\mathbf b\cdot\mathbf r} \\) links agrees with Resta's supercell formula (\\( \bmod 1 \\)) | \\( 10^{-10} \\) |
| Invariant | Gauge invariance of the Berry phase (arbitrary phases on eigenvectors) | \\( 10^{-12} \\) |
| Analytic | Thouless pump: the change in polarization over one cycle is 1 | \\( 10^{-6} \\) |
| Analytic | FHS Chern number of the QWZ model is \\( \pm1 \\), 0 | \\( 10^{-9} \\) |
| Analytic | 3-site nonorthogonal cluster: \\( -\partial^2\Omega/\partial\theta\_i\partial\theta\_j = 2J \\) (magnetic force theorem) | \\( 10^{-6} \\) |
| Invariant | The continued-fraction poles reproduce the Fermi function | \\( 10^{-12} \\) |
| Invariant | On a periodic chain, \\( \sum\_{\mathbf L}J\_{i\mathbf 0,j\mathbf L} \\) (pole sum) agrees with the periodic-sum formula | \\( 10^{-10} \\) |

Comparison with the existing implementation (fixtures, to be created):

| Target | Contents | Approximate tolerance |
|---|---|---|
| `.scfout` | Read the `.scfout` of the distributed input examples (methane Cluster, Si Band, bcc Fe collinear, GaAs noncollinear SOC); the rewritten file matches byte for byte (excluding content after the NUL in item 22 and the first element of 4-element vectors) | Exact |
| `.scfout` | \\( g\_s\sum\rho S = N\_{\mathrm e} \\) from the read \\( H, S, \rho \\). The relation for \\( \mathbf D \\): \\( \mathbf D\_{ab}(\mathbf L) + (\mathbf R\_I - \mathbf R\_J - \mathbf L)S\_{ab}(\mathbf L) = \mathbf D\_{ba}(-\mathbf L) \\) (indices transposed) | \\( 10^{-6} \\) (grid-integration error) |
| Mulliken | Per-atom and per-orbital values in `.out` | \\( 10^{-6} \\) |
| Voronoi | `.out` with `Voronoi.charge on` (methane, Si) | \\( 10^{-6} \\) |
| ESP | The manual example (methane, `-c 0 -s 1.4 2.0`): C \\( \approx -1.0696 \\), H \\( \approx 0.2674 \\) | \\( 10^{-4} \\) |
| DOS | DosMain tetrahedron and Gaussian output for Si (`Dos.Kgrid 12 12 12`) | \\( 10^{-5} \\) (relative) |
| Band | Eigenvalues in `.Band` | \\( 10^{-8} \\) Hartree |
| jx | Periodic sum and individual values of \\( J \\) for bcc Fe (varying `Num.Kgrid` and `Num.Poles`) | \\( 10^{-3} \\) meV |
| polB | Dipole (Debye) of NaCl (manual Section 44), and \\( Z^\*\approx 1.05 \\) from the difference when Na is displaced by ±0.05 Å | \\( 10^{-6} \\) Debye |
| CDDF | \\( \varepsilon\_{xx}(\omega) \\) of Si (manual Section 62) | \\( 10^{-4} \\) (relative) |

## [Open]

- Whether the \\( \tfrac12 \\) in the Cluster mode of jx is intended (Section 6.3). Confirm with the authors, including the correspondence with the prefactors in the public documents (Han–Ozaki–Yu 2004, manual Eq. 42).
- The sign of the off-diagonal components of the optical conductivity, the treatment of the intraband term (absent in the collinear case), and details of the matrix elements of approach 1 (treatment of the nonlocal term). Comparison with the equation numbers of the two technical notes remains to be done.
- Conditions under which the link \\( \mathbf b \\) of Z2FH (Section 8.2) affects the result.
- The Berry-phase path during SCF (Section 7.3) and its correspondence to keywords.
- Details of the contents of `.Dos.vec` for `DosGauss.fileout`, and the formats of `mode 5` and `mode 6` for DC and NEGF.
- Details of the formula for projecting the noncollinear DOS onto local spin axes (Section 3.2).
- Automatic generation of Band k-paths (manual Section 22.1), Fermi-surface output (`FermiSurfer.fileout`, Section 65), output for BoltzTraP (Section 63),
  band unfolding (Section 55), COHP, natural orbital analysis (Section 36), CWF/MLWF (Section 47), XPS (Section 60), and ionization potentials (Section 61) have not been investigated.
- Which state items 12 and 17 of `.scfout` contain when the SCF did not converge during geometry optimization or MD.
- Details of the columns of the output files of each kSpin mode (Section 9 is an overview only).
