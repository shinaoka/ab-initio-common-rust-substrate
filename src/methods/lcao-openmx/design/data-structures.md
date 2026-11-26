# Data structures for a parallel SCF [Design]

> Status: first draft. This chapter turns [Parallelization](parallelization.md) and the
> [sparse block format](sparse-block-format.md) into concrete data structures for the first Rust implementation.
> It does not describe the data layout of OpenMX or its GPU version.

## 1. Scope of the first implementation

The first target is an SCF that reproduces the total energy of the reference (the GPU version of OpenMX built for the CPU)
within \\( 10^{-7} \\) Ha on small systems:

- solvers `Cluster` and `Band`; non-spin-polarized and collinear spin;
- LDA (CA-PZ) and GGA-PBE; partial core correction;
- neutral-atom potential by projector expansion (the default) and, as a first step, on the grid;
- mixing `Simple`, `RMM-DIIS` (density matrix) and `RMM-DIISK` (density in reciprocal space).

Forces, stress, noncollinear spin, SOC, DFT+U, O(N) solvers, NEGF, EXX, ESM, MD and GPU offload are out of scope,
but the structures below must not block them (sections 9 and 10).

Although the first target is small, **every data structure is designed for distributed memory from the start**.
A serial run is a distributed run on one rank, through the same code.

## 2. Principles

1. **Replicated structure, distributed values.** Data whose size grows like the number of atoms times a small constant
   (species tables, positions, pair lists, grid shape, k mesh) are replicated on every rank and built deterministically by
   every rank from the same input. Data whose size grows like atoms times orbitals squared, or like grid points
   (block values, grid fields, orbital values on the grid, eigenvectors) are distributed.
2. **One owner per item.** Every distributed item has exactly one owning rank, given by a partition object.
   Partitions are explicit values, not implicit conventions.
3. **Pure kernels, thin drivers.** Physics lives in functions that act on local data only and never communicate
   (for example: one overlap block from two species and a displacement; XC at one grid point).
   Drivers loop over the items a rank owns, call kernels, and communicate through precomputed plans.
4. **Communication through a trait.** Drivers are generic over a communicator trait with three implementations:
   one rank, several simulated ranks inside one process, and MPI. Most tests use the first two.
5. **Plans are values.** Every irregular communication (grid halo gather, scatter-add to the grid, reductions of density
   matrices) is described by a plan built once per structure and reused in every SCF iteration.
   Plans are inspectable and testable without running the physics.
6. **Deterministic reductions.** No sum depends on message arrival order (section 7).

## 3. Index spaces

Each index space is a distinct newtype, so that an atom index cannot be passed where a pair index is expected.

| Type | Meaning | Range |
|---|---|---|
| `SpeciesId` | Species in the input order | `0..n_species` |
| `AtomId` | Atom in the input order | `0..n_atoms` |
| `Cell` | Lattice translation \\( \mathbf n \in \mathbb Z^3 \\) of the column-side atom | small integers |
| `PairId` | Entry \\( (I, J, \mathbf n) \\) of the pair graph | `0..n_pairs` |
| `OrbId` | Global orbital index (atom-major, then the order of [PAO basis](../modules/pao-basis.md)) | `0..n_orb` |
| `GridIndex` | Linear index of a grid point in the unit cell, row-major \\( (i\_1, i\_2, i\_3) \\) | `0..N1*N2*N3` |
| `KIndex` | Point of the irreducible k list | `0..n_k` |
| `Spin` | Spin channel | `0..n_spin` (1 or 2) |
| `KTask` | Pair `(KIndex, Spin)`; the unit of work of the eigenvalue stage | `0..n_k*n_spin` |

The row-major grid order matches pencil-rs, which is used for the distributed FFT.

## 4. Replicated objects

| Object | Content | Built from |
|---|---|---|
| `SpeciesData` | Radial tables of PAO, VPS projectors, neutral-atom potential and its projector expansion, Hartree potential of the atomic density, atomic density, PCC density (and the `scf.pcc.opencore` flag), cutoff radii, valence charge; their spherical Bessel transforms. **All** PAOs stored in the file are kept, not only the selected basis, because they seed the neutral-atom projectors | PAO and VPS files, basis specification |
| `GlobalTables` | Quantities that depend on all species together: the maximum angular momentum of the neutral-atom expansion (largest \\( l \\) of all species plus 6) and the tables built with it | `SpeciesData` |
| `Structure` | Lattice vectors, atomic positions, `AtomId -> SpeciesId`, initial up and down electron numbers per atom | Input |
| `OrbitalLayout` | Number of orbitals per atom and the offset of each atom in `OrbId` (prefix sum) | `Structure`, `SpeciesData` |
| `PairGraph` | Pair list in CSR form by row atom (section 4.1) | `Structure`, cutoff radii |
| `GridSpec` | Division numbers, lattice vectors, origin | Energy cutoff, `Structure` ([grid](../modules/grid-poisson-xc.md)) |
| `KMesh` | Fractional coordinates and integer weights of the k points | Input ([k points](../modules/kpoints-bloch.md)) |
| `Partitions` | `AtomPartition`, `GridPartition`, `KTaskPartition` (section 5) | All of the above, number of ranks |

All replicated objects are immutable during an SCF run. A new geometry produces new objects;
nothing is updated in place. Some choices made at the first geometry must survive later geometries and are passed
forward explicitly as `FrozenSettings`: the grid division numbers and origin (for clusters, the origin is the mean
atomic position of the first geometry; [grid](../modules/grid-poisson-xc.md), section 1.4) and the screening
parameters \\( s\_{IJ} \\) of the core-core energy ([total energy](../modules/energy-force.md)).
For clusters, `GridSpec` depends on cutoff radii and the energy cutoff as well as on `Structure`.

### 4.1 Pair graph

```text
row_ptr[I] .. row_ptr[I+1]          pairs whose row atom is I
pair_col[p]       : AtomId          column atom J
pair_cell[p]      : Cell            lattice translation n of J
pair_transpose[p] : PairId          the pair (J, I, -n)
block_offset[p]   : usize           start of block p in a value array (row-owner local)
```

- Inside one row, pairs are sorted by `(J, n)`. The order is a pure function of the structure, so every rank builds
  the same graph and the same `PairId`s.
- `pair_transpose` is stored because Hermiticity tests, the density matrix and mixing all need it.
- The neighbour search uses cell lists; its result is sorted before use, so the search algorithm does not affect `PairId`.
- One graph suffices. The third atoms \\( K \\) of the nonlocal and neutral-atom projector terms of pair \\( (I, J+\mathbf L) \\) are selected
  with the **PAO** cutoff radii: \\( K \\) must be a neighbour of \\( I \\) and of \\( J+\mathbf L \\) in this graph
  ([nonlocal pseudopotential](../modules/nonlocal-pseudopotential.md), section 3.3; [neutral atom](../modules/neutral-atom.md), section 3.6).
  The list is the intersection of two replicated rows, so it costs no communication.
- The damping function of the nonlocal and projector terms is applied per pair \\( (I, J+\mathbf L) \\) after the sum over \\( K \\),
  in the driver, not inside the per-\\( K \\) kernel. It is not applied to the same-atom neutral-atom blocks.
- `block_offset` depends on the rank (it indexes the local value array) and is therefore stored in the distributed
  container, not in the replicated graph.

## 5. Partitions

| Partition | Items | First rule | Later |
|---|---|---|---|
| `AtomPartition` | Row atoms \\( I \\) (and with them all pairs of row \\( I \\)) | Contiguous ranges of atoms after sorting by a space-filling curve, balanced by the cost \\( \sum\_{p \in \text{row } I} n\_I n\_J \\) | Recursive bisection, repartition on geometry change |
| `GridPartition` | Grid points | The pencil decomposition of pencil-rs | Matched to `AtomPartition` to reduce halos |
| `KTaskPartition` | `KTask`s | Round robin over k groups of size 1 | Groups of several ranks with a distributed eigensolver |

A partition answers `owner(item) -> rank` and `local_items(rank) -> range or list`.
Every distributed container holds a shared reference to its partition, so two containers can be checked for compatibility.

## 6. Distributed objects

### 6.1 `BlockMatrix<T>`

Values of a real-space matrix \\( O(\mathbf L) \\) for the rows owned by this rank.

```text
graph      : shared PairGraph
partition  : shared AtomPartition
n_spin     : 1 or 2   (collinear: one value array per spin)
values     : [n_spin][sum over local pairs of n_I * n_J] T   (each block column-major)
```

- Used for \\( S, T, V\_{\mathrm{NL}}, V\_{\mathrm{NA}} \\), \\( H\_0 \\), \\( H \\), the density matrix \\( \rho \\) and the energy-weighted density matrix.
- Operations: elementwise linear combination, inner product over all ranks (for mixing), Hermiticity check through
  `pair_transpose`, conversion to and from the sparse block format of the book (for fixtures and output).
- `T` is `f64` in this scope. Complex blocks for SOC and noncollinear spin are a later addition and do not change the layout.

### 6.2 `GridField`

A scalar field on the grid, distributed by `GridPartition`, stored as a pencil-rs array with an extra leading
dimension for spin. Real-space and reciprocal-space fields have distinct types (`RealField`, `RecipField`),
because they have different layouts after the FFT. `RecipField` fixes one convention for the index map of
\\( \mathbf G \\) and the normalization of the forward transform, taken from the [grid](../modules/grid-poisson-xc.md) chapter;
the mixing chapter uses the unnormalized forward transform, and the conversion is a named function, not a convention left to callers.
If a half-spectrum (real-to-complex) layout is used, inner products weight the conjugate pairs explicitly.

GGA needs gradients by a 27-point polynomial fit, applied twice (gradient of the density, then divergence of
\\( \partial e/\partial \nabla n \\)). A `Halo` plan exchanges a layer of width 1 **including edge and corner points**
(non-orthogonal cells use them), once per pass and spin. The 27 by 27 fit matrix is replicated in `GridSpec`.

Constant fields built once per structure by scatter-add through `AtomSpheres`:

| Field | Use |
|---|---|
| \\( n^a \\) | Superposition of atomic densities (renormalized); reference for \\( \delta n = n - n^a \\) and the neutral-atom decomposition |
| \\( n^{\mathrm{pcc}} \\) | Partial core correction, added half per spin in XC |
| \\( n^{\mathrm{init}} \\) | Initial density: atomic densities scaled by the input up and down electron numbers (differs from \\( n^a \\)) |
| \\( V\_{NA} \\) | Neutral-atom potential on the grid, only when `scf.ProExpn.VNA off` |

### 6.3 `AtomSpheres`: the bridge between atoms and grid

The grid points inside the cutoff sphere of each owned atom \\( I \\):

```text
atom_ptr[i] .. atom_ptr[i+1]   points of the i-th owned atom
point_grid[q] : GridIndex      the grid point (folded into the unit cell)
point_disp[q] : [f64; 3]       r_q - R_I in Cartesian coordinates (unfolded)
```

- All grid integrals of this scope reduce to sums over the sphere of the row atom:
  the matrix element \\( \langle \phi\_{I\alpha} \rvert V \lvert \phi\_{J\beta}(\cdot - \mathbf L) \rangle \\) is a sum over points of sphere \\( I \\) where
  \\( \phi\_{J\beta} \\) is also nonzero, and the density contribution of pair \\( (I, J, \mathbf L) \\) lives on the same points.
  A rank therefore needs grid values only on the spheres of its own atoms.
- Orbital values \\( \phi\_{I\alpha}(\mathbf r\_q) \\) are cached per owned atom. Values of the column atom \\( \phi\_{J\beta} \\) on sphere \\( I \\)
  are evaluated on the fly from the replicated radial tables (no communication). A cache for column atoms is an
  optional optimization behind the same interface.
- Both paths call one pure kernel: orbital values of one species at a list of displacements.

### 6.4 `EigenSet`

For each local `KTask`: eigenvalues, eigenvectors (dense, column-major, \\( n\_{\mathrm{orb}} \times n\_{\mathrm{st}} \\)) and later occupations.
Eigenvalues of all tasks are gathered to every rank (small) because the chemical potential needs all of them.

- The scalar type is a parameter: real for `Cluster` (and a 1 1 1 k grid, which switches to the real path), complex for `Band`.
- The number of computed states \\( n\_{\mathrm{st}} \\) may change between iterations
  ([generalized eigenproblem](../modules/generalized-eigen.md)), so the containers are resized, not fixed at construction.
- The regularized transformation \\( X = U \tilde\Lambda^{-1/2} \\) from the decomposition of \\( S(\mathbf k) \\) is cached per local k point
  for the whole SCF and shared by both spins. `KTaskPartition` therefore assigns **k points**, with both spins of one k point
  on the same group, so that \\( S(\mathbf k) \\) is decomposed once.

### 6.5 Mixing state

The mixed quantity is a `BlockMatrix` holding \\( \rho \\) (`Simple`, `RMM-DIIS`), a `BlockMatrix` holding \\( H \\)
(the first stage of `RMM-DIIS`, before the eigenproblem), or a `RecipField` holding \\( n(\mathbf G) \\) (`Kerker`, `RMM-DIISK`).
All implement one trait:

```text
trait MixVector:
    inner product with a diagonal metric (per-G weights for Kerker / RMM-DIISK, 1 for blocks), over all ranks
    residual norm with its own normalization (per method: sum over blocks, divided by N_atom for H, by N_g or N_g N_s for G)
    elementwise preconditioner (the Kerker factor), axpy, scale, clone
```

[density mixing](../modules/density-mixing.md) fixes the metrics and normalizations; they live in the implementations, not in the history code.
The history and the small Pulay system are generic over this trait. The Gram matrix is reduced over all ranks;
the small linear system is solved on rank 0 and broadcast, so all ranks use bitwise the same coefficients.

The scalar control state (current weight, previous residual norms, band energies, stage counter) is replicated
and updated identically on every rank from reduced scalars.

The output density matrix \\( \rho^{\mathrm{out}} \\) is kept separately from the mixed input, because the energy traces use it
even when the grid density is what is mixed. At the first iteration there is no input density matrix: the input is
\\( n^{\mathrm{init}} \\) on the grid.

## 7. Communication

### 7.1 Communicator trait

```text
trait Comm: rank, size, split, barrier,
            bcast, allgatherv, alltoallv, reduce_ordered, allreduce_ordered
```

| Implementation | Use |
|---|---|
| `SelfComm` | One rank, no communication |
| `SimComm` | \\( P \\) ranks as threads in one process, with shared mailboxes. For tests of partitioning and plans without `mpirun` |
| `MpiComm` | rsmpi-rt |

The distributed FFT of pencil-rs requires an MPI communicator. The Poisson solver is therefore a separate trait with a
local implementation (pencil-fft without MPI, used with `SelfComm`, and with `SimComm` by gathering to one rank)
and a distributed implementation (pencil-fft over `MpiComm`).

### 7.2 Plans

| Plan | From | To | Used for |
|---|---|---|---|
| `SphereGather` | Grid owners | Sphere points of owned atoms | \\( \delta V\_H + V\_{xc} \\) (+ \\( V\_{NA} \\) on the grid) for matrix elements |
| `SphereScatterAdd` | Sphere points of owned atoms | Grid owners | Electron density from \\( \rho \\) |
| `BlockAllgather` | Row owners | All ranks of a k group | \\( H(\mathbf L), S(\mathbf L) \\) for the Bloch sum |
| `BlockReduceScatter` | All k groups | Row owners | Density matrix and energy-weighted density matrix |

`SphereGather` and `SphereScatterAdd` are one plan used in two directions (an alltoallv with fixed send and receive lists).
In this scope the eigenvalue stage uses replicated dense matrices per k group, so `BlockAllgather` copies all blocks.
A distributed eigensolver later replaces it with a block-cyclic redistribution without changing the other stages.

### 7.3 Determinism

- With a fixed number of ranks, results are bitwise reproducible.
- Every reduction adds partial results in a fixed order: within a rank in increasing `AtomId`, `PairId` or `KTask`;
  across ranks in increasing rank. `reduce_ordered` gathers partial results and adds them in rank order;
  it never uses an MPI reduction whose order is unspecified.
- Threads inside a rank (threadteam) produce per-thread partial sums over fixed ranges, added in range order.
- Results with different numbers of ranks or threads agree within tolerance only.

## 8. One SCF iteration: who owns what

| Stage | Input (owner) | Output (owner) | Communication |
|---|---|---|---|
| Hartree and XC potential | \\( n^{\mathrm{in}} \\) (grid) | \\( \delta V\_H, V\_{xc} \\) (grid) | FFT transposes, halo for GGA gradients |
| Grid matrix elements | Potentials (grid), `AtomSpheres` | \\( H \\) blocks of owned rows | `SphereGather` |
| \\( H = H\_0 + H\_{\mathrm{grid}} \\) | Block matrices (rows) | \\( H \\) (rows) | none |
| Bloch sum, eigenproblem | \\( H, S \\) (rows) | `EigenSet` (k groups) | `BlockAllgather` |
| Chemical potential | Eigenvalues (all ranks) | \\( \mu \\), occupations (all ranks) | allgather of eigenvalues |
| Density matrix | `EigenSet` (k groups) | \\( \rho \\) (rows) | `BlockReduceScatter` |
| Density on the grid | \\( \rho \\) (rows), `AtomSpheres` | \\( n^{\mathrm{out}} \\) (grid) | `SphereScatterAdd` |
| Mixing | \\( \rho \\) (rows) or \\( n(\mathbf G) \\) (grid) | next input | inner products, broadcast of coefficients |
| Energy terms | Grid fields, block matrices | scalars (all ranks) | `allreduce_ordered` |

\\( H\_0 = T + V\_{\mathrm{NL}} + V\_{\mathrm{NA}} \\) (projector expansion) is built once per structure by the row owners, without communication:
every kernel needs only replicated species tables and the displacement. The projector overlaps of the column atom
\\( \langle \beta\_K \vert \phi\_J \rangle \\) are recomputed by the row owner (duplicated work, no messages); a cache keyed by `PairId` is an optional optimization.

The position of the mixing step depends on the method: `RMM-DIIS` mixes \\( H \\) between "\\( H = H\_0 + H\_{\mathrm{grid}} \\)" and the
eigenproblem in its first stage, and \\( \rho \\) after the density matrix later; `Kerker` and `RMM-DIISK` mix \\( n(\mathbf G) \\) after the
output density. The SCF driver is therefore a sequence of stages with the mixing stage inserted by the method.

### 8.1 Energy terms

| Term | Inputs | Owner of the work |
|---|---|---|
| \\( U\_{\mathrm{ele}} \\) (band energy) | eigenvalues, occupations | replicated |
| \\( U\_{\mathrm{kin}}, U\_{\mathrm{nl}}, U\_{\mathrm{na}} \\) (projector part) | \\( \mathrm{Tr}\\,\rho X \\) over blocks, with the spin factor | row owners, ordered reduction |
| \\( U\_{H1}, U\_{\mathrm{xc}} \\), \\( U\_{\mathrm{na}} \\) on the grid | \\( n\_\sigma, n^a, n^{\mathrm{pcc}}, \delta V\_H \\), \\( \epsilon\_{\mathrm{xc}} \\) per spin | grid owners, ordered reduction |
| \\( U\_{\mathrm{core}} + U\_{H0} \\) (screened atom-pair terms) | pair distances, frozen \\( s\_{IJ} \\), \\( Z \\), \\( V\_H^a \\) tables | row owners over the pair graph, once per structure |
| \\( \Delta E\_{\mathrm{xc}}^{(0)} \\) (atomic-density correction) | grid term on \\( n^a + n^{\mathrm{pcc}} \\); atom-centred quadrature over neighbours | grid owners and row owners, once per structure |

[Open] The [grid](../modules/grid-poisson-xc.md) and [total energy](../modules/energy-force.md) chapters must be checked against each other on
whether \\( n^{\mathrm{pcc}} \\) enters the atom-centred term of \\( \Delta E\_{\mathrm{xc}}^{(0)} \\) before this term is implemented.

## 9. Testability [Design]

| Level | What is tested | Communicator |
|---|---|---|
| Kernels | One block, one grid point, one radial transform: analytic values, invariants, fixtures | none |
| Replicated objects | Pair graph (symmetry, transpose, independence from search order), grid spec, k mesh | none |
| Partitions and plans | Every item has one owner; gather then scatter-add is the adjoint of each other; empty ranks | `SimComm` with \\( P = 1, 2, 3, 5 \\) |
| Drivers | \\( H\_0 \\), density, matrix elements, density matrix: identical within tolerance for every \\( P \\), bitwise for repeated runs | `SimComm` |
| Stage fixtures | Each stage accepts its inputs from a fixture (for example the reference \\( \rho \\)) so that it is checked in isolation | `SelfComm` |
| End to end | SCF on the target systems against the reference | `SelfComm`; a few cases with `MpiComm` under `mpirun` |

Rules:

- Every stage of section 8 is a function from explicit inputs to explicit outputs; there is no global state and no hidden I/O.
  This is what makes fixture injection possible.
- `SimComm` with \\( P \\) larger than the number of atoms is part of the default test set, so empty ranks are always exercised.
- Tests are small and fast: two to five atoms, coarse grids, fixed seeds, exact or tight assertions.
  Full SCF regressions run periodically, not in the default test set.

## 10. Crates

The implementation is called **sphaira-rs** (from Greek *sphaira*, sphere: the cutoff spheres of the atoms carry the whole design).
The repository is `sphaira-rs`, one Cargo workspace; every crate is named `sphaira-*`.

| Crate | Content | Depends on |
|---|---|---|
| `sphaira-comm` | `Comm` trait, `SelfComm`, `SimComm`, `MpiComm` (feature) | rsmpi-rt (optional) |
| `sphaira-atomdata` | PAO/VPS parsers, real spherical harmonics, radial grids, spherical Bessel transforms | — |
| `sphaira-structure` | `Structure`, `OrbitalLayout`, `PairGraph`, `KMesh`, partitions | `sphaira-comm` |
| `sphaira-blocks` | `BlockMatrix`, block plans, Bloch sums, density matrix | `structure`, tenferro-rs |
| `sphaira-twocenter` | \\( S, T, V\_{\mathrm{NL}}, V\_{\mathrm{NA}} \\) kernels and drivers | `atomdata`, `blocks` |
| `sphaira-grid` | `GridSpec`, fields, `AtomSpheres`, sphere plans, Poisson, XC | `structure`, pencil-rs |
| `sphaira-scf` | Eigenproblem with regularization, occupations, mixing, total energy, the SCF loop | all above |
| `sphaira-cli` | Input-file reader (subset of the `.dat` keywords), output | `sphaira-scf` |

The crate boundaries follow the dependency order of the SCF, so each crate is tested before the next exists.

## [Open]

- The agreement target. The reference stores orbital values on the grid in single precision (relative error of matrix elements
  about \\( 10^{-7} \\)) and adds tiny random numbers to \\( S \\); the chapter tolerances are about \\( 10^{-6} \\) Ha. Whether \\( 10^{-7} \\) Ha
  is reachable with a double-precision implementation, or a single-precision mode for comparison is needed, is to be settled with fixtures.
- The rule that matches `GridPartition` to `AtomPartition`; until then pencil-rs default pencils are used.
- Memory of the orbital cache for large systems, and when to switch to on-the-fly evaluation for row atoms too.
- The redistribution between row-distributed blocks and a block-cyclic distributed eigensolver.
- Whether `BlockReduceScatter` should avoid the temporary full density matrix on each k group for large systems.
