# Parallelization

> Status: first version (awaiting review).

## Policy [Design]

- Parallelization does not change the computed results. It is therefore not made a [Convention] item.
- Parallelization of the new implementation is designed from scratch.
- The internal MPI and GPU structure of the existing implementation is not described in this note.
- The goal is that changing the number of ranks or threads yields the same results within tolerance.
- The units of parallelization are determined based on the structure (the list of atom pairs) of the [sparse block format](sparse-block-format.md).

## Published methods [Public]

Only the ideas are summarized here in our own words. Refer to the original papers for details.

### Atom partitioning and domain decomposition

Key points of T. V. T. Duy and T. Ozaki, Comput. Phys. Commun. 185, 777 (2014) (arXiv:1209.4506).

- Atoms are sorted by their position along the principal axis and partitioned so that nearby atoms go to the same process.
- The partitioning of real-space grid points is matched to the partitioning of atoms. If the grid points belonging to an atom are on the same process, communication is reduced.
- Because the workload per atom is uneven, the load is balanced by recursive bisection.
- The 3D FFT of the Poisson equation is performed with a scheme that requires little communication.
- The paper reports a parallel efficiency of about 68% when going from 16,384 to 131,072 cores on the K computer.

### Minimum-communication decomposition for multidimensional FFT

Key points of T. V. T. Duy and T. Ozaki, Comput. Phys. Commun. 185, 153 (2014) (arXiv:1302.6189).

- Multidimensional data are treated as a 1D sequence and distributed evenly among processes.
- The number of decomposed dimensions can be changed according to the number of processes.
- The communication volume depends on the order of transposes. The orders that minimize the communication volume are given for 3 to 5 dimensions.

These serve as references for the design. However, the new implementation is not bound to these schemes.

## Parallelization design of the new implementation [Design]

### Levels of parallelization

| Level | What is parallelized | Nature of communication |
|---|---|---|
| 1 | Spin and k points | Nearly independent. The density matrix is summed at the end |
| 2 | Atom pairs (two-center integrals, assembly of matrix elements) | Only sharing of atomic positions and basis information |
| 3 | Domains of the real-space grid (density, potential, XC) | Exchange of boundary grid points |
| 4 | FFT (Poisson) | All-to-all transposes |
| 5 | Eigenvalue problem | Depends on the method (below) |
| 6 | Local problems of O(N) methods | Independent. Only boundary information is exchanged |

A separate communicator is used for each level.
A rule determines the assignment to each level according to the number of ranks.

### Eigenvalue problem

- With many k points: assign k points (and spins) to groups of ranks, and solve one eigenvalue problem within each group. A group size of 1 is the basic form.
- With few k points and large matrices: distribute the matrix and use a distributed dense eigensolver.
- Which to use is decided from the number of k points, the matrix dimension, and the number of ranks. The selection rule is [Open].
- In O(N) methods, the local problem of each atom is solved independently. Only the data of overlapping regions are received from neighboring ranks.

### Distribution of the sparse block format

- The unit of a block is an atom pair \\( (I, J, \mathbf L) \\). Blocks are not split.
- Atoms are assigned to ranks (atom partitioning); the rank that owns the atom of row \\( I \\) holds all blocks of that row.
- By Hermiticity the block \\( (J, I, -\mathbf L) \\) resides on another rank, but since the format stores duplicates, each rank builds its values itself without communication.
- The structure (list and layout of atom pairs) is built deterministically by all ranks from the same rule. Only the values are held per rank.
- Grid quantities are split by domain decomposition, chosen to overlap the atom partitioning as much as possible. This raises the fraction of grid points needed for integrating atom-pair matrix elements that are local.
- The load is estimated from the number of orbitals per atom and the number of neighboring atom pairs. Repartitioning is done only when the structure changes (geometry optimization and molecular dynamics steps).

### Threads within a rank

- Parallelization within a rank uses threads over atom-pair loops or blocks of grid points.
- Sums are formed by adding per-thread partial sums in a fixed order.
- Room is left for processing products of small dense blocks together (in batches) rather than with threads.

### GPU offload boundary

- Candidates are the dense eigenvalue problem, products of small dense blocks per atom pair, FFT on the grid, and pointwise evaluation of XC.
- The values of the sparse blocks and the density and potential on the grid stay on the device during SCF iterations.
- Only scalar quantities and the minimum data needed for mixing are copied between host and device at every iteration.
- A CPU path that gives the same results is always kept, even without a GPU.

### Determinism and reproducibility

- With the same numbers of ranks and threads, the goal is to obtain the same results (bitwise identical) on every run. No sums depend on the arrival order of messages.
- With different numbers of ranks or threads, the order of summation changes, so bitwise identity is not required. Agreement within tolerance suffices.
- Sums (reductions) are fixed to one of the following:
  - A tree-structured sum in a globally defined order (e.g. by atom index).
  - Adding partial sums in a fixed order.
- If random numbers are used, the random sequences are generated from per-element (atom, k point) seeds so that they do not depend on the number of ranks.
- Quantities that must be the same on all ranks, such as the mixing history, are computed on one rank and broadcast. They are not computed separately on each rank.

### MPI and numerical libraries

- Communication uses MPI through rsmpi. Collective communication between ranks is placed behind a thin abstraction layer.
- The abstraction layer has a single-rank implementation (no communication). This lets non-parallel tests go through the same code.
- Dense matrix operations use faer and tenferro. For the distributed eigensolver, an external library is placed behind this layer.

## Tests

The approach of the [test hierarchy](../../../testing/strategy.md) is followed.

| Check | Content | Tolerance |
|---|---|---|
| Independence of rank count | Compare total energy, forces, eigenvalues, and density matrix with 1, 2, 4, 8 ranks | Total energy about 1e-8 Ha (provisional) |
| Independence of thread count | Compare likewise with 1, 2, 4 threads | Same as above |
| Reproducibility | Run twice under the same conditions; results agree | Bitwise identical |
| Independence of partitioning | Compare with different methods of atom partitioning and k-point partitioning | Same as above |
| Single-rank path | Compare the single-rank implementation of the communication abstraction layer with the MPI implementation | Same as above |
| CPU vs. GPU | Compare results for the same input | Provisionally 1e-7 Ha |
| Load balance | Max/mean ratio of the numbers of atoms, atom pairs, and grid points per rank | Provisionally 1.2 or less |
| Communication volume | Number of grid points exchanged in domain decomposition and number of atom pairs | Recorded only (used to detect performance regressions) |
| Scaling | Time when increasing the number of ranks (strong scaling) | Recorded only |

Also check running a small system (a few atoms) on many ranks. Confirm that the program does not hang when a rank has no work (an empty rank).

## Open [Open]

- The rules for rank assignment at each level, and the rule for selecting the eigenvalue method.
- The external library used for the distributed eigensolver.
- The concrete atom-partitioning algorithm (a variant of recursive bisection, or a space-filling curve).
- The targets of GPU offload (which kernels to implement first) and the GPU libraries to use.
- The tolerance values (the table above is provisional).
- The effect on load of the mismatch between atom partitioning and grid partitioning in force and stress calculations.
