# H(R)/S(R) sparse block format [Design]

> This chapter is a design draft. It does not mirror the internal data layout of OpenMX;
> it defines an interface for the new implementation.

## Target

Real-space matrices in a localized basis,

$$
O\_{I\alpha, J\beta}(\mathbf L) = \langle \phi\_{I\alpha}(\mathbf r - \mathbf R\_I) | \hat O | \phi\_{J\beta}(\mathbf r - \mathbf R\_J - \mathbf L) \rangle
$$

are stored as dense blocks per atom pair and lattice vector. \\( \hat O \\) is the overlap (identity operator),
the kinetic energy, the Hamiltonian, etc.

## Separation of structure and values

- **Structure**: the list of atom pairs \\( (I, J, \mathbf L) \\), the number of orbitals of each atom, the block shapes, and the positions in the value array.
  It is invariant as long as the atomic configuration, basis, and cutoffs do not change.
- **Values**: the numerical values of the blocks. They are updated at every SCF iteration.

The structure is built once, and only the values can be replaced.

## Conventions to be decided

| Item | Draft |
|---|---|
| Range of atom pairs | All pairs satisfying \\( \lVert \mathbf R\_J + \mathbf L - \mathbf R\_I \rVert \le r\_{c,I} + r\_{c,J} \\) (including the atom itself; the boundary is included, as in the existing implementation; the values of pairs on the boundary are 0) |
| Duplicate pairs | \\( (J, I, -\mathbf L) \\) is also stored (it can be recovered by Hermiticity, but simplicity on the user side is prioritized). However, the density matrix of the O(N) solvers (`DC`, `Krylov`, `DC-LNO`) is not Hermitian in the above sense, \\( \rho\_{I\alpha,J\beta}(\mathbf L) \ne \rho\_{J\beta,I\alpha}(-\mathbf L) \\), so it cannot be recovered and **both directions must be stored** ([O(N) methods](../modules/on-methods.md)) |
| Representation of lattice vectors | Integer triple \\( (n\_1, n\_2, n\_3) \\) |
| Ordering within a block | Rows are orbitals of \\( I \\), columns are orbitals of \\( J \\). The orbital ordering follows the convention of the [PAO basis](../modules/pao-basis.md) |
| Memory layout | Undecided (column-major is the first candidate, considering exchange with faer and tenferro) |
| Spin | Collinear: a separate matrix per spin. Noncollinear: a complex matrix with 2×2 spin blocks |
| Bloch sum | \\( O(\mathbf k) = \sum\_{\mathbf L} e^{+i\mathbf k\cdot\mathbf L} O(\mathbf L) \\). The phase is attached to the lattice vector only and does not include atomic positions (same as the existing implementation; see [Convention] below) |
| Representation of k points | Fractional coordinates \\( \mathbf q \\) with respect to the reciprocal lattice vectors. \\( \mathbf k\cdot\mathbf L = 2\pi\\,\mathbf q\cdot\mathbf n \\) |
| Noncollinear k-space matrix | \\( 2N\times 2N \\), spin outer and orbital inner (\\( p = \sigma N + a \\)). The overlap is \\( I\_2\otimes S(\mathbf k) \\) |

## k-space conventions of the existing implementation [Convention]

To align design choices with the existing implementation, the verified facts are recorded here (details in [k points and Bloch sums](../modules/kpoints-bloch.md)).

- The phase is \\( e^{+2\pi i\\,\mathbf q\cdot\mathbf n} \\), where \\( \mathbf n \\) is the integer coordinate of the lattice vector pointing to the periodic image of the column-side (ket-side) atom \\( J \\).
  Atomic positions do not enter the phase (consistent with Eqs. 1 and 28–29 of technical note 1). Multiple periodic images of the same \\( (I, J) \\) are summed.
- For Cluster (Γ point), \\( \sum\_{\mathbf L} O(\mathbf L) \\) is built as a real matrix.
- In the noncollinear case, the three blocks \\( \uparrow\uparrow \\), \\( \downarrow\downarrow \\), \\( \uparrow\downarrow \\) are built by Bloch sums,
  and the \\( \downarrow\uparrow \\) block is filled with the Hermitian conjugate of the \\( \uparrow\downarrow \\) block.
  In real space, the real part of each spin block and the imaginary part due to spin-orbit coupling etc. are stored separately.
- The real-space density matrix is \\( \rho(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k}\sum\_\mu f\_{\mu\mathbf k}\\, e^{+i\mathbf k\cdot\mathbf L}\\, \overline{c\_{I\alpha,\mu}}\\, c\_{J\beta,\mu} \\)
  (the real part in the collinear case), with the same phase direction as \\( H(\mathbf L) \\).

## Open issues [Open]

- How to store the real and imaginary parts and the spin components in the noncollinear case (the existing implementation keeps the real and imaginary parts in separate arrays; whether the new format should use complex blocks).
- In collinear calculations, whether to leave the k-point reduction by time reversal to the user side (density-matrix assembly), or to provide a function on the format side that takes \\( \mathrm{Re} \\).
- Whether to provide, as a conversion function, a phase convention that includes atomic positions (in case post-processing needs it).
- Under MPI distribution, which rank holds which atom pairs.
