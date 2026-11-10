# Exact exchange, hybrid functionals, and dispersion correction

> Status: first draft (awaiting review). For the dispersion correction (DFT-D2, DFT-D3), the methods of Grimme et al.'s papers and Section 54 of the manual
> are rewritten in the notation of this note, and conventions absent from the public documents were confirmed against the existing implementation.
> For exact exchange and hybrid functionals, this chapter organizes the methods in the public documents and confirms that, in the existing implementation (a derivative of OpenMX 4.0),
> they cannot be reached from the input. Numerical checks of the formulas are in `tools/check_dispersion.py`.

## Purpose and scope

This chapter specifies the dispersion correction (output label `UvdW`), a correction that depends only on atomic positions and not on the electron density,
and the handling of `HSE` and `NLX` accepted by `scf.XcType`.

- [Public] The dispersion corrections are Grimme's DFT-D2 (JCC **27**, 1787 (2006), hereafter G06) and DFT-D3
  (JCP **132**, 154104 (2010), hereafter G10; BJ damping: Grimme–Ehrlich–Goerigk, JCC **32**, 1456 (2011), hereafter G11).
  Section 54 of the manual (54.1 DFT-D2, 54.2 DFT-D3) describes the input.
- [Convention] For exact exchange (Fock exchange) and hybrid functionals, the examined version has **no implementation reachable from the input**.
  The public documents (such as Toyoda–Ozaki's method for electron repulsion integrals) are organized here as a reference for future extensions.
- [Design] How the term enters the total energy, and the overall form of forces and stress, follow the [Total energy and forces](energy-force.md) chapter.
  This chapter specifies the contents of its \\( U\_{\mathrm{vdW}} \\) term.

## Notation [Design]

In addition to the symbols of [Notation](../../../conventions/notation.md) and [Total energy and forces](energy-force.md), the following are used.

| Symbol | Meaning |
|---|---|
| \\( Z\_I \\) (this chapter only) | Atomic number of atom \\( I \\) (not the valence charge of the pseudopotential) |
| \\( \mathbf d\_{IJ\mathbf L} = \mathbf R\_I - \mathbf R\_J - \mathbf L \\), \\( d = \lVert \mathbf d\_{IJ\mathbf L}\rVert \\) | Interatomic vector (the sign-reversed \\( \boldsymbol\tau\_{IJ\mathbf L} \\) of [Total energy and forces](energy-force.md)) and distance |
| \\( \mathbf a\_k \\) (\\( k = 1,2,3 \\)) | Lattice vectors. For isolated systems (`Cluster`), the vectors of the computational box |
| \\( h\_k = \Omega / \lVert \mathbf a\_i\times\mathbf a\_j\rVert \\) | Interplanar spacing along \\( \mathbf a\_k \\) (\\( (i,j,k) \\) cyclic) |
| \\( p\_I \in \\{0, 1\\} \\) | Per-atom periodicity flag (`DFTD.periodicity`) |
| \\( r\_{\mathrm{cut}} \\), \\( r\_{\mathrm{cn}} \\) | Cutoffs for the pair sum and the coordination number (`DFTD.rcut_dftD`, `DFTD.cncut_dftD`) |
| \\( \mathrm{CN}\_I \\) | DFT-D3 coordination number |
| \\( \Sigma\_{\gamma\eta} = \partial U / \partial\epsilon\_{\gamma\eta} \\) | Strain derivative (same as the definition of stress in [Total energy and forces](energy-force.md)) |

Energies are in Hartree and lengths in bohr.

## Input keywords and defaults [Convention] [Public]

[Public] Section 54 of the manual describes the meaning of the keywords and most of the defaults.
[Convention] The defaults in the table are the values confirmed in the existing implementation. Those that differ from the manual are annotated.

| Keyword | Value | Default | Meaning |
|---|---|---|---|
| `scf.dftD` | on, off | off | Enable the dispersion correction |
| `version.dftD` | 2, 3 | 2 | DFT-D2 or DFT-D3. Any other value is an input error (stop) |
| `DFTD3.damp` | `zero`, `bj` | **`zero`** (note 1) | DFT-D3 damping function |
| `DFTD.Unit` | `Ang`, `AU` | **`AU`** (note 2) | Unit of `DFTD.rcut_dftD` and `DFTD.cncut_dftD` |
| `DFTD.rcut_dftD` | real | 100 (in the unit of `DFTD.Unit`) | Pair-sum cutoff \\( r\_{\mathrm{cut}} \\) |
| `DFTD.cncut_dftD` | real | 40 (in the unit of `DFTD.Unit`) | Coordination-number cutoff \\( r\_{\mathrm{cn}} \\) (used by D3 only) |
| `DFTD.IntDirection` | 3 integers | `1 1 1` | Exclude periodic images along lattice directions set to 0 from the sum |
| `<DFTD.periodicity` … `DFTD.periodicity>` | each line "atom index, flag" | 1 for all atoms | Per-atom periodicity flag \\( p\_I \\) |
| `DFTD.d` | real | 20 | Steepness \\( d \\) of the D2 damping |
| `DFTD.scale6` | real | D2: 0.75, D3: 1.0 | \\( s\_6 \\) |
| `DFTD.scale8` | real | note 3 | D3 \\( s\_8 \\) |
| `DFTD.sr6` | real | note 3 | \\( s\_{r,6} \\) of D3 zero damping |
| `DFTD.a1`, `DFTD.a2` | real | note 3 | \\( a\_1 \\), \\( a\_2 \\) of D3 BJ damping (\\( a\_2 \\) is always in bohr; not affected by `DFTD.Unit`) |

- Note 1: Section 54.2 of the manual gives the default as `bj`, but the existing implementation uses `zero` when `DFTD3.damp` is absent.
  [Design] The reimplementation follows the existing implementation and defaults to `zero`; inputs used for comparison always specify it explicitly.
- Note 2: The manual does not state the default explicitly (the example in Section 54.1 uses `Ang`; Section 54.2 recommends `AU`).
  The existing implementation assumes `AU` when the keyword is absent. Hence the default cutoffs are 100 bohr and 40 bohr.
  With `Ang`, values are converted from Å to bohr (the conversion factor is the bohr–Å value of [Notation](../../../conventions/notation.md)).
- Note 3: The remaining D3 defaults depend on `scf.XcType`.
  For `GGA-PBE`, \\( s\_{r,6} = 1.217 \\), \\( a\_1 = 0.4289 \\), \\( a\_2 = 4.4407 \\), \\( s\_8 = 0.722 \\) (zero) or \\( 0.7875 \\) (bj),
  the PBE values of G10 and G11. For other functionals, \\( s\_8 = 0.5 \\), \\( s\_{r,6} = 1.0 \\), \\( a\_1 = 0.5 \\), \\( a\_2 = 5.0 \\);
  this is not a set published for any particular functional (the manual also requires the values to be specified for functionals other than PBE).
- The D2 values \\( s\_6 = 0.75 \\), \\( d = 20 \\) are the PBE values of G06. D2 does not change its defaults according to `scf.XcType`.
- In each line of `DFTD.periodicity`, the first column is the atom index in `Atoms.SpeciesAndCoordinates` and the second column is the flag (1: periodic, 0: nonperiodic).
  When the block is present, as many lines as there are atoms are read.
- The D2 correction uses only the D2 keywords (`DFTD.d`, the D2 `DFTD.scale6`), and D3 only the D3 ones.

## Pair sum and periodic images [Convention]

### Range of the sum

[Convention] The dispersion energy is a sum over triples of an atom \\( I \\) in the unit cell, an atom \\( J \\) (including \\( J = I \\)), and a lattice vector
\\( \mathbf L = n\_1\mathbf a\_1 + n\_2\mathbf a\_2 + n\_3\mathbf a\_3 \\). A triple \\( (I, J, \mathbf L) \\) is included when all of the following hold.

1. \\( 0.1 < d^2 \\) (bohr²). This excludes the self term \\( J = I \\), \\( \mathbf L = 0 \\);
   two atoms closer than 0.32 bohr are also silently excluded.
2. \\( d \le r\_{\mathrm{cut}} \\) (D2), \\( d < r\_{\mathrm{cut}} \\) (D3 pair sum). The cutoff is not smooth, and the energy is discontinuous there.
3. \\( \mathbf L = 0 \\), or \\( p\_J = 1 \\).
4. \\( \lvert n\_k\rvert \le N\_k \\). \\( N\_k \\) is determined by the rules below; if the \\( k \\)-th component of `DFTD.IntDirection` is 0, \\( N\_k = 0 \\).

The weight \\( w\_{IJ\mathbf L} \\) is 1 when \\( \mathbf L \ne 0 \\) and \\( p\_I = 0 \\), \\( p\_J = 1 \\), and \\( \tfrac12 \\) otherwise.
If all atoms are periodic, this is the usual \\( \tfrac12\sum\_I\sum\_{J,\mathbf L} \\).
When the flags are mixed, this weight makes the interaction of a nonperiodic atom \\( I \\) with the periodic images of a periodic atom \\( J \\) counted exactly once,
and interactions of nonperiodic atoms with each other's periodic images are not included.

### Range of periodic images \\( N\_k \\)

- [Convention] D3 (pair sum and coordination number): \\( N\_k = \lceil r / h\_k \rceil \\) (\\( r \\) is \\( r\_{\mathrm{cut}} \\) or \\( r\_{\mathrm{cn}} \\)).
  This is the smallest parallelepiped range containing a sphere of radius \\( r \\), and misses no image inside the sphere for any lattice.
- [Convention] D2: \\( N\_k = \lceil r\_{\mathrm{cut}} / \lVert\mathbf a\_k\rVert \rceil \\). For non-orthorhombic lattices \\( \lVert\mathbf a\_k\rVert > h\_k \\), so
  some images within \\( r\_{\mathrm{cut}} \\) may be missed (example in `tools/check_dispersion.py`: for an oblique lattice only 52 of the 92 required images are included).
  With the default \\( r\_{\mathrm{cut}} \\) (100 bohr) the missed images are distant and their effect on the energy is small.
- [Convention] When atomic positions change (MD, geometry optimization, cell optimization), \\( N\_k \\) is redetermined from the lattice vectors at that time.

### Isolated systems [Convention]

For isolated systems (`scf.EigenvalueSolver Cluster`) as well, the rules above are applied **as is**, with the vectors of the computational box taken as lattice vectors.
That is, with the defaults, dispersion interactions with the periodic images of the box are included. For an isolated molecule, the correct usage is to specify `DFTD.IntDirection 0 0 0`
(the distributed DFT-D3 examples do so).

[Design] The reimplementation reproduces the behavior of the existing implementation as the default, and issues a warning for isolated systems when `DFTD.IntDirection` is not `0 0 0`.

### Empty atoms and supported elements [Convention]

- Species with atomic number 0 (empty atoms that carry only basis functions) contribute nothing in D2 (treated as \\( C\_6 = 0 \\)) and are excluded from the sums (including the coordination number) in D3.
- D2 has values only for \\( Z \le 54 \\) (H to Xe); atoms with \\( 55 \le Z \le 103 \\) are treated as \\( C\_6 = 0 \\) and silently contribute nothing. \\( Z > 103 \\) is an error.
- D3 supports \\( 1 \le Z \le 94 \\) (H to Pu); other atomic numbers are an error (stop).

## DFT-D2 [Public] [Convention]

[Public] G06 Eq. 11–13, written with the sum range above, is

$$
U\_{\mathrm{vdW}}^{\mathrm{D2}} = -\sum\_{(I,J,\mathbf L)} w\_{IJ\mathbf L}\\, s\_6\\, \frac{C\_6^{IJ}}{d^6}\\, f(d),\qquad
f(d) = \frac{1}{1 + \exp\left[-d\_{\mathrm{D}}\left(d / R\_0^{IJ} - 1\right)\right]}
$$

Here \\( d\_{\mathrm{D}} \\) is `DFTD.d` (a subscript is added to distinguish it from the symbol \\( d \\)),
\\( C\_6^{IJ} = \sqrt{C\_6^{I} C\_6^{J}} \\) (G06 Eq. 13), and \\( R\_0^{IJ} = R\_0^{I} + R\_0^{J} \\).

[Public] [Convention] The per-element \\( C\_6 \\) and \\( R\_0 \\) are the values of G06 Table 1 (\\( C\_6 \\) in J nm⁶ mol⁻¹, \\( R\_0 \\) in Å) converted to atomic units.

- \\( 1\\ \mathrm{J\\,nm^6\\,mol^{-1}} = 17.3453\\ \mathrm{Hartree\\,bohr^6} = 34.6906\\ \mathrm{Ry\\,bohr^6} \\), \\( 1\\ \text{Å} = 1.889726\\ \mathrm{bohr} \\).
- The existing implementation holds values in Ry bohr⁶ rounded to three decimal places (a conversion table also used in other codes), and divides them by 2 on reading to obtain Hartree.
  For example, for H, \\( C\_6 = 4.857 \\) Ry bohr⁶ (0.14 J nm⁶ mol⁻¹ in G06) and \\( R\_0 = 1.892 \\) bohr (1.001 Å).
  The difference due to this rounding is at most \\( 2\times10^{-3} \\) relative in \\( C\_6 \\) (confirmed for several cases in `tools/check_dispersion.py`).
- [Design] The reimplementation takes the values from G06 Table 1 and converts them with the factors above. When transcribing the table, state the source;
  for comparison with the existing implementation, either allow a relative difference of about \\( 2\times10^{-3} \\) in \\( U\_{\mathrm{vdW}} \\), or provide the rounded table as a fixture.

## DFT-D3 [Public] [Convention]

### Energy

[Public] The two-body term of G10 Eq. 3 is written with the sum range above.

$$
U\_{\mathrm{vdW}}^{\mathrm{D3}} = -\sum\_{(I,J,\mathbf L)} w\_{IJ\mathbf L}\left[ s\_6\\, C\_6^{IJ}\\, g\_6(d) + s\_8\\, C\_8^{IJ}\\, g\_8(d)\right]
$$

The distance function \\( g\_n \\) including damping is one of the following.

- Zero damping (G10 Eq. 4–5):

$$
g\_n(d) = \frac{1}{d^n}\\, \frac{1}{1 + 6\left(d / (s\_{r,n} R\_0^{IJ})\right)^{-\alpha\_n}},\qquad
\alpha\_6 = 14,\quad \alpha\_8 = 16,\quad s\_{r,8} = 1
$$

  \\( R\_0^{IJ} \\) is the table of cutoff radii per element pair (defined in G10; in Å in the distributed data; the existing implementation converts it to bohr).
- BJ damping (G11 Eq. 5–6):

$$
g\_n(d) = \frac{1}{d^n + \left(a\_1 R\_0^{IJ,\mathrm{BJ}} + a\_2\right)^n},\qquad
R\_0^{IJ,\mathrm{BJ}} = \sqrt{C\_8^{IJ}/C\_6^{IJ}}
$$

  \\( a\_2 \\) is in bohr. The PBE value 4.4407 of G11 is used as a value in bohr (consistent with the reference implementation).
  Since \\( C\_8/C\_6 = 3 Q\_I Q\_J \\) does not depend on the coordination number, neither does \\( R\_0^{\mathrm{BJ}} \\).

[Convention] The three-body term (the Axilrod–Teller–Muto-type \\( E^{(3)} \\) of G10 Eq. 14–16) is **not implemented**.
What is sometimes called "three-body" in the output or comments of the existing implementation is the force term arising from the coordination-number dependence of \\( C\_6 \\) below, not a three-body energy term.

### Determination of \\( C\_6 \\) and \\( C\_8 \\)

[Public] Coordination number (G10 Eq. 15):

$$
\mathrm{CN}\_I = \sum\_{J,\mathbf L} \frac{1}{1 + \exp\left[-k\_1\left(R\_{\mathrm{cov}}^{IJ}/d - 1\right)\right]},\qquad
R\_{\mathrm{cov}}^{IJ} = k\_2\left(R\_{\mathrm{cov}}^{I} + R\_{\mathrm{cov}}^{J}\right),\quad k\_1 = 16,\quad k\_2 = 4/3
$$

[Convention] The sum runs over triples with \\( 0.1 < d^2 < r\_{\mathrm{cn}}^2 \\) and \\( \mathbf L = 0 \\) or \\( p\_J = 1 \\) (unweighted, including the periodic images of \\( J = I \\)).
The covalent radii \\( R\_{\mathrm{cov}} \\) are the values of Pyykkö–Atsumi (Chem. Eur. J. **15**, 188 (2009)), reduced by 10% for metals (the G10 convention).

[Public] \\( C\_6 \\) is a Gaussian-weighted average of reference values (G10 Eq. 16–17):

$$
C\_6^{IJ}(\mathrm{CN}\_I, \mathrm{CN}\_J) = \frac{\sum\_{a,b} C\_{6,\mathrm{ref}}^{IJ,ab}\\, L\_{ab}}{\sum\_{a,b} L\_{ab}},\qquad
L\_{ab} = \exp\left[-k\_3\left\\{(\mathrm{CN}\_I - \mathrm{CN}^{a}\_{I})^2 + (\mathrm{CN}\_J - \mathrm{CN}^{b}\_{J})^2\right\\}\right],\quad k\_3 = 4
$$

\\( a \\) (\\( b \\)) runs over the reference systems of element \\( I \\) (\\( J \\)) (at most 5 per element), and \\( \mathrm{CN}^{a}\_{I} \\) is the corresponding reference coordination number.

[Convention]

- Only reference pairs with \\( C\_{6,\mathrm{ref}} > 10^{-12} \\) enter the sum.
- The maximum exponent among the pairs in the sum is subtracted before taking \\( \exp \\) (a factor common to numerator and denominator, so the ratio is unchanged).
  This keeps the denominator nonzero even when \\( \mathrm{CN} \\) is far from all reference coordination numbers. In the limit, the result is the \\( C\_6 \\) of the nearest reference in coordination-number space,
  consistent with the behavior of the reference implementation.
  [Open] The distributed OpenMX 4.0 itself may have set \\( C\_6 = C\_8 = 0 \\) when the denominator fell below \\( 10^{-12} \\) (a change in the examined version).
  This makes no difference for ordinary molecules and solids, but \\( U\_{\mathrm{vdW}} \\) changes for atoms far from the reference coordination numbers (e.g. highly coordinated P).
- [Public] \\( C\_8^{IJ} = 3\\, C\_6^{IJ}\\, Q\_I Q\_J \\) (the form of G10 Eq. 6). \\( Q\_I \\) is a per-element coefficient,
  included in the distributed data as the value corresponding to \\( \sqrt{\tfrac12 \sqrt{Z\_I}\\, \langle r^4\rangle\_I/\langle r^2\rangle\_I} \\).

### Source and format of the coefficient tables [Public] [Convention]

[Convention] The existing implementation embeds the following tables in its source. All are data of the DFT-D3 reference program distributed by Grimme et al. (identical to the supplementary material of G10).

| Table | Format | Unit |
|---|---|---|
| \\( C\_{6,\mathrm{ref}} \\) | 5 numbers per record: \\( C\_{6,\mathrm{ref}} \\), encoded index of element A, encoded index of element B, reference coordination number of A, reference coordination number of B. The encoded index is \\( Z + 100\\,a \\) (\\( a \\) is the 0-based reference index). About 32,000 records | Hartree bohr⁶; coordination numbers dimensionless |
| \\( R\_0^{IJ} \\) | Lower triangle (\\( Z\_A \ge Z\_B \\)) for \\( Z \le 94 \\), laid out row by row | Å (divided by 0.52917726 on reading to obtain bohr) |
| \\( R\_{\mathrm{cov}} \\) | Per element | Å (multiplied by \\( k\_2 \\) and divided by 0.52917726 on reading) |
| \\( Q \\) | Per element | bohr |

- Each record is registered for both \\( (A,a,B,b) \\) and \\( (B,b,A,a) \\) (the order of the reference coordination numbers is also swapped).
- [Design] The reimplementation does not transcribe these tables into the source; it reads the distributed data (the coefficient files of the reference program) at run time,
  or generates them from that source at build time. The table values themselves are not reproduced in this note.
  For the bohr–Å conversion, 0.52917726 is used to match the existing implementation (the relative difference from the value in [Notation](../../../conventions/notation.md) is \\( 10^{-8} \\) and does not affect comparisons).

## Forces and stress [Public] [Convention]

### Two-body term

[Public] Let the pair function be \\( \phi\_{IJ}(d) \\) (in D2, \\( -s\_6 C\_6 f/d^6 \\); in D3, \\( -(s\_6 C\_6 g\_6 + s\_8 C\_8 g\_8) \\) with \\( C\_6 \\) held fixed). Then

$$
\frac{\partial U}{\partial \mathbf R\_K}\Big\rvert\_{C\_6} = \sum\_{(I,J,\mathbf L)} w\_{IJ\mathbf L}\\, \phi\_{IJ}'(d)\\, \frac{\mathbf d\_{IJ\mathbf L}}{d}\\,(\delta\_{KI} - \delta\_{KJ}),\qquad
\Sigma\_{\gamma\eta}\Big\rvert\_{C\_6} = \sum\_{(I,J,\mathbf L)} w\_{IJ\mathbf L}\\, \phi\_{IJ}'(d)\\, \frac{d\_\gamma\\, d\_\eta}{d}
$$

The derivatives are

$$
\frac{d}{dd}\left[\frac{f}{d^6}\right] = \frac{f}{d^6}\left[\frac{d\_{\mathrm{D}}}{R\_0}\\,\frac{e}{1+e} - \frac{6}{d}\right],\quad e = \exp\left[-d\_{\mathrm{D}}(d/R\_0 - 1)\right]
$$

$$
\frac{d g\_n}{dd} = \frac{g\_n}{d}\left[-n + \alpha\_n\\, 6 p\_n \frac{1}{1 + 6p\_n}\right],\quad p\_n = \left(\frac{d}{s\_{r,n}R\_0}\right)^{-\alpha\_n}\quad(\text{zero}),\qquad
\frac{d g\_n}{dd} = -\frac{n\\, d^{n-1}}{\left(d^n + R^n\right)^2}\quad(\text{BJ})
$$

(\\( R = a\_1 R\_0^{\mathrm{BJ}} + a\_2 \\)).

### Term through the coordination number (D3)

[Public] Since \\( C\_6 \\) depends on atomic positions through the coordination number, the chain rule adds the following terms.
When all atoms are periodic, with \\( G\_{IJ} = \sum\_{\mathbf L}\left[s\_6 g\_6 + s\_8\\, 3Q\_IQ\_J\\, g\_8\right] \\) (sum over \\( \mathbf L \\) in the pair-sum range, unweighted),

$$
\frac{\partial U}{\partial \mathrm{CN}\_I} = -\sum\_J \frac{\partial C\_6^{IJ}}{\partial \mathrm{CN}\_I}\\, G\_{IJ},\qquad
\frac{\partial U}{\partial \mathbf R\_K}\Big\rvert\_{\mathrm{CN}} = \sum\_I \frac{\partial U}{\partial\mathrm{CN}\_I}\sum\_{J,\mathbf L} c'(d)\\, \frac{\mathbf d\_{IJ\mathbf L}}{d}\\,(\delta\_{KI}-\delta\_{KJ})
$$

$$
\Sigma\_{\gamma\eta}\Big\rvert\_{\mathrm{CN}} = \sum\_I \frac{\partial U}{\partial\mathrm{CN}\_I}\sum\_{J,\mathbf L} c'(d)\\, \frac{d\_\gamma\\, d\_\eta}{d},\qquad
c'(d) = -\frac{e\_c\\, k\_1 R\_{\mathrm{cov}}^{IJ}}{(1+e\_c)^2 d^2},\quad e\_c = \exp\left[-k\_1(R\_{\mathrm{cov}}^{IJ}/d - 1)\right]
$$

(the sum over \\( J,\mathbf L \\) is over the coordination-number sum range). \\( \partial C\_6/\partial\mathrm{CN}\_I \\) is obtained from the derivative of the quotient of Gaussian weights using
\\( \partial L\_{ab}/\partial\mathrm{CN}\_I = -2k\_3(\mathrm{CN}\_I - \mathrm{CN}^a\_I) L\_{ab} \\).

> **Check**: For the D2 and D3-type formulas (D3 with a fictitious reference table generated with random numbers; both zero and BJ),
> on a periodic non-orthorhombic lattice, the analytic gradient and stress agree with central differences to within \\( 10^{-9} \\),
> as confirmed in `tools/check_dispersion.py`. The same holds for the D2 formula with mixed periodicity flags.

### Conventions of the existing implementation [Convention]

- The force is added to the gradients of the other terms as \\( -\partial U\_{\mathrm{vdW}}/\partial\mathbf R\_K \\). The stress \\( \Sigma \\) is also added to the other terms
  (units and output are in the [Total energy and forces](energy-force.md) chapter).
- Since \\( U\_{\mathrm{vdW}} \\) does not depend on the electron density, it has the same value in every SCF iteration. The value appears as `UvdW` in the output.
- With mixed periodicity flags, the energy of the existing implementation follows the weighted formula above, but the gradient drops the contribution to the periodic atom \\( J \\)
  of the term "nonperiodic atom \\( I \\) and the periodic images of \\( J \\)" (it is added only on the \\( I \\) side).
  Consequently the forces do not sum to zero, and the forces are not the exact gradient of the energy. The weight of that term in the stress also remains \\( \tfrac12 \\).
  In D3, in addition, the coordination number (\\( \mathrm{CN}\_J \\) does not count the periodic images of \\( I \\)) and \\( G\_{IJ} \\) become asymmetric in \\( I, J \\).
- [Design] The reimplementation implements the formulas above (the exact gradient including the weights), and comparison with the existing implementation is restricted to the case \\( p\_I = 1 \\) for all atoms.
  With mixed flags, only the energy is compared.

## Exact exchange and hybrid functionals

### Status of the existing implementation [Convention]

The following are facts for the examined version (a derivative of OpenMX 4.0).

- `scf.XcType` accepts `HSE` and `NLX` in addition to `LDA`, `LSDA-CA`, `LSDA-PW`, and `GGA-PBE` (case-insensitive).
  Section 11 of the manual 4.0 does not mention these two.
- The exchange-correlation evaluation on the grid handles only the four functionals above. When `HSE` or `NLX` is specified, no exchange-correlation formula is evaluated,
  and the input check does not stop either. The results obtained are not meaningful calculations.
- Only for `HSE`, the following four keywords are read and auxiliary arrays are allocated:
  `scf.xmin.Rec.Coulomb` (default −14.0), `scf.xmax.Rec.Coulomb` (default 14.0), `scf.Ng1.Rec.Coulomb` (default 11), `scf.Nrank.Rec.Coulomb` (default 900).
  These are the grid range, number of points, and rank of a feature under development that discretizes a reciprocal-space Coulomb-type kernel (the regularized \\( 1/(\lvert\mathbf q\_1 + \mathbf q\_2\rvert^2 + c) \\)) on a uniform grid
  and approximates it in the separable form \\( \sum\_i \lambda\_i A\_i(\mathbf q\_1) B\_i(\mathbf q\_2) \\) by singular value decomposition.
  That processing is not included in the build and is not called from anywhere.
- The source of Toyoda's exact-exchange module (which uses the electron-repulsion-integral method below to build the Fock matrix and the exchange energy) is also included,
  but it is not included in the build, and the switch enabling it is never set anywhere. For reference, its default settings are:
  for isolated systems, the bare Coulomb kernel without truncation; for periodic systems, the interaction truncated at 10 Å; no screening (\\( \omega = 0 \\)).
  It was designed to be unusable with solvers other than `Cluster` and `Band`.
- There is no exact-exchange contribution to forces or stress.
- No processing corresponding to `NLX` is found. From the name it can be guessed that nonlocal exchange (such as the method of PRA 2011 below) was intended, but this is [Open].

[Design] The reimplementation accepts `HSE` and `NLX` but stops as "not supported" (it does not silently compute without exchange-correlation).
If exact exchange is implemented in the future, the methods of the public documents below are the starting point, and the specification is written as a separate chapter.

### Matrix elements of Fock exchange [Public]

The following writes the general definitions in the notation of this note (not an implementation specification, but an organization for future extensions).

For real basis functions \\( \phi\_a \\) and the density matrix \\( \rho^\sigma\_{ab} \\) of spin \\( \sigma \\) ([Total energy and forces](energy-force.md)),
define the electron repulsion integral as

$$
(ac\vert bd) = \int d^3r\\, d^3r'\\, \frac{\phi\_a(\mathbf r)\phi\_c(\mathbf r)\\, \phi\_b(\mathbf r')\phi\_d(\mathbf r')}{\lVert\mathbf r - \mathbf r'\rVert}
$$

Then the Hartree–Fock exchange energy and the Fock exchange matrix are

$$
E\_x^{\mathrm{HF}} = -\frac12\sum\_\sigma\sum\_{abcd} \rho^\sigma\_{ab}\\, \rho^\sigma\_{cd}\\, (ac\vert bd),\qquad
K^\sigma\_{ab} = \frac{\partial E\_x^{\mathrm{HF}}}{\partial\rho^\sigma\_{ab}} = -\sum\_{cd} (ac\vert bd)\\, \rho^\sigma\_{cd}
$$

(for real symmetric \\( \rho \\); in periodic systems, the basis-function indices carry lattice vectors, and a sum over \\( \mathbf k \\) points appears).
For a hybrid functional with mixing fraction \\( \alpha \\), \\( \alpha K^\sigma \\) is added to the Kohn–Sham matrix, and the exchange of the semilocal functional is scaled by \\( (1-\alpha) \\).
In the HSE type (Heyd–Scuseria–Ernzerhof, JCP **118**, 8207 (2003)), only the short-range part, with \\( 1/r \\) replaced by \\( \mathrm{erfc}(\omega r)/r \\), is treated by Fock exchange.
Since \\( K \\) depends on \\( \rho \\), it is rebuilt in every SCF iteration.

### Method for computing electron repulsion integrals [Public]

The method of Toyoda–Ozaki (JCP **130**, 124114 (2009), technical note [exx/JCP2009](https://www.openmx-square.org/exx/JCP2009/index.html)) is rederived
with the [Fourier transform convention](../../../conventions/notation.md) of this note.

1. Regard the product of two basis functions as a single function centered at some point \\( \mathbf c \\), \\( F(\mathbf r - \mathbf c) = \phi\_a(\mathbf r - \mathbf R\_a)\phi\_c(\mathbf r - \mathbf R\_c) \\),
   and expand it in spherical harmonics: \\( F(\mathbf r) = \sum\_{lm} P\_{lm}(r)\\, Y\_{lm}(\hat{\mathbf r}) \\).
   For the expansion of a basis function centered away from \\( \mathbf c \\), Löwdin's α function (built with a spherical Bessel transform) is used.
   \\( G(\mathbf r - \mathbf c') \\) is built in the same way.
2. From the Fourier transform \\( 4\pi/k^2 \\) of \\( 1/r \\) and Parseval's identity, with \\( \mathbf R = \mathbf c' - \mathbf c \\),

$$
(F\vert G) = \frac{1}{2\pi^2}\int d^3k\\, \frac{1}{k^2}\\, \tilde F^{\ast}(\mathbf k)\\, \tilde G(\mathbf k)\\, e^{-i\mathbf k\cdot\mathbf R}
$$

3. Substituting \\( \tilde F(\mathbf k) = 4\pi\sum\_{lm}(-i)^l\\, \tilde P\_{lm}(k)\\, Y\_{lm}(\hat{\mathbf k}) \\) (\\( \tilde P\_{lm} \\) is the spherical Bessel transform) and
   the plane-wave expansion of \\( e^{-i\mathbf k\cdot\mathbf R} \\), and performing the angular integral with Gaunt coefficients \\( G(lm;\\, l'm',\\, \lambda\mu) \\)
   (same definition as in [Two-center integrals](two-center-overlap-kinetic.md)), gives

$$
(F\vert G) = 32\pi \sum\_{lm}\sum\_{l'm'}\sum\_{\lambda\mu} i^{\\,l - l' - \lambda}\\, G(lm;\\, l'm',\\, \lambda\mu)\\, Y^{\ast}\_{\lambda\mu}(\hat{\mathbf R})
\int\_0^\infty dk\\, \tilde P^{\ast}\_{lm}(k)\\, \tilde Q\_{l'm'}(k)\\, j\_\lambda(kR)
$$

   (\\( 32\pi = (4\pi)^3/(2\pi^2) \\), \\( Q \\) is the radial function of \\( G \\)). This has the same structure as the overlap formula, differing only by the kernel \\( k^{-2} \\).
4. For the screened kernel \\( \mathrm{erfc}(\omega r)/r \\), replace \\( 4\pi/k^2 \\) by \\( (4\pi/k^2)(1 - e^{-k^2/4\omega^2}) \\).
5. A fast method (Toyoda–Ozaki, CPC **181**, 277 (2010), technical note [exx/CPC2009](https://www.openmx-square.org/exx/CPC2009/index.html)) can be used for the spherical Bessel transform.
   JCP 2009 reports that the angular momentum of the product functions must be taken up to \\( l \lesssim 15 \\), that about 1024 radial grid points are used, that the final \\( k \\) integral
   is done by Gauss–Laguerre quadrature (about 60 points), and that derivatives with respect to atomic positions (forces) are computed with reduced cost using sum rules.

> **Check**: Formula 3 above was integrated numerically for s–s and p–s Gaussian product densities (bare kernel and erfc-screened kernel),
> and agreement with the closed-form expressions to within \\( 10^{-6} \\) relative is confirmed in `tools/check_dispersion.py`.

> **Note on the technical note JCP2009**: The final formula of that note appears to differ from this formula in the form of the phase factor, but
> that note differs from this note in the sign of the Fourier transform (\\( e^{+i\mathbf k\cdot\mathbf r} \\) type) and in where complex conjugates are placed, so the two cannot be compared directly.
> The reimplementation uses the formula in the conventions of this note (confirmed by the check above).

### Exchange functional based on a model exchange hole [Public]

Toyoda–Ozaki (PRA **83**, 032515 (2011), technical note [exx/PRA2011](https://www.openmx-square.org/exx/PRA2011/index.html))
propose an exchange functional in which the exchange hole is split into short- and long-range parts by \\( \mathrm{erfc}(\mu s) \\); the short-range part is treated exactly, and the long-range part is replaced by a model hole
built from a hydrogen-like 1s orbital (Becke–Roussel type). The two parameters of the model hole are determined at each point so that the zeroth moment of the hole is reproduced exactly
and the first moment as closely as possible; the potential is obtained analytically by the chain rule.
The calculations in that paper were done with the authors' finite-element code for atoms; no implementation or input in OpenMX is given in the public documents.
The examined version of OpenMX has no corresponding processing either (the relation to `NLX` is [Open]).

## Edge cases [Convention] [Design]

- [Convention] With `scf.dftD off` (default), \\( U\_{\mathrm{vdW}} = 0 \\), and there is no contribution to forces or stress.
- [Convention] For an isolated system with a single atom, \\( U\_{\mathrm{vdW}} = 0 \\) (when periodic images are excluded).
- [Convention] Because the energy is discontinuous at the cutoff, energy conservation in MD is slightly violated. With the default cutoff, the size of the discontinuity is
  \\( C\_6/r\_{\mathrm{cut}}^6 \\) per pair, i.e. of order \\( 10^{-11} \\) to \\( 10^{-10} \\) Hartree.
- [Convention] When two atoms come closer than 0.32 bohr, the pair is silently dropped from the sum.
- [Design] The reimplementation has the same discontinuity (for comparison), and a smooth cutoff is a separate option.

## Design [Design]

Requirements for a modular reimplementation.

- The **dispersion-correction module** is a pure function that takes as input the atomic numbers, atomic positions, lattice vectors (the box for isolated systems), periodicity flags, direction mask,
  and parameter set (version, damping function, \\( s\_6, s\_8, s\_{r,6}, a\_1, a\_2, d\_{\mathrm{D}} \\), cutoffs), and returns
  \\( U\_{\mathrm{vdW}} \\), per-atom gradients, and \\( \Sigma \\) (3×3). It does not depend on electronic-structure data.
- The part of the parameter defaults determined by `scf.XcType` (the D3 PBE values) is resolved in the input-interpretation layer before being passed to the module.
- The coefficient tables (G06 Table 1 for D2, reference data for D3) are read as external data, and their source, version, and units are recorded.
- Enumerating periodic images with the D3 interplanar-spacing rule is correct, but for comparison with the existing D2 implementation,
  the rule based on \\( \lVert\mathbf a\_k\rVert \\) is provided as a compatibility option.
- The three-body term (ATM term) is disabled by default for compatibility; if added, it is a separate option.
- For large systems, the sum over atom pairs is made \\( O(N) \\) with a neighbor list (cell partitioning). The existing implementation scans all atom pairs, \\( O(N^2) \\).
- `HSE` and `NLX` stop as not supported. If exact exchange is added, a separate layer builds the Fock matrix as atom-pair blocks of the [H(R), S(R) sparse block format](../design/sparse-block-format.md)
  and adds it to the SCF Hamiltonian assembly.

## Tests [Design]

### Analytic checks and invariants

| Kind | Content | Tolerance (guideline) |
|---|---|---|
| Analytic | D2 for 2 atoms (isolated, `IntDirection 0 0 0`): \\( U = -s\_6\sqrt{C\_6^IC\_6^J} f(d)/d^6 \\) | \\( 10^{-14} \\) (relative) |
| Analytic | D3 BJ for 2 atoms: finite even near \\( d \to 0 \\), approaching the \\( d = 0 \\) limit \\( -s\_6C\_6/R^6 - s\_8C\_8/R^8 \\) | \\( 10^{-12} \\) (relative) |
| Analytic | Unit conversion: \\( C\_6 \\), \\( R\_0 \\) of H, C, Si in G06 vs. the Ry bohr⁶ and bohr values of the existing implementation | relative \\( 2\times10^{-3} \\) |
| Finite difference | Gradient: central difference displacing an atom by \\( \pm10^{-4} \\) bohr | \\( 10^{-8} \\) Hartree/bohr |
| Finite difference | Stress: central difference with strain \\( \pm10^{-7} \\) (small strain to avoid pairs crossing the cutoff) | \\( 10^{-8} \\) Hartree |
| Invariant | All atoms periodic: gradients sum to 0 | rounding error |
| Invariant | \\( U \\) invariant under translation of all atoms and re-choice of the unit cell (a different primitive cell of the same crystal) (D3 interplanar-spacing rule, sufficiently large cutoff) | \\( 10^{-10} \\) (excluding differences in cutoff handling) |
| Invariant | \\( U \\) invariant under permutation of atoms; gradients permuted identically | rounding error |
| Invariant | \\( \Sigma \\) symmetric | rounding error |
| Invariant | Coordination number: for two isolated atoms, \\( \mathrm{CN} \\) is \\( \tfrac12 \\) when \\( d = R\_{\mathrm{cov}}^{IJ} \\) | \\( 10^{-14} \\) |
| Invariant | \\( C\_6 \\) interpolation: when \\( \mathrm{CN} \\) coincides with a pair of reference coordination numbers and is sufficiently far from the other references, it equals that reference value; far away, it approaches the nearest reference value | \\( 10^{-8} \\) |
| Analytic | Electron-repulsion-integral formula (for future use): closed form for Gaussian product densities (bare, erfc-screened) | relative \\( 10^{-6} \\) |

`tools/check_dispersion.py` implements the unit conversion, finite-difference checks of D2/D3-type gradients and stress, a comparison of the periodic-image range rules,
the far-field limit of \\( C\_6 \\), and the check of the electron-repulsion-integral formula.

### Comparison with fixtures

Compare with the output of the existing implementation (`UvdW`, forces, stress).

| Target | Setting | Quantities compared | Tolerance (guideline) |
|---|---|---|---|
| Benzene dimer (isolated) | D3, `zero` and `bj`, `DFTD.IntDirection 0 0 0`, several interplanar distances (same structure as the distributed example) | `UvdW`, vdW contribution to forces | \\( 10^{-10} \\) Hartree, \\( 10^{-9} \\) Hartree/bohr |
| Same as above | D2 | Same as above | \\( 10^{-10} \\) Hartree (when the same rounded table is used) |
| Graphite or h-BN layered crystal (`Band`) | D2, D3 (both dampings), default cutoffs, including oblique lattices | `UvdW`, forces, stress | Same as above. For D2, when the existing implementation's periodic-image rule is reproduced via the compatibility option |
| Molecule + surface model | Mixed flags | `UvdW` only | \\( 10^{-10} \\) Hartree |

- The vdW contributions to forces and stress are extracted as the difference between `scf.dftD` on and off (running SCF only once with the same density matrix, or from the same converged density).
- Since `UvdW` does not depend on the electronic state, the SCF need not be converged, and the tolerance can be made tight.

## [Open]

- In DFT-D3 of the distributed OpenMX 4.0 itself, the handling when the denominator of the Gaussian weights for \\( C\_6 \\) is small (set to 0, or subtract the maximum as in the examined version).
- For the defaults of `DFTD3.damp` and `DFTD.Unit`, whether the reimplementation's default should follow the manual or the existing implementation (this note follows the existing implementation).
- What value `NLX` was intended to be. Its relation to the method of PRA 2011.
- The final meaning and defaults of the four keywords read for `HSE` (separable approximation of the reciprocal-space kernel).
- Whether to provide a compatibility option that reproduces forces and stress with mixed periodicity flags in the same (incomplete) form as the existing implementation.
- The actual values of the fixtures (benzene dimer, layered crystal). They will be created once an environment for running the existing implementation is available.
- How to extract the vdW contributions to forces and stress individually from the output files (including energy-decomposition output).
