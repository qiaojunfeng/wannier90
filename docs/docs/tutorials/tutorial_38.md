# 38: Silicon &#151; Manifold-remixed Wannier functions

- Outline: *Split a combined valence + conduction Wannierisation of Si into
    two isolated manifolds using parallel transport, and maximally localise
    each manifold separately (manifold-remixed Wannier functions, MRWF).*

- Generation Details: *As tutorial 3 (from `pwscf`, norm-conserving
    pseudopotentials, 4$\times$4$\times$4 k-point grid, atom-centred sp$^3$
    starting guess).*

- Directory: `tutorials/tutorial38/`

- Input Files

    - `silicon.win` *The master input file*

    - `silicon.mmn` *The overlap matrices*

    - `silicon.amn` *The projections onto the trial orbitals*

    - `silicon.eig` *The Bloch eigenvalues at each k-point*

Tutorial 3 disentangles eight MLWFs spanning the four valence and four
low-lying conduction bands of silicon as a single, entangled set. Often one
would instead like *separate* sets of Wannier functions for the valence and
conduction manifolds &#151; for example to build manifold-specific models.
Diagonalising the Wannier Hamiltonian at each k-point and simply splitting
the eigenvectors gives a discontinuous (random) gauge within each manifold,
which cannot be maximally localised directly. The MRWF procedure fixes this
by *parallel-transporting* each manifold to a smooth gauge first.

## Running the calculation

The input file is the tutorial-3 silicon calculation with the MRWF keywords
appended:

```vi title="silicon.win"
mrwf            = .true.
mrwf_num_val    = 4
mrwf_run_maxloc = .true.
```

`mrwf_num_val = 4` requests a two-way split: the four lowest Wannier bands
form the valence manifold and the remaining four the conduction manifold.
Run `wannier90` as usual:

```bash
wannier90.x silicon
```

After the ordinary disentanglement and wannierisation of the combined
eight-band set, `wannier90` builds the Wannier Hamiltonian
$H^{\mathrm{W}}(\mathbf{k})$, diagonalises it, splits the eigenvectors into
the two manifolds, parallel-transports each to a smooth gauge, and &#151;
because `mrwf_run_maxloc = .true.` &#151; maximally localises each manifold.
The final spreads are reported in `silicon.wout`; the valence manifold is
noticeably more localised than the conduction one.

## Output

Each manifold is written as a self-contained Wannier problem in its own
subdirectory, `mrwf_g1/` (valence) and `mrwf_g2/` (conduction):

- `silicon.mmn`, `silicon.eig`, `silicon.win` *The overlaps, eigenvalues and
    a ready-to-run input file for that manifold*

- `silicon.amn` *The smooth (parallel-transported, then localised) gauge*

- `silicon_split.amn` *The block gauge mapping the original Bloch bands to
    the manifold's Wannier functions*

To plot the Wannier functions of a manifold, add `mrwf_write_unk = .true.`
(and provide the `UNK` files, as in the Wannier-function plotting tutorials);
the rotated `UNK` files are written into each `mrwf_g<i>/` directory.

## Further ideas

- Replace `mrwf_num_val` with an explicit `begin mrwf_manifolds ... end`
    block to split into more than two manifolds.

- Compare the per-manifold spreads with those of the combined eight-band
    Wannierisation of tutorial 3.
