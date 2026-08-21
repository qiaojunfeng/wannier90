# Parallel Transport and Manifold-Remixed Wannier Functions

`wannier90` can construct a smooth Bloch gauge across the
Monkhorst&#151;Pack $k$-mesh by *parallel transport*, and use it to split a
combined valence+conduction Wannierisation into separate, isolated band
manifolds &#151; the *manifold-remixed Wannier functions* (MRWF).

Both are post-processing steps that act on an existing Wannierisation. They
can be run in the same job (after the minimisation), or on a saved
checkpoint via the `restart` keyword.

## Parallel transport

Setting `parallel_transport = .true.` replaces the current gauge with a
smooth one, obtained by parallel-transporting a reference frame along each
of the three $k$-mesh directions and removing the residual obstructions
dimension by dimension. The construction follows the homotopy algorithm of
Gontier, Levitt and Siraj-dine [@gontier-jmp19]. It requires a full
Monkhorst&#151;Pack mesh whose finite-difference stencil contains the six
Cartesian neighbours ($\pm k_x, \pm k_y, \pm k_z$).

The smoothed gauge is written to a new checkpoint (`seedname.chk`) and to
`seedname.pt.amn`. The `wout` file reports the neighbour-overlap smoothness
error before and after the transport.

- `parallel_transport`: run the parallel-transport smoothing.
- `parallel_transport_use_gauge`: seed the transport from the current gauge
  instead of the identity.
- `parallel_transport_log_interp`: use logarithmic interpolation for the
  edge/surface obstructions.

To run it as a standalone step on a converged calculation:

```vi title="Input file"
restart = parallel_transport
parallel_transport = .true.
```

## Manifold-remixed Wannier functions (MRWF)

Setting `mrwf = .true.` splits a Wannierisation that spans several isolated
energy manifolds (typically valence + conduction) into one smooth,
maximally-localisable problem per manifold. For each $k$-point the
Wannier-gauge Hamiltonian
$H^{\mathrm{W}}(\mathbf{k}) = U(\mathbf{k})^\dagger\,
\mathrm{diag}[\epsilon_{n\mathbf{k}}]\, U(\mathbf{k})$ is diagonalised; its
eigenvectors are grouped into the requested manifolds, each manifold is
parallel-transported to a smooth gauge, and the result is written as a
self-contained Wannier problem in a subdirectory `mrwf_g<i>/` (containing
`seedname.mmn`, `seedname.eig`, `seedname.amn`, `seedname.win`, and the block
gauge `seedname_split.amn` that maps the original bands to the manifold).

The manifolds are specified either by a simple valence/conduction cut:

```vi title="Input file"
mrwf = .true.
mrwf_num_val = 4        ! bands 1..4 valence, the rest conduction
```

or explicitly by a block, one contiguous band range per manifold (the
number of manifolds is counted automatically):

```vi title="Input file"
mrwf = .true.
begin mrwf_manifolds
  1  4
  5  8
end mrwf_manifolds
```

Exactly one of `mrwf_num_val` or the `mrwf_manifolds` block must be given,
and the manifolds must tile the Wannier bands $1\ldots$`num_wann`
contiguously.

- `mrwf`: enable the manifold split.
- `mrwf_num_val`: two-way valence/conduction split at this band index.
- `mrwf_manifolds`: block of `first_band last_band` rows, one per manifold.
- `mrwf_run_maxloc`: run a full maximal localisation on each manifold after
  the transport.
- `mrwf_write_unk`: write per-manifold `UNK` files for plotting the Wannier
  functions.

With `mrwf_run_maxloc = .true.` each manifold is passed through the standard
maximal-localisation minimiser after parallel transport, and the localised
gauge and its spread are reported and written out. With
`mrwf_write_unk = .true.` the input `UNK` files are rotated into each
manifold and written to the corresponding `mrwf_g<i>/` directory, ready for
Wannier-function plotting.

To run MRWF as a standalone step on a converged calculation:

```vi title="Input file"
restart = mrwf
mrwf = .true.
mrwf_num_val = 4
```
