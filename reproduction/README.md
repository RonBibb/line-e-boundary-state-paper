# Reproduction inputs

The checkout bundles the frozen CBE1, CBE2i, CBE3, CBE3T and CBE3b records and sources under `inputs/`. The CBE3b inputs include the three original binary traces at resolutions 1024, 2048 and 4096. The CBE3c, Deep40 and LimFit directories preserve the result records used to assess the limits of the historical claims.

From this directory, verify the input bytes with `sha256sum -c SHA256SUMS`. On macOS use `shasum -a 256 -c SHA256SUMS`.

Install the pinned Python packages from `requirements.txt` in a Linux virtual environment. See the [continuation replay instructions](../research/2026-11-25_snapshot_observability/extensions/2026-11-25_nonlinear_snapshot_obstructions/continued_investigation/REPRODUCE.md) and the [historical E1 verification commands](../REPRODUCIBILITY_MANIFEST.md).

The historical producers in `inputs/` write alongside their sources. Run them from a scratch copy of the inputs if you want to regenerate the frozen outputs.
