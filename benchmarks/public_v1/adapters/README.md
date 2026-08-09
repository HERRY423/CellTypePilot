# Comparator adapters

These adapters receive a fold-train H5AD whose `obs.cell_type` is populated and a
truth-stripped test H5AD. They must emit exactly
`cell_id,predicted_label,confidence` for every test cell.

- `singler.json` targets the Bioconductor SingleR implementation. Its confidence is a
  bounded transformation of the maximum similarity score, not a calibrated probability.
  The isolated R library uses zellkonverter's native R HDF5 reader so it does not create
  an untracked Python environment.
- `popv.json` targets popV `retrain` mode on the fold reference. Its confidence is
  the fraction of executed experts agreeing with the consensus, not a probability. The
  Windows adapter launches the pinned `celltypepilot-popv:0.6.1` Linux image and mounts
  only the current fold workspace; the accepted image identity and package freeze live
  under `environment/`.
- CellTypist is implemented natively by `benchmark_runner.py` and is retrained per fold.
- Azimuth is intentionally not given a fake generic adapter. A valid Azimuth run requires
  a tissue-compatible reference containing `ref.Rds` and `idx.annoy`. Until a frozen
  reference has a training-study overlap audit, the release records Azimuth as unavailable.

Hosted or pretrained references may be useful as a separate external-reference track, but
they must not be described as fold-trained or study-independent without an overlap audit.
The adapter contracts remain unavailable until their dependency preflight and every locked
fold complete successfully; merely shipping these files does not validate either tool.
