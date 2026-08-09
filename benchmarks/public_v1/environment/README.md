# Frozen comparator environments

This directory separates three claims that must not be conflated:

1. **Environment available**: the package imports and its dependency checker passes.
2. **Fold executable**: the adapter completes a truth-stripped, fold-train-only run and
   writes a validated prediction row for every query cell.
3. **Benchmark evidence**: every locked fold completes and the release-level donor,
   study, batch, and QC analyses pass their predeclared checks.

Only the third state supports comparative benchmark reporting. It does not by itself
establish performance on a new tissue, platform, or clinical population.

`freeze_benchmark_environments.py` writes and hashes:

- `python-celltypist.freeze.txt`: the isolated CellTypePilot/CellTypist Python packages;
- `r-singler.packages.csv` and `r-session-info.txt`: the isolated SingleR R library;
- `popv-image.json`: the complete Docker image inspection record;
- `python-popv.freeze.txt`: packages inside the accepted popV image;
- `installation_attempts.jsonl`: accepted and rejected setup paths, including failures;
- `environment_manifest.json`: UTC freeze time and SHA-256 for every artifact above.

The accepted popV image uses the digest-pinned Python 3.12 base declared in
`popv.Dockerfile`, CPU-only PyTorch, and popV 0.6.1. The benchmark adapter always uses
`prediction_mode="retrain"` and the current fold's training cells. Its agreement score is
not a calibrated probability. SingleR likewise trains only on the current fold; its
reported confidence is an affine transformation of the maximum similarity score.

Azimuth references are frozen separately because their training-study provenance is an
outcome-relevant part of the method. A reference that includes an evaluation study is
ineligible for the primary independent track even when it runs successfully.
