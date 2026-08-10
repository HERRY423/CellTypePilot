# CellTypePilot release procedure

CellTypePilot has two coordinated release artifacts with the same semantic version:

1. `celltypepilot` on PyPI is the deterministic Python backend.
2. `celltypepilot-plugin-<version>.zip` on GitHub Releases is the complete Codex/Claude plugin.

The GitHub Release is created before the irreversible PyPI upload. A failed verification job
must not create either release.

## One-time Trusted Publisher setup

Create a `pypi` environment in the GitHub repository, then configure a PyPI Trusted Publisher
with exactly these values:

| Field | Value |
|---|---|
| PyPI project | `celltypepilot` |
| GitHub owner | `HERRY423` |
| Repository | `CellTypePilot` |
| Workflow | `release.yml` |
| Environment | `pypi` |

For the first upload, configure a pending publisher in PyPI account settings. Do not add a PyPI
API token to GitHub. The publish job receives only `id-token: write`; build and GitHub Release
jobs do not receive that permission.

## Release checklist

1. Confirm the intended commit is the current `main` HEAD and its CI is green.
2. Confirm `pyproject.toml`, `celltypepilot.__version__`, both plugin manifests, CHANGELOG, and
   the proposed `v<version>` tag agree.
3. Confirm the CHANGELOG heading contains the actual release date and retains the scientific
   validation boundary.
4. Build locally with `python -m build` and run `python -m twine check dist/*`.
5. Build the plugin bundle with
   `python scripts/build_plugin_bundle.py --tag v<version> --output-dir release-assets`.
6. Create the annotated tag at current `main`, then push only that tag. The workflow rejects a
   tag that does not point to current `origin/main`.
7. Verify the `Verify & Build Release Artifacts` and `Publish GitHub Plugin Release` jobs before
   allowing the protected `pypi` environment deployment.
8. After publication, verify the GitHub assets and hashes, install the wheel from PyPI in a clean
   environment, and confirm the PyPI/GitHub versions match.

## Non-negotiable failure behavior

- Never move or overwrite a published tag.
- Never use `skip-existing` for the production PyPI release.
- Never upload a rebuilt file under an already published PyPI version.
- Never describe a technical release as biological superiority or a completed validation release.
- If PyPI publishing fails after the GitHub Release exists, keep the evidence, correct the
  publisher configuration, and rerun only through the audited workflow; do not publish manually.
