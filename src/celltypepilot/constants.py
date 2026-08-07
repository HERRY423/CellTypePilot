"""Shared constants for CellTypePilot."""

from pathlib import Path

# Package root
PKG_ROOT = Path(__file__).resolve().parent

# Built-in marker atlas path
ATLAS_PATH = PKG_ROOT / "data" / "marker_atlas.json"

# Confidence tiers
CONFIDENCE_HIGH = "high"
CONFIDENCE_MEDIUM = "medium"
CONFIDENCE_LOW = "low"
CONFIDENCE_REVIEW = "needs_review"

# Consensus tiers (Phase 2)
TIER0 = 0  # marker score + LLM agree → adopt
TIER1 = 1  # disagreement → escalate to multi-model debate

# Minimum cells for a cluster to be annotated
MIN_CLUSTER_SIZE = 10

# Marker scoring thresholds
MARKER_PCT_THRESHOLD = 0.25  # marker expressed in ≥25% of cluster cells
MARKER_FC_THRESHOLD = 0.5   # log2 fold change ≥ 0.5
MARKER_SPECIFICITY_THRESHOLD = 0.3  # specificity score ≥ 0.3

# Critic thresholds
CRITIC_NEG_MARKER_PCT_THRESHOLD = 0.15  # negative marker expressed in >15% → flag
CRITIC_DOUBLET_COEXPR_THRESHOLD = 0.3   # two mutually exclusive lineages co-expressed >30% → flag
CRITIC_LOW_COVERAGE_THRESHOLD = 0.2     # <20% of expected markers detected → low confidence

# Default output filenames
OUTPUT_ANNOTATED = "data.annotated.h5ad"
OUTPUT_EVIDENCE = "evidence_table.csv"
OUTPUT_MANIFEST = "manifest.json"
OUTPUT_FIGURES_DIR = "figures"
OUTPUT_REPORT = "report_draft.html"

# Species
SPECIES_HUMAN = "human"
SPECIES_MOUSE = "mouse"

# Color-blind friendly palette (Wong palette)
CB_PALETTE = [
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # green
    "#CC79A7",  # pink
    "#56B4E9",  # light blue
    "#D55E00",  # vermillion
    "#F0E442",  # yellow
    "#000000",  # black
]
