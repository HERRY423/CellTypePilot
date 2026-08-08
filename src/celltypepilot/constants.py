"""Shared constants for CellTypePilot."""

from pathlib import Path

# Package root
PKG_ROOT = Path(__file__).resolve().parent

# Built-in marker atlas path
ATLAS_PATH = PKG_ROOT / "data" / "marker_atlas.json"
STATE_ATLAS_PATH = PKG_ROOT / "data" / "state_atlas.json"

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
MARKER_FC_THRESHOLD = 0.5  # log2 fold change ≥ 0.5
MARKER_SPECIFICITY_THRESHOLD = 0.3  # specificity score ≥ 0.3

MARKER_FDR_THRESHOLD = 0.05  # Benjamini-Hochberg adjusted p-value

# Critic thresholds
CRITIC_NEG_MARKER_PCT_THRESHOLD = (
    0.20  # negative marker expressed in >20% → flag (tolerates dropout noise)
)
CRITIC_DOUBLET_COEXPR_THRESHOLD = 0.3  # cross-lineage co-expression ≥30% → doublet candidate
CRITIC_DOUBLET_ACTIVE_COVERAGE = 0.4  # a signature counts as "active" at ≥40% coverage
CRITIC_DOUBLET_OVERLAP_JACCARD = 0.4  # marker Jaccard ≥0.4 → redundant signatures, not a doublet
CRITIC_LOW_COVERAGE_THRESHOLD = 0.2  # <20% of all expected markers expressed → low evidence

# Reference scoring thresholds
REF_MIN_SHARED_GENES = 100  # min shared genes for KNN/correlation
REF_CORR_MIN_GENES = 100  # min genes for correlation backend
REF_SCANVI_MIN_GENES = 200  # min genes for scANVI backend
REF_CELLTYPIST_MIN_GENES = 50  # min genes for CellTypist backend
REF_KNN_DEFAULT_K = 15  # default K for KNN label transfer
REF_KNN_MAX_K = 100  # max K for KNN label transfer

# Ensemble scoring thresholds
ENSEMBLE_AGREEMENT_THRESHOLD = 0.2  # score diff < 0.2 = agreement
ENSEMBLE_MARKER_HIGH = 0.6  # marker score ≥ this → marker-heavy weight
ENSEMBLE_MARKER_LOW = 0.3  # marker score ≤ this → reference-heavy weight
ENSEMBLE_REF_OVERRIDE = 0.5  # ref score ≥ this + marker low → ref override

# Default output filenames
OUTPUT_ANNOTATED = "data.annotated.h5ad"
OUTPUT_EVIDENCE = "evidence_table.csv"
OUTPUT_MANIFEST = "manifest.json"
OUTPUT_FIGURES_DIR = "figures"
OUTPUT_REPORT = "report_draft.html"

# Species
SPECIES_HUMAN = "human"
SPECIES_MOUSE = "mouse"
SPECIES_RAT = "rat"
SPECIES_ZEBRAFISH = "zebrafish"
SPECIES_CHICKEN = "chicken"
SPECIES_PIG = "pig"
SPECIES_COW = "cow"
SPECIES_MACAQUE = "macaque"
SPECIES_DOG = "dog"
SPECIES_UNKNOWN = "unknown"

# Ensembl gene ID prefix → species (longer prefixes first to avoid
# ambiguity, e.g. ENSGALG would otherwise match the ENSG prefix).
ENSEMBL_PREFIX_SPECIES = [
    ("ENSMUSG", SPECIES_MOUSE),
    ("ENSRNOG", SPECIES_RAT),
    ("ENSDARG", SPECIES_ZEBRAFISH),
    ("ENSGALG", SPECIES_CHICKEN),
    ("ENSSSCG", SPECIES_PIG),
    ("ENSBTAG", SPECIES_COW),
    ("ENSMMUG", SPECIES_MACAQUE),
    ("ENSCAFG", SPECIES_DOG),
    ("ENSG", SPECIES_HUMAN),
]

# Fraction of sampled genes that must match an Ensembl prefix / symbol
# convention before the species call is considered confident.
SPECIES_DOMINANCE_RATIO = 0.5
SPECIES_SYMBOL_RATIO = 2.0

# obs column names that may carry tissue context (matched case-insensitively)
TISSUE_COLUMN_SYNONYMS = [
    "tissue",
    "tissue_type",
    "tissue_origin",
    "tissue_source",
    "sample_tissue",
    "organ",
    "organ_system",
    "anatomy",
    "anatomical_site",
    "anatomy_site",
    "body_site",
    "body_part",
    "source",
    "sample_source",
    "sample_type",
    "location",
    "site",
    "origin",
]

# Fallback keywords scanned as substrings of obs column names
TISSUE_COLUMN_KEYWORDS = ["tissue", "organ", "anatom", "body_site"]

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
