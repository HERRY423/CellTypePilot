"""Gene Alias & Synonym Resolution Module.

Maps gene symbol aliases/synonyms to canonical names present in the AnnData matrix.
Prevents false-negative marker dropouts caused by gene naming variations (e.g. CD303 <-> CLEC4C).
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

# Common HGNC gene symbol aliases for single-cell marker panels
# Format: Alias/Synonym -> Canonical HGNC Symbol
GENE_ALIAS_MAP: dict[str, str] = {
    # Immune / Dendritic cell markers
    "CD303": "CLEC4C",
    "BDCA2": "CLEC4C",
    "BDCA-2": "CLEC4C",
    "CD304": "NRP1",
    "BDCA4": "NRP1",
    "BDCA-4": "NRP1",
    "CD141": "THBD",
    "BDCA3": "THBD",
    "BDCA-1": "CD1C",
    "CD123": "IL3RA",
    "CD11C": "ITGAX",
    "CD11B": "ITGAM",
    "CD16": "FCGR3A",
    "CD32": "FCGR2A",
    "CD64": "FCGR1A",
    "CD3E": "CD3E",
    "CD3D": "CD3D",
    "CD3G": "CD3G",
    "CD4": "CD4",
    "CD8A": "CD8A",
    "CD8B": "CD8B",
    "CD19": "CD19",
    "CD20": "MS4A1",
    "CD22": "CD22",
    "CD25": "IL2RA",
    "CD127": "IL7R",
    "CD278": "ICOS",
    "CD279": "PDCD1",
    "CD274": "CD274",
    "CTLA4": "CTLA4",
    "CD56": "NCAM1",
    "CD161": "KLRB1",
    "CD314": "KLRK1",
    "CD335": "NCR1",
    "CD337": "NCR3",
    "CD94": "KRC1",
    "TIGIT": "TIGIT",
    "LAG3": "LAG3",
    "TIM3": "HAVCR2",
    "CD366": "HAVCR2",
    "CD152": "CTLA4",
    "FOXP3": "FOXP3",
    "ROSA26": "Gt(ROSA)26Sor",

    # Brain / Neural markers
    "GFAP": "GFAP",
    "Gfap": "Gfap",
    "SLCA1A3": "SLC1A3",
    "GLAST": "SLC1A3",
    "GLT1": "SLC1A2",
    "GLT-1": "SLC1A2",
    "NEUN": "RBFOX3",
    "NeuN": "Rbfox3",
    "IONA": "AIF1",
    "IBA1": "AIF1",
    "Iba1": "Aif1",
    "CX3CR1": "CX3CR1",
    "TMEM119": "TMEM119",
    "Tmem119": "Tmem119",
    "OLIG2": "OLIG2",
    "Olig2": "Olig2",
    "MBP": "MBP",
    "Mbp": "Mbp",

    # Mouse equivalent aliases (lowercase/titlecase)
    "Cd303": "Clec4c",
    "Bdca2": "Clec4c",
    "Cd304": "Nrp1",
    "Bdca4": "Nrp1",
    "Cd141": "Thbd",
    "Cd123": "Il3ra",
    "Cd11c": "Itgax",
    "Cd11b": "Itgam",
    "Cd16": "Fcgr3a",
    "Cd56": "Ncam1",
    "Iba1": "Aif1",
}


def build_var_alias_index(var_names: list[str] | set[str]) -> dict[str, str]:
    """Build a mapping from any gene name (or alias) to the exact string present in var_names.

    Args:
        var_names: List or set of gene names present in the AnnData matrix.

    Returns:
        dict mapping: input_query_name -> exact_var_name_in_matrix
    """
    var_set = set(var_names)
    index: dict[str, str] = {}

    # Direct identity mapping
    for v in var_set:
        index[v] = v
        index[v.upper()] = v
        index[v.capitalize()] = v

    # Alias mappings
    for alias, canonical in GENE_ALIAS_MAP.items():
        if canonical in var_set:
            index[alias] = canonical
            index[alias.upper()] = canonical
            index[alias.capitalize()] = canonical
        elif alias in var_set:
            index[canonical] = alias

    return index


def resolve_marker_list(markers: list[str], alias_index: dict[str, str]) -> tuple[list[str], list[str]]:
    """Resolve a list of expected marker gene names against the AnnData var_names alias index.

    Args:
        markers: List of target marker gene symbols.
        alias_index: Precomputed alias index from build_var_alias_index.

    Returns:
        tuple (resolved_present_in_matrix, missing_from_matrix)
    """
    present = []
    missing = []

    for m in markers:
        resolved = alias_index.get(m) or alias_index.get(m.upper()) or alias_index.get(m.capitalize())
        if resolved and resolved not in present:
            present.append(resolved)
        elif not resolved:
            missing.append(m)

    return present, missing
