from scripts.build_public_label_maps import broad_class


def test_broad_class_avoids_substring_lineage_collisions() -> None:
    assert broad_class("club cell") == "epithelial"
    assert broad_class("goblet cell") == "epithelial"
    assert broad_class("intestinal tuft cell") == "epithelial"
    assert broad_class("B cell") == "b_cell"
    assert broad_class("T cell") == "t_cell"


def test_broad_class_preserves_abstention() -> None:
    assert broad_class("Unknown") == "Unknown"
