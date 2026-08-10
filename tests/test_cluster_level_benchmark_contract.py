import pandas as pd

from celltypepilot.benchmark import build_cluster_level_track, evaluate_holdout_predictions


def test_cluster_track_aggregates_every_method_to_same_locked_unit():
    cells = ["c1", "c2", "c3", "c4", "c5", "c6"]
    truth = pd.Series(["A", "A", "B", "B", "B", "A"], index=cells)
    assignments = pd.DataFrame(
        {
            "cell_id": cells,
            "fold_id": ["study=s1"] * 6,
            "role": ["test"] * 6,
            "held_out_study": ["s1"] * 6,
            "held_out_donor": ["s1::d1"] * 6,
        }
    )
    clusters = pd.Series(["0", "0", "0", "1", "1", "1"], index=cells)
    predictions = pd.DataFrame(
        [
            {"cell_id": cell, "fold_id": "study=s1", "method": method, "predicted_label": label}
            for method, labels in {
                "celltypepilot": ["A", "A", "A", "B", "B", "B"],
                "celltypist": ["A", "A", "B", "B", "B", "A"],
            }.items()
            for cell, label in zip(cells, labels, strict=True)
        ]
    )

    cluster_truth, cluster_assignments, cluster_predictions, diagnostics = (
        build_cluster_level_track(truth, assignments, predictions, clusters)
    )
    assert len(cluster_truth) == 2
    assert len(cluster_assignments) == 2
    assert set(cluster_predictions["method"]) == {"celltypepilot", "celltypist"}
    assert set(diagnostics["truth_purity"]) == {2 / 3}

    results, _ = evaluate_holdout_predictions(
        cluster_truth,
        cluster_assignments,
        cluster_predictions,
        expected_methods=("celltypepilot", "celltypist"),
        bootstrap_ci=False,
    )
    by_method = results.set_index("method")
    assert by_method.loc["celltypepilot", "accuracy"] == 1.0
    assert by_method.loc["celltypist", "accuracy"] == 1.0


def test_cluster_truth_tie_is_excluded_and_reported_not_silently_broken():
    truth = pd.Series(["A", "B"], index=["c1", "c2"])
    assignments = pd.DataFrame(
        {
            "cell_id": ["c1", "c2"],
            "fold_id": ["study=s1", "study=s1"],
            "role": ["test", "test"],
            "held_out_study": ["s1", "s1"],
            "held_out_donor": ["s1::d1", "s1::d1"],
        }
    )
    predictions = pd.DataFrame(
        {
            "cell_id": ["c1", "c2"],
            "fold_id": ["study=s1", "study=s1"],
            "method": ["celltypepilot", "celltypepilot"],
            "predicted_label": ["A", "A"],
        }
    )
    clusters = pd.Series(["0", "0"], index=["c1", "c2"])

    cluster_truth, assignments_out, predictions_out, diagnostics = build_cluster_level_track(
        truth, assignments, predictions, clusters
    )
    assert cluster_truth.empty
    assert assignments_out.empty
    assert predictions_out.empty
    assert diagnostics.iloc[0]["status"] == "ambiguous_truth_tie"
