import qiime2
import json
import pandas as pd
import pytest

from q2_ntcard._visualizers import _build_spec, _visualize


def test_visualizer_writes_self_contained_html(tmp_path):
    data = pd.DataFrame(
        {
            "kmer_length": [21, 21, 31, 31],
            "frequency": [1, 2, 1, 2],
            "estimated_distinct_kmers": [40.0, 10.0, 30.0, 8.0],
        },
        index=pd.Index(["k21-f1", "k21-f2", "k31-f1", "k31-f2"], name="id"),
    )
    summary = pd.DataFrame(
        {
            "kmer_length": [21, 31],
            "estimated_distinct_kmers": [50.0, 38.0],
        },
        index=pd.Index(["k21", "k31"], name="id"),
    )
    _visualize(tmp_path, qiime2.Metadata(data), qiime2.Metadata(summary), title="Test")
    observed = (tmp_path / "index.html").read_text()
    assert "Test" in observed
    assert "vegaEmbed" in observed
    spec = json.loads((tmp_path / "vega-lite-spec.json").read_text())
    assert spec["facet"]["field"] == "facet_label"
    assert spec["columns"] == 2
    assert {row["kmer_length"] for row in spec["data"]["values"]} == {21, 31}
    assert {row["F0"] for row in spec["data"]["values"]} == {50, 38}
    assert {row["facet_label"] for row in spec["data"]["values"]} == {
        "k = 21 · F0 = 50",
        "k = 31 · F0 = 38",
    }


def test_vega_lite_spec_uses_log_scales_and_tooltips():
    spec = _build_spec([], "Test")
    assert spec["$schema"].endswith("vega-lite/v6.json")
    assert spec["spec"]["encoding"]["x"]["scale"]["type"] == "log"
    assert spec["spec"]["encoding"]["y"]["scale"]["type"] == "log"
    assert len(spec["spec"]["encoding"]["tooltip"]) == 4


def test_visualizer_rejects_missing_f0(tmp_path):
    histogram = pd.DataFrame(
        {
            "kmer_length": [21],
            "frequency": [1],
            "estimated_distinct_kmers": [40.0],
        },
        index=pd.Index(["k21-f1"], name="id"),
    )
    summary = pd.DataFrame(
        {
            "kmer_length": [31],
            "estimated_distinct_kmers": [38.0],
        },
        index=pd.Index(["k31"], name="id"),
    )
    with pytest.raises(ValueError, match="missing k-mer size"):
        _visualize(tmp_path, qiime2.Metadata(histogram), qiime2.Metadata(summary))
