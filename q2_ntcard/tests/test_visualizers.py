import qiime2
import pandas as pd

from q2_ntcard._visualizers import visualize_kmer_histograms


def test_visualizer_writes_self_contained_html(tmp_path):
    data = pd.DataFrame({
        "kmer_length": [21, 21], "frequency": [1, 2],
        "estimated_distinct_kmers": [40.0, 10.0],
    }, index=pd.Index(["k21-f1", "k21-f2"], name="id"))
    visualize_kmer_histograms(tmp_path, qiime2.Metadata(data), title="Test")
    observed = (tmp_path / "index.html").read_text()
    assert "Test" in observed
    assert "estimated_distinct_kmers" in observed
