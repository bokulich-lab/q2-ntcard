import html
import json
from pathlib import Path

import qiime2


def visualize_kmer_histograms(
    output_dir: str, histogram: qiime2.Metadata,
    title: str = "ntCard k-mer coverage histograms",
) -> None:
    frame = histogram.to_dataframe().reset_index(drop=True)
    required = {"kmer_length", "frequency", "estimated_distinct_kmers"}
    if not required.issubset(frame.columns):
        raise ValueError(
            "Histogram metadata must contain kmer_length, frequency, and "
            "estimated_distinct_kmers columns."
        )
    data = frame[list(required)].to_dict(orient="records")
    template = Path(__file__).parent / "assets" / "index.html"
    document = template.read_text().replace("__TITLE__", html.escape(title))
    document = document.replace("__DATA__", json.dumps(data))
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "index.html").write_text(document)
