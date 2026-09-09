import html
import json
from pathlib import Path

import qiime2


def _build_spec(data, title):
    return {
        "$schema": "https://vega.github.io/schema/vega-lite/v6.json",
        "title": {
            "text": title,
            "subtitle": "Estimated distinct k-mers by multiplicity",
            "anchor": "start",
        },
        "data": {"values": data},
        "transform": [
            {"filter": "datum.estimated_distinct_kmers > 0"},
        ],
        "facet": {
            "field": "facet_label",
            "type": "ordinal",
            "title": None,
            "sort": {"field": "kmer_length", "op": "min"},
        },
        "columns": 2,
        "spec": {
            "width": 360,
            "height": 240,
            "mark": {
                "type": "line",
                "point": {"filled": True, "size": 35},
                "strokeWidth": 2,
                "color": "#3659a2",
            },
            "encoding": {
                "x": {
                    "field": "frequency",
                    "type": "quantitative",
                    "scale": {"type": "log"},
                    "axis": {"title": "k-mer frequency", "grid": True},
                },
                "y": {
                    "field": "estimated_distinct_kmers",
                    "type": "quantitative",
                    "scale": {"type": "log"},
                    "axis": {"title": "estimated distinct k-mers", "grid": True},
                },
                "order": {"field": "frequency", "type": "quantitative"},
                "tooltip": [
                    {
                        "field": "kmer_length",
                        "type": "ordinal",
                        "title": "k-mer length",
                    },
                    {
                        "field": "F0",
                        "type": "quantitative",
                        "title": "F0 (estimated distinct k-mers)",
                        "format": ",.3~g",
                    },
                    {
                        "field": "frequency",
                        "type": "quantitative",
                        "title": "frequency",
                        "format": ",",
                    },
                    {
                        "field": "estimated_distinct_kmers",
                        "type": "quantitative",
                        "title": "estimated distinct k-mers",
                        "format": ",.3~g",
                    },
                ],
            },
        },
        "resolve": {"scale": {"x": "independent", "y": "independent"}},
        "config": {
            "view": {"stroke": "#d8dee4"},
            "facet": {"spacing": 24},
            "axis": {"labelFontSize": 11, "titleFontSize": 12},
            "header": {"labelFontSize": 14, "labelFontWeight": "bold"},
        },
    }


def _visualize(
    output_dir: str,
    histogram: qiime2.Metadata,
    summary: qiime2.Metadata,
    title: str = "ntCard k-mer coverage histograms",
) -> None:
    frame = histogram.to_dataframe().reset_index(drop=True)
    required = {"kmer_length", "frequency", "estimated_distinct_kmers"}
    if not required.issubset(frame.columns):
        raise ValueError(
            "Histogram metadata must contain kmer_length, frequency, and "
            "estimated_distinct_kmers columns."
        )
    summary_frame = summary.to_dataframe().reset_index(drop=True)
    summary_required = {"kmer_length", "estimated_distinct_kmers"}
    if not summary_required.issubset(summary_frame.columns):
        raise ValueError(
            "Summary metadata must contain kmer_length and "
            "estimated_distinct_kmers columns."
        )
    f0_by_k = {
        int(row.kmer_length): float(row.estimated_distinct_kmers)
        for row in summary_frame.itertuples()
    }
    missing = set(frame["kmer_length"].astype(int)) - set(f0_by_k)
    if missing:
        raise ValueError(
            "Summary metadata is missing k-mer size(s): "
            + ", ".join(map(str, sorted(missing)))
        )
    frame["F0"] = frame["kmer_length"].map(lambda value: f0_by_k[int(value)])
    frame["facet_label"] = frame.apply(
        lambda row: f"k = {int(row.kmer_length)} · F0 = {row.F0:,.0f}", axis=1
    )
    columns = [
        "kmer_length",
        "frequency",
        "estimated_distinct_kmers",
        "F0",
        "facet_label",
    ]
    data = frame[columns].to_dict(orient="records")
    spec = _build_spec(data, title)
    template = Path(__file__).parent / "assets" / "index.html"
    document = template.read_text().replace("__TITLE__", html.escape(title))
    embedded_spec = json.dumps(spec).replace("<", "\\u003c")
    document = document.replace("__SPEC__", embedded_spec)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "index.html").write_text(document)
    (destination / "vega-lite-spec.json").write_text(json.dumps(spec, indent=2) + "\n")
