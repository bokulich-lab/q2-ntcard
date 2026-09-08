from pathlib import Path
import re
import subprocess
import tempfile

import pandas as pd
import qiime2

from ._run_command import run_command
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt


def _parse_kmer_sizes(value):
    try:
        sizes = [int(item.strip()) for item in value.split(",")]
    except ValueError as error:
        raise ValueError("kmer_sizes must be a comma-separated list of integers.") \
            from error
    if not sizes or any(size < 1 for size in sizes):
        raise ValueError("All k-mer sizes must be positive integers.")
    if len(set(sizes)) != len(sizes):
        raise ValueError("kmer_sizes must not contain duplicates.")
    return sizes


def _read_histograms(paths, expected_sizes):
    records = []
    reported = {}
    for path in sorted(paths):
        rows = pd.read_csv(path, sep=r"\s+", comment="#", header=None)
        if rows.empty:
            continue
        if rows.shape[1] == 3:
            frame = rows.set_axis(
                ["kmer_length", "frequency", "estimated_distinct_kmers"],
                axis="columns",
            ).apply(pd.to_numeric, errors="raise")
        elif rows.shape[1] == 2:
            match = re.search(r"_k(\d+)\.hist$", Path(path).name)
            if match is None:
                raise RuntimeError(
                    f"Cannot determine k-mer length from ntCard output {path}."
                )
            kmer_length = int(match.group(1))
            metrics = rows[rows[0].isin(["F0", "F1"])]
            reported[kmer_length] = {
                str(row[0]): float(row[1]) for _, row in metrics.iterrows()
            }
            frame = rows[~rows[0].isin(["F0", "F1"])].copy()
            frame.columns = ["frequency", "estimated_distinct_kmers"]
            frame.insert(0, "kmer_length", kmer_length)
            frame = frame.apply(pd.to_numeric, errors="raise")
        else:
            raise RuntimeError(
                f"Unexpected ntCard histogram layout in {path}: "
                f"expected two or three columns, observed {rows.shape[1]}."
            )
        records.append(frame)
    if not records:
        raise RuntimeError("ntCard completed but did not produce a histogram file.")

    histogram = pd.concat(records, ignore_index=True)
    histogram["kmer_length"] = histogram["kmer_length"].astype(int)
    histogram["frequency"] = histogram["frequency"].astype(int)
    observed = set(histogram["kmer_length"])
    missing = set(expected_sizes) - observed
    if missing:
        raise RuntimeError(
            "ntCard output is missing requested k-mer size(s): "
            + ", ".join(map(str, sorted(missing)))
        )
    histogram.index = [
        f"k{row.kmer_length}-f{row.frequency}"
        for row in histogram.itertuples()
    ]
    histogram.index.name = "id"
    histogram.attrs["ntcard_summary"] = reported
    return histogram


def _summarize(histogram):
    rows = []
    reported = histogram.attrs.get("ntcard_summary", {})
    for kmer_length, group in histogram.groupby("kmer_length", sort=True):
        kmer_length = int(kmer_length)
        counts = group["estimated_distinct_kmers"]
        frequencies = group["frequency"]
        distinct = reported.get(kmer_length, {}).get("F0", counts.sum())
        total = reported.get(kmer_length, {}).get(
            "F1", (frequencies * counts).sum()
        )
        singletons = counts[frequencies == 1].sum()
        occupied = group.loc[counts > 0, "frequency"]
        rows.append({
            "id": f"k{kmer_length}",
            "kmer_length": int(kmer_length),
            "estimated_distinct_kmers": float(distinct),
            "estimated_total_kmers": float(total),
            "estimated_singleton_kmers": float(singletons),
            "maximum_observed_frequency": (
                int(occupied.max()) if not occupied.empty else 0
            ),
        })
    return pd.DataFrame(rows).set_index("id")


def estimate_kmer_histograms(
    reads: CasavaOneEightSingleLanePerSampleDirFmt,
    kmer_sizes: str = "21,31",
    max_frequency: int = 1000,
    threads: int = 1,
) -> tuple[qiime2.Metadata, qiime2.Metadata]:
    """Estimate aggregate k-mer histograms over all FASTQ files in a dataset."""
    sizes = _parse_kmer_sizes(kmer_sizes)
    input_paths = sorted(Path(str(reads)).glob("*.fastq.gz"))
    if not input_paths:
        raise ValueError("The Casava directory does not contain any .fastq.gz files.")

    with tempfile.TemporaryDirectory(prefix="q2-ntcard-") as tmpdir:
        prefix = Path(tmpdir) / "ntcard"
        command = [
            "ntcard", f"-k{','.join(map(str, sizes))}",
            f"-c{max_frequency}", f"-t{threads}", "-p", str(prefix),
            *map(str, input_paths),
        ]
        try:
            run_command(command, capture_output=True, text=True)
        except subprocess.CalledProcessError as error:
            detail = (error.stderr or error.stdout or "no diagnostic output").strip()
            raise RuntimeError(
                f"ntCard failed with exit code {error.returncode}: {detail}"
            ) from error
        histogram = _read_histograms(Path(tmpdir).glob("ntcard*.hist"), sizes)

    return qiime2.Metadata(_summarize(histogram)), qiime2.Metadata(histogram)
