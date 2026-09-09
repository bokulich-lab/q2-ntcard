from pathlib import Path
import subprocess
from unittest.mock import patch

import pandas.testing as pdt
import pytest

from q2_ntcard._methods import (
    _parse_kmer_sizes,
    _read_histograms,
    _summarize,
    _estimate,
)

DATA = Path(__file__).parent / "data"


def test_parse_kmer_sizes():
    assert _parse_kmer_sizes("21, 31") == [21, 31]
    with pytest.raises(ValueError, match="duplicates"):
        _parse_kmer_sizes("21,21")


def test_parse_and_summarize_fixture():
    histogram = _read_histograms([DATA / "ntcard_k21.hist"], [21])
    summary = _summarize(histogram)
    assert summary.loc["k21", "estimated_distinct_kmers"] == 52
    assert summary.loc["k21", "estimated_total_kmers"] == 66
    assert summary.loc["k21", "estimated_singleton_kmers"] == 40


def test_parse_released_ntcard_format():
    histogram = _read_histograms([DATA / "ntcard_legacy_k21.hist"], [21])
    summary = _summarize(histogram)
    assert summary.loc["k21", "estimated_distinct_kmers"] == 52
    assert summary.loc["k21", "estimated_total_kmers"] == 66


def test_estimate_command_and_outputs(tmp_path):
    (tmp_path / "sample_S1_L001_R1_001.fastq.gz").touch()

    def fake_run(command, **kwargs):
        prefix = Path(command[command.index("-p") + 1])
        (prefix.parent / "ntcard_k21.hist").write_text(
            (DATA / "ntcard_k21.hist").read_text()
        )

    with patch("q2_ntcard._methods.run_command", side_effect=fake_run) as run:
        summary, histogram = _estimate(tmp_path, kmer_sizes="21", threads=2)
    assert run.call_args.args[0][0] == "ntcard"
    assert "-t2" in run.call_args.args[0]
    pdt.assert_frame_equal(
        summary.to_dataframe(),
        _summarize(histogram.to_dataframe()),
        check_dtype=False,
    )


def test_estimate_failure(tmp_path):
    (tmp_path / "sample_S1_L001_R1_001.fastq.gz").touch()
    error = subprocess.CalledProcessError(2, ["ntcard"], stderr="bad reads")
    with patch("q2_ntcard._methods.run_command", side_effect=error):
        with pytest.raises(RuntimeError, match="bad reads"):
            _estimate(tmp_path, kmer_sizes="21")
