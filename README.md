# q2-ntcard

`q2-ntcard` is a QIIME 2 community plugin that wraps
[ntCard](https://github.com/BirolLab/ntCard), a streaming estimator of k-mer
cardinality and multiplicity histograms.

The plugin exposes one public pipeline, backed by two private actions:

- `_estimate` runs ntCard on single- or paired-end
  `SampleData` backed by `CasavaOneEightSingleLanePerSampleDirFmt`. It returns
  summary and long-form histogram tables as `ImmutableMetadata` artifacts.
- `_visualize` renders all requested histograms as interactive Vega-Lite
  facets with logarithmic axes, F0 estimates, and hover details.
- `estimate-kmer-histograms` combines the method and visualizer in one pipeline.

## Installation

Create a fresh QIIME 2 tiny environment from one of the files in
`environment-files`, then install this repository:

```shell
conda env create -n q2-ntcard-dev \
  --file environment-files/q2-ntcard-rachis-tiny-release.yml
conda activate q2-ntcard-dev
python -m pip install -v .
qiime dev refresh-cache
qiime ntcard --help
```

## Example

```shell
qiime ntcard estimate-kmer-histograms \
  --i-reads demux.qza \
  --p-kmer-sizes 21,31,51 \
  --p-threads 4 \
  --o-summary summary.qza \
  --o-histogram histogram.qza \
  --o-visualization histogram.qzv
```

The result tables can be exported as metadata TSV files or supplied to actions
that accept QIIME 2 Metadata.
