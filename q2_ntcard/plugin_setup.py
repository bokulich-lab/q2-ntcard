from qiime2.plugin import Citations, Int, Plugin, Range, Str, Visualization
from q2_types.metadata import ImmutableMetadata
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)

from q2_ntcard import __version__
from q2_ntcard._methods import estimate_kmer_histograms
from q2_ntcard._pipelines import estimate_and_visualize
from q2_ntcard._visualizers import visualize_kmer_histograms

citations = Citations.load("citations.bib", package="q2_ntcard")
reads_type = SampleData[
    SequencesWithQuality | PairedEndSequencesWithQuality
]
parameters = {
    "kmer_sizes": Str,
    "max_frequency": Int % Range(1, None),
    "threads": Int % Range(1, None),
}
parameter_descriptions = {
    "kmer_sizes": "Comma-separated positive k-mer lengths (for example, 21,31).",
    "max_frequency": "Maximum k-mer frequency represented in each histogram.",
    "threads": "Number of OpenMP threads used by ntCard.",
}

plugin = Plugin(
    name="ntcard", version=__version__, package="q2_ntcard",
    website="https://github.com/qiime2/q2-ntcard",
    description=("Estimate k-mer cardinality and coverage histograms from "
                 "demultiplexed sequencing reads using ntCard."),
    short_description="QIIME 2 wrapper for ntCard.",
)

plugin.methods.register_function(
    function=estimate_kmer_histograms,
    inputs={"reads": reads_type}, parameters=parameters,
    outputs=[("summary", ImmutableMetadata),
             ("histogram", ImmutableMetadata)],
    input_descriptions={
        "reads": "Single- or paired-end reads in Casava 1.8 directory format."
    }, parameter_descriptions=parameter_descriptions,
    output_descriptions={
        "summary": "One row of aggregate k-mer statistics per k-mer length.",
        "histogram": "Long-form estimated k-mer multiplicity histograms.",
    }, name="Estimate k-mer histograms with ntCard",
    description=("Run ntCard over all reads and return analysis-ready "
                 "ImmutableMetadata tables."),
    citations=[citations["Mohamadi-Khan-Birol-2017"]],
)

plugin.visualizers.register_function(
    function=visualize_kmer_histograms,
    inputs={"histogram": ImmutableMetadata}, parameters={"title": Str},
    input_descriptions={"histogram": "Histogram table produced by ntCard."},
    parameter_descriptions={"title": "Title displayed above the chart."},
    name="Visualize ntCard k-mer histograms",
    description="Interactively inspect estimated k-mer multiplicity curves.",
    citations=[citations["Mohamadi-Khan-Birol-2017"]],
)

plugin.pipelines.register_function(
    function=estimate_and_visualize,
    inputs={"reads": reads_type},
    parameters={**parameters, "title": Str},
    outputs=[("summary", ImmutableMetadata),
             ("histogram", ImmutableMetadata),
             ("visualization", Visualization)],
    input_descriptions={"reads": "Single- or paired-end demultiplexed reads."},
    parameter_descriptions={
        **parameter_descriptions, "title": "Title displayed above the chart."
    }, output_descriptions={
        "summary": "Aggregate k-mer statistics.",
        "histogram": "Long-form k-mer multiplicity histograms.",
        "visualization": "Interactive k-mer histogram chart.",
    }, name="Estimate and visualize k-mer histograms",
    description="Run ntCard and immediately create an interactive report.",
)
