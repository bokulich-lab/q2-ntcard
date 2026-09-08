def estimate_and_visualize(
    ctx, reads, kmer_sizes="21,31", max_frequency=1000, threads=1,
    title="ntCard k-mer coverage histograms",
):
    estimate = ctx.get_action("ntcard", "estimate_kmer_histograms")
    visualize = ctx.get_action("ntcard", "visualize_kmer_histograms")
    summary, histogram = estimate(
        reads=reads, kmer_sizes=kmer_sizes,
        max_frequency=max_frequency, threads=threads,
    )
    visualization, = visualize(histogram=histogram, title=title)
    return summary, histogram, visualization
