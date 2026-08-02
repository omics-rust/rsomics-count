# rsomics-count

`rsomics-count` assigns short-read SAM or BAM alignments to annotated genomic
features and writes one count matrix plus a complete assignment summary. It
also normalizes strict gene-by-sample matrices to TPM, FPKM, or FPKM-UQ.

The counting contract follows featureCounts 2.1.1 for the supported slice:
GTF, GFF3, and SAF annotations; feature or meta-feature counting; multiple
libraries; stranded and paired-fragment modes; overlap, mapping, and pair
filters; and full or fractional multi-overlap and multi-mapping contributions.

```console
rsomics-count sample.bam -a genes.gtf -o counts.tsv --summary counts.summary.tsv
rsomics-count normalize counts.tsv -o tpm.tsv --mode tpm
```

Outputs are transactional. Parse, I/O, malformed alignment, annotation, pair,
and normalization errors leave existing outputs unchanged and exit non-zero.

This crate uses `rsomics-bamio` for borrowed BAM records, `rsomics-intervals`
for checked genomic geometry, `rsomics-common` for runtime and output
transactions, and `rsomics-help` for the shared CLI presentation. Counting
work is distributed through bounded worker queues; query-name-sorted pairs are
streamed without retaining the library in memory.

Licensed under MIT OR Apache-2.0.
