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
rsomics-count normalize counts.tsv -o fpkm.tsv --mode fpkm
rsomics-count normalize counts.tsv -o fpkm-uq.tsv --mode fpkm-uq
```

Outputs are transactional. Parse, I/O, malformed alignment, annotation, pair,
and normalization errors leave existing outputs unchanged and exit non-zero.

This crate uses `rsomics-bamio` for borrowed BAM records, `rsomics-intervals`
for checked genomic geometry, `rsomics-common` for runtime and output
transactions, and `rsomics-help` for the shared CLI presentation. Counting
work is distributed through bounded worker queues; query-name-sorted pairs are
streamed without retaining the library in memory.

The recorded featureCounts 2.1.1 gate has exact normalized counts and summary
columns for the benchmark cases. The current implementation is slower in wall
time and uses substantially less peak resident memory; the measurements and
provenance are in [`benchmarks/2026-08-02`](benchmarks/2026-08-02/README.md).

## Installation and delivery

The recorded registry release is `rsomics-count 0.1.0`:

```console
cargo install --locked rsomics-count --version 0.1.0
```

Native CI prepares Linux and macOS binaries on both x86_64 and aarch64,
rebuilt and installed from a verified Cargo source package. These are delivery
candidates until a GitHub release and its checksums are actually published;
this README does not claim that pending release is already available.
Each native carrier includes its exact source/lock/compiler/host identities,
installed count-to-normalization lifecycle, dependency graph and complete
notice texts. Binaries are not signed or notarized; minimum OS and external
dynamic-library compatibility are not certified by this source inventory.

A source package newly generated at the GitHub delivery commit includes its
new documentation, notices and VCS identity. It is not the byte-identical
immutable registry 0.1.0 archive and does not republish or replace that archive.
No version, algorithm, public API or locked dependency is changed for delivery.

The accepted historical performance implementation is `3a0154b1`; later
`5bc764d8` adds evidence/CI without changing its source or lockfile. Recorded
wall time is 1.67–2.56 times featureCounts in the declared cases, with about
87% lower RSS in the two reported four-worker cases. Those are historical
measurements, not fresh timings of the new native carriers or universal gains.

CI retains all four native debug/release and frozen compatibility suites.
Pinned official featureCounts 2.1.1 runs live on Linux x86_64 and both macOS
architectures. Linux aarch64 has frozen compatibility, not a fourth live
oracle; local checks with FEATURECOUNTS unset do not execute upstream software.

CRAM, transcript-abundance inference, read-position projection, junction counts,
read-group stratification, long-read mode, per-record assignments,
largest-overlap and non-overlap limits remain outside this released slice.
See [attribution](ATTRIBUTION.md) and [third-party notices](THIRD_PARTY_NOTICES.md).

Local development uses `bash scripts/guarded-cargo.sh <cargo arguments>` with
the pinned external Rust toolchain and physically checked Cargo, temporary and
Count-specific Clang/ccache/sccache directories. The guard rejects inherited
compiler wrappers and flags, including host/target-specific native overrides,
before invoking Cargo. Run with those overrides unset; do not change shared
user configuration to satisfy the guard. Its `--exec` mode is a trusted
environment helper, not an arbitrary-command sandbox.

Licensed under MIT OR Apache-2.0.
