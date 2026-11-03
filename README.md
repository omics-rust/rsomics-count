# rsomics-count

`rsomics-count` assigns short-read SAM or BAM alignments to annotated genomic
features and writes one count matrix plus a complete assignment summary. It
also normalizes strict gene-by-sample matrices to TPM, FPKM, or FPKM-UQ.

The counting contract follows featureCounts 2.1.1 for the supported slice:
GTF, GFF3, and SAF annotations; feature or meta-feature counting; multiple
libraries; stranded and paired-fragment modes; overlap, mapping, and pair
filters; and full or fractional multi-overlap and multi-mapping contributions.

The unreleased 0.2.0 candidate corrects deletion-gap coverage and read-fraction calculations
in the 0.1.0 distribution. `--frac-overlap` uses original CIGAR geometry,
including insertion and soft-clip contributions; overlapping paired reads
use fragment unions rather than summed mate lengths. Target overlap also
counts shared mate coverage once. These corrections are not yet released.
The [read-fraction correction measurements](benchmarks/2026-11-03/README.md) at
`1e63643` preserve exact
counts and summaries on representative inputs, with 82–85% lower peak RSS
and slower wall time than featureCounts.

Unreleased source also supports `--read-extension5`, `--read-extension3` and
`--read2pos {3,5}`. Each mate is extended in its own biological orientation,
clipped to the alignment-header chromosome bounds, then optionally reduced
to one reference base. Original CIGAR length and split status are retained
for filtering. These API additions are included in the 0.2.0 candidate; they are
not present in the immutable 0.1.0 packages or covered by the earlier timings.
Its [fresh projection measurements](benchmarks/2026-11-03/projection.md) retain
the lower-RSS/slower-wall result with exact complete output agreement.

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
provenance are in the [recorded benchmark evidence](https://github.com/omics-rust/rsomics-count/blob/8ae7a2b0294557d78edba30c0a1aa9cfef1c0350/benchmarks/2026-08-02/README.md).

## Installation and delivery

Version 0.2.0 is a delivery candidate, not a public release. The installation
links below remain pinned to the immutable 0.1.0 distributions, which do not
include the corrections or projection options described above.

The [GitHub v0.1.0 release](https://github.com/omics-rust/rsomics-count/releases/tag/v0.1.0)
provides native binaries for Linux and macOS on x86_64 and aarch64, pinned to
source commit `8ae7a2b0294557d78edba30c0a1aa9cfef1c0350`.
Download your archive and [SHA256SUMS](https://github.com/omics-rust/rsomics-count/releases/download/v0.1.0/SHA256SUMS),
verify its checksum, extract it and place `rsomics-count` in your chosen binary
directory. For example, after downloading the
[macOS aarch64 archive](https://github.com/omics-rust/rsomics-count/releases/download/v0.1.0/rsomics-count-0.1.0-macos-aarch64.tar.gz):

```console
grep '  rsomics-count-0.1.0-macos-aarch64.tar.gz$' SHA256SUMS | shasum -a 256 -c -
tar -xzf rsomics-count-0.1.0-macos-aarch64.tar.gz
./rsomics-count-0.1.0-macos-aarch64/rsomics-count --version
```

The [clean GitHub source package](https://github.com/omics-rust/rsomics-count/releases/download/v0.1.0/rsomics-count-0.1.0-source.crate)
can be installed from its extracted directory:

```console
grep '  rsomics-count-0.1.0-source.crate$' SHA256SUMS | shasum -a 256 -c -
tar -xzf rsomics-count-0.1.0-source.crate
cargo install --locked --path rsomics-count-0.1.0
```

The existing registry distribution is separately installed with:

```console
cargo install --locked rsomics-count --version 0.1.0
```

All four native jobs rebuilt and installed the verified clean source package;
complete originals were adopted before publication. A fresh public-download
Mac aarch64 installation passes the nine-call count-to-all-normalizations
lifecycle, including complete matrix/summary, stdin and failure preservation.
Each native carrier includes its exact source/lock/compiler/host identities,
installed count-to-normalization lifecycle, dependency graph and complete
notice texts. Binaries are not signed or notarized; minimum OS and external
dynamic-library compatibility are not certified by this source inventory.

The 0.1.0 source package generated at the GitHub delivery commit includes its
new documentation, notices and VCS identity. It is not the byte-identical
immutable registry 0.1.0 archive and does not republish or replace that archive.
That historical delivery changed no version, algorithm, public API or locked dependency.
Original published archives retain their construction-time README/BUILD
wording and relative benchmark link; the release notes and this README record
the later public delivery. The tag and original assets are not replaced by
the 0.2.0 candidate preparation.

The accepted historical performance implementation is `3a0154b1`; later
`5bc764d8` adds evidence/CI without changing its source or lockfile. Recorded
wall time is 1.67–2.56 times featureCounts in the declared cases, with about
87% lower RSS in the two reported four-worker cases. Those are historical
measurements, not fresh timings of the native carriers or universal gains.
The candidate's accepted implementation and projection performance source is
`2b5e450`; its [projection measurements](benchmarks/2026-11-03/projection.md)
show 82–85% lower peak RSS and 1.43–1.64 times featureCounts wall time in the
declared cases. Delivery binaries are not independently remeasured.

CI retains all four native debug/release and frozen compatibility suites.
Pinned official featureCounts 2.1.1 runs live on Linux x86_64 and both macOS
architectures. Linux aarch64 has frozen compatibility, not a fourth live
oracle; local checks with FEATURECOUNTS unset do not execute upstream software.

CRAM, transcript-abundance inference, junction counts,
read-group stratification, long-read mode, per-record assignments,
largest-overlap and non-overlap limits remain outside the candidate slice.
See [attribution](ATTRIBUTION.md) and [third-party notices](THIRD_PARTY_NOTICES.md).
The notices explicitly label the team-owned MIT/Apache supplement for the
historical `rsomics-bamio 0.5.0` license-text omission; it is not recovered
original bamio text or legal clearance.

Local development uses `bash scripts/guarded-cargo.sh <cargo arguments>` with
the pinned external Rust toolchain and physically checked Cargo, temporary and
Count-specific Clang/ccache/sccache directories. The guard rejects inherited
compiler wrappers and flags, including host/target-specific native overrides,
before invoking Cargo. Run with those overrides unset; do not change shared
user configuration to satisfy the guard. Its `--exec` mode is a trusted
environment helper, not an arbitrary-command sandbox.

Licensed under MIT OR Apache-2.0.
