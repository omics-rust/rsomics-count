# Original-CIGAR fraction correction

The current source repairs read/fragment fraction counting; it does not change
the public API, dependencies or package version. The previous 0.1.0 timing
evidence remains historical and is not reused as a measurement of this code.

Apple M2, 8 GiB, macOS27.2 arm64; Rust1.91.0; samtools/htslib1.24; official
featureCounts2.1.1 arm64. The existing `generate_fixture.rs` creates500 genes,
1500 exons and either2,000,000 single-end BAM records or1,000,000 query-name
sorted fragments. BAM compression identities differ from the historical
fixtures; the annotation hash is unchanged.

Both tools use four threads and write the complete matrix and summary during
timing. Each case has one warmup and five measured pairs, alternating tool
order. Paired commands add `--fracOverlap .45 -p --countReadPairs` upstream
and `--frac-overlap .45 --paired` natively. The single case uses defaults.
Wall time is the elapsed subprocess invocation, including `/usr/bin/time`;
peak RSS comes from its macOS `-l` report, in bytes.

| Case | Native wall median, s | Source wall median, s | Native RSS median, bytes | Source RSS median, bytes | RSS reduction |
|---|---:|---:|---:|---:|---:|
| single defaults |0.469265667|0.317002375|10,108,928|68,501,504|85.24%|
| paired fraction |0.571420584|0.341280792|12,107,776|69,091,328|82.48%|

This passes the resource-use gate, not a throughput gate: native wall time is
1.48× and1.67× source time. Every measured native RSS is below every source
RSS in its case. All12 pairs, including warmups, match every normalized
feature row and every assignment-summary category exactly. These synthetic
fixtures do not establish universal or real-cohort performance.

The [sample record](samples.json) retains all measured distributions and
input/binary/original-packet identities. Full argv, warmups, output matrices,
summaries, stdout, stderr and source-file identities are preserved externally
at `/Volumes/KIOXIA/Developments/tmp/count-read-fraction-perf-20261103.jmlsTpR5`.
An independent replay rehashes66 identities and rechecks all12 complete
output pairs and RSS reports. It does not infer performance of CI-built or
future published binaries.
