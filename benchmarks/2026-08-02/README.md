# featureCounts 2.1.1 compatibility and performance gate

The benchmarked source is commit
`3a0154b15c8a985b999a98c4dd61d07a7533dec3`. Its release binary SHA-256 is
`3ba07d69818c8c805649a724079e0989411177dfe1ff05fbd2d32e67e3c22f1d`.

The machine was an Apple M2 running macOS 26.6 build 25G72. The toolchain was
Homebrew Rust 1.97.1, hyperfine 1.20.0, and samtools 1.24. The official
featureCounts 2.1.1 arm64 binary SHA-256 was
`3249dea2e90721618f8af4fb547ff1e9bab987620762c1e24148b149487f8ea5`;
its downloaded archive SHA-256 was
`58f4d9a08d88334871fd89a1a87df1fc5560fe6ebad2e4527514cd7f8e3bf49a`.

The synthetic single-end fixture has 500 genes, 1,500 exons, and 2,000,000
BAM records. The paired fixture has the same annotation and 1,000,000
query-name-sorted fragments. Fixture SHA-256 values are:

| File | SHA-256 |
|---|---|
| `genes.gtf` | `2e43ca96b0be1d328b300c33aa386fb053fd912c2b2f93f46860dcc5b976c914` |
| `reads.bam` | `514e875b9337f24eb847034bb7575cb08c990c3ad717a9f1218e0f377d8d3596` |
| `paired-reads.bam` | `66c90e8e8e0dbf986be94768d450c4ac899c3f0aeb1dc14035d051da86e0dbc8` |

Hyperfine used two warmups and ten measured runs for the single and policy
cases. The paired and two-library cases used three warmups and fifteen runs.
Each command wrote its count matrix and summary during the timed interval.

| Case | Threads | featureCounts wall, ms | rsomics wall, ms | rsomics / featureCounts | rsomics CPU, ms | featureCounts CPU, ms |
|---|---:|---:|---:|---:|---:|---:|
| single | 1 | 705.4 ± 5.2 | 1,175.3 ± 46.9 | 1.67× | 737.2 | 693.5 |
| single | 4 | 322.7 ± 2.7 | 812.5 ± 50.3 | 2.52× | 730.5 | 868.3 |
| paired | 4 | 352.6 ± 5.4 | 903.8 ± 56.2 | 2.56× | 1,096.7 | 913.6 |
| overlap, multimap, fractional, sense | 4 | 335.7 ± 10.5 | 796.5 ± 36.3 | 2.37× | 703.4 | 898.1 |
| two libraries | 4 | 662.2 ± 33.0 | 1,160.2 ± 51.4 | 1.75× | 1,446.1 | 1,767.7 |

The current implementation does not have a throughput advantage. It used
less total CPU in the four-thread single, policy, and two-library cases, but
more CPU in the one-thread and paired cases. Its transactional outputs include
file and directory synchronization; the table remains the end-to-end result
without adjusting for that contract.

Five `/usr/bin/time -l` runs measured median peak resident sizes of 68,108,288
bytes for featureCounts and 8,503,296 bytes for rsomics in the single case, an
87.5% reduction. Paired medians were 68,222,976 and 8,634,368 bytes, an 87.3%
reduction. The per-run measurements are in `rss-final.tsv`. This material
memory advantage, rather than throughput, is the performance gate satisfied by
the release.

For every timed case, normalized feature rows and all assignment-summary
columns matched featureCounts exactly. `compatibility-final.tsv` records the
shared normalized hashes. The JSON files retain every timing sample.

The command families were:

```console
featureCounts -T THREADS [POLICY] -a genes.gtf -o counts.tsv reads.bam...
rsomics-count reads.bam... -a genes.gtf -o counts.tsv --summary summary.tsv --threads THREADS [POLICY]
```

Paired commands added `-p --countReadPairs` or `--paired`. The policy case used
`-O -M --fraction -s 1` and
`--multi-overlap fractional --multi-map fractional --strand sense`.
