# Attribution and compatibility sources

The counting compatibility profile is validated against featureCounts 2.1.1,
part of the Subread package, copyright the Subread authors and distributed
under the GNU General Public License version 3. The relevant public contracts
are documented by the Subread project, its featureCounts documentation, and
the Subread Users Guide:

- <https://subread.sourceforge.net/>
- <https://subread.sourceforge.net/featureCounts.html>
- <https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf>

`rsomics-count` is an independent Rust implementation. It does not link,
bundle, or copy Subread source code. CI downloads pinned official featureCounts
binaries only to run live compatibility tests. Repository fixtures are
synthetic and team-owned.
