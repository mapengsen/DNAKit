# DNAKit execution: ukb-ppp-region-fetch

The module ukb_ppp_region_fetch.py exports UKBPPPClient.fetch_region; inspect --help for the CLI. A SORT1/EUR GRCh38 slice is bundled for offline use. Arbitrary regions use Synapse and need a PAT with access to syn51364943; the public-looking S3 URL can return AccessDenied. Software is MIT (Aviv Madar); bundled summary-statistic derivatives are CC-BY-4.0 with Sun et al. 2023 attribution. Sources: https://registry.opendata.aws/ukbppp/ and https://doi.org/10.1038/s41586-023-06592-6 . This is summary-level pQTL data, not individual UK Biobank records.

Install this complete folder with `dnakit skills install ukb-ppp-region-fetch --destination ./analysis-skills`. Inspect optional dependencies with `dnakit skills doctor ukb-ppp-region-fetch`. Resolve script/resource paths against the installed folder.

Import validation covers metadata, file integrity and Python syntax. It does not establish scientific accuracy, service availability or a successful pipeline run. Cite the upstream authors and the actual data/methods used; see UPSTREAM.md.
