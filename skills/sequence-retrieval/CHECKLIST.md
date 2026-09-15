# Sequence Retrieval Checklist

- [ ] Confirm organism, gene/locus, molecule, strain and intended isoform.
- [ ] Verify accession/version from the provider; preserve explicit requested versions.
- [ ] Route all RefSeq prefixes, including NZ_, to NCBI; verify INSDC alternatives.
- [ ] Use the protein database for proteins and genome packages for GCA_/GCF_ assemblies.
- [ ] Report total hits separately from the retrieved page; handle empty results.
- [ ] Confirm identity from titles/annotations, not just All Fields strain matches.
- [ ] Report reference/prediction status and evidence without inferring quality from prefix alone.
- [ ] Include organism, strain, molecule, biological length, topology and a short sequence preview.
- [ ] Count sequence symbols rather than FASTA header/newline characters.
- [ ] Summarize CDS/gene/tRNA/rRNA/regulatory annotations where present; mark missing fields unknown.
- [ ] Verify RefSeq/INSDC, BioProject and BioSample links and explain relevant alternatives.
- [ ] Preserve saved files and manifests; include provider, retrieval date, accession and hashes.
- [ ] Identify format conversion explicitly; it does not establish cross-database identity.
- [ ] Respect output protection and response limits; distinguish service failure from no results.
- [ ] Provide FASTA, annotated-format options, direct database links and appropriate citations.
