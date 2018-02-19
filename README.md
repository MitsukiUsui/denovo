# What's this?
Detect de novo gene birth by overprinting.

## Workflow
These steps below are integrated by `./meta.sh`. 

1. download
Fetch genomes & annotations from RefSeq ftp server.<br>
Issue "orf_id" for each CDS, which is a unique identifier connecting `.gff`, `.fna`, and `.faa`.

2. represent (should be phylo in future)
Run PhyloPhlAn to reconstruct a phylogenetic tree.<br>
Select representative strains if necessary.

3. prodigal
Run Prodigal to re-annotate genomes.<br>
