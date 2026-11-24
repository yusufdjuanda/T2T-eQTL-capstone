# metadata

Small text files that describe the data. Written by people (or by the scripts named below) and committed.

## sources.tsv

Where every input comes from. One row per resource.

| Column | Meaning |
|--------|---------|
| id | short name used in the docs |
| resource | what it is |
| origin | who made it (paper or project) |
| url_or_accession | download URL, DOI or archive accession |
| version_or_retrieved | release, or the date it was downloaded |
| local_path | where it goes in the project (`data/raw/...`) |
| obtained_by | the script or pipeline process that downloads it |
| published_checksum | checksum given by the provider (or by us, for the release asset), where there is one |
| terms | licence or use conditions |

`processed-mage` is the input of `pixi run reproduce`: the release asset made by the cluster run. Its archive checksum
is filled in when the asset is published (see `docs/cluster.md`, "Publishing the bundle").

## checksums.sha256

SHA-256 of every file of the unpacked release asset (`data/raw/mage-processed-v1/`), written by the pipeline
(process EXPORT_BUNDLE, `results/release/checksums.sha256`). `pixi run download` checks the files against it.

## checksums-cluster-inputs.sha256

SHA-256 of every input of the cluster run (genotype panels, genomes, annotations, T2T consortium tracks, liftOver
chains), paths under `data/raw/`. Check with `pixi run check-inputs [raw folder]`.

## expected-results.tsv

The published values that `pixi run check-results` (the last step of `pixi run reproduce`) compares the results with:
`metric`, `expected`, `tolerance` (0 = exact) and `description`. `changed_genes_sha256` fixes the exact list of changed
genes and their tiers.

## samples.tsv

One row per MAGE individual; written by `scripts/make_sample_lists.sh`.

| Column | Meaning |
|--------|---------|
| sample | 1000 Genomes individual ID |
| mage | in MAGE (always `yes`) |
| t2t_panel | has genotypes in the phased 1000 Genomes T2T-CHM13 panel |
| included | used in this study (`yes` = in MAGE and in the panel); 731 individuals |

## samplesheet-mage.csv

One row per MAGE sequencing run of an included individual (779 runs; 24 individuals have 3 runs); written by
`scripts/make_samplesheet.py` from the ENA run report of PRJNA851328.

| Column | Meaning |
|--------|---------|
| sample | individual (runs of one individual are summed per gene) |
| run | SRA run accession |
| read_pairs | read pairs according to ENA; the pipeline checks the downloaded run has exactly this many |
| read_length | read length (150 for all of MAGE) |
| sex | from the ENA sample title; a covariate of the eQTL model |
| ancestry | 1000 Genomes population (super-population) |

## chrom-names-chm13v2.0.tsv

RefSeq (`NC_...`) and UCSC (`chr...`) names of the T2T-CHM13v2.0 chromosomes. The genotype panels use `chr1`; the
NCBI T2T genome uses `NC_060925.1`.
