# Does the reference genome change eQTLs? T2T-CHM13 vs GRCh38 in MAGE

**Question.** An eQTL is a genetic variant that changes how much a gene is expressed. Finding one starts by mapping
RNA-seq reads to a reference genome. GRCh38, the standard reference, has gaps, collapsed duplications and a few false
duplications. T2T-CHM13 is the first complete human genome. This project asks: **if the same reads and the same
people are analysed on both references, which eQTLs change, and where in the genome are they?**

**Data.** RNA-seq of lymphoblastoid cell lines from 731 people of the MAGE study (Taylor et al. 2024, 26 populations,
2 x 150 bp) and their 1000 Genomes genotypes, which were called separately on each reference.

**Approach.** The same reads go through the same pipeline twice, once per reference: STAR alignment and salmon
quantification, then cis-eQTL mapping with tensorQTL. Variants called on both references are paired up, so a gene's
eQTL effect can be compared **at the same variant** on T2T and on GRCh38. A gene counts as *changed* only when a
paired test (which accounts for the two estimates coming from the same people) is significant under both STAR and
salmon, in the same direction.

## Pipeline

| Stage | What it does | Entry point |
|---|---|---|
| 00 | set-up: storage folders, font, checksums, lint | `pixi run setup`, `pixi run fonts`, `pixi run -e lint lint` |
| 01 | genotypes: download both 1000 Genomes panels, keep the MAGE individuals | `scripts/01-genotypes/` |
| 02 | references: STAR references and indexes for T2T-CHM13 and GRCh38 (same RefSeq annotation release) | `scripts/02-references/` |
| 04 | expression: SRA → STAR + salmon on both references (Nextflow), then gene × individual matrices | `bash scripts/04-expression/submit-cohort.sh mage` |
| 05 | cis-eQTL mapping with tensorQTL, chr1–22, both references and both quantifiers | `bash scripts/05-eqtl/submit-eqtl-cohort.sh mage` |
| 06 | comparison: variant pairing, gene classes, eGene and effect comparison, paired test, changed-gene sets | `scripts/06-compare/` |
| 07 | figures | `pixi run figures` |

Stage numbers are shared by `scripts/`, `results/` and `logs/` (there is no stage 03). Every script describes its
inputs, outputs and usage in its header.

### Order of the steps

```bash
# 01 genotypes
bash scripts/01-genotypes/chm13v2.0/fetch-panel.sh
scripts/submit.sh scripts/01-genotypes/grch38/fetch-panel.sbatch
bash scripts/01-genotypes/make-sample-lists.sh
scripts/submit.sh scripts/01-genotypes/chm13v2.0/subset-genotypes.sbatch
scripts/submit.sh scripts/01-genotypes/grch38/subset-genotypes.sbatch
scripts/submit.sh scripts/06-compare/match-variants.sbatch          # pairs T2T and GRCh38 variants; stage 05 needs it

# 02 references
scripts/submit.sh scripts/02-references/prepare-star-references.sbatch
scripts/submit.sh scripts/02-references/build-star-index.sbatch 149   # MAGE reads are 150 bp
pixi run python scripts/02-references/make-region-beds.py

# 04 expression
pixi run python scripts/04-expression/make-cohort-sheets.py
bash scripts/04-expression/submit-cohort.sh mage

# 06a gene annotations used by the comparison
bash scripts/06-compare/compare-annotation-overlaps.sh
pixi run python scripts/06-compare/make-gene-classes.py

# 05 eQTL mapping and the first comparison steps (one chain of Slurm jobs with dependencies)
bash scripts/05-eqtl/submit-eqtl-cohort.sh mage

# 06b paired test and changed-gene sets, then 07 figures
pixi run python scripts/06-compare/paired-test.py --cohort mage
pixi run python scripts/06-compare/define-gene-sets.py --cohort mage
pixi run figures
```

To test stage 05 cheaply, `scripts/05-eqtl/make-simulated-cohort.py` builds a small simulated cohort with planted
eQTLs on chr22; a working pipeline finds most of them.

## Requirements

- Linux x86-64 and [pixi](https://pixi.sh) ≥ 0.67 (tested with 0.67.2). Every tool (Python 3.11, R, tensorQTL, STAR,
  salmon, samtools, bcftools, plink2, Nextflow, …) comes from `pixi.lock`.
- A Slurm cluster for stages 01–06. Expression takes about 3,000 CPU hours for MAGE on two references. The reads
  (~2.5 TB) pass through scratch one batch (~160 GB) at a time. STAR indexes need ~30 GB per reference.

## Quick start

```bash
git clone https://github.com/yusufdjuanda/T2T-eQTL-capstone.git && cd T2T-eQTL-capstone
pixi install --locked                          # default environment; `-e bio` for STAR, salmon & co.
cp config/site.env.example config/site.env     # edit: storage folders, scratch, Slurm partitions
pixi run setup                                 # data/, results/, logs/, work/ as links or plain folders
pixi run -e lint lint                          # static checks (also run by CI on every push)
```

`config/site.env` is the only file that knows about your machine. Scripts contain no absolute paths: they find the
project root and read `config/site.env` through `scripts/lib/common.sh`. Slurm jobs are submitted from the project
root with `scripts/submit.sh`, which adds the partition and account from `config/site.env`.

## Input data

Nothing is redistributed here. `config/data-sources.tsv` lists every input (origin, URL or accession, version, local
path) and `config/checksums/inputs.sha256` the sha256 of every downloaded file (`pixi run verify-data`):

- 1000 Genomes phased genotypes on T2T-CHM13 (Zenodo 10.5281/zenodo.7612953) and on GRCh38 (NYGC, 2022-04-22)
- MAGE RNA-seq (NCBI SRA PRJNA851328)
- NCBI RefSeq genomes and annotation release RS_2025_08 for both references
- T2T consortium segmental-duplication and new-sequence tracks; UCSC liftOver chains (used only to pair variants)

## Software environments (`pixi.toml`)

| Environment | Contents | Used by |
|---|---|---|
| `default` | Python 3.11 (pandas, scipy, PyTorch CPU, pgenlib), R (ggplot2, data.table), plink2, liftOver | Python and R scripts |
| `bio` | STAR 2.7.11b, salmon 1.11.4, samtools, bcftools, bedtools, sra-tools, minimap2, Nextflow | shell and Slurm scripts, Nextflow tasks |
| `lint` | Python, ShellCheck | `pixi run -e lint lint`, CI |

tensorQTL 1.0.10 and qtl 0.1.10 are vendored from their PyPI wheels in `envs/vendor/`.

## Conventions

- **Chromosome names.** The genotype panels use `chr1`; the NCBI T2T genome uses `NC_…`. Convert with
  `config/chrom-names-chm13v2.0.tsv`.
- **Alleles.** REF and ALT differ between T2T and GRCh38 at about 20% of common SNPs, because the references carry
  different bases there. Effects are always compared per T2T ALT allele, with the sign set by the correlation of the
  two variants' genotypes.
- **Temporary files.** Every job uses its own folder under `$T2T_TMP`, never `/tmp`.

## Licence

Code: MIT ([`LICENSE`](LICENSE)). Documentation and figures: CC BY 4.0 ([`LICENSE-docs.md`](LICENSE-docs.md)).
Citation: [`CITATION.cff`](CITATION.cff).
