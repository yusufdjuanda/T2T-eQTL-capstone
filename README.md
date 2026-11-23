# Does the reference genome change eQTLs? T2T-CHM13 vs GRCh38 in MAGE

[![Reproduce](https://github.com/yusufdjuanda/T2T-eQTL-capstone/actions/workflows/reproduce.yml/badge.svg)](https://github.com/yusufdjuanda/T2T-eQTL-capstone/actions/workflows/reproduce.yml)
[![Tests](https://github.com/yusufdjuanda/T2T-eQTL-capstone/actions/workflows/tests.yml/badge.svg)](https://github.com/yusufdjuanda/T2T-eQTL-capstone/actions/workflows/tests.yml)

Author: Yusuf Maulana. Capstone project of CBIO7501 (Reproducible and Open Science): one dataset (MAGE) and one
comparison (T2T-CHM13 vs GRCh38) from a larger ongoing study.

## What is this project about?

An eQTL is a genetic variant that changes how much a gene is expressed. Finding one starts by mapping RNA-seq reads
to a reference genome. GRCh38, the standard reference, has gaps, collapsed duplications and a few false duplications;
T2T-CHM13 is the first complete human genome. **If the same reads from the same people are analysed on both
references, which eQTLs change, and where in the genome are they?**

Data: RNA-seq of lymphoblastoid cell lines from 731 people of the MAGE study (26 populations, 2 x 150 bp) and their
1000 Genomes genotypes, called separately on each reference. The same reads go through the same pipeline twice, once
per reference: STAR alignment, STAR and salmon quantification, cis-eQTL mapping with tensorQTL. Variants called on both
references are paired, so a gene's eQTL effect is compared **at the same variant**. A gene counts as *changed* only
when a paired test (which accounts for the two estimates coming from the same people) is significant under both STAR
and salmon, in the same direction.

The project has two parts, one Nextflow pipeline (`workflow/nextflow/main.nf`):

| Part | From | To | Needs |
|------|------|----|-------|
| `--mode reproduce` (what `pixi run reproduce` runs) | the processed MAGE inputs of the paired test, a 150 MB release asset | paired test, changed-gene sets, figures, report | a laptop, minutes |
| `--mode full` | raw reads (SRA) and genotype panels | the release asset, then the same analysis | a Slurm cluster, ~3,000 CPU hours |

## Is there a publication?

Not yet. The data are from Taylor et al. 2024, *Nature* 632:122–130 (MAGE,
<https://doi.org/10.1038/s41586-024-07708-2>); the reference genome from Nurk et al. 2022, *Science* 376:44–53
(<https://doi.org/10.1126/science.abj6987>). Behera et al. 2023 (*Genome Biology* 24:31, FixItFelix) showed how
GRCh38's false duplications distort eQTLs on chr21.

## Reproduce the results

You need [Pixi](https://pixi.prefix.dev/latest/) (0.67 or later) and an internet connection; Linux or macOS on Apple
silicon (on Windows, use WSL2). Then:

```bash
git clone https://github.com/yusufdjuanda/T2T-eQTL-capstone.git
cd T2T-eQTL-capstone
pixi run reproduce
```

`pixi run reproduce` installs the exact software in `pixi.lock`, then:

1. `pixi run download` downloads the processed MAGE bundle (GitHub release `data-v1`), checks the archive against
   `metadata/sources.tsv` and every file against `metadata/checksums.sha256`, and makes it read-only;
2. `pixi run pipeline` runs the Nextflow pipeline in reproduce mode with the settings in `config/params.yaml`:
   paired test, changed-gene sets, figures (`results/analysis/`);
3. `pixi run report` renders the Quarto report, `reports/t2t-eqtl-report.html`.

The *Reproduce* GitHub Action runs the same command on a fresh machine on every push.

Tests (no data needed): `pixi run test` runs the unit tests of the scripts (pytest), the tests of the Nextflow
processes, workflows and both pipeline modes (nf-test), and a stub run of the full cluster pipeline that checks every
step from SRA download through alignment and eQTL mapping to the export. `pixi task list` shows every task.

The full analysis from raw reads, on a Slurm cluster: `docs/cluster.md`.

## What hardware does it need?

- **`pixi run reproduce`**: any recent laptop or the GitHub runner (1–4 CPUs, under 4 GB of memory). Disk: 2.1 GB for
  the software environment (measured on Linux, 2026-11-23), about 150 MB for the bundle, a few MB of results. The
  analysis takes minutes; the first `pixi install` takes longer.
- **Full pipeline** (`--mode full`): a Slurm cluster. About 3,000 CPU hours (alignment and quantification of 779 runs on
  two references), up to 16 CPUs and 48 GB per task (STAR), 64 GB for the genotype steps. Disk: ~60 GB of inputs, ~35 GB
  per STAR reference and index, the reads (~2.5 TB) pass through the work folder one run at a time, and ~40 GB of eQTL
  results. Details and timings: `docs/cluster.md`.

## What outputs should I expect?

| Path | Content |
|------|---------|
| `reports/t2t-eqtl-report.html` | the report: question, data, methods, results, calibration, software versions |
| `results/analysis/paired-test/` | per gene and method: effects on both references, paired z and p (`effects-star.tsv`, `effects-salmon.tsv`), the calls under both methods (`robust.tsv`), bootstrap calibration (`calibration.tsv`), `summary.tsv` |
| `results/analysis/gene-sets/changed-genes.tsv` | the changed genes, with their tier (core, main, extent) and region class |
| `results/analysis/figures/` | fig01 eGenes per reference, fig02 effects at the same variant, fig03 changed genes by region class (PNG and PDF) |

Expected, from the published run (731 individuals): 13,563 genes compared under both STAR and salmon; 121 genes in
the main changed set (63 of them core), 59 of them in segmental duplications; 246 in the wider extent set. The paired
standard error matches a bootstrap over individuals (ratio 0.96), while the independent one overstates it about
9-fold (ratio 0.11).

## What is in this repository?

The layout of the course template; every folder has its own README, and `docs/structure.md` describes them all.

```text
workflow/nextflow/   the pipeline: main.nf, modules/ (processes), subworkflows/, conf/ (profiles local, hpc, test), tests/ (nf-test)
scripts/             the code each step runs (Python, R, shell), one script per step
config/              params.yaml (analysis settings); site.config.example (cluster paths and queues)
metadata/            sources.tsv, checksums, samples.tsv, samplesheet-mage.csv
reports/             the Quarto report
tests/               pytest unit tests and small test data
docs/                structure, cluster guide, design decisions, software versions
data/  results/      downloaded data and outputs (not in Git)
containers/  envs/   Docker/Apptainer recipes from pixi.lock; Conda export; vendored tensorQTL
```

## Data

| Data | Source | Licence / terms |
|------|--------|-----------------|
| MAGE RNA-seq reads (779 runs, 731 individuals) | NCBI SRA BioProject PRJNA851328 | public, no use restrictions |
| Phased 1000 Genomes genotypes on T2T-CHM13v2.0 | <https://doi.org/10.5281/zenodo.7612953>, human-pangenomics S3 | IGSR data reuse policy |
| Phased 1000 Genomes genotypes on GRCh38 (NYGC, 2022-04-22) | 1000 Genomes FTP | IGSR data reuse policy |
| T2T-CHM13v2.0 and GRCh38 genomes, RefSeq annotation RS_2025_08 | NCBI RefSeq | public |
| T2T segmental-duplication and new-sequence tracks; UCSC liftOver chains | T2T consortium; UCSC | public |
| Processed MAGE bundle (input of `pixi run reproduce`) | this project, GitHub release `data-v1` | CC BY 4.0 |

Every input with its URL or accession, version, local path, download script and checksum: `metadata/sources.tsv` and
`metadata/README.md`. No data is stored in this repository.

**Data availability statement.** Raw RNA-seq reads are in the NCBI Sequence Read Archive under BioProject
PRJNA851328. Genotypes are from the 1000 Genomes Project (T2T-CHM13v2.0 calls: Zenodo
10.5281/zenodo.7612953; GRCh38 calls: IGSR, 20220422 release). The processed data needed to reproduce every result
(normalised expression, covariates and genotype dosages of the compared variants, and the per-gene effect tables) are
attached to release `data-v1` of this repository, with their sha256 checksums in `metadata/`. Code and environment:
this repository (MIT).

## How to cite

See `CITATION.cff` (GitHub: "Cite this repository").

## License

Code: MIT (`LICENSE`). Documentation, report and figures: CC BY 4.0 (`LICENSE-docs.md`). The processed data bundle:
CC BY 4.0. Input data keep the terms of their providers (table above).

## Contact

Questions and problems: open an issue on GitHub. Contributions: see `CONTRIBUTING.md`; everyone taking part follows
`CODE_OF_CONDUCT.md`.
