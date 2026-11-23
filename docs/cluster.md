# The full analysis on a cluster

`pixi run reproduce` starts from the processed bundle. This page is for rebuilding that bundle from the raw reads:
the `--mode full` pipeline on a Slurm cluster.

## 1. Software

```bash
pixi install --locked                  # default: Nextflow, Python, R
pixi install --locked -e bio -e qtl    # STAR, salmon, samtools, bcftools, ... ; tensorQTL, plink2 (Linux only)
```

Compute nodes without internet access are fine: install on the login node first. Tasks find the `bio` and `qtl` tools
through `PIXI_PROJECT_ROOT`, so always launch the pipeline through `pixi run` (see `workflow/nextflow/conf/base.config`).

## 2. Site settings

```bash
cp config/site.config.example config/site.config     # not committed
```

Set where the raw inputs are (`params.raw_dir`), where STAR references go (`params.processed_dir`), the output folder
(`outputDir`), the work folder on fast scratch (`workDir`) and the Slurm queues. Nothing else in the project knows
about your machine.

## 3. Inputs (about 60 GB)

Under `raw_dir`, with the layout of `metadata/sources.tsv` (column `local_path`, without the `data/raw/` prefix):

```bash
pixi run fetch-t2t-panel   <raw_dir>/external/1kgp-chm13v2.0-phased-shapeit5-v1.1
pixi run fetch-grch38-panel <raw_dir>/external/1kgp-grch38-nygc-phased-20220422      # hours; use tmux
```

Genomes, annotations, T2T tracks, liftOver chains and the ENA run report: download from the URLs in
`metadata/sources.tsv`. Then check every file:

```bash
pixi run check-inputs <raw_dir>
```

The sample list and samplesheet are committed; to rebuild them from the ENA report and the panel:

```bash
pixi run -e bio bash scripts/make_sample_lists.sh <ENA report> <T2T panel BCF>
pixi run python scripts/make_samplesheet.py --ena <ENA report> --samples metadata/samples.tsv --out metadata/samplesheet-mage.csv
```

## 4. Check the wiring, then run

```bash
pixi run stub-full        # seconds: every step with stubs, on placeholder inputs
pixi run pipeline-full    # nextflow run workflow/nextflow/main.nf -profile hpc -c config/site.config --mode full -resume
```

Run the Nextflow head in `tmux` (or as a small, long Slurm job). `-resume` reuses every finished task; STAR references
and indexes are kept in `processed_dir` (storeDir) and never rebuilt.

Resources (our run, 731 individuals, 779 runs; cluster limit 128 cores, 400 GB): STAR 16 CPUs and 48 GB per task, about
8 minutes per run and reference; salmon 8 CPUs; at most 6 Slurm jobs at a time and 3 SRA downloads. In total about
3,000 CPU hours. Each run's FASTQ (~10 GB) exists only while its tasks run; genome BAMs are deleted after counting, and
only alignments in the candidate regions (`keep-regions.bed`) are kept.

## 5. Publishing the bundle

The run ends with `results/release/`: `mage-processed-v1/`, `mage-processed-v1.tar.gz`, its `.sha256`, and
`checksums.sha256`. To make it the input of `pixi run reproduce`:

1. Create the GitHub release `data-v1` and attach `mage-processed-v1.tar.gz`.
2. Put the archive's sha256 in `metadata/sources.tsv` (row `processed-mage`, column `published_checksum`,
   `sha256:<hash>`), and copy the lines of `checksums.sha256` into `metadata/checksums.sha256`.
3. Commit, push, and check that the *Reproduce* action is green.

A new bundle gets a new name and release (`mage-processed-v2`, `data-v2`), so earlier versions stay reproducible.
