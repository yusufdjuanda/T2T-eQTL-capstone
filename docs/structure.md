# Project structure

The layout follows the course project template (cookiecutter-repro-project). Every folder has a README.

| Folder | What goes in it | Rule |
|--------|-----------------|------|
| `data/raw/` | downloaded inputs: the processed bundle (`pixi run download`); on a cluster, the genotype panels, genomes and annotations | read-only, not committed |
| `data/processed/` | derived inputs reused between runs: STAR references and indexes | rebuilt by the pipeline, not committed |
| `metadata/` | sources, checksums, sample list, samplesheet, chromosome names, data dictionary | committed |
| `config/` | `params.yaml` (analysis settings), `site.config.example` (cluster paths and queues) | committed; `site.config` is not |
| `workflow/nextflow/` | the pipeline: `main.nf`, `subworkflows/`, `conf/`, `modules/<module>/` (processes, the scripts they run, nf-tests), `tests/` | committed |
| `scripts/` | project scripts the pipeline does not run: downloads, metadata, lint | committed |
| `reports/` | Quarto report and bibliography | source committed, HTML not |
| `results/` | everything the pipeline produces | rebuilt, not committed |
| `tests/` | pytest unit tests (`python/`) and small test data (`data/`) | committed |
| `notebooks/` | exploration only | committed, outputs cleared |
| `docs/` | this file, the cluster guide, design decisions, software versions | committed |
| `envs/` | Conda export of the default environment; vendored tensorQTL | committed |
| `containers/` | Dockerfile and Apptainer recipe built from `pixi.lock` | committed |
| `.github/workflows/` | *Reproduce* (clean clone, `pixi run reproduce`) and *Tests* | committed |

## How the pipeline fits together

```text
                        --mode full (Slurm, -profile hpc)
  metadata/samplesheet-mage.csv ─┐
  data/raw/external (panels) ────┼─ GENOTYPES ─── subsets per chromosome ── variant pairs T2T <-> GRCh38 ─┐
  data/raw/references ───────────┼─ REFERENCES ── STAR references + indexes (data/processed) ──────────────┤
  NCBI SRA ──────────────────────┴─ EXPRESSION ── FASTQ -> STAR x2 -> salmon x2 -> matrices ──────────────┤
                                                                                                          ▼
                                    EQTL ── tensorQTL per reference x method (matched variants) ── COMPARE
                                                                                                          │
                                                       EXPORT_BUNDLE ── mage-processed-v1 (release asset) ┘
                                                                                │
  --mode reproduce (local)   data/raw/mage-processed-v1 ◄── pixi run download ──┘
                                        │
                                    ANALYSIS ── paired test ── gene sets ── figures ── reports/ (Quarto)
```

The same three ANALYSIS processes run in both modes, so `pixi run reproduce` repeats exactly the last part of the
cluster run.
