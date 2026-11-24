# workflow/nextflow

The pipeline (Nextflow DSL2, strict syntax, Nextflow 26). Run it through the pixi tasks: `pixi run pipeline`
(reproduce mode, local), `pixi run pipeline-full` (full mode, Slurm), `pixi run stub-full` (wiring check).

```text
main.nf               entry point: params, --mode reproduce|full, workflow outputs
nextflow.config       profiles local, hpc, test
conf/                 base.config (software of the cluster steps), local/hpc/test.config
subworkflows/         how the processes connect: GENOTYPES, REFERENCES, EXPRESSION, EQTL, COMPARE, ANALYSIS, FULL
modules/<module>/     main.nf (the processes), scripts/ (the code they run), tests/ (nf-test, where present)
tests/                nf-test of the ANALYSIS subworkflow and of the whole pipeline in both modes
```

Every process calls at most one script of its own module (`${moduleDir}/scripts/...`); commands of one or a few lines
are written in the process itself. Every process has a `stub:` block for `-stub-run`.

| Module | Processes | Scripts |
|--------|-----------|---------|
| `genotypes` | SUBSET_GENOTYPES, MATCH_VARIANTS, VARIANT_MATCHING_SUMMARY | `match_variants.py` |
| `references` | PREPARE_STAR_REFERENCE, STAR_INDEX, KEEP_REGIONS | `prepare_star_reference.sh`, `keep_regions.py` |
| `expression` | FETCH_SRA, FASTERQ_DUMP, STAR_ALIGN, SALMON_QUANT, MERGE_EXPRESSION | `fetch_sra.sh`, `align_star.sh`, `merge_expression.py` |
| `eqtl` | PREPARE_GENOTYPES, GENOTYPE_PCS, PREPARE_PHENOTYPES, MAKE_COVARIATES, TENSORQTL_CIS, COLLECT_EQTL | `prepare_genotypes.sh`, `prepare_phenotypes.py`, `make_covariates.py`, `run_tensorqtl.py`, `collect_eqtl.py` |
| `compare` | ANNOTATION_OVERLAPS, GENE_CLASSES, COMPARE_EQTL, COMPARE_EFFECTS, EXPORT_BUNDLE | `exon_overlap_pairs.sh`, `annotation_overlaps.py`, `gene_classes.py`, `compare_eqtl.py`, `compare_effects.py`, `export_bundle.py` |
| `analysis` | PAIRED_TEST, GENE_SETS, FIGURES | `paired_test.py`, `gene_sets.py`, `make_figures.R` |

Each script states its inputs, outputs and usage at the top and can be run on its own. Python scripts are importable
modules; `tests/python/` tests their functions with pytest.
