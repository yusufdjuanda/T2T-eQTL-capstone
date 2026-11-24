# scripts

The code each pipeline step runs. Every script states its inputs, outputs and usage at the top, takes its paths as
arguments (no hard-coded locations), and can be run on its own. Python scripts are importable modules, so the tests in
`tests/python/` call their functions directly.

| Step | Scripts |
|------|---------|
| data and metadata | `download_data.py` (release asset), `fetch_t2t_panel.sh`, `fetch_grch38_panel.sh`, `check_inputs.sh`, `make_sample_lists.sh`, `make_samplesheet.py` |
| genotypes | `match_variants.py` |
| references | `prepare_star_reference.sh`, `keep_regions.py` |
| expression | `fetch_sra.sh`, `align_star.sh`, `merge_expression.py` |
| eQTL mapping | `prepare_genotypes.sh`, `prepare_phenotypes.py`, `make_covariates.py`, `run_tensorqtl.py`, `collect_eqtl.py` |
| comparison | `exon_overlap_pairs.sh`, `annotation_overlaps.py`, `gene_classes.py`, `compare_eqtl.py`, `compare_effects.py`, `export_bundle.py` |
| analysis | `paired_test.py`, `gene_sets.py`, `make_figures.R` |
| project | `lint.sh` (static checks), `install_fonts.sh` (optional figure font) |

Which pipeline process calls which script: `workflow/nextflow/modules/`.
