# Design decisions

Choices a reader might otherwise have to guess, and why they were made.

## Study design

- **Same reads, same people, two references.** Every difference between T2T-CHM13 and GRCh38 results then comes from
  the reference (sequence and annotation), not from sampling.
- **Same annotation release on both references** (NCBI RefSeq RS_2025_08); genes are matched by NCBI GeneID.
- **T2T made comparable to GRCh38_noalt**: chrY PARs hard-masked and chrM/chrEBV added, as in the GRCh38 analysis set
  (`workflow/nextflow/modules/references/scripts/prepare_star_reference.sh`).
- **Genotypes called on each reference**, never lifted over. liftOver is used only to pair variants called
  independently on each reference (`workflow/nextflow/modules/genotypes/scripts/match_variants.py`).
- **Effects compared at the same variant**, not by eGene counts: genes such as PWP2 and CBS are eGenes on both
  references, but their GRCh38 effect is inflated by the chr21 false duplication.
- **Allele coding from the genotype correlation** of the two variants, not from REF/ALT labels: REF and ALT differ
  between the references at about 20% of common SNPs, and some indels are written differently by the two panels.
- **Two quantification methods** (STAR gene counts: unique reads only; salmon: multi-mapping reads shared by EM). A
  change counts only if both agree, which removes many annotation and multi-mapping artefacts.

## Statistics

- **Paired test instead of independent z.** The two effect estimates share individuals and most of their noise
  (median residual correlation 0.99), so the independent z is about 9 times too conservative. The covariance of the two
  OLS slopes is computed exactly and checked by a bootstrap over individuals (`workflow/nextflow/modules/analysis/scripts/paired_test.py`, tests in
  `tests/python/test_paired_test.py`).
- **Bonferroni** over the genes compared, per method; effect floor of 0.1 SD per allele for the main set.
- **Artefact flags** (exon-overlap partners differ between annotations, MHC, PARs) remove genes whose change has a known
  non-biological cause.
- **Autosomes only**; chrX needs sex-aware genotype handling.

## Reproducibility

- **Two tiers.** The full run needs a cluster and 2.5 TB of reads, which no reviewer can repeat. `pixi run reproduce`
  starts from a 150 MB bundle of the paired test's inputs (release asset, checksummed), and reruns the analysis, the
  figures and the report with the same code the cluster run uses (subworkflow ANALYSIS in both modes).
- **Plain genotype dosages in the bundle** instead of pgen files, so the reproducible part needs no tensorQTL, PyTorch
  or plink2 and runs on macOS too.
- **One Nextflow pipeline** with `local`, `hpc` and `test` profiles; cluster paths only in `config/site.config`.
- **Stubs for every process**, so the wiring of the full pipeline is tested on every push without data or tools.
- **tensorQTL vendored** (envs/vendor/), because pixi could not reach its conda-PyPI mapping service from the cluster.
