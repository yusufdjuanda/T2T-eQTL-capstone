# config

Settings of the pipeline. Change a setting here, not inside the code.

| File | What it holds |
|------|---------------|
| `params.yaml` | analysis settings behind the published results (bootstrap, effect floor, MAF, permutations, FDR, chromosomes); every pipeline task passes it with `-params-file` |
| `site.config.example` | template for `site.config`: where the cluster run finds its inputs and writes its outputs, and which Slurm queues to use |
| `site.config` | your copy, not committed (machine-specific paths) |

Nextflow profiles (`local`, `hpc`, `test`) are in `workflow/nextflow/nextflow.config`.
