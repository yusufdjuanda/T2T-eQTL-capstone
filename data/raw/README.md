# data/raw

Input data exactly as received. **Never edit these files**; they are read-only after download.

- `mage-processed-v1/` — the processed MAGE inputs of the paired test, downloaded and checked by
  `pixi run download` (`scripts/download_data.py`) from the project's GitHub release `data-v1`.
  This is all that `pixi run reproduce` needs.
- `external/`, `references/`, `metadata/` — the inputs of the full cluster pipeline (genotype panels, genomes and
  annotations, ENA run report). About 60 GB; see `docs/cluster.md`. On the cluster they stay on project storage
  and `config/site.config` points the pipeline at them.

Where every file comes from: `metadata/sources.tsv`. Checksums: `metadata/checksums.sha256` (release asset) and
`metadata/checksums-cluster-inputs.sha256` (cluster inputs). Not committed to Git.
