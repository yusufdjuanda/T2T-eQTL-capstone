# tests/data/stub

Placeholder inputs for the stub run of the full pipeline (`pixi run stub-full`, profile `test`). Every file under
`raw/` is empty: in a stub run no tool reads them, but the pipeline checks that they exist, so the run proves that
every input path and every channel connection is right. The two runs in `samplesheet.csv` are the first two MAGE runs
of `metadata/samplesheet-mage.csv`; chromosomes chr21 and chr22 only.
