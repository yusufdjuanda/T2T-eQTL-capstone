# tests

| What | Where | Run |
|------|-------|-----|
| Unit tests of the scripts (pytest) | `tests/python/` | `pixi run test-python` |
| Tests of processes (nf-test) | `workflow/nextflow/modules/<module>/tests/` | `pixi run test-nf` |
| Tests of the ANALYSIS subworkflow and the whole pipeline (nf-test) | `workflow/nextflow/tests/` | `pixi run test-nf` |
| Stub run of the full cluster pipeline | profile `test` | `pixi run stub-full` |
| All of them | | `pixi run test` |

Test data (`tests/data/`, small, committed):

- `bundle/` — a synthetic bundle in the release format (80 individuals, 40 genes, 6 planted effect changes), made by
  `make_test_bundle.py`; the paired test must find exactly the planted genes.
- `expression/` — two tiny runs in STAR and salmon output format, for MERGE_EXPRESSION.
- `stub/` — empty placeholder inputs of the full pipeline, for the stub run.
