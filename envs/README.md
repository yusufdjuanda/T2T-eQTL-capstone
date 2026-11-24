# envs

- `environment.yml` — the default pixi environment exported for Conda users (`pixi run export-conda`):
  `conda env create -f envs/environment.yml`. Pixi users do not need it; `pixi.lock` is the exact record.
- `vendor/` — tensorQTL 1.0.10 and qtl 0.1.10, unpacked from their PyPI wheels (kept next to them, sha256 in each
  wheel's RECORD). The `qtl` environment puts them on `PYTHONPATH`; they were vendored because pixi could not reach its
  conda-PyPI mapping service from the cluster this was developed on.
