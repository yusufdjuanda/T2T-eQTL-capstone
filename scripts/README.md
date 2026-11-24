# scripts

Project scripts that the pipeline does not run: getting the data, building the metadata, and project tools. Each
states its inputs, outputs and usage at the top, and has a pixi task where it is run often.

| Script | What it does | Run with |
|--------|--------------|----------|
| `download_data.py` | downloads the processed MAGE bundle (release asset) and checks every checksum | `pixi run download` (first step of `reproduce`) |
| `fetch_t2t_panel.sh` | downloads the 1000 Genomes T2T-CHM13 genotype panel (cluster run) | `pixi run fetch-t2t-panel` |
| `fetch_grch38_panel.sh` | downloads the 1000 Genomes GRCh38 genotype panel (cluster run) | `pixi run fetch-grch38-panel` |
| `check_inputs.sh` | checks the cluster inputs against `metadata/checksums-cluster-inputs.sha256` | `pixi run check-inputs` |
| `make_sample_lists.sh` | writes `metadata/samples.tsv` (MAGE individuals in the T2T panel) | see `docs/cluster.md` |
| `make_samplesheet.py` | writes `metadata/samplesheet-mage.csv` from the ENA run report | see `docs/cluster.md` |
| `lint.sh` | static checks of the whole project | `pixi run -e lint lint` |
| `install_fonts.sh` | installs the figure font (optional) | `pixi run fonts` |

The scripts each pipeline step runs live with their Nextflow module: `workflow/nextflow/modules/<module>/scripts/`
(see `workflow/nextflow/README.md`).
