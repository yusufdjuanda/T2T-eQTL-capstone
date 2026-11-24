# reports

`t2t-eqtl-report.qmd` is the Quarto report: question, data, methods, results and software versions, with every
number, table and figure computed from `results/analysis/` when it renders. `pixi run report` renders it to
`t2t-eqtl-report.html` (not committed; `pixi run reproduce` rebuilds it). References: `references.bib`.

To render it from the synthetic test bundle instead of the real data (for a quick check of the report itself):

    pixi run nextflow run workflow/nextflow/main.nf -profile local --bundle tests/data/bundle
    T2T_BUNDLE=tests/data/bundle pixi run quarto render reports/t2t-eqtl-report.qmd
