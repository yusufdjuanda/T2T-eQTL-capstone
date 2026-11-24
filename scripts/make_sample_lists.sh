#!/usr/bin/env bash
# Which MAGE individuals have genotypes in the phased 1KGP T2T-CHM13v2.0 panel, i.e. can be used for eQTL mapping.
#
#   pixi run -e bio bash scripts/make_sample_lists.sh <ENA run report> <T2T panel BCF> [out=metadata/samples.tsv]
#
# Inputs (metadata/sources.tsv): the ENA run report of PRJNA851328 (sample_alias = individual) and the panel BCF
# (sample names from its header; bcftools query -l reads the header only).
# Writes one row per MAGE individual:  sample  mage  t2t_panel  included   (included = in MAGE and in the panel).
# Read by scripts/make_samplesheet.py and by the pipeline (genotype subsets).
set -euo pipefail
MAGE=${1:?ENA run report}; PANEL=${2:?T2T panel BCF}; OUT=${3:-metadata/samples.tsv}
tmp=$(mktemp -d)
trap 'rm -rf "${tmp:?}"' EXIT

awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="sample_alias") c=i; next} {print $c}' "$MAGE" | sort -u > "$tmp/mage"
bcftools query -l "$PANEL" | sort -u > "$tmp/panel"
awk -v p="$tmp/panel" 'BEGIN{
      OFS = "\t"
      while ((getline x < p) > 0) P[x] = 1
      print "sample", "mage", "t2t_panel", "included" }
    function yn(v) { return v ? "yes" : "no" }
    { print $1, "yes", yn($1 in P), yn($1 in P) }' "$tmp/mage" > "$OUT"
awk -F'\t' 'NR>1{ m++; if ($4=="yes") n++ } END{ printf "MAGE %d individuals, %d in the panel (included)\n", m, n }' "$OUT"
