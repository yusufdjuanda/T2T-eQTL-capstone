#!/usr/bin/env bash
# Step 1: which MAGE individuals have genotypes in the phased 1KGP T2T-CHM13v2.0 panel, i.e. can be used for eQTL
# mapping.
#
#   bash scripts/01-genotypes/make-sample-lists.sh
#
# Sources (in data/metadata/ and data/external/; see config/data-sources.tsv):
#   MAGE   ENA run report for PRJNA851328 (mage/ena-PRJNA851328-2026-10-17.tsv): sample_alias = individual
#   panel  sample names in the panel BCF header (3,202; bcftools query -l reads the header only)
# Writes config/samples/expression-individuals.tsv, one row per MAGE individual:
#   sample  mage  t2t_panel  included      (yes/no; included = in MAGE and in the panel)
# Read by scripts/04-expression/make-cohort-sheets.py (which individuals' runs go into the samplesheet) and by the
# genotype subset scripts (01-genotypes/chm13v2.0/subset-genotypes.sbatch, 01-genotypes/grch38/subset-genotypes.sbatch).

set -euo pipefail
ROOT=${T2T_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)}
source "$ROOT/scripts/lib/common.sh"
MAGE=$ROOT/data/metadata/mage/ena-PRJNA851328-2026-10-17.tsv
PANEL=$ROOT/data/external/1kgp-chm13v2.0-phased-shapeit5-v1.1/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz
OUT=$ROOT/config/samples/expression-individuals.tsv
t2t_env bio   # bcftools from pixi.toml
mkdir -p "$T2T_TMP"
tmp=$(mktemp -d -p "$T2T_TMP"); trap 'rm -rf $tmp' EXIT

# individual IDs from the sample_alias column, found by name
alias_col() { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="sample_alias") c=i; next} {print $c}' "$1" | sort -u; }
alias_col $MAGE > $tmp/mage
bcftools query -l $PANEL | sort -u > $tmp/panel

awk -v p=$tmp/panel 'BEGIN{
      OFS = "\t"
      while ((getline x < p) > 0) P[x] = 1
      print "sample", "mage", "t2t_panel", "included" }
    function yn(v) { return v ? "yes" : "no" }
    { print $1, "yes", yn($1 in P), yn($1 in P) }' $tmp/mage > $OUT

awk -F'\t' 'NR>1{ m++; if ($4=="yes") n++ } END{ printf "MAGE %d individuals, %d in the panel (included)\n", m, n }' $OUT
