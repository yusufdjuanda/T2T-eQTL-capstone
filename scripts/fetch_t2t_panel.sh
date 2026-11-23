#!/usr/bin/env bash
# Download the phased 1000 Genomes T2T-CHM13v2.0 panel (SHAPEIT5 v1.1, 3,202 samples) from the human-pangenomics
# S3 bucket. About 30 GB; check it afterwards with `pixi run check-inputs`.
#
#   pixi run fetch-t2t-panel [destination=data/raw/external/1kgp-chm13v2.0-phased-shapeit5-v1.1]
#
# HTTPS to S3 from some clusters is reset during the TLS handshake on most attempts, so curl retries with a short
# fixed delay and resumes partial files (-C -). Plain HTTP is refused by the bucket. Rerunning is safe.
set -euo pipefail
dest=${1:-data/raw/external/1kgp-chm13v2.0-phased-shapeit5-v1.1}
base=https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1
bcf=1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz
mkdir -p "$dest"
cd "$dest"

fetch() {
  local f=$1
  [[ -f $f && ! -w $f ]] && { echo "$f already complete"; return 0; }
  for _ in $(seq 1 2000); do
    if curl -sSf -C - --connect-timeout 20 --retry 50 --retry-delay 1 --retry-all-errors -o "$f" "$base/$f"; then
      echo "$(date -Is) done $f ($(stat -c %s "$f") bytes)"
      chmod a-w "$f"
      return 0
    fi
    sleep 1
  done
  echo "$(date -Is) FAILED $f" >&2
  return 1
}
for f in README.txt "$bcf.csi" "$bcf"; do fetch "$f"; done
