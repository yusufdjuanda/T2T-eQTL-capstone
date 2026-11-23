#!/usr/bin/env bash
# Download one SRA run and validate it (called by the pipeline process FETCH_SRA).
#
#   bash scripts/fetch_sra.sh <run accession> <output folder>
#
# Writes <output folder>/<run>.sra after vdb-validate passes (NCBI checksums). prefetch can loop for hours without
# failing (2026-10-17: stuck at ~120 MB, restarting its own partial file), so each attempt has a time limit; after a
# timed-out attempt the partial files are removed and it starts again.
set -euo pipefail
RUN=${1:?run accession}
OUT=${2:?output folder}
mkdir -p "$OUT"
for i in $(seq 1 100); do
  rc=0
  timeout 30m prefetch --max-size 100G -O "$OUT" "$RUN" || rc=$?
  [[ $rc -eq 0 ]] && break
  echo "$(date -u +%T) $RUN prefetch attempt $i exit $rc" >&2
  if [[ $rc -eq 124 ]]; then
    rm -f -- "${OUT:?}/${RUN:?}/${RUN:?}.sra.tmp" "${OUT:?}/${RUN:?}/${RUN:?}.sra.prf" "${OUT:?}/${RUN:?}/${RUN:?}.sra.lock"
  fi
  if [[ $i -eq 100 ]]; then
    echo "$RUN: prefetch failed 100 times" >&2
    exit 1
  fi
  sleep 10
done
vdb-validate "$OUT/$RUN/$RUN.sra"
mv "$OUT/$RUN/$RUN.sra" "$OUT/$RUN.sra"
rmdir "$OUT/$RUN"
