#!/usr/bin/env bash
# Check the inputs of the cluster pipeline against metadata/checksums-cluster-inputs.sha256, or record them.
#
#   pixi run check-inputs [raw folder=data/raw]           check every listed file (prints failures, then a summary)
#   bash scripts/check_inputs.sh write [raw folder]       rewrite the list from the files present (after a download)
#
# Paths in the list start with data/raw/; a different raw folder (e.g. project storage on a cluster) replaces that
# prefix. Covered: genotype panels, T2T consortium BED files, genomes, annotations, liftOver chains.
set -euo pipefail
LIST=metadata/checksums-cluster-inputs.sha256
MODE=check
if [[ ${1:-} == write || ${1:-} == check ]]; then
  MODE=$1
  shift
fi
RAW=${1:-data/raw}

case $MODE in
  write)
    ( cd "$RAW" && find -L external references -type f -print0 | sort -z | xargs -0 -n 1 -P "$(nproc)" sha256sum ) \
      | sed 's#  #  data/raw/#' | sort -k 2 > "$LIST.tmp"
    mv "$LIST.tmp" "$LIST"
    echo "wrote $LIST: $(wc -l < "$LIST") files"
    ;;
  check)
    if sed "s#  data/raw/#  $RAW/#" "$LIST" | sha256sum -c --quiet -; then
      echo "all $(wc -l < "$LIST") inputs match"
    else
      echo "some inputs differ or are missing" >&2
      exit 1
    fi
    ;;
esac
