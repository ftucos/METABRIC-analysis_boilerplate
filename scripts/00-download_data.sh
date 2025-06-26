#!/usr/bin/env bash
# Download specific sub-directories from cBioPortal/datahub
set -euo pipefail

OWNER="cBioPortal"
REPO="datahub"
REF="37f40e90f4b66e6ef48630129482f8863a4935c2"      # tag / branch / full SHA

BASE_DEST="../data/raw_data"                         # destination base directory

SUBDIRS=(
  "public/brca_metabric"
)

# prepare local destination dirs
for subdir in "${SUBDIRS[@]}"; do
  mkdir -p "${BASE_DEST}/$(basename "$subdir")"
done

WORK=$(mktemp -d)
# ensure WORK is removed on exit
trap 'rm -rf "$WORK"' EXIT

# print temporary directory
echo "Using temporary directory: $WORK"

# shallow clone with sparse-checkout pre-configured
git -C "$WORK" init -q
git -C "$WORK" remote add origin "https://github.com/${OWNER}/${REPO}.git"
git -C "$WORK" config core.sparseCheckout true
git -C "$WORK" lfs install # initialize Git LFS
printf '%s\n' "${SUBDIRS[@]}" > "$WORK/.git/info/sparse-checkout"
git -C "$WORK" pull --depth=1 origin "$REF"

# copy each subdir to BASE_DEST/$(basename)
for subdir in "${SUBDIRS[@]}"; do
  dest="${BASE_DEST}/$(basename "$subdir")"
  rsync -a --delete "$WORK/${subdir}/" "${dest}/"
  echo "✅  download complete: ${dest}"
done