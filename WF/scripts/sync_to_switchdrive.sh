#!/bin/bash
# One-way sync of the shareable result tables to SWITCHdrive (rclone remote "switchdrive", see docs/colab_sharing.md).
# Run on the vscode node (compute nodes have no internet). Files deleted locally are also deleted on SWITCHdrive.
set -euo pipefail

TABLES=/scicore/home/zavolan/GROUP/pbody_RNA_spang_zavolan_v2/tables/gene_expression_analysis/bulkRNAseq/Sanity_utils
REMOTE=${REMOTE:-switchdrive:pbody_RNA_Spang_shared}
FILTER=$(dirname "$0")/sync_to_switchdrive.filter

module load rclone/1.68.2 2>/dev/null || true
rclone sync "$TABLES" "$REMOTE" --filter-from "$FILTER" --checksum --transfers 4 --log-level INFO "$@"
