#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
TRAIN_SCRIPT="${REPO_ROOT}/compared_models/scripts/train_nt_compared_model.py"
COLLECT_SCRIPT="${REPO_ROOT}/compared_models/scripts/collect_nt_grid_search_results.py"
LAUNCHER_NAME="$(basename "$0")"
LAUNCHER_MODEL_NAME="generator"
LAUNCHER_CONDA_ENV="visualdna"
LAUNCHER_DEFAULT_FINETUNE_METHOD="frozen_linear_probe"
LAUNCHER_DEFAULT_DEVICE="cuda"
LAUNCHER_DEFAULT_LOG_ROOT="/zengxiangxiang/mps/visualdna/outputs/visualdna_benchmark/nt_grid_generator"

source "${SCRIPT_DIR}/nt_compared_model_grid_search_lib.sh"
cd "$REPO_ROOT"
run_model_grid_search_launcher "$@"
