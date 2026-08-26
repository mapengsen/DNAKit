#!/usr/bin/env bash

set -Eeuo pipefail

DEFAULT_LR_VALUES_RAW="1e-3 1e-4 1e-5"
DEFAULT_BATCH_SIZE_VALUES_RAW="16 32 64"
DEFAULT_WEIGHT_DECAY="1e-4"
DEFAULT_OPTIMIZER_NAME="adamw"
DEFAULT_MAX_EPOCHS="50"
DEFAULT_EARLY_STOPPING_PATIENCE="15"
DEFAULT_WARMUP_RATIO="0.1"
DEFAULT_MIN_LR_RATIO="0.1"
DEFAULT_RANDOM_SEED="42"
DEFAULT_RAW_DOWNLOAD="/zengxiangxiang/mps/visualdna/data/raw_download/nucleotide_transformer_downstream_tasks_revised"
DEFAULT_PROCESSED_DOWNLOAD="/zengxiangxiang/mps/ood_imageDNA/data/low_similarity_sequence_csv_original_parquet/nt"
DEFAULT_GPU_ID_LIST="0 1 2 3 4 5 6 7"
DEFAULT_PROGRESS_SECONDS="5"
DEFAULT_EMBEDDING_CACHE_DTYPE="float32"
DEFAULT_FEATURE_TRAIN_DTYPE="float32"

EMBEDDING_CACHE_EXTRACT_SCRIPT="${EMBEDDING_CACHE_EXTRACT_SCRIPT:-${REPO_ROOT}/compared_models/scripts/extract_nt_embedding_cache.py}"
EMBEDDING_CACHE_TRAIN_SCRIPT="${EMBEDDING_CACHE_TRAIN_SCRIPT:-${REPO_ROOT}/compared_models/scripts/train_nt_cached_linear_probe.py}"
EMBEDDING_CACHE_JOB_SCRIPT="${EMBEDDING_CACHE_JOB_SCRIPT:-${REPO_ROOT}/compared_models/scripts/bash/run_nt_cached_linear_probe_job.sh}"

require_model_launcher_config() {
  local var_name
  for var_name in LAUNCHER_MODEL_NAME LAUNCHER_CONDA_ENV; do
    if [[ -z "${!var_name:-}" ]]; then
      echo "[launcher] Missing required launcher variable: ${var_name}" >&2
      return 1
    fi
  done
}

print_model_launcher_usage() {
  cat <<EOF
Usage:
  bash ${LAUNCHER_NAME} [launcher args] [train args...]

Launcher args:
  --task-name VALUE             NT task name; supports all or comma-separated names
  --processed-download PATH     processed low-similarity NT root
  --raw-download PATH           raw download path recorded in run_config.json
  --log-root PATH               output root for this model grid search
  --gpu-id-list "0 1 2 3"        GPU ids; default: ${DEFAULT_GPU_ID_LIST}
  --lr-values "1e-3 1e-4"       learning-rate grid
  --batch-size-values "16 32"   batch-size grid
  --finetune-method VALUE       frozen_linear_probe / full / ia3
  --conda-env VALUE             conda env used by train jobs
  --device VALUE                train device passed to Python, default: cuda
  --num-workers INT             DataLoader workers, default: 4
  --chunk-forward-batch-size INT
  --embedding-cache-root PATH   cache root for frozen_linear_probe; default: <log-root>/embedding_cache
  --no-embedding-cache          disable cache and use original per-epoch backbone forward
  --precompute-embedding-cache-first
                              precompute one frozen_linear_probe cache per task before grid training
  --overwrite-embedding-cache   regenerate cache even if it already exists
  --embedding-cache-dtype VALUE float32 / float16 for saved features; default: ${DEFAULT_EMBEDDING_CACHE_DTYPE}
  --feature-train-dtype VALUE   float32 / float16 for cached linear training; default: ${DEFAULT_FEATURE_TRAIN_DTYPE}
  --extract-batch-size INT      DataLoader batch size for cache extraction; default: train batch size
  --poll-seconds INT            scheduler poll interval, default: 5
  --progress-seconds INT        launcher progress log interval, default: ${DEFAULT_PROGRESS_SECONDS}; 0 disables
  --skip-existing               skip jobs whose metrics.json has status ok
  --dry-run                     print commands without launching jobs
  --stop-on-failure             stop scheduling new jobs after the first failure
  --no-collect                  do not collect metrics/profile/runtime into CSV/JSON
  --help                        show this help

Fixed defaults:
  model-name = ${LAUNCHER_MODEL_NAME}
  conda-env = ${LAUNCHER_CONDA_ENV}
  lr = ${DEFAULT_LR_VALUES_RAW}
  batch size = ${DEFAULT_BATCH_SIZE_VALUES_RAW}
  weight decay = ${DEFAULT_WEIGHT_DECAY}
  optimizer = ${DEFAULT_OPTIMIZER_NAME}
  epochs = ${DEFAULT_MAX_EPOCHS}
  early stopping patience = ${DEFAULT_EARLY_STOPPING_PATIENCE}
  warmup ratio = ${DEFAULT_WARMUP_RATIO}
  min lr ratio = ${DEFAULT_MIN_LR_RATIO}
  seed = ${DEFAULT_RANDOM_SEED}

Unrecognized args are forwarded to:
  python compared_models/scripts/train_nt_compared_model.py
EOF
}

parse_model_launcher_args() {
  TASK_NAME="${TASK_NAME:-all}"
  RAW_DOWNLOAD="${RAW_DOWNLOAD:-$DEFAULT_RAW_DOWNLOAD}"
  PROCESSED_DOWNLOAD="${PROCESSED_DOWNLOAD:-$DEFAULT_PROCESSED_DOWNLOAD}"
  LOG_ROOT="${LOG_ROOT:-${LAUNCHER_DEFAULT_LOG_ROOT:-/zengxiangxiang/mps/visualdna/outputs/visualdna_benchmark/nt_grid_${LAUNCHER_MODEL_NAME}}}"
  GPU_ID_LIST="${GPU_ID_LIST:-$DEFAULT_GPU_ID_LIST}"
  LR_VALUES_RAW="${LR_VALUES_RAW:-$DEFAULT_LR_VALUES_RAW}"
  BATCH_SIZE_VALUES_RAW="${BATCH_SIZE_VALUES_RAW:-$DEFAULT_BATCH_SIZE_VALUES_RAW}"
  WEIGHT_DECAY="${WEIGHT_DECAY:-$DEFAULT_WEIGHT_DECAY}"
  OPTIMIZER_NAME="${OPTIMIZER_NAME:-$DEFAULT_OPTIMIZER_NAME}"
  MAX_EPOCHS="${MAX_EPOCHS:-$DEFAULT_MAX_EPOCHS}"
  EARLY_STOPPING_PATIENCE="${EARLY_STOPPING_PATIENCE:-$DEFAULT_EARLY_STOPPING_PATIENCE}"
  WARMUP_RATIO="${WARMUP_RATIO:-$DEFAULT_WARMUP_RATIO}"
  MIN_LR_RATIO="${MIN_LR_RATIO:-$DEFAULT_MIN_LR_RATIO}"
  RANDOM_SEED="${RANDOM_SEED:-$DEFAULT_RANDOM_SEED}"
  FINETUNE_METHOD="${FINETUNE_METHOD:-${LAUNCHER_DEFAULT_FINETUNE_METHOD:-frozen_linear_probe}}"
  CONDA_ENV="${CONDA_ENV:-$LAUNCHER_CONDA_ENV}"
  TASK_RESOLVE_CONDA_ENV="${TASK_RESOLVE_CONDA_ENV:-visualdna}"
  COLLECT_CONDA_ENV="${COLLECT_CONDA_ENV:-visualdna}"
  DEVICE="${DEVICE:-${LAUNCHER_DEFAULT_DEVICE:-cuda}}"
  TOKEN_READOUT="${TOKEN_READOUT:-auto}"
  MONITOR="${MONITOR:-accuracy}"
  NUM_WORKERS="${NUM_WORKERS:-4}"
  CHUNK_FORWARD_BATCH_SIZE="${CHUNK_FORWARD_BATCH_SIZE:-8}"
  USE_EMBEDDING_CACHE="${USE_EMBEDDING_CACHE:-1}"
  EMBEDDING_CACHE_ROOT="${EMBEDDING_CACHE_ROOT:-}"
  PRECOMPUTE_EMBEDDING_CACHE_FIRST="${PRECOMPUTE_EMBEDDING_CACHE_FIRST:-0}"
  OVERWRITE_EMBEDDING_CACHE=0
  EMBEDDING_CACHE_DTYPE="${EMBEDDING_CACHE_DTYPE:-$DEFAULT_EMBEDDING_CACHE_DTYPE}"
  FEATURE_TRAIN_DTYPE="${FEATURE_TRAIN_DTYPE:-$DEFAULT_FEATURE_TRAIN_DTYPE}"
  EXTRACT_BATCH_SIZE="${EXTRACT_BATCH_SIZE:-}"
  POLL_SECONDS="${POLL_SECONDS:-5}"
  PROGRESS_SECONDS="${PROGRESS_SECONDS:-$DEFAULT_PROGRESS_SECONDS}"
  TENSORBOARD=1
  ALLOW_REMOTE_FALLBACK=1
  SKIP_EXISTING=0
  DRY_RUN=0
  STOP_ON_FAILURE=0
  COLLECT_RESULTS=1
  EXTRA_ARGS=()

  while (( $# > 0 )); do
    case "$1" in
      --task-name|--task_name)
        TASK_NAME="$2"
        shift 2
        ;;
      --task-name=*|--task_name=*)
        TASK_NAME="${1#*=}"
        shift
        ;;
      --raw-download|--raw_download)
        RAW_DOWNLOAD="$2"
        shift 2
        ;;
      --raw-download=*|--raw_download=*)
        RAW_DOWNLOAD="${1#*=}"
        shift
        ;;
      --processed-download|--processed_download)
        PROCESSED_DOWNLOAD="$2"
        shift 2
        ;;
      --processed-download=*|--processed_download=*)
        PROCESSED_DOWNLOAD="${1#*=}"
        shift
        ;;
      --log-root)
        LOG_ROOT="$2"
        shift 2
        ;;
      --log-root=*)
        LOG_ROOT="${1#*=}"
        shift
        ;;
      --gpu-id-list)
        GPU_ID_LIST="$2"
        shift 2
        ;;
      --gpu-id-list=*)
        GPU_ID_LIST="${1#*=}"
        shift
        ;;
      --lr-values)
        LR_VALUES_RAW="$2"
        shift 2
        ;;
      --lr-values=*)
        LR_VALUES_RAW="${1#*=}"
        shift
        ;;
      --batch-size-values)
        BATCH_SIZE_VALUES_RAW="$2"
        shift 2
        ;;
      --batch-size-values=*)
        BATCH_SIZE_VALUES_RAW="${1#*=}"
        shift
        ;;
      --weight-decay)
        WEIGHT_DECAY="$2"
        shift 2
        ;;
      --weight-decay=*)
        WEIGHT_DECAY="${1#*=}"
        shift
        ;;
      --optimizer-name)
        OPTIMIZER_NAME="$2"
        shift 2
        ;;
      --optimizer-name=*)
        OPTIMIZER_NAME="${1#*=}"
        shift
        ;;
      --epochs)
        MAX_EPOCHS="$2"
        shift 2
        ;;
      --epochs=*)
        MAX_EPOCHS="${1#*=}"
        shift
        ;;
      --early-stopping-patience)
        EARLY_STOPPING_PATIENCE="$2"
        shift 2
        ;;
      --early-stopping-patience=*)
        EARLY_STOPPING_PATIENCE="${1#*=}"
        shift
        ;;
      --warmup-ratio)
        WARMUP_RATIO="$2"
        shift 2
        ;;
      --warmup-ratio=*)
        WARMUP_RATIO="${1#*=}"
        shift
        ;;
      --min-lr-ratio)
        MIN_LR_RATIO="$2"
        shift 2
        ;;
      --min-lr-ratio=*)
        MIN_LR_RATIO="${1#*=}"
        shift
        ;;
      --seed)
        RANDOM_SEED="$2"
        shift 2
        ;;
      --seed=*)
        RANDOM_SEED="${1#*=}"
        shift
        ;;
      --finetune-method)
        FINETUNE_METHOD="$2"
        shift 2
        ;;
      --finetune-method=*)
        FINETUNE_METHOD="${1#*=}"
        shift
        ;;
      --conda-env)
        CONDA_ENV="$2"
        shift 2
        ;;
      --conda-env=*)
        CONDA_ENV="${1#*=}"
        shift
        ;;
      --task-resolve-conda-env)
        TASK_RESOLVE_CONDA_ENV="$2"
        shift 2
        ;;
      --task-resolve-conda-env=*)
        TASK_RESOLVE_CONDA_ENV="${1#*=}"
        shift
        ;;
      --collect-conda-env)
        COLLECT_CONDA_ENV="$2"
        shift 2
        ;;
      --collect-conda-env=*)
        COLLECT_CONDA_ENV="${1#*=}"
        shift
        ;;
      --device)
        DEVICE="$2"
        shift 2
        ;;
      --device=*)
        DEVICE="${1#*=}"
        shift
        ;;
      --token-readout)
        TOKEN_READOUT="$2"
        shift 2
        ;;
      --token-readout=*)
        TOKEN_READOUT="${1#*=}"
        shift
        ;;
      --monitor)
        MONITOR="$2"
        shift 2
        ;;
      --monitor=*)
        MONITOR="${1#*=}"
        shift
        ;;
      --num-workers)
        NUM_WORKERS="$2"
        shift 2
        ;;
      --num-workers=*)
        NUM_WORKERS="${1#*=}"
        shift
        ;;
      --chunk-forward-batch-size)
        CHUNK_FORWARD_BATCH_SIZE="$2"
        shift 2
        ;;
      --chunk-forward-batch-size=*)
        CHUNK_FORWARD_BATCH_SIZE="${1#*=}"
        shift
        ;;
      --embedding-cache-root)
        EMBEDDING_CACHE_ROOT="$2"
        shift 2
        ;;
      --embedding-cache-root=*)
        EMBEDDING_CACHE_ROOT="${1#*=}"
        shift
        ;;
      --use-embedding-cache)
        USE_EMBEDDING_CACHE=1
        shift
        ;;
      --no-embedding-cache)
        USE_EMBEDDING_CACHE=0
        shift
        ;;
      --precompute-embedding-cache-first)
        PRECOMPUTE_EMBEDDING_CACHE_FIRST=1
        shift
        ;;
      --no-precompute-embedding-cache-first)
        PRECOMPUTE_EMBEDDING_CACHE_FIRST=0
        shift
        ;;
      --overwrite-embedding-cache)
        OVERWRITE_EMBEDDING_CACHE=1
        shift
        ;;
      --embedding-cache-dtype)
        EMBEDDING_CACHE_DTYPE="$2"
        shift 2
        ;;
      --embedding-cache-dtype=*)
        EMBEDDING_CACHE_DTYPE="${1#*=}"
        shift
        ;;
      --feature-train-dtype)
        FEATURE_TRAIN_DTYPE="$2"
        shift 2
        ;;
      --feature-train-dtype=*)
        FEATURE_TRAIN_DTYPE="${1#*=}"
        shift
        ;;
      --extract-batch-size)
        EXTRACT_BATCH_SIZE="$2"
        shift 2
        ;;
      --extract-batch-size=*)
        EXTRACT_BATCH_SIZE="${1#*=}"
        shift
        ;;
      --poll-seconds)
        POLL_SECONDS="$2"
        shift 2
        ;;
      --poll-seconds=*)
        POLL_SECONDS="${1#*=}"
        shift
        ;;
      --progress-seconds)
        PROGRESS_SECONDS="$2"
        shift 2
        ;;
      --progress-seconds=*)
        PROGRESS_SECONDS="${1#*=}"
        shift
        ;;
      --tensorboard)
        TENSORBOARD=1
        shift
        ;;
      --no-tensorboard)
        TENSORBOARD=0
        shift
        ;;
      --allow-remote-fallback)
        ALLOW_REMOTE_FALLBACK=1
        shift
        ;;
      --no-allow-remote-fallback)
        ALLOW_REMOTE_FALLBACK=0
        shift
        ;;
      --skip-existing)
        SKIP_EXISTING=1
        shift
        ;;
      --dry-run)
        DRY_RUN=1
        shift
        ;;
      --stop-on-failure)
        STOP_ON_FAILURE=1
        shift
        ;;
      --no-collect)
        COLLECT_RESULTS=0
        shift
        ;;
      --help|-h)
        print_model_launcher_usage
        exit 0
        ;;
      *)
        EXTRA_ARGS+=("$1")
        shift
        ;;
    esac
  done

  GPU_ID_LIST="${GPU_ID_LIST//,/ }"
  LR_VALUES_RAW="${LR_VALUES_RAW//,/ }"
  BATCH_SIZE_VALUES_RAW="${BATCH_SIZE_VALUES_RAW//,/ }"
  if [[ -z "$EMBEDDING_CACHE_ROOT" ]]; then
    EMBEDDING_CACHE_ROOT="${LOG_ROOT}/embedding_cache"
  fi
  read -r -a GPU_IDS <<<"$GPU_ID_LIST"
  read -r -a LR_VALUES <<<"$LR_VALUES_RAW"
  read -r -a BATCH_SIZE_VALUES <<<"$BATCH_SIZE_VALUES_RAW"

  if (( ${#GPU_IDS[@]} == 0 )); then
    echo "[launcher] At least one GPU id is required." >&2
    return 1
  fi
  if (( ${#LR_VALUES[@]} == 0 || ${#BATCH_SIZE_VALUES[@]} == 0 )); then
    echo "[launcher] lr and batch-size grids must be non-empty." >&2
    return 1
  fi
  if ! [[ "$POLL_SECONDS" =~ ^[1-9][0-9]*$ ]]; then
    echo "[launcher] --poll-seconds must be a positive integer: ${POLL_SECONDS}" >&2
    return 1
  fi
  if ! [[ "$PROGRESS_SECONDS" =~ ^[0-9]+$ ]]; then
    echo "[launcher] --progress-seconds must be a non-negative integer: ${PROGRESS_SECONDS}" >&2
    return 1
  fi
  if [[ "$EMBEDDING_CACHE_DTYPE" != "float32" && "$EMBEDDING_CACHE_DTYPE" != "float16" ]]; then
    echo "[launcher] --embedding-cache-dtype must be float32 or float16: ${EMBEDDING_CACHE_DTYPE}" >&2
    return 1
  fi
  if [[ "$FEATURE_TRAIN_DTYPE" != "float32" && "$FEATURE_TRAIN_DTYPE" != "float16" ]]; then
    echo "[launcher] --feature-train-dtype must be float32 or float16: ${FEATURE_TRAIN_DTYPE}" >&2
    return 1
  fi
  if (( PRECOMPUTE_EMBEDDING_CACHE_FIRST )) && { (( ! USE_EMBEDDING_CACHE )) || [[ "$FINETUNE_METHOD" != "frozen_linear_probe" ]]; }; then
    echo "[launcher] --precompute-embedding-cache-first requires frozen_linear_probe with embedding cache enabled." >&2
    return 1
  fi
}

validate_passthrough_args() {
  local arg
  for arg in "${EXTRA_ARGS[@]}"; do
    case "$arg" in
      --task-name|--task-name=*|--task_name|--task_name=*|--processed-download|--processed-download=*|--processed_download|--processed_download=*|--raw-download|--raw-download=*|--raw_download|--raw_download=*|--log-root|--log-root=*|--model-name|--model-name=*|--lr|--lr=*|--batch-size|--batch-size=*|--epochs|--epochs=*|--weight-decay|--weight-decay=*|--optimizer-name|--optimizer-name=*|--warmup-ratio|--warmup-ratio=*|--min-lr-ratio|--min-lr-ratio=*|--early-stopping-patience|--early-stopping-patience=*|--seed|--seed=*|--device|--device=*)
        echo "[launcher] Do not pass controlled train arg through EXTRA_ARGS: ${arg}" >&2
        return 1
        ;;
    esac
  done
}

ensure_processed_download_ready() {
  if [[ -d "$PROCESSED_DOWNLOAD" ]]; then
    return 0
  fi
  cat >&2 <<EOF
[launcher] processed-download does not exist: ${PROCESSED_DOWNLOAD}
[launcher] Materialize NT first, for example:
  conda run -n visualdna python ood_imageDNA/scripts/data_prep/materialize_raw_nt_low_similarity_split_to_visualdna_raw.py \\
    --nt-root /zengxiangxiang/mps/visualdna/data/raw_download/nucleotide_transformer_downstream_tasks_revised \\
    --nt-low-split-root ood_imageDNA/data/all_nt_original_no_dedup_dashing_k6/low_similarity_split \\
    --output-root ood_imageDNA/data/low_similarity_sequence_csv_original_parquet/nt \\
    --output-format parquet \\
    --parquet-compression zstd \\
    --overwrite
EOF
  return 1
}

resolve_all_processed_tasks() {
  local task_dir
  find "$PROCESSED_DOWNLOAD" -mindepth 1 -maxdepth 1 -type d -print \
    | sort \
    | while IFS= read -r task_dir; do
        if [[ -f "${task_dir}/meta.json" || -f "${task_dir}/raw/meta.json" ]]; then
          basename "$task_dir"
        fi
      done
}

resolve_requested_tasks() {
  mapfile -t ALL_SUPPORTED_TASKS < <(resolve_all_processed_tasks)
  if (( ${#ALL_SUPPORTED_TASKS[@]} == 0 )); then
    echo "[launcher] No processed NT tasks found under ${PROCESSED_DOWNLOAD}." >&2
    return 1
  fi

  if [[ "$TASK_NAME" == "all" ]]; then
    RESOLVED_TASKS=("${ALL_SUPPORTED_TASKS[@]}")
    return 0
  fi

  local normalized_task_names="${TASK_NAME// /}"
  IFS=',' read -r -a RESOLVED_TASKS <<<"$normalized_task_names"
  if (( ${#RESOLVED_TASKS[@]} == 0 )); then
    echo "[launcher] --task-name cannot be empty." >&2
    return 1
  fi

  local requested_task
  local supported_task
  local matched
  for requested_task in "${RESOLVED_TASKS[@]}"; do
    matched=0
    for supported_task in "${ALL_SUPPORTED_TASKS[@]}"; do
      if [[ "$requested_task" == "$supported_task" ]]; then
        matched=1
        break
      fi
    done
    if (( ! matched )); then
      echo "[launcher] Unknown NT task: ${requested_task}" >&2
      echo "[launcher] Available tasks: ${ALL_SUPPORTED_TASKS[*]}" >&2
      return 1
    fi
  done
}

job_task_name() {
  printf '%s\n' "${1%%::*}"
}

job_lr() {
  local remainder="${1#*::}"
  printf '%s\n' "${remainder%%::*}"
}

job_batch_size() {
  local remainder="${1#*::}"
  printf '%s\n' "${remainder#*::}"
}

sanitize_job_value() {
  local value="$1"
  value="${value//./p}"
  value="${value//-/_}"
  printf '%s\n' "$value"
}

should_precompute_embedding_cache_first() {
  (( PRECOMPUTE_EMBEDDING_CACHE_FIRST )) && (( USE_EMBEDDING_CACHE )) && [[ "$FINETUNE_METHOD" == "frozen_linear_probe" ]]
}

embedding_cache_path_for_task() {
  local task_name="$1"
  printf '%s/%s/%s/features.pt\n' "$EMBEDDING_CACHE_ROOT" "$LAUNCHER_MODEL_NAME" "$task_name"
}

build_job_queue() {
  JOB_QUEUE=()
  local task
  local lr
  local batch_size
  if (( USE_EMBEDDING_CACHE )) && [[ "$FINETUNE_METHOD" == "frozen_linear_probe" ]]; then
    for lr in "${LR_VALUES[@]}"; do
      for batch_size in "${BATCH_SIZE_VALUES[@]}"; do
        for task in "${RESOLVED_TASKS[@]}"; do
          JOB_QUEUE+=("${task}::${lr}::${batch_size}")
        done
      done
    done
    return 0
  fi

  for task in "${RESOLVED_TASKS[@]}"; do
    for lr in "${LR_VALUES[@]}"; do
      for batch_size in "${BATCH_SIZE_VALUES[@]}"; do
        JOB_QUEUE+=("${task}::${lr}::${batch_size}")
      done
    done
  done
}

build_train_command() {
  local task_name="$1"
  local lr="$2"
  local batch_size="$3"
  local task_log_dir="$4"
  local cache_path
  local effective_extract_batch_size

  if (( USE_EMBEDDING_CACHE )) && [[ "$FINETUNE_METHOD" == "frozen_linear_probe" ]]; then
    if (( ${#EXTRA_ARGS[@]} > 0 )); then
      echo "[launcher] embedding cache mode does not support extra passthrough train args yet: ${EXTRA_ARGS[*]}" >&2
      return 1
    fi
    cache_path="${EMBEDDING_CACHE_ROOT}/${LAUNCHER_MODEL_NAME}/${task_name}/features.pt"
    effective_extract_batch_size="${EXTRACT_BATCH_SIZE:-$batch_size}"
    TRAIN_CMD=(
      bash "$EMBEDDING_CACHE_JOB_SCRIPT"
      --conda-env "$CONDA_ENV"
      --extract-script "$EMBEDDING_CACHE_EXTRACT_SCRIPT"
      --train-script "$EMBEDDING_CACHE_TRAIN_SCRIPT"
      --task-name "$task_name"
      --processed-download "$PROCESSED_DOWNLOAD"
      --raw-download "$RAW_DOWNLOAD"
      --log-root "$task_log_dir"
      --cache-path "$cache_path"
      --model-name "$LAUNCHER_MODEL_NAME"
      --token-readout "$TOKEN_READOUT"
      --device "$DEVICE"
      --epochs "$MAX_EPOCHS"
      --batch-size "$batch_size"
      --extract-batch-size "$effective_extract_batch_size"
      --num-workers "$NUM_WORKERS"
      --lr "$lr"
      --weight-decay "$WEIGHT_DECAY"
      --optimizer-name "$OPTIMIZER_NAME"
      --warmup-ratio "$WARMUP_RATIO"
      --min-lr-ratio "$MIN_LR_RATIO"
      --early-stopping-patience "$EARLY_STOPPING_PATIENCE"
      --seed "$RANDOM_SEED"
      --monitor "$MONITOR"
      --chunk-forward-batch-size "$CHUNK_FORWARD_BATCH_SIZE"
      --feature-dtype "$EMBEDDING_CACHE_DTYPE"
      --feature-train-dtype "$FEATURE_TRAIN_DTYPE"
    )
    if (( TENSORBOARD )); then
      TRAIN_CMD+=(--tensorboard)
    else
      TRAIN_CMD+=(--no-tensorboard)
    fi
    if (( ALLOW_REMOTE_FALLBACK )); then
      TRAIN_CMD+=(--allow-remote-fallback)
    else
      TRAIN_CMD+=(--no-allow-remote-fallback)
    fi
    if (( OVERWRITE_EMBEDDING_CACHE )); then
      TRAIN_CMD+=(--overwrite-cache)
    fi
    if should_precompute_embedding_cache_first; then
      TRAIN_CMD+=(--skip-cache-extract)
    fi
    return 0
  fi

  TRAIN_CMD=(
    conda run --no-capture-output -n "$CONDA_ENV" python -u "$TRAIN_SCRIPT"
    --task-name "$task_name"
    --processed-download "$PROCESSED_DOWNLOAD"
    --raw-download "$RAW_DOWNLOAD"
    --log-root "$task_log_dir"
    --model-name "$LAUNCHER_MODEL_NAME"
    --finetune-method "$FINETUNE_METHOD"
    --token-readout "$TOKEN_READOUT"
    --device "$DEVICE"
    --epochs "$MAX_EPOCHS"
    --batch-size "$batch_size"
    --num-workers "$NUM_WORKERS"
    --lr "$lr"
    --weight-decay "$WEIGHT_DECAY"
    --optimizer-name "$OPTIMIZER_NAME"
    --warmup-ratio "$WARMUP_RATIO"
    --min-lr-ratio "$MIN_LR_RATIO"
    --early-stopping-patience "$EARLY_STOPPING_PATIENCE"
    --seed "$RANDOM_SEED"
    --monitor "$MONITOR"
    --chunk-forward-batch-size "$CHUNK_FORWARD_BATCH_SIZE"
  )
  if (( TENSORBOARD )); then
    TRAIN_CMD+=(--tensorboard)
  else
    TRAIN_CMD+=(--no-tensorboard)
  fi
  if (( ALLOW_REMOTE_FALLBACK )); then
    TRAIN_CMD+=(--allow-remote-fallback)
  else
    TRAIN_CMD+=(--no-allow-remote-fallback)
  fi
  TRAIN_CMD+=("${EXTRA_ARGS[@]}")
}

print_shell_command() {
  local gpu="$1"
  shift
  printf 'CUDA_VISIBLE_DEVICES=%q ' "$gpu"
  printf '%q ' "$@"
  printf '\n'
}

prepare_child_python_startup_env() {
  local startup_dir
  if [[ "${LAUNCHER_ENABLE_NUMPY_BLAS_FPE_PATCH:-0}" != "1" ]]; then
    return 0
  fi

  startup_dir="${LAUNCHER_PYTHON_STARTUP_DIR:-${REPO_ROOT}/compared_models/scripts/python_startup}"
  export VISUALDNA_SKIP_NUMPY_BLAS_FPE_CHECK="${VISUALDNA_SKIP_NUMPY_BLAS_FPE_CHECK:-1}"
  case ":${PYTHONPATH:-}:" in
    *":${startup_dir}:"*) ;;
    *) export PYTHONPATH="${startup_dir}${PYTHONPATH:+:${PYTHONPATH}}" ;;
  esac
}

should_skip_job() {
  local task_log_dir="$1"
  if (( ! SKIP_EXISTING )); then
    return 1
  fi
  [[ -f "${task_log_dir}/metrics.json" ]] && grep -q '"status": "ok"' "${task_log_dir}/metrics.json"
}

launch_job() {
  local job_key="$1"
  local gpu="$2"
  local task_name
  local lr
  local batch_size
  local lr_safe
  local batch_safe
  local job_dir_name
  local task_log_dir
  local launch_log_path

  task_name="$(job_task_name "$job_key")"
  lr="$(job_lr "$job_key")"
  batch_size="$(job_batch_size "$job_key")"
  lr_safe="$(sanitize_job_value "$lr")"
  batch_safe="$(sanitize_job_value "$batch_size")"
  job_dir_name="lr_${lr_safe}_bs_${batch_safe}"
  task_log_dir="${LOG_ROOT}/${task_name}/${job_dir_name}"
  launch_log_path="${LAUNCH_LOG_DIR}/${task_name}_${job_dir_name}_gpu_${gpu}.log"

  if should_skip_job "$task_log_dir"; then
    echo "[launcher] Skip existing ok job: ${task_name} / ${job_dir_name}"
    COMPLETED_JOBS+=("${job_key}:0:skipped")
    return 2
  fi

  mkdir -p "$task_log_dir" "$LAUNCH_LOG_DIR"
  build_train_command "$task_name" "$lr" "$batch_size" "$task_log_dir"

  echo "[launcher] Launch ${LAUNCHER_MODEL_NAME} ${task_name} / ${job_dir_name} -> GPU ${gpu}"
  if (( DRY_RUN )); then
    print_shell_command "$gpu" "${TRAIN_CMD[@]}"
    COMPLETED_JOBS+=("${job_key}:0:dry_run")
    return 2
  fi

  (
    cd "$REPO_ROOT"
    prepare_child_python_startup_env
    CUDA_VISIBLE_DEVICES="$gpu" "${TRAIN_CMD[@]}"
  ) >"$launch_log_path" 2>&1 &

  local pid=$!
  GPU_BUSY["$gpu"]="$pid"
  PID_TO_GPU["$pid"]="$gpu"
  PID_TO_JOB["$pid"]="$job_key"
  PID_TO_LOG["$pid"]="$launch_log_path"
  echo "[launcher] PID=${pid}, log=${launch_log_path}"
  return 0
}

build_cache_queue() {
  CACHE_QUEUE=("${RESOLVED_TASKS[@]}")
}

should_skip_cache_task() {
  local task_name="$1"
  local cache_path
  if (( OVERWRITE_EMBEDDING_CACHE )); then
    return 1
  fi
  cache_path="$(embedding_cache_path_for_task "$task_name")"
  [[ -f "$cache_path" ]]
}

build_cache_command() {
  local task_name="$1"
  local cache_path="$2"
  local effective_extract_batch_size
  effective_extract_batch_size="${EXTRACT_BATCH_SIZE:-${BATCH_SIZE_VALUES[0]}}"

  CACHE_CMD=(
    conda run --no-capture-output -n "$CONDA_ENV" python -u "$EMBEDDING_CACHE_EXTRACT_SCRIPT"
    --task-name "$task_name"
    --processed-download "$PROCESSED_DOWNLOAD"
    --cache-path "$cache_path"
    --model-name "$LAUNCHER_MODEL_NAME"
    --raw-download "$RAW_DOWNLOAD"
    --token-readout "$TOKEN_READOUT"
    --finetune-method frozen_linear_probe
    --device "$DEVICE"
    --batch-size "$effective_extract_batch_size"
    --num-workers "$NUM_WORKERS"
    --seed "$RANDOM_SEED"
    --chunk-forward-batch-size "$CHUNK_FORWARD_BATCH_SIZE"
    --feature-dtype "$EMBEDDING_CACHE_DTYPE"
    --skip-existing
  )
  if (( ALLOW_REMOTE_FALLBACK )); then
    CACHE_CMD+=(--allow-remote-fallback)
  else
    CACHE_CMD+=(--no-allow-remote-fallback)
  fi
  if (( OVERWRITE_EMBEDDING_CACHE )); then
    CACHE_CMD+=(--overwrite)
  fi
}

launch_cache_task() {
  local task_name="$1"
  local gpu="$2"
  local cache_path
  local launch_log_path

  cache_path="$(embedding_cache_path_for_task "$task_name")"
  launch_log_path="${LAUNCH_LOG_DIR}/cache_${task_name}_gpu_${gpu}.log"

  if should_skip_cache_task "$task_name"; then
    echo "[launcher] Skip existing cache: ${task_name}"
    CACHE_COMPLETED_TASKS+=("${task_name}:0:skipped")
    return 2
  fi

  mkdir -p "$(dirname "$cache_path")" "$LAUNCH_LOG_DIR"
  build_cache_command "$task_name" "$cache_path"

  echo "[launcher] Precompute cache ${LAUNCHER_MODEL_NAME} ${task_name} -> GPU ${gpu}"
  if (( DRY_RUN )); then
    print_shell_command "$gpu" "${CACHE_CMD[@]}"
    CACHE_COMPLETED_TASKS+=("${task_name}:0:dry_run")
    return 2
  fi

  (
    cd "$REPO_ROOT"
    prepare_child_python_startup_env
    CUDA_VISIBLE_DEVICES="$gpu" "${CACHE_CMD[@]}"
  ) >"$launch_log_path" 2>&1 &

  local pid=$!
  GPU_BUSY["$gpu"]="$pid"
  PID_TO_GPU["$pid"]="$gpu"
  PID_TO_JOB["$pid"]="$task_name"
  PID_TO_LOG["$pid"]="$launch_log_path"
  echo "[launcher] CACHE_PID=${pid}, log=${launch_log_path}"
  return 0
}

print_cache_progress() {
  if (( PROGRESS_SECONDS <= 0 )); then
    return 0
  fi

  local completed_count="${#CACHE_COMPLETED_TASKS[@]}"
  local failed_count="${#CACHE_FAILED_TASKS[@]}"
  local running_count="${#GPU_BUSY[@]}"
  local pending_count=$(( CACHE_TOTAL_TASKS - completed_count - running_count ))
  local timestamp
  local message
  if (( pending_count < 0 )); then
    pending_count=0
  fi

  timestamp="$(date '+%Y-%m-%d %H:%M:%S')"
  message="[launcher] cache progress ${timestamp}: completed ${completed_count}/${CACHE_TOTAL_TASKS}, failed ${failed_count}, running ${running_count}, pending ${pending_count}"
  if [[ -n "${LAUNCH_PROGRESS_LOG:-}" ]]; then
    echo "$message" >>"$LAUNCH_PROGRESS_LOG"
  fi
}

maybe_print_cache_progress() {
  local force="${1:-0}"
  if (( PROGRESS_SECONDS <= 0 )); then
    return 0
  fi

  local now
  now="$(date +%s)"
  if (( force || LAST_PROGRESS_TS == 0 || now - LAST_PROGRESS_TS >= PROGRESS_SECONDS )); then
    print_cache_progress
    LAST_PROGRESS_TS="$now"
  fi
}

reap_one_finished_cache_task() {
  while true; do
    local gpu
    for gpu in "${GPU_IDS[@]}"; do
      local pid="${GPU_BUSY[$gpu]:-}"
      if [[ -z "$pid" ]]; then
        continue
      fi
      if ! kill -0 "$pid" 2>/dev/null; then
        local task_name="${PID_TO_JOB[$pid]}"
        local launch_log_path="${PID_TO_LOG[$pid]}"
        local exit_code
        set +e
        wait "$pid"
        exit_code=$?
        set -e

        unset "GPU_BUSY[$gpu]"
        unset "PID_TO_GPU[$pid]"
        unset "PID_TO_JOB[$pid]"
        unset "PID_TO_LOG[$pid]"

        CACHE_COMPLETED_TASKS+=("${task_name}:${exit_code}:cache")
        if (( exit_code != 0 )); then
          CACHE_FAILED_TASKS+=("${task_name}")
          echo "[launcher] Cache failed on GPU ${gpu}: ${task_name}, exit=${exit_code}, log=${launch_log_path}"
          if (( STOP_ON_FAILURE )); then
            CACHE_STOP_REQUESTED=1
          fi
        else
          echo "[launcher] Cache completed on GPU ${gpu}: ${task_name}"
        fi
        maybe_print_cache_progress 1
        CACHE_REAPED_GPU="$gpu"
        return 0
      fi
    done
    maybe_print_cache_progress 0
    sleep "$POLL_SECONDS"
  done
}

schedule_next_cache_task_on_gpu() {
  local gpu="$1"
  while (( CACHE_JOB_INDEX < CACHE_TOTAL_TASKS )); do
    if (( CACHE_STOP_REQUESTED )); then
      return 1
    fi
    if launch_cache_task "${CACHE_QUEUE[$CACHE_JOB_INDEX]}" "$gpu"; then
      ((CACHE_JOB_INDEX+=1))
      return 0
    fi
    ((CACHE_JOB_INDEX+=1))
  done
  return 1
}

run_cache_queue() {
  CACHE_JOB_INDEX=0
  CACHE_TOTAL_TASKS="${#CACHE_QUEUE[@]}"
  CACHE_STOP_REQUESTED=0
  local gpu

  for gpu in "${GPU_IDS[@]}"; do
    schedule_next_cache_task_on_gpu "$gpu" || true
    if (( CACHE_JOB_INDEX >= CACHE_TOTAL_TASKS || CACHE_STOP_REQUESTED )); then
      break
    fi
  done
  maybe_print_cache_progress 1

  while (( CACHE_JOB_INDEX < CACHE_TOTAL_TASKS && ! CACHE_STOP_REQUESTED )); do
    CACHE_REAPED_GPU=""
    reap_one_finished_cache_task
    schedule_next_cache_task_on_gpu "$CACHE_REAPED_GPU" || true
  done

  while (( ${#GPU_BUSY[@]} > 0 )); do
    CACHE_REAPED_GPU=""
    reap_one_finished_cache_task
  done
}

write_cache_completion_table() {
  local output_path="${LOG_ROOT}/launcher_cache_completed.tsv"
  local item
  local task_name
  local exit_code
  local mode
  {
    printf 'task_name\texit_code\tmode\n'
    for item in "${CACHE_COMPLETED_TASKS[@]}"; do
      IFS=':' read -r task_name exit_code mode <<<"$item"
      printf '%s\t%s\t%s\n' "$task_name" "$exit_code" "$mode"
    done
  } >"$output_path"
}

run_embedding_cache_precompute_stage() {
  build_cache_queue
  echo "[launcher] Precompute embedding caches before grid training: ${#CACHE_QUEUE[@]} tasks"
  LAST_PROGRESS_TS=0
  run_cache_queue
  write_cache_completion_table || true
  LAST_PROGRESS_TS=0

  if (( ${#CACHE_FAILED_TASKS[@]} > 0 )); then
    echo "[launcher] Failed cache tasks: ${CACHE_FAILED_TASKS[*]}" >&2
    return 1
  fi
  echo "[launcher] All embedding caches are ready."
}

print_launcher_progress() {
  if (( PROGRESS_SECONDS <= 0 )); then
    return 0
  fi

  local completed_count="${#COMPLETED_JOBS[@]}"
  local failed_count="${#FAILED_JOBS[@]}"
  local running_count="${#GPU_BUSY[@]}"
  local pending_count=$(( TOTAL_JOBS - completed_count - running_count ))
  local timestamp
  local message
  if (( pending_count < 0 )); then
    pending_count=0
  fi

  timestamp="$(date '+%Y-%m-%d %H:%M:%S')"
  message="[launcher] progress ${timestamp}: completed ${completed_count}/${TOTAL_JOBS}, failed ${failed_count}, running ${running_count}, pending ${pending_count}"
  if [[ -n "${LAUNCH_PROGRESS_LOG:-}" ]]; then
    echo "$message" >>"$LAUNCH_PROGRESS_LOG"
  fi
}

maybe_print_launcher_progress() {
  local force="${1:-0}"
  if (( PROGRESS_SECONDS <= 0 )); then
    return 0
  fi

  local now
  now="$(date +%s)"
  if (( force || LAST_PROGRESS_TS == 0 || now - LAST_PROGRESS_TS >= PROGRESS_SECONDS )); then
    print_launcher_progress
    LAST_PROGRESS_TS="$now"
  fi
}

cleanup_children() {
  local pid
  for pid in "${!PID_TO_GPU[@]}"; do
    kill "$pid" 2>/dev/null || true
  done
  wait || true
}

reap_one_finished_job() {
  while true; do
    local gpu
    for gpu in "${GPU_IDS[@]}"; do
      local pid="${GPU_BUSY[$gpu]:-}"
      if [[ -z "$pid" ]]; then
        continue
      fi
      if ! kill -0 "$pid" 2>/dev/null; then
        local job_key="${PID_TO_JOB[$pid]}"
        local launch_log_path="${PID_TO_LOG[$pid]}"
        local exit_code
        set +e
        wait "$pid"
        exit_code=$?
        set -e

        unset "GPU_BUSY[$gpu]"
        unset "PID_TO_GPU[$pid]"
        unset "PID_TO_JOB[$pid]"
        unset "PID_TO_LOG[$pid]"

        COMPLETED_JOBS+=("${job_key}:${exit_code}:run")
        if (( exit_code != 0 )); then
          FAILED_JOBS+=("${job_key}")
          echo "[launcher] Failed on GPU ${gpu}: ${job_key}, exit=${exit_code}, log=${launch_log_path}"
          if (( STOP_ON_FAILURE )); then
            STOP_REQUESTED=1
          fi
        else
          echo "[launcher] Completed on GPU ${gpu}: ${job_key}"
        fi
        maybe_print_launcher_progress 1
        REAPED_GPU="$gpu"
        return 0
      fi
    done
    maybe_print_launcher_progress 0
    sleep "$POLL_SECONDS"
  done
}

schedule_next_job_on_gpu() {
  local gpu="$1"
  while (( JOB_INDEX < TOTAL_JOBS )); do
    if (( STOP_REQUESTED )); then
      return 1
    fi
    if launch_job "${JOB_QUEUE[$JOB_INDEX]}" "$gpu"; then
      ((JOB_INDEX+=1))
      return 0
    fi
    ((JOB_INDEX+=1))
  done
  return 1
}

run_job_queue() {
  JOB_INDEX=0
  TOTAL_JOBS="${#JOB_QUEUE[@]}"
  STOP_REQUESTED=0
  local gpu

  for gpu in "${GPU_IDS[@]}"; do
    schedule_next_job_on_gpu "$gpu" || true
    if (( JOB_INDEX >= TOTAL_JOBS || STOP_REQUESTED )); then
      break
    fi
  done
  maybe_print_launcher_progress 1

  while (( JOB_INDEX < TOTAL_JOBS && ! STOP_REQUESTED )); do
    REAPED_GPU=""
    reap_one_finished_job
    schedule_next_job_on_gpu "$REAPED_GPU" || true
  done

  while (( ${#GPU_BUSY[@]} > 0 )); do
    REAPED_GPU=""
    reap_one_finished_job
  done
}

write_launcher_manifest() {
  LAUNCH_LOG_DIR="${LOG_ROOT}/launcher_logs"
  LAUNCH_PROGRESS_LOG="${LAUNCH_LOG_DIR}/launcher_progress.log"
  mkdir -p "$LOG_ROOT" "$LAUNCH_LOG_DIR"
  : >"$LAUNCH_PROGRESS_LOG"
  cat >"${LOG_ROOT}/launcher_config.txt" <<EOF
model_name=${LAUNCHER_MODEL_NAME}
conda_env=${CONDA_ENV}
finetune_method=${FINETUNE_METHOD}
task_name=${TASK_NAME}
resolved_tasks=${RESOLVED_TASKS[*]}
raw_download=${RAW_DOWNLOAD}
processed_download=${PROCESSED_DOWNLOAD}
log_root=${LOG_ROOT}
gpu_ids=${GPU_IDS[*]}
lr_values=${LR_VALUES[*]}
batch_size_values=${BATCH_SIZE_VALUES[*]}
weight_decay=${WEIGHT_DECAY}
optimizer_name=${OPTIMIZER_NAME}
epochs=${MAX_EPOCHS}
early_stopping_patience=${EARLY_STOPPING_PATIENCE}
warmup_ratio=${WARMUP_RATIO}
min_lr_ratio=${MIN_LR_RATIO}
seed=${RANDOM_SEED}
num_workers=${NUM_WORKERS}
chunk_forward_batch_size=${CHUNK_FORWARD_BATCH_SIZE}
use_embedding_cache=${USE_EMBEDDING_CACHE}
precompute_embedding_cache_first=${PRECOMPUTE_EMBEDDING_CACHE_FIRST}
embedding_cache_root=${EMBEDDING_CACHE_ROOT}
embedding_cache_dtype=${EMBEDDING_CACHE_DTYPE}
feature_train_dtype=${FEATURE_TRAIN_DTYPE}
extract_batch_size=${EXTRACT_BATCH_SIZE:-}
progress_seconds=${PROGRESS_SECONDS}
launcher_progress_log=${LAUNCH_PROGRESS_LOG}
extra_args=${EXTRA_ARGS[*]:-}
EOF
}

print_launcher_summary() {
  echo "[launcher] model_name: ${LAUNCHER_MODEL_NAME}"
  echo "[launcher] conda_env: ${CONDA_ENV}"
  echo "[launcher] finetune_method: ${FINETUNE_METHOD}"
  echo "[launcher] task_name: ${TASK_NAME}"
  echo "[launcher] resolved_tasks: ${RESOLVED_TASKS[*]}"
  echo "[launcher] processed_download: ${PROCESSED_DOWNLOAD}"
  echo "[launcher] log_root: ${LOG_ROOT}"
  echo "[launcher] gpu_ids: ${GPU_IDS[*]}"
  echo "[launcher] lr_values: ${LR_VALUES[*]}"
  echo "[launcher] batch_size_values: ${BATCH_SIZE_VALUES[*]}"
  echo "[launcher] total_jobs: ${#JOB_QUEUE[@]}"
  if (( USE_EMBEDDING_CACHE )) && [[ "$FINETUNE_METHOD" == "frozen_linear_probe" ]]; then
    echo "[launcher] embedding_cache: enabled"
    echo "[launcher] embedding_cache_root: ${EMBEDDING_CACHE_ROOT}"
    echo "[launcher] precompute_embedding_cache_first: ${PRECOMPUTE_EMBEDDING_CACHE_FIRST}"
    echo "[launcher] embedding_cache_dtype: ${EMBEDDING_CACHE_DTYPE}"
    echo "[launcher] feature_train_dtype: ${FEATURE_TRAIN_DTYPE}"
  else
    echo "[launcher] embedding_cache: disabled"
  fi
  if (( PROGRESS_SECONDS > 0 )); then
    echo "[launcher] progress_log: ${LAUNCH_PROGRESS_LOG}"
  fi
  if (( ${#EXTRA_ARGS[@]} > 0 )); then
    echo "[launcher] extra_train_args: ${EXTRA_ARGS[*]}"
  fi
}

collect_grid_results() {
  if (( ! COLLECT_RESULTS || DRY_RUN )); then
    return 0
  fi
  echo "[launcher] Collecting metrics/profile/runtime summaries."
  if ! conda run --no-capture-output -n "$COLLECT_CONDA_ENV" python -u "$COLLECT_SCRIPT" --log-root "$LOG_ROOT"; then
    echo "[launcher] Result collection failed; training outputs are still available under ${LOG_ROOT}." >&2
    return 1
  fi
}

write_completion_table() {
  local output_path="${LOG_ROOT}/launcher_completed.tsv"
  local item
  {
    printf 'job_key\texit_code\tmode\n'
    for item in "${COMPLETED_JOBS[@]}"; do
      printf '%s\n' "$item" | awk -F ':' '{print $1"::"$3"::"$5"\t"$6"\t"$7}'
    done
  } >"$output_path"
}

run_model_grid_search_launcher() {
  require_model_launcher_config || return 1
  parse_model_launcher_args "$@" || return 1
  validate_passthrough_args || return 1
  ensure_processed_download_ready || return 1
  resolve_requested_tasks || return 1
  build_job_queue
  write_launcher_manifest
  print_launcher_summary

  declare -gA GPU_BUSY=()
  declare -gA PID_TO_GPU=()
  declare -gA PID_TO_JOB=()
  declare -gA PID_TO_LOG=()
  declare -ga FAILED_JOBS=()
  declare -ga COMPLETED_JOBS=()
  declare -ga CACHE_FAILED_TASKS=()
  declare -ga CACHE_COMPLETED_TASKS=()
  declare -g LAST_PROGRESS_TS=0

  trap 'echo "[launcher] Interrupted; stopping running child jobs."; cleanup_children; exit 130' INT TERM

  if should_precompute_embedding_cache_first; then
    run_embedding_cache_precompute_stage || return 1
  fi

  run_job_queue
  write_completion_table || true
  COLLECT_FAILED=0
  collect_grid_results || COLLECT_FAILED=1

  echo "[launcher] All scheduled jobs finished."
  if (( ${#FAILED_JOBS[@]} > 0 )); then
    echo "[launcher] Failed jobs: ${FAILED_JOBS[*]}" >&2
    return 1
  fi
  if (( COLLECT_FAILED )); then
    echo "[launcher] Result collection failed." >&2
    return 1
  fi
  echo "[launcher] All jobs completed successfully."
}
