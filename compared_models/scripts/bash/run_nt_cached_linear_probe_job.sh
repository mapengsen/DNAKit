#!/usr/bin/env bash

set -Eeuo pipefail

CONDA_ENV=""
EXTRACT_SCRIPT=""
TRAIN_SCRIPT=""
TASK_NAME=""
PROCESSED_DOWNLOAD=""
RAW_DOWNLOAD=""
LOG_ROOT=""
CACHE_PATH=""
MODEL_NAME=""
TOKEN_READOUT="auto"
DEVICE="cuda"
EPOCHS="50"
BATCH_SIZE="32"
EXTRACT_BATCH_SIZE=""
NUM_WORKERS="4"
LR="1e-3"
WEIGHT_DECAY="1e-4"
OPTIMIZER_NAME="adamw"
WARMUP_RATIO="0.1"
MIN_LR_RATIO="0.1"
EARLY_STOPPING_PATIENCE="15"
SEED="42"
MONITOR="accuracy"
CHUNK_FORWARD_BATCH_SIZE="8"
TENSORBOARD=1
ALLOW_REMOTE_FALLBACK=1
OVERWRITE_CACHE=0
SKIP_CACHE_EXTRACT=0
FEATURE_DTYPE="float32"
FEATURE_TRAIN_DTYPE="float32"

while (( $# > 0 )); do
  case "$1" in
    --conda-env) CONDA_ENV="$2"; shift 2 ;;
    --extract-script) EXTRACT_SCRIPT="$2"; shift 2 ;;
    --train-script) TRAIN_SCRIPT="$2"; shift 2 ;;
    --task-name) TASK_NAME="$2"; shift 2 ;;
    --processed-download) PROCESSED_DOWNLOAD="$2"; shift 2 ;;
    --raw-download) RAW_DOWNLOAD="$2"; shift 2 ;;
    --log-root) LOG_ROOT="$2"; shift 2 ;;
    --cache-path) CACHE_PATH="$2"; shift 2 ;;
    --model-name) MODEL_NAME="$2"; shift 2 ;;
    --token-readout) TOKEN_READOUT="$2"; shift 2 ;;
    --device) DEVICE="$2"; shift 2 ;;
    --epochs) EPOCHS="$2"; shift 2 ;;
    --batch-size) BATCH_SIZE="$2"; shift 2 ;;
    --extract-batch-size) EXTRACT_BATCH_SIZE="$2"; shift 2 ;;
    --num-workers) NUM_WORKERS="$2"; shift 2 ;;
    --lr) LR="$2"; shift 2 ;;
    --weight-decay) WEIGHT_DECAY="$2"; shift 2 ;;
    --optimizer-name) OPTIMIZER_NAME="$2"; shift 2 ;;
    --warmup-ratio) WARMUP_RATIO="$2"; shift 2 ;;
    --min-lr-ratio) MIN_LR_RATIO="$2"; shift 2 ;;
    --early-stopping-patience) EARLY_STOPPING_PATIENCE="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --monitor) MONITOR="$2"; shift 2 ;;
    --chunk-forward-batch-size) CHUNK_FORWARD_BATCH_SIZE="$2"; shift 2 ;;
    --feature-dtype) FEATURE_DTYPE="$2"; shift 2 ;;
    --feature-train-dtype) FEATURE_TRAIN_DTYPE="$2"; shift 2 ;;
    --tensorboard) TENSORBOARD=1; shift ;;
    --no-tensorboard) TENSORBOARD=0; shift ;;
    --allow-remote-fallback) ALLOW_REMOTE_FALLBACK=1; shift ;;
    --no-allow-remote-fallback) ALLOW_REMOTE_FALLBACK=0; shift ;;
    --overwrite-cache) OVERWRITE_CACHE=1; shift ;;
    --skip-cache-extract) SKIP_CACHE_EXTRACT=1; shift ;;
    *)
      echo "[cached-job] Unknown arg: $1" >&2
      exit 2
      ;;
  esac
done

if [[ -z "$CONDA_ENV" || -z "$EXTRACT_SCRIPT" || -z "$TRAIN_SCRIPT" || -z "$TASK_NAME" || -z "$PROCESSED_DOWNLOAD" || -z "$LOG_ROOT" || -z "$CACHE_PATH" || -z "$MODEL_NAME" ]]; then
  echo "[cached-job] Missing required args." >&2
  exit 2
fi

if [[ -z "$EXTRACT_BATCH_SIZE" ]]; then
  EXTRACT_BATCH_SIZE="$BATCH_SIZE"
fi

EXTRACT_CMD=(
  conda run --no-capture-output -n "$CONDA_ENV" python -u "$EXTRACT_SCRIPT"
  --task-name "$TASK_NAME"
  --processed-download "$PROCESSED_DOWNLOAD"
  --cache-path "$CACHE_PATH"
  --model-name "$MODEL_NAME"
  --raw-download "$RAW_DOWNLOAD"
  --token-readout "$TOKEN_READOUT"
  --finetune-method frozen_linear_probe
  --device "$DEVICE"
  --batch-size "$EXTRACT_BATCH_SIZE"
  --num-workers "$NUM_WORKERS"
  --seed "$SEED"
  --chunk-forward-batch-size "$CHUNK_FORWARD_BATCH_SIZE"
  --feature-dtype "$FEATURE_DTYPE"
  --skip-existing
)
if (( ALLOW_REMOTE_FALLBACK )); then
  EXTRACT_CMD+=(--allow-remote-fallback)
else
  EXTRACT_CMD+=(--no-allow-remote-fallback)
fi
if (( OVERWRITE_CACHE )); then
  EXTRACT_CMD+=(--overwrite)
fi

TRAIN_CMD=(
  conda run --no-capture-output -n "$CONDA_ENV" python -u "$TRAIN_SCRIPT"
  --task-name "$TASK_NAME"
  --processed-download "$PROCESSED_DOWNLOAD"
  --raw-download "$RAW_DOWNLOAD"
  --log-root "$LOG_ROOT"
  --cache-path "$CACHE_PATH"
  --model-name "$MODEL_NAME"
  --device "$DEVICE"
  --epochs "$EPOCHS"
  --batch-size "$BATCH_SIZE"
  --num-workers "$NUM_WORKERS"
  --lr "$LR"
  --weight-decay "$WEIGHT_DECAY"
  --optimizer-name "$OPTIMIZER_NAME"
  --warmup-ratio "$WARMUP_RATIO"
  --min-lr-ratio "$MIN_LR_RATIO"
  --early-stopping-patience "$EARLY_STOPPING_PATIENCE"
  --seed "$SEED"
  --monitor "$MONITOR"
  --feature-train-dtype "$FEATURE_TRAIN_DTYPE"
)
if (( TENSORBOARD )); then
  TRAIN_CMD+=(--tensorboard)
else
  TRAIN_CMD+=(--no-tensorboard)
fi

if (( SKIP_CACHE_EXTRACT )); then
  if [[ ! -f "$CACHE_PATH" ]]; then
    echo "[cached-job] Missing precomputed cache: ${CACHE_PATH}" >&2
    exit 1
  fi
  echo "[cached-job] Reuse precomputed cache: ${CACHE_PATH}"
else
  echo "[cached-job] Ensure cache: ${CACHE_PATH}"
  "${EXTRACT_CMD[@]}"
fi
echo "[cached-job] Train cached linear probe: ${TASK_NAME}, lr=${LR}, batch_size=${BATCH_SIZE}"
"${TRAIN_CMD[@]}"
