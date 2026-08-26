#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
EXTRACT_SCRIPT="${REPO_ROOT}/compared_models/scripts/extract_nt_embedding_cache.py"

TASK_NAME="${TASK_NAME:-promoter_tata}"
MASK_RATIOS_RAW="${MASK_RATIOS:-2 5 10 20 40}"
MODEL_LIST_RAW="${MODEL_LIST:-alphagenome dnabert2 enformer grover hyenadna ntv2 janusdna caduceus lucaone}"
MASK_ROOT="${MASK_ROOT:-data/random_n_mask/nt/seed42}"
RAW_DOWNLOAD="${RAW_DOWNLOAD:-data/raw_download/nucleotide_transformer_downstream_tasks_revised}"
LOG_BASE="${LOG_BASE:-outputs/visualdna_benchmark/random_n_mask/seed42}"
EVAL_LOG_BASE="${EVAL_LOG_BASE:-outputs/visualdna_benchmark/random_n_mask_eval/seed42_from_clean}"
CLEAN_LOG_PREFIX="${CLEAN_LOG_PREFIX:-outputs/visualdna_benchmark/ood_nt_k6}"
CLEAN_RUN_NAME="${CLEAN_RUN_NAME:-lr_1e_4_bs_32}"
CHECKPOINT_MANIFEST="${CHECKPOINT_MANIFEST:-}"
GPU_ID_LIST_GLOBAL="${GPU_ID_LIST:-0 1 2 3 4}"
LR_VALUES="${LR_VALUES:-1e-3 1e-4 1e-5}"
BATCH_SIZE_VALUES="${BATCH_SIZE_VALUES:-16 32 64}"
EVAL_BATCH_SIZE="${EVAL_BATCH_SIZE:-64}"
EPOCHS="${EPOCHS:-50}"
EARLY_STOPPING_PATIENCE="${EARLY_STOPPING_PATIENCE:-15}"
WEIGHT_DECAY="${WEIGHT_DECAY:-1e-4}"
OPTIMIZER_NAME="${OPTIMIZER_NAME:-adamw}"
WARMUP_RATIO="${WARMUP_RATIO:-0.1}"
MIN_LR_RATIO="${MIN_LR_RATIO:-0.1}"
SEED="${SEED:-42}"
NUM_WORKERS="${NUM_WORKERS:-4}"
FEATURE_DTYPE="${FEATURE_DTYPE:-float32}"
FEATURE_TRAIN_DTYPE="${FEATURE_TRAIN_DTYPE:-float32}"
RUN_CACHE="${RUN_CACHE:-1}"
RUN_TRAIN="${RUN_TRAIN:-1}"
RUN_EVAL="${RUN_EVAL:-0}"
DRY_RUN=0
SKIP_EXISTING=1

export XLA_PYTHON_CLIENT_PREALLOCATE="${XLA_PYTHON_CLIENT_PREALLOCATE:-false}"

print_usage() {
  cat <<EOF
Usage:
  bash compared_models/scripts/bash/run_random_n_mask_compared_models_parallel_cache.sh [args]

Args:
  --task-name VALUE             default: ${TASK_NAME}
  --mask-ratios "2 5 10"        default: ${MASK_RATIOS_RAW}
  --models "dnabert2 ntv2"      default: ${MODEL_LIST_RAW}
  --mask-root PATH              default: ${MASK_ROOT}
  --raw-download PATH           default: ${RAW_DOWNLOAD}
  --log-base PATH               default: ${LOG_BASE}
  --eval-log-base PATH          default: ${EVAL_LOG_BASE}
  --clean-log-prefix PATH       default: ${CLEAN_LOG_PREFIX}
  --clean-run-name VALUE        default: ${CLEAN_RUN_NAME}
  --checkpoint-manifest PATH    optional TSV with compared/model/task best checkpoints
  --gpu-id-list "0 1 2 3 4"     default: ${GPU_ID_LIST_GLOBAL}
  --cache-only                  only precompute embedding caches
  --train-only                  only run cached linear-probe grids
  --eval-only                   precompute caches and evaluate clean linear probes on masked test
  --no-skip-existing            do not pass --skip-existing to grid training
  --dry-run                     print commands without running
  --help                        show this help

This script supports one NT task at a time. Use one mask ratio directory per
processed root, for example:
  data/random_n_mask/nt/seed42/mask_10/promoter_tata/raw/promoter_tata.parquet
EOF
}

while (( $# > 0 )); do
  case "$1" in
    --task-name)
      TASK_NAME="$2"
      shift 2
      ;;
    --task-name=*)
      TASK_NAME="${1#*=}"
      shift
      ;;
    --mask-ratios)
      MASK_RATIOS_RAW="$2"
      shift 2
      ;;
    --mask-ratios=*)
      MASK_RATIOS_RAW="${1#*=}"
      shift
      ;;
    --models)
      MODEL_LIST_RAW="$2"
      shift 2
      ;;
    --models=*)
      MODEL_LIST_RAW="${1#*=}"
      shift
      ;;
    --mask-root)
      MASK_ROOT="$2"
      shift 2
      ;;
    --mask-root=*)
      MASK_ROOT="${1#*=}"
      shift
      ;;
    --raw-download)
      RAW_DOWNLOAD="$2"
      shift 2
      ;;
    --raw-download=*)
      RAW_DOWNLOAD="${1#*=}"
      shift
      ;;
    --log-base)
      LOG_BASE="$2"
      shift 2
      ;;
    --log-base=*)
      LOG_BASE="${1#*=}"
      shift
      ;;
    --eval-log-base)
      EVAL_LOG_BASE="$2"
      shift 2
      ;;
    --eval-log-base=*)
      EVAL_LOG_BASE="${1#*=}"
      shift
      ;;
    --clean-log-prefix)
      CLEAN_LOG_PREFIX="$2"
      shift 2
      ;;
    --clean-log-prefix=*)
      CLEAN_LOG_PREFIX="${1#*=}"
      shift
      ;;
    --clean-run-name)
      CLEAN_RUN_NAME="$2"
      shift 2
      ;;
    --clean-run-name=*)
      CLEAN_RUN_NAME="${1#*=}"
      shift
      ;;
    --checkpoint-manifest)
      CHECKPOINT_MANIFEST="$2"
      shift 2
      ;;
    --checkpoint-manifest=*)
      CHECKPOINT_MANIFEST="${1#*=}"
      shift
      ;;
    --gpu-id-list)
      GPU_ID_LIST_GLOBAL="$2"
      shift 2
      ;;
    --gpu-id-list=*)
      GPU_ID_LIST_GLOBAL="${1#*=}"
      shift
      ;;
    --cache-only)
      RUN_CACHE=1
      RUN_TRAIN=0
      RUN_EVAL=0
      shift
      ;;
    --train-only)
      RUN_CACHE=0
      RUN_TRAIN=1
      RUN_EVAL=0
      shift
      ;;
    --eval-only)
      RUN_CACHE=1
      RUN_TRAIN=0
      RUN_EVAL=1
      shift
      ;;
    --no-skip-existing)
      SKIP_EXISTING=0
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --help|-h)
      print_usage
      exit 0
      ;;
    *)
      echo "[random-mask] Unknown arg: $1" >&2
      print_usage >&2
      exit 2
      ;;
  esac
done

MASK_RATIOS_RAW="${MASK_RATIOS_RAW//,/ }"
MODEL_LIST_RAW="${MODEL_LIST_RAW//,/ }"
GPU_ID_LIST_GLOBAL="${GPU_ID_LIST_GLOBAL//,/ }"
read -r -a MASK_RATIOS <<<"$MASK_RATIOS_RAW"
read -r -a MODEL_LIST <<<"$MODEL_LIST_RAW"

if [[ "$TASK_NAME" == "all" || "$TASK_NAME" == *","* ]]; then
  echo "[random-mask] This parallel-cache script currently supports one task at a time, got: ${TASK_NAME}" >&2
  exit 2
fi
if (( ${#MASK_RATIOS[@]} == 0 )); then
  echo "[random-mask] No mask ratios were provided." >&2
  exit 2
fi
if (( ${#MODEL_LIST[@]} == 0 )); then
  echo "[random-mask] No models were provided." >&2
  exit 2
fi

model_config() {
  local model="$1"
  case "$model" in
    alphagenome)
      SCRIPT="${SCRIPT_DIR}/run_grid_alphagenome.sh"
      CONDA_ENV="alphagenome"
      GPU_IDS_RAW="0 1 2 3"
      EXTRACT_BS="32"
      CHUNK_BS="1"
      LOG_NAME="grid_alphagenome_frozen_linear_probe"
      ;;
    dnabert2)
      SCRIPT="${SCRIPT_DIR}/run_grid_dnabert2.sh"
      CONDA_ENV="visualdna"
      GPU_IDS_RAW="0 1 2 3"
      EXTRACT_BS="32"
      CHUNK_BS="8"
      LOG_NAME="grid_dnabert2_frozen_linear_probe"
      ;;
    enformer)
      SCRIPT="${SCRIPT_DIR}/run_grid_enformer.sh"
      CONDA_ENV="visualdna"
      GPU_IDS_RAW="0 1 2 3"
      EXTRACT_BS="4"
      CHUNK_BS="8"
      LOG_NAME="grid_enformer_frozen_linear_probe"
      ;;
    grover)
      SCRIPT="${SCRIPT_DIR}/run_grid_grover.sh"
      CONDA_ENV="visualdna"
      GPU_IDS_RAW="0 1 2 3"
      EXTRACT_BS="32"
      CHUNK_BS="8"
      LOG_NAME="grid_grover_frozen_linear_probe"
      ;;
    hyenadna)
      SCRIPT="${SCRIPT_DIR}/run_grid_hyenadna.sh"
      CONDA_ENV="visualdna"
      GPU_IDS_RAW="0 1 2 3"
      EXTRACT_BS="32"
      CHUNK_BS="8"
      LOG_NAME="grid_hyenadna_frozen_linear_probe"
      ;;
    ntv2)
      SCRIPT="${SCRIPT_DIR}/run_grid_ntv2.sh"
      CONDA_ENV="visualdna"
      GPU_IDS_RAW="0 1"
      EXTRACT_BS="32"
      CHUNK_BS="8"
      LOG_NAME="grid_ntv2_frozen_linear_probe"
      ;;
    janusdna)
      SCRIPT="${SCRIPT_DIR}/run_grid_janusdna.sh"
      CONDA_ENV="janusdna"
      GPU_IDS_RAW="0 1 2 3"
      EXTRACT_BS="128"
      CHUNK_BS="128"
      LOG_NAME="grid_janusdna_batched_frozen_linear_probe"
      ;;
    caduceus)
      SCRIPT="${SCRIPT_DIR}/run_grid_caduceus.sh"
      CONDA_ENV="visualdna"
      GPU_IDS_RAW="0 1 2 3"
      EXTRACT_BS="256"
      CHUNK_BS="256"
      LOG_NAME="grid_caduceus_frozen_linear_probe"
      ;;
    lucaone)
      SCRIPT="${SCRIPT_DIR}/run_grid_lucaone.sh"
      CONDA_ENV="visualdna"
      GPU_IDS_RAW="0 1 3 4 5 6 7"
      EXTRACT_BS="32"
      CHUNK_BS="8"
      LOG_NAME="grid_lucaone_frozen_linear_probe"
      ;;
    *)
      echo "[random-mask] Unknown model: ${model}" >&2
      return 1
      ;;
  esac
  GPU_IDS_RAW="$GPU_ID_LIST_GLOBAL"
}

log_root_for() {
  local ratio="$1"
  printf '%s/mask_%s/%s\n' "$LOG_BASE" "$ratio" "$LOG_NAME"
}

processed_root_for() {
  local ratio="$1"
  printf '%s/mask_%s\n' "$MASK_ROOT" "$ratio"
}

cache_path_for() {
  local ratio="$1"
  printf '%s/embedding_cache/%s/%s/features.pt\n' "$(log_root_for "$ratio")" "$MODEL_NAME" "$TASK_NAME"
}

eval_log_root_for() {
  local ratio="$1"
  printf '%s/mask_%s/%s/%s/from_clean_%s\n' "$EVAL_LOG_BASE" "$ratio" "$LOG_NAME" "$TASK_NAME" "$CLEAN_RUN_NAME"
}

clean_checkpoint_path_for() {
  if [[ -n "$CHECKPOINT_MANIFEST" ]]; then
    awk -F '\t' \
      -v kind="compared" \
      -v model="$MODEL_NAME" \
      -v task="$TASK_NAME" \
      'NR > 1 && $1 == kind && $2 == model && $3 == task { print $5; found = 1; exit }
       END { if (!found) exit 1 }' \
      "$CHECKPOINT_MANIFEST"
    return
  fi
  printf '%s_%s/%s/%s/checkpoints/best.pt\n' "$CLEAN_LOG_PREFIX" "$LOG_NAME" "$TASK_NAME" "$CLEAN_RUN_NAME"
}

run_or_print() {
  if (( DRY_RUN )); then
    printf '%q ' "$@"
    printf '\n'
    return 0
  fi
  "$@"
}

launch_cache_one() {
  local ratio="$1"
  local gpu="$2"
  local processed_root
  local log_root
  local cache_path
  local launch_log

  processed_root="$(processed_root_for "$ratio")"
  log_root="$(log_root_for "$ratio")"
  cache_path="$(cache_path_for "$ratio")"
  launch_log="${log_root}/launcher_logs/cache_${TASK_NAME}_gpu_${gpu}.log"

  if [[ ! -d "$processed_root" ]]; then
    echo "[random-mask] Missing processed root: ${processed_root}" >&2
    return 1
  fi

  mkdir -p "$(dirname "$cache_path")" "$(dirname "$launch_log")"

  echo "[random-mask] Cache ${MODEL_NAME} mask_${ratio}/${TASK_NAME} -> GPU ${gpu}"
  if (( DRY_RUN )); then
    printf 'CUDA_VISIBLE_DEVICES=%q ' "$gpu"
    run_or_print \
      conda run --no-capture-output -n "$CONDA_ENV" python -u "$EXTRACT_SCRIPT" \
      --task-name "$TASK_NAME" \
      --processed-download "$processed_root" \
      --cache-path "$cache_path" \
      --model-name "$MODEL_NAME" \
      --raw-download "$RAW_DOWNLOAD" \
      --token-readout auto \
      --finetune-method frozen_linear_probe \
      --device cuda \
      --batch-size "$EXTRACT_BS" \
      --num-workers "$NUM_WORKERS" \
      --seed "$SEED" \
      --chunk-forward-batch-size "$CHUNK_BS" \
      --feature-dtype "$FEATURE_DTYPE" \
      --allow-remote-fallback \
      --skip-existing
    return 0
  fi

  (
    cd "$REPO_ROOT"
    CUDA_VISIBLE_DEVICES="$gpu" conda run --no-capture-output -n "$CONDA_ENV" python -u "$EXTRACT_SCRIPT" \
      --task-name "$TASK_NAME" \
      --processed-download "$processed_root" \
      --cache-path "$cache_path" \
      --model-name "$MODEL_NAME" \
      --raw-download "$RAW_DOWNLOAD" \
      --token-readout auto \
      --finetune-method frozen_linear_probe \
      --device cuda \
      --batch-size "$EXTRACT_BS" \
      --num-workers "$NUM_WORKERS" \
      --seed "$SEED" \
      --chunk-forward-batch-size "$CHUNK_BS" \
      --feature-dtype "$FEATURE_DTYPE" \
      --allow-remote-fallback \
      --skip-existing
  ) >"$launch_log" 2>&1
}

wait_cache_wave() {
  local failed=0
  local item
  local pid
  local ratio
  local gpu

  for item in "$@"; do
    IFS=':' read -r pid ratio gpu <<<"$item"
    if wait "$pid"; then
      echo "[random-mask] Cache done ${MODEL_NAME} mask_${ratio} on GPU ${gpu}"
    else
      echo "[random-mask] Cache failed ${MODEL_NAME} mask_${ratio} on GPU ${gpu}" >&2
      failed=1
    fi
  done
  return "$failed"
}

precompute_model_caches() {
  local gpu_ids=()
  local running=()
  local ratio
  local gpu
  local pid
  local index=0

  read -r -a gpu_ids <<<"${GPU_IDS_RAW//,/ }"
  if (( ${#gpu_ids[@]} == 0 )); then
    echo "[random-mask] ${MODEL_NAME} has no GPUs configured." >&2
    return 1
  fi

  echo "[random-mask] Precompute ${MODEL_NAME} caches with GPUs: ${gpu_ids[*]}"
  for ratio in "${MASK_RATIOS[@]}"; do
    gpu="${gpu_ids[$(( index % ${#gpu_ids[@]} ))]}"
    if (( DRY_RUN )); then
      launch_cache_one "$ratio" "$gpu"
    else
      launch_cache_one "$ratio" "$gpu" &
      pid="$!"
      running+=("${pid}:${ratio}:${gpu}")
    fi
    index=$(( index + 1 ))

    if (( ! DRY_RUN && ${#running[@]} >= ${#gpu_ids[@]} )); then
      wait_cache_wave "${running[@]}"
      running=()
    fi
  done

  if (( ! DRY_RUN && ${#running[@]} > 0 )); then
    wait_cache_wave "${running[@]}"
  fi
}

run_model_train_grid_for_ratio() {
  local ratio="$1"
  local processed_root
  local log_root
  local train_cmd=()

  processed_root="$(processed_root_for "$ratio")"
  log_root="$(log_root_for "$ratio")"

  train_cmd=(
    bash "$SCRIPT"
    --task-name "$TASK_NAME"
    --raw-download "$RAW_DOWNLOAD"
    --processed-download "$processed_root"
    --log-root "$log_root"
    --embedding-cache-root "${log_root}/embedding_cache"
    --gpu-id-list "$GPU_IDS_RAW"
    --lr-values "$LR_VALUES"
    --batch-size-values "$BATCH_SIZE_VALUES"
    --extract-batch-size "$EXTRACT_BS"
    --finetune-method frozen_linear_probe
    --chunk-forward-batch-size "$CHUNK_BS"
    --epochs "$EPOCHS"
    --early-stopping-patience "$EARLY_STOPPING_PATIENCE"
    --weight-decay "$WEIGHT_DECAY"
    --optimizer-name "$OPTIMIZER_NAME"
    --warmup-ratio "$WARMUP_RATIO"
    --min-lr-ratio "$MIN_LR_RATIO"
    --seed "$SEED"
    --feature-train-dtype "$FEATURE_TRAIN_DTYPE"
    --precompute-embedding-cache-first
  )
  if (( SKIP_EXISTING )); then
    train_cmd+=(--skip-existing)
  fi
  if (( DRY_RUN )); then
    train_cmd+=(--dry-run)
  fi

  echo "[random-mask] Train grid ${MODEL_NAME} mask_${ratio}"
  run_or_print "${train_cmd[@]}"
}

run_model_eval_for_ratio() {
  local ratio="$1"
  local gpu="$2"
  local processed_root
  local cache_path
  local eval_log_root
  local checkpoint_path

  processed_root="$(processed_root_for "$ratio")"
  cache_path="$(cache_path_for "$ratio")"
  eval_log_root="$(eval_log_root_for "$ratio")"
  checkpoint_path="$(clean_checkpoint_path_for)"

  if (( ! DRY_RUN )) && [[ ! -d "$processed_root" ]]; then
    echo "[random-mask] Missing processed root: ${processed_root}" >&2
    return 1
  fi
  if (( ! DRY_RUN )) && [[ ! -f "$cache_path" ]]; then
    echo "[random-mask] Missing embedding cache: ${cache_path}" >&2
    return 1
  fi
  if (( ! DRY_RUN )) && [[ ! -f "$checkpoint_path" ]]; then
    echo "[random-mask] Missing clean checkpoint: ${checkpoint_path}" >&2
    return 1
  fi

  mkdir -p "$eval_log_root"
  echo "[random-mask] Eval clean ${MODEL_NAME} checkpoint on mask_${ratio} -> GPU ${gpu}"
  if (( DRY_RUN )); then
    printf 'CUDA_VISIBLE_DEVICES=%q ' "$gpu"
    run_or_print \
      conda run --no-capture-output -n "$CONDA_ENV" python -u "${REPO_ROOT}/compared_models/scripts/train_nt_cached_linear_probe.py" \
      --task-name "$TASK_NAME" \
      --processed-download "$processed_root" \
      --raw-download "$RAW_DOWNLOAD" \
      --log-root "$eval_log_root" \
      --cache-path "$cache_path" \
      --model-name "$MODEL_NAME" \
      --device cuda \
      --batch-size "$EVAL_BATCH_SIZE" \
      --num-workers "$NUM_WORKERS" \
      --seed "$SEED" \
      --feature-train-dtype "$FEATURE_TRAIN_DTYPE" \
      --monitor mcc \
      --eval-only \
      --checkpoint-path "$checkpoint_path"
    return 0
  fi

  (
    cd "$REPO_ROOT"
    CUDA_VISIBLE_DEVICES="$gpu" conda run --no-capture-output -n "$CONDA_ENV" python -u "${REPO_ROOT}/compared_models/scripts/train_nt_cached_linear_probe.py" \
      --task-name "$TASK_NAME" \
      --processed-download "$processed_root" \
      --raw-download "$RAW_DOWNLOAD" \
      --log-root "$eval_log_root" \
      --cache-path "$cache_path" \
      --model-name "$MODEL_NAME" \
      --device cuda \
      --batch-size "$EVAL_BATCH_SIZE" \
      --num-workers "$NUM_WORKERS" \
      --seed "$SEED" \
      --feature-train-dtype "$FEATURE_TRAIN_DTYPE" \
      --monitor mcc \
      --eval-only \
      --checkpoint-path "$checkpoint_path"
  )
}

run_model_eval_all_ratios() {
  local gpu_ids=()
  local ratio
  local gpu
  local index=0

  read -r -a gpu_ids <<<"${GPU_IDS_RAW//,/ }"
  if (( ${#gpu_ids[@]} == 0 )); then
    echo "[random-mask] ${MODEL_NAME} has no GPUs configured." >&2
    return 1
  fi

  for ratio in "${MASK_RATIOS[@]}"; do
    gpu="${gpu_ids[$(( index % ${#gpu_ids[@]} ))]}"
    run_model_eval_for_ratio "$ratio" "$gpu"
    index=$(( index + 1 ))
  done
}

main() {
  local model
  local ratio

  cd "$REPO_ROOT"
  if [[ -n "$CHECKPOINT_MANIFEST" && ! -f "$CHECKPOINT_MANIFEST" ]]; then
    echo "[random-mask] Missing checkpoint manifest: ${CHECKPOINT_MANIFEST}" >&2
    exit 1
  fi
  for model in "${MODEL_LIST[@]}"; do
    MODEL_NAME="$model"
    model_config "$MODEL_NAME"
    echo "[random-mask] ===== ${MODEL_NAME} ====="
    if (( RUN_CACHE )); then
      precompute_model_caches
    fi
    if (( RUN_TRAIN )); then
      for ratio in "${MASK_RATIOS[@]}"; do
        run_model_train_grid_for_ratio "$ratio"
      done
    fi
    if (( RUN_EVAL )); then
      run_model_eval_all_ratios
    fi
  done
}

main "$@"
