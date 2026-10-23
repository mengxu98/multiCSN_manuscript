# Sourced by the top-level entry points; this file does not schedule a pipeline.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LC_ALL=C
export PYTHONDONTWRITEBYTECODE=1

task_rscript() {
  local environment=R_4.5
  if [[ -f "$HOME/miniconda3/envs/$environment/bin/Rscript" ]]; then
    printf '%s\n' "$HOME/miniconda3/envs/$environment/bin/Rscript"
  else
    command -v Rscript
  fi
}

task_python() {
  if [[ -f "$HOME/.omicos/env/.venv/bin/python" ]]; then
    printf '%s\n' "$HOME/.omicos/env/.venv/bin/python"
  else
    command -v python3
  fi
}

r_string() {
  local value=$1
  value=${value//\\/\\\\}
  value=${value//\"/\\\"}
  value=${value//$'\n'/\\n}
  printf '"%s"' "$value"
}

run_task() (
  set -euo pipefail
  local script=$1 stage='' timed=false
  local command
  if [[ "$script" == *.R ]]; then
    command=("$(task_rscript "$script")" --vanilla "$script")
    if [[ $# == 2 ]]; then
      [[ "$2" =~ ^[A-Za-z][A-Za-z0-9_]*$ ]] || exit 1
      command=("$(task_rscript "$script")" --vanilla -e "source($(r_string "$script")); $2()")
    fi
  else
    command=("$(task_python)" "$script")
  fi
  case "$script" in
    code/gse274113/prepare_control_candidates.R) stage=gse274113_multicsn_candidate_construction ;;
    code/gse274113/fit_control_timepoints.R) stage=gse274113_multicsn_network_inference ;;
    code/gse274113/run_pando.R) stage=gse274113_pando_network_inference ;;
    code/gse274113/prepare_scenicplus_inputs.py) stage=gse274113_scenicplus_input_construction ;;
    code/gse274113/run_pycistopic.py) stage=gse274113_scenicplus_pycistopic ;;
    code/gse274113/build_control_object.R|code/gse274113/fit_replicate_stability.R|code/gse274113/build_heldout_truth.R|code/gse274113/evaluate_tf_gene_comparators.R) timed=true ;;
  esac
  if [[ -n "$stage" ]]; then
    [[ -d /proc && -x /usr/bin/time ]] || { echo 'Resource monitoring requires Linux /proc and GNU time.' >&2; exit 1; }
    local directory task_pid='' sampler_pid='' status=0 sampler_status=0
    directory=$(mktemp -d)
    cleanup() {
      [[ -z "$task_pid" ]] || kill -TERM -- "-$task_pid" 2>/dev/null || true
      [[ -z "$sampler_pid" ]] || kill -TERM "$sampler_pid" 2>/dev/null || true
      rm -rf -- "$directory"
    }
    trap cleanup EXIT
    trap 'exit 130' INT TERM
    setsid bash -c 'while [[ ! -f "$1/go" ]]; do sleep 0.02; done; shift; exec /usr/bin/time -v "$@"' bash "$directory" "${command[@]}" &
    task_pid=$!
    "$(task_rscript "$script")" --vanilla -e "source(\"code/functions/sample_process_tree_rss.R\"); sample_rss($task_pid, $(r_string "$directory"))" &
    sampler_pid=$!
    while [[ ! -f "$directory/ready" ]]; do
      kill -0 "$sampler_pid" 2>/dev/null || { wait "$sampler_pid" || true; exit 1; }
      sleep 0.02
    done
    touch "$directory/go"
    wait "$task_pid" || status=$?
    wait "$sampler_pid" || sampler_status=$?
    cat "$directory/summary.tsv"
    printf 'resource_command_exit_status\t%s\nresource_sampler_exit_status\t%s\n' "$status" "$sampler_status"
    task_pid=''; sampler_pid=''
    [[ "$status" == 0 && "$sampler_status" == 0 ]] || exit "$((status != 0 ? status : sampler_status))"
  elif "$timed"; then
    /usr/bin/time -v "${command[@]}" || exit "$?"
  else
    "${command[@]}" || exit "$?"
  fi
)
