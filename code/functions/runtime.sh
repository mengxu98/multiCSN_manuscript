# Sourced by the top-level entry points; this file does not schedule a pipeline.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1
export LC_ALL=C

task_rscript() {
  local environment=R_4.4
  [[ "$1" == code/beeline/* ]] && environment=beeline_r453
  if [[ -f "$HOME/miniconda3/envs/$environment/bin/Rscript" ]]; then
    printf '%s\n' "$HOME/miniconda3/envs/$environment/bin/Rscript"
  else
    command -v Rscript
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
  local script=$1 stamp log='' stage='' timed=false
  local command
  if [[ "$script" == *.R ]]; then
    command=("$(task_rscript "$script")" --vanilla "$script")
    if [[ $# == 2 ]]; then
      [[ "$2" =~ ^[A-Za-z][A-Za-z0-9_]*$ ]] || exit 1
      command=("$(task_rscript "$script")" --vanilla -e "source($(r_string "$script")); $2()")
    fi
  else
    command=(python3 "$script")
  fi
  stamp="$(date -u +%Y%m%dT%H%M%SZ)_${BASHPID:-$$}"
  if [[ "$script" == code/scmultisim/* ]]; then
    log="logs/scmultisim/$(basename "$script" .R)_$stamp.log"
  elif [[ "$script" == code/gse274113/* ]]; then
    local folder='' name=''
    case "$(basename "$script")" in
      build_control_object.R) folder=preparation; name=build_gse274113_control_object ;;
      prepare_control_candidates.R) folder=candidate_preparation; name=prepare_gse274113_control_candidates; stage=gse274113_multicsn_candidate_construction ;;
      fit_control_timepoints.R) folder=inference; name=run_gse274113_control_timepoint_multicsn; stage=gse274113_multicsn_network_inference ;;
      fit_replicate_stability.R) folder=replicate_stability; name=fit_gse274113_replicate_stability ;;
      build_heldout_truth.R) folder=causal_truth; name=build_gse274113_heldout_pseudobulk_truth ;;
      run_pando.R) folder=comparators; name=run_pando_gse274113_control_full; stage=gse274113_pando_network_inference ;;
      prepare_scenicplus_inputs.py) folder=comparators; name=prepare_scenicplus_gse274113_inputs; stage=gse274113_scenicplus_input_construction ;;
      run_pycistopic.py) folder=comparators; name=run_pycistopic_gse274113_full; stage=gse274113_scenicplus_pycistopic ;;
      run_scenicplus.py) folder=comparators; name=run_scenicplus_gse274113 ;;
      evaluate_tf_gene_comparators.R) folder=causal_evaluation; name=evaluate_gse274113_causal_comparators_full ;;
    esac
    if [[ -n "$folder" ]]; then
      log="logs/gse274113/$folder/${name}_$stamp.log"
      [[ "$script" == *run_scenicplus.py ]] || timed=true
    fi
  fi
  if [[ -n "$log" ]]; then
    mkdir -p "$(dirname "$log")"
    # All following output, including GNU time and RSS provenance, shares the log.
    exec > >(tee "$log") 2>&1
  fi
  printf 'stage_started\t%s\t%s\n' "$(basename "$script" .R)" "$(date -u +%FT%TZ)"
  if [[ -n "$stage" ]]; then
    [[ -d /proc && -x /usr/bin/time ]] || { echo 'Resource monitoring requires Linux /proc and GNU time.' >&2; exit 1; }
    local directory task_pid='' sampler_pid='' status=0 sampler_status=0 command_text
    directory=$(mktemp -d)
    cleanup() {
      [[ -z "$task_pid" ]] || kill -TERM -- "-$task_pid" 2>/dev/null || true
      [[ -z "$sampler_pid" ]] || kill -TERM "$sampler_pid" 2>/dev/null || true
      rm -rf -- "$directory"
    }
    trap cleanup EXIT
    trap 'exit 130' INT TERM
    printf -v command_text '%q ' "${command[@]}"
    # The gate lets the sampler record provenance before the timed task begins.
    setsid bash -c 'while [[ ! -f "$1/go" ]]; do sleep 0.02; done; shift; exec /usr/bin/time -v "$@"' bash "$directory" "${command[@]}" &
    task_pid=$!
    "$(task_rscript "$script")" --vanilla -e "source(\"code/functions/sample_process_tree_rss.R\"); sample_rss($task_pid, $(r_string "$directory"), $(r_string "$script"), $(r_string "$stage"), $(r_string "$command_text"))" &
    sampler_pid=$!
    while [[ ! -f "$directory/ready" ]]; do
      kill -0 "$sampler_pid" 2>/dev/null || { wait "$sampler_pid" || true; exit 1; }
      sleep 0.02
    done
    cat "$directory/header.tsv"
    touch "$directory/go"
    wait "$task_pid" || status=$?
    wait "$sampler_pid" || sampler_status=$?
    cat "$directory/summary.tsv"
    printf 'resource_command_exit_status\t%s\nresource_sampler_exit_status\t%s\n' "$status" "$sampler_status"
    printf 'resource_block_completed_utc\t%s\nresource_block_completed\t%s\n' "$(date -u +%FT%TZ)" "$(cat "$directory/run_id")"
    task_pid=''; sampler_pid=''
    [[ "$status" == 0 && "$sampler_status" == 0 ]] || exit "$((status != 0 ? status : sampler_status))"
  elif "$timed"; then
    /usr/bin/time -v "${command[@]}" || exit "$?"
  else
    "${command[@]}" || exit "$?"
  fi
  printf 'stage_completed\t%s\t%s\n' "$(basename "$script" .R)" "$(date -u +%FT%TZ)"
)
