#!/usr/bin/env bash
set -euo pipefail
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
test_root=$(mktemp -d)
trap 'rm -rf -- "$test_root"' EXIT
mkdir -p "$test_root/code/functions" "$test_root/code/scmultisim" "$test_root/code/gse274113"
cp "$repo_root/code/functions/runtime.sh" "$repo_root/code/functions/sample_process_tree_rss.R" "$test_root/code/functions/"
cd "$test_root"
source code/functions/runtime.sh
cat > code/scmultisim/fixture.R <<'RS'
stopifnot(Sys.getenv("OMP_NUM_THREADS")=="1", Sys.getenv("LC_ALL")=="C")
cat("stdout fixture\n"); message("stderr fixture")
RS
run_task code/scmultisim/fixture.R > output.log 2>&1
log=$(find logs/scmultisim -name '*.log')
for text in stage_started 'stdout fixture' 'stderr fixture' stage_completed; do
  grep -q "$text" "$log"
done
cat > code/scmultisim/fail.R <<'RS'
quit(status=7L)
RS
if run_task code/scmultisim/fail.R > failure.log 2>&1; then
  echo 'Expected task failure' >&2; exit 1
else
  [[ "$?" == 7 ]]
fi
! grep -q stage_completed failure.log
echo 'PASS direct R execution, thread settings, combined logs and failure status'
# A named batch task is a function call; no configuration reaches commandArgs.
cat > code/scmultisim/named.R <<'RS'
run_fixture <- function() cat("named task fixture\n")
if (sys.nframe() == 0L) stop("Unexpected direct task entry")
RS
run_task code/scmultisim/named.R run_fixture > named.log 2>&1
grep -q 'named task fixture' named.log

if [[ -d /proc ]]; then
  cat > code/gse274113/fit_control_timepoints.R <<'RS'
stopifnot(Sys.getenv("OMP_NUM_THREADS")=="1")
p <- processx::process$new("sleep", "2")
cat("fixture inference\n"); Sys.sleep(2); p$wait()
RS
  run_task code/gse274113/fit_control_timepoints.R > monitored.log 2>&1
  for field in resource_block_started resource_producer_sha256 resource_wrapper_sha256 resource_sampler_sha256 \
    resource_host_fingerprint resource_command_started_epoch_ns resource_command_completed_epoch_ns resource_block_completed; do
    grep -q "^$field" monitored.log
  done
  awk -F '\t' '$1=="process_tree_rss_sum_peak_kb" {rss=$2} $1=="process_tree_processes_at_peak" {n=$2}
    $1=="process_tree_rss_samples" {samples=$2} END {exit !(rss>0 && n>=3 && samples>1)}' monitored.log
  grep -q $'resource_command_exit_status\t0' monitored.log
  grep -q 'Maximum resident set size' monitored.log

  # Python uses the same entry monitoring without a Python launcher process.
  cat > code/gse274113/run_pycistopic.py <<'PYTHON'
import os, subprocess, time
assert os.environ["OMP_NUM_THREADS"] == "1"
child = subprocess.Popen(["sleep", "2"])
print("fixture topic model", flush=True)
time.sleep(2)
child.wait()
PYTHON
  run_task code/gse274113/run_pycistopic.py > python_monitored.log 2>&1
  grep -q 'fixture topic model' python_monitored.log
  grep -q $'resource_sampler_exit_status\t0' python_monitored.log
  awk -F '\t' '$1=="process_tree_rss_sum_peak_kb" {rss=$2} $1=="process_tree_processes_at_peak" {n=$2}
    $1=="process_tree_rss_samples" {samples=$2} END {exit !(rss>0 && n>=3 && samples>1)}' python_monitored.log
  echo 'PASS direct Python task, shared logging, GNU time and process-tree RSS'
  # Failed monitored tasks still finish the provenance block and preserve status.
  printf 'quit(status=7L)\n' > code/gse274113/fit_control_timepoints.R
  if run_task code/gse274113/fit_control_timepoints.R > monitored_failure.log 2>&1; then
    exit 1
  else
    [[ "$?" == 7 ]]
  fi
  grep -q $'resource_command_exit_status\t7' monitored_failure.log
  grep -q '^resource_block_completed' monitored_failure.log
  ! grep -q stage_completed monitored_failure.log
  echo 'PASS Linux process-tree RSS, GNU time, provenance and monitored failure status'
  cat > code/gse274113/fit_control_timepoints.R <<'RS'
child <- parallel::mcparallel(Sys.sleep(30))
writeLines(as.character(c(Sys.getpid(), child$pid)), "worker_pids.txt")
Sys.sleep(30)
RS
  # Send the signal to the driver group, as a terminal cancellation does.
  setsid bash -c 'source code/functions/runtime.sh; run_task code/gse274113/fit_control_timepoints.R' > cancellation.log 2>&1 &
  driver=$!
  for attempt in {1..200}; do
    [[ ! -f worker_pids.txt ]] || break
    sleep 0.02
  done
  [[ -f worker_pids.txt ]]
  kill -TERM -- "-$driver"
  wait "$driver" && exit 1
  sleep 0.2
  while read -r pid; do
    if kill -0 "$pid" 2>/dev/null; then
      # A terminated orphan can briefly remain a zombie until init reaps it.
      state=$(ps -o stat= -p "$pid")
      [[ "$state" == Z* ]]
    fi
  done < worker_pids.txt
  echo 'PASS monitored task cancellation stops the R process and its forked worker'
fi
