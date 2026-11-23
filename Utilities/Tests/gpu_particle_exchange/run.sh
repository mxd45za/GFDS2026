#!/usr/bin/env bash
set -euo pipefail

# Runs serially. Each result directory is new, so prior baselines are never overwritten.
test_source=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
fds_binary=$(realpath -- "${1:?Usage: run.sh /path/to/fds [new-result-directory] [CPU-reference-root]}")
reference_root=${3:-}
if [[ $# -ge 2 ]]; then
    mkdir -- "$2"
    test_output=$(realpath -- "$2")
else
    test_output=$(mktemp -d -t gfds-particle-exchange.XXXXXXXX)
fi
printf 'Results: %s\n' "$test_output"

for case_name in one_way bidirectional reinsertion removal periodic linked many_particles multidest count_devices; do
    ranks=2
    [[ $case_name != bidirectional ]] || ranks=3
    for mode in C GC G1 G; do
        result_dir="$test_output/${case_name}_${mode}"
        if [[ $mode == C && -n $reference_root ]] &&
            cmp -s -- "$test_source/$case_name.fds" "$reference_root/${case_name}_C/$case_name.fds"; then
            ln -s -- "$(realpath -- "$reference_root/${case_name}_C")" "$result_dir"
            python3 "$test_source/check_particles.py" "$result_dir" --case "$case_name" > /dev/null
            printf '%s C: reused verified baseline\n' "$case_name"
            continue
        fi
        mkdir -- "$result_dir"
        cp -- "$test_source/$case_name.fds" "$result_dir/"
        mode_flag=$mode
        run_ranks=$ranks
        if [[ $mode == G1 ]]; then
            mode_flag=G
            run_ranks=1
        fi
        (
            cd -- "$result_dir"
            export OMP_NUM_THREADS=1
            if [[ $mode == GC ]]; then export CUDA_VISIBLE_DEVICES=-1; fi
            mpirun -np "$run_ranks" "$fds_binary" "$case_name.fds" "--$mode_flag" > run.log 2>&1
        )
        check_args=("$result_dir" --case "$case_name")
        [[ $mode == C ]] || check_args+=(--reference "$test_output/${case_name}_C")
        python3 "$test_source/check_particles.py" "${check_args[@]}" > "$result_dir/check.json"
        printf '%s %s: passed\n' "$case_name" "$mode"
    done
done
