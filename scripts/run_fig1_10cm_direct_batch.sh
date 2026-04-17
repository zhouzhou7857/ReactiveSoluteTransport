#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
input_dir="$repo_root/Input"
output_dir="$repo_root/Output/fig1_10cm_direct_generated"
binary="$repo_root/Code/Release/ReactiveTransportPart"
file_names="$input_dir/File_names.txt"
backup_file="$(mktemp)"

cp "$file_names" "$backup_file"
trap 'cp "$backup_file" "$file_names"; rm -f "$backup_file"' EXIT

mkdir -p "$output_dir"

cases=(
  "9 1.5"
  "9 2.0"
  "9 2.5"
  "12 1.5"
  "12 2.0"
  "12 2.5"
  "16 1.5"
  "16 2.0"
  "16 2.5"
)

for case in "${cases[@]}"; do
  read -r p a <<<"$case"
  a_tag="${a/./p}"
  dfn_file="DFN_direct10cm_p${p}_a${a_tag}.txt"
  run_dir="$output_dir/p${p}_a${a_tag}"
  log_file="$run_dir/run.log"
  mkdir -p "$run_dir/input" "$run_dir/output"

  cat >"$file_names" <<EOF
Domain_10cm_square.txt
Simulation_dfn_only_fast.txt
$dfn_file
Chemistry_decay_reactivity_case_3.txt
EOF

  (
    cd "$repo_root/Code/Release"
    "$binary"
  ) >"$log_file" 2>&1

  cp "$repo_root/Input/Domain_files/Domain_10cm_square.txt" "$run_dir/input/"
  cp "$repo_root/Input/DFN_files/$dfn_file" "$run_dir/input/"
  cp "$repo_root/Output/DFN_raw.txt" "$run_dir/output/DFN_raw.txt"
  [[ -f "$repo_root/Output/DFN_init.txt" ]] && cp "$repo_root/Output/DFN_init.txt" "$run_dir/output/DFN_init.txt"
  [[ -f "$repo_root/Output/DFN.txt" ]] && cp "$repo_root/Output/DFN.txt" "$run_dir/output/DFN.txt"
done

echo "Saved direct-generated 10 cm outputs to $output_dir"
