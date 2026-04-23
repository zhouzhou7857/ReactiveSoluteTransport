#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
input_dir="$repo_root/Input"
output_root="$repo_root/Output/song_selected_10cm_rtm_vp5e-6_np1e5_dp100Pa"
binary="$repo_root/Code/Release/ReactiveTransportPart"

mkdir -p \
  "$output_root/generated_inputs/Domain_files" \
  "$output_root/generated_inputs/Simulation_files" \
  "$output_root/file_names"

pressure_pa="100.0"
rho="1000.0"
g="9.8"
left_head="$(awk -v p="$pressure_pa" -v rho="$rho" -v g="$g" 'BEGIN{printf "%.15e", p/(rho*g)}')"

domain_name="song_selected_10cm_dp100Pa_domain.txt"
domain_path="$input_dir/Domain_files/$domain_name"
cat >"$domain_path" <<EOF
0.1 0.1
1e-10 0.05
$left_head 0.0
EOF
cp "$domain_path" "$output_root/generated_inputs/Domain_files/$domain_name"

simulation_name="song_selected_10cm_rtm_vp5e-6_np1e5.txt"
simulation_path="$input_dir/Simulation_files/$simulation_name"
cat >"$simulation_path" <<'EOF'
100000
0.05
0
1e-1
1e6
200
1
1e6
1e5
100

5.000000000000e-06
1.000000000000e+00
EOF
cp "$simulation_path" "$output_root/generated_inputs/Simulation_files/$simulation_name"

cases=(
  "p9_a1p5 song_selected_10cm/p9_a1p5/p9_a1p5_raw_filemode.txt"
  "p12_a2p0 song_selected_10cm/p12_a2p0/p12_a2p0_raw_filemode.txt"
  "p16_a2p5 song_selected_10cm/p16_a2p5/p16_a2p5_raw_filemode.txt"
)

clean_main_output() {
  local pattern
  for pattern in \
    "DFN_step*.txt" \
    "particle_positions_t*.csv" \
    "DFN.txt" \
    "DFN_aperture_delta.txt" \
    "DFN_init.txt" \
    "DFN_raw.txt" \
    "cdf.txt" \
    "pdf.txt"; do
    find "$repo_root/Output" -maxdepth 1 -type f -name "$pattern" -delete
  done
}

copy_case_outputs() {
  local target_dir="$1"
  local pattern
  mkdir -p "$target_dir"
  for pattern in \
    "DFN_step*.txt" \
    "particle_positions_t*.csv" \
    "DFN.txt" \
    "DFN_aperture_delta.txt" \
    "DFN_init.txt" \
    "DFN_raw.txt" \
    "cdf.txt" \
    "pdf.txt"; do
    find "$repo_root/Output" -maxdepth 1 -type f -name "$pattern" -exec cp {} "$target_dir/" \;
  done
}

summary_csv="$output_root/summary.csv"
echo "case,status,log_file" >"$summary_csv"

for entry in "${cases[@]}"; do
  read -r case_name dfn_name <<<"$entry"
  case_dir="$output_root/$case_name"
  input_snapshot_dir="$case_dir/input"
  output_snapshot_dir="$case_dir/output"
  log_file="$case_dir/run.log"
  file_names_path="$output_root/file_names/${case_name}.txt"

  rm -rf "$case_dir"
  mkdir -p "$input_snapshot_dir" "$output_snapshot_dir"

  cat >"$file_names_path" <<EOF
$domain_name
$simulation_name
$dfn_name
EOF

  cp "$domain_path" "$input_snapshot_dir/"
  cp "$simulation_path" "$input_snapshot_dir/"
  cp "$input_dir/DFN_files/$dfn_name" "$input_snapshot_dir/"
  cp "$file_names_path" "$input_snapshot_dir/File_names.txt"

  clean_main_output

  (
    cd "$repo_root/Code/Release"
    RST_USE_VP_WIDTH_CORRECTION=1 "$binary" "$file_names_path"
  ) >"$log_file" 2>&1

  copy_case_outputs "$output_snapshot_dir"
  echo "$case_name,completed,$log_file" >>"$summary_csv"
done

echo "Saved batch inputs and outputs to $output_root"
