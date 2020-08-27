#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  output_path
  log_path
  analysis_id
  normal_genotype_file
  tumor_genotype_file
USAGE
}

required_args=5
. $(dirname $0)/parse_args_bash.sh

analysis_name='spia'

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"

[[ -d "$output_prefix" ]] || mkdir -p "$output_prefix"

normal_genotype="$(abs_path.py $4)"
tumor_genotype="$(abs_path.py $5)"

export spia_functions_path="$(dirname $(abs_path.py $(which SPIA.R)))/SPIAfunctions.R"

cfg_dir="${project_dir}/cfg/${analysis_name}/$analysis_id"
[[ -d "$cfg_dir" ]] || mkdir -p "$cfg_dir"

export vcf_files_list="${cfg_dir}/sample_info_file.tsv"
spia_cfg="${cfg_dir}/config_file.R"

cat > "$vcf_files_list" <<EOF
$tumor_genotype
$normal_genotype
EOF

export output_spia_file="${output_prefix}/spia.tsv"
export save_spia_plot=T
export output_plot_file="${output_prefix}/spia_plot.pdf"
export save_genotype=T

mo "${tool_res_dir}/cfg/spia_config_file_template.tpl" > "$spia_cfg"

spia "$spia_cfg" &> "${log_prefix}.spia.log"

