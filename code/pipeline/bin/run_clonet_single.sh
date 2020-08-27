#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [-t|--threads NUM_THREADS]
  output_path
  log_path
  analysis_id
  normal_analysis_id
  normal_sample_file
  tumor_analysis_id
  tumor_sample_file
  informative_snps_dir
  sample_segments_file
USAGE
}

threads=1

required_args=9
declare -A longoptspec
longoptspec=( [threads]=1 )
optspec="t:"
set_args () {
  local handled=0
  case "${1}" in
    t|threads)
        threads="${2}"
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh


analysis_name='clonet'

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"

normal_sample_name="$4"
normal_sample="$(abs_path.py $5)"
tumor_sample_name="$6"
tumor_sample="$(abs_path.py $7)"

informative_snps_dir="$8"
export segments_file="$9"

export num_threads="${threads}"


cfg_dir="${project_dir}/cfg"
clonet_path="/usr/share/clonet/src"


cfg_dir="${project_dir}/cfg/${analysis_name}/$analysis_id"
[[ -d "$cfg_dir" ]] || mkdir -p "$cfg_dir"

export sample_info_file="${cfg_dir}/sample_info_file.tsv"
clonet_cfg="${cfg_dir}/config_file.R"

cat > "$sample_info_file" <<EOF
Tumor.Array.Name	Tumor.Bam.Name	Normal.Array.Name	Normal.Bam.Name
$tumor_sample_name	$tumor_sample	$normal_sample_name	$normal_sample
EOF

export clonet_functions_path="$clonet_path/"
export output_dir="${output_prefix}/"
export informative_snps_path="${informative_snps_dir}/"
export error_table_file="${tool_res_dir}/error_table.tsv"
export adm_method='2D'
export stages='1,2,3,4,5,6'

mo "${tool_res_dir}/cfg/clonet_config_file_template.tpl" \
  > "$clonet_cfg"

clonet "$clonet_cfg"            \
  >& "${log_prefix}.clonet.log"
