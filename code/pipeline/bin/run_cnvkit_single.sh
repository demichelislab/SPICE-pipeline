#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [-s|--sex m|f]
  output_path
  log_path
  analysis_id
  normal_analysis_id
  normal_sample_file
  tumor_analysis_id
  tumor_sample_file
  reference_genome_name
  kit_name
USAGE
}

sex="auto"

required_args=9
declare -A longoptspec
longoptspec=( [sex]=1 )
optspec="s:"
set_args () {
  local handled=0
  case "${1}" in
    s|sex)
        sex="${2}"
        if [[ "$sex" != 'm' && "${sex}" != 'f' ]]
        then
          echo "Sex option set to '${sex}'. Valid values are 'm' or 'f'"
          exit 1
        fi
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

analysis_name='cnvkit'

reference_genome_name="$8"
kit_name="$9"

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"
normal_sample_name="$4"
normal_sample_file="$5"
tumor_sample_name="$6"
tumor_sample_file="$7"

sex_option=''

if [[ "$sex" != 'auto' ]]
then
  sex_option=( )
  if [[ "$sex" == 'm' ]]
  then
    sex_option+=( '--male-reference' )
  fi
fi


cnvkit batch "$tumor_sample_file"                                                                  \
  "${sex_option[@]}"                                                                               \
  --normal "$normal_sample_file"                                                                   \
  --targets "${kit_bait_bed}"                                                                      \
  --annotate "${tool_res_dir}/${reference_genome_name}/ucsc_ref_flat-${reference_genome_name}.txt" \
  --fasta "${reference_genome_fasta}"                                                              \
  --access "${tool_res_dir}/${reference_genome_name}/access_5kb-${reference_genome_name}.bed"      \
  --output-dir "${output_prefix}/"                                                                 \
  >& "${log_prefix}.log"

