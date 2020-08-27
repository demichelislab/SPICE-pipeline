#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [--downsample_to_coverage COVERAGE]
  output_path
  log_path
  analysis_id
  normal_sample_file
  tumor_sample_file
  reference_genome_name
  kit_name
  dbsnp_version
  cosmic_version
USAGE
}

downsample_to_coverage=10000

required_args=9
declare -A longoptspec
longoptspec=( [downsample_to_coverage]=1 )
set_args () {
  local handled=0
  case "${1}" in
    downsample_to_coverage)
      downsample_to_coverage="${2}"
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

analysis_name='mutect2'

reference_genome_name="$6"
kit_name="$7"
dbsnp_version="$8"
cosmic_version="$9"

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"
normal_sample_file="$4"
tumor_sample_file="$5"

  # --downsample_to_coverage ${downsample_to_coverage}                                                                                        \
gatk -T MuTect2                                                                                                                             \
  -U ALLOW_SEQ_DICT_INCOMPATIBILITY                                                                                                         \
  --read_group_black_list "${tool_res_dir}/read_group_blacklist.txt"                                                                        \
  --intervals "${kit_target_bed}"                                                                                                           \
  -I:normal "$normal_sample_file"                                                                                                           \
  -I:tumor "$tumor_sample_file"                                                                                                             \
  --reference_sequence "$reference_genome_fasta"                                                                                            \
  --dbsnp "${dbsnp_restricted_kit_target_vcf}"                                                                                              \
  --cosmic "${res_dir}/cosmic/${reference_genome_name}/${cosmic_version}/cosmic-coding_muts-${cosmic_version}-${reference_genome_name}.vcf" \
  --normal_panel "${tool_res_dir}/${reference_genome_name}/normal_panel-${reference_genome_name}.vcf"                                       \
  --out "${output_prefix}.vcf"                                                                                                              \
  --log_to_file "${log_prefix}.log"                                                                                                         \
  &> "${log_prefix}.cli.log"

