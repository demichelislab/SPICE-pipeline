#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [-t|--threads NUM_THREADS]
  [--flag_pick]
  output_path
  log_path
  analysis_id
  input_file
  reference_genome_name
USAGE
}

threads=1
pick_option='--pick'

required_args=5
declare -A longoptspec
longoptspec=( [threads]=1 )
optspec="t:"
set_args () {
  local handled=0
  case "${1}" in
    t|threads)
      threads="${2}"
        ;;
    flag_pick)
      pick_option='--flag_pick'
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

analysis_name='vep'
reference_genome_name="$5"
. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
input_file="$4"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"

reference_genome_name_vep=$(
awk -F '\t' -v reference_genome_name="${reference_genome_name}" '
NR > 1 {
  m[$1] = $2;
}
END {
  if (reference_genome_name in m) {
    print m[reference_genome_name];
  } else {
    print "Cannot find", reference_genome_name, "in genome_reference_name_map.tsv" > "/dev/stderr";
    exit 1;
  }
}
' ${tool_res_dir}/genome_reference_name_map.tsv)

vep                                        \
 --input_file "$input_file"                \
 --format 'vcf'                            \
 --output_file "${output_prefix}.tsv"      \
 --assembly "${reference_genome_name_vep}" \
 --dir "${tool_res_dir}/vep_data"          \
 --everything                              \
 --fork "$threads"                         \
 --tab                                     \
 --offline                                 \
 --dont_skip                               \
 "$pick_option"                            \
 --stats_text                              \
 --warning_file "${log_prefix}.log"        \
 &> "${log_prefix}.cli.log"


