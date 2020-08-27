#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [--include_indels]
  output_path
  log_path
  analysis_id
  input_file
USAGE
}

include_indels=0

required_args=4
declare -A longoptspec
longoptspec=( [include_indels]=0 )
set_args () {
  local handled=0
  case "${1}" in
    include_indels)
      include_indels=1
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

analysis_name='generate_bed_from_vcf'

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"

input_file="$4"

exec >& "${log_prefix}.log"

awk -F'\t' -v OFS='\t' -v include_indels=${include_indels} '
$0 !~ /^#/ && (include_indels || (tolower($4) ~ /^[acgt]$/ && tolower($5) ~ /^[acgt]$/)) {
  print $1, $2 - 1, $2;
}' "$input_file" |
  uniq \
> "${output_prefix}.bed"

