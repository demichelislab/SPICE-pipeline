#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [-m|--method METHOD]
  [--drop_low_coverage]
  output_path
  log_path
  analysis_id
  cnvkit_coverage_file
USAGE
}

drop_low_coverage=0
method='cbs'

required_args=4
declare -A longoptspec
longoptspec=( [method]=1 [drop_low_coverage]=0 )
optspec="m:"
set_args () {
  local handled=0
  case "${1}" in
    m|method)
        method="${2}"
        ;;
    drop_low_coverage)
      drop_low_coverage=1
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

analysis_name='cnvkit_segment'

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"
cnvkit_coverage_file="$4"

dlc_option=''

if (( drop_low_coverage ))
then
  dlc_option="--drop-low-coverage"
fi


cnvkit segment              \
  ${dlc_option}             \
  --method "${method}"      \
  ${cnvkit_coverage_file}   \
  -o "${output_prefix}.cns" \
  >& "${log_prefix}.log"

