#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  output_path
  log_path
  analysis_id
  cnvkit_segment_file
USAGE
}

required_args=4
. $(dirname $0)/parse_args_bash.sh

analysis_name='cnvkit_export_segfile'

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"
cnvkit_segment_file="$4"

cnvkit export seg           \
  ${cnvkit_segment_file}    \
  -o "${output_prefix}.seg" \
  >& "${log_prefix}.log"

