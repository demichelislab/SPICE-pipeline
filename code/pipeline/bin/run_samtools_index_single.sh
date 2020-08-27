#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [-f|--force_rebuild]
  output_path
  log_path
  analysis_id
  bam_file
USAGE
}

force_rebuild=0

required_args=4
declare -A longoptspec
longoptspec=( [force_rebuild]=0 )
optspec="f"
set_args () {
  local handled=0
  case "${1}" in
    f|force_rebuild)
        force_rebuild=1
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

analysis_name='samtools_index'

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"

bam_file="$4"
bam_filename=$(sed -re 's/\.bam$//' <<<"$bam_file")

exec &> "${log_prefix}.${analysis_name}.log"


if [[ (! ( -e "${bam_filename}.bam.bai" ||  -e "${bam_filename}.bai")) || (${force_rebuild} -eq 1) ]]
then
  samtools index \
    "${bam_file}"
else
  echo "Index exists for sample ${bam_file} and force_rebuild is not specified. Exiting"
fi
