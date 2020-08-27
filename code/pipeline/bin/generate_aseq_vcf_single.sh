#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  output_path
  log_path
  analysis_id
  input_file
USAGE
}

required_args=4
. $(dirname $0)/parse_args_bash.sh

analysis_name='generate_aseq_vcf'

. $(dirname $0)/get_common_paths.sh


shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"

input_file="$4"

exec >& "$log_path"

(
echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tPASS";
awk -F'\t' -v OFS='\t' '$0 !~ /^#/ { print $1, $2, $3, $4, $5, ".", ".", "." }' "$input_file" |
sort -k1,1V -k2,2n |
uniq
) > "${output_prefix}.vcf"

