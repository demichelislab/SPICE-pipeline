#!/usr/bin/env bash

usage() {
cat <<USAGE
USAGE: $0
  [-s|--sample_name SAMPLE_NAME]
  output_path
  log_path
  analysis_id
  pacbam_snps
USAGE
}

required_args=4
declare -A longoptspec
longoptspec=( [sample_name]=1 )
optspec="s:"
set_args () {
  local handled=0
  case "${1}" in
    s|sample_name)
        sample_name="${2}"
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

analysis_name='generate_read_count'

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
pacbam_snps="${4}"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"

sample_name="${sample_name-$analysis_id}"

awk -v sample_name="$sample_name" '
BEGIN {
  OFS = "\t";
  print "Gene.id", "chr", "start", "end", "sample", "ref.count", "alt.count";
}{
  if (NR == 1) {
    for (i = 1; i <= NF; i++) {
      f[$i] = i;
    }
  } else {
    print "", $f["chr"], $f["pos"], $f["pos"], sample_name, $f[$f["ref"]], $f[$f["alt"]];
  }
}' "$pacbam_snps" \
  > "${output_prefix}.txt"

