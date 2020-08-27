#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [--cfg_file_name CFG_FILE_NAME]
  input_file
  cfg_base_dir
  output_file
USAGE
}

cfg_file_name="config_file.R"

required_args=3
declare -A longoptspec
longoptspec=( [cfg_file_name]=1 )
optspec=""
set_args () {
  local handled=0
  case "${1}" in
    cfg_file_name)
        cfg_file_name="${2}"
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

input_file="$1"
cfg_base_dir="$2"
output_file="$3"

cfg_abs_dir="$( abs_path.py ${cfg_base_dir} )"

(echo -e "sampleName\tpath";
awk -v OFS="\t" -v cfg_dir="${cfg_abs_dir}" \
    -v cfg_file_name="${cfg_file_name}" '
BEGIN {
  sep = "/";
}
{
  if (NR == 1) {
    for (i = 1; i <= NF; i++) {
      fn[$i] = i;
    }
  } else {
    id = $fn["analysis_id"];
    print id, cfg_dir sep id sep cfg_file_name;
  }
}' "$input_file" |
   sort -k1,1 ) \
  > "$output_file"
