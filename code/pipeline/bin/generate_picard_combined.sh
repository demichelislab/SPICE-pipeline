#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [--postfix_to_remove POSTFIX_TO_REMOVE]
  output_path
  log_path
  input_folder
  file_name_pattern
  output_file_name
USAGE
}

postfix_to_remove=''

required_args=5
declare -A longoptspec
longoptspec=( [postfix_to_remove]=1 )
set_args () {
  local handled=0
  case "${1}" in
    postfix_to_remove)
        postfix_to_remove="${2}"
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

project_dir="$( abs_path.py $(dirname $0)/../ )"
dbs_dir="${project_dir}/dbs"
res_dir="${project_dir}/resources"
analysis_name='generate_file_combined'

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
output_prefix="$output_path"
input_folder="${3%%+(/)}"
file_name_pattern="$4"
output_file_name="$5"

exec >& "$log_path"

input_files=( $(find "$input_folder" -iname "$file_name_pattern" | sort) )

echo "${input_files[@]}"

if [[ ${#input_files[@]} -gt 0 ]]
then
  awk -v postfix_to_remove="$postfix_to_remove" -f <( cat <<-'EOP'
                     FNR == 1 {
                       sample_name = gensub("^(.*/)?([^/]+)"postfix_to_remove"\\..+$", "\\2", "g", FILENAME);
                       files++;
                     }
                     out_section {
                       out_sample = (files == 1 && out_section == 1) ?  "SAMPLE" : sample_name;
                       if (!/^\s*$/) {
                         printf "%s\t%s\n", out_sample, $0;
                         out_section++;
                       } else {
                         out_section = 0;
                       }
                     }
                     /METRICS CLASS\tpicard.analysis.directed.HsMetrics/ {
                       out_section = 1;
                       if (files != 1) {
                         getline;
                       }
                     }
EOP
) "${input_files[@]}" > "$output_prefix/$output_file_name"
else
  echo "Cannot find any file matching the pattern ${file_name_pattern} in the path " \
       "${output_path}"
  exit 1
fi
