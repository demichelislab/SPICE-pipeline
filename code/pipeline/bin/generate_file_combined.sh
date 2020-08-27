#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [--num_header_rows NUM_HEADER_ROWS]
  [--fields_to_output FIELDS]
  [--header_line HEADER_LINE]
  output_path
  log_path
  input_folder
  file_name_pattern
  output_file_name
USAGE
}

num_header_rows=1
fields_to_output=""
header_line=0

required_args=5
declare -A longoptspec
longoptspec=( [num_header_rows]=1 [fields_to_output]=1 [header_line]=1 )
set_args () {
  local handled=0
  case "${1}" in
    num_header_rows)
        num_header_rows="${2}"
        ;;
    fields_to_output)
        fields_to_output="${2}"
        ;;
    header_line)
        header_line="${2}"
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

analysis_name='generate_file_combined'

. $(dirname $0)/get_common_paths.sh


shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
output_prefix="$output_path"
log_prefix="$log_path"
input_folder="${3%%+(/)}"
file_name_pattern="$4"
output_file_name="$5"

exec &> "$log_path"

input_files=( "$(find $input_folder -iname $file_name_pattern)" )

if [[ ${#input_files[@]} -gt 0 ]]
then

  awk -F'\t' \
    -v fields="$fields_to_output"
    -v num_header_rows="$num_header_rows"
    -v header_line="$header_line" '
  function print_selected () {
    out = \$0;
    if (fields) {
      out = \"\";
      for (i = 1; i <= f_len; i++) {
        if (!fields_check) {
          cur_f = f[i];
          if (cur_f < 1 || cur_f > NF) {
            printf(\"ERROR: %i is outside of the field range [1, %i]\", cur_f, NF);
            exit(1);
          }
        }
        out = out\$f[i];
        if (i < f_len) {
          out = out\"\t\";
        } else if (!fields_check) {
            fields_check = 1;
        }
      }
    }
    print out;
  }
  BEGIN {
    split(fields, f, \",\");
    f_len = length(f);
    fields_check = 0;
  }
  {
    if (FNR > num_header_rows) {
      print_selected();
    }
    if (NR <= num_header_rows) {
      if (FNR == header_line) {
        print_selected();
      } else {
        print;
      }
    }
  }' ${input_files[@]} \
    > "$output_prefix/$output_file_name"

else
  echo "Cannot find any file matching the pattern ${file_name_pattern} in the path " \
       "${output_path}"
  exit 1
fi
