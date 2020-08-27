#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  input_file
  output_file
USAGE
}

required_args=2
. $(dirname $0)/parse_args_bash.sh

input_file="$1"
output_file="$2"

(echo -e "pair_id\tanalysis_id\tfile\tsex\treference_genome\tkit_name\tdbsnp_version\tcosmic_version\ttype";
awk '{
  if (NR == 1) {
    OFS = "\t";
    for (i = 1; i <= NF; i++) {
      head[$i] = i
    }
  } else {
    print $head["analysis_id"],
      $head["analysis_id"]"-n",
      $head["file_normal"],
      $head["sex"],
      $head["reference_genome"],
      $head["kit_name"],
      $head["dbsnp_version"],
      $head["cosmic_version"],
      "normal";
    print  $head["analysis_id"],
      $head["analysis_id"],
      $head["file_tumor"],
      $head["sex"],
      $head["reference_genome"],
      $head["kit_name"],
      $head["dbsnp_version"],
      $head["cosmic_version"],
      "tumor";
  }
}' "$input_file" |
   sort |
   uniq ) \
  > "$output_file"
