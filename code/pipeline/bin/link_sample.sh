#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  sample_path
  destination_path
USAGE
}

required_args=2
. $(dirname $0)/parse_args_bash.sh

sample="$1"
sample_name=$(sed -re 's/\.bam$//' <<<"$sample")
destination="$2"
destination_name=$(sed -re 's/\.bam$//' <<<"$destination")

have_errors=0

for ext in .bam .bam.bai .bai
do
  cur_sample=$(abs_path.py "${sample_name}${ext}")
  if [[ -e "$cur_sample" ]]
  then
    ln -sf "$cur_sample" "$destination_name$ext"
  elif [[ -z "$ext" ]]
  then
    echo "ERROR: Sample $cur_sample doesn't exist"
    have_errors=1
  fi
done

exit ${have_errors}

