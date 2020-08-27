#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  output_path
  log_path
  analysis_id
  bam_file
  reference_genome_name
  kit_name
USAGE
}


required_args=6
. $(dirname $0)/parse_args_bash.sh

analysis_name='picard_hsmetrics'

reference_genome_name="$5"
kit_name="$6"

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
bam_file="$4"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"

picard CollectHsMetrics                             \
  VALIDATION_STRINGENCY="LENIENT"                   \
  I="$bam_file"                                     \
  R="$reference_genome_fasta"                       \
  BAIT_INTERVALS="$kit_bait_interval_list"          \
  TARGET_INTERVALS="$kit_target_interval_list"      \
  O="${output_prefix}-hsmetrics.txt"                \
  PER_TARGET_COVERAGE="${output_prefix}-target.tsv" \
  &> "${log_prefix}.log"

