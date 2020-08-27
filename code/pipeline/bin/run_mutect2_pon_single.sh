#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  output_path
  log_path
  analysis_id
  normal_sample_file
  reference_genome_name
  kit_name
  dbsnp_version
  cosmic_version
USAGE
}


required_args=8
. $(dirname $0)/parse_args_bash.sh

analysis_name='mutect2_pon'

reference_genome_name="$5"
kit_name="$6"
dbsnp_version="$7"
cosmic_version="$8"

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"
sample_file="$4"

  # --downsample_to_coverage ${downsample_to_coverage}                                                                                        \
gatk -T MuTect2                                                                                                                             \
  --artifact_detection_mode                                                                                                                 \
  --intervals "${kit_target_bed}"                                                                                                           \
  -I:tumor "$sample_file"                                                                                                                   \
  --reference_sequence "$reference_genome_fasta"                                                                                            \
  --dbsnp "${dbsnp_restricted_kit_target_vcf}"                                                                                              \
  --cosmic "${res_dir}/cosmic/${reference_genome_name}/${cosmic_version}/cosmic-coding_muts-${cosmic_version}-${reference_genome_name}.vcf" \
  --out "${output_prefix}.vcf"                                                                                                              \
  --log_to_file "${log_prefix}.log"                                                                                                         \
  &> "${log_prefix}.cli.log"

