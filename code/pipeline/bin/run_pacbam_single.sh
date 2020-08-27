#!/usr/bin/env bash

usage() {
cat <<USAGE
usage: $0
  [-b|--bed BED_FILE]
  [-v|--vcf VCF_FILE]
  [-k|--kit_name KIT_NAME]
  [-d|--dbsnp_version DBSNP_VERSION]
  [-t|--threads THREADS]
  [--min_base_quality MIN_BASE_QUALITY]
  [--min_read_quality MIN_READ_QUALITY]
  [--region_perc REGION_PERC]
  output_path
  log_path
  analysis_id
  bam_file
  reference_genome_name
USAGE
}

min_read_quality=20
min_base_quality=20
region_perc=0.5
threads=1

required_args=5
declare -A longoptspec
longoptspec=( [bed]=1 [vcf]=1 [kit_name]=1 [dbsnp_version]=1 [threads]=1 [min_read_quality]=1 [min_base_quality]=1 [region_perc]=1 )
optspec="b:v:k:d:t:"
set_args () {
  local handled=0
  case "${1}" in
    b|bed)
        bed_file="${2}"
        ;;
    v|vcf)
        vcf_file="${2}"
        ;;
    k|kit_name)
        kit_name="${2}"
        ;;
    d|dbsnp_version)
        dbsnp_version="${2}"
        ;;
    t|threads)
        threads="${2}"
        ;;
    min_base_quality)
        min_base_quality="${2}"
        ;;
    min_read_quality)
        min_read_quality="${2}"
        ;;
    region_perc)
        region_perc="${2}"
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

analysis_name='pacbam'

reference_genome_name="${5}"

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path=${1%%+(/)}
log_path=${2%%+(/)}
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"

bam_file="${4}"
bed_file="${bed_file-${kit_target_bed}}"
vcf_file="${vcf_file-${dbsnp_restricted_kit_target_vcf}}"

pacbam                                    \
  bam="${bam_file}"                       \
  bed="${bed_file}"                       \
  vcf="${vcf_file}"                       \
  fasta="${reference_genome_fasta}"       \
  mode=2                                  \
  threads=${threads}                      \
  mbq=${min_base_quality}                 \
  mrq=${min_read_quality}                 \
  regionperc=${region_perc}               \
  out="${output_prefix}"                  \
  >& "${log_prefix}.${analysis_name}.log"

