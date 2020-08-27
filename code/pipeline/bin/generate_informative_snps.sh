#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [-m|--min_depth_of_coverage MIN_DEPTH_OF_COVERAGE]
  [--min_depth_of_coverage_tumor MIN_DEPTH_OF_COVERAGE_TUMOR]
  [--heterozygous_perc HETEROZYGOUS_PERC]
  output_path
  log_path
  analysis_id
  normal_snps
  tumor_snps
USAGE
}

min_depth_of_coverage=20
heterozygous_perc=0.2

required_args=5
declare -A longoptspec
longoptspec=( [min_depth_of_coverage]=1 [min_depth_of_coverage_tumor]=1 [heterozygous_perc]=1 )
optspec="m"
set_args () {
  local handled=0
  case "${1}" in
    m|min_depth_of_coverage)
        min_depth_of_coverage="${2}"
        ;;
    min_depth_of_coverage_tumor)
        min_depth_of_coverage_tumor="${2}"
        ;;
    heterozygous_perc)
        heterozygous_perc="${2}"
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

analysis_name='generate_informative_snps'

. $(dirname $0)/get_common_paths.sh


shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"

normal_snps="$4"
tumor_snps="$5"

min_depth_of_coverage_tumor=${min_depth_of_coverage_tumor-${min_depth_of_coverage}}

[[ -e $output_prefix ]] || mkdir -p $output_prefix/{Normal,Tumor}Pileup

paste -d'\n' ${normal_snps} ${tumor_snps} |
awk -v OFS=","                                                      \
    -v heterozygous_perc="${heterozygous_perc}"                     \
    -v output_prefix="${output_prefix}"                             \
    -v analysis_id="${analysis_id}"                                 \
    -v min_depth_of_coverage="${min_depth_of_coverage}"             \
    -v min_depth_of_coverage_tumor="${min_depth_of_coverage_tumor}" \
    'BEGIN {
      heterozygous_perc_top = 1 - heterozygous_perc;
      norm_file  = output_prefix"/NormalPileup/"analysis_id"-n.snps";
      tumor_file = output_prefix"/TumorPileup/"analysis_id".snps";
      printf("") > norm_file;
      printf("") > tumor_file;
    }{
      if (NR == 1) {
        for (i = 1; i <= NF; i++) {
          fn[$i] = i;
        }
      }
      if ((length($fn["ref"]) == 1 && length($fn["alt"]) == 1) &&
          (NR % 2 == 1 && $fn["af"] >= heterozygous_perc && $fn["af"] <= heterozygous_perc_top) || NR == 1) {
        norm_snp_id = "";
        if ($fn["cov"] >= min_depth_of_coverage) {
          print >> norm_file;
          norm_snp_id = $fn["rsid"];
        }
        getline;
        if (norm_snp_id && $fn["cov"] >= min_depth_of_coverage_tumor) {
          if (norm_snp_id != $fn["rsid"]) {
            print "Tumor snp "$fn["rsid"]" is not present or not in the same order as in normal snps";
            exit 1;
          }
          print >> tumor_file;
        }
      }
    }' \
&> "${log_prefix}.${analysis_name}.log"
