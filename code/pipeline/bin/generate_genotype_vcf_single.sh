#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [--min_depth_of_coverage MIN_DEPTH_OF_COVERAGE]
  [--heterozygous_perc HETEROZYGOUS_PERC]
  [-s|--snps_list_file SNPS_LIST_FILE]
  output_path
  log_path
  analysis_id
  input_file
USAGE
}

min_depth_of_coverage=10
heterozygous_perc=0.2
snps_list_file="/dev/null"

required_args=4
declare -A longoptspec
longoptspec=( [min_depth_of_coverage]=1 [heterozygous_perc]=1 [snps_list_file]=1 )
optspec="s:"
set_args () {
  local handled=0
  case "${1}" in
    snps_list_file)
        snps_list_file="${2}"
        ;;
    min_depth_of_coverage)
        min_depth_of_coverage="${2}"
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

analysis_name='generate_genotype_vcf'

. $(dirname $0)/get_common_paths.sh


shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"

input_file="$4"

(
echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${analysis_id}";
awk -v htperc="${heterozygous_perc}" -v min_depth_of_coverage="${min_depth_of_coverage}" \
  -v snps_list_file="${snps_list_file}" '
BEGIN {
  OFS          = "\t";
  s            = ":";
  htperc_top   = 1 - htperc;
  print_all    = 1;
  fixed_fields = "." OFS "." OFS "." OFS "GT";
}
{
  if (FILENAME == snps_list_file && $0 !~ /^(#|CHR)/) {
    ids[$1 s $2 s $3 s $4 s $5]   = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS fixed_fields OFS "./.";
    print_all = 0;
  } else {
    key  = $fn["chr"] s $fn["pos"] s $fn["rsid"] s $fn["ref"] s $fn["alt"];
    if (FNR == 1) {
      for (i = 1; i <= NF; i++) {
        fn[$i] = i;
      }
    } else if (print_all || (key in ids && !(key in snp_counts))) {
      snp_counts[key]++;
      af = $fn["af"];
      gt = "./.";
      if ($fn["cov"] >= min_depth_of_coverage) {
        if (af < htperc) {
          gt = "0/0";
        } else if (af > htperc_top) {
          gt = "1/1";
        } else {
          gt = "0/1";
        }
      }
      print $fn["chr"], $fn["pos"], $fn["rsid"], $fn["ref"], $fn["alt"], fixed_fields, gt;
    }
  }
}
END {
  if (!print_all) {
    for (id in ids) {
      if (!(id in snp_counts)) {
        print ids[id];
      }
    }
  }
}' "${snps_list_file}" "$input_file" | sort -k1,1V -k2,2n
) > "${output_prefix}.vcf"

