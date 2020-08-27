#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [-n|--normal_cov NORMAL_COV]
  [-t|--tumor_cov TUMOR_COV]
  output_path
  log_path
  analysis_id
  normal_snps
  tumor_snps
USAGE
}

normal_cov=10

required_args=5
declare -A longoptspec
longoptspec=( [normal_cov]=1 [tumor_cov]=1 )
optspec="n:t:"
set_args () {
  local handled=0
  case "${1}" in
    n|normal_cov)
        normal_cov="${2}"
        ;;
    t|tumor_cov)
        tumor_cov="${2}"
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

analysis_name='generate_facets_input'

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"

normal_snps="$4"
tumor_snps="$5"

tumor_cov="${tumor_cov-$normal_cov}"

(
echo -e "Chromosome,Position,Ref,Alt,File1R,File1A,File1E,File1D,File2R,File2A,File2E,File2D";
paste ${normal_snps} ${tumor_snps} |
awk -v OFS="," -v normal_cov="${normal_cov}" -v tumor_cov="${tumor_cov}" '{
       if (NR == 1) {
         st = "_n";
         for (i = 1; i <= NF; i++) {
           fn[$i st] = i;
           if (i == (NF / 2)) {
             st = "_t"
           }
         }
       } else {
         if ((length($fn["ref_n"]) == 1 && length($fn["alt_n"]) == 1 && $fn["alt_n"] != ".") &&
              $fn["cov_n"] >= normal_cov && $fn["cov_t"] >= tumor_cov) {
           tot_n = $fn["A_n"] + $fn["C_n"] + $fn["G_n"] + $fn["T_n"];
           tot_t = $fn["A_t"] + $fn["C_t"] + $fn["G_t"] + $fn["T_t"];
           ref_n = $fn[$fn["ref_n"]"_n"];
           alt_n = $fn[$fn["alt_n"]"_n"];
           ref_t = $fn[$fn["ref_n"]"_t"];
           alt_t = $fn[$fn["alt_n"]"_t"];
           print $fn["chr_n"], $fn["pos_n"], $fn["ref_n"], $fn["alt_n"],
                 ref_n, alt_n, tot_n - (ref_n + alt_n), 0,
                 ref_t, alt_t, tot_t - (ref_t + alt_t), 0;
         }
       }
     }'
) > "${output_prefix}.txt"

