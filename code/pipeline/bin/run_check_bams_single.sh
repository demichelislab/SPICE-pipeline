#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  [-n|--no_chr]
  output_path
  log_path
  analysis_id
  bam_file
  reference_genome_name
USAGE
}

no_chr=0

required_args=5
declare -A longoptspec
longoptspec=( [no_chr]=0 )
optspec="n"
set_args () {
  local handled=0
  case "${1}" in
    n|no_chr)
      no_chr=1
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

analysis_name='check_bams'

reference_genome_name="$5"

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
bam_file="$4"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"

chromosome_size_file="${output_prefix}-chromosome_size.tsv"

exec >& "${log_path}.log"

samtools view -H ${bam_file} |
  awk -F'\t' -v OFS='\t' '
BEGIN {
  label_re = "^[^:]+:";
  print "chromosome", "size";
}
/^@SQ/ {
  sub(label_re, "", $2);
  sub(label_re, "", $3);
  print $2, $3;
}' > "${chromosome_size_file}"

awk -F '\t' -v no_chr=$no_chr '
BEGIN {
  pref = "";
  if (!no_chr) {
    pref = "chr";
  }
  for (i = 1; i < 23) {
    main_chromosomes[i] = 1;
  }
  main_chromosomes["X"] = 1;
  main_chromosomes["Y"] = 1;
  for (i = 1; i < ARGC; i++) {
    AP[ARGV[i]] = i;
  }
}
AP[$FILENAME] == 1 {
  chr_sizes[] = ;
}
AP[$FILENAME] == 2 {
}
' "$reference_genome_chromosome_size" "$chromosome_size_file"

# compare chr sizes
# check if chr/not
# chr option
