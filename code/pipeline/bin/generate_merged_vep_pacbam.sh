#!/usr/bin/env bash

usage() {
cat <<USAGE
USAGE: $0
  output_path
  log_path
  analysis_id
  vep_output
  pacbam_normal_output
  pacbam_tumor_output
USAGE
}

required_args=6
. $(dirname $0)/parse_args_bash.sh

analysis_name='generate_merged_vep_pacbam'

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
analysis_id="$3"
output_prefix="${output_path}/$analysis_id"
log_prefix="${log_path}/$analysis_id"
vep_output="$4"
pacbam_normal_output="$5"
pacbam_tumor_output="$6"
output_file_path="${output_prefix}.txt"

exec &> "${log_prefix}.log"

Rscript <(
cat <<'EOF'
  cmd_args              <- commandArgs(TRUE)
  vep_output            <- cmd_args[[1]]
  pacbam_normal_output  <- cmd_args[[2]]
  pacbam_tumor_output   <- cmd_args[[3]]
  output_file_path      <- cmd_args[[4]]
  pacbam_cmp_fields     <- c('chr', 'pos', 'alt')
  pacbam_extract_fields <- c('rc_ref', 'rc_alt', 'af', 'cov')
  vep_cmp_fields        <- c('Location', 'Allele')

  vep_lines      <- readLines(vep_output)
  vep_is_header  <- grepl('^#', vep_lines)
  vep_header     <- vep_lines[vep_is_header]
  vep_content    <- vep_lines[!vep_is_header]
  vep_colnames   <- scan(text = sub('^#', '', tail(vep_header, 1)), what = 'chr',
                         quiet = TRUE)
  n_fields       <- length(gregexpr('\t', vep_content[[1]], fixed = TRUE)[[1]]) + 1
  vep <- local({
    dd <- scan(text = vep_content, what = setNames(as.list(rep('char', n_fields)),
               vep_colnames), quiet = TRUE)
    as.data.frame(dd, stringsAsFactors = FALSE)
  })

  prepare_keys   <- function (x) apply(x, 1, function (dd) paste(dd, collapse = ':'))
  extract_counts <- function (x, cn) mapply(function (ii, cc) getElement(x[ii, ], cc),
                                                  seq_len(nrow(x)), getElement(x, cn))
  read_pacbam_add_cols <- function (fn) {
    dd <- read.table(fn, sep = '\t', quote = '', header = TRUE, comment.char = '',
                     stringsAsFactors = FALSE, colClasses = 'character')
    cbind(dd,
          rc_ref = as.integer(extract_counts(dd, 'ref')),
          rc_alt = as.integer(extract_counts(dd, 'alt')),
          key    = prepare_keys(dd[, pacbam_cmp_fields]),
          stringsAsFactors = FALSE)
  }

  pacbam_normal  <- read_pacbam_add_cols(pacbam_normal_output)
  pacbam_tumor   <- read_pacbam_add_cols(pacbam_tumor_output)

  vep_keys       <- prepare_keys(vep[, vep_cmp_fields])

  normal_cols    <- setNames(pacbam_normal[match(vep_keys, pacbam_normal$key), pacbam_extract_fields],
                             paste0(pacbam_extract_fields, '_normal'))
  tumor_cols     <- setNames(pacbam_tumor[match(vep_keys, pacbam_tumor$key), pacbam_extract_fields],
                             paste0(pacbam_extract_fields, '_tumor'))

  merged         <- cbind(vep, normal_cols, tumor_cols, stringsAsFactors = FALSE)

  writeLines(c(head(vep_header, -1), paste0('#', paste0(colnames(merged), collapse = '\t'))), output_file_path)
  suppressWarnings(write.table(merged, output_file_path, sep = '\t', quote = FALSE,
                               col.names = FALSE, row.names = FALSE, na = '-',
                               append = TRUE))
EOF
) "$vep_output" "$pacbam_normal_output" "$pacbam_tumor_output" "$output_file_path"
