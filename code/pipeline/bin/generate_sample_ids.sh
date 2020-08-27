#!/usr/bin/env bash

usage() {
cat <<USAGE
USAGE: $0
  [--id_separator ID_SEPARATOR]
  [--id_internal_separator ID_INTERNAL_SEPARATOR]
  input_file
  output_file
  log_file
USAGE
}

id_separator="-"
id_internal_separator="_"

required_args=3
declare -A longoptspec
longoptspec=( [id_separator]=1 [id_internal_separator]=1 )
set_args () {
  local handled=0
  case "${1}" in
    id_separator)
        id_separator="${2}"
        ;;
    id_internal_separator)
        id_internal_separator="${2}"
        ;;
    *)
      handled=1
        ;;
  esac
  return ${handled}
}
. $(dirname $0)/parse_args_bash.sh

analysis_name='generate_sample_info_files'

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
input_file="$1"
output_path=$2
log_path="$3"

exec &> "$log_path"

Rscript <(
cat <<'EOF'
  cmd_args              <- commandArgs(TRUE)
  input_file            <- cmd_args[[1]]
  id_separator          <- cmd_args[[2]]
  id_internal_separator <- cmd_args[[3]]
  output_path           <- cmd_args[[4]]

  check_cols     <- function (tc, cc) length(intersect(tc, cc)) == length(cc)
  normalize      <- Vectorize(function (nn, sep = '_') {
                              regex <- paste('^', '+|', '+$', sep = sep)
                              gsub(regex, '', gsub('[^[:alnum:]]+', sep, nn))
                            }, vectorize.args = 'nn')
  per_line_list  <- function (ll, indent = 1) {
    format <- sprintf('%s%%s\n', do.call(paste0, as.list(rep('\t', indent))))
    do.call(paste0, lapply(ll, function (l) sprintf(format, l)))
  }
  paste_df_lines <- function (dd) apply(dd, 1, function (ll) do.call(paste, c(as.list(ll), sep = '\t')))
  sif            <- read.table(input_file, sep = '\t', quote = '', header = TRUE,
                               stringsAsFactors = F, comment.char = '')
  cols           <- colnames(sif)
  col_types      <- c('tumor', 'normal')
  file_cols      <- paste('file', col_types, sep = '_')
  id_cols        <- paste('analysis_id', col_types, sep = '_')
  have_file      <- check_cols(cols, file_cols)
  have_ids       <- check_cols(cols, id_cols)
  have_id        <- check_cols(cols, 'analysis_id')
  cols_to_write  <- c('analysis_id', file_cols, 'sex', 'reference_genome', 'kit_name', 'dbsnp_version', 'cosmic_version')
  id_sep         <- id_separator
  id_int_sep     <- id_internal_separator
  bad_samples    <- intersect(sif$file_normal, sif$file_tumor)
  dup_rows       <- sif[duplicated(sif), ]
  if (!have_file) {
    stop(sprintf('%s columns are needed', paste(file_cols)))
  }
  if (length(bad_samples)) {
    stop(sprintf('the following samples are specified as both normal and tumor:\n%s',
                 per_line_list(bad_samples)))
  }
  if (nrow(dup_rows)) {
    dup_rows_str <- paste_df_lines(dup_rows)
    warning(sprintf('the following rows in the input file are duplicated and will be collapsed:\n%s',
                    per_line_list(dup_rows_str)))
    sif <- unique(sif)
  }
  for (fc in file_cols) {
    cc <- getElement(sif, fc)
    dup <- duplicated(cc)
    if (any(dup)) {
      warning(sprintf('the following samples are duplicated in the %s column:\n%s',
              fc, per_line_list(cc[dup])))
    }
  }
  if (!have_id) {
    if (!have_ids) {
      id_normal <- gsub('\\.bam$', '', basename(sif$file_normal))
      id_tumor  <- gsub('\\.bam$', '', basename(sif$file_tumor))
    } else {
      id_normal <- sif$analysis_id_normal
      id_tumor  <- sif$analysis_id_tumor
    }
    sif$analysis_id <- paste(normalize(id_normal, id_int_sep),
                             normalize(id_tumor, id_int_sep), sep = id_sep)
  } else {
    sif$analysis_id <- normalize(sif$analysis_id)
  }
  dup_ids <- (function (x) which(x %in% x[duplicated(x)]))(sif$analysis_id)
  if (length(dup_ids)) {
    if (!have_id) {
      id_source <- if (have_ids) do.call(paste, c(as.list(id_cols), 'columns')) else 'file names'
      cat(sprintf(paste("IMPORTANT:",
                  "Automatically generated IDs are obtained by substituting '%s'\n",
                  "to all runs of non-alphanumeric characters and removing '%s'\n",
                  "occurrencies at the end or at the begninning of the id.\n",
                  "The IDs are generated starting from %s.\n",
                  "-----------------------------------------\n\n"), id_int_sep, id_int_sep,
                  id_source))
    }
    stop(sprintf('the ids at the following lines are conflicting:\n%s',
                 per_line_list(paste_df_lines(sif[dup_ids, ]))))
  }
  write.table(sif[, cols_to_write], output_path, sep = '\t',
              quote = FALSE, row.names = FALSE, na = '')
EOF
) "${input_file}" "${id_separator}" "${id_internal_separator}" "${output_path}"
