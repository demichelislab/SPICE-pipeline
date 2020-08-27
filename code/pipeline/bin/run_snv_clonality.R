#!/usr/bin/env Rscript

file_pref   <- '^--file='
script_path <- dirname(gsub(file_pref, '',
                            grep(file_pref, commandArgs(trailingOnly = FALSE),
                                 value = TRUE)))
source(file.path(script_path, '..', 'R', 'utils.R'))

docstring <- '
Usage:
  run_snv_clonality.R [options] CLONET_FOLDER NORMAL_SAMPLE_NAME TUMOR_SAMPLE_NAME ALLELIC_IMBALANCE_TABLE VEP_TABLE_OUTPUT OUTPUT_FOLDER

-h --help                Shows help informations
CLONET_FOLDER            The folder containing clonet functions
NORMAL_SAMPLE_NAME       The name of the normal sample
TUMOR_SAMPLE_NAME        The name of the tumor sample
ALLELIC_IMBALANCE_TABLE  The path of the file containing the allelic imbalance table for the sample.
VEP_TABLE_OUTPUT         The path of the file containing VEP output in the tabular format
OUTPUT_FOLDER            The folder where the output will be created.

options:
  --log_file PATH  The path to the output file
'


my_read_table <- function (...)  read.table(..., stringsAsFactors = FALSE,
                                            header = TRUE, sep = '\t',
                                            comment.char = '', quote = '')

read_vep_table <- function (filename) {
  vep_lines       <- readLines(filename)
  is_heading_line <- grepl('^#', vep_lines)
  dl              <- vep_lines[!is_heading_line]
  lr              <- c(sub('^#', '', tail(vep_lines[is_heading_line], 1)), dl)
  list(vep_data      = my_read_table(text = lr, na.strings = c('-', '.')),
       raw_lines     = vep_lines,
       heading_lines = which(is_heading_line))
}

get_snv_containing_lines <- function (vep_data) {
  vep_cmp_fields    <- c('Location', 'Allele')
  vep_keys          <- apply(vep_data[, vep_cmp_fields], 1, function (dd) paste(dd, collapse = ':'))
  vep_only_snv_keys <- vep_keys[vep_data$VARIANT_CLASS == 'SNV']
  match(vep_keys, vep_only_snv_keys)
}

load_clonet_functions <- function (clonet_path) {
  files_to_load  <- paste('CLONET', c('BasicFunctions', 'clonality_snv_functions'),
                          'R', sep = '.')
  functions_dir  <- file.path(clonet_path, 'Functions')
  files_paths <- file.path(functions_dir, files_to_load)
  if (!dir.exists(functions_dir)) {
      die(paste('Folder', functions_dir, 'does not exist or is not a folder'))
  }
  invisible(lapply(files_paths, source))
}

snv_clonality_cols <- c('t_cov', 't_af', 'n_admReads', 't_ref_count_corr',
                        't_af_corr', 'cn.int', 'CN_SNVmut', 'VAFexp',
                        'SNV.clonality', 'SNV.clonality.int')

main_function <- function () {
  if (!is.null(cli_args$log_file)) {
    sink(cli_args$log_file)
  }
  load_clonet_functions(cli_args$CLONET_FOLDER)
  output_name             <- sub('(\\.[^.]+)?$', '', basename(cli_args$VEP_TABLE_OUTPUT))
  output_path             <- file.path(cli_args$OUTPUT_FOLDER, paste0(output_name, '.tsv'))
  allelic_imbalance_table <- my_read_table(cli_args$ALLELIC_IMBALANCE_TABLE)
  vep_data                <- read_vep_table(cli_args$VEP_TABLE_OUTPUT)
  preprocessed_vep_data   <- preprocess_snvtable_vep(vep_data$vep_data, cli_args$NORMAL_SAMPLE_NAME,
                                                     cli_args$TUMOR_SAMPLE_NAME)
  snv_extended            <- extend_SNVt_with_bt(preprocessed_vep_data,
                                                 allelic_imbalance_table,
                                                 Ncores = 1)
  snvs_with_clonality     <- do.call(rbind, lapply(1:nrow(snv_extended),
                                                   computeSNVclonality,
                                                   snv_extended))
  vep_heading_lines       <- with(vep_data, raw_lines[heading_lines])
  writeLines(c(head(vep_heading_lines, -1),
               paste0(c(tail(vep_heading_lines, 1), snv_clonality_cols), collapse = '\t')),
             output_path)
  output_data <- data.frame(with(vep_data, raw_lines[setdiff(1:length(raw_lines), heading_lines)]),
                            snvs_with_clonality[get_snv_containing_lines(vep_data$vep_data),
                                                snv_clonality_cols])
  write.table(output_data, output_path, append = TRUE, sep = '\t',
              quote = FALSE, col.names = FALSE, row.names = FALSE, na = '-')
}

(make_main(docstring, c('parallel'), main_function))()
