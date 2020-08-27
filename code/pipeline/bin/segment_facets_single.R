#!/usr/bin/env Rscript

file_pref   <- '^--file='
script_path <- dirname(gsub(file_pref, '',
                            grep(file_pref, commandArgs(trailingOnly = FALSE),
                                 value = TRUE)))
source(file.path(script_path, '..', 'R', 'utils.R'))

docstring <- '
Usage:
  segment_facets_single.R [options] FACETS_INPUT [OUTPUT_SEGFILE]

-h --help           Shows help informations
FACETS_INPUT        The file containing the input for FACETS segmentation.
OUTPUT_SEGFILE      Output file containing the segments.

options:
  --reference_genome GENOME_NAME  Specifies the reference genome to use [default: hg19]
  --sample_name NAME              Specifies the name for the sample to be used in the
                                    output [default: AUTO]
  --cval CVAL                     Changes the critical value for segmentation [default: 150]
  --err_thresh THRESH             Ignores positions that exceed the error threshold [default: Inf]
  --del_thresh THRESH             Ignores positions that exceed the deletion threshold [default: Inf]
  --min_depth DEPTH               Specifies the minimum depth in normal for positions to be
                                    considered [default: 35]
  --max_depth DEPTH               Specifies the maximum depth in normal for position to be
                                    considered [default: 1000]
  --old_pileup                    If set the input file is read as the old format.
  --log_file PATH                 The path to the output file
'

get_segments_from_results <- function (x, sample_name, have_chr = FALSE) {
  chrs_map <- c(as.character(1:22), 'X', 'Y')
  chr_pref <- if (have_chr) 'chr' else ''
  ok_segs  <- x$jointseg[!is.na(x$jointseg$cnlr), ]
  n_mark   <- x$out$num.mark
  mark_sum <- cumsum(n_mark)
  data.frame(
    id               = sample_name,
    chrom            = paste0(chr_pref, chrs_map[x$out$chr]),
    loc_start        = ok_segs$maploc[c(1, head(mark_sum, -1) + 1)],
    loc_end          = ok_segs$maploc[mark_sum],
    num_mark         = n_mark,
    seg_mean         = x$out$cnlr.median,
    stringsAsFactors = FALSE
  )
}

run_facets <- function (dat, cval = 150, min_depth = 35, max_depth = 1000, reference_genome = 'hg19') {
  preproc_dat <- preProcSample(dat, ndepth = min_depth, ndepthmax = max_depth,
                               gbuild = sub('_nochr$', '', reference_genome))
  procSample(preproc_dat, cval = cval)
}

is_chr_sample <- function (ff) {
  ll <- readLines(ff, 1000)
  any(grepl('^chr', ll[-1], ignore.case = TRUE))
}

main_function <- function () {
  if (!is.null(cli_args$log_file)) {
    sink(cli_args$log_file)
  }
  available_genomes <- c('hg18', 'hg19', 'hg19_nochr', 'hg38')
  greater_than_0    <- 'The value must be greater than 0'
  cval              <- parse_arg_pred(cli_args, 'cval', check_if_positive,
                                      greater_than_0, converter = conv_to_numeric)
  err_thresh        <- parse_arg_pred(cli_args, 'err_thresh', check_if_positive,
                                      greater_than_0, converter = conv_to_numeric)
  del_thresh        <- parse_arg_pred(cli_args, 'del_thresh', check_if_positive,
                                      greater_than_0, converter = conv_to_numeric)
  min_depth         <- parse_arg_pred(cli_args, 'min_depth', check_if_positive,
                                      greater_than_0, converter = conv_to_numeric)
  max_depth         <- parse_arg_pred(cli_args, 'max_depth', check_if_positive,
                                      greater_than_0, converter = conv_to_numeric)
  reference_genome  <- parse_arg_pred(cli_args, 'reference_genome',
                                      do.call(check_arg_in_list ,as.list(available_genomes)),
                                      paste('Valid genomes are:', paste0(available_genomes, collapse = ', ')))

  out_file          <- cli_args$OUTPUT_SEGFILE
  sample_name       <- cli_args$sample_name

  if (sample_name == 'AUTO') {
    sample_name <- gsub('\\..+$', '', basename(cli_args$FACETS_INPUT))
  }

  if (is.null(out_file)) {
    out_file    <- stdout()
  }

  have_chr      <- is_chr_sample(cli_args$FACETS_INPUT)
  dat           <- readSnpMatrix(cli_args$FACETS_INPUT, err.thresh = err_thresh,
                                 del.thresh = del_thresh, perl.pileup = cli_args$old_pileup)
  facets_output <- run_facets(dat, cval, min_depth, max_depth, reference_genome = reference_genome)
  segfile       <- get_segments_from_results(facets_output, sample_name, have_chr)
  write.table(segfile, out_file, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

(make_main(docstring, c('facets'), main_function))()
