#!/usr/bin/env Rscript

file_pref   <- '^--file='
script_path <- dirname(gsub(file_pref, '',
                            grep(file_pref, commandArgs(trailingOnly = FALSE),
                                 value = TRUE)))
source(file.path(script_path, '..', 'R', 'utils.R'))

docstring <- '
Usage:
  generate_cna.R [options] GENE_MODEL SEGFILE [OUTPUT_PATH]

-h --help    Shows help informations
GENE_MODEL   BED file containing the ids of the genes and the respective positions
SEGFILE      Input file containing the segments
OUTPUT_PATH  Where to put the outputs of the program (raw per-gene SCNAs and
                                                      discretized per-gene SCNAs)

options:
  --scna_breaks BREAKS                 The array of values to use to discretize the
                                       SCNAs [default: -Inf,-1.4,-0.4,0.4,1.2,Inf]
  --parallel_samples N                 How many samples in parallel to analyze [default: 1]
  --parallel_genes N                   How many genes to analyze in parallel [default: 1]
  --gene_id_col_name NAME              The name of the column containing the gene id
                                       in the gene model [default: HUGO]
  --segments_scna_col_name NAME        The column containing the value to use as segment
                                       copy number alteration [default: log2.corr]
  --scna_summarization_method METHOD   The method to use to summarize multiple segments
                                       in one gene (max: max absolute CNA, mean, median,
                                       w.mean) [default: max]
  --wild_type_value N                  The value that will be subtracted from the quantized
                                       SCNAs so the wild type version of the gene is 0
                                       [default: 3]
  --raw_scnas_file_name NAME           The name of the file containing raw SCNAs
                                       [default: genes_scnas_raw.tsv]
  --discrete_scnas_file_name NAME      The name of the file containing discrete SCNAs
                                       [default: genes_scnas_discrete.tsv]
  --log_file PATH                      The path to the output file
  --transpose                          Write the output transposed (Genes as rows,
                                         samples as columns)
'

generate_gene_based_scna_matrix <- function (x, gene_model, method = 'max',
                                        segments_scna_col_name = 'log2.corr',
                                        gene_id_col_name = 'HUGO',
                                        parallel_samples = 1, parallel_genes = 1) {
  samples        <- unique(x$sample)
  n_samples      <- length(samples)
  gene_ids       <- getElement(gene_model, gene_id_col_name)
  res            <- matrix(rep(NA, nrow(gene_model) * n_samples), ncol = n_samples)
  colnames(res)  <- samples
  rownames(res)  <- gene_ids
  check_all_na   <- function (f) function (dd) if (any(complete.cases(dd))) f(dd) else NA
  fun            <- switch(method,
                            mean   = check_all_na(function (dd) mean(dd$v, na.rm = T)),
                            median = check_all_na(function (dd) median(dd$v, na.rm = T)),
                            max    = check_all_na(function (dd) dd[which.max(abs(dd$v)), , drop = FALSE]$v),
                            w.mean = check_all_na(function (dd) weighted.mean(dd$v, dd$s, na.rm = T))
                          )
  segments       <- with(x,          GRanges(chr, IRanges(start, end)))
  gene_intervals <- with(gene_model, GRanges(chr, IRanges(start, end)))

  matches        <- findOverlaps(segments, gene_intervals)
  samples_groups <- split(seq_len(nrow(x)), x$sample)

  res_apply <- mclapply(setNames(samples, samples), function (sn) {
    sample_matches   <- matches[queryHits(matches) %in% samples_groups[[sn]]]
    dd               <- data.frame(v = getElement(x[queryHits(sample_matches), , drop = FALSE],
                                                  segments_scna_col_name),
                                   s = width(segments[queryHits(sample_matches)]))
    segments_by_gene <- split(dd, subjectHits(sample_matches))
    present_genes    <- gene_ids[as.integer(names(segments_by_gene))]
    res_by_gene      <- mclapply(setNames(segments_by_gene, present_genes), fun,
                                 mc.cores = parallel_genes)
  }, mc.cores = parallel_samples, mc.preschedule = FALSE)
  for (ss in samples) {
    dd                 <- res_apply[[ss]]
    res[names(dd), ss] <- unlist(dd)
  }
  t(res[order(rownames(res)), order(colnames(res)), drop = FALSE])
}

quantize_scna_matrix <- function (x, breaks, wt_val = 3) {
  res <- x
  if (!all(is.na(x))) {
    res <- matrix(cut(x, breaks, labels = F) - wt_val, nrow = nrow(x),
          dimnames = list(rownames(x), colnames(x)))
  }
  res
}

parse_scna_breaks <- function (arg) {
  res <- as.numeric(unlist(strsplit(arg, ',')))
  if (any(is.na(res) | is.nan(res))) {
    stop('All arguments must be numeric')
  }
  if (!identical(sort(res), res)) {
    stop('Breaks must be in increasing order')
  }
  if (!identical(unique(res), res)) {
    stop('Breaks must be distinct')
  }
  if (length(res) < 2) {
    stop('At least two values are needed')
  }
  if (length(intersect(res, c(-Inf, Inf))) < 2) {
    print('WARNING: the breaks could not cover the entire possible domain of CNA values')
  }
  res
}

main_function <- function () {
    if (!is.null(cli_args$log_file)) {
      sink(cli_args$log_file)
    }
    scna_breaks            <- parse_arg(cli_args, 'scna_breaks', parse_scna_breaks)
    parallel_samples       <- parse_arg(cli_args, 'parallel_samples', conv_to_integer)
    parallel_genes         <- parse_arg(cli_args, 'parallel_genes', conv_to_integer)
    wild_type_value        <- parse_arg(cli_args, 'wild_type_value', conv_to_integer)
    method                 <- parse_arg(cli_args, 'scna_summarization_method', function (arg) {
                                match.arg(arg, c('max', 'mean', 'median', 'w.mean'))
                              })
    is_transposed          <- cli_args$transpose
    object_name            <- if (is_transposed) 'hugo' else 'sample'
    flip_if_transposed     <- if (is_transposed) t else identity

    gene_model             <- read.table(cli_args$GENE_MODEL, header = TRUE, sep = '\t', quote = '', stringsAsFactors = F, comment.char = '')
    segfile                <- local({
      dd <- read.table(cli_args$SEGFILE, header = TRUE, sep = '\t', quote = '', stringsAsFactors = F, comment.char = '')
      setNames(dd, c('sample', 'chr', 'start', 'end', tail(colnames(dd), -4)))
    })
    gene_id_col_name       <- parse_arg(cli_args, 'gene_id_col_name', function (arg) {
                                match.arg(arg, colnames(gene_model))
                              })
    segments_scna_col_name <- parse_arg(cli_args, 'segments_scna_col_name', function (arg) {
                                match.arg(arg, colnames(segfile))
                              })
    genes_scnas            <- flip_if_transposed(generate_gene_based_scna_matrix(segfile, gene_model, method = method,
                                                              segments_scna_col_name = segments_scna_col_name,
                                                              gene_id_col_name = gene_id_col_name,
                                                              parallel_samples = parallel_samples,
                                                              parallel_genes = parallel_genes))
    genes_scnas_quantized  <- quantize_scna_matrix(genes_scnas, scna_breaks,
                                                   wt_val = wild_type_value)

    output_path            <- cli_args$OUTPUT_PATH

    if (is.null(output_path)) {
      output_path <- getwd()
    }
    if (!file.exists(output_path)) {
      dir.create(output_path, recursive = TRUE)
    }
    write.table(setNames(data.frame(rownames(genes_scnas), genes_scnas, stringsAsFactors = FALSE),
                         c(object_name, colnames(genes_scnas))),
                file.path(output_path, cli_args$raw_scnas_file_name),
                row.names = FALSE, sep = '\t', quote = FALSE, na = '')
    write.table(setNames(data.frame(rownames(genes_scnas_quantized), genes_scnas_quantized, stringsAsFactors = FALSE),
                         c(object_name, colnames(genes_scnas_quantized))),
                      file.path(output_path, cli_args$discrete_scnas_file_name),
                      row.names = FALSE, sep = '\t', quote = FALSE, na = '')
}

(make_main(docstring, c('parallel', 'IRanges', 'GenomicRanges'), main_function))()
