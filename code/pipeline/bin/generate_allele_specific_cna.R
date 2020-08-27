#!/usr/bin/env Rscript

file_pref   <- '^--file='
script_path <- dirname(gsub(file_pref, '',
                            grep(file_pref, commandArgs(trailingOnly = FALSE),
                                 value = TRUE)))
source(file.path(script_path, '..', 'R', 'utils.R'))

docstring <- '
Usage:
  generate_cna_ase.R [options] GENE_MODEL AI_TABLE [OUTPUT_PATH]

-h --help    Shows help informations
GENE_MODEL   BED file containing the ids of the genes and the respective positions
AI_TABLE     Input file containing the allelic imbalance table
OUTPUT_PATH  Where to put the outputs of the program (raw per-gene SCNAs and
                                                      discretized per-gene SCNAs)

options:
  --homo_log2_threshold HOMO_LOG2_THRESHOLD    Log2 corrected threshold to add (replace)
                                                 the status of genes with small homo del
  --max_cn_dist MAX_CN_DIST                    Maximum distance of the segment from the nearest
                                                 integer copy number allowed (to filter too noisy
                                                 segments)
  --gene_id_colname GENE_ID_COLNAME            The name of the column containing the gene id
                                                 in the gene model [default: HUGO]
  --segments_log2_colname LOG2_COLNAME         The name of the column containing the value to
                                                 use as segment copy number alteration [default: log2.corr]
  --segements_cna_colname CNA_COLNAME          The name of the column containing the allele specific copy number
                                                 A [default: cnA]
  --segements_cnb_colname CNB_COLNAME          The name of the column containing the allele specific copy number
                                                 B [default: cnB]
  --cn_disc_colname CN_DISC_COLNAME            The name of the column that will contain the discretized allele
                                                 specific number in the updated allelic imbalance
                                                 table [default: as_cn_disc]
  --segment_id_colname SEGMENT_ID_COLNAME      The name of the column that will contain the id of the segment
                                                 in the updated allelic imbalance table [default: segment_id]
  --updated_ai_table_filename FILENAME         The path for the updated allelic imbalance table. The file is
                                                 written only when this parameter is set.
  --genes_segment_id_matrix_filename FILENAME  The file containing the matrix that identifies which segment is
                                                 assigned to each gene. The file is written only when this parameter is set.
  --log_file PATH                              The path to the output file
  --gene_as_cn_file_name FILENAME              The name of the file containing raw SCNAs [default: genes_as_cn.tsv]
'

# get_cn_ab_matrix <- function (){
#   allele_values  <- get_allele_specific_values()
#   acast(allele_values, cn_a ~ cn_b, value.var = 'name_as_factor')
# }
#
# get_allele_specific_values <- function () {
#   res <- data.frame(
#            name             = c('N.D.', 'Homo Del', 'Hemi Del', 'CNNL', 'Del | Gain', 'WT', 'Gain', 'Amp'),
#            priority         = c(     8,          1,          3,      4,            5,    7,      6,     2),
#            cn_a             = c(     0,          0,          1,      2,            3,    1,      2,     3),
#            cn_b             = c(     1,          0,          0,      0,            0,    1,      1,     1),
#            stringsAsFactors = FALSE
#          )
#   within(res, {
#     name_as_factor <- factor(name, name[order(-priority)])
#   })
# }

get_nd_string <- function () {
  'nd'
}

get_allele_specific_order <- function () {
  c(
    get_nd_string(),
    'wt',
    'gain',
    'gain_unb',
    'gain_del',
    'cnnl',
    'hemi_del',
    'amp',
    'amp_unb',
    'amp_del',
    'homo_del'
  )
}

get_as_matrix <- function() {
  c(
      'homo_del',
      'hemi_del',       'wt',
          'cnnl',     'gain',    'gain',
      'gain_del', 'gain_unb',     'amp',     'amp',
      'gain_del',  'amp_unb', 'amp_unb',     'amp', 'amp',
       'amp_del',  'amp_unb', 'amp_unb', 'amp_unb', 'amp', 'amp'
  )
}

get_cn_ab_matrix <- function (){
  matrix_size               <- 6
  res                       <- factor(rep(get_nd_string(), matrix_size^2), get_allele_specific_order())
  dim(res)                  <- rep(matrix_size, 2)
  res[upper.tri(res, TRUE)] <- get_as_matrix()
  t(res)
}

cn_ab_to_named_value <- function (cn_a, cn_b) {
  cn_a_offset <- ifelse(cn_b > 3, cn_b - 3, 0)
  cn_a_pos    <- pmin(cn_a - cn_a_offset, 5) + 1L
  cn_b_pos    <- pmin(cn_b, 3) + 1L
  cn_matrix   <- get_cn_ab_matrix()
  cn_matrix[cbind(cn_a_pos, cn_b_pos)]
}

preprocess_cna <- function (cn_a,
                            cn_b,
                            log2_val,
                            homo_log2_threshold = NULL,
                            max_cn_dist         = NULL) {
  cn_a_res       <- cn_a
  cn_b_res       <- cn_b
  cn_a_round     <- round(cn_a)
  cn_b_round     <- round(cn_b)
  is_too_distant <- rep(FALSE, length(cn_a))
  if (!is.null(homo_log2_threshold)) {
    is_homo_deletion       <- !is.na(log2_val) & log2_val < homo_log2_threshold
    if (any(is_homo_deletion)) {
      cn_a_res[is_homo_deletion] <- 0
      cn_b_res[is_homo_deletion] <- 0
    }
  }
  is_na          <- is.na(cn_a_res) | is.na(cn_b_res)
  if (!is.null(max_cn_dist)) {
    is_too_distant       <- !is_na & (abs(cn_a_round - cn_a) > max_cn_dist | abs(cn_b_round - cn_b) > max_cn_dist)
  }
  is_negative            <- !is_na & (cn_a_round < 0 | cn_b_round < 0)
  to_set_to_na           <- is_na | is_too_distant | is_negative
  if (any(to_set_to_na)) {
    cn_a_res[to_set_to_na] <- 0
    cn_b_res[to_set_to_na] <- 1
  }
  list(
    a = cn_a_res,
    b = cn_b_res
  )
}

get_ai_table_with_discrete_allele_specific_levels <- function (ai_table,
                                                               segments_scna_colname = 'log2.corr',
                                                               gene_id_colname       = 'HUGO',
                                                               cn_a_colname          = 'cnA',
                                                               cn_b_colname          = 'cnB',
                                                               cn_disc_colname       = 'as_cn_disc',
                                                               segment_id_colname    = 'segment_id',
                                                               homo_log2_threshold   = NULL,
                                                               max_cn_dist           = NULL) {
  ai_table                    <- ai_table[stri_order(ai_table$sample), ]
  cn                          <- preprocess_cna(getElement(ai_table, cn_a_colname),
                                                getElement(ai_table, cn_b_colname),
                                                getElement(ai_table, segments_scna_colname),
                                                homo_log2_threshold,
                                                max_cn_dist)
  ai_table[[cn_disc_colname]]    <- cn_ab_to_named_value(round(cn$a), round(cn$b))
  ai_table[[segment_id_colname]] <- with(ai_table, paste(sample, chr, start, end, sep = ':'))
  ai_table
}

generate_gene_based_allele_specific_cna_matrix <- function (ai_table_with_cn_disc,
                                                            gene_model,
                                                            cn_disc_colname,
                                                            gene_id_colname,
                                                            segment_id_colname,
                                                            generate_segment_id_matrix = FALSE) {
  gene_model      <- gene_model[stri_order(getElement(gene_model, gene_id_colname)), ]
  sample_names    <- unique(ai_table_with_cn_disc$sample)
  segs_ranges     <- with(ai_table_with_cn_disc,   GRanges(chr, IRanges(start, end)))
  gene_ranges     <- with(gene_model, GRanges(chr, IRanges(start, end)))
  hugo_id         <- getElement(gene_model, gene_id_colname)
  range_matches   <- findOverlaps(segs_ranges, gene_ranges)
  genes_as_cns    <- matrix(get_nd_string(), length(hugo_id), length(sample_names), dimnames = list(hugo_id, sample_names))
  genes_seg_ids   <- if (generate_segment_id_matrix) {
    matrix(NA_character_, length(hugo_id), length(sample_names), dimnames = list(hugo_id, sample_names))
  }
  segs_hits       <- queryHits(range_matches)
  gene_hits       <- subjectHits(range_matches)
  all_genes_as_cn <- local({
    rr <- data.frame(
            sample_id        = match(ai_table_with_cn_disc$sample[segs_hits], sample_names),
            gene_id          = gene_hits,
            cn_value         = getElement(ai_table_with_cn_disc, cn_disc_colname)[segs_hits],
            segment_id       = getElement(ai_table_with_cn_disc, segment_id_colname)[segs_hits],
            stringsAsFactors = FALSE
          )
    rr[order(rr$cn_value), ]
  })
  idxes                <- as.matrix(all_genes_as_cn[, c('gene_id', 'sample_id')])
  genes_as_cns[idxes]  <- as.character(all_genes_as_cn$cn_value)
  if (generate_segment_id_matrix) {
    genes_seg_ids[idxes] <- all_genes_as_cn$segment_id
  }
  list(
    allele_specific_cns = genes_as_cns,
    segments_ids        = genes_seg_ids
  )
}

# generate_gene_based_allele_specific_cna_matrix <- function (ai_table,
#                                                             gene_model,
#                                                             segments_scna_colname = 'log2.corr',
#                                                             gene_id_colname       = 'HUGO',
#                                                             cn_a_colname          = 'cnA',
#                                                             cn_b_colname          = 'cnB',
#                                                             homo_log2_threshold   = NULL,
#                                                             max_cn_dist           = NULL) {
#   ai_table         <- ai_table[order(ai_table$sample), ]
#   gene_model       <- gene_model[order(getElement(gene_model, gene_id_colname)), ]
#   gm_ord           <- sort(getElement(gene_model, gene_id_colname))
#   sample_names     <- unique(ai_table$sample)
#   cn               <- preprocess_cna(getElement(ai_table, cn_a_colname),
#                                      getElement(ai_table, cn_b_colname),
#                                      getElement(ai_table, segments_scna_colname),
#                                      homo_log2_threshold,
#                                      max_cn_dist)
#   cn_named         <- cn_ab_to_named_value(round(cn$a), round(cn$b))
#   segs_ranges      <- with(ai_table,   GRanges(chr, IRanges(start, end)))
#   gene_ranges      <- with(gene_model, GRanges(chr, IRanges(start, end)))
#   hugo_id          <- getElement(gene_model, gene_id_colname)
#   range_matches    <- findOverlaps(segs_ranges, gene_ranges)
#   res              <- matrix(get_nd_string(), length(hugo_id), length(sample_names), dimnames = list(hugo_id, sample_names))
#   segs_hits        <- queryHits(range_matches)
#   gene_hits        <- subjectHits(range_matches)
#   all_genes_as_cn  <- local({
#     rr <- data.frame(
#             sample_id        = match(ai_table$sample[segs_hits], sample_names),
#             gene_id          = gene_hits,
#             cn_value         = cn_named[segs_hits],
#             stringsAsFactors = FALSE
#           )
#     rr[order(rr$cn_value), ]
#   })
#   res[as.matrix(all_genes_as_cn[, c('gene_id', 'sample_id')])] <- as.character(all_genes_as_cn$cn_value)
#   res
# }

# testing_code <- function () {
#   library(data.table)
#   library(GenomicRanges)
#   gm       <- fread('/shares/CIBIO-Storage/CO/SPICE/resources/gene_models/hg38/biomart.ensembl.hg38.20171221.bed', data.table = FALSE)
#   ai       <- fread('/shares/CIBIO-Storage/CO/SPICE/analysis/pipeline/data/aggregation/blca_all_manual_corrected-tcga_180116.tsv', data.table = FALSE)
#   # ai       <- fread('pixz -d < /shares/CIBIO-Storage/CO/SPICE/analysis/pipeline/data/aggregation/clonet_all_segs_finished_ok_first_run.tsv.xz', data.table = FALSE)
#   test_ddd <- fread('/shares/CIBIO-Storage/CO/SPICE/analysis/discretize_clonet_tables/output/discrete_cn_blca_tcga_ai_tb_new_pipeline.txt')
#   test_ddd <- as.matrix(test_ddd[, -1], rownames = test_ddd[[1]])
#
#   start_time <- Sys.time()
#   rrr      <- generate_gene_based_allele_specific_cna_matrix(ai, gm, homo_log2_threshold = -1.3)
#   end_time <- Sys.time()
#
#   end_time - start_time
#
#   dt_rrr   <- as.data.table(rrr)
#   dt_rrr[, N := .N, by = c('sample', 'gene')]
#   dt_rrr[N > 1]
# }

null_or_numeric <- function (x) if (!is.null(x)) suppressWarnings(as.numeric(x))

main_function <- function () {
    if (!is.null(cli_args$log_file)) {
      sink(cli_args$log_file)
    }
    homo_log2_threshold    <- parse_arg(cli_args, 'homo_log2_threshold', null_or_numeric)
    max_cn_dist            <- parse_arg(cli_args, 'max_cn_dist', null_or_numeric)

    gene_model             <- read.table(cli_args$GENE_MODEL, header = TRUE, sep = '\t', quote = '', stringsAsFactors = F, comment.char = '')
    ai_table               <- local({
      dd <- read.table(cli_args$AI_TABLE,
                       header           = TRUE,
                       sep              = '\t',
                       quote            = '',
                       stringsAsFactors = FALSE,
                       comment.char     = '')
      setNames(dd, c('sample', 'chr', 'start', 'end', tail(colnames(dd), -4)))
    })
    ai_table_colnames      <- colnames(ai_table)
    gene_id_colname        <- parse_arg(cli_args, 'gene_id_colname', function (arg) {
                                match.arg(arg, colnames(gene_model))
                              })
    segments_scna_col_name <- parse_arg(cli_args, 'segments_log2_colname', function (arg) {
                                match.arg(arg, ai_table_colnames)
                              })
    segements_cna_colname  <- parse_arg(cli_args, 'segements_cna_colname', function (arg) {
                                match.arg(arg, ai_table_colnames)
                              })
    segements_cnb_colname  <- parse_arg(cli_args, 'segements_cnb_colname', function (arg) {
                                match.arg(arg, ai_table_colnames)
                              })
    cn_disc_colname        <- parse_arg(cli_args, 'cn_disc_colname', function (arg) {
                                if (any(grepl('[[:space:]]', arg))) {
                                  die('Spaces cannot be part of the cn_disc_colname column name')
                                }
                                arg
                              })
    segment_id_colname     <- parse_arg(cli_args, 'segment_id_colname', function (arg) {
                                if (any(grepl('[[:space:]]', arg))) {
                                  die('Spaces cannot be part of the segment_id_colname column name')
                                }
                                arg
                              })

    generate_segment_id_matrix <- !is.null(cli_args$genes_segment_id_matrix_filename)

    output_path            <- cli_args$OUTPUT_PATH

    ai_table_with_cn_disc  <- get_ai_table_with_discrete_allele_specific_levels(ai_table,
                                                                                segments_scna_col_name,
                                                                                gene_id_colname,
                                                                                segements_cna_colname,
                                                                                segements_cnb_colname,
                                                                                cn_disc_colname,
                                                                                segment_id_colname,
                                                                                homo_log2_threshold,
                                                                                max_cn_dist)

    result_matrices        <- generate_gene_based_allele_specific_cna_matrix(
                                ai_table_with_cn_disc,
                                gene_model,
                                cn_disc_colname,
                                gene_id_colname,
                                segment_id_colname,
                                generate_segment_id_matrix
                              )

    if (is.null(output_path)) {
      output_path <- getwd()
    }
    if (!file.exists(output_path)) {
      dir.create(output_path, recursive = TRUE)
    }
    write.table(data.frame(hugo = rownames(result_matrices$allele_specific_cns),
                           result_matrices$allele_specific_cns,
                           stringsAsFactors = FALSE),
                file.path(output_path, cli_args$gene_as_cn_file_name),
                row.names = FALSE, sep = '\t', quote = FALSE, na = '')
    if (generate_segment_id_matrix) {
      write.table(data.frame(hugo = rownames(result_matrices$segments_ids),
                             result_matrices$segments_ids,
                             stringsAsFactors = FALSE),
                  file.path(output_path, cli_args$genes_segment_id_matrix_filename),
                  row.names = FALSE, sep = '\t', quote = FALSE, na = '')
    }
    if (!is.null(cli_args$updated_ai_table_filename)){
      write.table(ai_table_with_cn_disc,
                  file.path(output_path, cli_args$updated_ai_table_filename),
                  row.names = FALSE, sep = '\t', quote = FALSE, na = '')
    }


}

(make_main(docstring, c('stringi', 'IRanges', 'GenomicRanges'), main_function))()
