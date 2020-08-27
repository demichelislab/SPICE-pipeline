#!/usr/bin/env bash

usage () {
cat <<USAGE
USAGE: $0
  output_path
  log_path
  oncotator_output_file
  gene_model_file
  output_file_name
USAGE
}

required_args=5
. $(dirname $0)/parse_args_bash.sh

analysis_name='generate_mutation_matrix'

. $(dirname $0)/get_common_paths.sh

shopt -s extglob
output_path="${1%%+(/)}"
log_path="${2%%+(/)}"
output_prefix="$output_path"
log_prefix="$log_path"
oncotator_output="$3"
gene_model="$4"
output_file_name="$5"
output_file_path="${output_prefix}/$output_file_name"

exec &> "$log_path"

Rscript <(
cat <<EOF
  prepare_changes                <- function (cc) {
    res <- gsub('^p\\\.', '', unique(cc))
    res <- as.list(setdiff(res, if (length(res) == 1) '' else 'wt'))
    do.call(paste, c(res, sep = ','))
  }
  oncotator                        <- read.table('$oncotator_output', skip = 3, stringsAsFactors = F, sep = '\t', quote = '', header = T, comment.char = '')
  gene_model                       <- read.table('$gene_model', stringsAsFactors = F, sep = '\t', quote = '', header = T, comment.char = '')
  common_genes                     <- oncotator\$Hugo_Symbol %in% gene_model\$HUGO
  dat                              <- setNames(oncotator[common_genes, c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification')], c('sample', 'gene', 'change'))
  dat_unique_idx                   <- unique(dat[, -3])
  dat[dat\$change == '', 'change'] <- 'wt'
  samples_names                    <- unique(dat\$sample)
  res                              <- matrix(rep('wt', length(samples_names) * nrow(gene_model)), nrow = length(unique(samples_names)),
                                             dimnames = list(samples_names, gene_model\$HUGO))
  muts_combined                    <- by(dat\$change, paste(dat\$sample, dat\$gene, sep = ':'), prepare_changes)
  invisible(lapply(seq_len(nrow(dat_unique_idx)), function (nn) {
    row                         <- dat_unique_idx[nn, ]
    res[row\$sample, row\$gene] <<- muts_combined[paste(row\$sample, row\$gene, sep = ':')]
  }))
  res <- res[order(rownames(res)), order(colnames(res))]
  write.table(cbind(sample = sort(samples_names), res), '$output_file_path', sep = '\t', quote = F, row.names = F, na = '')
EOF
)
