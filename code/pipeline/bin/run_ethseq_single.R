#!/usr/bin/env Rscript

file_pref   <- '^--file='
script_path <- dirname(gsub(file_pref, '',
                            grep(file_pref, commandArgs(trailingOnly = FALSE),
                                 value = TRUE)))
source(file.path(script_path, '..', 'R', 'utils.R'))

docstring <- '
Usage:
  run_ethseq_single.R [options] MODEL_FILE VCF_INPUT OUTPUT_FOLDER

-h --help          Shows help informations
MODEL_FILE         The file containing the model to use for ethnicity annotation.
VCF_INPUT          The file containing the genotypes of the samples to analyze.
OUTPUT_FOLDER      The folder where the output will be created.

options:
  --num_threads N   Specifies the name for the sample to be used in the
                      output [default: 1]
  --log_file PATH   The path to the output file
'

run_ethseq <- function (vcf_dat, output_path, model_path, threads = 1) {
  res <- ethseq.Analysis(target.vcf = vcf_dat,
                         out.dir = output_path,
                         model.gds = model_path,
                         space = '3D',
                         composite.model.call.rate = 0.99,
                         run.genotype = FALSE,
                         cores = threads)
  for (ff in file.path(output_path, paste0(c('Aggregated', 'Target'), '.gds'))) {
    if (file.exists(ff)) {
      file.remove(ff)
    }
  }
  res
}

main_function <- function () {
  if (!is.null(cli_args$log_file)) {
    sink(cli_args$log_file)
  }
  greater_than_0    <- 'The value must be greater than 0'
  num_threads       <- parse_arg_pred(cli_args, 'num_threads', check_if_positive,
                                      greater_than_0, converter = conv_to_integer)

  res <- run_ethseq(cli_args$VCF_INPUT, cli_args$OUTPUT_FOLDER, cli_args$MODEL_FILE,
             num_threads)
  quit('no', !res)
}

(make_main(docstring, c('EthSEQ'), main_function))()
