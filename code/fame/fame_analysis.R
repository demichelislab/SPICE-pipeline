# Installs or imports the required libraries
librarian::shelf(
  stringi,
  parallel,
  pipeR,
  data.table,
  dplyr,
  tidyr,
  purrr,
  juliangehring/HighSpeedStats
)

# Moves to the root folder of the repo (all paths are relative to the repo)
while (!dir.exists('code')) {
  setwd('..')
  if (getwd() == '/') {
    stop('Cannot find the root folder of the project')
  }
}

# Combines two binary matrices with a specific operator (defaults to '|')
combine_mats <- function (mm1, mm2, FUN = `|`) {
  res <- FUN(mm1, mm2)
  storage.mode(res) <- 'integer'
  res
}

matrix_to_vector <- function (mm) {
  dim(mm) <- NULL
  mm
}

# Construct a binary aberration matrix based on the aberration passed as input.
#   Undefined samples are taken in consideration
build_aberration_matrix <- function (dat, aberration_to_extract, na_aberration = 'nd') {
  mm                       <- dat[, -1]
  res                      <- mm == aberration_to_extract
  res[mm == na_aberration] <- NA
  rownames(res)            <- dat[[1]]
  storage.mode(res)        <- 'integer'
  res
}

# Computes the contingency tables for all the genes in the matrices passed in
#   in input. Expected format of the matrices is genes on the rows and samples
#   on the columns
compute_counts <- function (mat1, mat2) {
  are_matrices_identical <- identical(mat1, mat2)
  mms <- list(
    m1 = mat1[[1]],
    m2 = mat2[[1]]
  )
  nns <- pipeline({
    setNames(mms, c('row', 'col'))
    map(rownames)
  })
  ll <- map(nns, length)
  def_elems <- pipeline({
    mms
    map(~ {
      xx               <- !is.na(.x)
      storage.mode(xx) <- 'integer'
      xx
    })
  })
  mms_na <- pipeline({
    mms
    map2(
      names(mms),
      function(mm, nn) pipeline({
        mm
        replace_na(0L)
        list(1L - .)
        setNames(paste0(nn, c('', '_n')))
        map(~ .x * def_elems[[nn]])
      })
    )
    flatten
  })
  res <- with(
    mms_na,
    data.table(
      g1   = rep(nns$row, ll$row),
      g2   = rep(nns$col, each = ll$col),
      n_11 = matrix_to_vector(tcrossprod(  m1,   m2)),
      n_10 = matrix_to_vector(tcrossprod(  m1, m2_n)),
      n_01 = matrix_to_vector(tcrossprod(m1_n,   m2)),
      n_00 = matrix_to_vector(tcrossprod(m1_n, m2_n))
    )
  )
  if (are_matrices_identical) {
    res <- pipeline({
      ll$row
      diag
      upper.tri
      matrix_to_vector
      (res[.])
    })
  }
  res
}

# Runs the fisher tests on the count data
run_tests <- function (dat, filter_freq = 0, fisher_fun = ultrafastfet) {
  pipeline({
    dat
    mutate(
      n_tot      = n_11 + n_10 + n_01 + n_00,
      n_1        = (n_11 + n_10),
      n_2        = (n_11 + n_01),
      f_1        = n_1 / n_tot,
      f_2        = n_2 / n_tot,
      f_min      = pmin(f_1, f_2),
      odds_ratio = (n_11 / n_01) / (n_10 / n_00),
      p_value    = 1L,
      fdr        = 1L
    )
    filter(f_min > filter_freq)
    (dd ~ {
      if (nrow(dd) > 0) {
        pipeline({
          dd
          mutate(
            p_value = fisher_fun(n_11, n_10, n_01, n_00),
            fdr     = p.adjust(p_value, method = 'BH')
          )
          arrange(fdr, p_value, f_min)
          select(
            g1,
            g2,
            p_value,
            fdr,
            odds_ratio,
            f_min,
            everything()
          )
        })
      } else dd
    })
    as.data.table
  })
}

num_cores <- detectCores()

# Reads all the files with genomic data from the specified folder
data_per_gene <- pipeline({
  list.files('data/genomic', full.names = TRUE)
  map_dfr(fread)
  group_by(sample_id)
  filter(any(as_cn_disc != 'nd'))
  ungroup
  as.data.table
})

# All the possible combinations of these genes will be tested
interest_genes <- pipeline({
  fread('data/resources/gene_sets/cancer_and_druggable_genes.tsv.gz')
  pull(gene)
})

# The genomic data of the interest genes
interest_genes_data <- pipeline({
  data_per_gene[hugo %in% interest_genes]
  nest(data = -dataset)
  arrange(stri_order(dataset))
  (setNames(.$data, .$dataset))
  map( ~ pipeline({
      .x
      arrange(hugo, sample_id)
      as.data.table
    })
  )
})

# This matrix contains the allele specific data of the loaded datasets
as_cn_matrix <- pipeline({
  interest_genes_data
  map(~ dcast(.x, hugo ~ sample_id, value.var = 'as_cn_disc', fill = NA_character_))
})

# This builds the binary aberration matrix for a specific allele specific copy
#   number aberration. This matrix will have a value of 1 if a sample have
#   an Homozygous deletion at a specific gene 0 otherwise.
homo_del_matrix <- pipeline({
  as_cn_matrix
  map(~ build_aberration_matrix(.x, 'homo_del'))
})

# This matrix will have a value of 1 if a sample have an Hemizygous deletion at
#   a specific gene 0 otherwise.
hemi_del_matrix <- pipeline({
  as_cn_matrix
  map(~ build_aberration_matrix(.x, 'hemi_del'))
})

# This matrix will have a value of 1 if a sample have an Copy Neutral LOH at a
# specific gene 0 otherwise.
cnnl_matrix <- pipeline({
  as_cn_matrix
  map(~ build_aberration_matrix(.x, 'cnnl'))
})

# This matrix will have a value of 1 if a sample have an amplification at a
# specific gene 0 otherwise.
amp_matrix <- pipeline({
  as_cn_matrix
  map(~ build_aberration_matrix(.x, 'amp'))
})

# This matrix will have a value of 1 if a sample have an SNV at a
# specific gene 0 otherwise.
snv_matrix <- pipeline({
  interest_genes_data
  map(~ {
    pipeline({
      dcast(.x, hugo ~ sample_id, value.var = 'count_snvs_deleterious', fill = 0L, fun.aggregate = function (x) as.integer(sign(x)))
    })
  })
  map(~ as.matrix(.x[, -1], rownames.value = .x[[1]]))
})

all_datasets <- pipeline({
  data_per_gene
  (dataset)
  unique
  stri_sort
  setNames(.)
})

# Named list with the simple aberrations to test
aberrations_bases <- list(
  'hemi_del' = hemi_del_matrix,
  'cnnl'     = cnnl_matrix,
  'snv'      = snv_matrix
)

# Named list with the combined aberrations to test. In this case the aberrations
#   are combined with an OR operation.
aberration_combinations <- pipeline({
  list(
    hemi_cnnl     = c('hemi_del', 'cnnl'),
    hemi_cnnl_snv = c('hemi_del', 'cnnl', 'snv')
  )
  map(~ {
    selected_abs <- aberrations_bases[.x]
    map(all_datasets, function (nn) {
      cc <- map(selected_abs, ~ .x[[nn]])
      reduce(cc, combine_mats)
    })
  })
})

aberrations <- c(
  aberrations_bases,
  aberration_combinations
)

# This combines all the matrices for all the dataset in an unique matrix in
#   order to test all the dataset in a pan-cancer fashion.
pancancer_aberrations <- pipeline({
  aberrations
  map(~ list(pancancer = reduce(.x, cbind)))
})

# Generates all the possible combinations of aberrations. FAME efficiency enable
#   the possibility to test many aberrations combinations.
aberration_comparisons <- pipeline({
  c(
    'hemi_del',
    'cnnl',
    'hemi_cnnl',
    'hemi_cnnl_snv',
    'snv'
  )
  (expand.grid(
    a2               = .,
    a1               = .,
    stringsAsFactors = FALSE
  ))
  filter(
    a1 <= a2
  )
  select(a1, a2)
  as.data.table
})

# Tests all the aberrations combinations on all the pairs of genes separately
#   for each dataset.
results_per_project <- pipeline({
  aberration_comparisons
  pmap(c)
  map_dfr(~ {
    cat(.x, '\n')
    pipeline({
      aberrations[.x]
      map(~ { .x[all_datasets] })
      pmap(list)
      map(~ compute_counts(.x[1], .x[2]))
      mcmapply(names(.), FUN = function (dd, nn) {
        print(nn)
        res <- run_tests(dd)
        gc()
        res
      }, mc.cores = num_cores, SIMPLIFY = FALSE)
      bind_rows(.id = 'dataset')
      mutate(
        a1 = .x[[1]],
        a2 = .x[[2]]
      )
      select(a1, a2, everything())
      as.data.table
    })
  })
  (? gc())
})

fwrite(results_per_project, 'data/result_pairs.tsv', sep = '\t')

# Tests all the aberrations combinations on all the pairs of genes on the
#   pancancer dataset (all the datasets combined).
results_pancancer <- pipeline({
  aberration_comparisons
  pmap(c)
  map_dfr(~ {
    cat(.x, '\n')
    pipeline({
      pancancer_aberrations[.x]
      pmap(list)
      map(~ compute_counts(.x[1], .x[2]))
      mcmapply(names(.), FUN = function (dd, nn) {
        print(nn)
        res <- run_tests(dd)
        gc()
        res
      }, mc.cores = num_cores, SIMPLIFY = FALSE)
      bind_rows(.id = 'dataset')
      mutate(
        a1 = .x[[1]],
        a2 = .x[[2]]
      )
      select(a1, a2, everything())
      # filter(odds_ratio < 1)
      as.data.table
      (? gc())
    })
  })
  (? gc())
})

fwrite(results_pancancer, 'data/result_pairs_pancancer.tsv', sep = '\t')
