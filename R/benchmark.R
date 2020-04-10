#' Test meta-analysis against contrasts for a drug in CMAP.
#' @export
benchmark_ma <- function(dprimes, anals, drug_info = NULL, suffix = "", data_dir = getwd(), pvals = NULL, pval_breaks = c(seq(0.01, 0.3, 0.01), seq(0.3, 1, 0.05))) {


  if (is.null(drug_info)) {
    data('cmap_es', package = 'ccdata')
    drug_info <- cmap_es
    rm(cmap_es)
  }


  # get contrast results
  cons_res <- lapply(dprimes$contrasts, function(con) {
    crossmeta::query_drugs(con, drug_info)
  })

  # get contrast results using pvals to choose ngenes
  cons_pres <- list()
  cons_ngenes <- list()

  for (p in pval_breaks) {
    pc <- as.character(p)
    ngenes <- sapply(pvals$contrasts, function(con) sum(con <= p))

    pres <- lapply(seq_along(dprimes$contrasts), function(i) {
      con <- dprimes$contrasts[[i]]
      crossmeta::query_drugs(con, drug_info, ngenes = ngenes[i])
    })
    names(pres) <- names(dprimes$contrasts)
    cons_pres[[pc]] <- pres
    cons_ngenes[[pc]] <- ngenes
  }

  # get sample sizes
  nsamples  <- c()
  contrasts <- names(cons_res)

  for (con in contrasts) {
    # gse_name and groups
    gse_name <- gsub('_.+?$', '', con)
    groups <- gsub('^.+?_', '', con)
    groups <- strsplit(groups, '-')

    test <- groups[[1]][1]
    ctrl <- groups[[1]][2]

    # sample size for contrasts
    design  <- anals[[gse_name]]$ebayes_sv$design
    nsamples <- c(nsamples, sum(design[, c(test, ctrl)]))
  }

  names(nsamples) <- contrasts

  # get meta-analysis results
  meta_res <- crossmeta::query_drugs(dprimes$meta, drug_info)

  # get meta-analysis results using pvals to choose ngenes
  meta_ngenes <- c()
  meta_pres <- list()
  for (p in pval_breaks) {
    pc <- as.character(p)
    ngenes <- sum(pvals$meta < p)

    meta_pres[[pc]] <- crossmeta::query_drugs(dprimes$meta, drug_info, ngenes = ngenes)
    meta_ngenes[pc] <- ngenes
  }

  res <- list(meta = meta_res, cons = cons_res, nsamples = nsamples,
              meta_pres = meta_pres, meta_ngenes = meta_ngenes,
              cons_pres = cons_pres, cons_ngenes = cons_ngenes)

  saveRDS(res, file.path(data_dir, paste0("ma_res", suffix, ".rds")))

  return(res)
}

#' Compare similarity metrics
#' @export
benchmark_ccmap <- function(dprimes, drug_info = NULL, drug_prl = NULL, suffix = "") {

  if (is.null(drug_info)) {
    library(ccdata)
    data(cmap_es)
    drug_info <- cmap_es
    rm(cmap_es)
  }

  if (is.null(drug_prl)) {
    drug_names <- colnames(drug_info)
    drug_prl <- drug_info

    for (drug in drug_names) {
      rank <- order(drug_prl[, drug])
      drug_prl[rank, drug] <- nrow(drug_prl):1
    }
  }

  # get meta ccmap results
  ccmap_meta <- crossmeta::query_drugs(dprimes$meta, drug_info)

  # meta xcos results
  xcos_meta <- xcos(dprimes$meta, drug_info, 100)

  # meta xsum results
  xsum_meta <- xsum(dprimes$meta, drug_info, 100)

  # meta ks results
  ks_meta <- ks(dprimes$meta, drug_prl, 100)

  # meta euclidean results
  euc_meta <- euclidean(dprimes$meta, drug_info, 100)

  # pearson correlation
  pear_meta <- pear(dprimes$meta, drug_info, 100)

  #spearman correlation
  spear_meta <- spear(dprimes$meta, drug_info, 100)


  # get contrast ccmap results
  ccmap_cons <- lapply(dprimes$contrasts, function(con) {
    crossmeta::query_drugs(con, drug_info)
  })

  # get contrast xcos results
  xcos_cons <- lapply(dprimes$contrasts, function(con) {
    xcos(con, drug_info, 100)
  })

  # get contrast xsum results
  xsum_cons <- lapply(dprimes$contrasts, function(con) {
    xsum(con, drug_info, 100)
  })

  # get contrast ks results
  ks_cons <- lapply(dprimes$contrasts, function(con) {
    ks(con, drug_prl, 100)
  })

  # get contrast euclidean results
  euc_cons <- lapply(dprimes$contrasts, function(con) {
    euclidean(con, drug_info, 100)
  })

  # get contrast pearson results
  pear_cons <- lapply(dprimes$contrasts, function(con) {
    pear(con, drug_info, 100)
  })

  # get contrast euclidean results
  spear_cons <- lapply(dprimes$contrasts, function(con) {
    spear(con, drug_info, 100)
  })


  res <- list(ccmap_meta = ccmap_meta, xcos_meta = xcos_meta,
              xsum_meta = xsum_meta, ks_meta = ks_meta, euc_meta = euc_meta,
              pear_meta = pear_meta, spear_meta = spear_meta,
              ccmap_cons = ccmap_cons, xcos_cons = xcos_cons,
              xsum_cons = xsum_cons, ks_cons = ks_cons, euc_cons = euc_cons,
              pear_cons = pear_cons, spear_cons = spear_cons
  )

  saveRDS(res, paste0("ccmap_res", suffix, ".rds"))

  return(res)
}

#' Compare different number of selected genes
#' @export
benchmark_ngenes <- function(dprimes, drug_info, data_dir = getwd(), suffix = '_ngenes') {

  if (is.null(drug_info)) {
    library(ccdata)
    data(cmap_es)
    drug_info <- cmap_es
    rm(cmap_es)
  }

  max_genes <- max(sapply(dprimes$contrasts, length))

  # get contrast results
  cons_res <- list()

  for (n in seq(50, max_genes, 50)) {
    nc <- as.character(n)

    res <- lapply(dprimes$contrasts, function(con) {
      crossmeta::query_drugs(con, drug_info, ngenes = n)
    })

    cons_res[[nc]] <- res

  }


  # get meta-analysis results
  meta_res <- list()

  for (n in seq(50, max_genes, 50)) {
    nc <- as.character(n)
    res <- crossmeta::query_drugs(dprimes$meta, drug_info, ngenes = n)
    meta_res[[nc]] <- res
  }

  saveRDS(res, file.path(data_dir, paste0("ma_res", suffix, ".rds")))

  return(res)
}


#' Test characteristic direction authors original processing of l1000 data
#' @export
benchmark_chdir <- function(dprimes, l1000_es = NULL, chdir_es = NULL) {

  # make common
  common   <- intersect(colnames(l1000_es), colnames(chdir_es))
  l1000_es <- l1000_es[, common]
  chdir_es <- chdir_es[, common]


  # get contrast ccmap (xcos) results
  ccmap_res <- lapply(dprimes$contrasts, function(con) crossmeta::query_drugs(con, l1000_es))

  # get contrast CD (cosine distance) results
  cd_res <- lapply(dprimes$contrasts, function(con) query_cd(con, chdir_es))

  # save results
  saveRDS(ccmap_res, "chdir_res_ccmap.rds")
  saveRDS(cd_res,    "chdir_res.rds")

}

#' Test my characteristic direction processing of l1000 data
#' @export
benchmark_cd <- function(dprimes, cd_es = NULL, l1000_es = NULL) {

  # make common
  common   <- intersect(colnames(cd_es), colnames(l1000_es))
  cd_es    <- cd_es[, common]
  l1000_es <- l1000_es[, common]


  # get contrast ccmap (cos dist) results
  cd_res <- lapply(dprimes$contrasts, function(con) crossmeta::query_drugs(con, cd_es))

  # get contrast CD results
  l1000_res <- lapply(dprimes$contrasts, function(con) query_cd(con, l1000_es))

  # save results
  saveRDS(l1000_res, "cd_res_ccmap.rds")
  saveRDS(cd_res,   "cd_res.rds")

}


#' Compare diff_cds to diff_expr processed query signatures
#' @export
benchmark_diff_cds <- function(anals_cd, dprimes, cd_es, l1000_es) {

  # make anals_cd and dprimes common
  dprimes$contrasts <- dprimes$contrasts[names(anals_cd$contrasts)]

  # get contrast CD (cos dist) results
  cd_res <- lapply(anals_cd$contrasts, function(con) query_cd(con, cd_es))

  # get contrast ccmap results
  l1000_res <- lapply(dprimes$contrasts, function(con) crossmeta::query_drugs(con, l1000_es))

  # save results
  saveRDS(l1000_res, "diff_cd_res_ccmap.rds")
  saveRDS(cd_res,   "diff_cd_res.rds")
}
