
#' Run Correlation Query
#'
#' @param query_genes Numeric vector of effect size values with genes as names.
#' @param drug_info Matrix of effect size values to query against. Rows names are genes, column names are signature names.
#' @param ngenes Number of top upregulated and number of top downregulated \code{query_genes} that are in
#'   common with genes in \code{drug_info}. 
#' @inheritParams stats::cor
#'
#' @return Vector of query values with names of corresponding signature in \code{drug_info} sorted by decreasing similarity.
#' @export
#' 
query_cor <- function(query_genes, drug_info, ngenes = 100, method = 'pearson') {
  
  # use only common genes
  query_genes <- query_genes[names(query_genes) %in% row.names(drug_info)]
  
  # top 100 up/down genes
  query_genes <- sort(query_genes, TRUE)
  query_genes <- c(head(query_genes, ngenes), tail(query_genes, ngenes))
  drug_info   <- drug_info[names(query_genes), ,drop = FALSE]
  
  # correlation
  res <- cor(query_genes, drug_info, method = method)
  res <- structure(c(res), names=colnames(res))
  
  return(sort(res, decreasing = TRUE))
}


#' Run XCOS Query
#'
#' @inheritParams query_cor
#' @inherit query_cor return
#' @export
#'
query_xcos <- function(query_genes, drug_info, ngenes = 100) {
  
  # use only common genes
  query_genes <- query_genes[names(query_genes) %in% row.names(drug_info)]
  
  # top 100 up/down genes
  query_genes <- sort(query_genes, TRUE)
  query_genes <- c(head(query_genes, ngenes), tail(query_genes, ngenes))
  drug_info   <- drug_info[names(query_genes), ,drop = FALSE]
  
  # cosine similarity
  drug_info   <- sweep(drug_info, 1, query_genes, `*`)
  xcos        <- colSums(drug_info)
  names(xcos) <- colnames(drug_info)
  
  return(sort(xcos, decreasing = TRUE))
  
}

#' Run Xsum Query
#'
#' @inheritParams query_cor
#' @inherit query_cor return
#' @export
#'
query_xsum <- function(query_genes, drug_info, ngenes = 100) {
  
  # same genes
  drug_info <- drug_info[row.names(drug_info) %in% names(query_genes), ]
  query_genes <- query_genes[names(query_genes) %in% row.names(drug_info)]
  
  query_genes <- sort(query_genes, TRUE)
  query_up <- head(query_genes, ngenes)
  query_dn <- tail(query_genes, ngenes)
  
  
  res <- apply(drug_info, 2, function(x) {
    x <- sort(x, TRUE)
    x <- c(head(x, n), tail(x, n))
    
    x_up <- intersect(names(x), names(query_up))
    x_dn <- intersect(names(x), names(query_dn))
    
    sum(x[x_up]) - sum(x[x_dn])
  })
  
  res <- sort(res, TRUE)
  return(res)
}


#' Run KS Query
#'
#' @inheritParams query_cor
#' @inherit query_cor return
#' @param drug_prl Matrix of PRLs to query against. Rows names are genes, column names are signature names.
#' @export
#'
query_ks <- function(query_genes, drug_prl, ngenes = 100) {
  
  
  # same genes
  query_genes <- query_genes[names(query_genes) %in% row.names(drug_prl)]
  
  # Vj's
  Vj_up <- sort(query_genes, TRUE)[1:ngenes]
  Vj_up <- drug_prl[names(Vj_up), ]
  
  Vj_dn <- sort(query_genes)[1:ngenes]
  Vj_dn <- drug_prl[names(Vj_dn), ]
  
  n <- nrow(drug_prl)
  
  # a's
  a_up <- apply(Vj_up, 2, function(Dr_up) {
    max((1:ngenes/ngenes) - (sort(Dr_up)/n))
  })
  
  a_dn <- apply(Vj_dn, 2, function(Dr_dn) {
    max((1:ngenes/ngenes) - (sort(Dr_dn)/n))
  })
  
  # b's
  b_up <- apply(Vj_up, 2, function(Dr_up) {
    max((sort(Dr_up)/n) - ((0:(ngenes-1))/ngenes))
  })
  
  b_dn <- apply(Vj_dn, 2, function(Dr_dn) {
    max((sort(Dr_dn)/n) - ((0:(ngenes-1))/ngenes))
  })
  
  #ks's
  ks_up <- ifelse(a_up > b_up, a_up, -b_up)
  ks_dn <- ifelse(a_dn > b_dn, a_dn, -b_dn)
  
  si <- ks_up - ks_dn
  p <- max(si)
  q <- min(si)
  
  ks_opp <- sign(ks_up) != sign(ks_dn)
  
  # S
  S <- numeric(length(si))
  names(S) <- colnames(drug_prl)
  
  S[ks_opp & si > 0] <- si[ks_opp & si > 0]/p
  S[ks_opp & si < 0] <- -(si[ks_opp & si < 0]/q)
  
  S <- S[order(S, ks_up, decreasing=TRUE)]
  
  return(S)
}


#' Run Euclidean Distance Query
#'
#' @inheritParams query_cor
#' @inherit query_cor return
#' @export
#'
query_euclidean <- function(query_genes, drug_info, ngenes = 100) {
  
  drug_info <- drug_info[row.names(drug_info) %in% names(query_genes), ]
  query_genes <- query_genes[names(query_genes) %in% row.names(drug_info)]
  
  query_genes <- sort(query_genes, TRUE)
  query_genes <- c(head(query_genes, ngenes), tail(query_genes, ngenes))
  
  drug_info <- drug_info[names(query_genes), ]
  
  res <- apply(drug_info, 2, function(x) dist(rbind(x, query_genes)))
  
  res <- sort(res)
  return(res)
}


#' Run Cosine Distance Query for Characteristic Direction Method
#'
#' @inheritParams query_cor
#' @inherit query_cor return
#' @export
#'
query_cd <- function(query_genes, drug_info, ngenes = 100) {
  
  
  # use only common genes
  query_genes <- query_genes[names(query_genes) %in% row.names(drug_info)]
  
  # top 100 up/down genes
  query_genes <- sort(query_genes, TRUE)
  query_genes <- c(head(query_genes, ngenes), tail(query_genes, ngenes))
  drug_info   <- drug_info[names(query_genes), ,drop = FALSE]
  
  # cosine distance
  dist   <- 1 - apply(drug_info, 2, function(col) lsa::cosine(query_genes, col))
  
  return(sort(dist))
}

