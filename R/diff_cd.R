#' Characteristic Direction Differential Expression
#' @export
diff_cds <- function(esets, data_dir = getwd(), annot = "SYMBOL",
                     prev_anals = list(NULL)) {
  
  esets <- esets[names(esets) %in% names(prev_anals)]
  prev_anals <- prev_anals[names(esets)]
  anals <- list()
  
  # check for annot column
  chk <- sapply(esets, function(x) annot %in% colnames(Biobase::fData(x)))
  
  if (FALSE %in% chk) {
    stop(annot, " column in fData missing for esets: ",
         paste(names(which(!chk)), collapse = ", "))
  }
  
  for (i in seq_along(esets)) {
    
    eset <- esets[[i]]
    gse_name <- names(esets)[i]
    prev_anal <- prev_anals[[i]]
    
    gse_folder <- strsplit(gse_name, "\\.")[[1]][1]  # name can be "GSE.GPL"
    gse_dir <- paste(data_dir, gse_folder, sep = "/")
    
    cat('Working on', paste0(gse_name, ':'), '\n')
    
    # select contrasts
    cons <- crossmeta:::add_contrasts(eset, gse_name, prev_anal)
    contrasts <- cons$contrasts
    
    # remove replicates and annotate with human SYMBOL
    dups <- crossmeta:::iqr_replicates(cons$eset, annot=annot)
    eset <- dups$eset
    
    pdata <- Biobase::pData(eset)
    groups <- pdata$group
    
    diff_cd <- list()
    # CD for each contrast
    for (con in contrasts) {
      
      cat('    ', paste0(con, '...', '\n'))
      
      # contrast levels
      ctrl <- gsub('^.+?-', '', con)
      test <- gsub('-.+?$', '', con)
      
      # control and experimental matrices
      ctrl <- Biobase::exprs(eset)[, groups == ctrl, drop=FALSE]
      expm <- Biobase::exprs(eset)[, groups == test, drop=FALSE]
      
      genes <- row.names(ctrl)
      
      #run chdir analysis
      con_name <- paste(gse_name, con, sep="_")
      diff_cd[[con_name]] <- tryCatch(chdir(ctrl, expm, genes)[, 1],
                                      error = function(e) return(NULL))
    }
    
    # save results for eset
    save_name <- paste(gse_name, "diff_cd", tolower(annot), sep = "_")
    save_name <- paste0(save_name, ".rds")
    saveRDS(diff_cd, file = paste(gse_dir, save_name, sep = "/"))
    
    # store results for eset
    anals <- c(anals, diff_cd)
  }
  
  return(list(contrasts = anals))
}

# used by diff_cds
chdir <- function(ctrl,expm,genes,r=1) {
  # This function caclulates the characteristic direction for a gene expression dataset.
  #  	ctrl: control gene expressoion data, a matrix object
  #  	expm: experiment gene expression data, a matrix object
  #  	b: return value, a vector of n-components, representing the characteristic
  #          direction of the gene expression dataset. n equals to the number of genes in the
  #          expression dataset. b is also a matrix object. b is sorted by its components'
  #          absolute values in descending order.
  #  	r: regularized term. A parameter that smooths the covariance matrix and reduces
  #          potential noise in the dataset. The default value for r is 1, no regularization.
  #
  #       For the input matrix rows are genes and columns are gene expression profiles.
  #       r is the regulization term ranging [0,1]. b is the characteristic direction.
  #       ctrl(control) and expm(experiment) matrices should have the same number
  #       of genes(rows).
  #
  #       Author: Qiaonan Duan
  #       Ma'ayan Lab, Icahn School of Medicine at Mount Sinai
  #       Jan.13, 2014
  #
  #		Add gene symbols to results. Apr. 4, 2014
  
  
  
  if(dim(ctrl)[1]!=dim(expm)[1]){
    stop('Control expression data must have equal number of genes as experiment expression data!')
  }
  
  if(any(is.na(ctrl))||any(is.na(expm))){
    stop('Control expression data and experiment expression data have to be real numbers. NA was found!')
  }
  
  
  # There should be variance in expression values of each gene. If
  # gene expression values of a gene are constant, it would dramatically
  # affect the LDA caculation and results in a wrong answer.
  constantThreshold <- 1e-5;
  ctrlConstantGenes <- diag(var(t(ctrl))) < constantThreshold
  expmConstantGenes <- diag(var(t(expm))) < constantThreshold
  
  # if (any(ctrlConstantGenes)){
  #     errMes <- sprintf('%s row(s) in control expression data are constant. Consider Removing the row(s).',paste(as.character(which(ctrlConstantGenes)),collapse=','))
  #     stop(errMes)
  # }else if(any(expmConstantGenes)){
  #     errMes <- sprintf('%s row(s) in experiment expression data are constant. Consider Removing the row(s).',paste(as.character(which(expmConstantGenes)),collapse=','))
  #     stop(errMes)
  # }
  
  # place control gene expression data and experiment gene expression data into
  # one matrix
  combinedData <- cbind(ctrl,expm)
  
  # get the number of samples, namely, the total number of replicates in  control
  # and experiment.
  dims <- dim(combinedData)
  samplesCount <- dims[2]
  
  # the number of output components desired from PCA. We only want to calculate
  # the chdir in a subspace that capture most variance in order to save computation
  # workload. The number is set 20 because considering the number of genes usually
  # present in an expression matrix 20 components would capture most of the variance.
  componentsCount <- min(c(samplesCount-1,20))
  
  
  # use the nipals PCA algorithm to calculate R, V, and pcvars. nipals algorithm
  # has better performance than the algorithm used by R's builtin PCA function.
  # R are scores and V are coefficients or loadings. pcvars are the variances
  # captured by each component
  pcaRes <- nipals(t(combinedData),componentsCount,1e5,1e-4)
  R <- pcaRes$T
  V <- pcaRes$P
  pcvars <- pcaRes$pcvar
  
  
  # we only want components that cpature 95% of the total variance or a little above.
  # cutIdx is the index of the compoenent, within which the variance is just equal
  # to or a little greater than 95% of the total.
  cutIdx <- which(cumsum(pcvars)>0.95)
  if(length(cutIdx)==0){
    cutIdx <- componentsCount
  }else{
    cutIdx <- cutIdx[1]
  }
  
  # slice R and V to only that number of components.
  R <- R[,1:cutIdx]
  V <- V[,1:cutIdx]
  
  # the difference between experiment mean and control mean.
  meanvec <- rowMeans(expm) - rowMeans(ctrl)
  
  
  # all the following steps calculate shrunkMats. Refer to the ChrDir paper for detail.
  # ShrunkenMats are the covariance matrix that is placed as denominator
  # in LDA formula. Notice the shrunkMats here is in the subspace of those components
  # that capture about 95% of total variance.
  Dd <- t(R)%*%R/samplesCount
  Dd <- diag(diag(Dd))
  sigma <- mean(diag(Dd))
  shrunkMats <- r*Dd + sigma*(1-r)*diag(dim(R)[2])
  
  # The LDA formula.
  #  V%*%solve(shrunkMats)%*%t(V) transforms the covariance matrix from the subspace to full space.
  b <- V%*%solve(shrunkMats)%*%t(V)%*%meanvec
  
  # normlize b to unit vector
  b <- b*as.vector(sqrt(1/t(b)%*%b))
  
  # sort b to by its components' absolute value in decreasing order and get the
  # sort index
  sortRes <- sort(abs(b),decreasing=TRUE,index.return=TRUE)
  
  # sort b by the sort index
  bSorted <- as.matrix(b[sortRes$ix])
  # sort genes by the sort index
  genesSorted <- genes[sortRes$ix]
  # assign genesSorted as the row names of bSorted
  rownames(bSorted) <- genesSorted
  
  # return bSorted
  bSorted <- bSorted
}
nipals <- function(X,a,it=10,tol=1e-4) {
  #fct nipals calculates the principal components of a given data matrix X according to
  #the NIPALS algorithm (Wold).
  #X...data matrix, a...number of components,
  #it...maximal number of iterations per component,
  #tol...precision tolerance for calculation of components
  
  Xh <- scale(X,center=TRUE,scale=FALSE)		#mean-centering of data matrix X
  nr <- 0
  T <- NULL
  P <- NULL
  pcvar <- NULL
  varTotal <- sum(diag(var(Xh)))
  currVar <- varTotal
  
  for (h in 1:a){
    th <- Xh[,1]		#starting value for th is 1st column of Xh
    ende <- FALSE
    #3 inner steps of NIPALS algorithm
    while (!ende){
      nr <- nr+1
      
      # the result of matrix multiplication operation (%*%) is a matrix of a single
      # valule. A matrix cannot multiply another using scalar multiplication (*).
      # as.vector convert a value of class matrix to a value of class double.
      # (A'*B)' = B'*A
      ph <- t((t(th)%*%Xh) * as.vector(1/(t(th)%*%th)))	#LS regression for ph
      ph <- ph * as.vector(1/sqrt(t(ph)%*%ph))		#normalization of ph
      thnew <- t(t(ph)%*%t(Xh) * as.vector(1/(t(ph)%*%ph)))	#LS regression for th
      prec <- t(th-thnew)%*%(th-thnew)	#calculate precision
      # cat("actual precision: ",sqrt(prec),"\n")
      th <- thnew	#refresh th in any case
      #check convergence of th
      if (prec <= (tol^2)) {
        ende <- TRUE
      }
      else if (it <= nr) {	#too many iterations
        ende <- TRUE
        cat("\nWARNING! Iteration stop in h=",h," without convergence!\n\n")
      }
    }
    Xh <- Xh-(th%*%t(ph))	#calculate new Xh
    T <- cbind(T,th)	#build matrix T
    P <- cbind(P,ph)	#build matrix P
    oldVar <- currVar
    currVar <- sum(diag(var(Xh)))
    pcvar <- c(pcvar,(oldVar-currVar)/varTotal)
    nr <- 0
  }
  list(T=T,P=P,pcvar=pcvar)
}

#' Load Characteristic Direction Differential Expression Signatures
#' @export
load_diff_cds <- function(gse_names, data_dir = getwd(), annot = "SYMBOL") {
  
  anals <- list()
  for (gse_name in gse_names) {
    gse_dir <- file.path (data_dir, gse_name)
    
    # get paths
    pattern <- paste0(paste(".*_diff_cd", tolower(annot), sep="_"), ".rds")
    anal_paths <- list.files(gse_dir, pattern, full.names = TRUE)
    
    # load each diff path
    # multiple if more than one platform per GSE)
    for (path in anal_paths) {
      anals <- c(anals, readRDS(path))
    }
  }
  return (list(contrasts = anals))
}
