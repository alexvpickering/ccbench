library(crossmeta)
library(ccmap)

cmap_es <- readRDS('data/raw/cmap_es_ind.rds')
l1000_es <- readRDS("data/raw/l1000_es.rds")
l1000_fdr <- readRDS('data/raw/l1000_fdr.rds')


# common
cmap_drugs <- gsub('^([^_]+)_.+?$', '\\1', colnames(cmap_es))
l1000_drugs <- gsub('^([^_]+)_.+?$', '\\1', colnames(l1000_es))

cmap_drugs <- tolower(cmap_drugs)
cmap_drugs <- gsub('-| ', '', cmap_drugs)

l1000_drugs <- tolower(l1000_drugs)
l1000_drugs <- gsub('-| ', '', l1000_drugs)

cmap_unique <- unique(cmap_drugs)
l1000_unique <- unique(l1000_drugs)

common <- intersect(cmap_unique, l1000_unique)

# rename
colnames(l1000_es) <- colnames(l1000_fdr) <- l1000_drugs
colnames(cmap_es) <- cmap_drugs

# get query signatures/pvals
dprimes <- lapply(common, function(drug) l1000_es[, drug])
pvals <- lapply(common, function(drug) l1000_fdr[, drug])
names(dprimes) <- names(pvals) <- common

# run queries based on pvals and fixed number of genes
pval_res  <- query_pval_ngenes(dprimes, cmap_es, pvals)
fixed_res <- query_fixed_ngenes(dprimes, cmap_es)

# save results
saveRDS(pval_res, 'data/processed/l1000_vs_cmap_pvals.rds')
saveRDS(fixed_res, 'data/processed/l1000_vs_cmap_fixed.rds')
