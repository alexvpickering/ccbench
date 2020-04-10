library(ggplot2)

# MRR vs FDR cutoff for L1000 queries against CMAP02 ----

pval_results <- readRDS('data/processed/l1000_vs_cmap_pvals.rds')$resl


# p-value cutoffs and using all genes but queries excluded as for p-value cutoffs
ranks <- lapply(pval_results, get_ranks)
ranks_excluded <- lapply(ranks, function(x) {tmp <- ranks$`1`; tmp[is.na(x)] <- NA; tmp})

df <- data.frame(mrrs = sapply(ranks, mrr),
                 mrrs_excluded = sapply(ranks_excluded, mrr),
                 breaks = as.numeric(names(ranks)))


jpeg('output/l1000_vs_cmap_pvals_mrr.jpeg')
plot(df$breaks,
     df$mrrs,
     col = 'red',
     xlab = 'FDR Cutoff',
     ylab = 'MRR',
     at = seq(0, 1, by = .05),
     ylim = c(min(df$mrrs),max(df$mrrs_excluded)))

points(df$breaks, df$mrrs_excluded, col = 'blue')

mtext(side=3, line=3, at=-0.07, adj=0, cex=1, 'MRR vs FDR Cutoff for L1000 Queries Against CMAP02:')
mtext(side=3, line=2, at=-0.07, adj=0, cex=1, col = 'red', 'Genes and Queries Excluded')
mtext(side=3, line=1, at=-0.07, adj=0, cex=1,  col = 'blue', 'Only Queries Excluded')
dev.off()
