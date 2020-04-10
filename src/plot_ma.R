# Plot comparison of meta-analysis to contrasts

library(ggplot2)
library(gridExtra)
library(RColorBrewer)

setwd("~/Documents/Batcave/GEO/1-meta")

# paths_sva <- list.files(pattern = "ma_res_sva.rds", full.names = TRUE, recursive = TRUE)
paths  <- list.files(pattern = "ma_res.rds",  full.names = TRUE, recursive = TRUE)

# names(paths_sva) <- stringr::str_extract(paths_sva, "[A-Z]+")
names(paths)  <- stringr::str_extract(paths,  "[A-Z]+")

drugs <- c("^tretinoin", "dexamethasone", "doxorubicin", "LY-294002", "metformin",
           "progesterone", "quercetin", "sirolimus", "resveratrol","vorinostat")

# load results
# results_sva <- lapply(paths_sva, function(path) {readRDS(path)})
results  <- lapply(paths,  function(path) {readRDS(path)})

# Contrasts ------------------------------

# overlap from individual contrasts
consov <- mapply(function(res, drug) {

     lapply(res$cons, function(con) {
       con[grep(drug, names(con))]
    })

}, results, drugs)

consov <- unlist(consov)


# results from individual contrasts
cons <- mapply(function(res, drug) {

    lapply(res$cons, function(con) {
        grep(drug, names(con))
    })

}, results, drugs)

cons <- unlist(cons)


# Meta-Analysis ------------------------------


# overlap from meta-analyses
metaov <- mapply(function(res, drug) {
  res$meta[grep(drug, names(res$meta))]
}, results, drugs)


# results from meta-analyses
meta <- mapply(function(res, drug) {
    grep(drug, names(res$meta))
}, results, drugs)

meta <- unlist(meta)



# Violin ------------------------------

covdf <- data.frame(overlap = consov,
                    category = rep("Contrasts", length(consov)),
                    class = rep("no sva", length(consov)))

movdf <- data.frame(overlap = metaov,
                    category = rep("Meta-Analyses", length(metaov)),
                    class = rep("no sva", length(metaov)))

# covdf_sva <- data.frame(overlap = consov_sva,
#                         category = rep("Contrasts (sva)", length(consov_sva)),
#                         class = rep("sva", length(consov_sva)))

# movdf_sva <- data.frame(overlap = metaov_sva,
#                         category = rep("Meta-analyses (sva)", length(metaov_sva)),
#                         class = rep("sva", length(metaov_sva)))


ovdf <- rbind(movdf, covdf)


# save json
whiskerdata <- data.frame(type = as.character(ovdf$category), correlation = ovdf$overlap, stringsAsFactors = FALSE)
jsonlite::write_json(whiskerdata, "whiskerdata.json")


ovpl <- ggplot(ovdf, aes(category, overlap)) +
    geom_violin(aes(colour = category),
                draw_quantiles = 0.5,
                adjust = 0.5, size=1) +
    geom_jitter(width = 0.30, alpha = 0.3) +
    labs(x="", y= "Similarity") +
    scale_x_discrete(breaks=c("", "")) +
    scale_colour_manual(values=c('#4582ec', '#d9534f')) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_line(colour = "#dddddd", size = 0.6),
          plot.background = element_rect(fill = "transparent",colour = NA),
          plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
          axis.text = element_text(size=14),
          axis.text.x = element_text(margin=margin(5,5,10,5,"pt")),
          axis.text.y = element_text(margin=margin(5,5,10,5,"pt")),
          axis.title = element_text(size=16),
          legend.background = element_blank(),
          legend.text = element_text(size=16),
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.position = 'bottom')


ovpl


# ------------------------------ ROCR

get_rates <- function(res) {

    # number of true positives (one per con/ma)
    tp <- length(res)

    # number of false positives (1308 per con/ma)
    fp <- tp * 1308

    # for each posible position, add to tpr and fpr
    tpr <- c(0)
    fpr <- c(0)

    for (i in 1:1309) {

        # add to tpr for results that are correct at this position
        tpr <- c(tpr, tail(tpr, 1) + (sum(res == i) / tp))

        # add to fpr for results that are incorrect at this position
        fpr <- c(fpr, tail(fpr, 1) + (sum(res != i) / fp))
    }
    return(list(tpr = tpr, fpr = fpr))
}

# contrasts: (auc = 0.812)
crt <- get_rates(cons)
crtdf <- data.frame(crt, Approach = "Contrasts")
MESS::auc(x = crtdf$fpr, y = crtdf$tpr)

# meta analyses: (auc = 0.956)
mrt <- get_rates(meta)
mrtdf <- data.frame(mrt, Approach = "Meta Analyses")
MESS::auc(x = mrtdf$fpr, y = mrtdf$tpr)


# contrasts: (auc = 0.730)
# crt_sva <- get_rates(cons_sva)
# crtdf_sva <- data.frame(crt_sva, Approach = "Contrasts (sva)")
# MESS::auc(x = crtdf_sva$fpr, y = crtdf_sva$tpr)

# meta analyses: (auc = 0.922)
# mrt_sva <- get_rates(meta_sva)
# mrtdf_sva <- data.frame(mrt_sva, Approach = "Meta-analyses (sva)")
# MESS::auc(x = mrtdf_sva$fpr, y = mrtdf_sva$tpr)

scaleFUN <- function(x) sprintf("%.1f", x)

# plot together
rtdf <- rbind(mrtdf, crtdf)

# save json
rocdata <- list(cons = as.matrix(crtdf[, c('fpr', 'tpr')]),
                meta = as.matrix(mrtdf[, c('fpr', 'tpr')]))

jsonlite::write_json(rocdata, "rocdata.json")

rtpl <- ggplot(rtdf) +
    geom_line(aes(y = tpr, x = fpr, colour = Approach),
              size=1, show.legend = FALSE) +
    scale_colour_manual(values=c('#4582ec', '#d9534f')) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
  scale_y_continuous(labels=scaleFUN) +
  scale_x_continuous(labels=scaleFUN) +
    theme(
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent",colour = NA),
      panel.grid.major = element_line(colour = "#dddddd", size = 0.6),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
      axis.text.x = element_text(margin=margin(5,5,10,5,"pt")),
      axis.text.y = element_text(margin=margin(5,5,10,5,"pt")),
      axis.title = element_text(size=16),
      axis.text = element_text(size=14))


rtpl

# Largest sample size ----

# number of queries (contrasts + 1)
num_cons <- sapply(results, function(res) length(res$cons))

# rank data frame
cons_rank <- data.frame(rank = cons, special = 'none', drug = rep(names(meta), num_cons), stringsAsFactors = FALSE)



# largest sample size
maxn <- sapply(results, function(res) names(which.max(res$nsamples)))
maxn <- paste(names(maxn), unname(maxn), sep='.')

cons_rank[maxn, 'special'] <- 'Largest Study'

# add meta-analysis
ma_rank <- data.frame(rank = meta, special = 'Meta Analysis  ', drug = names(meta))
cons_rank <- rbind(ma_rank, cons_rank)

# plot
rkpl <- ggplot(cons_rank) +
    geom_jitter(aes(y = rank, x = drug, colour = special, alpha= special, size = special), width=0.15) +
    scale_alpha_manual(breaks=c('Meta Analysis  ', 'Largest Study'), values=c(1, 0.2, 1)) +
    scale_size_manual(breaks=c('Meta Analysis  ', 'Largest Study'), values=c(3, 1, 3)) +
    scale_y_log10(breaks=c(1, 10, 100, 1000)) +
    scale_colour_manual(breaks=c('Meta Analysis  ', 'Largest Study'), values=c('#4582ec', '#000000', '#d9534f')) +
    ylab("Rank") +
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#dddddd", size = 0.6),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
        axis.text.x = element_text(margin=margin(5,5,10,5,"pt")),
        axis.text.y = element_text(margin=margin(5,5,10,5,"pt")),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        axis.title.x = element_blank(),
        legend.title=element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size=16),
        legend.key=element_blank(),
        legend.position = 'bottom')


rkpl




# Legend Function ---------------------------------

get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

# Plot Together -------------------------

save_for_web <- function(fname, pl, rotx = FALSE) {

  ggsave(paste0(fname, '-xl.svg'), pl, width = 908/72, height = 475/72, bg="white", dpi = 72)
  ggsave(paste0(fname, '-lg.svg'), pl, width = 773/72, height = 405/72, bg="white", dpi = 72)
  ggsave(paste0(fname, '-md.svg'), pl, width = 721/72, height = 377/72, bg="white", dpi = 72)
  ggsave(paste0(fname, '-sm.svg'), pl, width = 599/72, height = 313/72, bg="white", dpi = 72)
  
  if (rotx)
    pl <- pl + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(paste0(fname, '-xs.svg'), pl, width = 400/72, height = 250/72, bg="white", dpi = 72)
    
}

# get legend
# legend <- get_legend(rtpl)

# remove legend
# rtpl  <- rtpl  + theme(legend.position="none")

mult <- arrangeGrob(ovpl, rtpl, ncol=2, widths=c(1, 1))

save_for_web('ma', mult)
save_for_web('malg', rkpl, rotx = TRUE)


# save at 1400 width svg
ggsave("plot_ma.svg", mult, width = 10, height = 5.233, bg="transparent")

ggsave("plot_ma_lg.svg", rkpl, width = 10, height = 5.233, bg="transparent")

# save at 1260 width svg
ggsave("plot_ma_thumb.svg", mult, width = 10, height = 5.233, bg="#f8f8f8")
