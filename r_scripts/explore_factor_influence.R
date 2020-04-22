library(vegan)

source('parser.R')

fit_environment <- function(data, factors, nperms=1e5, k=3, parallel=60) {
    abd.all <- as(t(otu_table(data)), 'matrix')
    nmds.bc.all <- metaMDS(abd.all, k=3, trymax=500, parallel=60)

    meta.all <- as(sample_data(data), 'data.frame')[, factors]

    ef <- envfit(ord=nmds.bc.all,
                 env=meta.all,
                 permutations=nperms,
                 na.rm=T)

    plot_env(ef, outdir)

    ## Check dispersion in each group
    dists <- vegdist(abd.all)
    model.bdisper <- betadisper(dists, meta.all$Year)
    plot(model.bdisper, label=T, ellipse=TRUE, hull=FALSE, label.cex=.5)
}
