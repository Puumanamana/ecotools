## library(configr)
library(vegan)
library(pairwiseAdonis)
library(dplyr)
library(stringr)

library(usedist)
library(Hmisc)

library(doParallel)
registerDoParallel(cores=50)

source('parser.R')
source('loader.R')
source('explore_factor_influence.R')
source('display.R')
source('stats.R')

args <- parse_args()

raw <- load_data_h5(args$h5, args$tree)
raw <- filter_data(raw, group=args$groupby)
data <- transform_sample_counts(raw, args$norm_fn)

if (args$distance == 'unifrac') {
    distances <- UniFrac(data, weighted=FALSE, normalized=TRUE, parallel=TRUE)
} else if (args$distance == 'wunifrac') {
    distances <- UniFrac(data, weighted=TRUE, normalized=TRUE, parallel=TRUE)
} else {
    distances <- distance(data, 'bray')
}

#################################
### Step 1: environmental fit ###
#################################

abd.all <- as(t(otu_table(data)), 'matrix')
nmds.bc.all <- metaMDS(abd.all, k=3, trymax=500, parallel=60)
stressplot(nmds.bc.all)

## cols <- c('Year', 'Month', 'Lifestyle', 'Aquifer', 'Flowpath', 'Divide')
meta.all <- as(sample_data(data), 'data.frame')[, args.factors]

ef <- envfit(ord=nmds.bc.all,
             env=meta.all, 
             permutations=args$nperms, 
             na.rm=T)
plot_env(ef, outdir)

## Check dispersion in each group
dists <- vegdist(abd.all)
model.bdisper <- betadisper(dists, meta.all$Year)
plot(model.bdisper, label=T, ellipse=TRUE, hull=FALSE, label.cex=.5)

####################################################
### Step 2: Visualize factor influence with NMDS ###
####################################################
## We define a condition as a tuple (Site_type, Year, Month). 
## We repeat the analysis for each condition

## Number of samples per condition
strata.fact.names <- c('Year', 'Month')
meta <- sample_data(data)
strata <- interaction(meta[, strata.fact.names], drop=TRUE)
strata.lvls <- sort(levels(strata))
print(table(strata))

## Process all conditions separately ##
nmds.bc <- as.data.frame(matrix(0, nrow=length(strata), ncol=2))
rownames(nmds.bc) <- rownames(meta)
colnames(nmds.bc) <- c('NMDS1', 'NMDS2')

for(stratum in levels(strata)) {
    print(sprintf("Processing %s", stratum))
    ## Subset the phyloseq data with the corresponding samples
    ps.cond <- subset_samples(data, strata==stratum)
    ps.cond <- prune_taxa(taxa_sums(ps.cond) > 0, ps.cond)
  
    if(nsamples(ps.cond) < 10) {
      print(sprintf("Not enough samples for %s (%s samples). Skipping.", stratum, nsamples(ps.cond)))
      next
    }
    ## Compute the NMDS components
    ## dist.cond <- distance(ps.cond, 'bray')
    ## dist.cond <- UniFrac(ps.cond, weighted=TRUE, normalized=TRUE, parallel=TRUE)
    components.cond <- scores(
        metaMDS(t(otu_table(ps.cond)), k=2, trymax=200, parallel=50)
    )

    nmds.bc[rownames(components.cond), 'NMDS1'] <- components.cond[, 1]
    nmds.bc[rownames(components.cond), 'NMDS2'] <- components.cond[, 2]
}

## Visualize
nmds.bc <- nmds.bc[rownames(meta), ]
nmds.bc <- cbind(nmds.bc, meta[!(colnames(meta) %in% colnames(nmds.bc))])

write.csv(nmds.bc, sprintf('%s/nmds_components_splitby-%s.csv', table_dir, paste0(strata.fact.names, collapse='-')))

for (factor in c('Flowpath', 'Aquifer')){
  output <- sprintf("%s/NMDSbray_by-%s_splitby-%s.pdf", 
                    fig_dir, factor, paste0(strata.fact.names, collapse='-'))

  plot_nmds(nmds.bc, size=2, text.size=10,
            color.fact=factor, 
            row.fact='Year',
            col.fact='Month',
            output=output)
}

################################################
### Step 3: Statistical analysis (PERMANOVA) ###
################################################

empty_table <- function(rows, ref=NULL) {
  data <- as.data.frame(matrix(nrow=length(rows), ncol=5))
  rownames(data) <- rows
  colnames(data) <- c('R2', 'p.adjusted', 'F.Model', 'n1', 'n2')
  
  return(data)
}

for (factor.name in c('Flowpath', 'Aquifer')) {
  f.lvls <- levels(meta[[factor.name]])
  labels <- apply(combn(f.lvls, 2), 2, paste0, collapse=" vs ")
  
  results <- lapply(labels, function(x) empty_table(strata.lvls, ref=x))
  names(results) <- labels
  
  for(stratum in levels(strata)) {
      print(sprintf("Processing %s (factor: %s)", stratum, factor.name))
      ## Subset the phyloseq data with the corresponding samples
      ps.cond <- subset_samples(data, strata==stratum)
      ps.cond <- prune_taxa(taxa_sums(ps.cond) > 0, ps.cond)
    
      abd.cond <- t(as(otu_table(ps.cond), 'matrix'))
      meta.cond <- as(sample_data(ps.cond), 'data.frame')[[factor.name]]
      
      freqs <- table(meta.cond)
    
      if (length(unique(as.character(meta.cond))) > 1) {
          model <- pairwise.adonis(abd.cond[!is.na(meta.cond),], 
                                   meta.cond[!is.na(meta.cond)], 
                                   perm=N_PERMS)
          rownames(model) <- as.character(model$pairs)
          print(model)
      
          for (comp in rownames(model)) {
              res <- c(model[comp, c('R2', 'p.adjusted', 'F.Model')])
            
              if (!(comp %in% names(results))) {
                comp <- gsub('(.*) vs (.*)', '\\2 vs \\1', comp)
              }
            
              lvls <- unlist(strsplit(comp, ' vs '))
              res$n1 <- freqs[lvls[1]]
              res$n2 <- freqs[lvls[2]]

              results[[comp]][stratum, ] <- res
          }
      } else {
          print(sprintf('Only one level for %s', stratum))
      }
  }
  
  for (label in labels) {
      write.csv(results[[label]][complete.cases(results[[label]]), ], 
                sprintf('%s/adonis_%s_%s.csv', table_dir, factor.name, gsub(' ', '-', label)), 
                quote=F, row.names=T)
  }
}

