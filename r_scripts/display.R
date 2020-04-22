library(ggplot2)
library(reshape2)

plot_nmds <- function(data, metadata=NULL, 
                      color.fact=NULL, 
                      shape.fact=NULL, 
                      size=2, text.size=10,
                      alpha=0.8,
                      row.fact=NULL, row.fact2=NULL,
                      col.fact=NULL, col.fact2=NULL,
                      output=NULL) {
  
  if (is.null(shape.fact)) {
    esth <- aes(x=NMDS1,y=NMDS2, colour=data[[color.fact]])
  } else {
    esth <- aes(x=NMDS1,y=NMDS2, shape=data[[shape.fact]], colour=data[[color.fact]])
  }
  
  p <- ggplot(data, esth) + 
    geom_point(size=size, alpha=alpha) +
    coord_equal() +
    ## stat_ellipse(level=0.9, show.legend=FALSE, linetype=4, size=0.8) +
    labs(color=factor(color.fact), shape=factor(shape.fact)) +
    theme_linedraw(base_size=text.size)
    ## theme(legend.title = element_text(size=text.size+2, face=2))
  
  if (color.fact == 'Flowpath') {
    colors <- c("Blue"="#3399FF", 
                "Green"="#009900", 
                "Orange"="#FFAE41", 
                "Purple"="#660099", 
                "Red"="#FF0000",
                "White"="#808080")
    p <- p + scale_color_manual(breaks=names(colors), values=unname(colors))
  } else {
    p <- p + scale_color_discrete(color.fact)
  }
  
  formula <- NULL
  if(!is.null(row.fact2)) {
    formula <- paste(row.fact, '+', row.fact2, '~', col.fact2, '+', col.fact)
  } else if(!is.null(col.fact2)) {
    formula <- paste(row.fact, '~', col.fact2, '+', col.fact)
  } else if(!is.null(row.fact) && !is.null(col.fact)) {
    formula <- paste(row.fact, '~', col.fact)
  } else if(!is.null(row.fact)) {
    formula <- paste(row.fact, '~', '.')
  } else if(!is.null(col.fact)) {
    formula <- paste('.', '~', col.fact)
  }
  
  p <- p + facet_grid(as.formula(formula))
  
  ggsave(output, scale=2)
  
  return(p)
}


plot_env <- function(ef, outdir) {
  ef_df <- melt(cbind(p.value=ef$factors$pvals, r2=ef$factors$r))
  colnames(ef_df) <- c('Factor', 'Metric', 'Score')
  levels(ef_df$Metric) <- c('p-value', 'r2 correlation score')
  
  p <- ggplot(ef_df, aes(x=Factor, y=Score, fill=Factor)) +
              facet_grid(Metric ~ ., scales='free') +
              geom_bar(position="dodge", stat="identity", color='black', alpha=0.7, width=0.5) +
      scale_fill_brewer(palette="Dark2") + 
      theme_linedraw()

  ggsave(sprintf('%s/figures/bray_env-fit.pdf', outdir))
  
  return(p)
}
