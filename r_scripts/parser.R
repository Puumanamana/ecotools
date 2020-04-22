library("argparse")

parse_args <- function() {
    ## create parser object
    parser <- ArgumentParser()

    parser$add_argument("--h5", type="character",
                        help="Path to h5 formatted data file")
    parser$add_argument("--tree", type="character",
                        help="Path to phylogenetic tree")    
    parser$add_argument("--otu-thresh", type="integer",
                        default=100,
                        help="OTU similarity threshold")
    parser$add_argument("--norm-fn", type="character",
                        help="Abundance normalization method",
                        choices=c("NULL", "assr", "ratio"),
                        default="NULL")
    parser$add_argument("--distance", type="character",
                        help="Metric or Dissimilarity index for analysis",
                        choices=c("bray", "unifrac", "wunifrac", "jaccard"),
                        default="bray")
    parser$add_argument("--factors", type="character",
                        help="Factor names in metadata to consider for analysis")
    parser$add_argument("--groupby", type="character",
                        help="Group OTU by some taxonomic level",
                        choices=c("NULL", "Genus", "Family", "Order", "Class", "Phylum"),
                        default="Genus")
    parser$add_argument("--nperms", type="integer", default=1e5,
                        help="Number of permutation for Permanova tests")
    parser$add_argument("--fig-dir", type="character", default="./figures",
                        help="Location to save figures")
    parser$add_argument("--table-dir", type="character", default="./tables",
                        help="Location to save tables")
    
    ## get command line options, if help option encountered print help and exit,
    ## otherwise if options not found on command line then set defaults, 
    args <- parser$parse_args()

    cat("\n")

    return(args)

    
