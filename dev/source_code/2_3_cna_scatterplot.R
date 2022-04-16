#' @title Visualise Malignant and Non-Malignant Subsets
#' @description Visualise Malignant and Non-Malignant Subsets of cells. This is achieved by plotting, for each cell, its CNA signal over its CNA correlation. Please see `cna_signal` and `cna_cor` for details.
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param cor.method character string indicating the method to use for the pairwise correlations. E.g. 'pearson', 'spearman'. Default: 'pearson'
#' @param samples a character vector of sample names list of character vectors (cell/column names). Can be provided if the cna matrix contains multiple samples of cells, i.e. multiple tumours, such that the cell-groups correlations are calcualted for each groups/tumour in turn. Default: FALSE
#' @param gene.quantile calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include. Default: NULL
#' @param gene.quantile.for.corr as above but for CNA correlations specifically. Default: gene.quantile 
#' @param gene.quantile.for.signal as above but for CNA signal specifically. Default: gene.quantile
#' @param groups character vector of cell IDs to colour. Default: NULL
#' @param groups.col colour to colour groups by. Default: 'magenta'
#' @param border.col colour for point border. Default: 'black'
#' @param alpha colour transparency. Default: 0.3
#' @param hline y-intercept line. Default: NULL
#' @param vline x-intercept line. Default: NULL
#' @param refCells a character vector of cell ids to exclude from average CNA profile that each cell is correlated to. You can pass reference normal cell ids to this argument if these are known. Default: NULL
#' @param samples if CNA correlations should be calculated within cell subgroups, provide i) a list of cell id groups, ii) a character vector of sample names to groups cells by, iii) TRUE to extract sample names from cell ids and subsequently groups. Default: NULL
#' @param ... other arguments passed to base plot 
#' @return a base R plot. If return value is saved to a variable, instead returns data points for cna correlations and cna signal in list form.
#' @rdname cna_scatterplot
#' @export 
cna_scatterplot = function(cna,
                           cor.method = 'pearson',
                           cor.cell.quantile = NULL,
                           cor.gene.quantile.for.cells = NULL,
                           cor.gene.quantile = NULL,
                           cor.cell.quantile.for.genes = NULL,
                           signal.gene.quantile = NULL,
                           signal.cell.quantile.for.genes = NULL,
                           groups = NULL,
                           groups.col = discrete_colours[1:length(groups)],
                           border.col = 'black', 
                           alpha = 0.3,
                           cex = 0.8,
                           pch = 20,
                           hline = NULL,
                           vline = NULL,
                           refCells = NULL,
                           samples = NULL,
                           ...) {
  cors = cna_cor(cna,
                 cell.quantile = cor.cell.quantile,
                 gene.quantile.for.cells = cor.gene.quantile.for.cells,
                 gene.quantile = cor.gene.quantile,
                 cell.quantile.for.genes = cor.cell.quantile.for.genes,
                 refCells = refCells,
                 samples = samples)
  signals = cna_signal(cna,
                       gene.quantile = signal.gene.quantile,
                       cell.quantile.for.genes = signal.cell.quantile.for.genes,
                       refCells = refCells,
                       samples = samples)
  plot(cors,
       signals,
       xlab = 'CNA Correlation',
       ylab = 'CNA Signal',
       pch = 1,
       col = border.col,
       cex = cex,
       ...)
  
  if (!is.null(groups)) {
    if (!is.list(groups)) {
      groups = list(groups)
      groups = sapply(groups, function(x) x[x %in% colnames(cna)], simplify = F)
    }
    groups.col = scales::alpha(groups.col, alpha)
    
    .Points = function(cors, signals, groups, groups.col) {
      points(cors[groups],
             signals[groups],
             pch = pch,
             col = groups.col,
             cex = cex)
    }
    
    Map(.Points,
        groups = groups,
        groups.col = groups.col,
        MoreArgs = list(cors = cors, signals = signals))
  }
  
  if (!is.null(vline)) abline(v = vline, lty = 2)
  if (!is.null(hline)) abline(h = hline, lty = 2)
  
  return(invisible(list(cna.cor = cors, cna.signal = signals)))
}

