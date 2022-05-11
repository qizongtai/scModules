#' @title Find the Cells with the highest CNA signal
#' @description Find Cells in the top nth quantile for CNA Signal values
#' @param cna a matrix of cell rows by cell columns containing CNA values
#' @param cell.quantile calculate CNA measures including only top / "hotspot" cells according to their squared CNA values across all genes. Value between 0 and 1 denoting the quantile of cells to include. 
#' @return cell names in the top nth quantile, where n is specified via <cell.quantile>
#' @rdname cnaHotspotCells
#' @export 
cnaHotspotCells = function(cna, cell.quantile, gene.quantile = NULL) {
  cna = as.matrix(cna)
  stopifnot(is.numeric(cell.quantile))
  if (!is.null(gene.quantile)) {
    stopifnot(is.numeric(gene.quantile))
    cna = cna[cnaHotspotGenes(cna, gene.quantile = gene.quantile), ]
  }
  msq = colMeans(cna^2)
  names(msq)[msq >= quantile(msq, cell.quantile)]
}

#' @title Find the Genes with the highest CNA signal
#' @description Find Genes in the top nth quantile for CNA Signal values
#' @param cna a matrix of gene rows by cell columns containing CNA values
#' @param gene.quantile calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include. 
#' @return gene names in the top nth quantile, where n is specified via <gene.quantile>
#' @rdname cnaHotspotGenes
#' @export 
cnaHotspotGenes = function(cna, gene.quantile, cell.quantile = NULL) {
  cna = as.matrix(cna)
  stopifnot(is.numeric(gene.quantile))
  if (!is.null(cell.quantile)) {
    stopifnot(is.numeric(cell.quantile))
    cna = cna[, cnaHotspotCells(cna, cell.quantile = cell.quantile)]
  }
  msq = rowMeans(cna^2)
  names(msq)[msq >= quantile(msq, gene.quantile)]
}

#' @title Get unique sample names
#' @description Extract unique sample names from cell ids
#' @param x character vector of cell ids
#' @param sep separator Default: '-|_'
#' @param pos position of sample name given separator. First position is 1. Default: 1
#' @param max.nchar maximum number of characters to take as sample name (starting from the first). Default: NULL
#' @param replace a list of character vectors to replace. Takes the form list(c(old, new), c(old, new)). Default: NULL
#' @return character vector of unique sample names
#' @rdname unique_sample_names
#' @export 
#' @importFrom stringr str_split
unique_sample_names = function(x, sep = "-|_", pos = 1, max.nchar = NULL, replace = NULL) {
  #samples = sapply(stringr::str_split(x, sep), `[`, pos)
  samples = stringr::str_split(x, sep)
  samples = sapply(samples, `[`, pos, simplify = F)
  sep = sapply(stringr::str_split(sep, "\\|"), `[`, 1)
  samples = sapply(samples, paste0, collapse = sep, simplify = F)
  samples = as.character(unlist(samples))
  
  if (!is.null(max.nchar)) {
    samples = sapply(samples, function(sample) {
      if (nchar(sample) <= max.nchar) sample
      else substr(sample, 1, max.nchar)})
  }
  
  samples = unique(samples)
  
  if (!is.null(replace)) {
    for (i in 1:length(replace)) {
      r.old = paste0(replace[[i]], collapse = '|')
      r.new = replace[[i]][[2]]
      samples = stringr::str_replace(samples, r.old, r.new)
    }
  }
  
  unique(samples)
}

#' @title Unlist, keeping original list or element names
#' @description Unlist, keeping original list or element names
#' @param L list to flatten
#' @param nested.names logical; keep nested list names rather than expanding list names. Default: FALSE
#' @return a vector of the same length as the combined lengths of list elements. Names will either be the list names replicated, or, if nested.names is TRUE, the original list element names will be kept.
#' @seealso 
#'  \code{\link[stats]{setNames}}
#' @rdname Unlist
#' @export 
#' @importFrom stats setNames
Unlist = function(L, nested.names = FALSE) {
  if (nested.names) {
    Names = unlist(sapply(L, names), use.names = F)
  } else {
    Names = rep(names(L), lengths(L))
  }
  stats::setNames(unlist(L), Names)
}


