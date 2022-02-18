.cna_signal = function(cna) {
  sqmat = cna^2
  colMeans(sqmat)
}

#' @title Calculate the Means of Squared CNA Values
#' @description Calculate the Mean of Squared CNA Values
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param gene.quantile calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include. Default: NULL
#' @param refCells a character vector of cell ids (e.g. normal reference cell ids) to exclude from calculation of cnaHotspotGenes. Only relevant if gene.quantile is not NULL. Default: NULL
#' @param samples if cnaHotspotGenes should be calculated within cell subgroups, provide i) a list of cell id groups, ii) a character vector of sample names to group cells by, iii) TRUE to extract sample names from cell ids and subsequently group. Default: NULL
#' @param ... other arguments passed to unique_sample_names if samples = TRUE.
#' @return a numeric vector of CNA signal values or the Mean of Squared CNA values
#' @rdname cnaSignal
#' @export 
cna_signal = function(cna, 
                      gene.quantile = NULL, 
                      cell.quantile.for.genes = NULL, 
                      refCells = NULL, 
                      samples = NULL,
                      verbose = TRUE) {
  if(verbose) {message("running 'cna_signal'... return a numeric vector of CNA signal(the Mean of Squared CNA values)")}
  
  cna = as.matrix(cna)
  
  if (is.null(gene.quantile)) {return(.cna_signal(cna))}
  
  #cell filter
  if (is.null(refCells)) tmpcna = cna
  else tmpcna = cna[, !colnames(cna) %in% unlist(refCells)]
  
  #---single sample/group---#
  if (is.null(samples)) {
    genes = cnaHotspotGenes(tmpcna, gene.quantile = gene.quantile, cell.quantile = cell.quantile.for.genes)
    return(.cna_signal(cna[genes, ]))
  }
  
  #---multiple samples/groups---#
  if (isTRUE(samples)) {
    samples = unique_sample_names(colnames(cna), ...)
    message('Samples identified:\n', paste0(samples, collapse = '\n'))
  }
  if (is.character(samples)) {
    samples = split_by_sample_names(colnames(cna), samples = samples)
    samples = samples[lengths(samples) != 0]
  }
  stopifnot(is.list(samples))
  #gene filter
  if (!is.null(gene.quantile)) {
    tmpcnalist = sapply(samples, function(cells) tmpcna[, tmpcna %in% cells, drop = FALSE], simplify = F)
    genelist = sapply(tmpcnalist, cnaHotspotGenes, gene.quantile = gene.quantile, cell.quantile = cell.quantile.for.genes)
    rm(tmpcnalist)
  } else {
    genelist = replicate(length(samples), rownames(cna), simplify = F)
  }
  
  cnalist = Map(function(m, x, y) m[x, y],
                x = genelist,
                y = samples,
                MoreArgs = list(m = cna))
  
  Unlist(sapply(cnalist, .cna_signal, simplify = F), nested.names = T)
}


