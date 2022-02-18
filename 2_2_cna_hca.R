hca_cells = function(cna, 
                     groups = NULL, 
                     subset.genes = NULL, 
                     cor.method = 'pearson',
                     dist.method = 'euclidean', 
                     cluster.method = 'average', 
                     verbose = TRUE,
                     ...) {
  #Note: '...' argument is for hca done at the single sample
  if(verbose) {message("running 'hca_cells'... return a list of ordered mtx or ordered groups" )}
  cna.orig = cna
  if (!is.null(subset.genes)) {
    subset.genes = subset.genes[subset.genes %in% rownames(cna)]
    cna = cna[subset.genes, ]
    cat('Using', nrow(cna), 'genes for ordering...')
  }
  # if (!is.null(subset.genes)) {
  #     atts = sapply(subset.genes, .getAttributeName)
  #     Args = split(subset.genes, atts)
  #     cna = filterGenes(x = cna, value = Args, attribute = names(Args))
  # }
  
  if (!is.null(groups)) {
    res = hca_groups(m = cna, groups = groups, interorder = F, intraorder = T, 
                     cor.method = cor.method, dist.method = dist.method, cluster.method = cluster.method)
    ordered.id = colnames(res$ordered.mtx) 
    ordered.groups = res$ordered.groups; rm(res)
  } else {
    res = hca(x = center_row(cna), cor.method = cor.method, dist.method = dist.method, 
              cluster.method = cluster.method, ...)
    ordered.id = res$order
    ordered.groups = res$groups; rm(res)
  }
  
  cat('done')
  return(list(mtx=cna.orig[,ordered.id], groups=ordered.groups))
  
}

#' @title Grouped hca order and Ungroup 
#' @description Keep desired groups together but reorder members within each group. E.g. If cells from multiple samples are being reordered, samples can be kept together and cells reordered within. Option to also reorder the order of groups or samples, such that groups with similar average profiles are placed next to one another. 
#' @param m matrix to be reordered
#' @param groups groups within which reordering should take place.
#' @param interorder if TRUE, group order itself is reordered such that groups with similar profiles are placed near one another. E.g. Samples with similar average CNA profiles. Default: FALSE
#' @param cor.method desired correlation metric. Default: 'pearson'
#' @param dist.method desired distance metric to be used on top of correlation matrix. Default: 'euclidean'
#' @param cluster.method desired agglomeration method. Default: 'average'
#' @param Names return the vector of ordered IDs instead of the reordered matrix. Default: FALSE
#' @return reordered matrix (same dimensions as input) or a character vector of ordered column names if Names = T. 
#' @rdname hca_groups
#' @export 
hca_groups = function(m,
                      groups,
                      interorder = FALSE,
                      intraorder = TRUE,
                      cor.method = 'pearson',
                      dist.method = 'euclidean',
                      cluster.method = 'average',
                      return.m = TRUE,
                      return.id = TRUE,
                      verbose = T
                      ) {
  if(verbose) {message("running 'hca_groups'... return ordered mtx or ordered groups or a list of both" )}
  if (interorder) {
    groupAvgs = sapply( groups, function(x) rowMeans(m[, x, drop = F]) )
    groups = groups[hca(center_row(groupAvgs),
                        cor.method = cor.method,
                        dist.method = dist.method,
                        cluster.method = cluster.method)$order]
  }
  if (intraorder) {
    m.list = sapply( groups, function(x) m[, x], simplify = F )
    groups = sapply( m.list, function(m) hca(center_row(m),
                                             cor.method = cor.method,
                                             dist.method = dist.method,
                                             cluster.method = cluster.method)$order )
  }
  ord = as.character(unlist(groups, use.names = F))
  m = m[, ord]
  if (return.m && return.id) { return(list(ordered.mtx=m, ordered.groups=groups)) }
  if (return.m) { return(m) } 
  if (return.id) { return(groups) }
  
}

#note this function doesn't specify dis.method and cluster.method
#probably for line 24 groups = groups[hca_order(center_rwo(groupAvgs))], no parameters needed in there.
# hca_order = function(x, return.steps = F, cor.method = 'pearson', ...) {
#   if (return.steps) return(.hca(x, hclust.end = T, cor.method = cor.method, ...))
#   .hca(x, hclust.end = T, cor.method = cor.method, ...)$order
# }

# #this may be a updated version for hca_order
# hca_order = function(x, return.steps = F, cor.method = 'pearson', dist.method = NULL, cluster.method = NULL) {
#   if (return.steps) return(.hca(x, hclust.end = T, cor.method = cor.method, dist.method = dist.method, cluster.method = cluster.method))
#   .hca(x, hclust.end = T, cor.method = cor.method, dist.method = dist.method, cluster.method = cluster.method)$order
# }

# hca_reorder = function(x, col = T, row = T,cor.method = 'none', ...) {
#   if (col) x = x[, .hca(x = x, hclust.end = T, cor.method = cor.method, ...)$order]
#   if (row) x = x[.hca(x = t(x), hclust.end = T, cor.method = cor.method, ...)$order, ]
#   x
# }

# hca = function(x,
#                cor.method = 'pearson', #scalop::cor.methods
#                dist.method = 'euclidean', #scalop::dist.methods
#                cluster.method = 'average', #scalop::cluster.methods
#                max.dist = 1,
#                h = NULL,
#                k = NULL,
#                min.size = 5,
#                max.size = 0.5,
#                verbose = verbose,
#                ...) {
#   if(verbose){message("running 'hca', return a list of results from steps..." )}
#   .hca(x,
#        cor.method = cor.method,
#        dist.method = dist.method,
#        cluster.method = cluster.method,
#        max.dist = max.dist,
#        h = h,
#        k = k,
#        min.size = min.size,
#        max.size = max.size,
#        ...)
# }

###=============================================
load("/scratch/splab/zqi/htan/code/data/cor.methods.rda")
load("/scratch/splab/zqi/htan/code/data/cluster.methods.rda")
load("/scratch/splab/zqi/htan/code/data/dist.methods.rda")
#.hca = function(x,
hca = function(x,
               cor.method = 'pearson', #scalop::cor.methods
               dist.method = 'euclidean', #scalop::dist.methods
               cluster.method = 'average', #scalop::cluster.methods
               max.dist = 1,
               cutree.h = NULL,
               cutree.k = NULL,
               min.size = 5,
               max.size = 0.5, 
               cor.end = F,
               dist.end = F,
               hclust.end = F,
               verbose = T) {
  if(verbose){message("running 'hca'... return a list of results from steps" )}
  res = c()
  
  #if (class(x) == 'matrix') { #class(matrix) ouput: "matrix" "array" 
  if (any(class(x) %in% 'matrix')) {
    x = .hca_cor(x, method = cor.method)
    res = c(res, list(cr = x))
    if (cor.end) return(res)
    x = .hca_dist(x, method = dist.method, max.dist = max.dist)
    res = c(res, list(dist = x))
    if (dist.end) return(res)
  }
  
  if (class(x) == 'dist') {
    x = .hca_tree(x, method = cluster.method) #stats::hclust(d, method = method)
    res = c(res, list(tree = x, order = x$labels[x$order]))
    if (hclust.end) return(res)
  }
  
  if (class(x) == 'hclust') {
    if (is.null(cutree.h)) cutree.h = res$tree$height
    x = .hca_cutree(x,
                    k = cutree.k,
                    h = cutree.h,
                    min.size = min.size, 
                    max.size = max.size)
    res = c(res, list(groups = x))
  }
  
  res
}

.hca_cor = function(m, method) {
  method = match.arg(method, cor.methods)
  if (method == 'none') return(m)
  if (is_cor(m)) {
    warning('\nComputing correlations over a correlation matrix...\n',
            'Set `cor.method = "none"` to skip correlation step.')
  }
  stats::cor(m, method = method)
}

.hca_dist = function(m, method, max.dist = NULL) {
  method = match.arg(method, dist.methods) #scalop::
  if (!is_square(m)) m = t(m)
  if (!is.null(max.dist)) {
    if (is_square(m) && unique(diag(m)) != max.dist) {
      warning("<max.dist> = ", max.dist, " but <m> diagonal = ", unique(diag(m)), "...")
    }
    m = max.dist - m
  }
  if (method == 'none') return(stats::as.dist(m))
  stats::dist(m, method = method)
}

.hca_tree = function(d, method) {
  method = match.arg(method, cluster.methods) #scalop::
  stats::hclust(d, method = method)
}

.hca_cutree = function(tree, k, h, min.size, max.size) {
  groups = .hca_cutree_as_list(tree = tree, k = k, h = h)
  ncells = length(tree$labels)
  min.size = .cluster_size(min.size, ncells)
  max.size = .cluster_size(max.size, ncells)
  lens = lengths(groups)
  groups[lens >= min.size & lens <= max.size]
}

.hca_cutree_as_list = function(tree, k, h) {
  stopifnot(!is.null(k) || !is.null(h))
  groups = stats::cutree(tree = tree, k = k, h = h)
  if (!has_dim(groups)) return(split(names(groups), groups))
  
  .clusterNames = function(cutreeOutput) {
    # change df colnames (tree heights) to have 4 decimal places followed by "_"
    colnames(cutreeOutput) = paste0(round(as.numeric(colnames(cutreeOutput)), 4), "_")
    # new clusterNames
    names(unlist(apply(cutreeOutput, 2, function(col) 1:length(table(col)))))
  }
  
  clusterNames = .clusterNames(cutreeOutput = groups)
  labels = rownames(groups)
  groups = as.list(as.data.frame(groups))
  groups = sapply(groups, function(ID) split(labels, ID), simplify =F)
  groups = unlist(groups, recursive = F, use.names = F)
  groups = stats::setNames(groups, clusterNames)
  unique(groups)
}

.cluster_size = function(cluster.size, ncells) {
  stopifnot(cluster.size >= 0 && cluster.size <= ncells)
  if (cluster.size > 1) return(cluster.size)
  cluster.size * ncells
}
