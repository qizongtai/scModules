#v2: v1 + multiple k + multiple clustering methods
#KNN graph-based community clustering
knn_cluster <- function(mtx,transposed=FALSE,k=c(30,50),methods=c('walktrap', 'louvain'),all.methods=FALSE,verbose=TRUE) {
  if (verbose) message("running 'knn_cluster'... will return a data frame: cells x cluster.numbers")
  funCalls = list(
    #leading_eigen = quote(igraph::cluster_leading_eigen(g)$membership),
    #label_prop = quote(igraph::cluster_label_prop(g)$membership),
    #fast_greedy = quote(igraph::cluster_fast_greedy(g)$membership), #get an error; see Onenote(scRNA: cluster)
    infomap = quote(igraph::cluster_infomap(g)$membership),
    louvain = quote(igraph::cluster_louvain(g)$membership),
    walktrap = quote(igraph::cluster_walktrap(g)$membership)
  )
  if (all.methods) methods = names(funCalls)
  if (!all(methods %in% names(funCalls))) stop()
  funCalls = funCalls[names(funCalls) %in% methods]
  
  if (!transposed) {mtx <- t(mtx)}
  #make sure k (nearest neighbours/cells) is NOT bigger than number of possible neighbours/cells
  k = k[k < nrow(mtx)] #rows are cells if transponsed
  k = as.list(k)
  
  out = list()
  for (ki in k) {
    #find neighbours by KNN; input: rows=objects(cells) x columns=dimentions(genes)
    knn <- FNN::get.knn(mtx, k = ki)
    #convert KNN to a weighted edge list (data frame)
    knn.df = data.frame(from=rep(1:nrow(knn$nn.index),ki), to=as.vector(knn$nn.index),
                        weight=1/(1+as.vector(knn$nn.dist)))
    #build KNN graph and do community clustering
    g <- igraph::graph_from_data_frame(knn.df, directed = FALSE)
    g <- igraph::simplify(g) #remove self loops
    if (verbose) message("vertex/node number ", igraph::vcount(g), ";", "edge number ", igraph::ecount(g))
    dat = sapply(funCalls, eval, envir = environment(), simplify = F)
    dat = do.call(cbind.data.frame, dat)
    dat$k = rep(ki, nrow(dat))
    dat$id = rownames(mtx)
    out = c(out, list(dat))
  }
  out = suppressMessages(Reduce(dplyr::full_join, out))
  out = out[, ncol(out):1] #reverse the column order
  print( lapply(split(out, out$k), function(x){lapply(as.data.frame(x[,methods]),table)}) )
  return(out)
}

#v2: v1 + added transposed argument
#SNN graph-based community clustering
#use PCA to reduce dimensions(genes) if argument ndim < number of genes(dimensions)
#default: ndim=all_genes(nrow(mtx)/ncol(mtx)); so no dimension reduction
snn_cluster <- function(mtx,transposed=FALSE,k=c(30,50),methods=c('walktrap','louvain'),all.methods=FALSE,
                        ndim=ifelse(transposed==FALSE,nrow(mtx),ncol(mtx)),verbose=TRUE, ...) {
  if (verbose) message("running 'snn_cluster'... will return a data frame: cells x cluster.numbers")
  funCalls = list(
    #leading_eigen = quote(igraph::cluster_leading_eigen(g)$membership),
    #label_prop = quote(igraph::cluster_label_prop(g)$membership),
    #fast_greedy = quote(igraph::cluster_fast_greedy(g)$membership), #get an error; see Onenote(scRNA:cluster)
    infomap = quote(igraph::cluster_infomap(g)$membership),
    louvain = quote(igraph::cluster_louvain(g)$membership),
    walktrap = quote(igraph::cluster_walktrap(g)$membership)
  )
  if (all.methods) methods = names(funCalls)
  if (!all(methods %in% names(funCalls))) stop()
  funCalls = funCalls[names(funCalls) %in% methods]
  
  if (!transposed) {mtx <- t(mtx)}
  #make sure k (nearest neighbours/cells) is NOT bigger than number of possible neighbours/cells
  k = k[k < nrow(mtx)]
  k = as.list(k)
  #make sure that dimensions(genes) less or equal to those in <mtx>
  ndim = ndim[ndim <= ncol(mtx)]
  ndim = as.list(ndim)
  
  out = list()
  for (d in ndim) {
    for (ki in k) {
      # build SNN graph
      g = scran::buildSNNGraph(x = mtx, k = ki, d = d, transposed=TRUE, ...)
      if (verbose) message(paste("vertex/node number", igraph::vcount(g), ";", "edge number", igraph::ecount(g)))
      dat = sapply(funCalls, eval, envir = environment(), simplify = F)
      dat = do.call(cbind.data.frame, dat)
      dat$ndim = rep(d, nrow(dat))
      dat$k = rep(ki, nrow(dat))
      dat$id = rownames(mtx) #cell ids are on the row after transposed
      out = c(out, list(dat))
    }
  }
  out = suppressMessages(Reduce(dplyr::full_join, out))
  out = out[, ncol(out):1]
  print( lapply(split(out, out$k), function(x){lapply(as.data.frame(x[,methods]),table)}) )
  return(out)
}

wide_gcluster <- function(df.gcluster, prefix=NULL, verbose=TRUE) {
  if (verbose) message("running 'wide_cluster'... will return a data frame: cells x each cluster id")
  df.by.k <- split(df.gcluster, df.gcluster$k)
  df.list <- list()
  for (df in df.by.k) {
    #delete $ndim as current version of snn_cluster does not use it to reduce dimensions
    if (!is.null(df$ndim)) { df$ndim <- NULL }
    #move $id to row.names and delete $id
    row.names(df) <- df$id; df$id <- NULL
    #delete $k and add k info to other columns' names
    k <- unique(df$k); df$k <- NULL
    colnames(df) <- paste0(prefix, ".", "k", k, ".", colnames(df))
    df.list <- c(df.list, list(df))
  } 
  return(do.call(data.frame, df.list))
}



