#identify the highly varible genes
hvg_row <- function(mtx.norm,method="std",top.n=3000,quantile=NULL,verbose=TRUE) {
  if (verbose) message("running 'hvg_generow'... return a matrix of hvg.norm.")
  if (method == "std"){ var = apply(mtx.norm, 1, sd)
  } else if (method == "mad"){ var = apply(mtx.norm, 1, mad)
  } else {stop("invalid argument <method>; choose std or mad")}
  if (!is.NULL(top.n)) { 
    mtx.norm.hvg = mtx.norm[order(var, decreasing = T),][1:top.n,]
  } else if (!is.NULL(quantile)) {
    mtx.norm.hvg = mtx.norm[var > quantile(var, prob = quantile),]
  }
  sketch_mtx(mtx.norm.hvg); return(mtx.norm.hvg)
}

#pca
run_pca = function(mtx.hvg,npcs=100,retx=TRUE,center=TRUE,scale=FALSE,verbose=TRUE) {
  if (verbose) message("running 'pca'... will return a pca obj = matrix: cells x PCs")
  #pca input: rows are objects(cells); columns are dimentions(genes); 
  pca.obj <- irlba::prcomp_irlba(t(mtx.hvg),n=npcs,retx=retx,center=center,scale.=scale)
  rownames(pca.obj$x) <- colnames(mtx.hvg) #assign cell id
  return(pca.obj)
}

#extract PC embedings
extract_pc <- function(pca.obj, npc=2, verbose=TRUE) {
  if (verbose) message("running 'extract_pc'... will return a matrix: cells x extracted.PCs ")
  #checking requested n.pc is larger than existing PCs
  if(npc > ncol(pca.obj$x)){stop("Error: the requested PCs(n.pc) is larger than the existing ones")}
  return(pca.obj$x[,c(1:npc)])
}

#extract top loadings genes


#umap
run_umap <- function(mtx.hvg,n.neighbors=100,min.dist=0.0001,spread=5,metric="correlation",verbose=TRUE){
  if (verbose) message("running 'umap'... will return an umap obj = matrix: cells x umap.dimensions")
  #umap input: rows are objects(cells); columns are dimentions(genes)
  umap.obj <- uwot::umap(t(as.matrix(mtx.hvg)), n_neighbors=n.neighbors, min_dist=min.dist, spread=spread, metric=metric)
  colnames(umap.obj) <- c("umap.x","umap.y")
  row.names(umap.obj) <- colnames(mtx.hvg)
  return(umap.obj)
}
