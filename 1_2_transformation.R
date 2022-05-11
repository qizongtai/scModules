#remove functions
#rm(list=lsf.str())

#count to CPM
.cpm <- function(mtx.count,verbose=FALSE){
  if (verbose) message("running 'cpm'... return a cpm mtx")
  #if (has_dim(mtx.count)) mtx.count = as.matrix(mtx.count)
  count.sum <- Matrix::colSums(mtx.count) #apply(mtx.count,2,sum)
  #Note: NaN can be introduced in cpm transformation if colSums = 0 (cells have 0 counts).
  if (any(count.sum==0)) { message ("Note: some cell(s) have 0 counts") }
  mtx.cpm <- t(t(mtx.count)/count.sum)*1000000
  #mtx.cpm <- scale(mtx, center = F, scale = colSums(counts)/1e6L)
  sketch_mtx(mtx.cpm); return(mtx.cpm)
}

cpm <- function(mtx.count,cells.cut=NULL,parallel=FALSE,verbose=FALSE){
  if(is.null(cells.cut)){
    .cpm(mtx.count=mtx.count,verbose=verbose)
  } else {
    cut = ntile(1:ncol(mtx.count), cells.cut)
    cut.uniq =  unique(cut)
    if(isTRUE(parallel)) { lst = parallel::mclapply(cut.uniq, mc.cores=cells.cut, function(x){.cpm(mtx.count[,cut==x], verbose=FALSE)} ) 
    } else { lst = lapply(cut.uniq, function(x){.cpm(mtx.count[,cut==x], verbose=FALSE)}) }
    mtx.cpm = do.call(cbind, lst); rm(lst)
    sketch_mtx(mtx.cpm); return(mtx.cpm)
  }
}

#cpm to log2cpm; sc: scale=10; bulk: scale=1
.log2cpm <- function(mtx.cpm,scale=10,verbose=FALSE) {
  if (verbose) message("running 'log2cpm'... return a log2-transformed mtx")
  #if (has_dim(mtx.cpm)) mtx.cpm = as.matrix(mtx.cpm)
  mtx.log2cpm <- log2((mtx.cpm/scale)+1)
  sketch_mtx(mtx.log2cpm); return(mtx.log2cpm)
}

log2cpm <- function(mtx.cpm,scale=10,cells.cut=NULL,parallel=FALSE,verbose=FALSE){
  if(is.null(cells.cut)){
    .log2cpm(mtx.cpm=mtx.cpm,scale=scale,verbose=verbose)
  } else {
    cut = ntile(1:ncol(mtx.cpm), cells.cut)
    cut.uniq =  unique(cut)
    if(isTRUE(parallel)) { lst = parallel::mclapply(cut.uniq, mc.cores=cells.cut, function(x){.log2cpm(mtx.cpm[,cut==x], verbose=FALSE)} ) 
    } else { lst = lapply(cut.uniq, function(x){.log2cpm(mtx.cpm[,cut==x], verbose=FALSE)}) }
    mtx.log2cpm = do.call(cbind, lst); rm(lst)
    sketch_mtx(mtx.log2cpm); return(mtx.log2cpm)
  }
}

#unlog2cpm; sc: scale=10; bulk: scale=1
.unlog2cpm <- function(mtx.log2cpm,scale=10,verbose=FALSE) {
  if (verbose) message("running 'unlog2cpm'... return a unlog2 cpm mtx")
  #if (has_dim(mtx.log2cpm)) mtx.log2cpm = as.matrix(mtx.log2cpm)
  mtx.cpm <- scale*(2^(mtx.log2cpm)-1)
  sketch_mtx(mtx.cpm); return(mtx.cpm)
}

unlog2cpm <- function(mtx.log2cpm,scale=10,cells.cut=NULL,parallel=FALSE,verbose=FALSE){
  if(is.null(cells.cut)){
    .unlog2cpm(mtx.log2cpm=mtx.log2cpm,scale=scale,verbose=verbose)
  } else {
    cut = ntile(1:ncol(mtx.log2cpm), cells.cut)
    cut.uniq =  unique(cut)
    if(isTRUE(parallel)) { lst = parallel::mclapply(cut.uniq, mc.cores=cells.cut, function(x){.unlog2cpm(mtx.log2cpm[,cut==x], verbose=FALSE)} ) 
    } else { lst = lapply(cut.uniq, function(x){.unlog2cpm(mtx.log2cpm[,cut==x], verbose=FALSE)}) }
    mtx.unlog2cpm = do.call(cbind, lst); rm(lst)
    sketch_mtx(mtx.unlog2cpm); return(mtx.unlog2cpm)
  }
}

# another version of unlog2cpm (from https://rdrr.io/github)
# unlog2cpm = function(m, bulk = F) {
#   # wrapper around scalop::tpm since scalop::tpm is confusing.
#   # in that it does not generate tpm from counts, but rather removes log, scaling and pseudocount
#   if (has_dim(m)) m = as.matrix(m)
#   if (bulk) x = 1
#   else x = 10
#   (2^m) * x - 1 
#   #x * (2^(m) - 1)
# }

#normalization: center and scale
norm_row <- function(mtx,center=TRUE,centerby="mean",scale=FALSE,scaleby=NULL,verbose=TRUE){
  if (verbose) message("running 'norm_row'... return a centered and/or scaled mtx")
  if(center==TRUE){
    if (centerby=="mean") {
      if (scale==TRUE) {
        if (scaleby=="std") {
          mtx.norm = t(scale(x=t(mtx), center=T, scale=T))
        } else if (scaleby=="mad") {
          mtx.norm = t(scale(x=t(mtx), center=T, scale=apply(mtx,1,mad)))
        } else {message("unrecognized logical for argument <scaleby>; default 'NULL'; options are 'std' and 'mad")}
      } else if (scale==FALSE) {
        mtx.norm = t(scale(x=t(mtx), center=T, scale=F))
      } else {message("unrecognized logical for argument <scale>; default 'FALSE'; alternative 'TRUE'")}
    } else if (centerby=="median") {
      if (scale==TRUE) {
        if (scaleby=="std") {
          mtx.norm = t(scale(x=t(mtx), center=apply(mtx,1,median), scale=T))
        } else if (scaleby=="mad") {
          mtx.norm = t(scale(x=t(mtx), center=apply(mtx,1,median), scale=apply(mtx,1,mad)))
        } else {message("unrecognized logical for argument <scaleby>; default 'NULL'; options are 'std' and 'mad")}
      } else if (scale==FALSE) {
        mtx.norm = t(scale(x=t(mtx), center=apply(mtx,1,median), scale=F))
      } else {message("unrecognized logical for argument <scale>; default 'FALSE'; alternative 'TRUE'")}
    } else {message("unrecognized logical for argument <centerby>; default 'mean'; alternative 'median'")}
  } else if (center==FALSE) {
    mtx.norm = mtx
  } else {message("unrecognized char for argument <center>; default 'TURE'; alternative 'FALSE")}
  sketch_mtx(mtx.norm); return(mtx.norm)
}

center_col = function(mtx,by="mean") {
  mtx = as.matrix(mtx)
  if (by == "mean")  by = Matrix::colMeans(mtx, na.rm = T)
  else if (by == "median") by = matrixStats::colMedians(mtx, na.rm = T)
  else stopifnot(is.numeric(by) & length(by) == ncol(mtx))
  scale(mtx, center = by, scale = F)
}

center_row = function(mtx,by="mean") {
  #mtx = as.matrix(mtx)
  if (by == "mean")  by = Matrix::rowMeans(mtx, na.rm = T)
  else if (by == "median") by = matrixStats::rowMedians(mtx, na.rm = T)
  else stopifnot(is.numeric(by) & length(by) == nrow(mtx))
  t(scale(t(mtx), center = by, scale = F))
}
