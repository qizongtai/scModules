###======infercna utils======###

##########=====================================================================
#' @title Load Example Data
#' @description Load example dataset of human single-cell RNA seq data. Cells and genes are columns and rows, respectively, and all passed QC. Matrix is row-centered.
#' @param rowcenter center the expression matrix row-wise. Default: F
#' @return centered expression matrix of genes X cells 
#' @rdname useData
#' @export
useData = function(x = NULL, rowcenter = F) {
  if (is.null(x)) x = cbind(as.matrix(bt771), as.matrix(mgh125))
  else x = as.matrix(x)
  x = x[, !duplicated(colnames(x))]
  if (rowcenter) x = center_row(x, by = 'mean')
  x
}

##########======================================================================
clip <- function(mtx, range=c(-3, 3),verbose=TRUE) {
  message("running 'clip'... will return a clipped matrix (gene x cell)")
  mtx = as.matrix(mtx)
  mtx[mtx < range[[1]]] <- range[[1]]
  mtx[mtx > range[[2]]] <- range[[2]]
  return(mtx)
}

##########======================================================================
# #this function is replaced by top_expr_genes
# .mostExpressedGenes = function(m, ngenes) {
#     if (ngenes < 0) stop('<ngenes> cannot be negative.')
#     if (ngenes > nrow(m)) stop('<ngenes> cannot be larger than nrow(m).')
#     if (ngenes >= 0 & ngenes <= 1) ngenes = ngenes * nrow(m)
#     else if (ngenes > nrow(m)) ngenes = nrow(m)
#     rom = rowMeans(m)
#     names(sort(rowMeans(m), decreasing = T)[1:ngenes])
# }

# #to replace .mostExpressedGenes
exp_topgenes <- function(mtx, ngenes, verbose=TRUE) {
    if (verbose) message("running 'top_expr_genes'... will return a character vector of gene names")
    if (ngenes < 0) stop('<ngenes> cannot be negative.')
    else if (ngenes > nrow(mtx)) stop('<ngenes> cannot be larger than nrow(mtx).')
    else if (ngenes >= 0 & ngenes <= 1) ngenes = ngenes * nrow(m)
    return(names(sort(rowMeans(mtx), decreasing = T)[1:ngenes]))
}

##########======================================================================
.match_Genv_v1 <- function(genes) { match(Genv$gene, genes, nomatch = 0) }

#' @title Order Genes by their Genomic Positions
#' @description Order Genes by their Genomic Positions
#' @param x a matrix with gene names as rows or a character vector of gene names.
#' @return ordered matrix or character vector
#' @examples 
#' m = infercna::useData()
#' orderGenes(m) %>% rownames %>% head
#' orderGenes(rownames(m)) %>% head
#' @rdname orderGenes 
#' @export

order_genes <- function(x) {
  #return a ordered vector or matrix same as input format 
  if (is.null(dim(x))) return(x[.match_Genv_v1(genes = x)])
  x[.match_Genv_v1(rownames(x)), , drop = FALSE]
}

##########=======================================================================
#return position indices of Genv$gene for the matched genes
.match_Genv_v2 = function(genes) {match(genes, Genv$gene, nomatch = 0)}

#filter Genv based on matched genes and return a filtered Genv (list)
.glob_Genv = function(genes = NULL) {
  #return a filtered Genv (list)
  L = as.list(Genv)
  L = L[sapply(L, function(l) is.atomic(l) && length(l) > 1)] #filter1: class is atomic and length > 1
  if (!is.null(genes)) {L = sapply(L, `[`, .match_Genv_v2(genes), simplify = F)} #filter2: mathed genes
  return(L)
}

# old function; use if new one not work
# .order_chrom = function(chrom.names=NULL, return.value="chr") {
#     #retrun a sorted vector of character(chrom) or index of chrom.names
#     xnum = sort(as.numeric(chrom.names[chrom.names %in% as.character(1:100)]))
#     xchar = sort(chrom.names[!chrom.names %in% xnum]) #sort(chrom.names[!chrom.names %in% as.character(1:100)])
#     if (return.value=="chr") { return (c(xnum, xchar))
#     } else if (return.value=="index") { return (match(c(xnum, xchar), chrom.name)) }
# }

.order_chrom = function(chrom.names=NULL, return.value="chr", exclude.y=TRUE) {
  #retrun a sorted vector of character(chrom) or index of chrom.names
  if (isTRUE(exclude.y)) { chrom.ref = c(1:22, "X")
  } else { chrom.ref = c(1:22, "X", "Y") }
  if (return.value=="chr") { return ( chrom.ref[sort(match(chrom.names, chrom.ref))] )
  } else if (return.value=="index") { return ( na.omit(order(match(chrom.names, chrom.ref))) ) 
  } else ( stop( ' use a valid character for <return.value>: "chr" or "index" ') )
}

# old function; use if new one not work
# .order_arm = function(chrom.arms=NULL, return.value="chr", ) {
#   #retrun a sorted vector of char or index of chrom.names
#   #make arm reference
#   arm.ref = paste0(rep(c(1:22, "X", "Y"), each=2), c("p","q"))
#   if (return.value=="chr") { return (chrom.arms[order(match(chrom.arms, arm.ref))])
#   } else if (return.value=="index") { return (order(match(chrom.arms, arm.ref))) }
# }

.order_arm = function(chrom.arms=NULL, return.value="chr", exclude.y=TRUE) {
  #retrun a sorted vector of character(chrom) or index of chrom.arms
  if (isTRUE(exclude.y)) { arm.ref = paste0(rep(c(1:22, "X"), each=2), c("p","q"))
  } else { arm.ref = paste0(rep(c(1:22, "X", "Y"), each=2), c("p","q")) }
  if (return.value=="chr") { return ( arm.ref[sort(match(chrom.arms, arm.ref))] )
  } else if (return.value=="index") { return ( na.omit(order(match(chrom.arms, arm.ref))) ) 
  } else ( stop( ' use a valid character for <return.value>: "chr" or "index" ') )
}

#' @title Split Genes By Chromosome (Arm)
#' @description Split Genes By Chromosome (Arm). Input can be one of a matrix to be split into several matrices or a character vector of gene names to be split into several character vectors. 
#' @param x a matrix with gene names as rows or a character vector of gene names.
#' @param by a string; one of 'chr' or 'arm', determining how the matrix or character vector should be split. Default: 'chr'
#' @return if a matrix was provided, a list of matrices. If a character vector was provided, a list of character vectors. Both lists will be of length # of chromosome (arms) in the genome in question. Ie. length 24 for H.sapiens. 
#' @examples 
#' m = infercna::useData()
#' a = splitGenes(x = m)
#' b = splitGenes(x = rownames(m))
#' all(sapply(a, nrow) == lengths(b))
#' names(a)
#' names(splitGenes(x = m, by = 'arm'))
#' @rdname splitGenes
#' @export 
split_genes = function(x, by = 'chr', exclude.y=TRUE) {
  #message("running 'split_genes'... return a named (chrom) list of chr vectors (gene symbols) or matrix (gene x cell)")
  if (is.null(dim(x))) { genes = x } else { genes = rownames(x) }
  L = .glob_Genv(genes = genes) #filter Genv based on matched gene and return a list
  stopifnot(by %in% names(L)) #check <by='chr'> is in the name of the list
  splut = split(L$gene, L[[by]]) #return a list of char vectors (gene symbols)
  if (by=="chr") {
    splut = splut[.order_chrom(chrom.names = names(splut), exclude.y=exclude.y)] #order by the names(chromosomes) of splut
  } else if (by=="arm") {
    splut = splut[.order_arm(chrom.arms = names(splut), exclude.y=exclude.y)] #order by the names(arms) of splut
  }
  if (is.null(dim(x))) return(splut) #retun a list of char vectors
  sapply(splut, function(group) x[group, ], simplify = F) #return a list of matrixes
  # sapply(splut, "[", group, simplify = F)
}


##########=====================================================================
run_mean <- function(mtx,k=100,endrule='mean',align='center',verbose=TRUE) {
  if (verbose) message("running 'run_mean'... will return a data.frame or a named numeric vector")
  if (!is.null(dim(mtx))) {mtx = as.matrix(mtx)}
  if (is.null(dim(mtx))) return(FALSE)
  if (nrow(mtx) == 0) return(FALSE)
  if (nrow(mtx) < k) {
    k = nrow(mtx)
    if (verbose) message('Adjusting window to the max. number of genes in chromosome (', k, ')')
  }
  mout = caTools::runmean(mtx, k = k, endrule = endrule, align = align)
  if (!is.null(dim(mtx))) {
    colnames(mout) = colnames(mtx)
    rownames(mout) = rownames(mtx)
    mout = as.data.frame(mout)
  } else {
    names(mout) = names(mtx)
  }
  return(mout)
}

###########======================================================================
#return the min and max cpm for each gene of the refCells
#.refrange = function(cna, ...) # this version does calculations based on log2 format
.refrange = function (cna, isLog = FALSE, verbose = TRUE, ...) {
  #dots is refCell itself; it is a list of ref cells id
  #str(refCell)
  #List of 2
  #$ oligodendrocytes: chr [1:219] "MGH102-P1-A06" "MGH102-P1-A08" "MGH102-P1-D11" "MGH102-P1-E09" ...
  #$ macrophages     : chr [1:707] "MGH102-P1-B02" "MGH102-P1-B03" "MGH102-P1-B07" "MGH102-P1-B12" ...
  dots = list(...) #print(str(dots)); print(names(dots))
  if (isLog) {
    if (verbose) {message("converting log2-based CNA values to cpm-based for ref correction")}
    cna = unlog2cpm(cna) #old unlog2tpm function has a bulk=F argument
  } 
  v = sapply(dots, function(ref) rowMeans(cna[, ref, drop = F]), simplify = F)
  # v is a list (refcell) of named(gene) numeric vectors(gene mean)
  # str(v) #Note this results are based on log2 not unlog2 format
  # List of 2
  # $ fib: Named num [1:2949] 0.00874 0.00674 0.01002 0.00541 0.0103 ...
  # ..- attr(*, "names")= chr [1:2949] "HES4" "ISG15" "SDF4" "MXRA8" ...
  # $ end: Named num [1:2949] -0.0947 -0.0893 -0.0848 -0.0849 -0.0866 ...
  # ..- attr(*, "names")= chr [1:2949] "HES4" "ISG15" "SDF4" "MXRA8" ...
  L_min_max = list(Min = do.call(pmin, v), Max = do.call(pmax, v))
  #return(L_min_max) #this is unlog2cpm cna
  # str(L_min_max)
  # List of 2 #Note this results are based on log2 not unlog2 format
  # $ Min: Named num [1:2949] -0.0947 -0.0893 -0.0848 -0.0849 -0.0866 ...
  # ..- attr(*, "names")= chr [1:2949] "HES4" "ISG15" "SDF4" "MXRA8" ...
  # $ Max: Named num [1:2949] 0.00874 0.00674 0.01002 0.00541 0.0103 ...
  # ..- attr(*, "names")= chr [1:2949] "HES4" "ISG15" "SDF4" "MXRA8" ...
  L_min_max_log2cpm = sapply(L_min_max, log2cpm, simplify = F)
  return(L_min_max_log2cpm)
}
#a <- .refrange(cna = cna, isLog = FALSE, refCells[[1]], refCells[[2]])
#a <- .refrange(cna = cna, isLog = FALSE, a=refCells[[1]], b=refCells[[2]])

#return a numeric vector for the gene(row); input is uncorrected numeric vector across all cells for the gene(row)
.refcenter = function(gene, Min, Max, ref.noise = NULL, ref.scale = NULL) {
  #Expand reference boundaries by fixed noise factor or scaling percentage
  if (!is.null(ref.noise)) {
    Min = Min - ref.noise
    Max = Max + ref.noise
  } else if (!is.null(ref.scale)) {
    Min = Min - ref.scale*abs(Min)
    Max = Max + ref.scale*abs(Max)
  }
  above = gene > Max
  below = gene < Min
  normal = !above & !below
  gene[above] <- gene[above] - Max
  gene[below] <- gene[below] - Min
  gene[normal] <- 0
  return(gene)
}

#' @title Convert Relative CNA Values To Absolute
#' @description If the identities of normal cells are known, the expected CNA values of these cells should be 0. Thus, CNA values within the range of normal cell CNA values are corrected to 0; CNA values below the range have the minimum substracted; CNA values above the range have the maximum subtracted.
#' @param cna a matrix of rows (genes) by columns (cells) of (log2-based??) CNA values. 
#' @param ref.noise a numeric value, which if given, increases the boundaries within which CNA values are considered 0. Increases by <ref.noise> i.e. the bounds become Minimum - ref.noise and Maximum + ref.noise. Default: NULL
#' @param ... normal cell column IDs. Expects at least two character vectors. 
#' @rdname correct_ref
#' @export 
correct_ref = function(cna, ref.noise = NULL, ref.scale=NULL, isLog = FALSE, ...) {
  #dots are refCell=refCell, it is a list of ref cells id
  dots = list(...) #print(str(dots)); print(names(dots))
  genes = rownames(cna)
  Args = c(list(cna = cna, isLog = isLog), dots) #print(str(Args)); print(names(Args))
  rg = do.call(.refrange, Args)
  # str(rg)
  # List of 2
  # $ Min: Named num [1:2949] -0.493 -0.463 -0.434 -0.441 -0.45 ...
  # ..- attr(*, "names")= chr [1:2949] "HES4" "ISG15" "SDF4" "MXRA8" ...
  # $ Max: Named num [1:2949] 0.163 0.146 0.167 0.132 0.169 ...
  # ..- attr(*, "names")= chr [1:2949] "HES4" "ISG15" "SDF4" "MXRA8" ...
  n = nrow(cna)
  cna = t(sapply(1:n, function(i) {
          .refcenter(cna[i, ],
                     Min = rg$Min[i],
                     Max = rg$Max[i],
                     ref.noise = ref.noise,
                     ref.scale = ref.scale) }))
  rownames(cna) = genes
  cna
}
#Args = c(list(cna = cna, ref.noise = ref.noise, ref.scale = ref.scale, isLog = TRUE), refCells) #Args needs to be a list
#cna = do.call(refCorrect, Args)
