#remove functions
#rm(list=lsf.str())

#filter gene_cell matrix by one variable in meta (eg.pct_MT, pct_RB, pct_HB)
#return a lost of 2 elements (mtx and meta)
# filter_cell_bymeta <- function(mtx, meta, var="pct_MT", max=20){
#   message("running 'filter_cell_bymeta' ... return a list(mtx, meta)")
#   #check if cells match among mtx and meta
#   n_cell_meta <- nrow(meta)
#   n_cell_mtx <- ncol(mtx)
#   name_cell_meta <- row.names(meta)
#   name_cell_mtx <- colnames(mtx)
#   logical_name <- name_cell_mtx %in% name_cell_meta
#   message("In metadata, the sum cell(s): ", n_cell_meta)
#   message("In gene-cell mtx, the sum cell(s): ", n_cell_mtx)
#   message("The overlapped cells between meta and gene-cell mtx: ", sum(logical_name))
#   if ( identical (name_cell_meta, name_cell_mtx) ) {
#     mtx <- mtx[,meta[[var]]<max] #in mtx, cell is at column
#     meta <- meta[meta[[var]]<max,] #in meta, cell is at row
#     n_removed <- n_cell_mtx - sum(meta[[var]]<max)
#     message(n_removed, " cells are removed from matrix and meta")
#   } else { stop("the overlapped cells between metadata and mtx dose NOT match")}
#   sketch_meta(meta); sketch_mtx(mtx); return(list(mtx=mtx, meta=meta))
# }

filter_cells_bymeta <- function(mtx,meta,meta.col=NULL,max=NULL,min=NULL,character=NULL,verbose=TRUE){
  if (verbose) message("running 'filter_cell_bymeta'... return a list(mtx, meta)")
  #check if cells match among mtx and meta
  n_cell_meta <- nrow(meta)
  n_cell_mtx <- ncol(mtx)
  name_cell_meta <- row.names(meta)
  name_cell_mtx <- colnames(mtx)
  logical_name <- name_cell_mtx %in% name_cell_meta
  if (verbose) message("In metadata, the sum cell(s): ", n_cell_meta)
  if (verbose) message("In gene-cell mtx, the sum cell(s): ", n_cell_mtx)
  if (verbose) message("The overlapped cells between meta and gene-cell mtx: ", sum(logical_name))
  if ( check_mtx_meta(mtx=mtx, meta=meta, verbose = F) ) {
    if (!is.null(max)) {
      mtx <- mtx[,meta[[meta.col]]<max] #in mtx, cell is at column
      meta <- meta[meta[[meta.col]]<max,] #in meta, cell is at row
      n_removed <- n_cell_mtx - sum(meta[[meta.col]]<max)
    } else if (!is.null(min)) {
      mtx <- mtx[,meta[[meta.col]]>min] #in mtx, cell is at column
      meta <- meta[meta[[meta.col]]>min,] #in meta, cell is at row
      n_removed <- n_cell_mtx - sum(meta[[meta.col]]>min)
    } else if (!is.null(character)) {
      mtx <- mtx[,meta[[meta.col]]==character] #in mtx, cell is at column
      meta <- meta[meta[[meta.col]]==character,] #in meta, cell is at row
      n_removed <- n_cell_mtx - sum(meta[[meta.col]]==character)
    } else {stop ("Need to provide values for arguments <max> or <min> or <character>") }
    message(n_removed, " cells are removed from matrix and meta")
  } else { stop ("Note!!! The overlapped cells between metadata and mtx dose NOT match")}
  return(list(mtx=mtx, meta=meta))
}

#geneset percentage
geneset_percent <- function(mtx.count,gene=NULL,pattern=NULL,verbose=TRUE){
  if (verbose) message("running 'geneset_percent'... return a named(cell ID) vector(percentage)")
  #if (format != "count") {warning("input mtx: gene-cell UMI count.")}
  #check arguments <gene> and <pattern> #stopifnot(!is.null(h) | !is.null(k)) 
  if (!is.null(gene) && !is.null(pattern)) {
    warning("Both <pattern> and <gene> provided. <pattern> is being ignored.")
    gene <- gene
  } else if (is.null(gene) && !is.null(pattern)) {
    gene <- grep(pattern = pattern, x = rownames(mtx.count), value = TRUE)
  } else if (!is.null(gene) && is.null(pattern)) {
    gene <- gene
  } else {
    stop("Neither <gene> or <pattern> is provided")
  }
  #check if the row names exits as subseting row names gives NAs and no error is generated
  logical.gene <- gene %in% rownames(mtx.count) 
  if(all(logical.gene)) {
    message("all gene(s) are found in the UMI matrix(gene:row)"); print(gene)
  } else if (!all(logical.gene)) {
    stop("None of the gene(s) are found in the UMI matrix(gene:row)") ; print(gene)
  } else { 
    warning("some gene(s) are missing")
    miss.gene <- row.names(mtx.count)[!logical.gene]
    print(miss.gene)
  }
  #calculate the percentage #
  percent.gene <- colSums(mtx.count[gene, ,drop = FALSE])*100/colSums(mtx.count)
  return(percent.gene)
}

###filter1 count matrix by cells (vetical)
###Minimal 1000 expressed genes (>0 counts) to filter cells
###Note: run before filter2_gene; all genes are required to calcuate umi count 
# non-parallel version
# filter1_cells <- function(mtx.count,countcutoff=0,mingene=1000,
#                           maxgene=NULL,mincount=NULL,maxcount=NULL,verbose=TRUE){
#   if (verbose) message("running 'filter1_cell'... return a cell-filtered mtx.count")
#   if (verbose) {message("Before cell filter:"); sketch_mtx(mtx.count)}
#   #default filter: mingene
#   genespercell <- apply(mtx.count,2,function(x)sum(x>countcutoff))
#   mtx.count <- mtx.count[,genespercell > mingene]
#   #optional filter: maxgene, mincount, maxcount
#   if(!is.null(maxgene)){
#     mtx.count <- mtx.count[,genespercell < maxgene]}
#   if(!is.null(mincount)){
#     countpercell < colSums(mtx.count)
#     mtx.count <- mtx.count[,countpercell > maxcount]}
#   if(!is.null(maxcount)){ 
#     countpercell < colSums(mtx.count)
#     mtx.count <- mtx.count[,countpercell < maxcount]}
#   if (verbose) {message("After cell filter:"); sketch_mtx(mtx.count)} 
#   return(mtx.count)
# }

filter1_cells <- function(mtx.count,countcutoff=0,mingene=1000,
                          maxgene=NULL,mincount=NULL,maxcount=NULL,verbose=TRUE,cells.cut=NULL){
  if (verbose) message("running 'filter1_cell'... return a cell-filtered mtx.count")
  if (verbose) {message("Before cell filter:"); sketch_mtx(mtx.count)}
  #default filter: mingene
  if(is.null(cells.cut)){
    genespercell <- apply(mtx.count,2,function(m) sum(m>countcutoff))
  } else {
    cut = dplyr::ntile(1:ncol(mtx.count), cells.cut) #1 1 1 2 2 2 3 3 3
    cut.uniq =  unique(cut) #1 2 3 
    lst = parallel::mclapply(cut.uniq, mc.cores=cells.cut, function(x){
      apply(mtx.count[,cut==x],2,function(m) sum(m>countcutoff)) }
    )
    genespercell = unlist(lst, use.names = F)
  }
  mtx.count <- mtx.count[,genespercell > mingene]
  #optional filter: maxgene, mincount, maxcount
  if(!is.null(maxgene)){
    mtx.count <- mtx.count[,genespercell < maxgene]}
  if(!is.null(mincount)){
    countpercell < colSums(mtx.count)
    mtx.count <- mtx.count[,countpercell > maxcount]}
  if(!is.null(maxcount)){ 
    countpercell < colSums(mtx.count)
    mtx.count <- mtx.count[,countpercell < maxcount]}
  if (verbose) {message("After cell filter:"); sketch_mtx(mtx.count)} 
  return(mtx.count)
}

#filter2 count matrix by genes (horizontal); 
#Note: cpm needs to be calculated without filtering genes in the function; 
#Note: cpm can not be calculated after runing this gene-filtering function because cpm by definition required a full gene list.
#Note: this is why this function returns a gene-filtered cpm matrix; 
filter2_genes <- function(mtx.count,countcutoff=5,countcells=20,minavlog2=4,verbose=TRUE,whitelist=NULL,genes.cut=NULL){
  if (verbose) message("running 'filter2_gene'... return a gene-filtered mtx.cpm")
  #filter1: countcutoff and countcells
  if(is.null(genes.cut)){
    cellspergene = apply(mtx.count,1,function(m) sum(m>=countcutoff))
  } else {
    cut = dplyr::ntile(1:nrow(mtx.count), genes.cut) #1 1 1 2 2 2 3 3 3
    cut.uniq =  unique(cut) #1 2 3
    lst = parallel::mclapply(cut.uniq, mc.cores=genes.cut, function(x){
      apply(mtx.count[cut==x,],1,function(m) sum(m>=countcutoff)) }
    )
    #cellspergene = do.call(c,lst)
    cellspergene = unlist(lst, use.names = F)
  }
  countcells.threshold = cellspergene >= countcells
  
  #filter2: minavlog in cpm
  mtx.cpm <- cpm(mtx.count=mtx.count, verbose=T)
  #Note: NaN can be introduced in cpm transformation if colsums = 0 (cells have 0 counts)
  cpm.gene.av <- Matrix::rowMeans(mtx.cpm) #apply(mtx.cpm,1,mean)
  minavlog2.threshold <- log2(cpm.gene.av+1) >= minavlog2
  
  #filter3: whitelist
  whitelist.logi <- rownames(mtx.cpm) %in% whitelist
  
  #combine filter 1+2+3
  mtx.cpm <- mtx.cpm[countcells.threshold|minavlog2.threshold|whitelist.logi,]
  
  if (verbose) {message("Before gene filter:"); sketch_mtx(mtx.count)}
  if (verbose) {message("After gene filter:"); sketch_mtx(mtx.cpm)}
  is_gene_duplicated(mtx.cpm)
  return(mtx.cpm)
}

