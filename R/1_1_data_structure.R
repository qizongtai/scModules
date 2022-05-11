#create mtx from the cellranger output
read_mtx <- function(indir,is.gz=TRUE, remove=FALSE, verbose=TRUE) {
  message("running 'read_mtx'... return a UMI count matrix (gene x cell)")
  if(is.gz){ lapply(paste0(indir,"/",list.files(indir, pattern="\\.gz$")), GEOquery::gunzip, remove=remove) }
  mtx <- Matrix::readMM(paste0(indir,"/matrix.mtx"))
  b <- read.table(paste0(indir,"/barcodes.tsv"), sep='\t')
  g <- read.table(paste0(indir, "/features.tsv"), sep='\t')
  colnames(mtx) <- b[,1]
  row.names(mtx) <- g[,2] #col1:Ensumble ID; col2:gene symbol
  row.names(mtx) <- sub("GRCh38_", "", rownames(mtx))
  row.names(mtx) <- sub("hpv", "HPV", rownames(mtx))
  sketch_mtx(mtx); return(mtx)
  #check the duplicated gene symbols
  is_gene_duplicated(mtx)
  # #combine expression values for the duplicated rows
  # #method1:group_by and summarise
  # mtx %>% tibble::rownames_to_column() %>% dplyr::group_by(rowname) %>% 
  # dplyr::summarise(across(where(is.numeric), sum)) %>% dplyr::select(-"rowname")
  # #method2:aggregate()
}

#creat metadata
create_metadata <- function(mtx.count,verbose=TRUE){
  if (verbose) message("running 'create_metadata'... return a meta df (rownames = cell ID)")
  #if (format != "count") {warning("input mtx: gene-cell UMI count.")}
  n_counts <- colSums(mtx.count)
  n_features <- colSums(mtx.count>0)
  meta <- data.frame(n_counts = n_counts, n_features = n_features, row.names = colnames(mtx.count))
  sketch_meta(meta); return(meta)
}

#add metadata (addmeta is a named vector (name = cell barcode))
add_metadata <- function(meta, addmeta, addname=colnames(as.data.frame(addmeta)),prefix=NULL,verbose=TRUE){
  if (verbose) message("running 'add_metadata'... return a new meta df with added cell info.")
  #check if addname exits in the meta
  if (!is.null(prefix)) {addname <- paste0(prefix, ".", addname)}
  if (any(addname %in% colnames(meta))) {stop("The addname(s) already exit in the metadata")}
  addmeta <- as.data.frame(addmeta) #transform to df if vector or matrix
  #addmeta <- as.data.frame(lapply(addmeta, as.factor)) #transform class
  #check if addname matches addmeta
  if (ncol(addmeta) != length(addname)) {stop("The number of addname(s) do not match with features in the added metadata")}
  #check the rownames(cells) of the addmeta with meta
  colnames(addmeta) <- addname
  n_addmeta <- nrow(addmeta)
  n_meta <- nrow(meta)
  logical.cells <- row.names(addmeta) %in% row.names(meta)
  overlap <- row.names(addmeta)[logical.cells]
  unfound <- row.names(addmeta)[!logical.cells]
  n_overlap <- length(overlap)
  n_unfound <- length(unfound)
  if (verbose) message("In the original metadata, the sum cells (barcodes): ", n_meta)
  if (verbose) message("In the added metadata, the sum cells (barcodes): ", n_addmeta)
  if (verbose) message("For the added metadata, the overlaped cell(s): ", n_overlap, "; the unfound cell(s): ", n_unfound)
  #merge meta with addmeta by overlapped row.names
  if (all(logical.cells) && n_addmeta == n_meta) {
    message("All cells in metadata are remained")
    merged <- merge(x=meta, y=addmeta, by="row.names")
    #move column $Row.names to row.names
    row.names(merged) <- merged$Row.names; merged$Row.names <- NULL
    sketch_meta(merged); return(merged)
  } else {
    message("Only the overlapped cells in metadata are remained.")
    pause <- function () {
      var <- readline(prompt = "proceed adding to the metadata? (y/n): ")
      if (toupper(var)=='Y') {
        message("entered 'Y/y'... proceed adding ...")
        merged <- merge(x=meta, y=addmeta, by="row.names")
        #move column $Row.names to row.names
        row.names(merged) <- merged$Row.names; merged$Row.names <- NULL
        sketch_meta(merged); return(merged)
      } else if (toupper(var)=='N') {
        message("entered 'N/n'... stop adding.")
        sketch_meta(meta); return(meta)
      } else { message("please enter valid letter 'y' or 'n'"); pause() }
    }
    pause()
  }
}
# #check the format (named vector or df) for addmeta
# if (is.vector(addmeta) || !is.null(names(addmeta))) {addmeta <- as.data.frame(addmeta)
# } else if (is.data.frame(addmeta)) { addmeta <- addmeta
# } else {stop("The addmeta is either named vector or data frame")}

#get cell and gene number
sketch_mtx <- function(mtx) {
  if (!has_dim(mtx)) { 
    len <- length(mtx) 
    message("mtx does not have dimension and its length is ", len, ".")
  } else {
  genes <- dim(mtx)[1]
  cells <- dim(mtx)[2]
  message(genes, " genes", "; ", cells, " cells")
  }
}

#get cell and feature information
sketch_meta <- function(meta) {
  if (has_dim(meta)) mtx.cpm = as.data.frame(meta)
  cells <- dim(meta)[1]
  features <- dim(meta)[2]
  #print(paste(cells, "cells", ";", features, "features"))
  #print(paste(colnames(meta), collapse = "; "))
  message(cells, " cells", "; ", features, " features")
  message(paste(colnames(meta), collapse = "; "))
}

#check the order of metadata and mtx
check_mtx_meta <- function(mtx, meta, verbose=TRUE) {
  if (verbose) {message("running 'check_mtx_meta'... return a logical value")}
  if (verbose) {sketch_mtx(mtx); sketch_meta(meta)}
  if ( identical(colnames(mtx), rownames(meta)) ) {
    message("Cell IDs are 1)matched and 2)in order in both mtx and meta"); return(TRUE)
  } else if ( all(colnames(mtx) %in% rownames(meta)) ) {
    message("Cell IDs are 1)matched but NOT 2)in order! To order: meta = meta[colnames(mtx),]"); return(FALSE)
  } else {message("Cell IDs are NOT matched!"); return(FALSE)}
}

#check the order of mtx1 and mtx2
check_mtx1_mtx2 <- function(mtx1, mtx2, verbose=TRUE) {
  if (verbose) {message("running 'check_mtx1_mtx2'... return a logical value")}
  if (verbose) {sketch_mtx(mtx1); sketch_mtx(mtx2)}
  if ( identical(colnames(mtx1), colnames(mtx2)) ) {
    message("Cell IDs are 1)matched and 2)in order in both mtx and meta"); return(TRUE)
  } else if ( all(colnames(mtx1) %in% colnames(mtx2)) ) {
    message("Cell IDs are 1)matched but NOT 2)in order! To order: mtx1 = mtx1[,colnames(mtx2)]"); return(FALSE)
  } else {message("Cell IDs are NOT matched!"); return(FALSE)}
}

#check the order of mtx1 and mtx2
checkgenes_mtx1_mtx2 <- function(mtx1, mtx2, verbose=TRUE) {
  if (verbose) {message("running 'check_mtx1_mtx2'... return a logical value")}
  if (verbose) {sketch_mtx(mtx1); sketch_mtx(mtx2)}
  if ( identical(rownames(mtx1), rownames(mtx2)) ) {
    message("Gene IDs are 1)matched and 2)in order in both mtx and meta"); return(TRUE)
  } else if ( all(rownames(mtx1) %in% rownames(mtx2)) ) {
    message("Gene IDs are 1)matched but NOT 2)in order! To order: mtx1 = mtx1[rownames(mtx2),]"); return(FALSE)
  } else {message("Gene IDs are NOT matched!"); return(FALSE)}
}

#check the order of meta1 and meta2
check_meta1_meta2 <- function(meta1, meta2, verbose=TRUE) {
  if (verbose) {message("running 'check_meta1_meta2'... return a logical value")}
  if (verbose) {sketch_meta(meta1); sketch_meta(meta2)}
  if ( identical(rownames(meta1), rownames(meta2)) ) {
    message("Cell IDs are 1)matched and 2)in order in both mtx and meta"); return(TRUE)
  } else if ( all(rownames(meta1) %in% rownames(meta2)) ) {
    message("Cell IDs are 2)matched but NOT 2)in order! To order: meta1 = meta1[rownames(meta2),]"); return(FALSE)
  } else {message("Cell IDs are NOT matched!"); return(FALSE)}
}

#sub-function
.pause <- function (m, format,logical.cells) {
  var <- readline(prompt = paste("proceed matching and return", format, "? (y/n): ") )
  if (toupper(var)=='Y') {
    message("entered 'Y/y'... proceed matching ...")
    if (format=="mtx") { 
      mtx.match <- m[,logical.cells]
      message("The final mtx:"); sketch_mtx(mtx.match); return(mtx.match)
    } else if (format=="meta") {
      meta.match <- m[logical.cells,]
      message("The final metadata:"); sketch_meta(meta.match); return(meta.match)
    } else { stop("please enter a valid character for <format>: mtx or meta") }
  } else if (toupper(var)=='N') {
    stop("entered 'N/n'... stop matching.")
  } else { message("please enter valid letter 'y' or 'n'"); .pause(m, format, logical.cells) }
}

#match the metadata and mtx
match_mtx_meta <- function(mtx, meta, return="mtx", force=FALSE, verbose=TRUE) {
  if (verbose) {message("running 'match_mtx_bymeta'... return a matched mtx or metadata (the same format as input)")}
  if (verbose) {message("In the mtx:"); sketch_mtx(mtx); message("In the metadata:"); sketch_meta(meta)}
  #check if mtx and meta are already matched using check_mtx_meta
  if (check_mtx_meta(mtx=mtx, meta=meta, verbose=FALSE)) {stop("###===The mtx and meta are matched already!===###")}
  #check the rownames(cells) of the addmeta with meta
  if (return=="mtx") {
    logical.cells <- colnames(mtx) %in% rownames(meta)
    overlap <- colnames(mtx)[logical.cells]
    unfound <- colnames(mtx)[!logical.cells]
    n_overlap <- length(overlap)
    n_unfound <- length(unfound)
    if (verbose) message("In the mtx: the overlaped cell(s): ", n_overlap, "; the unfound cell(s): ", n_unfound)
    if (force) {
      mtx.match <- mtx[,logical.cells]
      message("The final mtx:"); sketch_mtx(mtx.match); return(mtx.match)
    } else {.pause(m=mtx, format=return, logical.cells=logical.cells) }
  } else if (return=="meta") {
    logical.cells <- rownames(meta) %in% colnames(mtx)
    overlap <- rownames(meta)[logical.cells]
    unfound <- rownames(meta)[!logical.cells]
    n_overlap <- length(overlap)
    n_unfound <- length(unfound)
    if (verbose) message("In the metadata: the overlaped cell(s): ", n_overlap, "; the unfound cell(s): ", n_unfound)
    if (force) {
      meta.match <- meta[logical.cells,]
      message("The final metadata:"); sketch_meta(meta.match); return(meta.match)
    } else { .pause(m=meta, format=return, logical.cells=logical.cells) }
  } else { stop("please enter a valid character for <return>: mtx or meta") }
}

# old version
# match_mtx_meta <- function(mtx, meta, return="mtx", verbose=TRUE) {
#   if (verbose) {message("running 'match_mtx_bymeta'... return a matched mtx or metadata (the same format as input)")}
#   if (verbose) {message("In the mtx:"); sketch_mtx(mtx); message("In the metadata:"); sketch_meta(meta)}
#   #check the rownames(cells) of the addmeta with meta
#   logical_cells <- colnames(mtx) %in% rownames(meta)
#   overlap <- colnames(mtx)[logical_cells]
#   unfound <- colnames(mtx)[!logical_cells]
#   n_overlap <- length(overlap)
#   n_unfound <- length(unfound)
#   if (verbose) message("The overlaped cell(s): ", n_overlap, "; the unfound cell(s): ", n_unfound)
#   pause <- function () {
#     var <- readline(prompt = paste("proceed matching and return", return, "? (y/n): ") )
#     if (toupper(var)=='Y') {
#       message("entered 'Y/y'... proceed matching ...")
#       if (return=="mtx") { 
#         mtx.match <- mtx[,logical_cells]
#         message("The final mtx:"); sketch_mtx(mtx.match); return(mtx.match)
#       } else if (return=="meta") {
#         meta.match <- meta[logical_cells,]
#         message("The final metadata:"); sketch_meta(meta.match); return(meta.match)
#       } else { stop("please enter a valid character for <return>: mtx or meta") }
#     } else if (toupper(var)=='N') {
#       stop("entered 'N/n'... stop matching.")
#     } else { message("please enter valid letter 'y' or 'n'"); pause() }
#   }
#   pause()
# }

#match mtx1 by mtx2
match_mtx1_bymtx2 <- function(mtx1, mtx2, force=FALSE, verbose=TRUE) {
  if (verbose) {message("running 'match_mtx1_mtx2'... return a matched mtx1 by mtx2 (the same format mtx)")}
  if (verbose) {message("In the mtx1:");sketch_mtx(mtx1); message("In the mtx2:");sketch_mtx(mtx2)}
  #check if mtx1 and mtx2 are already matched using check_mtx1_mtx2
  if (check_mtx1_mtx2(mtx1=mtx1, mtx2=mtx2, verbose=FALSE)) {stop("###===The mtx and meta are matched already!===###")}
  #check the rownames(cells) of the addmeta with meta
  logical.cells <- colnames(mtx1) %in% colnames(mtx2)
  overlap <- colnames(mtx1)[logical.cells]
  unfound <- colnames(mtx1)[!logical.cells]
  n_overlap <- length(overlap)
  n_unfound <- length(unfound)
  if (verbose) message("The overlaped cell(s): ", n_overlap, "; the unfound cell(s): ", n_unfound)
  pause <- function () {
    var <- readline(prompt = "proceed matching mtx1 by mtx2? (y/n): ")
    if (toupper(var)=='Y') {
      message("entered 'Y/y'... proceed matching ...")
      mtx1.match <- mtx1[,logical.cells]
      sketch_mtx(mtx1.match); return(mtx1.match)
    } else if (toupper(var)=='N') {
      stop("entered 'N/n'... stop matching.")
    } else { message("please enter valid letter 'y' or 'n'"); pause() }
  }
  if(force) {
    mtx1.match <- mtx1[,logical.cells]
    sketch_mtx(mtx1.match); return(mtx1.match)
  } else { pause() }
}

#match meta1 by meta2
match_meta1_bymeta2 <- function(meta1, meta2, force=FALSE, verbose=TRUE) {
  if (verbose) {message("running 'match_meta1_meta2'... return a matched meta1 by meta2 (the same format as meta1)")}
  if (verbose) {sketch_meta(meta1); sketch_meta(meta2)}
  #check if mtx1 and mtx2 are already matched using check_mtx1_mtx2
  if (check_meta1_meta2(meta1=meta1, meta2=meta2, verbose=FALSE)) {stop("###===The mtx and meta are matched already!===###")}
  #check the rownames(cells) of the addmeta with meta
  logical.cells <- colnames(meta1) %in% colnames(meta2)
  overlap <- colnames(meta1)[logical.cells]
  unfound <- colnames(meta1)[!logical.cells]
  n_overlap <- length(overlap)
  n_unfound <- length(unfound)
  if (verbose) message("The overlaped cell(s): ", n_overlap, "; the unfound cell(s): ", n_unfound)
  pause <- function () {
    var <- readline(prompt = "proceed matching meta1 by meta2? (y/n): ")
    if (toupper(var)=='Y') {
      message("entered 'Y/y'... proceed matching ...")
      meta1.match <- meta1[logical.cells,]
      sketch_meta(meta1.match); return(meta1.match)
    } else if (toupper(var)=='N') {
      stop("entered 'N/n'... stop matching.")
    } else { message("please enter valid letter 'y' or 'n'"); pause() }
  }
  if(force) {
    meta1.match <- meta1[logical.cells,]
    sketch_meta(meta1.match); return(meta1.match)
  } else { pause() }
}

