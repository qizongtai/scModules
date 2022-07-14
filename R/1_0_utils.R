#' Create directories if not exist
#'
#' @param dir a character vector of path names;
#' @param remove.dir a logical value. If FALSE, only the path not exist be created.
#'
#' @return a path
#' @export
#' @concept utils
#'
#' @examples
#' create_dirs("tmp_dir")
create_dirs <- function(dir, remove.dir=FALSE) {
  if (dir.exists(dir)) {
      message ("The directory ", dir, " already exists!")
      if (isFALSE(remove.dir)) { message ("...skip creating...") } 
      if (isTRUE(remove.dir)) {
        pause <- function () {
          var <- readline(prompt = "proceed removing the existing directory and all files? (y/n): ")
          if (toupper(var)=='Y') {
            message("entered 'Y/y'... proceed removing ...")
            unlink(dir, recursive = TRUE)
            res=dir.create(dir); 
            if(res) {message ("The directory ", dir, " is created!")
            } else {stop("error: The directory ", dir, " is NOT created!")}
          } else if (toupper(var)=='N') {
            message("entered 'N/n'... stop removing.")
            message ("...skip creating...") 
          } else { message("please enter valid letter 'y' or 'n'"); pause() }
        }
        pause()
      }
  } else {
    res=dir.create(dir); 
    if(res) {message ("The directory ", dir, " is created!")
    } else {stop ("error: The directory ", dir, " is created!")}
  }
}

#' A format time output as a string character 
#' 
#' if only ymd set TRUE, year-month-day
#' if only hms set TRUE, hour minute second
#' if all TRUE, year_month_day_hour minute second
#'
#' @param ymd a logical value
#' @param hms a logical value
#'
#' @return
#' a format time
#' @export
#' @concept utils
#'
#' @examples
#' time_str()
time_str = function(ymd=T, hms=T) {
  #stopifnot(ymd || hms) #output empty string instead of stop
  if (ymd==F && hms==F) return ("")
  if (ymd==T && hms==F) return (format(Sys.time(), "%Y-%m-%d"))
  if (ymd==F && hms==T) return (format(Sys.time(), "%H%M%S"))
  if (ymd==T && hms==T) return (format(Sys.time(), "%Y_%m_%d_%H%M%S"))
}

#NA == NA -> TRUE (return TRUE if both are NAs)
`%!=%` <- function(v1, v2) { (v1 != v2 | (is.na(v1) & !is.na(v2)) | (is.na(v2) & !is.na(v1))) & !(is.na(v1) & is.na(v2)) }
`%==%` <- function(v1, v2) { (v1 == v2 | (is.na(v1) & is.na(v2)) ) & !( (is.na(v1) & !is.na(v2)) | (is.na(v2) & !is.na(v1)) ) }
#NA == NA -> FALSE (return FALSE if both are NAs. In all cases, FALSE whenever seeing NA)
`%!=na%` <- function(v1, v2) { v1 != v2 | (is.na(v1) & !is.na(v2)) | (is.na(v2) & !is.na(v1)) | (is.na(v1) & is.na(v2)) }
`%==na%` <- function(v1, v2) { (v1 == v2) & !( (is.na(v1) & !is.na(v2)) | (is.na(v2) & !is.na(v1)) | (is.na(v1) & is.na(v2)) ) }

#shead a two dimential data (matrix or dataframe)
shead <- function(data,n=6){print(data[1:n,1:n])}
stail <- function(data,n=6){print(data[(nrow(data)-n):nrow(data),1:n])}

#check the duplicated genes at rows
is_gene_duplicated <- function(mtx, verbose=TRUE) {
  if (verbose) message("running 'is_gene_duplicated' ... return a logical value and message duplicated genes")
  dup.rownames <- row.names(duplicated(rownames(mtx)))
  n.dup.rownames <- length(dup.rownames)
  message(n.dup.rownames, " duplicated gene symbols")
  message(dup.rownames)
  if (n.dup.rownames==0) {return(FALSE)} else {return(TRUE)}
}

discrete16=c("#e5192c","#3a77b7","#3cac4c","#813c93",
             "#f36c24","#37b8c3","#a54922","#6b7627",
             "#28996b","#965b6a","#e9148f","#595b5e",
             "#76c3ad","#80d08a","#d29099","#f2e010")
#show_col(discrete16)

#library(RColorBrewer)
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#brewer.discrete74 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
brewer.discrete74= c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", 
                     "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", 
                     "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", 
                     "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", 
                     "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", 
                     "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", 
                     "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", 
                     "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", 
                     "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", 
                     "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", 
                     "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", 
                     "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", 
                     "#CCEBC5", "#FFED6F")

if(F){
  mytheme <- theme(plot.title = element_text(size = 14,color="black",hjust = 0.5),
                   axis.title = element_text(size = 14,color ="black"), 
                   axis.text = element_text(size= 14,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  )
}

###====== https://rdrr.io/github/ ======###
#' @export
# old version
# has_dim <- function(x) {
#   if (is.data.frame(x)) x = as.matrix(x)
#   if (class(mtx.count) == "dgTMatrix" | class(mtx.count) == "dgCMatrix") x = Matrix::as.matrix(x)
#   !is.null(attr(x, "dim"))
# }
#attr(df, "dim") = NULL; #so need convert df to matrix
#but dim(df) returns dimensions

#' @export
# new version
has_dim <- function(x) {!is.null(dim(x))}

#' @title <dim> for many matrices
#' @description Returns the result of dim for every matrix in a list
#'
#' @param mats a list of matrices (or a single matrix)
#'
#' @return dim for each matrix provided.
#' @concept utils
#' @rdname dims
#' @export 
dims <- function(mats) {
  # if mats is a single matrix:
  if (!is.null(dim(mats))) {
    return(dim(mats))
  }
  # if mats is a list of matrices:
  sapply(mats, dim, simplify = T)
}

#' @title <ncol> for many matrices
#' @description Returns the result of ncol for every matrix in a list
#' @param mats a list of matrices (or a single matrix)
#' @return ncol for each matrix provided.
#' @concept utils
#' @rdname ncols
#' @export 
ncols <- function(mats) {
  # if mats is a single matrix:
  if (!is.null(dim(mats))) {
    return(ncol(mats))
  }
  # if mats is a list of matrices:
  sapply(mats, ncol, simplify = T)
}


#' @title <nrow> for many matrices
#' @description Returns the result of nrow for every matrix in a list
#' @param mats a list of matrices (or a single matrix)
#' @return nrow for each matrix provided.
#' @concept utils
#' @rdname nrows
#' @export 
nrows <- function(mats) {
  # if mats is a single matrix:
  if (!is.null(dim(mats))) {
    return(nrow(mats))
  }
  # if mats is a list of matrices:
  sapply(mats, nrow, simplify = T)
}

#' split matrix
#' 
#' @param m a matrix
#' @param by a col names
#'
#' @export
#' @concept utils
split_matrix = function(m, by) {
  stopifnot(has_dim(m))
  stopifnot(is.character(by))
  stopifnot(all(by %in% colnames(m)))
  list(x = m[, by, drop = F], y = m[, !colnames(m) %in% by, drop = F])
}

#' @concept utils
#' @export
have_equal_nrows = function(m1, m2) {
  nrow(m1) == nrow(m2)
}

#' @concept utils
#' @export
have_equal_rownames = function(m1, m2) {
  all(rownames(m1) == rownames(m2))
}

#' @concept utils
#' @export
is_square = function(m) {
  nrow(m) == ncol(m)
}

#' @concept utils
#' @export
have_equal_dims = function(m1, m2) {
  identical(dim(m1), dim(m2))
}

#' @concept utils
#' @export
is_cor = function(m) {
  rg = range(m)
  if ((is_square(m)) & (rg[1] >= -1) & (rg[2] <= 1)) {
    dg = unique(diag(m))
    return(length(dg) == 1 & dg == 1)
  }
  FALSE
}

#' @concept utils
#' @export
is_symm = function(m) {
  (is_square(m)) && (sum(m == t(m)) == nrow(m)^2)
}



