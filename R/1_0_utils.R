#remove functions
#rm(list=lsf.str())

#creat dirtories if not exist
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

#time ouput as a string character 
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
#NA == NA -> FALSE (return FALSE if both are NAs. In this case, whenever seeing NA)
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
#' @param mats a list of matrices (or a single matrix)
#' @return dim for each matrix provided.
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

#' @export
split_matrix = function(m, by) {
  stopifnot(has_dim(m))
  stopifnot(is.character(by))
  stopifnot(all(by %in% colnames(m)))
  list(x = m[, by, drop = F], y = m[, !colnames(m) %in% by, drop = F])
}

#' @export
have_equal_nrows = function(m1, m2) {
  nrow(m1) == nrow(m2)
}

#' @export
have_equal_rownames = function(m1, m2) {
  all(rownames(m1) == rownames(m2))
}

#' @export
is_square = function(m) {
  nrow(m) == ncol(m)
}

#' @export
have_equal_dims = function(m1, m2) {
  identical(dim(m1), dim(m2))
}

#' @export
is_cor = function(m) {
  rg = range(m)
  if ((is_square(m)) & (rg[1] >= -1) & (rg[2] <= 1)) {
    dg = unique(diag(m))
    return(length(dg) == 1 & dg == 1)
  }
  FALSE
}

#' @export
is_symm = function(m) {
  (is_square(m)) && (sum(m == t(m)) == nrow(m)^2)
}



