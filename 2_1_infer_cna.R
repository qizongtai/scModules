#' @title Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @description Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @param m a matrix of genes X cells containing scRNA-seq expression data. The matrix should NOT be row-centered.
#' @param ref.cells a list of two or more character vectors, each containg cell IDs of normal cell types (one cell type per list element). Since these cell types are presumed normal cells, their CNA values can be used to correct the remaining CNA values. Note that all cell IDs should correspond to column names in <m>. Default: NULL
#' @param window the size of the window to use for the rolling mean. Units are in number of genes. Default: 100
#' @param range values in <m> above and below this range will be set to values in range. Default: c(-3, 3)
#' @param n a numeric value indicating if only top genes should be used in the CNA calculation and how many. if n is a fraction between 0 and 1, the number of genes included will n * nrow(m). Else the number of genes included will be n. Default: NULL
#' @param noise a numeric value, which if given, increases the boundaries within which CNA values are considered 0. Increases by <noise> i.e. the bounds become Minimum - noise and Maximum + noise. Default: 0.1
#' @param center.method method by which to center the cells after calculating CNA values. One of 'median', 'mean', etc.... Default: 'median'
#' @param isLog boolean value indicating whether the input expression matrix <m> is in log2 form. Default: FALSE
#' @param verbose print progress messages. Default: FALSE
#' @return a matrix of genes X cells of inferred CNA values. Note that n = (window - 1)/2 genes will be lost from either extremity of the genome (ie. n genes lost at the start of chromosome 1 and at the end of chromosome Y, if the genome in question is H.sapiens.)
#' @details Correction with reference cells' <ref.cells> CNAs: the boundaries of reference CNA values are the boundaries for what should be considered a CNA of 0. Thus, if the boundary is -0.1 and 0.1, then a cell with CNA = -0.5 will be corrected to -0.4, a cell with CNA value of 1 will be corrected to 0.9 and a cell with CNA value of 0.05 will be corrected to 0.
#' @examples 
#'  m = infercna::useData()
#'  cna = infercna::infer_cna(m = m, ref.cells = ref.lst)
#' @rdname infer_cna
#' @export 
infer_cna = function(m,
                     ref.cells = NULL,
                     min.ref.cells = 10,
                     window = 100,
                     range = c(-3, 3),
                     topgenes.n = NULL,
                     ref.noise = 0.2,
                     ref.scale = NULL,
                     cna.noise = 0.15,
                     center.method = 'median',
                     isLog = FALSE,
                     verbose = TRUE) {

    #QC1: m should NOT be row(gene)-centered
    if (all(round(range(rowMeans(m)), 3) == 0)) {
        stop('Matrix is row-centered. Please provide non-centered data.')
    }
    
    #QC2: remove small ref.cells
    pass.logical = lengths(ref.cells) > min.ref.cells
    if (all(pass.logical)) {
        message("All cell types in ref.cells have more than ", min.ref.cells, " cells.")
        message("The orignal genes x cells matrix is kept."); sketch_mtx(m)
    } else {
        failed.ref.cells = ref.cells[!pass.logical]
        failed.ref.cells.names = names(failed.ref.cells)
        message("Note: ", failed.ref.cells.names, " in ref.cells has/have less than ", min.ref.cells, " cells.")
        message("Remove :", failed.ref.cells.names, " from the matrix.")
        message("The oringal matrix stats:"); sketch_mtx(m)
        m = m[ , !colnames(m) %in% unlist(failed.ref.cells, use.names = F)]
        message("After removing, the matrix stats:"); sketch_mtx(m)
    }
   
    # ===calcualte top expr genes by average cpm per gene
    if (isLog) {
        if (verbose) message('Converte log2cpm to cpm...')
        m = unlog2cpm(m)
    }
    if (!is.null(topgenes.n)) {
        if (verbose) message('Filter the expression matrix to include only top ', topgenes.n, 'genes...')
        m = m[exp_topgenes(m, ngenes = topgenes.n), ]
    }
    # ===center, order and clip expression range in log2 space 
    if (verbose) message('Converte cpm to log2cpm...')
    m = log2cpm(m)
    if (verbose) message('Perfore mean-centering of the genes...')
    m = norm_row(m, centerby = 'mean')
    if (verbose) message('Order the genes by their genomic position...')
    m = order_genes(m); sketch_mtx(m)
    if (verbose) message('Cap expression matrix values to between ', range[[1]], ' and ', range[[2]], '..')
    m = clip(m, range = range); sketch_mtx(m)
    # ===calculate CNA by moving average of cpm
    if (verbose) message('Converte log2cpm to cpm for CNA calculation...')
    m = unlog2cpm(m)
    if (verbose) message('Splite genes by chromosome...')
    ms = split_genes(m, by = 'chr')
    if (verbose) message('Calculate rolling means with a window size of ', window, ' genes...')
    cna = sapply(ms, function(m) {try(run_mean(m, k = window, verbose = verbose))}, simplify = F)
    cna = cna[sapply(cna, class) != 'try-error' | !sapply(class, isFALSE)]
    cna = Reduce(rbind, cna) #a = do.call (rbind, cna0)
    # ===convert cpm CNA to log2cpm CNA and center by median
    if (verbose) message('Converte cpm CNA to log2 space...')
    cna = log2cpm(cna)
    if (verbose) message('Performe ', center.method, '-centering of the cells...')
    cna = center_col(cna, by = center.method) # note: using median centering here
    # ===correct log2 CNA by ref log2 CNA
    if (!is.null(ref.cells)) {
        if (verbose) message('Correcte log2 CNA profiles using log2 CNA from ref.cells...')
        Args = c(list(cna = cna, ref.noise = ref.noise, ref.scale=ref.scale, isLog = TRUE), ref.cells) #Args needs to be a list
        cna = do.call(correct_ref, Args)}
    # ===clean noise: this usually coverts > 50% of the low values to 0
    if (!is.null(cna.noise)){
        if (verbose) message('Set log2 CNA values within ', cna.noise, ' and -', cna.noise, ' to 0 ...')
        cna[cna<cna.noise & cna>(-cna.noise)] = 0 ; sketch_mtx(cna)
    } 
    if (verbose) message('Done!')
    cna
}

# load("data/infercna/bt771.rda")
# load("data/infercna/mgh125.rda")
# load("data/infercna/refCells.rda")
# m <- useData()
# cna = infer_cna(m = m, ref.cells = refCells)
