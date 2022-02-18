.getAttributeName = function(val) {
  if (val %in% levels(Genv$chr)) return('chr')
  if (val %in% levels(Genv$arm)) return('arm')
  if (val %in% Genv$symbol) return('symbol')
  else stop('Value ', val, ' not found in genome attributes.')
}

.axisSpacer = function(breaks, labels, limits, levels = NULL) {
  if (!is.null(labels) & !is.null(levels)) {
    breaks = levels %in% labels
    labels = levels[breaks]
  }
  if (is.null(breaks)) breaks = seq(limits[[1]], limits[[2]], limits[[2]])
  if (is.null(labels)) labels = breaks
  list(breaks = breaks, labels = labels)
}

#<group> a list of cell ids
.cellBreaks = function(groups, halfway = F, hide = NULL) {
  #message(running 'cellBreaks'...return a named(group) num(position) vector)
  n = length(unlist(groups)) # total cell number
  chrsum = cumsum(lengths(groups)) # cell number in each group
  Breaks = chrsum/max(chrsum) * n # cell percentage length in each group
  # halfway useful for placing x axis chromosome text labels
  # halfway = F for placing x axis chromosome lines
  if (halfway) {
    b = Breaks
    Breaks = c(b[1]/2, b[2:length(b)] - diff(b)/2)
  }
  if (!is.null(hide)) {
    names(Breaks)[which(names(Breaks) %in% hide)] = rep('', length(hide))
  }
  Breaks
}

.chromosomeBreaks = function(genes = NULL, halfway = F, hide = NULL) {
  #message("running 'chromosomeBreaks'...return a named(chromosome) num(position) vector")
  n = length(genes) # total gene number
  chrsum = cumsum(lengths(split_genes(genes, by = 'chr'))) #gene number in each chr
  Breaks = (chrsum/max(chrsum)) * n # gene percentage length in each chr
  # n is the original gene number; genes that don't match Genv are removed from the output of split_genes  
  # halfway useful for placing x axis chromosome text labels
  # halfway = F for placing x axis chromosome lines
  if (halfway) {
    b = Breaks
    Breaks = c(b[1]/2, b[2:length(b)] - diff(b)/2)
  }
  if (!is.null(hide)) {
    names(Breaks)[which(names(Breaks) %in% hide)] = rep('', length(hide))
  }
  Breaks
}
# dat = reshape2::melt(as.matrix(cna2))
# colnames(dat) = c('Gene', 'Cell', 'CNA')
# genes=levels(dat$Gene)
# str(levels(dat$Gene))
# split_genes(genes, by = 'chr')
# xLineBreaks = chromosomeBreaks(levels(dat$Gene))

# ord = hca_order(x = center_row(cna2[,1:100]))
# ord = hca_order(x = center_row(cna2[,1:500]))
# #ord = colnames(hca_order(x = center_row(cna2[,1:100])))
# 
# head(ord)
# shead(cna2)

# if (group.ord) {
#   # breaks for horizontal lines denoting order.cells
#   yLineBreaks = .cellBreaks(order.cells)
#   # and breaks for their labels
#   yTextBreaks = .cellBreaks(order.cells, halfway = TRUE, hide = x.hide)
# }
# 
# # breaks for vertical lines denoting chromosomes
# xLineBreaks = .chromosomeBreaks(levels(dat$Gene))
# # and breaks for their labels
# xTextBreaks = .chromosomeBreaks(levels(dat$Gene), halfway = TRUE, hide = x.hide)
# yLineBreaks = ggplot2::waiver()
# yTextBreaks = ggplot2::waiver()
# 
# legend = .axisSpacer(breaks = legend.breaks, labels = legend.labels, limits = limits)
# legend.breaks = legend$breaks
# legend.labels = legend$labels

.cna_heatmap = function(...) {
  list2env(list(...), envir = environment())
  
  #setting up geomtile
  geomtile = ggplot2::geom_tile(color = tile.col, size = tile.size)
  if (any(sapply(list(tile.size, tile.col), is.null))) { geomtile = ggplot2::geom_tile() }
  
  #major plot
  G = ggplot2::ggplot(dat, ggplot2::aes(x = as.numeric(Gene),
                                        y = as.numeric(Cell),
                                        fill = CNA))  +
    geomtile +
    # ggplot2::scale_fill_gradientn(colors = cols,
    #                               limits = limits,
    #                               oob = scales::squish,
    #                               breaks = legend.breaks,
    #                               labels = legend.breaks,
    #                               name = legend.title,
    #                               guide = ggplot2::guide_colorbar(frame.colour = 'black',
    #                                                               ticks = F,
    #                                                               barheight = grid::unit(legend.height, "cm"),
    #                                                               barwidth = grid::unit(legend.width, "cm"),
    #                                                               title.position = legend.title.position)) +
    ggplot2::scale_fill_gradient2(low= "steelblue", 
                                  mid = "white", 
                                  high= "darkred", 
                                  midpoint= 0,
                                  limits = limits,
                                  oob = scales::squish,
                                  breaks = legend.breaks,
                                  labels = legend.breaks,
                                  name = legend.title,
                                  guide = ggplot2::guide_colorbar(frame.colour = 'black',
                                                                  ticks = F,
                                                                  barheight = grid::unit(legend.height, "cm"),
                                                                  barwidth = grid::unit(legend.width, "cm"),
                                                                  title.position = legend.title.position)) +
    ggplot2::labs(x = x.name,
                  y = y.name,
                  title = title,
                  subtitle = subtitle,
                  caption = caption) +
    hrbrthemes::theme_ipsum(base_size = base.size,
                            axis_title_size = axis.title.size,
                            axis_text_size = axis.text.size) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = xTextBreaks, labels = names(xTextBreaks)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = yTextBreaks, labels = names(yTextBreaks)) +
    #y.angle +
    #x.angle +
    ggplot2::geom_vline(xintercept = xLineBreaks, size = 0.08, linetype = 2) +
    ggplot2::geom_hline(yintercept = yLineBreaks, size = 0.08, linetype = 1) +
    ggplot2::theme(aspect.ratio = ratio,
                   title = ggplot2::element_text(colour = base.col),
                   #title = ggplot2::element_text(colour = "red"),
                   rect = ggplot2::element_rect(colour = base.col),
                   #rect = ggplot2::element_rect(colour = "red"),
                   #line = ggplot2::element_line(colour = base.col),
                   line = ggplot2::element_line(colour = "red"),
                   panel.border = ggplot2::element_rect(colour = 'black', fill = NA),
                   #panel.border = ggplot2::element_rect(colour = 'red', fill = NA),
                   legend.text = ggplot2::element_text(size = ggplot2::rel(legend.rel)), #0.9 *
                   legend.title = ggplot2::element_text(size = ggplot2::rel(legend.title.rel), angle = legend.title.angle),
                   plot.margin = grid::unit(c(1, 0, 0.3, 0.3), "cm"),
                   legend.margin = ggplot2::margin(0,0.3,0,-0.3,'cm'),
                   legend.justification = legend.justification)
  G
}


#heatCols <- colorRampPalette(c("blue", "white", "red"))(100)

#' @title Plot a Heatmap of CNA Values
#' @description Plot a heatmap of CNA values.
#' @param cna matrix of CNA values (genes by cells).
#' @param limits colour range. Cells >= upper limit will be the same colour, and likewise for cells <= lower limit. Default: c(-1, 1)
#' @param ratio aspect ratio of the panel. Default: NULL
#' @param cols character vector of colours to use. Default: heatCols
#' @param order.cells a boolean or a list of groups if ordering should be performed within each only; this is useful if plotting multiple samples in the same panel. Cells are ordered by hierarchical clusrering. Default: F
#' @param x.name x axis label. Default: 'Chromosome'
#' @param y.name y axis label. Default: 'Cell'
#' @param angle angle of axes tick labels. x.angle and y.angle inherit from angle. Default: NULL
#' @param x.angle angle of x axis tick labels. If left, will inherit from angle. Default: NULL
#' @param y.angle angle of y axis tick labels. If left, will inherit from angle. Default: 0
#' @param axis.rel relative size of axes labels. Default: 1
#' @param base.size base text size. Default: 12
#' @param axis.title.size axes titles text size. Default: 12
#' @param axis.text.size axes labels text size. Default: 11
#' @param base.col base text colour. Default: '#073642'
#' @param title title of plot. Default: NULL
#' @param subtitle subtitle of plot. Default: NULL
#' @param caption caption of plot. Default: NULL
#' @param text.size base text size. Default: 12
#' @param x.hide x axis labels to hide from plot. Default: c("13", "18", "21", "Y")
#' @param y.hide y axis labels to hide from plot. Default: NULL
#' @param tile.size size of tile borders. No border if is equal to NULL. Default: 0.1
#' @param tile.col colour of tile borders. No border if is equal to NULL. Default: NULL
#' @param legend.position legend position on plot. Default: 'right'
#' @param legend.height legend height. Default: 2
#' @param legend.width legend width. Default: 0.6
#' @param legend.rel relative size of legend text. Default: 0.9
#' @param legend.colour legend text colour. Default: 'black'
#' @param legend.breaks breaks to subset.genes on colour key. Default: NULL
#' @param legend.labels specify labels on colour key. Default: NULL
#' @param legend.title title of legend. Default: 'Inferred CNA `[log2 ratio`]'
#' @param legend.justification legend justification on plot. Default: 'top'
#' @param legend.title.position position of legend title on plot. Default: 'bottom'
#' @param legend.title.angle angle of legend title. Default: NULL
#' @param legend.title.rel relative size of legend title text. Default: 0.9
#' @return a list containing $p, the ggplot plot, and $data, the data (in dataframe format).
#' @seealso 
#'  \code{\link[reshape2]{melt}}
#'  \code{\link[ggpubr]{rotate_axis_text}}
#' @rdname cna_heatmap
#' @export 
#' @importFrom reshape2 melt
#' @importFrom ggpubr rotate_x_text rotate_y_text
cna_heatmap = function(cna,
                   limits = c(-1, 1),
                   ratio = 0.5,
                   #cols = heatCols, 
                   x.name = 'Chromosome',
                   y.name = 'Cell',
                   legend.title = 'Inferred CNA\n[log2 ratio]',
                   x.hide = c('Y'),
                   clones = T,
                   order.cells = F,
                   cutree.k=NULL,
                   cutree.h=NULL,
                   subset.genes = NULL,
                   cor.method = 'pearson',
                   dist.method = 'euclidean',
                   cluster.method = 'average',
                   angle = NULL,
                   x.angle = NULL,
                   y.angle = 0,
                   axis.rel = 1,
                   base.size = 12,
                   axis.title.size = 12,
                   axis.text.size = 11,
                   base.col = "#073642",
                   title = NULL,
                   subtitle = NULL,
                   caption = NULL,
                   text.size = 12,
                   y.hide = NULL,
                   tile.size = 0.1,
                   tile.col = NULL,
                   legend.position = 'right',
                   legend.height = 2,
                   legend.width = 0.6,
                   legend.rel = 0.9,
                   legend.colour = 'black',
                   legend.breaks = NULL,
                   legend.labels = NULL,
                   legend.justification = 'top',
                   legend.title.position = 'bottom',
                   legend.title.angle = NULL,
                   legend.title.rel = 0.9,
                   verbose=T) {
  
  # order genes by genomic position
  cna = order_genes(cna)
  # order cells and set y axis break and text
  group.ord = !is.logical(order.cells) && !is.null(order.cells)
  uni.ord = is.logical(order.cells) && isTRUE(order.cells)
  yLineBreaks = ggplot2::waiver()
  yTextBreaks = ggplot2::waiver()
  if (group.ord) {
    res.lst = hca_cells(cna,
                        subset.genes = subset.genes,
                        groups = order.cells,
                        cor.method = cor.method,
                        dist.method = dist.method,
                        cluster.method = cluster.method)
    cna = res.lst$mtx
    # breaks for horizontal lines denoting order.cells
    yLineBreaks = .cellBreaks(res.lst$groups)
    # and breaks for their labels
    yTextBreaks = .cellBreaks(res.lst$groups, halfway = TRUE, hide = x.hide)
  } else if (uni.ord) {
    res.lst = hca_cells(cna,
                        subset.genes = subset.genes,
                        groups = NULL,
                        cor.method = cor.method,
                        dist.method = dist.method,
                        cluster.method = cluster.method,
                        cutree.k=cutree.k, cutree.h=cutree.h, verbose=T)
    cna = res.lst$mtx
    if( !is.null(cutree.k) || !is.null(cutree.h) ) { 
    # res.lst$groups is a named(name labels) vector of memberships(1,2,3...) from stats:cutree
    yLineBreaks = .cellBreaks(res.lst$groups) #breaks for positions
    yTextBreaks = .cellBreaks(res.lst$groups, halfway = TRUE, hide = x.hide) #breaks for labels
    }
  } else if (clones) {
    res.lst = find_clones(cna, 
                          gene.quantile=0.5,
                          umap.n.neighbors=10, 
                          umap.spread=5, 
                          umap.min.dist=0.01,
                          cluster.methods='louvain',
                          cluster.k=20,
                          exclude.y=TRUE,
                          amp=0.2, 
                          del=-0.2,
                          verbose=TRUE)
    cna=res.lst$cna
    yLineBreaks = .cellBreaks(res.lst$groups) #breaks for positions
    yTextBreaks = .cellBreaks(res.lst$groups, halfway = TRUE, hide = x.hide) #breaks for labels
  }
  rm(res.lst)
  
  # prepare dataframe
  dat = reshape2::melt(as.matrix(cna))
  colnames(dat) = c('Gene', 'Cell', 'CNA')

  # breaks for vertical lines denoting chromosomes
  xLineBreaks = .chromosomeBreaks(levels(dat$Gene))
  # and breaks for their labels
  xTextBreaks = .chromosomeBreaks(levels(dat$Gene), halfway = TRUE, hide = x.hide)
  
  if (!is.null(angle)) {
    x.angle = angle
    y.angle = angle
  }
  if (isTRUE(x.angle)) x.angle = 45
  if (!is.null(x.angle)) x.angle = ggpubr::rotate_x_text(angle = x.angle)
  if (isTRUE(y.angle)) y.angle = 90
  if (!is.null(y.angle)) y.angle = ggpubr::rotate_y_text(angle = y.angle)
  legend = .axisSpacer(breaks = legend.breaks, labels = legend.labels, limits = limits)
  legend.breaks = legend$breaks
  legend.labels = legend$labels
  # arguments to be passed to .cna_heatmap()
  Args = mget(ls()[-1 * which(ls() %in% c('cna', 'x.hide', 'order.cells'))])
  print(str(Args)); print(names(Args))
  # make ggplot object
  G = do.call(.cna_heatmap, Args)
  list(p = G, data = tibble::as_tibble(dat))
}


