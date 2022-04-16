arm_mean = function(cna, exclude.y=TRUE) {
  #return a named num
  cna.arm.mean = sapply(split_genes(cna, by="arm"), mean, simplify = F)
  cna.arm.mean[sapply(cna.arm.mean, is.nan)] = 0
  if(isTRUE(exclude.y)) { temp = unlist(cna.arm.mean); return ( temp[!names(temp) %in% c("Yp","Yq")] ) 
  } else { return ( unlist(cna.arm.mean) ) }
}

arm_class = function(cna, amp=0.2, del=-0.2) {
  #return a matrix; same as origianl input matrix
  cna[cna > amp] = 1
  cna[cna < del] = -1
  cna[cna <= amp & cna >= del] = 0
  return(cna)
}

identical_col_lst <- function(df) {
  #input data frame
  #return a list of col names that are duplciated
  cols_with_duplicates <- names(which(duplicated(t(df))))
  duplicate_column_groups <- lapply(cols_with_duplicates, 
                                    function(column) { names(which(apply(df, 2, identical, df[, column])))})
  unique(duplicate_column_groups)
}

# #remove duplicated clusters (columns that have the same chromosome amp/del)
# cna.cluster.arm.class[,!duplicated(as.list(as.data.frame(cna.cluster.arm.class)))]

find_clones = function(cna, 
                       gene.quantile=0.5,
                       umap.n.neighbors=10, 
                       umap.spread=5, 
                       umap.min.dist=0.01,
                       cluster.methods='louvain',
                       cluster.k=15,
                       exclude.y=TRUE,
                       amp=0.2, 
                       del=-0.2,
                       merge.duplicated=TRUE,
                       verbose=TRUE) {
  if (verbose) message('running find_clones...return a list of 2: #1 cell-ordered cna matrix; #2 sorted cluster id')
  
  #filter genes
  if (verbose) message('selecting ', gene.quantile, ' quantile of genes...')
  cna.topgene = cna[cnaHotspotGenes(cna, gene.quantile = gene.quantile),]
  
  #order genes by chromosomes
  if (isTRUE(exclude.y)) {
    if (verbose) message('excluding y chromosome genes...(note: the returned cna matrix is not affected)')
    genes = split_genes(x = rownames(cna.topgene), exclude.y = exclude.y) %>% unlist() %>% unname()
    cna.topgene = cna.topgene[genes,]
  }

  if (verbose) message('dimension reduction by umap...')
  cna.topgene.umap <- run_umap(mtx.hvg=cna.topgene, 
                               n.neighbors=umap.n.neighbors, 
                               spread=umap.spread, 
                               min.dist=umap.min.dist)
  if (verbose) message('graph-based clustering by ', cluster.methods, ' k=', cluster.k, '...')
  cna.topgene.umap.snn.cluster = snn_cluster(mtx = cna.topgene.umap, methods=cluster.methods, k = cluster.k, transposed = T)
  cna.topgene.umap.snn.cluster.wide = wide_gcluster(cna.topgene.umap.snn.cluster, prefix = "umap.snn")
  
  #build a list (cluster id: cell id)
  cluster.name = paste0('umap.snn.k', cluster.k, '.', cluster.methods)
  cluster.uniq.id = sort(unique(cna.topgene.umap.snn.cluster.wide[[cluster.name]]))
  names(cluster.uniq.id) <- cluster.uniq.id
  cluster.cell.id.lst = lapply(cluster.uniq.id, 
                                function(x) row.names(cna.topgene.umap.snn.cluster.wide)[cna.topgene.umap.snn.cluster.wide[[cluster.name]]==x] )
  
  if (verbose) message('calculating the mean of cna by chromosome arms for ', length(cluster.cell.id.lst), ' clusters...')
  cna.cluster.arm.mean = sapply(cluster.cell.id.lst, 
                                function(x) arm_mean(cna=cna.topgene[,x], exclude.y=exclude.y) )
  if (verbose) message('classifying chromosome arms: amp(1) > ', amp, '; del(-1) < ', del, '; no change(0)')
  cna.cluster.arm.class = arm_class(cna=cna.cluster.arm.mean, amp=amp, del=del)
  
  #return(cna.cluster.arm.mean)
  #return(cna.cluster.arm.class)
  
  if (isTRUE(merge.duplicated)) {
  if (verbose) message('finding duplicated clusters in the amp/del matrix.')
  dup.lst = identical_col_lst(cna.cluster.arm.class)
  print(dup.lst)
  if (length(dup.lst)==0) {
    message("no duplicated clusters were found in the amp/del matrix.")
  } else if (length(dup.lst)>0) {
    message("found ", length(dup.lst), " duplicated clusters in the amp/del matrix.")
    sapply(dup.lst, 
           function(x, name){cna.topgene.umap.snn.cluster.wide[[name]][ cna.topgene.umap.snn.cluster.wide[[name]] %in% x[-1] ] <<- x[1]}, 
           name=cluster.name)
    if (verbose) message('repeating the previous two steps')
    #build a list (cluster id: cell id)
    cluster.name = paste0('umap.snn.k', cluster.k, '.', cluster.methods)
    cluster.uniq.id = sort(unique(cna.topgene.umap.snn.cluster.wide[[cluster.name]]))
    #cluster.uniq.id = seq_along(unique(cna.topgene.umap.snn.cluster.wide[[cluster.name]]))
    names(cluster.uniq.id) <- cluster.uniq.id
    cluster.cell.id.lst = lapply(cluster.uniq.id, 
                                 function(x) row.names(cna.topgene.umap.snn.cluster.wide)[cna.topgene.umap.snn.cluster.wide[[cluster.name]]==x] )
    
    if (verbose) message('RE-calculating the mean of cna by chromosome arms for ', length(cluster.cell.id.lst), ' clusters...')
    cna.cluster.arm.mean = sapply(cluster.cell.id.lst, 
                                  function(x) arm_mean(cna=cna.topgene[,x], exclude.y=exclude.y) )
    if (verbose) message('RE-classifying chromosome arms: amp(1) > ', amp, '; del(-1) < ', del, '; no change(0)')
    cna.cluster.arm.class = arm_class(cna=cna.cluster.arm.mean, amp=amp, del=del)
  } else { stop("error: duplicated list!") }
  }
  
  #return(cna.cluster.arm.mean)
  #return(cna.cluster.arm.class)
  #return(list(mean=cna.cluster.arm.mean,wide=cna.topgene.umap.snn.cluster.wide))
  
  if (verbose) message('hc on the cluster level by mean of the chromosome amp/del...')
  res = hca(x = cna.cluster.arm.mean, hclust.end = T)
  #options1: use factor; option2: use match
  groups = factor(cna.topgene.umap.snn.cluster.wide[[cluster.name]], levels=res$order, ordered = T) %>% sort()
  groups = split(groups, groups) # groups: make a list of (Ord.factors)
  ordered.index = factor(cna.topgene.umap.snn.cluster.wide[[cluster.name]], levels=res$order, ordered = T) %>% order()
  
  return(list(cna=cna[,ordered.index], groups=groups))
  
}

#tried to subset the arm-splut matrix to avoid repeated running of split_genes but have errors below
#ms = split_genes(cna.hta69.epi.topgene, by="arm")
#ms is a list of matrixes; each matrix correspone to one chrom arm
#Note1: if no gene in the arm; the element is num[0,1:cellnumber]; class: matrix; dim: [0, cellnumber]
#Note2: if one gene in the arm; the element is named num [1:cellnumber]; class: numeric; dim: NULL; create errors if subset by []
#Note3: if > one gene in the arm; the element is named num [1:genenumber,1:cellnumber]; class: matrix; dim: [genenumber, cellnumber]
#both below 2 codes have errors due to Note2
#ms.1 <- lapply(ms, "[", ,cluster.list[[1]])
#ms.1 <- lapply(ms, function(x) {x[ ,cluster.list[[1]] ]})

# a = split_genes(x = rownames(cna.hta69.epi), exclude.y = TRUE)
# l = find_clones(cna = cna.hta69.epi)
# cna.cluster.arm.class = find_clones(cna = cna.hta69.epi)
# cna.cluster.arm.mean = find_clones(cna = cna.hta69.epi)
# cluster.cell.id.lst = find_clones(cna = cna.hta69.epi)
# 
# cna.cluster.arm.class = arm_class(cna=cna.cluster.arm.mean, amp=0.2, del=-0.2)
