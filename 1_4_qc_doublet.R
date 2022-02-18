detect_db <- function(mtx.count,methods=c("sce", "srt.dbfinder"), srt.pc.num=30, verbose=TRUE){
if (verbose) message("running 'db_detection'... will return a meta df with doublet info.")
###====doublet rate=====###
#method1
cell.num <- ncol(mtx.count)
dbr <- 0.004*cell.num/500
#method2 from scDblFinder https://bioconductor.org/packages/devel/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
#dbr <- 0.01*cell.num/1000

###==============================sce based doublet detection====================================###
if (verbose) message("running sce beased db detection: scran, scds and scDblFinder")
library(SingleCellExperiment)
sce <- SingleCellExperiment::SingleCellExperiment(assay=list(counts=mtx.count))

###===scDblFinder(scran;scdbfinder)===###
#sce <- scDblFinder::scDblFinder(sce) #default dbr=0.01*cell/1000; higher than 0.004*cell/500
sce <- scDblFinder::scDblFinder(sce, dbr=dbr)
#colData(sce); table(sce$scDblFinder.class)

###======scran======### Note: don't assign the result to sce
#new version of the function is maintained in scDblFinder
sce$scran.db <- scDblFinder::computeDoubletDensity(sce)
#old verison is deprecated. 
#sce$scran.db <- scran::doubletCells(sce)

###======scds======###
# #Annotate doublet using co-expression based doublet scoring:
# sce <- scds::cxds(sce)
# #Annotate doublet using binary classification based doublet scoring:
# sce <- scds::bcds(sce)
#Combine both annotations into a hybrid annotation
sce <- scds::cxds_bcds_hybrid(sce)
#check
#cd <- colData(sce); head(cbind(cd$cxds_score,cd$bcds_score, cd$hybrid_score))

###=======summarize sce based doublet detection=======###
db.num <- round(cell.num*dbr, digits = 0)
db.vec <- c(rep("doublet", db.num), rep("singlet", (cell.num-db.num)))
ref.colnames <- c("scran.db", "hybrid_score")
add.colnames <- c("scran.class", "hybrid.class")
#sort the colData by doublet socres and then add a new column that has doublets and singlets definition
sort_addcol <- function(colData, ref.colnames, add.colnames, db.vec) {
  for (i in seq_along(ref.colnames)) {
    colData <- colData[order(colData[[i]], decreasing = T),]
    colData[[add.colnames[i]]] <- db.vec
  }
  return(colData)
}
sce.colData <- sort_addcol(colData=SingleCellExperiment::colData(sce), ref.colnames=ref.colnames, add.colnames=add.colnames, db.vec=db.vec)

if(!"srt.dbfinder" %in% methods) {
###=====combine sce doublet results to one column(singlet/doublet)=====###
if (verbose) message("combining db results from scran, scds and scDblFinder")
# #combine three sce-based methods: scran, scds, scDblFinder
# #old code 
# sce.colData$sce.comb <- ifelse(sce.colData$scDblFinder.class=="doublet", 
#                                ifelse(sce.colData$scran.class=="doublet",
#                                       ifelse(sce.colData$hybrid.class=="doublet", "doublet", "singlet"), 
#                                       "singlet"), 
#                                "singlet")

#recode "singlet" and "doublet" as 1 and 2
scDblFinder.class.num = ifelse(sce.colData$scDblFinder.class=="doublet", 0, 1)
scran.class.num = ifelse(sce.colData$scran.class=="doublet", 0, 1)
hybrid.class.num = ifelse(sce.colData$hybrid.class=="doublet", 0, 1)
sce.comb.num = scDblFinder.class.num + scran.class.num + hybrid.class.num
sce.colData$sce.comb = ifelse(sce.comb.num>=2, "singlet", "doublet")

print(table(sce.colData$sce.comb))
message("total cell number: ", dim(sce.colData)[1])
return(as.data.frame(sce.colData))
#return(tibble::column_to_rownames(sce.colData, var="Row.names"))

} else {
###==============================seurat based doubletfinder====================================###
if (verbose) message("running seurat based doubletfinder")
library(Seurat)
library(DoubletFinder)
#setting parameters
pc.num=1:srt.pc.num
#seurat preprocessing
srt <- CreateSeuratObject(counts=mtx.count)
if (T) { srt <- SCTransform(srt) 
} else { srt <- NormalizeData(srt) %>% FindVariableFeatures() %>% ScaleData() }
srt <- RunPCA(srt) %>% FindNeighbors(dims = pc.num) %>% FindClusters(resolution=0.8) %>% RunUMAP(dims=pc.num)
#Doublet checking
#find the optimal pK
sweep.res.list <- paramSweep_v3(srt, PCs = pc.num, sct = T) #different combinations of pN and pK.
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

#excluding the doublets from the same cell type and optimazing the expected doublets rate
DoubletRate = dbr
if (T) {homotypic.prop <- modelHomotypic(srt$SCT_snn_res.0.8) #or use seurat_clusters
} else {homotypic.prop <- modelHomotypic(srt$SingleR_BE)} #better with celltype annotation
nExp_poi <- round(DoubletRate*ncol(srt))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

##use optmized paramers for doubletFinder
srt <- doubletFinder_v3(srt, PCs = pc.num, pN = 0.25, pK = pK_bcmvn,
                        nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
## plotting and results are in srt@meta.data
srt.dbfider.colname <- colnames(srt@meta.data)[grep("DF.classifications", colnames(srt@meta.data))]
p <- DimPlot(srt, reduction = "umap", group.by = srt.dbfider.colname); p
#saving plots
#if (!dir.exists(outdir[["doublet"]])) {dir.create(outdir[["doublet"]], recursive = TRUE)}
#ggsave(paste0(outdir[["doublet"]], "doublet.pdf"), p, width = 10, height = 10)

#export doublet and singlet results
#doublet.tb1 <- table(srt[[doublet.colname]]); print(doublet.tb1)
#write.table(doublet.tb1, paste0(outdir[["doublet"]], "doulbet.txt"), sep="\t", col.names = T, row.names = F, quote=F)

#delete doublet
#need to rename the original colname to a defined one as the original name changes as different parameters are used
#didn't figure out how to use "xxx" character to reference the column name. tired as.name. didn't work.
#srt.1 <- subset(srt, subset = "DF.classifications_0.25_0.29_593" == "Singlet")
#srt.1 <- subset(srt, subset = as.name("DF.classifications_0.25_0.29_593") == "Singlet")
if(F) {
  colnames(srt@meta.data)[grep("DF.classifications", colnames(srt@meta.data))] <- "doublet_finder"
  srt.1 <- subset(srt, subset = doublet_finder == "Singlet")
}

###=====combine sce and srt doublet results to one column(singlet/doublet)=====###
if (verbose) message("combining sce and srt doublet results")
#extract srt doublet info
srt.meta.df <- srt@meta.data[srt.dbfider.colname]
db.sce.srt.comb <- merge(as.data.frame(sce.colData), y=srt.meta.df, by=0)
#combine three sce-based methods: scran, scds, scDblFinder
db.sce.srt.comb$sce.comb <- ifelse(db.sce.srt.comb$scDblFinder.class=="doublet", 
                                   ifelse(db.sce.srt.comb$scran.class=="doublet",
                                          ifelse(db.sce.srt.comb$hybrid.class=="doublet", "doublet", "singlet"), 
                                          "singlet"), 
                                   "singlet")
#combine thre sce-based and seurat-based DoubletFinder
db.sce.srt.comb$all.comb <- ifelse(db.sce.srt.comb$scDblFinder.class=="doublet", 
                                   ifelse(db.sce.srt.comb$scran.class=="doublet",
                                          ifelse(db.sce.srt.comb$hybrid.class=="doublet", 
                                                 ifelse(db.sce.srt.comb[srt.dbfider.colname]=="Doublet", "doublet", "singlet"), 
                                                 "singlet"), 
                                          "singlet"), 
                                   "singlet")
message("sce.combo x srt.dbfinder")
print(table(db.sce.srt.comb$sce.comb, db.sce.srt.comb$all.comb))
#return(tibble::column_to_rownames(db.sce.srt.comb, var="Row.names"))
return(as.data.frame(sce.colData))
}
}

