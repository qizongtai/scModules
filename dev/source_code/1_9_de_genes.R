#Differentially expressed genes
de_genes<-function(matrix,vector,no=50,just="BH",max_pv=0.05,volcano=FALSE){
  clusterlist<-unique(vector)  
  comparisons <- parallel::mclapply(clusterlist,mc.cores=min(10,length(clusterlist)), function(x){
  #comparisons <- lapply(clusterlist, function(x){
    logfc<- rowMeans(matrix[,vector==x]) - rowMeans(matrix[,vector!=x])
    pval  <- apply(matrix,1,function(y) t.test(y[vector==x],y[vector!=x])$p.value)
    padj<- p.adjust(pval, method = just)
    a<-cbind.data.frame(logfc,padj,x)
    b<-names(logfc)
    c<-cbind.data.frame(b,a)
    d<-subset(c, c$logfc > 0 & c$padj < max_pv)
    if(!is.null(no)){
      e<-head(d[order(d$logfc, decreasing = T),], n=no)}else{
        e<-d[order(d$logfc, decreasing = T),]
      }
    colnames(e)<-c("Gene","logFC","adj.pval","Cluster")
    e
    
  })
  deg<- do.call(rbind,comparisons)
  rownames(deg)<-seq(1,nrow(deg))
  if(isTRUE(volcano)){
    deg[deg$Cluster==clusterlist[2],"logFC"]<- -1*deg[deg$Cluster==clusterlist[2],"logFC"]
  }
  return(deg)
}
