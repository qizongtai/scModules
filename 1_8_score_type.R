#cell score (input log2 transformed)
geneset_score <- function(mtx.log2,genes,bin=100,center=TRUE,no_bins=30,verbose=TRUE){
  if (verbose) message("running 'geneset_score'... will return a named(cell ID) vector(score)")
  #check if the signature genes exist in the mtx.log2
  logical_genes <- genes %in% rownames(mtx.log2)
  overlap <- sum(logical_genes)
  unfound <- length(logical_genes) - overlap 
  if (verbose) message(overlap, " signature genes are matched with genes(rows) in mtx.log2")
  if (verbose) message(unfound, " signature genes are NOT matched with genes(rows) in mtx.log2")
  if (verbose) message("Note: only use the matched genes to calculate the score...")
  #extract genes by matched genes
  genes<-genes[logical_genes]
  #extract mtx.log2 by the matched genes
  a1<-mtx.log2[rownames(mtx.log2) %in% genes,] #a1<-mtx.log2[genes,]
  if(is.null(bin)){score<-colMeans(a1)}
  if(!is.null(bin)){
    aa<-rowMeans(mtx.log2)
    bb<-dplyr::ntile(aa,no_bins) #named int
    print(table(bb))
    names(bb)<-names(aa)
    fc<-parallel::mclapply(genes,function(x){
      a2<-bb[x] #get bin number; a2: named int
      b2<-sample(names(bb[bb==a2]),size = bin,replace = F) #sample 100 genes; b2: chr[1:100]
      a3<-mtx.log2[x,] #for signature gene, get the expression value of all cells; a3: named numeric
      b3<-colMeans(mtx.log2[b2,]) #for sampled 100 genes, get mean expression values within each cells; b3: named numeric
      a3 - b3 # (sig1-ref1(=mean(100ref))) for all cells
    })
    fc2<-do.call(rbind,fc)
    score<-colMeans(fc2) # colMeans(fc2) = mean[(sig1-ref1)+(sig2-ref2)+(sig3-ref3)+(sign-refn)] = mean(sig1+2+3+n) - mean(ref1+2+3+n) 
  }
  if(isTRUE(center)){score <- score-mean(score)}
  return(score)
}

#Classify cells by using celltype
type_score<-function(mtx.log2,typemat,gv='Gene',tv='Celltype',size=2,conf_int=0,
                     bin=NULL,center=TRUE,conf="scale",no_bins=30,verbose=TRUE){
  if (verbose) message("running 'type_score'... will return a list(type,typescorelist)")
  #Only genes existing in tested mtx.log2
  typemat<-typemat[typemat[,gv] %in% rownames(mtx.log2),]
  #Filter cell types with few genes
  tl<- typemat %>% dplyr::group_by(typemat[,tv]) %>% dplyr::summarise_all(length)
  # `typemat[, tv]`    Gene Program
  # <chr>             <int>   <int>
  # 1 EpiSen_A           46      46
  # 2 EpiSen_B           63      63
  # 3 G1S                86      86
  # 4 G2M                79      79
  # 5 Immune             37      37
  tl<-tl[tl$Gene>size,]
  tl2<-as.data.frame(tl[,1])
  c1<-as.character(tl2[,1])
  typemat<-typemat[typemat[,tv] %in% c1,]
  #Score cells
  t1<-unique(typemat[,tv]) #unique cell types: "EpiSen_A" "EpiSen_B" "OxPhos" "Immune" "pEMT_LQ" "Stress" "pEMT"
  #option1: use mclapply
  typescorelist<-parallel::mclapply(t1,function(x){
    typescore<-geneset_score(mtx.log2,typemat[typemat[,tv]==x,gv],bin=bin,center=center,no_bins = no_bins)
  })
  # List of 10
  # $ : Named num [1:100] -0.138 -0.2 -0.685 0.282 -0.144 ...
  # ..- attr(*, "names")= chr [1:100] "AAACCTGAGCAGACTG-1-OCSCC1-post-3" "AAACCTGCATACTACG-1-OCSCC1-post-3" 
  # $ : Named num [1:100] -0.0786 0.037 -0.1851 -0.0557 -0.1798 ...
  # ..- attr(*, "names")= chr [1:100] "AAACCTGAGCAGACTG-1-OCSCC1-post-3" "AAACCTGCATACTACG-1-OCSCC1-post-3" 
  typescores<-do.call(cbind,typescorelist)
  colnames(typescores)<-t1
  #                                     EpiSen_A      EpiSen_B      OxPhos        Immune      pEMT_LQ       Stress         pEMT
  # AAACCTGAGCAGACTG-1-OCSCC1-post-3 -0.13796066 -0.0785503091 -0.16322768 -0.0007544539 -0.462373240  0.785821734  0.374618397
  # AAACCTGCATACTACG-1-OCSCC1-post-3 -0.20028954  0.0370149013  0.45271896 -1.2621032554  0.080715449  1.270509287  1.470107801
  # AAACCTGGTCTGATCA-1-OCSCC1-post-3 -0.68496917 -0.1850738054 -0.30376915  0.9488453587 -0.664185762  0.410816362 -0.409109653
  # AAACCTGGTGGCGAAT-1-OCSCC1-post-3  0.28207583 -0.0556764897 -0.45954217 -0.8667186604  0.146998490 -0.202995813 -0.382384352
  # AAACCTGTCCACGTGG-1-OCSCC1-post-3 -0.14410518 -0.1797775145  0.29231518 -0.4140463928  0.177249945 -0.256873284 -0.357516895
  names(typescorelist)<-t1
  #Classify cells by type with the highest expression
  Type<-character(length=nrow(typescores))
  names(Type)<-rownames(typescores)
  for(i in 1:nrow(typescores)){
    a<- which.max(typescores[i,])
    b<- max(typescores[i,-a])
    c<- max(typescores[i,])
    if(conf_int!=0){
      d<-ifelse(conf=="scale",conf_int*c,conf_int)
      Type[i]<-ifelse( c-b > d, names(a),"Unresolved")}else{
        Type[i]<-names(a)}
  }
  return(list(type=Type,typescorelist=typescorelist))
}