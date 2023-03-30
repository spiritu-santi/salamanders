library(ape)
library(foreach)
library(taxize)
library(tidyverse)
library(magrittr)
library(phyloch)
library(strap)
library(phytools)
library(MetBrewer)
library(here)
library(sp)
library(rgdal)
library(patchwork)
library(ggtree)
library(ggnewscale)


#### Process alingment ####
fasta_seqs <- read.FASTA("Caudata_outgrp.outaln_v1.fas")
names(fasta_seqs)
input.tre = "data/RAxML_bestTree.Caudata_v5.tre"
tre <- read.tree(input.tre)
match(names(fasta_seqs),tre$tip.label) -> iden
class <- classification(tre$tip.label,db="ncbi")
new_name <- lapply(class, function(x) 
{ if(length(which(x$rank == "subspecies")) == 0) paste(x[x$rank %in% c("family","species"),"name"],collapse = "_")
  else paste(x[x$rank %in% c("family","subspecies"),"name"],collapse = "_")
} )  %>% unlist(.)  %>% gsub(" ","_",.)
names(new_name) <- NULL
tre$tip.label <- new_name
tre$tip.label[iden] -> names(fasta_seqs)
names(fasta_seqs) <- unlist(lapply(lapply(strsplit(names(fasta_seqs),"_"),"[",1:3),paste,collapse="_"))
dupas <- unique(names(fasta_seqs)[which(duplicated(names(fasta_seqs)))])
fasta_seqs_bu <- fasta_seqs
cat("Number of sequences",length(fasta_seqs),"\n")
for (dp in 1:length(dupas)){
dupi <- del.gaps(fasta_seqs[which(names(fasta_seqs) == dupas[dp])])
ll <- sapply(dupi,length)
gon <- which(names(fasta_seqs) == dupas[dp])[which(ll!=max(ll))]
fasta_seqs <- fasta_seqs[-gon]
cat("Number of sequences",length(fasta_seqs),"\n")
}
dupas <- unique(names(fasta_seqs)[which(duplicated(names(fasta_seqs)))])
cat("Number of sequences",length(fasta_seqs),"\n")
for (dp in 1:length(dupas)){
  dupi <- del.gaps(fasta_seqs[which(names(fasta_seqs) == dupas[dp])])
  ll <- sapply(dupi,length)
  gon <- which(names(fasta_seqs) == dupas[dp])[sample(1:2,1)]
  fasta_seqs <- fasta_seqs[-gon]
  cat("Number of sequences",length(fasta_seqs),"\n")
}
write.FASTA(fasta_seqs,"Caudata_outgrp.outaln_v2.fas")
unique(lapply(fasta_seqs,length))
#####
#### Read original alignment and index duplicate 16S clusters ####
fasta_seqs <- read.dna("Caudata_outgrp.outaln_v2.fas",format = "fasta")
second_16S <- fasta_seqs[, 144973 : 148798] ## positions come from the partitions file
second_16S <- second_16S[-which(unlist((lapply(del.gaps(second_16S),length))==0)),]
first_16S <- fasta_seqs[, 12052 : 16406] ## positions come from the partitions file
first_16S <- first_16S[-which(unlist((lapply(del.gaps(first_16S),length))==0)),]
nombres_16S <- c(rownames(first_16S),rownames(second_16S))
dups_16S <- nombres_16S[duplicated(nombres_16S)]

second_16S <- fasta_seqs[, 144973 : 148798]
first_16S <- fasta_seqs[, 12052 : 16406]
ll_first <- lapply(del.gaps(first_16S[which(rownames(first_16S)%in%dups_16S),]),length)
ll_second <- lapply(del.gaps(second_16S[which(rownames(second_16S)%in%dups_16S),]),length)

### Check which sequences are longer in first clusters relative to second cluster ##
### and discard them from second cluster. Write as unaligned fasta files ##
names((unlist(ll_first) > unlist(ll_second))[unlist(ll_first) > unlist(ll_second)]) -> dos
second_16S <- second_16S[-which(rownames(second_16S)%in% dos),]
second_16S <- second_16S[-which(unlist((lapply(del.gaps(second_16S),length))==0)),]
write.FASTA(del.gaps(second_16S),"data/Second16_toalign.fasta")

### Check which sequences are longer in second clusters relative to first cluster ##
### and discard them from first cluster. Write as aligned fasta files (without all-gap sequences ##
names((unlist(ll_first) > unlist(ll_second))[!unlist(ll_first) > unlist(ll_second)]) -> uno
first_16S <- first_16S[-which(rownames(first_16S)%in% uno),]
first_16S <- first_16S[-which(unlist((lapply(del.gaps(first_16S),length))==0)),]
write.FASTA(first_16S,"data/First16.fasta")

#### Algin in clustal-omega: first cluster as profile, second cluster unaligned ####
#### Before: add empty sequences (all gaps) for the two extra species
#### Before: add extra 16S sequences to second cluster
#### Incorporate 16S alingment into all-cluster alingment ####
fasta_news <- read.dna("data/output_16S.fasta",format = "fasta")
fasta_news <- fasta_news[,-1]
fasta_seqs <- read.dna("Caudata_outgrp.outaln_v2.fas",format = "fasta")
first_16S <- fasta_seqs[, 12052 : 16406]
which(is.na(match(rownames(fasta_news),rownames(fasta_seqs))))
fasta_seqs[match(rownames(fasta_news),rownames(fasta_seqs)),12052 : 16406]  <- fasta_news
first_16S[-which(unlist((lapply(del.gaps(first_16S),length))==0)),] ### before
first_16S <- fasta_seqs[, 12052 : 16406]
first_16S[-which(unlist((lapply(del.gaps(first_16S),length))==0)),] ### after!!
write.FASTA(fasta_seqs,"data/Caudata_outgrp.outaln_v2.1.fas")
#####
### Repeat for COI ####
fasta_seqs <- read.dna("data/Caudata_outgrp.outaln_v2.1.fas",format = "fasta")
second_16S <- fasta_seqs[,98360 : 98960] ## positions come from the partitions file
second_16S <- second_16S[-which(unlist((lapply(del.gaps(second_16S),length))==0)),]
first_16S <- fasta_seqs[, 29863 : 33980] ## positions come from the partitions file
first_16S <- first_16S[-which(unlist((lapply(del.gaps(first_16S),length))==0)),]
nombres_16S <- c(rownames(first_16S),rownames(second_16S))
dups_16S <- nombres_16S[duplicated(nombres_16S)]
dups_16S
rownames(second_16S)
#### Just eliminate second COI cluster !!!
first_16S <- first_16S[-which(rownames(first_16S)=="Plethodontidae_Chiropterotriton_nubilus"),]
write.FASTA(first_16S,"data/FirstCOI.fasta")
#### Algin in clustal-omega: first cluster as profile, extra sequences unaligned ####
#### Incorporate 16S alingment into all-cluster alingment ####
fasta_news <- read.dna("data/output_COI.fasta",format = "fasta")
fasta_news
fasta_seqs <- read.dna("data/Caudata_outgrp.outaln_v2.1.fas",format = "fasta")
first_16S <- fasta_seqs[, 29863 : 33980]
first_16S
which(is.na(match(rownames(fasta_news),rownames(fasta_seqs))))
fasta_seqs[match(rownames(fasta_news),rownames(fasta_seqs)), 29863 : 33980]  <- fasta_news
first_16S[-which(unlist((lapply(del.gaps(first_16S),length))==0)),] ### before
first_16S <- fasta_seqs[, 29863 : 33980]
first_16S[-which(unlist((lapply(del.gaps(first_16S),length))==0)),] ### after!!
write.FASTA(fasta_seqs,"data/Caudata_outgrp.outaln_v2.2.fas")

write.table(sort(rownames(fasta_seqs)),"data/constraints.tre",row.names = F,quote = F)


grep("Plethodontidae_Pseudoeurycea_lynchi",rownames(fasta_seqs[,5657:7597]))

fasta_seqs[325,5657:7597]


###### delete shorter duplicates from second cluster -> delete gaps
###### delete shorter duplicates from first cluster -> keep gaps
###### in theory clustal's profile alignment will keep the length, so I can easily replace sequences.......

##### Utilities ####
rename_tree <- function(input.tre) {
  tre <- read.tree(input.tre)
  class <- classification(tre$tip.label,db="ncbi")
 new_name <- lapply(class, function(x) 
  { if(length(which(x$rank == "subspecies")) == 0) paste(x[x$rank %in% c("family","species"),"name"],collapse = "_")
    else paste(x[x$rank %in% c("family","subspecies"),"name"],collapse = "_")
  } )  %>% unlist(.)  %>% gsub(" ","_",.)
  names(new_name) <- NULL
  tre$tip.label <- new_name
  write.tree(tre,file="RAxML_names.Caudata_v1.tre")
}
getphylo_x <- function(tree, node) {
  if(is.character(node)) {
    node <- which(c(tree$tip.label, tree$node.label) == node)
  }
  pi <- tree$edge[tree$edge[,2]==node, 1]
  if (length(pi)) {
    ei<-which(tree$edge[,1]==pi & tree$edge[,2]==node)
    tree$edge.length[ei] + Recall(tree, pi)
  } else {
    if(!is.null(tree$root.edge)) {
      tree$root.edge
    } else {
      0
    }
  }
}
getphylo_y <- function(tree, node) {
  if(is.character(node)) {
    node <- which(c(tree$tip.label, tree$node.label)==node)
  }
  ci <- tree$edge[tree$edge[,1]==node, 2]
  if (length(ci)==2) {
    mean(c(Recall(tree, ci[1]), Recall(tree, ci[2])))
  } else if (length(ci)==0) {
    Ntip <- length(tree$tip.label)
    which(tree$edge[tree$edge[, 2] <= Ntip, 2] == node)
  } else {
    stop(paste("error", length(ci)))
  }
}

mantel.correlog_2 <-
  function(D.eco, D.geo=NULL, XY=NULL, n.class=0, break.pts=NULL,
           cutoff=TRUE, r.type="pearson", nperm=999, mult="holm",
           progressive=TRUE) {
    r.type <- match.arg(r.type, c("pearson", "spearman", "kendall"))
    mult   <- match.arg(mult, c("sidak", p.adjust.methods))
    
    epsilon <- .Machine$double.eps
    D.eco <- as.matrix(D.eco)
    
    ## Geographic distance matrix
    if(!is.null(D.geo)) {
      if(!is.null(XY))
        stop("you provided both a geographic distance matrix and a list of site coordinates:\nwhich one should the function use?")
      D.geo <- as.matrix(D.geo)
    } else {
      if(is.null(XY)) {
        stop("you did not provide a geographic distance matrix nor a list of site coordinates")
      } else {
        D.geo <- as.matrix(dist(XY))
      }
    }
    
    n <- nrow(D.geo)
    if(n != nrow(D.eco))
      stop("numbers of objects in D.eco and D.geo are not equal")
    n.dist <- n*(n-1)/2
    vec.D <- as.vector(as.dist(D.geo))
    vec.DD <- as.vector(D.geo)
    
    ## Number of classes and breakpoints
    
    if(!is.null(break.pts)) {
      ## Use the list of break points
      if(n.class > 0)
        stop("you provided both a number of classes and a list of break points:\nwhich one should the function use?")
      n.class = length(break.pts) - 1
      
    } else {
      ## No breakpoints have been provided: equal-width classes
      if(n.class == 0) {
        ## Use Sturges rule to determine the number of classes
        n.class <- ceiling(1 + log(n.dist, base=2))
      }
      ## Compute the breakpoints from n.class
      start.pt <- min(vec.D)
      end.pt <- max(vec.D)
      width <- (end.pt - start.pt)/n.class
      break.pts <- vector(length=n.class+1)
      break.pts[n.class+1] <- end.pt
      for(i in 1:n.class)
        break.pts[i] <- start.pt + width*(i-1)
    }
    
    half.cl <- n.class %/% 2
    
    ## Move the first breakpoint a little bit to the left
    break.pts[1] <- break.pts[1] - epsilon
    
    ## Find the break points and the class indices
    class.ind <- break.pts[1:n.class] +
      (0.5*(break.pts[2:(n.class+1)]-break.pts[1:n.class]))
    
    ## Create the matrix of distance classes
    vec2 <- vector(length=n^2)
    for(i in 1:n^2)
      vec2[i] <- min( which(break.pts >= vec.DD[i]) ) - 1
    
    ## Start assembling the vectors of results
    class.index <- NA
    n.dist <- NA
    mantel.r <- NA
    mantel.p <- NA
    ## check.sums = matrix(NA,n.class,1)
    
    ## Create a model-matrix for each distance class, then compute a Mantel test
    for(k in 1:n.class) {
      class.index <- c(class.index, class.ind[k])
      vec3 <- rep(0, n*n)
      sel <- which(vec2 == k)
      vec3[sel] <- 1
      mat.D2 <- matrix(vec3,n,n)
      diag(mat.D2) <- 0
      n.dis <- sum(mat.D2)
      n.dist <- c(n.dist, n.dis)
      if(n.dis == 0) {
        mantel.r <- c(mantel.r, NA)
        mantel.p <- c(mantel.p, NA)
      } else {
        row.sums <- rowSums(mat.D2)
        ## check.sums[k,1] = length(which(row.sums == 0))
        if((cutoff==FALSE) ||
           !(cutoff==TRUE && k > half.cl && any(row.sums == 0))) {
          temp <- vegan::mantel(mat.D2, D.eco, method=r.type, permutations=nperm)
          mantel.r <- c(mantel.r, -temp$statistic)
          temp.p <- temp$signif
          
          ## The mantel() function produces a one-tailed p-value
          ## (H1: r>0) Here, compute a one-tailed p-value in
          ## direction of the sign
          if(temp$statistic < 0) {
            temp.p <- (sum(temp$perm <= temp$statistic)+1)/(nperm+1)
          }
          mantel.p <- c(mantel.p, temp.p)
        } else {
          mantel.r <- c(mantel.r, NA)
          mantel.p <- c(mantel.p, NA)
        }
      }
    }
    
    mantel.res <- cbind(class.index, n.dist, mantel.r, mantel.p)
    mantel.res <- mantel.res[-1,]
    
    ## Note: vector 'mantel.p' starts with a NA value
    mantel.p <- mantel.p[-1]
    n.tests <- length(mantel.p)  #length(which(!is.na(mantel.p))) 
    
    if(mult == "none") {
      colnames(mantel.res) <-
        c("class.index", "n.dist", "Mantel.cor", "Pr(Mantel)")
    } else {
      ## Correct P-values for multiple testing
      if(progressive) {
        p.corr <- mantel.p[1]
        if(mult == "sidak") {
          for(j in 2:n.tests)
            p.corr <- c(p.corr, 1-(1-mantel.p[j])^j)
        } else {
          for(j in 2:n.tests) {
            temp <- p.adjust(mantel.p[1:j], method=mult)
            p.corr <- c(p.corr, temp[j])
          }
        }
      } else {
        ## Correct all p-values for 'n.tests' simultaneous tests
        if(mult == "sidak") {
          p.corr <- 1 - (1 - mantel.p[1:n.tests])^n.tests
        } else {
          p.corr <- p.adjust(mantel.p[1:n.tests], method=mult)
        }
      }
      temp <- c(p.corr, rep(NA,(n.class-n.tests)))
      mantel.res <- cbind(mantel.res, temp)
      colnames(mantel.res) <-
        c("class.index", "n.dist", "Mantel.cor", "Pr(Mantel)", "Pr(corrected)")
    }
    rownames(mantel.res) <-
      rownames(mantel.res,do.NULL = FALSE, prefix = "D.cl.")
    
    ## Output the results
    res <- list(mantel.res=mantel.res, n.class=n.class, break.pts=break.pts,
                mult=mult, n.tests=n.tests, call=match.call() )
    class(res) <- "mantel.correlog"
    return(res)
  }
mpmcorrelogram_2 <-function (xdis, geodis, zdis=NULL, method="pearson",
                             alfa=0.05, nclass=NULL, breaks=NULL,
                             permutations=999, strata, simil=FALSE, plot=TRUE,
                             print=TRUE,method_adj="fdr"){
  
  
  xdis <- as.dist(xdis)
  ydis <- as.dist(geodis)
  if(!is.null(zdis)) zdis <- as.dist(zdis)
  
  # manage breaks and numbre of distace classes
  nclas <- NULL
  if(is.null(breaks)){
    # set the number of distance-classes
    if(!is.null(nclass)){
      nclas <- nclass
    }
    else {
      ## Sturge's law
      m <- length(ydis)
      nclas <- round(1 + 3.3*log10(m))
    }
    
    cmin <- range(ydis)[1]
    cmax <- range (ydis)[2]
    breaks <- seq(cmin, cmax, le = nclas+1)
    breaks[nclas+1] <- breaks[nclas+1]+1
  }
  if(is.null(nclas)) nclas <- length(breaks)-1 
  
  # inner function to compute p-values 
  pepe <- function(x) ifelse(x[1]>=0, length(which(x>=x[1]))/ length(x),
                             length(which(x<=x[1]))/ length(x))
  ydis <- as.matrix(ydis)
  geoclas <- NULL
  mantelr <- NULL
  estadistico <- NULL
  significatividad <- NULL
  pvalues <- NULL
  clases <- NULL
  
  # PArtial MAntel correlogram
  if(!is.null(zdis)){
    cat("evaluating distance class ")
    for(i in 1:nclas){
      cat(paste(i,", ",sep=""))
      geoclas[[i]] <- as.dist((ydis>=breaks[i] & ydis < breaks[i+1])*1)
      mantelr[[i]] <- vegan::mantel.partial(xdis, geoclas[[i]], zdis,
                                            method,permutations, strata)
      estadistico <- c(estadistico, 
                       ifelse(simil==FALSE, mantelr[[i]]$statistic*-1,
                              mantelr[[i]]$statistic))
      significatividad <- c(significatividad, mantelr[[i]]$signif)
      pvalues <- c(pvalues, pepe(c(mantelr[[i]]$statistic,
                                   mantelr[[i]]$perm)))
      clases <- c(clases, paste(signif(breaks[i],3)," - ",
                                signif(breaks[i+1],3)))
    }
    cat("\n\n")
  }
  else{
    # Mantel correlogram
    cat("evaluating distance class ")
    for(i in 1:nclas){
      cat(paste(i,", ",sep=""))
      geoclas[[i]] <- as.dist((ydis>=breaks[i] & ydis < breaks[i+1])*1)
      mantelr[[i]] <- mantel(xdis, geoclas[[i]], method, permutations, strata)
      estadistico <- c(estadistico,
                       ifelse(simil==FALSE, mantelr[[i]]$statistic*-1,
                              mantelr[[i]]$statistic))
      
      significatividad <- c(significatividad, mantelr[[i]]$signif)
      pvalues <- c(pvalues, pepe (c(mantelr[[i]]$statistic,
                                    mantelr[[i]]$perm)))
      clases <- c(clases, paste(round(breaks[i],3)," - ",round(breaks[i+1],3)))
    }
    cat("\n\n")
  }
  pvalues.adj <- p.adjust(pvalues,method=method_adj) #pvalues*(1:length(pvalues))
  
  if(plot==TRUE) {
    colores <- as.numeric(pvalues.adj < alfa)
    plot(1:nclas, estadistico, cex=1.5, pch=22-colores*7, type="b",
         ylab=expression(r[M]), xlab="distance classes", xlim= c(0,nclas +1),
         ylim=c(-1,1))
  }
  if(print==TRUE){
    print(data.frame(class = 1:nclas, distance.range=clases, rM = estadistico,
                     p =pvalues, p.Bonferroni=pvalues.adj))
    cat("\n\n")
  }
  result <- list(breaks=breaks, rM=estadistico, signif=significatividad,
                 pvalues=pvalues, pval.Bonferroni= pvalues.adj,
                 clases=clases)
  class(result) <- c("mpmcorrelogram", class(result))
  return(result)
}


Mode <- function(x) {
  ux <- unique(na.omit(x))
  if(length(ux) > 0) {names(which.max(table(x)))} else(NA)
}

### FUNCTIONS
mono.check<-function(input.tre="data/RAxML/best_tree/RAxML_bestTree.Caudata_v6.tre",datos_microbiomas = "data/Super_table_27_07_22_V2_60.csv",plot.PDF=TRUE,probabilidad = 80) {
  names <- strsplit(input.tre,"bipartitions.") %>% unlist(.) %>% .[2] %>% sub(".tre","",.)
  output = "FamilyGenusTable_"
  tre <- ape::read.tree(input.tre)
  tre <- ladderize(tre, TRUE)
  list <- tre$tip.label
  list_fams_un <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][1]} %>% unique(.)
  list_gens_un <- foreach(i=1:length(list),.combine=c) %do% {strsplit(list,"_")[[i]][2]} %>% unique(.)
  fams <- as.data.frame(matrix(nrow=length(list_fams_un),ncol=4))
  rownames(fams)<-list_fams_un; colnames(fams)<-c("MRCA_tips","Expected_tips","Monophyletic","Support")
  res <- as.data.frame(matrix(nrow=length(list_gens_un),ncol=4))
  rownames(res)<-list_gens_un; colnames(res)<-c("MRCA_tips","Expected_tips","Monophyletic","Support")
  for(i in 1:length(list_fams_un)){
    iden<-grep(list_fams_un[i],x = tre$tip.label)
    cat(i,"----",list_fams_un[i],"\n")
    if(length(iden) == 1) {
      fams[i,2]<-length(iden);fams[i,1]<-length(iden)
      fams[i,3]<-"YES"
      fams[i,4]<-NA
    }
    else{ 
      fams[i,2]<-length(iden)
      subtre <- extract.clade(tre,node=getMRCA(tre,tre$tip.label[iden]))
      cat("-----------NUMBER OF TIPS:",length(subtre$tip.label),"\n")
      fams[i,1]<-length(subtre$tip.label)
      fams[i,4]<-as.numeric(subtre$node.label[1])
      if(length(subtre$tip.label)-length(iden)!=0) {color="red"} else {color="black"}
      if(length(subtre$tip.label)-length(iden)!=0) {fams[i,3]<-"NO"} else {fams[i,3]<-"YES"}
    }}
  for(i in 1:length(list_gens_un)){
    iden<-grep(list_gens_un[i],x = tre$tip.label)
    cat(i,"----",list_gens_un[i],"\n")
    if(length(iden) == 1) {
      res[i,2]<-length(iden);res[i,1]<-length(iden)
      res[i,3]<-"YES"
      res[i,4]<-NA
    }
    else{ 
      res[i,2]<-length(iden)
      subtre <- extract.clade(tre,node=getMRCA(tre,tre$tip.label[iden]))
      cat("-----------NUMBER OF TIPS:",length(subtre$tip.label),"\n")
      res[i,1]<-length(subtre$tip.label)
      res[i,4]<-as.numeric(subtre$node.label[1])
      if(length(subtre$tip.label)-length(iden)!=0) {color="red"} else {color="black"}
      if(length(subtre$tip.label)-length(iden)!=0) {res[i,3]<-"NO"} else {res[i,3]<-"YES"}
    }}
  write.table(rbind(fams,res),file = paste("output/",output,names,".csv",sep=""), sep = ",",quote = F,row.names = T,col.names = T)
if(plot.PDF){
  cat("Making PDF of annotated tree","\t")
a <- ape::read.tree(input.tre) %>% ladderize(.,FALSE)
familias <- foreach(i=1:length(a$tip.label),.combine=c) %do% {strsplit(a$tip.label,"_")[[i]][1]} %>% unique(.) %>% as.data.frame(.)
generos <- foreach(i=1:length(a$tip.label),.combine=c) %do% {strsplit(a$tip.label,"_")[[i]][2]} %>% unique(.)
familias$Color <-rev(viridis::viridis(dim(familias)[1],end = 0.7))
familias$Color <- rainbow(11)
nodas <- lapply(familias$.,function(x) extract.clade(a, node=getMRCA(a,tip=a$tip.label[grep(paste(x,"_",sep=""),a$tip.label)]))$tip.label)
tcol <- tip.color(a,nodas,col=as.character(familias$Color))
ecol <- edge.color(a,nodas,col = as.character(familias$Color),what="stem")
w <- which(as.numeric(a$node.label) >= probabilidad)
w <- w + length(a$tip.label)
child <- (a$edge[which(a$edge[,2] %in% w),2])
parent <- (a$edge[which(a$edge[,2] %in% w),1])
x_child <- unlist(lapply(child,function(xx) getphylo_x(tree=a,node=xx)))
y_child <- unlist(lapply(child,function(xx) getphylo_y(tree=a,node=xx)))
x_parent <- unlist(lapply(parent,function(xx) getphylo_x(tree=a,node=xx)))
nodas <- data.frame("Child_node"=child,"Parent_node"=parent,"X_child"=x_child,"Y_child"=y_child,"X_parent"=x_parent)
pdf(file=paste("output/",names,".pdf",sep=""),width=20,height=180,onefile=F,useDingbats=FALSE)
plot(a,edge.width=2,tip.color="white",edge.color="black",show.tip.label=T,cex=1.5)
for (i in 1:length(x_child)) {
  lines(x = c(x_child[i],x_parent[i]),y=c(y_child[i],y_child[i]),lwd=11,col="black",lend="butt") }
tiplabels(text=a$tip.label,cex=1.5,col="black",frame = "rect",bg=tcol,adj=-.1)
micro <- readxl::read_xlsx(datos_microbiomas,sheet = 7)
a_micro <- unique(micro$Organism) %>% str_squish(.) %>% strsplit(.," ") %>% 
lapply(.,function(x)paste(x[1],x[2],sep="_")) %>% unlist(.) %>% unique(.)
a_macro <- sub(".*?_","",a$tip.label)
a_macro[grep("sp$",a_macro)] <- "Pseudoeurycea_sp."
which(is.na(match(a_micro,a_macro))) -> nana
a_micro[nana]
colores <- rep("white",length(a$tip.label))
cexa <- rep(0.1,length(a$tip.label))
colores[match(a_micro,a_macro)] <- "red"
cexa[match(a_micro,a_macro)] <- 3
tiplabels(pch=15,offset=0.02,cex=cexa,col=colores,frame="none")
dev.off()
}
}

### Process Data to map ####
do_map <- function(
  richness = "data/maps/Richness_10km_AMPHIBIANS_dec2017_EckertIV_Caudata.tif",
  data="data/Super_table_27_07_22_V2_60.csv",output="map.pdf"){
  raster::raster(here(richness)) -> anura
  raster::projectRaster(anura,crs="+proj=eck4 +lon_0=-90 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") -> anura
  raster::aggregate(anura,fact=5,fun=max)  %>% as(., "SpatialPixelsDataFrame") -> caudata
  projection = "+proj=eck4 +lon_0=-90 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  roads <- rnaturalearth::ne_coastline(scale = 110,returnclass = "sf") %>% sf::st_transform(., crs = projection)
  counts <- rnaturalearth::ne_countries(scale = 110,returnclass = "sf") %>%  sf::st_transform(., crs = projection)
  geos <- data.table::fread(here(data)) %>% group_by(long,lat,Organism) %>% mutate(Organism=sub(" ","_",Organism)) %>% summarise(.,n=n(),Family=first(Family),origin=first(origin)) %>% ungroup()
  g <- raster::aggregate(anura,fact=5,fun=max) %>% as(., 'SpatialPixels')
  geos %>% filter(!is.na(long)) %>% select(long,lat) %>% 
    sf::st_as_sf(x = ., coords = c("long", "lat"),crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% 
    sf::st_transform(.,crs = projection) %>%  sf::st_coordinates() %>% SpatialPoints(.,proj4string= CRS(projection)) %>%    sp::over(.,g) %>% enframe(.,name="name") %>% rename(., CellID=value) -> to_go
  
  geos %>% filter(!is.na(long)) %>% bind_cols(.,CellID=to_go %>% pull(CellID)) -> geos
  geos %>% group_by(Family,CellID) %>% summarize(n=sum(n),long=first(long),lat=first(lat)) %>% ungroup() %>% select(long,lat,n) %>% sf::st_as_sf(x = ., coords = c("long", "lat"),agr="n", 
   crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% sf::st_transform(.,crs =projection) -> geos_proj
  geos %>% group_by(Family, CellID) %>% summarise(Family=first(Family)) -> famas
  tibble(x=sf::st_coordinates(geos_proj)[,1],y=sf::st_coordinates(geos_proj)[,2],Family=famas$Family,n=geos_proj$n) -> ptos
  set.seed(1000)
  ptos %>% mutate(x_jit=jitter(x,factor=1000),y_jit=jitter(y,factor=1000)) -> ptos
  ggplot() + geom_sf(data=counts,colour="grey80",fill="grey80") +
    geom_tile(data=as.data.frame(caudata),aes(y=y, x=x,fill = Richness_10km_AMPHIBIANS_dec2017_EckertIV_Caudata)) +
    scale_fill_viridis(option = "mako",direction=-1,begin=0,end=1,name="SR")+
    theme(legend.position = "none") + geom_sf(data=roads,size=0.1) +
    geom_point(data=ptos,aes(x=x,y=y),size=0.2,color="black") +
    geom_segment(data=ptos,aes(x=x,y=y,xend=x_jit,yend=y_jit)) +
    geom_point(data=ptos,aes(x=x_jit,y=y_jit,color=Family,size=n)) +
    scale_color_manual(values = wesanderson::wes_palette("Darjeeling1",n = 5) )  +
    theme(legend.key.width = unit(2, "cm"),
          panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
          axis.text.x = element_blank(),axis.text.y = element_blank(),
          panel.background=element_rect(colour="black",fill=NA),
          panel.border=element_rect(colour="black",fill=NA))
ggsave(last_plot(),file=here("figures",output))
}

plot_alphas <- function(data="data/Super_table_27_07_22_V2_60.csv",output="figures_alpha.pdf"){ 
  a <- data.table::fread(here(data)) %>% as_tibble()
  a <- a %>% mutate(Organism=as.factor(Organism),Family=as.factor(Family))
  a$Family <- factor(a$Family,levels=c("Hynobiidae","Cryptobranchidae","Salamandridae","Ambystomatidae","Plethodontidae"))
  a %<>% mutate(origin = sub("-caught","",origin)) %>% 
    mutate(origin = ifelse(origin=="wild",origin,"captive"))
  wilx <- wilcox.test(shannon_entropy ~ morphotype,data=a,alternative="two")
  capt <- a %>% ggplot(aes(x=morphotype,y=shannon_entropy,fill=morphotype)) + 
    geom_boxplot(scale="count",trim=F,outlier.colour = "NA",width=0.3) +
    ggdist::stat_dots(aes(group=origin),side="left",justification=1.25,scale=0.7,size=0.1,dotsize=4,layout="bin")
  d <- ggplot_build(capt)$data[[2]] 
  p1 <- a %>% ggplot(aes(x=morphotype,y=shannon_entropy,fill=morphotype)) + 
    geom_boxplot(scale="count",trim=F,outlier.colour = "NA",width=0.3) +
    scale_fill_manual(values=wesanderson::wes_palette("Moonrise3")[1:2]) +
    ggdist::stat_dots(side="left",justification=1.25,scale=0.7,size=0.1,dotsize=4,layout="bin") +
    ggnewscale::new_scale_fill() +
    ggdist::stat_dots(aes(alpha=origin),fill="black",side="left",justification=1.25,scale=0.7,size=0.1,dotsize=4,layout="bin") +
    scale_alpha_manual(values=c(0.8,0)) +
    #stat_summary(geom="crossbar",fun=median,size=1,width=1) +
    theme(axis.title.x = element_text(size=14),axis.text.y = element_text(size=12),axis.text.x = element_text(size=10),axis.line = element_line(),legend.position = "none",panel.background = element_blank(),panel.grid = element_blank()) + scale_fill_manual(values=c(wesanderson::wes_palette("GrandBudapest2")[c(2,1,3)])) + 
    #geom_jitter(width=0.15,alpha=0.2,shape=19,size=2) + 
    labs(x="",y="Alpha diversity",title="Wilcoxon Test",subtitle = paste("W = ",wilx$statistic,", p-val = ",sprintf("%.4f",wilx$p.value),sep=""))
  wilx <- kruskal.test(shannon_entropy ~ Family,data=a)
  p2 <- a %>%
    ggplot(aes(x=Family,y=shannon_entropy,fill=Family)) + 
    geom_boxplot(scale="count",trim=F,outlier.colour = "NA",width=0.3) +
    ggdist::stat_dots(side="left",justification=1.25,scale=0.7,size=0.1,dotsize=4,layout="bin") +
    scale_fill_manual(values=c(wesanderson::wes_palette("Darjeeling1"))) +
    ggnewscale::new_scale_fill() +
    ggdist::stat_dots(aes(alpha=origin),fill="black",side="left",justification=1.25,scale=0.7,size=0.1,dotsize=4,layout="bin") +
    scale_alpha_manual(values=c(0.8,0)) +
    #stat_summary(geom="crossbar",fun=median,width=.5) + 
    theme(axis.title.x = element_text(size=14),axis.text.y = element_text(size=12),axis.text.x = element_text(size=8),axis.line = element_line(),legend.position = "none",panel.background = element_blank(),panel.grid = element_blank()) + 
    scale_x_discrete(guide = guide_axis(n.dodge=2)) +
    labs(x="",y="Alpha diversity",title="Kruskal-Wallis Test",subtitle = paste("Chi-2 = ",round(wilx$statistic,1),", p-val = ",sprintf("%.4f",wilx$p.value),sep="")) +
    NULL
  
  wilx <- kruskal.test(shannon_entropy ~ Organism,data=a)
  p3 <- a %>% ggplot(aes(x=fct_reorder(.f=Organism,.x=Family,.fun=first),y=shannon_entropy,fill=Family)) + geom_boxplot(width=0.3,outlier.colour = "NA",size=0.1) + 
    theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=9),axis.line = element_line(),legend.position = "none",panel.background = element_blank(),panel.grid = element_blank(),axis.text.y = element_text(size=8,face="italic")) + 
    scale_fill_manual(values=c(wesanderson::wes_palette("Darjeeling1"))) + 
    ggdist::stat_dots(fill="black",side="left",justification=1.25,scale=0.7,size=0.4,dotsize=4,layout="bin") +
    #geom_jitter(width=0.05,alpha=0.2,shape=19,size=1) + labs(x="",y="Faith's phylogenetic diversity") + 
    coord_flip() +
    labs(x="",y="Faith's phylogenetic diversity",title="Kruskal-Wallis Test",subtitle = paste("Chi-2 = ",round(wilx$statistic,1),", p-val = ",sprintf("%.4f",wilx$p.value),sep=""))
  
  layout="
  AACC
  BBCC"
  
  pp <- wrap_plots(list(p1,p2,p3)) +
    plot_layout(design=layout,byrow=T,nrow = 2,ncol=1,guides = "collect",tag_level = 'new') + 
    plot_annotation(title = "Alpha diversity",subtitle="By habitat and family") & theme(legend.position = '')
  ggsave(filename = here("figures",output),plot =  pp)
  
}

plot_alpha_phylo <- function(tree="data/trees/treePL_RAxML_bestTree.Caudata_v6.tre",
                             data="data/Super_table_27_07_22_V2_60.csv",
                             output="pruned_tree.pdf"){ 
  tree <- ape::read.tree(tree)
  all_data <- data.table::fread(data)%>% as_tibble()
  geos <- all_data %>% group_by(Organism) %>% mutate(V1=as.factor(V1)) %>% 
    mutate(Organism=sub(" ","_",Organism)) %>% summarise(across(where(is.numeric),mean,na.rm=T),Family=first(Family))
  tipas <- unlist(lapply(lapply(strsplit(tree$tip.label,"_"),"[",2:3),paste,collapse="_"))
  stopifnot("Error: tree and database not fully matched!"=geos$Organism[which(is.na(match(geos$Organism,tipas)))]==0)
  fin_tip <- tipas %>% match(geos$Organism,.) %>% .[!duplicated(.)]
  subtree <- keep.tip(tree,tree$tip.label[fin_tip])
  subtree <- ladderize(subtree,FALSE)
  subtree$root.time <- max(branching.times(subtree),na.rm=T)
  geos %>% distinct(Family) %>% pull() -> familias
  nodas <- lapply(familias[familias!="Hynobiidae"],function(x) extract.clade(subtree,node = getMRCA(subtree,tip=grep(paste(x,"_",sep=""),subtree$tip.label) %>% subtree$tip.label[.]))) %>% lapply(.,function(x) x$tip.label)
nodas[[5]] <- nodas[[4]]
nodas[[4]] <- grep(paste(familias[familias=="Hynobiidae"],"_",sep=""),subtree$tip.label) %>% subtree$tip.label[.]
names(nodas) <- familias
Color <- wesanderson::wes_palette("Royal2")[c(4,2,3,1,5)]
tcol <- tip.color(subtree,nodas,col=Color)
ecol <- edge.color(subtree,nodas,col=Color,what="stem")
  
geos <- geos %>% mutate(Color = tcol) %>% mutate(Tips = paste(Family,"_",Organism,sep=""),.before=Organism)
  nodes <- lapply(familias,function(x) getMRCA(subtree,tip=grep(paste(x,"_",sep=""),subtree$tip.label)))
  nodes[[4]] <- subtree$edge %>% as_tibble() %>% filter(V2== grep(paste(familias[familias=="Hynobiidae"],"_",sep=""),subtree$tip.label)) %>% pull(V1)
subtree <- groupClade(subtree,nodes)
subtree$tip.label %>% as_tibble() %>% separate(value,into=c("fam","sp"),sep="_",extra = "merge") %>% pull(sp) -> subtree$tip.label
df <- geos %>% mutate(Tips=Organism)
p1 <- ggtree(subtree,as.Date=F, right=F,branch.length = "time",root.position = -max(branching.times(subtree))) %<+% df + geom_tippoint(aes(color=Family),size=3,shape=19) + 
  geom_tiplab(size = 3,offset=3) + scale_color_manual(values=Color) + theme_tree2(legend.position="bottom") +
  guides(color=guide_legend(nrow=2)) +
  geom_point(aes(x=-100,y=4.7),shape=21,fill=Color[1],size=4,stroke=0.2) +
  geom_point(aes(x=-150,y=2.6),shape=21,fill=Color[2],size=4,stroke=0.2) +
  geom_point(aes(x=-165,y=1),shape=21,fill=Color[3],size=4,stroke=0.2) +
  geom_point(aes(x=-130,y=22.5),shape=21,fill=Color[4],size=4,stroke=0.2) +
  geom_point(aes(x=-105,y=8.6),shape=21,fill=Color[5],size=4,stroke=0.2) +
  NULL

ggsave(filename = here("figures",output),plot =  p1) 
  
#df <- geos %>% mutate(Tips=Organism) %>% as.data.frame()
#  row.names(df) <- df$Tips
#pp <- gheatmap(p1,df[,"faith_pd",drop=F],offset=100,width=.1,colnames = F,colnames_position = "top",colnames_offset_y=1) + 
#    scale_fill_fermenter(n.breaks=8,type = "seq",palette="OrRd",direction = 1,name="PD")  +
#    theme(legend.position="bottom")
}

plot_habitat_relative <- function(data="data/Super_table_27_07_22_V2_60.csv",
                                  output="habitat_relative.pdf"){
  
  all_data <- data.table::fread(data)%>% as_tibble()
  all_data %>% mutate(V1=as.factor(V1)) %>% 
    mutate(Organism=sub(" ","_",Organism)) %>% group_by(Organism,morphotype) %>% mutate(Organism=sub(" ","_",Organism)) %>% summarise(Sample=n()) %>% pivot_wider(names_from = morphotype,values_from = Sample) %>% mutate(across(where(is.numeric),replace_na,0)) %>% mutate(Total=Aquatic+Terrestrial) %>% mutate(RelAquatic=Aquatic/Total,RelTerrestrial=Terrestrial/Total) %>% dplyr::select(Organism,RelAquatic,RelTerrestrial) %>% pivot_longer(cols = 2:3,names_to = "Habitat",values_to = "Relative") -> habitat 
  p1 <- habitat %>% ggplot(aes(fill=Habitat, y=Relative, x=Organism)) + 
    geom_bar(position="stack", stat="identity") + theme(legend.position = "right") + 
    coord_flip()
  ggsave(filename = here("figures",output),plot =  p1) 
}

plot_betas <- function(data="data/Super_table_27_07_22_V2_60.csv",output="beta_diversity.pdf"){
distance_matrices <- list.files(path="data",pattern = "*_distance*",full.names = T)
data.table::fread(data) %>% as_tibble() -> geos
matrices <- lapply(distance_matrices,FUN = function(x){data.table::fread(x) %>% as.data.frame()})
for (i in 1:length(matrices)){
  rownames(matrices[[i]]) <- matrices[[i]][,1]
  matrices[[i]] <- matrices[[i]][,-1]
}
names(matrices) <- sub("-matrix.tsv","",distance_matrices) %>% substring(.,6)
nmds <- list()
pcoas <- list()

a <- data.table::fread(here(data)) %>% as_tibble()
#a <- a %>% filter(region_16S=="V4")
for (i in 1:length(matrices)){
  iden <- which(is.na(match(rownames(matrices[[i]]),a$Sample_id)))
vegan::metaMDS(dist(matrices[[i]][-iden,-iden]),sfgrmin = 1e-9,k=6,maxit=2000,try=20,engine="monoMDS",weakties = TRUE, noshare=0.2) -> nmds[[i]]
dist(matrices[[i]][-iden,-iden]) %>% pcoa(.,correction="none") -> pcoas[[i]]
}

i=2
iden <- which(is.na(match(rownames(matrices[[i]]),a$Sample_id)))
wuf_hab <- #nmds[[i]]$points %>% as_tibble(rownames = "Sample_id") %>% rename(Axis.1 = MDS1, Axis.2 = MDS2) %>% 
  pcoas[[i]]$vectors %>% as_tibble() %>% mutate(Sample_id=rownames(matrices[[i]][-iden,-iden]),.before=Axis.1) %>% dplyr::select(Sample_id:Axis.5) %>% 
  left_join(.,a,by="Sample_id") %>%
  ggplot(aes(x=Axis.1,y=Axis.2,color=morphotype,fill=morphotype)) +
  geom_point()

wuf_fam <- #nmds[[i]]$points %>% as_tibble(rownames = "Sample_id") %>% rename(Axis.1 = MDS1, Axis.2 = MDS2) %>% 
  pcoas[[i]]$vectors %>% as_tibble() %>% mutate(Sample_id=rownames(matrices[[i]][-iden,-iden]),.before=Axis.1) %>% dplyr::select(Sample_id:Axis.5) %>% 
  left_join(.,a,by="Sample_id") %>%
  ggplot(aes(x=Axis.1,y=Axis.2,color=Family,fill=Family)) +
  geom_point()

i=3
iden <- which(is.na(match(rownames(matrices[[i]]),a$Sample_id)))
uf_hab <- #nmds[[i]]$points %>% as_tibble(rownames = "Sample_id") %>% rename(Axis.1 = MDS1, Axis.2 = MDS2) %>% 
  pcoas[[i]]$vectors %>% as_tibble() %>% mutate(Sample_id=rownames(matrices[[i]][-iden,-iden]),.before=Axis.1) %>% dplyr::select(Sample_id:Axis.5) %>% 
  left_join(.,a,by="Sample_id") %>%
  ggplot(aes(x=Axis.1,y=Axis.2,color=morphotype,fill=morphotype)) +
  geom_point()

uf_fam <- #nmds[[i]]$points %>% as_tibble(rownames = "Sample_id") %>% rename(Axis.1 = MDS1, Axis.2 = MDS2) %>% 
  pcoas[[i]]$vectors %>% as_tibble() %>% mutate(Sample_id=rownames(matrices[[i]][-iden,-iden]),.before=Axis.1) %>% dplyr::select(Sample_id:Axis.5) %>% 
  left_join(.,a,by="Sample_id") %>%
  ggplot(aes(x=Axis.1,y=Axis.2,color=Family,fill=Family)) +
  geom_point()


layout="
  AACC
  BBDD"

pp <- wrap_plots(list(wuf_fam,wuf_hab,uf_fam,uf_hab)) +
  plot_layout(design=layout,byrow=T,nrow = 2,ncol=2,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = "Beta diversity",subtitle="By habitat and family") & theme(legend.position = 'bottom')
pp

# tabla con coordenadas del centroide y puntos ordenacion y hacer geom_segments

ggsave(filename = here("figures",output),plot =  pp)

sink(here("permanovas_beta.txt"))
i=2
iden <- which(is.na(match(rownames(matrices[[i]]),a$Sample_id)))
w_unifrac <- dist(matrices[[i]][-iden,-iden])
vegan::adonis2(w_unifrac ~ morphotype,data=as.data.frame(a))
vegan::adonis2(w_unifrac ~ Family,data=as.data.frame(a))

i=3
iden <- which(is.na(match(rownames(matrices[[i]]),a$Sample_id)))
uw_unifrac <- dist(matrices[[i]][-iden,-iden])
vegan::adonis2(uw_unifrac ~ morphotype,data=as.data.frame(a))
vegan::adonis2(uw_unifrac ~ Family,data=as.data.frame(a))
sink()

}

plot_betas_bray <- function(
  data="data/Super_table_27_07_22_V2_60.csv",
  output="beta_BC_diversity.pdf"){
  distance_matrices <- list.files(path="data/bc_samples",pattern = ".csv$*",full.names = T)
  data.table::fread(data) %>% as_tibble() -> geos
  matrices <- lapply(distance_matrices,FUN = function(x){data.table::fread(x) %>% as.data.frame()})
  for (i in 1:length(matrices)){
    rownames(matrices[[i]]) <- matrices[[i]][,1]
    matrices[[i]] <- matrices[[i]][,-1]
  }
  names(matrices) <- sub("data/bc_samples/","",distance_matrices) %>% substring(.,1,7)
  nmds <- list()
  pcoas <- list()
  
  a <- data.table::fread(here(data)) %>% as_tibble()
  #a <- a %>% filter(!Study%in%c("PRJNA830991","PRJNA505069"))
  for (i in 1:length(matrices)){
    cat(i,"\n")
    iden <- which(is.na(match(rownames(matrices[[i]]),a$Sample_id)))
  #vegan::metaMDS(dist(matrices[[i]][-iden,-iden]),sfgrmin = 1e-7,k=4,maxit=2000,engine="monoMDS",weakties = TRUE, noshare=0.2) -> nmds[[i]]
dist(matrices[[i]][-iden,-iden]) %>% pcoa(.,correction="lingoes") -> pcoas[[i]]
  }
  names(matrices)
  i=3
  iden <- which(is.na(match(rownames(matrices[[i]]),a$Sample_id)))
  #nmds[[i]]$points %>% as_tibble(rownames = "Sample_id") %>%
    pcoas[[i]]$vectors %>% as_tibble() %>% mutate(Sample_id=rownames(matrices[[i]][-iden,-iden]),.before=Axis.1) %>% dplyr::select(Sample_id:Axis.5) %>% 
    left_join(.,a,by="Sample_id") %>%
    ggplot(aes(x=Axis.1,y=Axis.2,color=morphotype,fill=morphotype)) +
    geom_point() +
      #scale_fill_stepsn(colors=MetBrewer::met.brewer("VanGogh3"),n.breaks=10)+
      #scale_color_stepsn(colors=MetBrewer::met.brewer("VanGogh3"),n.breaks=10) +
      NULL
  

    #nmds[[i]]$points %>% as_tibble(rownames = "Sample_id") %>% rename(Axis.1 = MDS1, Axis.2 = MDS2) %>% 
    pcoas[[i]]$vectors %>% as_tibble() %>% mutate(Sample_id=rownames(matrices[[i]][-iden,-iden]),.before=Axis.1) %>% dplyr::select(Sample_id:Axis.5) %>% 
    left_join(.,a,by="Sample_id") %>%
    ggplot(aes(x=Axis.1,y=Axis.2,color=Family,fill=Family)) +
    geom_point()
  
 
  
  
  
  
  
  layout="
  AACC
  BBDD"
  
  pp <- wrap_plots(list(wuf_fam,wuf_hab,uf_fam,uf_hab)) +
    plot_layout(design=layout,byrow=T,nrow = 2,ncol=2,guides = "collect",tag_level = 'new') + 
    plot_annotation(title = "Beta diversity",subtitle="By habitat and family") & theme(legend.position = 'bottom')
  
  ggsave(filename = here("figures",output),plot =  pp)
  
  sink(here("permanovas_beta.txt"))
  i=3
  iden <- which(is.na(match(rownames(matrices[[i]]),a$Sample_id)))
  mata <- dist(matrices[[i]][-iden,-iden])
  names(matrices)[[i]]
  vegan::adonis2(mata ~ morphotype,data=as.data.frame(a))
  vegan::adonis2(mata ~ Family,data=as.data.frame(a))
  
  i=3
  iden <- which(is.na(match(rownames(matrices[[i]]),a$Sample_id)))
  uw_unifrac <- dist(matrices[[i]][-iden,-iden])
  vegan::adonis2(uw_unifrac ~ morphotype,data=as.data.frame(a))
  vegan::adonis2(uw_unifrac ~ Family,data=as.data.frame(a))
  sink()
  
}

core_plot <- function(core="core_4.csv",
                      persistence = "Persistencia_lvl_4.csv",
                      abundance="Abundancia_lvl_4.csv",
                      presence="Persistencia_aparicion_lvl_4.csv",
                      tree="data/trees/treePL_RAxML_bestTree.Caudata_v6.tre",
                      data="data/Super_table_27_07_22_V2_60.csv"){ 
data.table::fread(here("data/core",core)) %>% as_tibble()  -> pam
data.table::fread(here("data/core",persistence)) %>% as_tibble()  -> pers
data.table::fread(here("data/core",abundance)) %>% as_tibble()  -> abund
data.table::fread(here("data/core",presence)) %>% as_tibble()  -> presencia

presencia %>% rename(Taxa = V1) %>% rowwise() %>% mutate(Total = 1 - sum(c_across(-Taxa)/1091),.before=Taxa) %>% #filter(Total>=0.80) %>% 
  distinct(Taxa,.keep_all = T) %>% select(Total:Taxa) -> core_80
  tree <- ape::read.tree(tree)
all_data <- data.table::fread(data)%>% as_tibble()
  geos <- all_data %>% group_by(Organism) %>% mutate(V1=as.factor(V1)) %>% 
    mutate(Organism=sub(" ","_",Organism)) %>% summarise(across(where(is.numeric),mean,na.rm=T),Family=first(Family))
  tipas <- unlist(lapply(lapply(strsplit(tree$tip.label,"_"),"[",2:3),paste,collapse="_"))
  stopifnot("Error: tree and database not fully matched!"=geos$Organism[which(is.na(match(geos$Organism,tipas)))]==0)
  fin_tip <- tipas %>% match(geos$Organism,.) %>% .[!duplicated(.)]
  subtree <- keep.tip(tree,tree$tip.label[fin_tip])
  subtree <- ladderize(subtree,FALSE)
  subtree$root.time <- max(branching.times(subtree),na.rm=T)

  abund %>% rename(Taxa = V1) %>% 
    pivot_longer(-1,names_to = "Species",values_to = "Abundances") %>% 
    mutate(Abundances=ifelse(Abundances!=0,1,0)) %>% group_by(Taxa) %>% summarise(n=sum(Abundances)) %>% mutate(Prop = n/41*100) %>% arrange(desc(Prop)) %>% left_join(.,pers %>% rename(Taxa = V1),by="Taxa") %>% filter(Prop==100) %>% pivot_longer(-c(1:3),names_to = "Salamander",values_to = "Persistence") %>%  left_join(.,core_80,by="Taxa") %>% relocate(Total,.after=Taxa) %>% filter(!grepl(";_",Taxa)) %>% group_by(Taxa) %>% 
    summarize(median=median(Persistence),min=min(Persistence),low = quantile(Persistence,0.1)) %>% arrange(min)
  
  
abund %>% rename(Taxa = V1) %>% 
  pivot_longer(-1,names_to = "Species",values_to = "Abundances") %>% 
   mutate(Abundances=ifelse(Abundances!=0,1,0)) %>% group_by(Taxa) %>% summarise(n=sum(Abundances)) %>% mutate(Prop = n/41*100) %>% arrange(desc(Prop)) %>% left_join(.,pers %>% rename(Taxa = V1),by="Taxa") %>% filter(Prop==100) %>% pivot_longer(-c(1:3),names_to = "Salamander",values_to = "Persistence") %>% 
  left_join(.,core_80,by="Taxa") %>% relocate(Total,.after=Taxa) %>% filter(!grepl(";_",Taxa)) %>%
  ggplot(aes(y=Persistence,x=fct_reorder(.f=Taxa,.x=Persistence,.fun=median),fill=Total)) +
  geom_boxplot(width=0.3,outlier.colour = "NA") + 
  theme(axis.title.x = element_text(size=12),
        axis.text.y = element_blank(),
        axis.line = element_line(),
        legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.key.width = unit(0.7,"in")) + 
  scale_fill_stepsn(colors=MetBrewer::met.brewer("Hiroshige"),n.breaks=10,name="Prevalence across\nall samples") + 
  #scale_color_stepsn(colors=colorspace::darken(MetBrewer::met.brewer("Hiroshige"),amount=0.5),n.breaks=10,name="Prevalence across\nall samples") + 
  ggdist::stat_dots(fill="black",side="left",justification=1.5,scale=0.5,size=1,dotsize=2,layout="bin",color="black") +
  coord_flip() + 
  scale_y_continuous(breaks=c(0,30,60,80,90,95,97,99,99.5,99.9)) +
  labs(y="Prevalence within species",x="Bacterial taxa") + 
  NULL

aver <- abund %>% rename(Taxa = V1) %>% 
  pivot_longer(-1,names_to = "Species",values_to = "Abundances") %>% 
  mutate(Abundances=ifelse(Abundances!=0,1,0)) %>% group_by(Taxa) %>% summarise(n=sum(Abundances)) %>% mutate(Prop = n/41*100) %>% arrange(desc(Prop)) %>% left_join(.,pers%>% rename(Taxa = V1),by="Taxa") %>% filter(Prop==100) %>% 
  pivot_longer(-c(1:3),names_to = "Salamander",values_to = "Persistence") %>% 
  left_join(.,core_80,by="Taxa") %>% relocate(Total,.after=Taxa) %>% filter(!grepl(";_",Taxa))
targets <- aver %>% distinct(Taxa) %>% pull()

dosa <- abund %>% rename(Taxa = V1) %>% 
  pivot_longer(-1,names_to = "Species",values_to = "Abundances") %>% filter(Taxa%in%all_of(targets)) %>% 
  left_join(.,core_80,by="Taxa") %>% relocate(Total,.after=Taxa) %>% filter(!grepl(";_",Taxa)) %>% rename(Salamander=Species)

geos %>% mutate(Salamander=sub("_"," ",Organism)) %>% left_join(dosa,.,by="Salamander") %>% 
  filter(Taxa=="Bacteria;Proteobacteria;Alphaproteobacteria;Caulobacterales") %>% 
  ggplot(aes(x=pre,y=Abundances,color=Taxa)) +
  geom_point() +
  geom_smooth(method="lm",se=FALSE)+
  #geom_smooth(method="gam",method.args=list(family="ziplss"),se = FALSE) +
theme(axis.title.x = element_text(size=12),
        axis.text.y = element_blank(),
        axis.line = element_line(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.key.width = unit(0.7,"in")) + 
  scale_fill_stepsn(colors=MetBrewer::met.brewer("Hiroshige"),n.breaks=10,name="Prevalence across\nall samples") +
  NULL

geos %>%mutate(Salamander=sub("_"," ",Organism)) %>% left_join(dosa,.,by="Salamander") %>% distinct(Taxa)
geos %>%mutate(Salamander=sub("_"," ",Organism)) %>% left_join(dosa,.,by="Salamander") %>% 
  #filter(Taxa=="Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales") %>% 
  mutate(Depth = 2500) %>% filter(!is.na(bio17))%>% mutate(Abundances= round(Abundances*2500/100)) -> test
library(phyloseq)
library(corncob)
test %>% select(Taxa,Salamander,Abundances) %>% pivot_wider(names_from = Salamander,values_from = Abundances) %>%
  select(-1) %>% as.matrix(.) %>% otu_table(.,taxa_are_rows = T) -> otus
test %>% select(Taxa,Salamander,Abundances) %>% pivot_wider(names_from = Salamander,values_from = Abundances) %>% select(1) %>% as.matrix(.) %>% tax_table(.) -> taxa
test %>% select(Salamander,Family,tm_max:bio19) %>% distinct(Salamander,.keep_all = TRUE) %>% as.data.frame(.) -> samples
row.names(samples) <- samples$Salamander
samples <- sample_data(samples)
wow <- merge_phyloseq(otus,taxa,samples)
bbdml(formula = sp9 ~ pre, phi.formula = ~ pre,data=wow)
wow_da <- differentialTest(formula= ~ pre, phi.formula = ~ pre,
                 formula_null = ~ 1,
                 phi.formula_null = ~ pre,
                 test="Wald",fdr_cutoff = 0.05,
                 data=wow)

otu_to_taxonomy(wow_da$significant_taxa,wow)
plot(wow_da)


data.table::fread(here("data/collapse/table-L4_8_8_22.txt")) %>% as_tibble()

###






plot(subtree)
tips  <- tibble(Species=subtree$tip.label,order= get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[1:41])  %>% mutate(Species = sub(".*dae_","",Species)) %>% mutate(Species = sub("_"," ",Species))

abund %>% rename(Taxa = V1) %>% 
  pivot_longer(-1,names_to = "Species",values_to = "Abundances") %>% 
  mutate(Abundances=ifelse(Abundances!=0,1,0)) %>% group_by(Taxa) %>% summarise(n=sum(Abundances)) %>% mutate(Prop = n/41*100) %>% arrange(desc(Prop)) %>% left_join(.,pers %>% rename(Taxa=V1),by="Taxa") %>% filter(Prop==100) %>% pivot_longer(-c(1:3),names_to = "Species",values_to = "Persistence") %>% 
  left_join(.,core_80,by="Taxa") %>% 
  relocate(Total,.after=Taxa) %>% filter(!grepl(";_",Taxa)) %>% left_join(.,tips,by="Species") -> core

colores = met.brewer("Hiroshige",n = core %>% distinct(Taxa) %>% nrow(),type="continuous")

#set.seed(32221) # for families
set.seed(33121) # for orders

core %>% left_join(.,abund %>% rename(Taxa = V1) %>% pivot_longer(-1,names_to = "Species",values_to = "Abundances"),by=c("Taxa","Species")) %>% group_by(Species) %>% summarise(n=sum(Abundances)) %>% summarize(range(n))


core %>% left_join(.,abund %>% rename(Taxa = V1) %>% pivot_longer(-1,names_to = "Species",values_to = "Abundances"),by=c("Taxa","Species")) %>% 
  ggplot(aes(fill=Taxa, y = Abundances, x = fct_reorder(Species,order))) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_discrete(type = sample(colores,length(colores),replace = F)) +
  #scale_fill_discrete(type = colores) +
  theme(legend.position = "",panel.background = element_blank(),panel.grid=element_blank(),
        axis.line = element_line()) + labs(x="",y="Relative Abundance") + 
  coord_flip() + 
  NULL
  
ggsave(filename = here("figures/core_relabund.pdf"),plot = last_plot())
  
  core %>% left_join(.,abund %>% rename(Taxa = V1) %>% pivot_longer(-1,names_to = "Species",values_to = "Abundances"),by=c("Taxa","Species")) %>% 
    ggplot(aes(fill=Taxa, y = Abundances, x = fct_reorder(Species,order))) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_discrete(type = sample(colores,length(colores),replace = F)) +
    #scale_fill_discrete(type = colores) +
    theme(legend.position = "left",panel.background = element_blank(),panel.grid=element_blank(),
          axis.line = element_line()) + labs(x="",y="Relative Abundance") + 
    coord_flip() + 
    NULL
  
  ggsave(filename = here("figures/core_relabund_legend.pdf"),plot = last_plot())
  
  }

plot_bioclims_locality <- function(data="data/Super_table_27_07_22_V2_60.csv") { 
  a <- data.table::fread(here(data)) %>% as_tibble()
  a <- a %>% mutate(Organism=as.factor(Organism),Family=as.factor(Family))
  a$Family <- factor(a$Family,levels=c("Hynobiidae","Cryptobranchidae","Salamandridae","Ambystomatidae","Plethodontidae"))
  a %<>% mutate(origin = sub("-caught","",origin)) %>% 
    mutate(origin = ifelse(origin=="wild",origin,"captive"))
  
p1 <- a %>% filter(origin=="wild") %>% distinct(lat,long,Family,.keep_all = TRUE) %>% 
   ggplot(aes(x=bio1,fill=bio1)) + 
    ggdist::geom_dots(side="right",color="NA",scale=1,size=0,dotsize=1,layout="bin",binwidth=0.5,shape=22,stackratio = 1.1) + scale_fill_stepsn(colors = wesanderson::wes_palette("Zissou1")[3:5],n.breaks=20) +
   theme(axis.title.x = element_text(size=14),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size=10),
         axis.ticks.y = element_blank(),
         axis.line.y = element_blank(),
         axis.line = element_line(),
         
         legend.position = "none",panel.background = element_blank(),panel.grid = element_blank()) +
  ggdist::stat_dotsinterval(dotsize=0)+
   NULL
  
 p2 <- a %>% filter(origin=="wild") %>% distinct(lat,long,Family,.keep_all = TRUE) %>% 
   ggplot(aes(x=bio12,fill=bio12)) + 
   ggdist::geom_dots(side="right",color="NA",scale=1,size=0.1,dotsize=1,layout="bin",binwidth=50,shape=22,stackratio = 1.1) + scale_fill_stepsn(colors=wesanderson::wes_palette("Darjeeling2")[c(4,2)],n.breaks=15) +
   theme(axis.title.x = element_text(size=14),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size=10),
         axis.ticks.y = element_blank(),
         axis.line.y = element_blank(),
         axis.line = element_line(),legend.position = "none",panel.background = element_blank(),panel.grid = element_blank()) +
   ggdist::stat_dotsinterval(dotsize=0)+
   NULL
 pp <- wrap_plots(list(p1,p2))
 ggsave(filename = here("figures","inset_maps.pdf"),plot =  pp)
 
}


mixed_model <- function(data = "data/Super_table_27_07_22_V2_60.csv",
                     to_drop = c("tm_min", "bio3", "bio4", "bio5","bio6","bio13","bio14"),
                     output="output_data/model_data.RData"){
  
  infossm <- read.table(here(data),sep=",",header=T) %>% as_tibble() %>% filter(origin!="captive") %>% filter(!is.na(elevation)) %>% mutate(across(16:38, ~scale(.x,scale=T,center=T),.names = "{.col}")) %>% 
    rename("Habitat"=morphotype) %>% mutate(whicha=faith_pd) #shannon_entropy
  infossm %>% select(c(16:38,47)) %>% drop_na() -> infossm_step 
  plot(infossm$whicha)
  infossm_step %>% colnames()
  #Step two
  #Deleting variables with pearson coefficient >= .7
  infossm %>% select(c(16:38,47)) %>% select(-all_of(to_drop)) %>% colnames() -> one
  one
  lm(data= infossm_step %>% select(-all_of(to_drop)), whicha ~ 1) -> FitStart
  lm(data= infossm_step %>% select(-all_of(to_drop)), whicha ~ .) -> FitAll
  step(FitStart,direction = "both",scope = formula(FitAll)) -> best_step
  one %>% .[which(one %in% names(best_step$coefficients))] -> two
  infossm %>% select(all_of(two),"Habitat","Family") %>% colnames() %>% c(.,"(1|Study)") %>% reformulate(.,response="whicha") -> lmm_formula
  infossm %>% drop_na() %>% mutate(Habitat = as.numeric(as.factor(Habitat))-1) -> infossm_lmm
  lmm2.1_c <- lme4::lmer(formula = lmm_formula, REML=F, data = infossm_lmm)
  #Cheking collinearity and deleting variables with VIF > 11 
  performance::check_collinearity(lmm2.1_c) -> vifs
  if(any(vifs$VIF>=10)) { 
    tibble(term=vifs$Term,VIF=vifs$VIF,SE=vifs$SE_factor) %>% filter(VIF == max(VIF)) %>% pull(term) -> high_vif
    terms(lmm_formula) -> terminos
    drop.terms(terminos,dropx = which(attr(terminos,"term.labels") %in% high_vif),keep.response = F) -> simplified
    attr(simplified,"term.labels") %>% sub(" + 1 | Study","",.,fixed = T) %>% paste(.," + (1|Study)",sep="") %>% reformulate(.,response="whicha") -> simplified
    lmm2.1_c <- lme4::lmer(formula = simplified, REML=F, data = infossm_lmm)
    performance::check_collinearity(lmm2.1_c) -> vifs
    if(any(vifs$VIF>=10)) {
      tibble(term=vifs$Term,VIF=vifs$VIF,SE=vifs$SE_factor) %>% filter(VIF == max(VIF)) %>% pull(term) -> high_vif
      terms(simplified) -> terminos
      drop.terms(terminos,dropx = which(attr(terminos,"term.labels")%in%high_vif),keep.response = F) -> simplified
      attr(simplified,"term.labels") %>% sub(" + 1 | Study","",.,fixed = T) %>% paste(.," + (1|Study)",sep="") %>% reformulate(.,response="whicha") -> simplified
      lmm2.1_c <- lme4::lmer(formula = simplified, REML=F, data = infossm_lmm)
      performance::check_collinearity(lmm2.1_c) -> vifs
      if(any(vifs$VIF>=10)) {
        tibble(term=vifs$Term,VIF=vifs$VIF,SE=vifs$SE_factor) %>% filter(VIF == max(VIF)) %>% pull(term) -> high_vif
        terms(simplified) -> terminos
        drop.terms(terminos,dropx = which(attr(terminos,"term.labels")%in%high_vif),keep.response = F) -> simplified
        attr(simplified,"term.labels") %>% sub(" + 1 | Study","",.,fixed = T) %>% paste(.," + (1|Study)",sep="") %>% reformulate(.,response="whicha") -> simplified
        lmm2.1_c <- lme4::lmer(formula = simplified, REML=F, data = infossm_lmm)
        performance::check_collinearity(lmm2.1_c)
      }
    }
  }
  summary(lmm2.1_c)
  MuMIn::r.squaredGLMM(lmm2.1_c)
  results <- list(data=infossm_lmm,model_s=lmm2.1_c,
                  R_squared_c=MuMIn::r.squaredGLMM(lmm2.1_c))
  
  #save(results, file = here("output_data",output))
  
  #broom.mixed::tidy(lmm2.1_c,conf.int=T) %>% mutate(across(is.numeric,round,2)) %>% flextable::flextable() %>% flextable::save_as_docx(path = here("figures/table_S2.docx"))
  lmm_fix <- broom.mixed::tidy(lmm2.1_c,conf.int=T) %>% slice(-1) %>% filter(!effect=="ran_pars") %>% mutate(nonzero = !(conf.low <=0 & conf.high >= 0)) %>% 
    mutate(sign = ifelse(estimate < 0,1,2)) %>% mutate(ok = nonzero*sign) %>% 
    arrange(estimate) %>% 
    ggplot(aes(x=fct_reorder(term,estimate),y=estimate)) +
    geom_hline(yintercept = 0,color="black",linetype="dashed") +
    geom_segment(aes(x=fct_reorder(term,estimate),xend=fct_reorder(term,estimate),y=conf.low,yend=conf.high,color=factor(sign,levels=c(1,0))),size=1.5) + 
    theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line(),legend.position = "bottom") + 
    labs(y="Fixed effect estimate",x="") +
    geom_point(aes(fill=factor(sign,levels=c(1,0))),size=3,shape=21) + coord_flip() +
    scale_color_manual(values= wesanderson::wes_palette("Moonrise3")[c(1:2)],name="fixed effects",labels=c("Negative","Positive"))+
    scale_fill_manual(values= wesanderson::wes_palette("Moonrise3")[c(1:2)],name="fixed effects",labels=c("Negative","Positive")) +
    NULL
  lmm_fix
  
  
  ggsave(filename = "figures/lmm_fixed.pdf", plot =  lmm_fix)

  data="data/Super_table_27_07_22_V2_60.csv"
  a <- data.table::fread(here(data)) %>% as_tibble()
  a <- a %>% mutate(Organism=as.factor(Organism),Family=as.factor(Family))
  a$Family <- factor(a$Family,levels=c("Hynobiidae","Cryptobranchidae","Salamandridae","Ambystomatidae","Plethodontidae"))
  a %<>% mutate(origin = sub("-caught","",origin)) %>% 
    mutate(origin = ifelse(origin=="wild",origin,"captive"))
  a %>% mutate(across(16:38, ~scale(.x,scale=T,center=T),.names = "{.col}")) %>% #group_by(Organism) %>% #summarise(across(.cols = c(15:37,42),mean,na.rm=T),Habitat=first(Habitat),Family=first(Family)) %>% 
    ggplot(aes(x=pre,y=shannon_entropy)) + geom_point(shape=21,size=2,alpha=0.6) + 
    theme(axis.title.x = element_text(size=14),axis.text.y = element_text(size=14),axis.text.x = element_text(size=10),axis.line = element_line(),legend.position = "none",panel.background = element_blank(),panel.grid = element_blank()) + geom_smooth(color="grey50",method = "lm",se = T) +
    #scale_fill_fermenter(n.breaks=9,na.value = "black",direction = 1,type = "seq",palette="GnBu")+
    labs(x="Monthly precipitation (mm)",y="Shannon's Entropy") + 
    NULL
  
  aver <- infossm_lmm %>% bind_cols(.,fit=predict(lmm2.1_s) )
  cofs <- coef(lmm2.1_s)$Study %>% as_tibble(rownames = "Study") %>%
    rename(Intercept=`(Intercept)`) 
  names(cofs)[-c(1:2)] <- paste0("B_",names(cofs)[-c(1:2)])
  
  aver %>% right_join(cofs,.,by="Study") %>% group_by(Study) %>% distinct(lat,long) %>% count(Study) %>% filter(n>3) %>% mutate(ok=TRUE) -> targets
  aver %>% right_join(cofs,.,by="Study") %>% 
    #right_join(targets,.,by="Study") %>% filter(ok) %>% 
  ggplot(aes(x=bio17,y=shannon_entropy,fill=bio17,group=Study)) + 
    geom_point(shape=21,size=3,stroke=0.2,alpha=0.7) + 
    theme(axis.title.x = element_text(size=14),axis.text.y = element_text(size=14),axis.text.x = element_text(size=10),axis.line = element_line(),legend.position = "none",panel.background = element_blank(),panel.grid = element_blank()) + 
    #geom_abline(aes(slope = B_pre,intercept = Intercept)) +
    geom_smooth(inherit.aes = FALSE,
                aes(x=bio17,y=shannon_entropy), method="lm",se = FALSE,color="black",size=3) +
    geom_smooth(method="lm",se = FALSE,color="black",size=0.5) +
    scale_fill_gradient(low=wesanderson::wes_palette("Darjeeling2")[c(4)],
                        high=wesanderson::wes_palette("Darjeeling2")[c(2)])  +
    NULL
  
  
  
}


mantel_multi_phylo <- function(data = "data/Super_table_27_07_22_V2_60.csv",
                               phylo="data/trees/treePL_RAxML_bootstrap.Caudata_v6.trees"){ 
  geos <- data.table::fread(data)
  geos2 <- geos
  geos <- geos %>% as_tibble() %>% #filter(origin!="captive") %>% 
    filter(Organism!="Echinotriton andersoni") %>% 
    group_by(Organism) %>% mutate(Organism=sub(" ","_",Organism)) %>% 
    summarise(across(where(is.numeric),mean,na.rm=T),Family=first(Family))
  
  tree_BP <- ape::read.tree(phylo)
  tree_BP  <- tree_BP[sample(1:length(tree_BP),100,replace=FALSE)]
  tree_BP[[length(tree_BP)+1]] <- ape::read.tree("data/trees/treePL_RAxML_bestTree.Caudata_v6.tre")
  
  level = list.files(path=here("data/bc_species"),pattern = "Specie_bc_vegan")
  nivel = stringr::str_sub(level,-9,-5)
  all_plot <- list()
  for(j in 1:length(level)){ 
    cat("Performing mantels for level",j,"\n")
    r.mantel <- list()
    sig.mantel <- list()
    r.mantel.partial <- list()
    sig.mantel.partial <- list()
    mantel.correlogram <- list()
    mantel.correlogram.partial <- list()
    which_tree <- list()
    for (i in 1:length(tree_BP)){
      cat("Tree:",i,"\n")
      quien = tree_BP[[i]]
      mat_dist <- adephylo::distTips(quien,tips = "all",method = "patristic")
      tipas <- unlist(lapply(lapply(strsplit(quien$tip.label,"_"),"[",2:3),paste,collapse="_"))
      stopifnot(length(geos$Organism[which(is.na(match(geos$Organism,tipas)))])==0) 
      fin_tip <- tipas %>% match(geos$Organism,.) %>% .[!duplicated(.)]
      mat_dist <- as.matrix(mat_dist)
      mat_dist[fin_tip,fin_tip] -> mat_dist
      micro <- read.table(here("data/bc_species/",level[j]),sep=",",header=T)
      colnames(micro) <- sub("\\.","_",colnames(micro))
      micro <- micro[,-1]
      row.names(mat_dist) %>% strsplit(.,"_") %>% lapply(.,"[",2:3) %>% lapply(.,paste,collapse="_") %>% unlist () -> wa
      which(is.na(match(colnames(micro),wa))) -> captive
      if(length(captive!=0)) micro[-captive,-captive] -> micro
      stopifnot(length(which(is.na(match(colnames(micro),wa))))==0)
      micro[match(wa,colnames(micro)),match(wa,colnames(micro))] -> micro
      rownames(micro) <- colnames(micro)
      micro <- as.matrix(micro)
      diag(mat_dist) <- NA
      diag(micro)<- NA
      
      rownames(mat_dist) <- rownames(mat_dist) %>% strsplit(.,"_") %>% lapply(.,"[",2:3) %>% lapply(.,paste,collapse="_") %>% unlist () 
      colnames(mat_dist) <- rownames(mat_dist)
      if(i==length(tree_BP)) {save(mat_dist,file=here("output_data/patristic_distances.RData"))
        save(micro,file=here("output_data/micro_distances.RData"))
        }
      a <- vegan::mantel(mat_dist,micro,method = "pearson",permutation=999)
      r.mantel[[i]] <- a$statistic
      sig.mantel[[i]] <- a$signif
      bios_mat <- geos2 %>% as_tibble() %>%  filter(Organism!="Echinotriton andersoni") %>%  
        distinct(long,lat,Organism,.keep_all = T) %>% group_by(Organism) %>% dplyr::select(starts_with("bio",ignore.case = F))
      bios_mat %>% ungroup() %>% filter(!is.na(bio1)) %>% data.frame(.) -> bios_df
      bios_pca <- ade4::dudi.pca(bios_df[,-1],scannf = FALSE, nf = 5)
      tibble(Species=bios_df[,1],bios_pca$li) %>% group_by(Species) %>% summarise(across(starts_with("Axis"), list(mean=mean))) -> pcas_bios
      tibble(Species=bios_df[,1],bios_pca$li) %>% ggplot(aes(x=Axis1,y=Axis2,color=Species,fill=Species)) + geom_point() + theme(legend.position = "none") + 
        geom_point(data=pcas_bios,aes(x=Axis1_mean,y=Axis2_mean,fill=Species),size=5,shape=21,color="black")
      pcas_bios %>% ungroup() %>% dplyr::select(-1) %>% dist(.) %>% as.matrix() -> dist_bio
      colnames(dist_bio) <- rownames(dist_bio) <- sub(" ","_",pcas_bios$Species)
      stopifnot(length(which(is.na(match(colnames(mat_dist),colnames(dist_bio)))))==0)
      dist_bio <- dist_bio[match(colnames(mat_dist),colnames(dist_bio)),match(colnames(mat_dist),colnames(dist_bio))]
      x_mat <- mat_dist / 2
      y_mat <- micro
      z_mat <- dist_bio
      test.mantel <- vegan::mantel.partial(x_mat,y_mat,z_mat)
      r.mantel.partial[[i]] <- test.mantel$statistic
      sig.mantel.partial[[i]] <- test.mantel$signif
      
      mantel.correlog_2(y_mat,x_mat,cutoff = F,mult="fdr",r.type="pearson",nperm = 999,progressive = T) -> x.mant
      mantel.correlogram[[i]] <- x.mant$mantel.res
      mpmcorrelogram_2(y_mat,x_mat,zdis = z_mat,method="pearson",permutation=999,breaks=x.mant$break.pts,plot=F,print=F,method_adj = "fdr") -> y.mant
      mantel.correlogram.partial[[i]] <- y.mant
      which_tree[[i]] <- i
      
    }
all_plot[[j]] <- tibble(level=nivel[j],tree=which_tree %>% unlist(),r.mantel = r.mantel %>% unlist(),
           sig.mantel = sig.mantel%>% unlist(),
           r.mantel.partial = r.mantel.partial %>% unlist(),
           sig.mantel.partial = sig.mantel.partial %>% unlist(),
           mantel.correlogram = mantel.correlogram,
           mantel.correlogram.partial = mantel.correlogram.partial
    ) 
  }
  res.mantels <- all_plot %>% do.call("rbind",.)
  save(res.mantels,file = "output_data/Mantel_results.Rdata")
}


plot_mantels <- function(data="Mantel_results.Rdata",
                         distance="patristic_distances.RData",
                         micro="micro_distances.RData"){
  load(here("output_data",data))
  load(here("output_data",distance))
  load(here("output_data",micro))
  quien = "lvl_4"
  res.mantels %<>% filter(level==quien) 
    res.mantels$r.mantel %>% quantile(probs=c(0.01,0.5,0.99)) -> aqui
    statistic <- paste("rM (BP trees) = " ,round(aqui[1],3)," \u2013 ",round(aqui[3],3),sep="")
    ML.stat <- paste("rM (ML tree) = ",round(res.mantels$r.mantel[nrow(res.mantels)],3))
    
    res.mantels$r.mantel.partial %>% quantile(probs=c(0.01,0.5,0.99)) -> aqui
    statistic.p <- paste("p-rM (BP trees) = " ,round(aqui[1],3)," \u2013 ",round(aqui[3],3),sep="")
    ML.stat.p <- paste("p-rM (ML tree) = ",round(res.mantels$r.mantel.partial[nrow(res.mantels)],3))
    
tibble(Phylo=as.vector(as.dist(mat_dist))/2,Microbiome=as.vector(as.dist(micro))) %>% ggplot(aes(x=Phylo,y=Microbiome)) + geom_jitter(height=0,width=1.5,shape=21,stroke=0.01,size=2,color="black",fill="grey30",alpha=0.8) + geom_smooth(method=lm, fill="grey50",color="grey10", se=TRUE)  +
      theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line.x = element_line(),axis.line.y = element_line()) + labs(x="Phylogenetic distance (Mya)",y="Microbiome distance (Bray Curtis)",title="Host phylogeny vs. microbiome",subtitle=quien) + 
      annotate("text",x=130,y=0.105,label=statistic) +
      annotate("text",x=130,y=0.15,label=ML.stat) +
      annotate("text",x=130,y=0.0135,label=statistic.p) +
      annotate("text",x=130,y=0.046,label=ML.stat.p) +
      NULL
nombre = paste0("Mantel_plot_",quien,".pdf")
ggsave(file = here("figures",nombre),last_plot())
    
mantel_plot <- res.mantels %>% mutate(mantel.correlogram=(lapply(res.mantels$mantel.correlogram, function(x){ x %>% as_tibble(rownames = "Dist.Class") %>% mutate(Dist.Class=sub("D.cl.","",Dist.Class))})))
    
    mantel_plot_partial <- res.mantels %>% 
      mutate(mantel.correlogram.partial=(lapply(res.mantels$mantel.correlogram.partial, function(x){ x$clases %>% strsplit(.,split="  -  ") -> xx
        x$clases <- xx %>% lapply(.,as.numeric) %>% lapply(.,mean) %>% unlist()
        tibble(Dist.Class=1:11,class.index=x$clases,Mantel.cor=x$rM,`Pr(Mantel)`=x$pvalues,
               `Pr(corrected)`=x$pval.Bonferroni) })))
    
    mantel_plot %>% slice(nrow(res.mantels)) %>% unnest(mantel.correlogram) %>% 
      mutate(alphas = case_when(`Pr(Mantel)` < 0.05 ~0.8, `Pr(Mantel)` >= 0.05 ~ 0.2)) %>% 
      mutate(alphas=case_when(!is.na(`Pr(Mantel)`) ~ alphas ,is.na(`Pr(Mantel)`) ~ 0)) %>% 
      mutate(cats=cut(`Pr(Mantel)`,breaks=c(0,0.01,0.025,0.05,0.06,0.08,0.1,0.3,0.5,0.8,1.1),include.lowest = T,labels=c(0.01,0.025,0.05,0.06,0.08,0.1,0.3,0.5,0.8,1))) -> ML
    
    mantel_plot %>% unnest(mantel.correlogram)  %>% 
      mutate(alphas = case_when(`Pr(Mantel)` < 0.05 ~0.8, `Pr(Mantel)` >= 0.05 ~ 0.2)) %>% 
      mutate(alphas=case_when(!is.na(`Pr(Mantel)`) ~ alphas ,is.na(`Pr(Mantel)`) ~ 0)) %>% 
      #mutate(`Pr(corrected)` = case_when(is.na(`Pr(corrected)`) ~ 1,!is.na(`Pr(corrected)`) ~ `Pr(corrected)`)) %>%   
      mutate(cats=cut(`Pr(Mantel)`,breaks=c(0,0.01,0.025,0.05,0.06,0.08,0.1,0.3,0.5,0.8,1.1),include.lowest = T,labels=c(0.01,0.025,0.05,0.06,0.08,0.1,0.3,0.5,0.8,1))) -> mantel_plot
    mantel_plot %>% group_by(Dist.Class) %>% summarize(mean(class.index))
    pallete <- c(MetBrewer::met.brewer("Demuth",15,type="continuous")[c(3,5,6,8,9,10,12,14)])
  ML <- mantel_plot %>% 
    group_by(level,Dist.Class) %>% 
    filter(row_number() ==  n()) %>% 
    mutate(alphas = case_when(`Pr(Mantel)` < 0.05 ~ 1.0, `Pr(Mantel)` >= 0.05 ~ 0)) %>% 
    mutate(alphas=case_when(!is.na(`Pr(Mantel)`) ~ alphas ,is.na(`Pr(Mantel)`) ~ 0)) %>% 
    mutate(cats=cut(`Pr(Mantel)`,breaks=c(0,0.01,0.025,0.05,0.06,0.08,0.1,0.3,0.5,0.8,1.1),include.lowest = T,labels=c(0.01,0.025,0.05,0.06,0.08,0.1,0.3,0.5,0.8,1))) #%>% filter(!Level%in%c("lvl_2","lvl_7"))
  
  set.seed(1234)
  mantel_plot %>% filter(!is.na(Mantel.cor)) %>% filter(!Dist.Class %in% c(10,11)) %>% 
    ggplot(aes(x=class.index,y=Mantel.cor,fill=`Pr(corrected)`,color=`Pr(corrected)`)) +
    geom_hline(yintercept = 0,alpha=0.6) +
    geom_line(data=  ML %>% filter(!is.na(Mantel.cor)) %>% filter(!Dist.Class %in% c(10,11)),aes(x=class.index,y=Mantel.cor),size=0.9,alpha=0.9,color="black") +
    geom_jitter(shape=21,height=0.001,width=2,alpha=0.8,stroke=0,size=2) + 
    geom_point(data=ML %>% filter(!is.na(Mantel.cor)) %>% filter(!Dist.Class %in% c(10,11)),aes(x=class.index,y=Mantel.cor,group=Dist.Class,fill=`Pr(corrected)`),shape=21,color="black",size=3,stroke=1,inherit.aes = F) +
    scale_fill_stepsn(colors=MetBrewer::met.brewer("Demuth",direction=1),breaks=c(0.01,0.05,0.1,1.0),name="p-value",limits=c(0,0.2)) +
    scale_color_stepsn(colors=MetBrewer::met.brewer("Demuth",direction=1),breaks=c(0.01,0.05,0.1,1.0),name="p-value",limits=c(0,0.2)) +
    coord_cartesian(ylim=c(-0.25,0.25)) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          plot.subtitle=element_text(size=10, hjust=0, color="black"),
          plot.title=element_text(size=14, hjust=0, face="bold", color="black")) + 
    labs(x="Host phylogenetic distance (Mya)",y="Mantel correlation coefficient",subtitle= paste0("Mantel Correlogram of phylogenetic distance vs. Bray Curtis dissimilarity ","\U0028",quien,"\U0029"),title="The influence of host phylogeny on microbiome beta diversity") +
    scale_alpha(range = c(0.5,1)) +
    scale_size(range=c(2,3)) +
    guides(alpha="none",size="none") +
    NULL
  
  nombre = paste0("Mantel_Corre_",quien,".pdf")
  ggsave(file = here("figures",nombre),last_plot())

  }  
  

core_abundances <- function(core="core_4.csv",
                      persistence = "Persistencia_lvl_4.csv",
                      abundance="Abundancia_lvl_4.csv",
                      presence="Persistencia_aparicion_lvl_4.csv",
                      tree="data/trees/treePL_RAxML_bestTree.Caudata_v6.tre",
                      data="data/Super_table_27_07_22_V2_60.csv"){ 
  data.table::fread(here("data/core",core)) %>% as_tibble()  -> pam
  data.table::fread(here("data/core",persistence)) %>% as_tibble()  -> pers
  data.table::fread(here("data/core",abundance)) %>% as_tibble()  -> abund

geos_sample <- data.table::fread(here("data/collapse/table-L4_8_8_22.txt")) %>% as_tibble()
all_data <- data.table::fread(data)%>% as_tibble()
tree <- ape::read.tree(tree)
geos <- all_data %>% group_by(Organism) %>% mutate(V1=as.factor(V1)) %>% 
    mutate(Organism=sub(" ","_",Organism)) %>% summarise(across(where(is.numeric),mean,na.rm=T),Family=first(Family))
  tipas <- unlist(lapply(lapply(strsplit(tree$tip.label,"_"),"[",2:3),paste,collapse="_"))
  stopifnot("Error: tree and database not fully matched!"=geos$Organism[which(is.na(match(geos$Organism,tipas)))]==0)
  fin_tip <- tipas %>% match(geos$Organism,.) %>% .[!duplicated(.)]
  subtree <- keep.tip(tree,tree$tip.label[fin_tip])
  subtree <- ladderize(subtree,FALSE)
  subtree$root.time <- max(branching.times(subtree),na.rm=T)
  
targets <- abund %>% rename(Taxa = V1) %>% 
    pivot_longer(-1,names_to = "Species",values_to = "Abundances") %>% 
    mutate(Abundances=ifelse(Abundances!=0,1,0)) %>% group_by(Taxa) %>% summarise(n=sum(Abundances)) %>% mutate(Prop = n/41*100) %>% arrange(desc(Prop)) %>% left_join(.,pers%>% rename(Taxa = V1),by="Taxa") %>% filter(Prop==100) %>% 
    pivot_longer(-c(1:3),names_to = "Salamander",values_to = "Persistence") %>% 
    filter(!grepl(";_",Taxa)) %>% distinct(Taxa) %>% pull()

  
all_data %<>% group_by(Organism) %>% mutate(V1=as.factor(V1)) %>% 
  mutate(Organism=sub(" ","_",Organism))
geos_sample %<>% pivot_longer(-1,names_to = "Sample_id",values_to = "Raw_reads") %>% left_join(.,all_data,by="Sample_id") %>% filter(!is.na(V1)) %>% filter(`#OTU ID` %in% all_of(targets)) %>% rename(Taxa=`#OTU ID`,Salamander=Organism)

library(phyloseq)
library(corncob)

geos_sample %>% filter(!is.na(tm_max)) %>% 
  select(Taxa,Sample_id,Raw_reads) %>% pivot_wider(names_from = Sample_id,values_from = Raw_reads) %>%
  select(-1) %>% as.matrix(.) %>% otu_table(.,taxa_are_rows = T) -> otus

geos_sample %>% filter(!is.na(tm_max)) %>% select(Taxa,Sample_id,Raw_reads) %>% pivot_wider(names_from = Sample_id,values_from = Raw_reads) %>% select(1) %>% as.matrix(.) %>% tax_table(.) -> taxa

geos_sample %>% filter(!is.na(tm_max)) %>% select(Sample_id,Salamander,tm_max:bio19) %>% distinct(Sample_id,.keep_all = TRUE) %>% as.data.frame(.) -> samples
row.names(samples) <- samples$Sample_id
samples <- sample_data(samples)
wow <- merge_phyloseq(otus,taxa,samples)
wow_da <- differentialTest(formula= ~ tm_max, phi.formula = ~ tm_max,
                           formula_null = ~ 1,
                           phi.formula_null = ~ tm_max,
                           test="LRT",fdr_cutoff = 0.05,
                           data=wow)

otu_to_taxonomy(wow_da$significant_taxa,wow)
plot(wow_da)


geos_sample %>% 
  select(Taxa,Sample_id,Raw_reads,lat,long,Salamander) %>% group_by(Taxa,Salamander,lat,long) %>% 
  summarise(Raw_reads=sum(Raw_reads),Sample_id=first(Sample_id)) %>% ungroup() %>%  select(-lat,-long,-Salamander) %>% pivot_wider(names_from = Sample_id,values_from = Raw_reads) %>%
  select(-1) %>% as.matrix(.) %>% otu_table(.,taxa_are_rows = T) -> otus

geos_sample %>%  select(Taxa,Sample_id,Raw_reads,lat,long,Salamander) %>% group_by(Taxa,Salamander,lat,long) %>% 
  summarise(Raw_reads=sum(Raw_reads),Sample_id=first(Sample_id)) %>% ungroup() %>%  select(-lat,-long,-Salamander) %>% pivot_wider(names_from = Sample_id,values_from = Raw_reads) %>%select(1) %>% as.matrix(.) %>% tax_table(.) -> taxa

geos_sample %>% 
  select(Taxa,Sample_id,Raw_reads,lat,long,Salamander) %>% group_by(Taxa,Salamander,lat,long) %>% 
  summarise(Raw_reads=sum(Raw_reads),Sample_id=first(Sample_id)) %>% ungroup() %>%  select(-lat,-long,-Salamander) %>% pivot_wider(names_from = Sample_id,values_from = Raw_reads) -> target_samples

geos_sample %>% select(Sample_id,Salamander,tm_max:bio19) %>% distinct(Sample_id,.keep_all = TRUE) %>% 
  filter(Sample_id %in% all_of(colnames(target_samples)[-1])) %>% as.data.frame(.) -> samples
row.names(samples) <- samples$Sample_id





geos %>% mutate(Salamander=sub("_"," ",Organism)) %>% left_join(dosa,.,by="Salamander") %>% 
    filter(Taxa=="Bacteria;Proteobacteria;Alphaproteobacteria;Caulobacterales") %>% 
    ggplot(aes(x=pre,y=Abundances,color=Taxa)) +
    geom_point() +
    geom_smooth(method="lm",se=FALSE)+
    #geom_smooth(method="gam",method.args=list(family="ziplss"),se = FALSE) +
    theme(axis.title.x = element_text(size=12),
          axis.text.y = element_blank(),
          axis.line = element_line(),
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.key.width = unit(0.7,"in")) + 
    scale_fill_stepsn(colors=MetBrewer::met.brewer("Hiroshige"),n.breaks=10,name="Prevalence across\nall samples") +
    NULL
  
  geos %>%mutate(Salamander=sub("_"," ",Organism)) %>% left_join(dosa,.,by="Salamander") %>% distinct(Taxa)
  geos %>%mutate(Salamander=sub("_"," ",Organism)) %>% left_join(dosa,.,by="Salamander") %>% 
    #filter(Taxa=="Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales") %>% 
    mutate(Depth = 2500) %>% filter(!is.na(bio17))%>% mutate(Abundances= round(Abundances*2500/100)) -> test
 
  

  ###
  
  
  
  
  
  
  plot(subtree)
  tips  <- tibble(Species=subtree$tip.label,order= get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[1:41])  %>% mutate(Species = sub(".*dae_","",Species)) %>% mutate(Species = sub("_"," ",Species))
  
  abund %>% rename(Taxa = V1) %>% 
    pivot_longer(-1,names_to = "Species",values_to = "Abundances") %>% 
    mutate(Abundances=ifelse(Abundances!=0,1,0)) %>% group_by(Taxa) %>% summarise(n=sum(Abundances)) %>% mutate(Prop = n/41*100) %>% arrange(desc(Prop)) %>% left_join(.,pers %>% rename(Taxa=V1),by="Taxa") %>% filter(Prop==100) %>% pivot_longer(-c(1:3),names_to = "Species",values_to = "Persistence") %>% 
    left_join(.,core_80,by="Taxa") %>% 
    relocate(Total,.after=Taxa) %>% filter(!grepl(";_",Taxa)) %>% left_join(.,tips,by="Species") -> core
  
  colores = met.brewer("Hiroshige",n = core %>% distinct(Taxa) %>% nrow(),type="continuous")
  
  #set.seed(32221) # for families
  set.seed(33121) # for orders
  
  core %>% left_join(.,abund %>% rename(Taxa = V1) %>% pivot_longer(-1,names_to = "Species",values_to = "Abundances"),by=c("Taxa","Species")) %>% group_by(Species) %>% summarise(n=sum(Abundances)) %>% summarize(range(n))
  
  
  core %>% left_join(.,abund %>% rename(Taxa = V1) %>% pivot_longer(-1,names_to = "Species",values_to = "Abundances"),by=c("Taxa","Species")) %>% 
    ggplot(aes(fill=Taxa, y = Abundances, x = fct_reorder(Species,order))) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_discrete(type = sample(colores,length(colores),replace = F)) +
    #scale_fill_discrete(type = colores) +
    theme(legend.position = "",panel.background = element_blank(),panel.grid=element_blank(),
          axis.line = element_line()) + labs(x="",y="Relative Abundance") + 
    coord_flip() + 
    NULL
  
  ggsave(filename = here("figures/core_relabund.pdf"),plot = last_plot())
  
  core %>% left_join(.,abund %>% rename(Taxa = V1) %>% pivot_longer(-1,names_to = "Species",values_to = "Abundances"),by=c("Taxa","Species")) %>% 
    ggplot(aes(fill=Taxa, y = Abundances, x = fct_reorder(Species,order))) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_discrete(type = sample(colores,length(colores),replace = F)) +
    #scale_fill_discrete(type = colores) +
    theme(legend.position = "left",panel.background = element_blank(),panel.grid=element_blank(),
          axis.line = element_line()) + labs(x="",y="Relative Abundance") + 
    coord_flip() + 
    NULL
  
  ggsave(filename = here("figures/core_relabund_legend.pdf"),plot = last_plot())
  
}


  #####

as_tibble(nmds$points) %>% mutate(Sample_id=rownames(nmds$points),.before=1) %>% inner_join(.,geos,by="Sample_id") -> geos 

geos %>% ggplot(aes(x=MDS1,y=MDS2,fill=Family,color=Family)) + geom_point(shape=21,alpha=0.8) + theme(panel.background = element_blank(),panel.grid = element_blank(),legend.position = "right",axis.line = element_line()) + #stat_ellipse() +
  scale_fill_manual(values=c(wesanderson::wes_palette("Darjeeling1"))) + scale_color_manual(values=c(wesanderson::wes_palette("Darjeeling1"))) + 
  annotate(geom="text", x=0.6, y=-.5, label=paste("Stress:",round(nmds$stress,3),sep=" "))

geos %>% ggplot(aes(x=MDS1,y=MDS2,fill=morphotype,color=morphotype)) + geom_point(shape=21,alpha=0.8) + theme(panel.background = element_blank(),panel.grid = element_blank(),legend.position = "right",axis.line = element_line()) + #stat_ellipse() +
  scale_fill_manual(values=c(wesanderson::wes_palette("Darjeeling2")[c(2,1)])) +scale_color_manual(values=c(wesanderson::wes_palette("Darjeeling2")[c(2,1)])) + annotate(geom="text", x=0.6, y=-.5, label=paste("Stress:",round(nmds$stress,3),sep=" "))

p2
p3 <- geos %>% ggplot(aes(x=MDS1,y=MDS2,fill=bio17)) + geom_point(shape=21,alpha=0.7,size=2,stroke=0.1) + theme(panel.background = element_blank(),panel.grid = element_blank(),legend.position = "right",axis.line = element_line()) + scale_fill_fermenter(n.breaks=9,na.value = "black",direction = 1,type = "seq",palette="GnBu")
#pdf("figure_1.pdf")
p1
p2
p3
#dev.off()

do_map_deprecated <- function(input.tre = "data/RAxML_Caudata_v5.tre",data="data/Super_table_v2.csv") { 
  geos <- data.table::fread(here(data)) %>% group_by(long,lat,Organism) %>% mutate(Organism=sub(" ","_",Organism))
  a <- ape::read.tree(here(input.tre)) %>% ladderize(.,FALSE)
  a_macro <- sub(".*?_","",a$tip.label)
  a_macro[grep("sp$",a_macro)] <- "Pseudoeurycea_sp"
  a_micro <- unique(geos$Organism) %>% str_squish(.)
  fams <- a$tip.label[match(a_micro,a_macro)] %>% strsplit(.,"_") %>% lapply(.,"[",1) %>% unlist(.) %>% cbind(.,a_micro)
  geos <- geos %>% summarise(.,n=n())
  geos <- geos %>% ungroup() %>% mutate(Family=fams[match(sub(" ","_",geos$Organism),fams[,2]),1])
  families <- unique(geos$Family)
  val_limits <- range(geos$n)
  plotas <- list()
  for (sc in 1:length(families)){
    plotas[[sc]] <- geos %>% filter(Family==families[sc]) %>%  
      ggplot() + #xlim(xlims) + ylim(ylims) +
      geom_sf(data=rnaturalearth::ne_countries(scale=110,type="countries",returnclass = "sf"),colour="grey80",fill="grey80")+
      geom_point(aes(y=jitter(lat,amount = 3), x=jitter(long,amount=3),colour = n)) +
      scale_colour_viridis(option = "viridis",direction=1,begin=0.2,end=0.95,discrete = F,limits=val_limits) +
      theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
            panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
            axis.text.x = element_blank(),axis.text.y = element_blank(),
            panel.background=element_rect(colour="black",fill="azure2"),panel.border=element_rect(colour="black",fill=NA)) +
      #labs(subtitle=families[sc]) +
      NULL
  }
  wrap_plots(plotas) +
    plot_layout(ncol = 2,byrow=T,guides="collect") +
    plot_annotation(title ="Microbiome data",
                    caption='Rebollar et al.') & theme(legend.position = 'bottom')
}

dated_phylo_data <- function(output = "MS/Master_Figure_S2.pdf",input.tre = "data/treePL/treePL_Caudata.v.3.tre",data= "data/Super_table_ss_10_06_21.csv") {
tree <- ape::read.tree(here(input.tre))
geos <- data.table::fread(here(data)) %>% group_by(long,lat,Organism)
tipas <- unlist(lapply(lapply(strsplit(tree$tip.label,"_"),"[",2:3),paste,collapse="_"))
geos$Organism <- sub(" ","_",geos$Organism)
geos$Organism[grep("Mozotal",geos$Organism)] <- "Pseudoeurycea_sp"
geos$Organism[grep("Cryptobranchus_alleganiensis",geos$Organism)] <- "Cryptobranchus_alleganiensis"
geos$Organism[which(is.na(match(geos$Organism,tipas)))]
fin_tip <- tipas %>% match(geos$Organism,.) %>% .[!duplicated(.)]
subtree <- keep.tip(tree,tree$tip.label[fin_tip])
mat_branch <- cophenetic(subtree)
diag(mat_branch) <- NA
tree <- ape::ladderize(tree)
lapply(strsplit(tree$tip.label,"_"),"[",1) %>% unlist() %>% as_tibble() -> familias
familias %>% rename("Familia"=value) %>% mutate(Colores=as_factor(Familia))-> familias
familias %>% mutate(Familia=as_factor(Familia))-> familias
levels(familias$Colores) <- viridis::viridis(length(levels(familias$Colores)),option = "B",direction = -1)
familias <- familias %>% mutate(Tips=tree$tip.label)
familias <- familias %>% arrange(factor(Familia, levels = levels(familias$Familia)[c(11,3,2,1,6,5,4,7,8,10,9)]))
nombres <- levels(familias$Familia)[c(11,3,2,1,6,5,4,7,8,10,9)]
noda <- c()
for (i in 1:length(nombres)) {
  c(noda,getMRCA(tree,grep(nombres[i],tree$tip.label))) -> noda
}

col_groups <- c("black",viridis::viridis(length(levels(familias$Familia)),option = "B",direction = -1))
g <- tree %>% groupClade(.,noda) %>% ggtree(layout="circular") + NULL
familias %>% dplyr::select(Familia) %>%  as.data.frame() -> famis
rownames(famis) <- familias$Tips
famis$Familia <- as.factor(famis$Familia)
gg <- gheatmap(g, famis, offset=-.1, width=0.1,colnames=FALSE,
         colnames_angle=0, colnames_offset_y = .25,legend_title = "Family") +
scale_fill_viridis_d(option = "F",direction = -1)
micro <- readxl::read_xlsx(here(data),sheet = 7)
a_micro <- unique(micro$Organism) %>% str_squish(.) %>% strsplit(.," ") %>% 
  lapply(.,function(x)paste(x[1],x[2],sep="_")) %>% unlist(.) %>% unique(.)
a_macro <- sub(".*?_","",tree$tip.label)
a_macro[grep("sp$",a_macro)] <- "Pseudoeurycea_sp."
which(is.na(match(a_micro,a_macro))) -> nana
familias %>% dplyr::select(Familia) %>%  as.data.frame() -> famis_2
famis_2 %>% mutate(Familia="No Data") %>% rename("Microbioma"=Familia)-> famis_2
rownames(famis_2) <- familias$Tips
famis_2$Microbioma[match(a_micro,a_macro)] <- "With Data"
gg <- gg + new_scale_fill()
ggg <- gheatmap(gg, famis_2, offset=.01, width=.04,colnames=FALSE,
               colnames_angle=0, colnames_offset_y = .25,legend_title = "Skin microbiome") + 
  scale_fill_manual(values=c("white","black")) + theme(legend.position = "")
ggsave(ggg,file=here(output),width=10,height=10)
}

no_idea_what_this_is <- function(){ 
tree <- ape::read.nexus("data/RAxML_bipartitions.Caudata_v5.tre")
geos <- readxl::read_xlsx("data/Possible_SSM_v2.1.xlsx",sheet = 8,na="NA") %>% group_by(long,lat,Organism)
tree$tip.label %>% strsplit(.,"_") %>% lapply(.,"[",2:3) %>% 
  lapply(.,paste,collapse = "_") %>% unlist() -> tipas
geos$Organism <- sub(" ","_",geos$Organism)
geos$Organism[which(is.na(match(geos$Organism,tipas)))]
fin_tip <- tipas %>% match(geos$Organism,.) %>% .[!duplicated(.)]
tree$tip.label <- tipas
subtree <- ape::keep.tip(tree,tree$tip.label[fin_tip])
plot(subtree,cex=0.6)

data <- data.table::fread("data/GBIF.csv")
unique(data$species) %>% .[order(.)]
taxonomy <- data.table::fread("data/redlist_species_data_318f534f-7b61-4bc6-826f-3c69ad7cba8b/synonyms.csv")
taxonomy <- taxonomy %>% mutate(SynonymName = taxonomy$name %>% stringr::str_squish(.) %>% strsplit(.," ") %>% lapply(.,"[",1:2) %>%  lapply(.,paste,collapse=" ") %>% unlist())
taxonomy <- taxonomy %>% filter(infraType=="") 
data <- data %>% mutate(RESOLVED = species)
unique(data$RESOLVED) %>% .[order(.)]
ids <- taxonomy$scientificName[match(data$species,taxonomy$SynonymName)]
non <- which(is.na(ids))
data$RESOLVED[which(!is.na(ids))] <- ids[which(!is.na(ids))]
unique(data$RESOLVED) %>% .[order(.)]
#### COORDINATE CLEANER
names(data)[51] <- "BINOMIAL"
f <- CoordinateCleaner::clean_coordinates(as.data.frame(data), lon = "decimalLongitude",lat = "decimalLatitude",centroids_rad = 1000,centroids_detail = "both",tests=c("centroids","institutions","seas", "equal", "gbif","capitals", "zeros"), species=NULL,value="flagged",seas_scale = 10)
#g <- CoordinateCleaner::cc_iucn(as.data.frame(data),lon = "decimalLongitude",lat = "decimalLatitude",species="BINOMIAL",value="flagged",range=shapes)
#table(g)
table(f)
f[which(data$BINOMIAL=="Plethodon glutinosus")]
#g[which(data$BINOMIAL=="Plethodon glutinosus")]
#f_g <- as.logical(f*g)
#table(f_g)
data$Valid <- f
data %>% filter(Valid) -> data_clean
g <- raster::raster(nrows=180*12,ncols=360*12,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1,crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% as(., 'SpatialPixels')
data_clean %>% dplyr::select(.,decimalLongitude,decimalLatitude) %>% SpatialPoints(.,proj4string= CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% 
  sp::over(.,g) %>% enframe(.,name="name") %>% rename(., CellID=value) %>% bind_cols(data_clean,.) %>% as_tibble() -> data_clean
data_clean
data_clean %>% add_row(CellID=tibble(decimalLongitude= -92.34113, decimalLatitude=15.42507) %>% SpatialPoints(.,proj4string= CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% sp::over(.,g) ,BINOMIAL="Pseudoeurycea_sp",decimalLongitude = -92.34113,decimalLatitude = 15.42507,datasetKey="Robito_PersComm",class="Amphibia",order="Caudata",taxonRank="SPECIES",family="Plethodontidae",.before = T) -> data_clean
data_clean %>% add_row(CellID=tibble(decimalLongitude= -97.080, decimalLatitude=19.433) %>% SpatialPoints(.,proj4string= CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% sp::over(.,g) ,BINOMIAL="Chiropterotriton_nubilus",decimalLongitude = -97.080,decimalLatitude = 19.433,datasetKey="GarciaCastillo_etal-2018",class="Amphibia",order="Caudata",taxonRank="SPECIES",family="Plethodontidae",.before = T) -> data_clean
data_clean %>% add_row(CellID=tibble(decimalLongitude= -97.043, decimalLatitude=19.439) %>% SpatialPoints(.,proj4string= CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% sp::over(.,g) ,BINOMIAL="Chiropterotriton_nubilus",decimalLongitude = -97.043,decimalLatitude = 19.439,datasetKey="GarciaCastillo_etal-2018",class="Amphibia",order="Caudata",taxonRank="SPECIES",family="Plethodontidae",.before = T) -> data_clean
data_clean %>% add_row(CellID=tibble(decimalLongitude= -96.946, decimalLatitude=19.586) %>% SpatialPoints(.,proj4string= CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% sp::over(.,g) ,BINOMIAL="Chiropterotriton_nubilus",decimalLongitude = -96.946,decimalLatitude = 19.586,datasetKey="GarciaCastillo_etal-2018",class="Amphibia",order="Caudata",taxonRank="SPECIES",family="Plethodontidae",.before = T) -> data_clean
data_clean %>% add_row(CellID=tibble(decimalLongitude= -96.984, decimalLatitude=19.521) %>% SpatialPoints(.,proj4string= CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% sp::over(.,g) ,BINOMIAL="Chiropterotriton_nubilus",decimalLongitude = -96.984,decimalLatitude = 19.521,datasetKey="GarciaCastillo_etal-2018",class="Amphibia",order="Caudata",taxonRank="SPECIES",family="Plethodontidae",.before = T) -> data_clean
data_clean %>% add_row(CellID=tibble(decimalLongitude= -96.945, decimalLatitude=19.582) %>% SpatialPoints(.,proj4string= CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% sp::over(.,g) ,BINOMIAL="Chiropterotriton_nubilus",decimalLongitude = -96.945,decimalLatitude = 19.582,datasetKey="GarciaCastillo_etal-2018",class="Amphibia",order="Caudata",taxonRank="SPECIES",family="Plethodontidae",.before = T) -> data_clean
data_clean %>% distinct(BINOMIAL,CellID,.keep_all = T)  %>% mutate(SpeciesRESOLVED=sub(" ","_",BINOMIAL)) %>% 
  dplyr::filter(SpeciesRESOLVED %in% all_of(subtree$tip.label)) %>% 
  dplyr::select(SpeciesRESOLVED) %>% table()

lista_ras <- list.files(path="~/Documents/7.MAPOTECA/WorldClim_2/wc2.0_10m_bio",".tif$",full.names = T)
stacka <- raster::stack(lista_ras)
data_clean %>% distinct(BINOMIAL,CellID,.keep_all = T)  %>% mutate(SpeciesRESOLVED=sub(" ","_",BINOMIAL)) %>% 
  dplyr::filter(SpeciesRESOLVED %in% all_of(subtree$tip.label)) -> data_final
names(stacka) <- names(stacka) %>% sub("wc2.0_","",.) %>% sub("_10m","",.)
data_final %>% dplyr::select(decimalLongitude,decimalLatitude,BINOMIAL,CellID) -> data_omi
data_omi %>% mutate(Presence=1) %>% pivot_wider(names_from = BINOMIAL,values_from = Presence) -> data_omi
names(data_omi)[-c(1:3)] %>% sub("_"," ",.) -> Scientific
Scientific %>% strsplit(.," ") %>% do.call(rbind,.) %>% apply(.,2,substr,1,3) %>% apply(.,1,paste,collapse="") -> codes
omi_object <- list()
data_omi %>% dplyr::select(decimalLongitude,decimalLatitude) %>% raster::extract(stacka,.) %>% as.data.frame() -> omi_object$bioclims
omi_object$xy <- data_omi %>% dplyr::select(decimalLongitude,decimalLatitude) %>% as.data.frame()
omi_object
data_omi %>% dplyr::select(-c(1:3)) %>% rename_with(~all_of(codes)) %>% mutate_if(is.numeric , replace_na, replace = 0) %>% as.data.frame() -> omi_object$presence
bind_cols(Scientific=Scientific,Code=codes) %>% as.data.frame() -> omi_object$species
omi_object
omi_object$bioclims %>% apply(.,2,function(x) which(is.na(x)))

library(subniche)
dudi1 <- dudi.pca(omi_object$bioclims, scale = TRUE, scan = TRUE, nf = 3)
nic1 <- niche(dudi1, omi_object$presence, scann = TRUE)
rownames(nic1$li) <- omi_object$species$Scientific
#rownames(nic1$li)[which(rownames(nic1$li)=="Pseudoeurycea sp")] <- "Pseudoeurycea sp."
rownames(nic1$li) <- rownames(nic1$li) %>% sub(" ","_",.)
rownames(nic1$l1) <- omi_object$species$Scientific
#rownames(nic1$l1)[which(rownames(nic1$l1)=="Pseudoeurycea sp")] <- "Pseudoeurycea sp."
rownames(nic1$l1) <- rownames(nic1$l1) %>% sub(" ","_",.)
dist(nic1$l1,method = "euclidean" ) %>% as.matrix() -> omi_dist
#which(is.na(match(colnames(omi_dist),wa)))
#omi_dist[match(wa,colnames(omi_dist)),match(wa,colnames(omi_dist))] -> omi_dist
omi_dist
plot(nic1)
### phylogenetic distance vs. microbial distance
plot(as.dist(mat_dist),as.dist(micro),xlab="Patristic distance",ylab="Microbial dissimilarity")
vegan::mantel(mat_dist,micro,method = "pearson")

### ecological distance vs. phylogenetic distance
plot(y=as.dist(omi_dist),x=as.dist(mat_dist),ylab="Ecological distance",xlab="Patristic distance")
vegan::mantel(as.dist(omi_dist),mat_dist,method = "pearson")

### ecological distance vs. microbial distance
plot(as.dist(omi_dist),as.dist(micro),xlab="Ecological distance",ylab="Microbial dissimilarity")
vegan::mantel(as.dist(omi_dist),micro,method = "pearson")



micro_dist
phylo_dist
micro_dist %>% labels()
phylo_dist %>% labels()
match(rownames(nic1$li),omi_object$species$Code)

hclust(dist(nic1$l1,method = "euclidean" ),method="average") %>% as.dendrogram() -> eco_dist
eco_dist %>% labels()

dl <- dendlist(
  eco_dist %>% 
    set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3),
  phylo_dist %>% 
    set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3),
  micro_dist %>% 
    set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3)
)

names(dl) <- c("eco_dist", "phylo_dist","micro_dist")
tanglegram(dl$eco_dist,dl$phylo_dist,
           common_subtrees_color_lines = F, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd = FALSE, 
           margin_inner=7,
           lwd=2,sort=F,lab.cex=0.6,
           match_order_by_labels=T
)
}

plot_core_climate <- function(climate_core="data/Climatic_taxa_corre_v2.xlsx",abundance = "data/table_abund_4.csv",data="data/Super_table_v2.csv",output="MS/Master_Figure_S1.pdf")
{
  readxl::read_xlsx(here(climate_core),na="NA") %>% rename_with(~c("taxa","core","rho_17","p_17","rho_19","p_19","rho_15","p_15","rho_2","p_2")) -> rhos
  data.table::fread(here(abundance),sep=",",header = T) %>% as_tibble() -> abunds
  rhos %>% filter(core=="core",p_17<= 0.05 | p_15<= 0.05 |p_19<= 0.05 |p_2<= 0.05) %>% pull(taxa) ->  target
  data.table::fread(here(data),sep=",",header = T) %>% as_tibble() %>% rename(Sample="Sample_id") -> master
  master %>% select(Sample, bio15) %>% right_join(.,abunds,by="Sample") %>% relocate(bio15,.after=Bio2) %>% rename(Bio15="bio15") -> abunds
  p1 <- abunds %>% pivot_longer(cols = -c(1:9),names_to = "taxa",values_to = "abund") %>% filter(taxa %in% all_of(target)) %>% ggplot(aes(x=Bio17,y=abund,fill=taxa,color=taxa)) + geom_point(shape=21,alpha=0.6) + scale_fill_viridis_d(name="bacteria",option="F",end=0.8) + scale_color_viridis_d(name="bacteria",option="F",end=0.8)  + geom_smooth(method="lm") + theme(axis.title.x = element_text(size=14),axis.text.y = element_text(size=14),axis.text.x = element_text(size=10),axis.line = element_line(),legend.position = "bottom",panel.background = element_blank(),panel.grid = element_blank()) + scale_y_log10()
  
  p2 <- abunds %>% pivot_longer(cols = -c(1:9),names_to = "taxa",values_to = "abund") %>% filter(taxa %in% all_of(target)) %>% ggplot(aes(x=Bio15,y=abund,fill=taxa,color=taxa))  + geom_point(shape=21,alpha=0.6)+ scale_fill_viridis_d(name="bacteria",option="F",end=0.8) + scale_color_viridis_d(name="bacteria",option="F",end=0.8) + geom_smooth(method="lm") +theme(axis.title.x = element_text(size=14),axis.text.y = element_text(size=14),axis.text.x = element_text(size=10),axis.line = element_line(),legend.position = "bottom",panel.background = element_blank(),panel.grid = element_blank()) + scale_y_log10()
  
  layout <- "
AAAA
CCCC
"
  pp <- wrap_plots(list(p1,p2)) +
    plot_layout(byrow=T,nrow = 1,ncol=2,guides = "collect",tag_level = 'new') + 
    plot_annotation(title = "Abundance correlation with climate",subtitle="Bacterial Core", caption= "Rebollar et al.") & theme(legend.position = 'bottom')
  pp
  
  ggsave(filename = here(output),plot = pp)
}

plot_relative_abund <- function(data="data/level-4.csv",core_data="data/Climatic_taxa_corre_v2.xlsx",pruned_tree="data/prunned.tree"){ 
  read.table(here(data),sep=",",header=T) %>% as_tibble() -> barras
  barras %>% summarise(across(starts_with("Bacteria"), sum)) %>% pivot_longer(cols = starts_with("Bacteria"),names_to = "Bacteria",values_to = "TotAbund") %>% arrange(desc(TotAbund)) %>% mutate(CumAbund = cumsum(TotAbund)/sum(TotAbund)) %>% filter(CumAbund < 0.9) %>% pull(Bacteria) -> lista_bacs
  barras %>% group_by(Organism) %>% summarise(across(starts_with("Bacteria"), mean)) %>%  mutate_if(is.numeric,~.x/2400) %>% rowwise() %>% mutate(Suma=sum(c_across(cols = starts_with("Bacteria"))),.after=Organism) -> barras
  
  readxl::read_xlsx(here(core_data),na="NA") %>% rename_with(~c("taxa","core","rho_17","p_17","rho_19","p_19","rho_15","p_15","rho_2","p_2")) -> rhos
  rhos %>% filter(core=="core") %>% select(taxa) %>% mutate(taxa=gsub("\\|","\\.",taxa)) %>% pull(taxa) -> core
  lista_bacs[c(1,2,4,5,7)] -> core
  
  colores = met.brewer("Java",n=5,type="continuous")
  read.tree(here(pruned_tree)) %>% ladderize(.,TRUE) -> tree
  tree$tip.label %>% strsplit(.,split="_") %>% lapply(.,"[",-1) %>% lapply(.,function(x) paste(x,collapse=" ")) %>% unlist() -> tipas
  
  barras$Organism[24] <- "Pseudoeurycea sp"
  match(tipas,barras$Organism)
  barras$Organism[match(tipas,barras$Organism)] -> ordenado
  
  barras %>% ungroup() %>% mutate(Organism=factor(Organism, levels=all_of(ordenado))) -> barras
  barras %>% select(-Suma) %>% select(Organism,all_of(core)) %>%
    pivot_longer(cols = -c(1),names_to = "Bacteria",values_to = "RelAbund") %>% 
    ggplot(aes(fill=Bacteria, y=RelAbund, x = Organism)) + 
    geom_bar(position="stack", stat="identity") + 
    ylim(c(0,1)) +
    scale_fill_discrete(type = rev(colores)) +
    theme(legend.position = "",panel.background = element_blank(),panel.grid=element_blank(),
          axis.line = element_line()) + labs(x="",y="Relative Abundance") + 
    coord_flip()+ 
    NULL
}

data.table::fread("data/lefse/feature-table_4_h_01_09_22.res") %>% as_tibble() %>% rename_with(~c("Taxa","log","Habitat","LDA","p_val")) -> lefses

lefses %<>% separate(Taxa,into=c("lvl1","lvl2","lvl3","lvl4"),sep="\\.") %>% 
  filter(!is.na(lvl4)) %>% filter(lvl4!="") %>% filter(LDA>=2)

save(lefses, file=here("output_data/lefse_resultsHabitat.Rdata"))

data.table::fread("data/lefse/feature-table_family_01_09_2022.res") %>% as_tibble() %>% rename_with(~c("Taxa","log","Family","LDA","p_val")) -> lefses_fam
lefses_fam %<>% separate(Taxa,into=c("lvl1","lvl2","lvl3","lvl4"),sep="\\.") %>% 
  filter(!is.na(lvl4)) %>% filter(lvl4!="") %>% filter(LDA>=2)

save(lefses_fam, file=here("output_data/lefse_resultsFamily.Rdata"))


lefses_fam %>% count(Family) %>% summarise(sum(n))
aqui <- lefses %>% full_join(lefses_fam,.,by="lvl4") %>% #filter(!is.na(Family),!is.na(Habitat)) %>% 
  group_by(Habitat,Family) %>% summarise(n=n()) %>% 
  mutate(Habitat=ifelse(is.na(Habitat),"Neither",Habitat)) %>%  mutate(n=ifelse(Habitat=="Aquatic",n*-1,n)) %>% filter(Habitat=="Neither") %>% mutate(n=n/2)

lefses %>% full_join(lefses_fam,.,by="lvl4") %>% #filter(!is.na(Family),!is.na(Habitat)) %>% 
  group_by(Habitat,Family) %>% summarise(n=n()) %>% 
  mutate(Habitat=ifelse(is.na(Habitat),"Neither",Habitat)) %>% 
 # mutate(n=ifelse(Habitat=="Aquatic",n*-1,n)) %>% 
 #mutate(n=ifelse(Habitat=="Neither",-n/2,n)) %>% 
  #bind_rows(.,aqui) %>% 
  #mutate(n=n/sum(n)) %>% 
  #filter(Habitat!="Neither") %>% 
  ggplot(aes(y=Family,x=n,fill=factor(Habitat,levels=c("Aquatic","Terrestrial","Neither")))) +
  geom_vline(xintercept = 0) +
  geom_bar(stat = "identity",position="stack",alpha=1/1.2)+
  scale_fill_manual(values=MetBrewer::met.brewer("Hiroshige")[c(9,2,5)]) +
    theme(legend.position = "none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line=element_line(),axis.title = element_blank()) +
  NULL






data.table::fread(here("data/core/Abundancia_lvl_4.csv")) %>% as_tibble()  -> abund

lefses %>% distinct(Taxa) %>% arrange(Taxa) %>% view


abund %>% filter(grepl("Spirochaetes",V1)) %>% distinct(V1)
lefses %>% separate(Taxa,into=c("lvl1","lvl2","lvl3","lvl4"),sep="\\.") %>% 
  filter(!is.na(lvl4)) %>% filter(lvl4!="") %>% filter(LDA>=2)




abund <- abund %>% rename(Taxa=V1) %>% mutate(Taxa = gsub("__$","",Taxa)) %>% mutate(Taxa = gsub("__;","",Taxa))
lefses <- lefses %>% filter(!is.na(LDA)) %>% mutate(Taxa = gsub("\\.",";",Taxa)) %>%
  left_join(.,abund,by="Taxa") %>% filter(!is.na(`Ambystoma altamirani`))

data.table::fread(here("data/lefse/feature-table_habiatat_01_09_2022_4.txt")) %>% as_tibble(.name_repair = "unique") -> abund
abund <- abund %>% rename(Taxa=Habitat) %>%  mutate(Taxa = gsub("\\|",";",Taxa))
lefses <- lefses %>% filter(!is.na(LDA)) %>% mutate(Taxa = gsub("\\.",";",Taxa)) %>%
  left_join(.,abund,by="Taxa") %>% filter(!is.na(Aquatic...2))


lefses


lefses %>% mutate(LDA = ifelse(Class == "Aquatic",-1*LDA,LDA)) %>% 
  ggplot(aes(y=fct_reorder(Taxa,LDA),x=LDA,fill=Class)) + 
  geom_col() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line()) +
  geom_vline(xintercept = 0,size=0.2) +
  NULL

lefses %>% select(-log,-LDA,-p_val) %>% mutate(across(where(is.numeric),~.x*2500)) %>% 
  pivot_longer(-c(Taxa,Class),names_to = "Species",values_to = "abund") %>% 
  #group_by(Taxa) %>% 
  #mutate(abund = abund/sum(abund)) %>% 
  mutate(abund=ifelse(abund==0,NA,abund)) %>% 
  ggplot(aes(x=fct_reorder(Taxa,Class),y=Species,fill=abund,color=Class)) +
  geom_tile() +
  scale_fill_stepsn(colours = MetBrewer::met.brewer("Greek",direction=-1),n.breaks=10,na.value="white")  +
  #facet_grid(~V3,drop=TRUE)+
  theme(axis.text.x = element_blank())





data.table::fread(here("data/lefse/feature-table_habiatat_01_09_2022_4.txt")) %>% as_tibble(.name_repair = "unique")  -> abund

lefses %>%  mutate(V1=gsub("\\.$","",V1)) %>%  mutate(V1=gsub("\\.","|",V1)) %>% rename(Taxa=V1) %>% left_join(.,abund %>% rename(Taxa=Habitat),by="Taxa") %>% filter(!is.na(Aquatic...2)) %>% filter(V4>=2) %>%
  select(-c(V2,V4:V5)) %>% pivot_longer(-c(Taxa:V3),names_to="Sample",values_to="abund") %>% 
  mutate(Habitat = ifelse(grepl("Aquatic",Sample),"Aquatic","Terrestrial")) %>% 
  mutate(colores = paste0(V3,Habitat)) %>% 
  ggplot(aes(y=fct_reorder(Sample,Habitat),x=fct_reorder(Taxa,V3),fill=(abund/2500))) +
  geom_tile() +
  scale_fill_stepsn(colours = MetBrewer::met.brewer("Greek",direction=-1),n.breaks=10,na.value="white")  +
  theme(axis.text.x = element_blank()) +
  #coord_flip() +
  NULL
  





data.table::fread("data/lefse/feature-table_family_01_09_2022.res") %>% as_tibble() -> lefses
lefses %>% filter(V5!="-")

lefses %>% filter(!is.na(V4)) %>% #mutate(V4 = ifelse(V3 == "Aquatic",-1*V4,V4)) %>% 
  arrange(V4) %>% 
  ggplot(aes(y=fct_reorder(fct_reorder(V1,V4),V3),x=V4,fill=V3)) + geom_col() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line()) +
  geom_vline(xintercept = 0,size=0.2) +
  NULL


data = "data/Super_table_ss_10_06_21.csv"
read.table(here(data),sep=",",header=T) %>% as_tibble() %>% distinct(Organism)

infossm <- read.table(here(data),sep=",",header=T) %>% as_tibble() %>% mutate(across(15:37, ~scale(.x,scale=T,center=T),.names = "{.col}")) %>% rename("Habitat"=morphotype)  %>% mutate(faith_pd=log(faith_pd))
infossm %>% filter(Study=="PRJNA632638") %>% count(Organism)

#### RDA
dferia <- readRDS(here("data/RDA/sites.rds")) %>% as_tibble()
Habcentroids1 <- readRDS(here("data/RDA/centroids.rds")) %>% as_tibble()
continuous_arrows <-  readRDS(here("data/RDA/arrows.rds")) %>% as_tibble()
dferia %>% rename("Class"=site) %>% mutate(Class="sample") -> dferia
Habcentroids1 %>% rename("Class"=Centroids) -> Habcentroids1
continuous_arrows %>% rename("Class"=class) -> continuous_arrows
mult = 7
#continuous_arrows %>% filter(Class %in% all_of(c("sample","FamilyAmbystomatidae","FamilyCryptobranchidae","FamilyHynobiidae","FamilyPlethodontidae","FamilySalamandridae","bio1","bio2","bio16","bio18","bio19")))  -> continuous_arrows
pp <- ggplot(dferia,aes(x = CAP1, y = CAP2,fill=Family,color=Family,group=Family)) +
  geom_hline(yintercept=0,color="grey20") +
  geom_vline(xintercept=0,color="grey20") +
  geom_point(shape=19,stroke=0.3,size=2.5,alpha=0.95) +  
  stat_ellipse(size=1) +
  scale_color_manual(values= rev(MetBrewer::met.brewer("Hokusai3",5)),name="Sample") +
  scale_fill_manual(values= rev(MetBrewer::met.brewer("Hokusai3",5)),name="Sample") +
  geom_segment(data = continuous_arrows,aes(x = 0, xend = mult * CAP1,y = 0, yend = mult * CAP2,group=Class,fill=NULL), arrow = arrow(length = unit(0.2, "cm")), colour = "black") + geom_text(data = continuous_arrows,aes(x= (mult + mult/5) * CAP1, y = (mult + mult/5) * CAP2, 
label = Class,group=Class,fill=NULL), 
size = 3,parse = TRUE,color="black") +
  xlim(-5, 5) + ylim(-5, 5) +
  labs(x="CAP1 (8.6%)",y="CAP2 (6.4%)") + theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line()) + 
  NULL
pp
ggsave(filename = "MS/dbRDA.1.pdf",plot = pp)





####### FINAL ANALYSES #####




# NMDS


##### DESDE AQUI!!!!





#rownames(mat_dist)[23] <- "Pseudoeurycea_sp."
subtree2 <- subtree 
subtree2$tip.label %>% strsplit(.,"_") %>% lapply(.,"[",2:3) %>% lapply(.,paste,collapse="_") %>% unlist () -> subtree2$tip.label
#subtree2$tip.label[which(subtree2$tip.label=="Pseudoeurycea_sp")] <- "Pseudoeurycea_sp."
library(dendextend)
micro
hclust(as.dist(micro), method = "average") %>% as.dendrogram() -> micro_dist
#hclust(as.dist(mat_dist),method = "average") %>% as.dendrogram() -> d2
subtree2$edge.length <- subtree2$edge.length/max(phytools::nodeHeights(subtree2)[,2]) * 0.5
subtree2 %>% phytools::force.ultrametric(.,method = "extend") %>% as.dendrogram() -> phylo_dist



dl <- dendlist(
  phylo_dist %>% 
    set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3),
  micro_dist %>% 
    set("branches_lty", 1)
)

names(dl) <- c("phylo_dist","micro_dist")
tanglegram(dl$micro_dist,dl$phylo_dist,
           common_subtrees_color_lines = F, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd = FALSE, 
           margin_inner=7,
           lwd=2,sort=F,lab.cex=0.8,
           match_order_by_labels=T
)

a <- readxl::read_xlsx("~/Desktop/datos_campo.xlsx",sheet = 1) %>% as_tibble()
a
a %>% filter(Especie=="fulva") %>% lm(as.numeric(Circinios)~log(as.numeric(Altura)),data=.) %>% summary() -> test_lm
statistic <- paste("beta"," = ",round(test_lm$coefficients[2,1],2),", pval < 0.001"," (R^2 = ",round(test_lm$r.squared,2),")",sep="")

a %>% filter(Especie=="fulva") %>% ggplot(aes(x=log(as.numeric(Altura)),y=as.numeric(Circinios),fill=as.numeric(D_base))) +  geom_smooth(method=lm,color="black") +geom_point(size=2,shape=21,stroke=0.5,alpha=1) +theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line(),legend.position = "bottom") + labs(x="Height (log)",y="No. Croizers",subtitle="Cyathea fulva (Summer 2021)",title="Emerging fronds as a function of plant height",caption=statistic) + scale_fill_fermenter(type = "seq",palette = "Greens",direction = 1,name="Diameter")


rnaturalearth::ne_coastline(scale = 110,returnclass = "sf") -> coasr
rgdal::readOGR("~/Downloads/ne_10m_time_zones/ne_10m_time_zones.shp")%>% sf::st_as_sf() -> time_zones

ggplot() + geom_sf(data=time_zones, aes(fill=time_zone),colour="grey50") + 
  scale_fill_viridis_d(option="H") +scale_color_viridis_d(option="G") + geom_sf(data=coasr, colour="grey10",fill="NA") + theme(legend.position="")
