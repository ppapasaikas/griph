###### Unsupervised graph-based clustering of SC data
###### DM: Count data (n genes x k cells) or directly a correlation k x k matrix. No Default
###### is.cor: Is the supplied matrix already a correlation matrix? Default: FALSE
###### impute: Perform imputation prior to clustering. Default: TRUE
###### filter: T/F Filter genes according to cv=f( mean ) fitting. Default: TRUE
###### diffuse.iter: Number of similarity difusion iterations -> [2,5]/ Default=3
###### rho: Inverse covariance matrix regularization (graph sparsification) parameter -> [0,1]. Default=0.3.  The parameter is then automatically scaled to account for number of variables 
###### or converted into a matrix and adjusted according to the batch.penalty factor to account for BatchAssignment (if given).
###### batch.penalty: [0,1] rho scaling factor for enforcing topological constraints variables according to BatchAssignment. For penalty p  -> rho_same_batch=rho^(1-p), rho_diff_batch=rho^(1+p)  
###### It is ignored if BatchAssignment==NULL 
###### ClassAssignment: If available a numeric vector of length k with numeric class labels (e.g-> c(1,2,1,1,1,2,3,3,3,1,2))
###### BatchAssignment: If available a numeric vector of length k with numeric batch labels (e.g-> c(3,1,1,2,1,3,2,2,2,3,3))
###### plotG:   if TRUE plots the resulting graph
###### plotSP: if TRUE plots the speactral projection of the graph in 2D
###### fsuffix: A suffix added to the output plots. If not given a random 5 character string.
require(igraph)       # Various methods for manipulating/analyzing graph objects
require(biotools)     # ::confusionmatrix
require(mclust)       #Optimal class label mapping and
require(doParallel)   # Parallelization
require(foreach)      # Parallelization

source ("/work/gbioinfo/papapana/FMI_groups/SingleCell_Analysis/SC_ParallelM.R")


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

RandString <- function(n=1, lenght=5){
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    lenght, replace=TRUE),collapse="")
  }
  return(randomString)
}

mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}





SC_clusterP <- function(DM,is.cor=FALSE,impute=TRUE,filter=TRUE,diffuse.iter=3,rho=0.2,batch.penalty=0.5,ClassAssignment=rep(1,ncol(DM)),BatchAssignment=NULL,plotG=TRUE,plotSP=TRUE,fsuffix=RandString() ){
  
#######Internal parameters for testing puproses only:  
comm.method=cluster_infomap  # Community detection algorithm. See igraph "communities" 
ncom=NULL  #Forces the community detection algorithm to a fixed number of communities. Only possible with hierarchical methods (like fast_greedy),
qnt=8 #Max gene expression decile for imputation (8->bottom 70% of the genes are imputed) 
rho=rho+(((ncol(DM)^1.2) /1e6)^0.25) #Scale rho for number of cells. 2000 adds ~0.3, 1e5 adds ~0.5,  
#####################################################  


if(!is.null(BatchAssignment)){
  if (length(BatchAssignment) != ncol(DM) ) {stop ("!Length of BatchAssignment vector should have length == ncol(DM)")}
}

ptm=proc.time()

  
if (!isTRUE(is.cor)) {  

if (impute==TRUE){
  message("Imputing...","\r")
  flush.console() 
  DMimp=FlashRNMF(DM)
  GF=rowSums(DM)/sum(rowSums(DM))
  QNT=quantile(GF,probs=seq(0,1,0.1))
  Low=which(GF <= QNT[qnt]  )  
  DM[Low,]=DMimp[Low,]
  DMimp=NULL
}
  
message("Preprocessing...","\r")
flush.console()    

AllZeroRows=which  ( rowSums(DM)<1e-9 )
if(length(AllZeroRows)>0){
  DM=DM[-AllZeroRows , ] 
}

corMethod="pearson"
if (  median(  apply(DM,1,function(x) sd(log2(x+2)) )) >1      ) {
  corMethod="spearman"
} 

##########  Remove invariant genes:
meanDM=mean(DM)
nSD=apply(DM,1,function(x) sd(x)/meanDM)
ConstRows=which   ( nSD < 0.25 )
if(length(ConstRows)>0){
DM=DM[-ConstRows , ]
}

#############CV=f(mean) -based filtering:
CellCounts=colSums(DM)
nDM=sweep(DM,2,CellCounts,FUN="/")
nDM=DM*10000 #Counts per 10K  
if (filter){
X1=log2(rowMeans(nDM+1/ncol(nDM)))
Y1=apply(nDM,1,function(x) log2(sd(x)/mean(x+1/length(x))+1/length(x) )  )
m=nls(Y1 ~ a*X1+b, start=list(a=-5,b=-10)  )
Yhat=predict(m)
DM=DM[which(Y1 > Yhat),]
}

NoData=which(colSums(DM)==0)
if (length(NoData>0)){
  DM=DM[,-NoData]
  ClassAssignment=ClassAssignment[,-NoData]
  BatchAssignment=BatchAssignment[,-NoData]
}
nDM<-NULL

C=list()

message("Calculating Pairwise Similarities...","\r")
flush.console() 

if (corMethod=="pearson") {C[[1]]=FlashPearsonCor(log2(DM+1))}
else {C[[1]]=FlashSpearmanCor(log2(DM+1))}

p=max(2,ceiling (0.5/ IQR((C[[1]])[upper.tri(C[[1]])])/mean(C[[1]][upper.tri(C[[1]])])) )
p=min(8,p)
C[[1]]=FlashGKRNL( (C[[1]]^p),sigma=1)  #Kernelized distance metric 

for(order in 2:diffuse.iter){
message(paste ("Calculating Diffused Similarities...",(order-1),"/",(diffuse.iter-1)),"\r")
flush.console() 
C[[order]]=FlashSpearmanCor( C[[order-1]])  
C[[order-1]]=NA
}

}

  
else {
flush.console() 
C=list()
C[[1]]=DM

p=max(2,ceiling (0.5/ IQR((C[[1]])[upper.tri(C[[1]])])/mean(C[[1]][upper.tri(C[[1]])])) )
p=min(8,p)
C[[1]]=FlashGKRNL( (C[[1]]^p),sigma=1)  #Kernelized distance metric 

for(order in 2:diffuse.iter){
  message(paste ("Calculating Diffused Similarities...",(order-1),"/",(diffuse.iter-1)),"\r")
  flush.console() 
C[[order]]=FlashSpearmanCor( C[[order-1]])  
C[[order-1]]=NA
} 
}

ClassAssignment=rev(as.numeric(factor(ClassAssignment)))
  
############### glasso-based graph structure estimation: #####################
message("Estimating Graph Structure...","\r")
flush.console() 

Cuse=C[[diffuse.iter]]
if(!is.null(BatchAssignment)){
  rL=min(1,rho^(1-batch.penalty))
  rS=rho^(1+batch.penalty)
  rho=mapply (  function(r,c) {if (BatchAssignment[r]==BatchAssignment[c]) {rL} else { rS } },row(Cuse),col(Cuse) )
  rho=matrix(rho,nrow=nrow(Cuse),ncol=ncol(Cuse) )
}


C=NULL
#Cuse=as.matrix(nearPD(Cuse,corr = TRUE)$mat) ##### Force Positive Definite Matrix

ADJ=-FlashGlasso(Cuse,rho=rho)
ADJ[which(ADJ>0)]=(Cuse[which(ADJ>0)])^2
ADJ[which(ADJ < 0)]=0
diag(ADJ)=1

GRAO<-graph.adjacency(ADJ[rev(1:nrow(ADJ)),rev(1:nrow(ADJ))],mode=c("max"),weighted=TRUE,diag=FALSE)
GRAO<-set.vertex.attribute(GRAO,"labels",value=rev(rownames(Cuse)) )
GRAO<-set.vertex.attribute(GRAO,"ids",value=V(GRAO))
GRAO<-set.vertex.attribute(GRAO,"class",value=ClassAssignment  )
Cuse<-NULL

######## COMMUNITY DETECTION #########
message("Detecting Graph Communities...","\r")
flush.console() 
memb=comm.method(GRAO)
if (!is.null(ncom)) {
  memb$membership=cut_at(memb,no=ncom)  
}
csize=table(memb$membership)

######### Optimal Mapping between true and estimated class assignments: ##########
misclErr=classError(ClassAssignment,memb$membership) #should be the same as classError(ClassAssignment,memb$membership) before mapping
mapping=as.vector(unlist(mapClass(memb$membership,ClassAssignment)[[2]]))
newmemb=memb$membership

################# Graph Layout and Plotting ####################==
comps <-memb$membership
colbar  <- gg_color_hue(length(table(comps)))
V(GRAO)$color <- colbar[comps]
V(GRAO)$shape<-"circle"
V(GRAO)$size<-10 /(length(memb$membership)/60 )^0.3
V(GRAO)$cex<-V(GRAO)$size / 3
V(GRAO)$frame.width=2 /(length(memb$membership)/60 )^0.3
V(GRAO)$classcolor<-c("gold","maroon","green","blue","red","black","indigo","darkorange","darkslategray","olive")[V(GRAO)$class]
E(GRAO)$width<-E(GRAO)$weight*2.0 / sqrt((length(memb$membership)/60 ))

pruneE <- function (x,pct=0.1){
  edges <- E(GRAO) [from(x)]  
  stopifnot (length(edges)>0)
  k=ceiling(pct * length(edges))
  e_top <- order(edges$weight, decreasing=TRUE)[1:k]
  keep_edg=edges[e_top]
}

if (plotG) {
  message("Computing Graph Layout and Rendering...","\r")
  flush.console() 
  ncores = detectCores()
  cl<-makeCluster(  min(48,ncores,ceiling(length(V(GRAO))/100))  )
  registerDoParallel(cl)
  pct=min(1, 4e3/(length(V(GRAO))^1.5) )
  KEEP=foreach(v=1:length(V(GRAO)), .combine=c, .packages=('igraph')  ) %dopar% {
  return(pruneE(v,pct)) 
  }
  stopCluster(cl)
  
  g=delete_edges(GRAO,E(GRAO)[-KEEP])
  g=delete_vertices(g,which(ego_size(g,4) < 6))

  l<-layout_with_fr(g)

  add.vertex.shape("fcircle", clip=igraph.shape.noclip, plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))
  fname=paste('graph_',fsuffix,'.pdf',sep="")
  pdf(fname,12,10)
  par(mar=c(5.1, 4.1, 4.1, 14.1), xpd=TRUE)
  plot(g, layout=l,asp=0,vertex.label=NA,vertex.frame.color=V(g)$classcolor, vertex.shape="fcircle",vertex.frame.width=V(g)$frame.width )
  legend("topright",inset=c(-0.35,0),title=" Pred.      True        ",
         legend= c(sort(unique(comps)),"" ,sort(unique(ClassAssignment))   ), 
         col=c(colbar[sort(unique(comps))],"white", c("gold","maroon","green","blue","red","black","purple","darkorange","darkslategray","olive")[1:length(unique(ClassAssignment))] ), 
         pch=c(rep(20 ,length(unique(comps)) ),1, rep(21 ,length(unique(ClassAssignment))  )   ),
         bty="n", border=F, ncol=2, text.width=0.02)
  dev.off()
  g<-NULL
}

########## SPECTRAL PROJECTION ########
if (plotSP) {
message("Computing Spectral Projection and Rendering...","\r")
flush.console() 
rcomps=rev(comps)
rclass=rev(ClassAssignment)
symbar <- c(21,24,22,25,23,c(0:14))
A=FlashGKRNL(ADJ^0.2,sigma=1+ 1/(bw.nrd0 ( as.vector( (ADJ)[upper.tri(ADJ)]) ) +1e-2)  )
D <- diag(apply(A, 1, sum)) # sum rows
U <- D-A
L <- diag(nrow(A)) - solve(D) %*% A
k   <- 2
evL <- eigen(U, symmetric=TRUE)
Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
Z=apply(Z,2,function(x) (x-min(x))/ (max(x)-min(x)) )
Z=apply(Z,2,function(x) jitter(x,amount=IQR(x)/10 ))
A<-NULL; D<-NULL; U<-NULL; L<-NULL; evL<-NULL


fname=paste('specProj_',fsuffix,'.pdf',sep="")
pdf(fname,12,10)
par(mar=c(5.1, 4.1, 4.1, 14.1), xpd=TRUE)

plot(Z, col=colbar[rcomps], pch=symbar[rclass],cex=10/(nrow(Z)^0.35 ),bty='L',xlab="SD1",ylab="SD2",
     main="Spectral Projection")

legend("topright",inset=c(-0.35,0),title=" Pred.      True        ",
       legend= c(sort(unique(rcomps)),"" ,sort(unique(rclass))   ), 
       col=c(colbar[sort(unique(rcomps))],"white", rep("black",length(unique(rclass)) ) ), 
       pch=c(rep(20 ,length(unique(rcomps)) ),1, symbar[sort(unique(rclass))]),
       bty="n", border=F, ncol=2, text.width=0.02)
dev.off()
}
#########################################
Te=(proc.time() - ptm)[3]
Te=signif(Te,digits=6)
message("Done...","\r",appendLF=FALSE)
message(paste("Elapsed Time: ", Te),"\r")
flush.console() 

ConfMatrix=confusionmatrix(newmemb,ClassAssignment)
colnames(ConfMatrix)=paste(rownames(ConfMatrix),".true",sep="")
rownames(ConfMatrix)=paste(rownames(ConfMatrix),".pred",sep="")
return(list(MEMB=rev(newmemb),DISTM=ADJ, specp=Z, ConfMatrix=ConfMatrix,miscl=misclErr,GRAO=GRAO,rho=rho  ))
}
















