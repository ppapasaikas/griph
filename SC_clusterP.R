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



sparsify <- function (x,pct=0.1){
  L=length(x)
  if (sum(x>0)<2 ){return (x)}
  else {
    cutoff=sort(x[which(x>0)],decreasing=TRUE)[ ceiling(pct*L)  ]  
    x[which(x<cutoff)]=0  
    return(x)
  }
  
}



pruneE <- function (x,pct=0.1){
  edges <- E(GRAOp) [from(x)]  
  stopifnot (length(edges)>0)
  k=ceiling(pct * length(edges))
  e_top <- order(edges$weight, decreasing=TRUE)[1:k]
  keep_edg=edges[e_top]
}




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





SC_clusterP <- function(DM,is.cor=FALSE,impute=FALSE,filter=FALSE,diffuse.iter=NULL,p=NULL,rho=0.25,batch.penalty=0.5,ClassAssignment=rep(1,ncol(DM)),BatchAssignment=NULL,plotG=TRUE,maxG=1000,plotSP=TRUE,fsuffix=RandString() ){
  
  #######Internal parameters for testing puproses only:  
  comm.method=cluster_infomap  # Community detection algorithm. See igraph "communities" 
  ncom=NULL  #Forces the community detection algorithm to a fixed number of communities. Only possible with hierarchical methods (like fast_greedy),
  qnt=8 #Max gene expression decile for imputation (8->bottom 70% of the genes are imputed) 
  rho=rho+( (ncol(DM)/1e9)^0.15) #Scale rho for number of cells. 
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
    
  }
  
  
  
  else {
    flush.console() 
    C=list()
    C[[1]]=DM
  }
  
  
  
  
  if( is.null(p) || is.null(diffuse.iter) ){
    message("Estimating Data Separability...","\r")
    if (!isTRUE(is.cor)) { 
      HQCells=which(CellCounts >  quantile(CellCounts,probs=seq(0,1,0.02))[2] )
      SMPL=sample(HQCells,min(length(HQCells),400) )
    }
    else {
      SMPL=sample(1:nrow(C[[1]]),min(nrow(C[[1]]),400) )
    }
    Ncells=length(SMPL)
    DIST=1-WScor( C[[1]][SMPL,SMPL],WM=C[[1]][SMPL,SMPL] )  
    DIST=DIST/median(  DIST[upper.tri(DIST)]  )
    MaxD=apply(DIST,2,function(x) sort(x,decreasing=TRUE)[5]  )
    MinD=apply(DIST,2,function(x) sort(x)[5]  ) #Essentially here k=4 since minimum is the diagonal element.
    midRange=(MaxD+MinD)/2
    #midRWind=sapply(1:nrow(DIST),function(x) mad(DIST[x,-(x)] ) )/2
    midRwind=(MaxD-MinD)/13
    sep.try=mean(sapply(1:nrow(DIST),function(x) sum(DIST[x,] > 5*MinD[x] ) ) / (nrow(DIST)-1))
    separ.index=sapply(1:nrow(DIST),function(x) sum( abs(DIST[x,] - midRange[x] ) < midRwind[x] ) ) / (Ncells-1)  
    separ.index=mean(separ.index+1e-6)
    #message("separ.index1 = ",separ.index,"  Index values < 0.1 indicate datasets with sharp divisions","\r")
    #message("separ.index  = ",sep.try,"  Index values > 0.1 indicate datasets with sharp divisions","\r")
    
    if (is.null(p)) {
      p1=1.1/(0.01 + abs(log10(separ.index)) ) ; p1=max(p1,0.5)   ;  
      p2=-(log10(sep.try+1e-4)) ; p2=max(p2,0.5)   ;  
      p=(p1+p2)/2
      message("Estimated p = ",p, "  p1=",p1,"  p2=",p2, "\r")
      message("Estimated p = ",p,"\r")  
    }
    
    if (is.null(diffuse.iter)) {diffuse.iter= ifelse(separ.index < 0.2,2,1)    ;  message("Estimated diffuse.iter = ",diffuse.iter,"\r")      }
    flush.console() 
  }
  
  
  
  
  if (diffuse.iter > 1) {
    for(order in 2:diffuse.iter){
      message(paste ("Calculating Diffused Similarities...",(order-1),"/",(diffuse.iter-1)),"\r")
      flush.console() 
      C[[order]]=FlashWSpearmanCor ( C[[order-1]],WM=C[[1]],O=order-1,p=p )  
      if (order >2) {C[[order-1]]=NA}
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
  
  
  
  
  
  
  GRAO<-graph.adjacency(ADJ [rev(1:nrow(ADJ)),rev(1:nrow(ADJ))],mode=c("max"),weighted=TRUE,diag=FALSE)
  
  
  
  for (i in 1:2){
    ADJ=FlashPR(GRAO)
    ADJ[which(ADJ < (0.01/(ncol(ADJ))  ) ) ]=0
    GRAO<-graph.adjacency(ADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
  }
  
  
  pct=1
  medDEG=median(degree(GRAO))
  if (medDEG > 8 ) {
    pct=min(1,0.8/sqrt(medDEG) )
    ADJ=apply(ADJ,1,function(x) sparsify(x,pct)  )
  }
  
  
  
  
  
  
  
  
  GRAO<-graph.adjacency(ADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
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
  V(GRAO)$classcolor<-c("gold","maroon","green","blue","red","black","purple","darkorange","darkslategray","olive")[V(GRAO)$class]
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
    
    if (length(V(GRAO) ) > 1.25*maxG ) {
      message("WARNING: Graph too large (>",maxG, " vertices). A sampled subgraph of ~", maxG, " vertices will be plotted","\r")
      message("$GRAO element of the results list will contain the complete graph object","\r")
      flush.console() 
      
      DEG=degree(GRAO)
      snowball_seeds=c()
      for (i in 1:length(csize)){
        members= which(memb$membership==i)
        minDegree=quantile(DEG[members],seq(0,1,0.1)  )[5]-1
        maxDegree=quantile(DEG[members],seq(0,1,0.1)  )[6]+1
        seedn=1
        module_seeds=sample (which(memb$membership==i & DEG>minDegree & DEG<maxDegree ),seedn)
        snowball_seeds=unique(c(snowball_seeds,module_seeds))
      }
      
      snowball=c()
      seed.ego_size=0
      while (length(snowball) < maxG/2){
        seed.ego_size=seed.ego_size+1  
        snowball=unique(unlist(ego(GRAO,seed.ego_size,snowball_seeds)))
      }
      if (length(snowball) > 1.25*maxG && seed.ego_size > 1 ) {
        seed.ego_size=seed.ego_size-1
        snowball=unique(unlist(ego(GRAO,seed.ego_size,snowball_seeds))) 
      }
      
      GRAOp=induced.subgraph(GRAO,sort(snowball) )
      #edgelist=E(GRAO) [from(snowball)]
      #GRAOp=subgraph.edges(GRAO,edgelist)
      
      message("Used vertices: ", length(V(GRAOp)),"  seed_size: ",seed.ego_size , "\r")
      flush.console() 
    }
    
    
    else{GRAOp=GRAO}
    
    
    
    ncores = detectCores()
    cl<-makeCluster(  min(48,ncores,ceiling(length(V(GRAOp))/100))  )
    registerDoParallel(cl)
    #pct=min(1, 4e3/(length(V(GRAOp))^1.5) )
    pct=1
    if (median(degree(GRAOp)) > 4 ) {
      #pct=min(1,0.1/log10(median(degree(GRAOp))) )
      pct=min(1,0.8/sqrt(median(degree(GRAOp))) )
      
    }
    message("percentage of kept edges: ",pct,"\r")
    flush.console()
    KEEP=foreach(v=1:length(V(GRAOp)), .combine=c, .packages=('igraph'), .export=c('pruneE')   ) %dopar% {
      if (degree(GRAOp,v)>1){
        return(pruneE(v,pct)) 
      }
      else {return (E(GRAOp) [from(v)])  }
    }
    stopCluster(cl)
    
    GRAOp=delete_edges(GRAOp,E(GRAOp)[-KEEP])
    GRAOp=delete_vertices(GRAOp,which(ego_size(GRAOp,3) < 6))
    
    l<-layout_with_fr(GRAOp)
    add.vertex.shape("fcircle", clip=igraph.shape.noclip, plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))
    fname=paste('graph_',fsuffix,'.pdf',sep="")
    pdf(fname,12,10)
    par(mar=c(5.1, 4.1, 4.1, 14.1), xpd=TRUE)
    plot(GRAOp, layout=l,asp=0,vertex.label=NA,vertex.frame.color=V(GRAOp)$classcolor, vertex.shape="fcircle",vertex.frame.width=V(GRAOp)$frame.width )
    legend("topright",inset=c(-0.35,0),title=" Pred.      True        ",
           legend= c(sort(unique(comps)),"" ,sort(unique(ClassAssignment))   ), 
           col=c(colbar[sort(unique(comps))],"white", c("gold","maroon","green","blue","red","black","purple","darkorange","darkslategray","olive")[1:length(unique(ClassAssignment))] ), 
           pch=c(rep(20 ,length(unique(comps)) ),1, rep(21 ,length(unique(ClassAssignment))  )   ),
           bty="n", border=F, ncol=2, text.width=0.02)
    dev.off()
    GRAOp<-NULL
  }
  
  
  
  
  
  
  ADJ [rev(1:nrow(ADJ)),rev(1:nrow(ADJ))]
  
  
  
  
  
  
  ########## SPECTRAL PROJECTION ########
  Z=NULL
  if (plotSP) {
    message("Computing Spectral Projection and Rendering...","\r")
    flush.console() 
    rcomps=rev(comps)
    rclass=rev(ClassAssignment)
    symbar <- c(21,24,22,25,23,c(0:14))
    #A=FlashGKRNL(ADJ^0.2,sigma=1+ 1/(bw.nrd0 ( as.vector( (ADJ)[upper.tri(ADJ)]) ) +1e-2)  )
    A=ADJ
    A=abs(A)^p
    message("p:",p,"\r")
    A[which(ADJ<0)]=-A[which(ADJ<0)]
    A=A+1
    
    "%^%" <- function(M, power)
      with(eigen(M), vectors %*% (values^power * solve(vectors)))
    
    D <- diag(apply(A, 1, sum)) # sum rows
    U <- D-A
    #L <- diag(nrow(A)) - solve(D) %*% A #simple Laplacian
    #L <- (D %^% (-1/2)) %*% U %*% (D %^% (-1/2))  # normalized Laplacian
    L <-solve(D) %*% U #Generalized Laplacian
    
    k   <- 2
    evL <- eigen(U, symmetric=TRUE)
    Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
    Z=apply(Z,2,function(x) (x-min(x))/ (max(x)-min(x)) )
    Z=apply(Z,2,function(x) jitter(x,amount=max(IQR(x),0.25) /10 ))
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
















