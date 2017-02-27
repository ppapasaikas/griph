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
###### BatchAssignment: If available a numeric vector of length k with numeric batch labels (e.g-> c(1,2,1,1,1,2,3,3,3,1,2))
###### plotG:   if TRUE plots the resulting graph
###### plotSP: if TRUE plots the speactral projection of the graph in 2D
###### fsuffix: A suffix added to the output plots. If not given a random 5 character string.
require(igraph)       # Various methods for manipulating/analyzing graph objects
require(QUIC)         # Graphical Lasso using Newton's method for coordinate descent
require(biotools)     # ::confusionmatrix
require(coop)         # (pcor) Fast correlation calculations using the BLAS libraries
require(corpcor)      #(cor.shrink)  Incorporates weights and Shrinkage
require(mclust)       # Optimal class label mapping and
require(rNMF)         # rnMF based imputation
require(KRLS)         # gaussian kernel distance
require(doParallel)   # Parallelization
require(foreach)      # Parallelization



Spcor <- function (M) { #Spearman's correlation using coop
  R=vapply(c(1:ncol(M)),function (x) rank(M[,x]),FUN.VALUE=double(length=nrow(M) ) )
  R=pcor(R)
  dimnames(R)=list(colnames(M),colnames(M))
  return(R)
}



WScor <- function (M, WM=matrix(1,ncol(M),ncol(M) ),O=1,l=0,p=1  ) { #Spearman's weighted correlation using cor.shrink
  R=vapply(c(1:ncol(M)),function (x) rank(M[,x]),FUN.VALUE=double(length=nrow(M) ) )
  s=10^(2*O-2)/2
  s=max(3,s)
  s=min(100,s)
  GK=gausskernel((WM),sigma=s ) #0.5-1
  W=apply(GK ,1,sum )^(3) #
  minW=quantile(W,probs=seq(0,1,0.1))[[2]]
  W=((minW)/(W+1*minW))  #
  W=W/max(W)
  R1=cor.shrink( R^(2.5*p),verbose=FALSE,lambda=0,w=W ) 
  R2=cor.shrink( R^(12*p),verbose=FALSE,lambda=0,w=W ) 
  R=(R1+R2)/2
  dimnames(R)=list(colnames(M),colnames(M))
  return(R)
}



PPR <- function (G,df=0.75){
if (!isTRUE(is.igraph(G))) {  
if (!isSymmetric(G) )  {stop ("!! G should be either a graph object or a symmetric matrix !!")}   
G=graph.adjacency(G [rev(1:nrow(G)),rev(1:nrow(G))],mode=c("max"),weighted=TRUE,diag=FALSE)
}
  
L=length(V(G))
PR=diag(nrow=L)
vals=sapply(1:L, function(x) page_rank(G, vids=c(1:L),personalized=PR[,x],damping=df )$vector   )
PR[1:length(PR)]=log2(unlist(vals)+(0.01/L))
PR[lower.tri(PR)]=(PR[lower.tri(PR)]+t(PR)[lower.tri(t(PR))])/2
PR[upper.tri(PR)]=t(PR)[upper.tri(t(PR))]
PR=2^PR
return(PR)
}



sparsify <- function (x,pct=0.1){
  N=which(x>0)
  L=length(N)
  if (L<2 ){return (x)}
  else {
    cutoff=sort(x[N],decreasing=TRUE)[ 1+ceiling(pct*L)  ]  
    x[which(x<cutoff)]=0  
    return(x)
  }
}



pruneE <- function (x,pct=0.1){
  edges <- E(GRAOp) [from(x)]  
  stopifnot (length(edges)>2)
  k=max(3,floor(pct * length(edges)))
  k=min(100,k)
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





SC_cluster <- function(DM,is.cor=FALSE,impute=FALSE,filter=FALSE,rho=0.3,diffuse.iter=2,p=1,batch.penalty=0.5,ClassAssignment=rep(1,ncol(DM)),BatchAssignment=NULL,plotG=TRUE,maxG=1000,plotSP=TRUE,fsuffix=RandString() ){
  
  #######Internal parameters for testing puproses only:  
  comm.method=cluster_infomap  # Community detection algorithm. See igraph "communities" 
  ncom=NULL  #Forces the community detection algorithm to a fixed number of communities. Only possible with hierarchical methods (like fast_greedy),
  qnt=8 #Max gened expression decile for imputation (8->bottom 70% of the genes are imputed) 
  rho=rho+( (ncol(DM)/1e9)^0.2) #Scale rho for number of cells. 
  
  
  if(!is.null(BatchAssignment)){
    if (length(BatchAssignment) != ncol(DM) ) {stop ("!Length of BatchAssignment vector should have length == ncol(DM)")}
  }
  
  ptm=proc.time()
  
  
  if (!isTRUE(is.cor)) {  
    
    if (impute==TRUE){
      message("Imputing...","\r")
      flush.console() 
      DMimp=rnmf(DM,k =6,alpha = 0.15,tol=1e-2,maxit=10,showprogress=FALSE,quiet=TRUE)$fit
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
    if (  median(  apply(DM,1,function(x) sd(log2(x+2)) )) >1      ) { #Switch to spearman correlation in presence of severe outliers
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
      CellCounts=CellCounts[-NoData]
    }
    nDM<-NULL
    
    C=list()
    
    message("Calculating Pairwise Similarities...","\r")
    flush.console() 
    
    if (corMethod=="pearson") {C[[1]]=pcor(log2(DM+1))} #Switch to spearman correlation in presence of severe outliers
    else {C[[1]]=Spcor(log2(DM+1))}
    
  }
  
  
  else {
    flush.console() 
    C=list()
    C[[1]]=DM
  }
  
  
  
  if (diffuse.iter > 1) {
  for(order in 2:diffuse.iter){
    message(paste ("Calculating Diffused Similarities...",(order-1),"/",(diffuse.iter-1)),"\r")
    flush.console() 
    C[[order]]=WScor( C[[order-1]],WM=C[[1]],O=order-1,p=p )
    #C[[order]]=Spcor( C[[order-1]] )  
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
  tol=2e-02
  maxIter=50
  if (ncol(Cuse)<200) {tol=5e-03;maxIter=200}
  if (ncol(Cuse)>800) {tol=8e-02;maxIter=25}
  #message("Rho: ",summary(rho),"\r")
  G<-QUIC(Cuse,rho=rho,tol=tol,maxIter=maxIter,msg=0)  #0.1 for C4/ 0.3 for C3 / 0.5 for C2
  
  ADJ= -G$X
  ADJ[which(ADJ>0)]=(Cuse[which(ADJ>0)])^2
  ADJ[which(ADJ>0)]=exp(-  (1-Cuse[which(ADJ>0)])^2 / (1-mean(Cuse[upper.tri((Cuse))]))^2)   #From Haren and Koren 2001 section 3
  ADJ[which(ADJ < 0)]=0
  G<-NULL
  diag(ADJ)=1
  

  GRAO<-graph.adjacency(ADJ [rev(1:nrow(ADJ)),rev(1:nrow(ADJ))],mode=c("max"),weighted=TRUE,diag=FALSE)
  
  
  
  niter=2
  for (i in 1:niter){
  df=0.75
  ADJ=PPR(GRAO,df=df)
  ADJ[which(ADJ < (0.01/(ncol(ADJ))  ) ) ]=0
  medDEG=median(apply(ADJ,1,function(x) (sum(x>0)-1)  ))
  thr=1000
  if (i==niter) {thr=5}
    if (medDEG > thr ) {
    pct=min(1,1/sqrt(medDEG/thr) )
    ADJ=apply(ADJ,1,function(x) sparsify(x,pct)  )
  }
  GRAO<-graph.adjacency(ADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
  }

  
  
  
  
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
        minDegree=quantile(DEG[members])[3]-1
        maxDegree=quantile(DEG[members])[4]+1
        seedn=1
        module_seeds=sample (which(memb$membership==i & DEG >= minDegree & DEG <= maxDegree ),seedn)
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
      
      message("Used vertices: ", length(V(GRAOp)),"  seed_size: ",seed.ego_size , "\r")
      flush.console() 
    }
    
    else{GRAOp=GRAO}
    
    
    
    ncores = detectCores()
    cl<-makeCluster(  min(16,ncores,ceiling(length(V(GRAOp))/100))  )
    registerDoParallel(cl)
    pct=1
    if (median(degree(GRAOp)) > 4 ) {
      pct=min(1,1/sqrt(1*median(degree(GRAOp)) ) )
      
    }
    message("percentage of displayed edges: ",pct,"\r")
    flush.console()
    KEEP=foreach(v=1:length(V(GRAOp)), .combine=c, .packages=('igraph'), .export=c('pruneE')  ) %dopar% {
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
  
  
  
  
  
  ADJ=ADJ [rev(1:nrow(ADJ)),rev(1:nrow(ADJ))]
  
  
  
  
  ########## SPECTRAL PROJECTION ########
  Z=NULL
  if (plotSP) {
    message("Computing Spectral Projection and Rendering...","\r")
    flush.console() 
    rcomps=rev(comps)
    rclass=rev(ClassAssignment)
    symbar <- c(21,24,22,25,23,c(0:14))
    A=ADJ
    A=log2(A + 0.1/ncol(A))
    A=A-min(A)
    A=A/max(A)
    #n=rnorm(length(A),0,0.1)
    #A=A+n
    #A=(A+1)
    #diag(A)=100
    

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
    #Z=apply(Z,2,function(x) jitter(x,amount=max(IQR(x),0.25) /10 ))
    D<-NULL; U<-NULL; L<-NULL; evL<-NULL

    mZ1=mean(Z[,1])
    mZ2=mean(Z[,2])
    outl=which(  abs(Z[,1]-mZ1) > 0.9 | abs(Z[,2]-mZ2)>0.9    )
    message(length(outl),"\r")
    flush.console() 
    
    
    if (length(outl)>0 && length(outl) < 8 ) {
    A=A[-outl,-outl]
    rcomps=rcomps[-outl]
    rclass=rclass[-outl]
    D <- diag(apply(A, 1, sum)) # sum rows
    U <- D-A
    L <-solve(D) %*% U #Generalized Laplacian
    evL <- eigen(U, symmetric=TRUE)
    Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
    Z=apply(Z,2,function(x) (x-min(x))/ (max(x)-min(x)) )
    Z=apply(Z,2,function(x) jitter(x,amount=max(IQR(x),0.25) /10 ))
    D<-NULL; U<-NULL; L<-NULL; evL<-NULL
    }
    A<-NULL;
    
    #PCA=prcomp(x=Z,scale=T) 
    
    
    fname=paste('specProj_',fsuffix,'.pdf',sep="")
    pdf(fname,12,10)
    par(mar=c(5.1, 4.1, 4.1, 14.1), xpd=TRUE)
    
    #plot(PCA$x[,1:2], col=colbar[rcomps], pch=symbar[rclass],cex=10/(nrow(Z)^0.35 ),bty='L',xlab="SD1",ylab="SD2",
    #     main="Spectral Projection")
    
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
  
  return(list(MEMB=rev(newmemb), DISTM=ADJ, specp=Z, ConfMatrix=ConfMatrix,miscl=misclErr,GRAO=GRAO  ))
}
















