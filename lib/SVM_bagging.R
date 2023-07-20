
library(SummarizedExperiment)
library(igraph)
library(bmrm)


T18.ontology <- function(T18,level=2) {
  dag <- graph_from_edgelist(unique(cbind(as.character(T18$class),"root")))
  if (level>=2) {
    dag <- dag + graph_from_edgelist(unique(cbind(as.character(T18$subclass),as.character(T18$class))))
    label <- T18$subclass[drop=TRUE]
  }
  if (level>=3) {
    dag <- dag + graph_from_edgelist(unique(cbind(as.character(T18$cluster),as.character(T18$subclass))))
    label <- T18$cluster[drop=TRUE]
  }
  list(dag=dag,label=label)
}


T18.train.bag <- function(T18,nGene=1000,nBag=16,nCell=100,mc.cores=6,level=3,LAMBDA=1e6,...) {
  dag <- T18.ontology(T18,level=level)
  
  # make the training set with log2(RPM+1)
  X <- t(assay(T18)) / T18$nUMI*1e6
  X@x <- log2(X@x+1)
  y <- dag$label
  
  
  DAG <- is.finite(distances(dag$dag,mode = "out"))+0
  DAG <- DAG[levels(y),setdiff(colnames(DAG),"root")]
  l <- distances(dag$dag)[levels(y),levels(y)]/2

  # train a bag of ontology-based multiclass SVM on (X,y,DAG,l)
  bag <- mclapply(seq(nBag),mc.cores=mc.cores,function(seed) {
    set.seed(seed)
    trn.cell <- unlist(lapply(levels(y),function(v) {sample(which(y %in% v),nCell,replace = TRUE)}))
    trn.gene <- sample(seq(ncol(X)),nGene)
    
    trn.X <- cbind(intercept=1000,scale(as.matrix(X[trn.cell,trn.gene]),center=TRUE,scale=FALSE))
    trn.y <- y[trn.cell]
    w <- nrbm(ontologyLoss(trn.X,trn.y,dag=DAG,l=l[trn.y,]),LAMBDA=LAMBDA,...)
    
    dim(w) <- attr(w, "model.dim")
    dimnames(w) <- attr(w, "model.dimnames")
    gradient(w) <- attr(w, "model.dim") <- attr(w, "model.dimnames") <- NULL
    return(w)
  })
  
  list(bag = bag,DAG = DAG,dag = dag)
}


T18.bag.predict <- function(o,X,mc.cores=6,DAG=o$DAG) {
  X <- cbind(intercept=1000,X)
  P <- mclapply(o$bag,mc.cores=mc.cores,function(W) {
    fn <- intersect(rownames(W),colnames(X))
    f <- as.matrix(X[,fn] %*% W[fn,])
    if (!is.null(DAG)) f <- tcrossprod(f,DAG)
    f
  })
  simplify2array(P)
}


