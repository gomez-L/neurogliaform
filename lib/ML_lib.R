

library(bmrm)
library(Biobase)


# given an object K of class dist, find n nearest neighbors
knn.index <- function(K,k=30) {
  K <- as.matrix(K)
  diag(K) <- -Inf
  t(apply(K,1,function(v) order(v)[1:k]))
}


# compute rowmeans of knn neigbors
# m: matrix to correct
# K: distance matrix between columns to consider for KNN
knn.rowmean <- function(m,K,k) {
  m <- as.matrix(m)
  nn <- knn.index(K,k)
  M <- t(apply(m,1,function(v) rowMeans(array(v[nn],dim(nn)))))
  colnames(M) <- colnames(m)
  M
}


# train multiclass SVM
feature.selected.msvm <- function(x.trn,y.trn,x.tst=x.trn,class.balancing=TRUE,n=25,LAMBDA=1,LAMBDA2=LAMBDA,...) {
  y.trn <- as.factor(y.trn)
  l <- (1 - table(seq_along(y.trn), y.trn))
  if (class.balancing) l <- balanced.loss.weights(y.trn) * l
  w <- nrbm(softMarginVectorLoss(x.trn,y.trn,l=l),LAMBDA=LAMBDA,...)
  
  W <- array(w,attr(w,"model.dim"))
  RK <- apply(W,2,function(v) rank.linear.weights(v)$rk)
  #RK <- array(rank.linear.weights(w)$rk,attr(w,"model.dim"))
  rk <- rowMin(RK)

  x.trn.reduced <- x.trn[,rk<=n]
  w.reduced <- nrbm(softMarginVectorLoss(x.trn.reduced,y.trn,l=l),LAMBDA=LAMBDA2,...)
  
  p <- predict(w.reduced,x.tst[,rk<=n])
  attr(p,"w") <- w.reduced
  attr(p,"w.global") <- w
  
  p
}


feature.selected.hinge <- function(x.trn,y.trn,x.tst=x.trn,class.balancing=TRUE,n=25,LAMBDA=1,LAMBDA2=LAMBDA,...) {
  if (is.logical(class.balancing)) {
  	if (class.balancing) {
    	lw <- balanced.loss.weights(y.trn)	
  	} else {
  	  lw <- 1
  	}
  } else {
  	lw <- class.balancing
  } 
  w <- nrbm(hingeLoss(x.trn,y.trn,loss.weights=lw),LAMBDA=LAMBDA,...)
  names(w) <- colnames(x.trn)
  rk <- rank.linear.weights(w)
  
  x.trn.reduced <- x.trn[,rk$rk<=n]
  w.reduced <- nrbm(hingeLoss(x.trn.reduced,y.trn,loss.weights=lw),LAMBDA=LAMBDA2,...)
  names(w.reduced) <- colnames(x.trn.reduced)
  
  p <- predict(w.reduced,x.tst[,rk$rk<=n])
  attr(p,"w") <- w.reduced
  attr(p,"rank") <- rk
  p
}


feature.selected.ordinal <- function(x.trn,y.trn,x.tst=x.trn,n=25,LAMBDA=1,LAMBDA2=LAMBDA,balanced=FALSE,...) {
  
  C <- costMatrix(y.trn)
  if (balanced) {
    # Adjust costs for number of sample per ordinal class
    C <- C / tabulate(y.trn)[col(C)] / tabulate(y.trn)
    C <- C / max(C)
  }

  w <- nrbm(ordinalRegressionLoss(x.trn,y.trn,C=C),LAMBDA=LAMBDA,...)
  names(w) <- colnames(x.trn)
  rk <- rank.linear.weights(w)
  
  x.trn.reduced <- x.trn[,rk$rk<=n]
  w.reduced <- nrbm(ordinalRegressionLoss(x.trn.reduced,y.trn,C=C),LAMBDA=LAMBDA2,...)
  names(w.reduced) <- colnames(x.trn.reduced)
  
  p <- predict(w.reduced,x.tst[,rk$rk<=n])
  attr(p,"w") <- w.reduced
  attr(p,"rank") <- rk
  p
}




balancedOrdinalRegressionLoss <- function(x,y) {
  C <- costMatrix(y)
  C <- C / tabulate(y.trn)[col(C)] / tabulate(y.trn)
  C <- C / max(C)
  ordinalRegressionLoss(x,y,C=C)
}


balancedOrdinalRegressionPairLoss <- function(x1,y1,x2,y2) {
  # compute the 2 ordinal losses
  L1 <- balancedOrdinalRegressionLoss(x1,y1)
  L2 <- balancedOrdinalRegressionLoss(x2,y2)
  function(w,...) {
    w1 <- L1(w,...)
    w2 <- L2(w,...)
    w <- w1
    lvalue(w) <- (lvalue(w1) + lvalue(w2))/2
    gradient(w) <- (gradient(w1) + gradient(w2))/2
    return(w)
  }
}

feature.selected.balanced.ordinal.regression.pair <- function(x1.trn,y1.trn,x2.trn,y2.trn,x.tst,n=30,keep=NULL,...) {
  if (any(colnames(x1.trn)!=colnames(x2.trn))) stop("colnames mismatch")
  
  w <- nrbm(balancedOrdinalRegressionPairLoss(x1.trn,y1.trn,x2.trn,y2.trn),...)
  names(w) <- colnames(x1.trn)
  
  rk <- rank.linear.weights(w)$rk
  rk[colnames(x.trn) %in% keep] <- 0
  
  x1.trn.reduced <- x1.trn[,rk<=n]
  x2.trn.reduced <- x2.trn[,rk<=n]
  w <- nrbm(balancedOrdinalRegressionPairLoss(x1.trn.reduced,y1.trn,x2.trn.reduced,y2.trn),...)
  names(w) <- colnames(x1.trn.reduced)
  
  f <- predict(w,x.tst[,colnames(x1.trn.reduced)])
  attr(f,"w") <- sort(w[w!=0])
  
  f
}



feature.selected.tsvm <- function(x.trn,y.trn,x.tst=x.trn,class.balancing=TRUE,n=25,LAMBDA=1,LAMBDA2=LAMBDA,...) {
  if (is.logical(class.balancing)) {
    if (class.balancing) {
      lw <- ifelse(is.na(y),1,balanced.loss.weights(y))*balanced.loss.weights(is.na(y))
    } else {
      lw <- 1
    }
  } else {
    lw <- class.balancing
  } 
  w <- nrbm(tsvmLoss(x.trn,y.trn,loss.weights=lw),LAMBDA=LAMBDA,...)
  names(w) <- colnames(x.trn)
  rk <- rank.linear.weights(w)
  
  x.trn.reduced <- x.trn[,rk$rk<=n]
  w.reduced <- nrbm(tsvmLoss(x.trn.reduced,y.trn,loss.weights=lw),LAMBDA=LAMBDA2,...)
  names(w.reduced) <- colnames(x.trn.reduced)
  
  p <- predict(w.reduced,x.tst[,rk$rk<=n])
  attr(p,"w") <- w.reduced
  attr(p,"rank") <- rk
  p
}





