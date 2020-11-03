
rr <- function(y, Z, K=NULL, X=NULL, method="REML"){
  stopifnot(is.matrix(Z))
  out <- rrBLUP::mixed.solve(y=y, Z=Z, K=K, X=X, method=method)
  return(structure(out, class="rr"))
}
predict.rr <- function(object, newZ){
  stopifnot(is.matrix(newZ))
  out <- as.vector(newZ %*% object$u)
  if(! is.null(rownames(newZ)))
    names(out) <- rownames(newZ)
  return(out)
}

### Partitions


folds <- cvFolds(n=nrow(X), K=5, R=10)


### Validation

callRR <- call("rr", y=y, Z=X)
system.time(
    out.cv <- cvTool(call=callRR, x=X, y=y, names=c("Z", "y"),
                     cost=cor, folds=folds))
out.cv # one row per replicate
mean(out.cv[,"CV"])
sd(out.cv[,"CV"])
