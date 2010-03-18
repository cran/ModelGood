MgPlans <- function(plan,B,N,M,k){
  ReName <- NULL
  k <- as.numeric(substring(grep("^cv[0-9]+$",plan,val=TRUE,ignore.case=TRUE),3))
  if (length(k)==0) k <- NULL
  if (is.null(k)){
    re.noinf <- length(grep("noinf",plan,val=FALSE,ignore.case=TRUE))>0
    re.boot <- length(grep("boot|plain",plan,val=FALSE,ignore.case=TRUE))>0
    re.bootcv <- length(grep("bootcv|out.of.bag|out-of-bag|b0|bootcv",plan,val=FALSE,ignore.case=TRUE))>0
    re.632 <- length(grep("632",plan,val=FALSE,ignore.case=TRUE))>0
    re.plus <- length(grep("plus|\\+",plan,val=FALSE,ignore.case=TRUE))>0
    ## plan <- match.arg(plan,c("none","plain","bootcv","boot632","boot.632","boot632plus","boot.632plus","noinf"))
    if (re.noinf==TRUE){plan <- "noinf"; ReName <- "no information error"}
    else if (re.bootcv==TRUE){plan <- "bootcv"; ReName <- "bootcv error"}
    else
      if (re.boot==TRUE){
        if (re.632==TRUE){
          if (re.plus==TRUE){plan <- "boot632plus"; ReName <- ".632+"}
          else{plan <- "boot632";ReName <- ".632"}
        }
        else{stop("Plan boot not supported.");plan <- "plain"; ReName <- "bootstrap"}
      }
    if (is.null(ReName)) {plan <- "noPlan"; ReName <- "no plan"}
  }
  else{
    plan <- "crossval"
    ReName <- paste(k,"fold cross-validation",sep="-")
  }
  if (missing(M)) M <- N
  stopifnot(M>0 && M<=N) 
  subsampling <- M!=N
  ## if (missing(na.accept)) na.accept <- M/10
  if (plan=="noPlan"|| plan=="noinf") {
    B <- 0
  }
  else{
    if (missing(B)){
      if (length(k)>0) B <- 1 # repeat k-fold CrossVal ones
      else B <- 100  # either `plain' or `Bootcv'
    }
    else if (B==0) stop("No. of crossvals must be a positive integer.")
  }
  if (length(k)>0){
    CrossvalIndex <- do.call("cbind",lapply(1:B,function(b){sample(rep(1:k,length.out=N))}))
  }
  else{
    if (plan %in% c("boot632plus","bootcv","boot632","plain")){
      CrossvalIndex <- matrix(sapply(1:B,function(b){sort(sample(1:N,size=M,replace=!subsampling))}),nrow=M,ncol=B)
      colnames(CrossvalIndex) <- paste("Boot",1:B,sep=".")
    }
    else{
      CrossvalIndex <- NULL
    }
  }
  out <- list(name=ReName,internal.name=plan,index=CrossvalIndex,k=k,B=B,M=M,N=N)
  class(out) <- "MgPlans"
  out
}

