# {{{ UseMethod

Roc <- function(object,...){
  UseMethod("Roc",object=object)
}

# }}}
# {{{ Roc.default, Roc.glm,etc.

Roc.default <- function(object,
                        y,
                        breaks,
                        crRatio=1,
                        ## pv=FALSE,
                        ## confint=FALSE,
                        ## confint.method="exact",
                        keep.tables=FALSE,
                        keep.breaks=FALSE){
  N <- length(object)
  if(length(y)!=N) stop("Arguments must have the same length")
  if(length(unique(y))!=2) stop("y must be binary")
  Disease <- as.integer(as.character(factor(y,labels=c("0","1"))))
  count.DiseasePos <- sum(Disease==1)
  count.DiseaseNeg <- sum(Disease==0)
  if (missing(breaks))
    breaks <- sort(unique(object))
  else
    breaks <- sort(unique(breaks))
  if (length(breaks)>1 & !is.factor(object)){
    eval.times <- breaks- min(diff(breaks))/2
    count.TestPos <- N-MgSindex2(jump.times=object,eval.times=eval.times)
    count.TestNeg <- N-count.TestPos
    a <- count.DiseasePos-MgSindex2(jump.times=object[Disease==1],eval.times=eval.times)
    b <- count.DiseaseNeg-MgSindex2(jump.times=object[Disease==0],eval.times=eval.times)
  }
  else{
    tabx <- table(object)
    count.TestPos <- tabx
    count.TestNeg <- tabx
    tabxpos <- table(object[Disease==1])
    tabxneg <- table(object[Disease==0])
    a <- count.DiseasePos-tabxpos
    b <- count.DiseaseNeg-tabxneg
  }
  c <- count.DiseasePos-a
  d <- count.DiseaseNeg-b
  sens <- crRatio*a/count.DiseasePos
  spec <- d/count.DiseaseNeg
  out <- list("Sensitivity"=c(sens,0),"Specificity"=c(spec,1))
  ##   if (confint==TRUE){
  ##     tmp.sens <- binconf(x=a,n=count.DiseasePos,method=confint.method)
  ##     sens <- tmp.sens[,"PointEst"]
  ##     ci.sens <- tmp.sens[,c("Lower","Upper")]
  ##     tmp.spec <- binconf(x=d,n=count.DiseaseNeg,method=confint.method)
  ##     spec <- tmp.spec[,"PointEst"]
  ##     ci.spec <- tmp.spec[,c("Lower","Upper")]
  ##     out <- list("Sensitivity"=c(sens,0),"Specificity"=c(spec,1),
  ##                 "CI.Sens"=rbind(ci.sens,c(0,0)),"CI.Spec"=rbind(ci.spec,c(1,1)))
  ##   }
  ##   if (pv==TRUE){
  ##     if (confint==TRUE){
  ##       tmp.ppv <- binconf(x=a,n=count.TestPos,method=confint.method)
  ##       ppv <- tmp.ppv[,"PointEst"]
  ##       ci.ppv <- tmp.ppv[,c("Lower","Upper")]
  ##       tmp.npv <- binconf(x=d,n=count.TestNeg,method=confint.method)
  ##       npv <- tmp.npv[,"PointEst"]
  ##       ci.npv <- tmp.npv[,c("Lower","Upper")]
  ##       out <- c(out,list("PPV"=c(ppv,1),"NPV"=c(npv,npv[length(npv)]),
  ##                         list("CI.PPV"=rbind(ci.ppv,c(1,1)),
  ##                              "CI.NPV"=rbind(ci.npv,ci.npv[NCOL(ci.npv),]))))
  ##     }
  ##     else{
  ##   ppv <- a/count.TestPos
  ##   npv <- d/count.TestNeg
  ## }
  ##   out <- c(out,list("PPV"=c(ppv,1),"NPV"=c(npv,npv[length(npv)])))
  ## }
  if (keep.breaks==TRUE)
    out <- c(out,list(breaks=breaks))
  if (keep.tables==TRUE){
    out <- c(out,list(tables=data.frame(a,b,c,d)))
  }
  ##   if (confint==TRUE)
  ##     out$confint.method <- confint.method
  ## class(out) <- "Roc"
  out
}


Roc.formula <- function(object,formula,data,crRatio=1,...){
  call <- match.call()
  m <- match.call(expand = FALSE)
  m <- m[match(c("","formula","data","subset","na.action"), names(m), nomatch = 0)]
  if (missing(data)) Terms <- terms(formula)
  else Terms <- terms(formula,data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m1 <- eval(m, parent.frame())
  if (NROW(m1) == 0) stop("No (non-missing) status values")
  m$formula <- delete.response(Terms)
  m2 <- eval(m, parent.frame())
  if (NROW(m2) == 0) stop("No (non-missing) test values")
  y <- model.extract(m1,"response")
  x <- m2[,1,drop=TRUE]
  if(any(match(c("Surv","Hist"),all.names(formula),nomatch=0)>0)){
    stop("Survival models not supported. Refer to the R package `pec'.")
  }
  else{
    out <- list("Roc"=list(threshModel=Roc.default(object,y, crRatio=crRatio,...)),
                call=match.call(),
                models=c("threshModel"=formula))
    class(out) <- "Roc"
    out
  }
}

Roc.glm <- function(object,formula,data,...){
  stopifnot(object$family$family=="binomial")
  Roc.list(object=list(object),formula,data,...)
}

Roc.lrm <- function(object,formula,data,...){
  Roc.list(object=list(object),formula,data,...)
}


Roc.rpart <- function(object,formula,data,...){
  Roc.list(object=list(object),formula,data,...)
}

Roc.randomForest <- function(object,formula,data,...){
  Roc.list(object=list(object),formula,data,...)
}

# }}}
# {{{ average Roc curves
avRoc <- function(list,grid,method="vertical"){
  if (missing(grid))
    grid <- switch(method,"vertical"=seq(0,1,.01),"horizontal"=seq(1,0,-.01))
  if (method=="vertical"){
    meanSens <- rowMeans(do.call("cbind",lapply(list,function(Roc){
      approx(x=Roc$Specificity,y=Roc$Sensitivity,xout=grid,ties=median,yleft=0,yright=1)$y
    })))
    meanRoc <- list(Sensitivity=c(meanSens,0),Specificity=c(grid,1))
  }
  else
    if (method=="horizontal"){
      meanSpec <- rowMeans(do.call("cbind",lapply(list,function(Roc){
        approx(x=Roc$Sensitivity,y=Roc$Specificity,xout=grid,ties=median,yleft=0,yright=1)$y
      })))
      meanRoc <- list(Sensitivity=c(grid,0),Specificity=c(meanSpec,1))
    }
  meanRoc
}
# }}}
# {{{ Roc.list

Roc.list <- function(object,
                     formula,
                     data,
                     plan="noPlan",
                     noinf.method=c("simulate"),
                     simulate="reeval",
                     B,
                     M,
                     breaks,
                     crRatio=1,
                     RocAverageMethod="vertical",
                     RocAverageGrid=switch(RocAverageMethod,,"vertical"=seq(0,1,.01),"horizontal"=seq(1,0,-.01)),
                     model.args=NULL,
                     model.parms=NULL,
                     keepModels=FALSE,
                     keepSampleIndex=FALSE,
                     keepCrossValRes=FALSE,
                     keepNoInfSimu,
                     na.accept=0,
                     verbose=TRUE,
                     ...){

# }}}
# {{{ models
  NF <- length(object) 
  if (is.null(names(object)))names(object) <- sapply(object,function(o)class(o)[1])
  else{names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])}
  names(object) <- make.names(names(object),unique=TRUE)

  # }}}
# {{{ formula
  if (missing(formula)){
    formula <- eval(object[[1]]$call$formula)
    if (class(formula)!="formula")
      stop("Argument formula is missing.")
    else
      if (verbose)
        warning("Argument formula is missing. I use the formula from the call to the first model instead.")
  }

  # }}}
# {{{ data
  if (missing(data)){
    data <- eval(object[[1]]$call$data)
    if (class(data)!="data.frame")
      stop("Argument data is missing.")
    else
      if (verbose)
        warning("Argument data is missing. I have (ab)used the data from the call\n of the first model instead.")
  
  }
  # }}}
# {{{ response

  m <- model.frame(formula,data,na.action=na.fail)
  Y <- model.response(m)
  
  if (is.factor(Y) && (length(levels(Y))==2) || length(unique(Y))==2) {
    Y <- factor(Y)
    Y <- as.numeric(Y==levels(Y)[2])
  }
  N <- length(Y)
  
  # }}}
# {{{ break points for the ROC
  if (missing(breaks))
  breaks <- seq(0,1,.01)
  
# }}}
# {{{ Plan

  Plan <- MgPlans(plan=plan,B=B,N=N,M=M,k=k)
  B <- Plan$B
  CrossvalIndex <- Plan$index
  if (!keepSampleIndex) Plan$index <- NULL
  k <- Plan$k
  do.crossval <- !(is.null(CrossvalIndex))
  if (missing(keepCrossValRes)) keepCrossValRes <- do.crossval
  if (missing(keepNoInfSimu)) keepNoInfSimu <- FALSE

  # }}}
# {{{ checking the models for compatibility with cross-validation
  if (do.crossval){
    cm <- MgCheck(object=object,model.args=model.args,model.parms=model.parms,Plan=Plan,verbose=verbose)
    model.args <- cm$model.args
    model.parms <- cm$model.parms
  }

  # }}}
# {{{ computation of ROC curves in a loop over the models 
  
  list.out <- lapply(1:NF,function(f){
   
    if (verbose && NF>1) message("\n",names(object)[f],"\n")
    fit <- object[[f]]
    extract <- model.parms[[f]]
    
  # }}}
# {{{ apparent ROC (use the same data for fitting and validation)

    pred <- do.call("predictStatusProb",c(list(object=fit,newdata=data),model.args[[f]]))
    AppRoc <- Roc.default(object=pred,y=Y,breaks=breaks,crRatio=crRatio)
    AppAuc <- Auc.default(object=AppRoc$Sensitivity,Spec=AppRoc$Specificity)
    AppBS <- Brier.default(object=pred,y=Y,crRatio=crRatio)

  # }}}
# {{{ No information error  

    if (Plan$internal.name %in% c("boot632plus","noinf")){
      if (noinf.method=="simulate"){
        if (verbose)
          cat("\nSimulate no information performance\n")
        compute.NoInfRocList <- lapply(1:B,function(b){
          if (verbose) MgTalk(b,B)
          data.b <- data
          ## permute the response variable
          responseName <- all.vars(formula)[1]
          data.b[,responseName] <- sample(factor(Y),replace=FALSE)
          if (simulate=="reeval")
            fit.b <- MgRefit(object=fit,data=data.b,step=b,silent=na.accept>0,verbose=verbose)
          ## evaluate the model in data with reeallocated responses
          else
            fit.b <- fit
          pred.b <- do.call("predictStatusProb",
                            c(list(object=fit.b,newdata=data.b),
                              model.args[[f]]))
          innerNoInfRoc <- Roc.default(object=pred.b,y=data.b[,responseName],breaks=breaks,crRatio=crRatio)
          innerNoInfBS <- Brier.default(object=pred.b,y=data.b[,responseName],crRatio=crRatio)
          list("innerNoInfRoc"=innerNoInfRoc,"innerNoInfBS"=innerNoInfBS)
        })
        if (verbose) cat("\n")
        NoInfRocList <- lapply(compute.NoInfRocList,function(x)x$innerNoInfRoc)
        NoInfRoc <- avRoc(list=NoInfRocList,grid=RocAverageGrid,method=RocAverageMethod)
        NoInfBS <- mean(sapply(compute.NoInfRocList,function(x){x$innerNoInfBS}))
        NoInfAuc <- mean(sapply(NoInfRocList,function(nil){Auc.default(object=nil$Sensitivity,nil$Specificity)}))
      }
      else{         
        NoInfRoc <- list(Sensitivity=c(breaks,0),Specificity=c(1-breaks,1))
        NoInfAuc <- 0.5
        NoInfBS <- .C("brier_noinf",bs=double(1),as.double(Y),as.double(pred),as.integer(N),NAOK=TRUE,PACKAGE="ModelGood")$bs
      }
    }
    if (Plan$internal.name %in% c("boot632plus","bootcv","boot632")){

      # }}}
      # {{{ Bootcv aka BootstrapCrossValidation
      if (verbose)
        cat("\nBootstrap cross-validation performance\n")
      compute.BootcvRocList <- lapply(1:B,function(b){
        if (verbose) MgTalk(b,B)
        vindex.b <- match(1:N,CrossvalIndex[,b],nomatch=0)==0
        val.b <- data[vindex.b,,drop=FALSE]
        train.b <- data[CrossvalIndex[,b],,drop=FALSE]
        fit.b <- MgRefit(object=fit,data=train.b,step=b,silent=na.accept>0,verbose=verbose)
        if (!is.null(extract)) {
          fit.parms.b <- fit.b[extract]
          names(fit.parms.b) <- paste(extract,paste("sample",b,sep="."),sep=":")
        }
        else fit.parms.b <- NULL
        if (is.null(fit.b)){
          failed <- "fit"
          innerBootcvRoc <- list(Sensitivity=NA,Specificity=NA)
          innerBCVBS <- NA
        }
        else{
          try2predict <- try(pred.b <- do.call("predictStatusProb",c(list(object=fit.b,newdata=val.b),model.args[[f]])),silent=na.accept>0)
          if (inherits(try2predict,"try-error")==TRUE){
            if (verbose) warning(paste("During bootstrapping: prediction for model ",class(fit.b)," failed in step ",b),immediate.=TRUE)
            failed <- "prediction"
            innerBootcvRoc <- list(Sensitivity=NA,Specificity=NA)
            innerBCVBS <- NA
          }
          else{
            failed <- NA
            innerBootcvRoc <- Roc.default(y=Y[vindex.b],pred.b,breaks=breaks,crRatio=crRatio)
            innerBCVBS <- Brier.default(object=pred.b,y=Y[vindex.b],crRatio=crRatio)
          }
        }
        list("innerBootcvRoc"=innerBootcvRoc,
             "fit.parms"=fit.parms.b,
             "failed"=failed,
             "innerBCVBS"=innerBCVBS)
      })
      if (verbose) cat("\n")
      if (!is.null(extract)) fitParms <- sapply(compute.BootcvRocList,function(x)x$fit.parms)
      failed <- na.omit(sapply(compute.BootcvRocList,function(x)x$failed))
      BootcvRocList <- lapply(compute.BootcvRocList,function(x)x$innerBootcvRoc)
      BootcvRoc <- avRoc(BootcvRocList,grid=RocAverageGrid,method=RocAverageMethod)
      BCVSens <- BootcvRoc$Sensitivity
      BCVSpec <- BootcvRoc$Specificity
      BCVAuc <- mean(sapply(BootcvRocList,function(ool){Auc.default(object=ool$Sensitivity,ool$Specificity)}))
      BCVBS <- mean(sapply(compute.BootcvRocList,function(x){x$innerBCVBS}))
    }

# }}}
# {{{ Bootstrap .632

    if (Plan$internal.name=="boot632"){
      B632Roc <- list(Sensitivity=.368 * AppRoc$Sensitivity + .632 * BCVSens,
                      Specificity=.368 * AppRoc$Specificity + .632 * BCVSpec)
      B632BS <- .368 * AppBS + .632 * BCVBS
      B632Auc <- .368 * AppAuc + .632 * BCVAuc
    }
  # }}}
# {{{ Bootstrap .632+
    if (Plan$internal.name=="boot632plus"){
      ## first we have to prepare the averaging
      if (RocAverageMethod=="vertical"){
        AppSens <- c(approx(AppRoc$Sensitivity,AppRoc$Specificity,xout=RocAverageGrid,yleft=0,yright=1,ties=median)$y,0)
        AppRoc <- list(Sensitivity=AppSens,Specificity=c(RocAverageGrid,1))
        AppAuc <- Auc.default(object=AppRoc$Sensitivity,Spec=AppRoc$Specificity)
        
        R632Plus <- MgFormule632(App=AppSens,BCV=BCVSens,NoInf=NoInfRoc$Sensitivity,SmallerBetter=FALSE)
        B632plusRoc <- list(Sensitivity=R632Plus$B632Plus,Specificity=c(RocAverageGrid,1))
      }
      else{ ## RocAverageMethod horizontal
        AppSpec <- c(approx(AppRoc$Sensitivity,AppRoc$Specificity,xout=RocAverageGrid,yleft=1,yright=0,ties=median)$y,1)
        AppRoc <- list(Sensitivity=c(RocAverageGrid,0),Specificity=AppSpec)
        R632Plus <- MgFormule632(App=AppSpec,BCV=BCVSpec,NoInf=NoInfRoc$Specificity,SmallerBetter=FALSE)
        B632plusRoc <- list(Sensitivity=c(RocAverageGrid,0),
                            Specificity=R632Plus$B632Plus)
      }
      ## add the real .632+ AUC
      B632PlusAuc <- MgFormule632(App=AppAuc,BCV=BCVAuc,NoInf=NoInfAuc,SmallerBetter=FALSE)
      B632PlusBS <- MgFormule632(App=AppBS,BCV=BCVBS,NoInf=NoInfBS,SmallerBetter=TRUE)
    }
    # }}}
# {{{ preparing the output
    out <- switch(Plan$internal.name,
                  "noPlan"=list("Roc"=AppRoc),
                  ## "plain"=list("Roc"=BootRoc,"AppRoc"=AppRoc),
                  "boot632"=list("Roc"=B632Roc,"AppRoc"=AppRoc,"BootcvRoc"=BootcvRoc),
                  "boot632plus"=list("Roc"=B632plusRoc,
                    "AppRoc"=AppRoc,
                    "BootcvRoc"=BootcvRoc,
                    "NoInfRoc"=NoInfRoc,
                    "weight"=R632Plus$weight,
                    "overfit"=R632Plus$overfit),
                  "bootcv"=list("Roc"=BootcvRoc,"AppRoc"=AppRoc),
                  "noinf"=list("AppRoc"=AppRoc,"NoInfRoc"=NoInfRoc))
    
    out$Auc <- switch(Plan$internal.name,
                      "noPlan"=list("Auc"=AppAuc),
                      ## "plain"=list("Auc"=BootAuc,"AppAuc"=AppAuc),
                      "boot632"= list("Auc"=B632Auc,"AppAuc"=AppAuc,"AucBCV"=BCVAuc),
                      "boot632plus"=list("Auc"=B632PlusAuc$B632Plus,
                        "AppAuc"=AppAuc,
                        "BootcvAuc"=BCVAuc,
                        "NoInfAuc"=NoInfAuc,
                        "weight"=B632PlusAuc$weight,
                        "overfit"=B632PlusAuc$overfit),
                      "bootcv"=list("Auc"=BCVAuc,"AppAuc"=AppAuc),
                      "noinf"=list("NoInfAuc"=NoInfAuc,"AppAuc"=AppAuc))
    out$Brier <- switch(Plan$internal.name,
                        "noPlan"=list("BS"=AppBS),
                        ## "plain"=list("BS"=BootBS,"AppBS"=AppBS),
                        "boot632"= list("BS"=B632BS,"AppBS"=AppBS,"BSBCV"=BCVBS),
                        "boot632plus"=list("BS"=B632PlusBS$B632Plus,
                          "AppBS"=AppBS,
                          "BootcvBS"=BCVBS,
                          "NoInfBS"=NoInfBS,
                          "weight"=B632PlusBS$weight,
                          "overfit"=B632PlusBS$overfit),
                        "bootcv"=list("BS"=BCVBS,"AppBS"=AppBS),
                        "noinf"=list("AppBS"=AppBS,"NoInfBS"=NoInfBS))

    ##     if (keepCrossValRes==TRUE && Plan$internal.name!="noPlan"){
    
    if (keepCrossValRes==TRUE && class(try(is.null(BootcvRocList),silent=TRUE))!="try-error"){
      if (Plan$internal.name!="noinf")
        out <- c(out,list("BootcvRocList"=BootcvRocList))
    }
    if (keepNoInfSimu==TRUE && Plan$internal.name!="noPlan"){
      if (Plan$internal.name!="noinf")
        out <- c(out,list("NoInfRocList"=NoInfRocList))
    }
    if (!is.null(extract)) out <- c(out,list("fitParms"=fitParms))
    if (na.accept>0) out <- c(out,list("failed"=failed))
    out
  })
  ## it might be that the first model has no extracted parameters
  ## but one of the other has
  if (length(model.parms)>0)
    names.lout <- c(names(list.out[[1]]),"fitParms")
  else
    names.lout <- names(list.out[[1]])
  out <- lapply(names.lout,function(w){
    e <- lapply(list.out,function(x){x[[w]]})
    names(e) <- names(object)
    e
  })
  names(out) <- names.lout
  if(keepModels==TRUE)
    outmodels <- object
  else if (keepModels=="Call"){
    outmodels <- lapply(object,function(o)o$call)
    names(outmodels) <- names(object)
  }
  else{
    outmodels <- names(object)
    names(outmodels) <- names(object)
  }
  out <- c(out,
           list(call=match.call(),
                Response=Y,
                models=outmodels,
                method=Plan,
                breaks=breaks,crRatio=crRatio))
  if (verbose) cat("\n")
  # }}}
  class(out) <- "Roc"
  out
}

