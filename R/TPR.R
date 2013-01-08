# {{{ UseMethod

TPR <- function(x,...){
  UseMethod("TPR",x)
}

TPR.numeric <- function(x,...){
  ci.tpr=binconf(x[4],(x[3]+x[4]),...)
  colnames(ci.tpr)[1]="TPR (Sensitivity)"
  ci.tpr
}

TPR.table <- function(x,...){
  ci.tpr=binconf(x[4],(x[3]+x[4]),...)
  colnames(ci.tpr)[1]="TPR (Sensitivity)"
  ci.tpr
}

TNR <- function(x,...){
  UseMethod("TNR",x)
}

TNR.table <- function(x,...){
  ci.tnr=binconf(x[1],(x[1]+x[2]),...)
  colnames(ci.tnr)[1]="TNR (Specificity)"
  ci.tnr
}

# }}}
