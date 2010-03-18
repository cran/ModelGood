.First.lib <- function(lib,pkg) 
  library.dynam("ModelGood",pkg,lib)

.Last.lib <- function(lib)
  library.dynam.unload("ModelGood",lib)
