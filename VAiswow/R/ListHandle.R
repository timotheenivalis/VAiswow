meltlists <- function(x,y)
{
  stopifnot(exprs = {
      is.list(x)
      is.list(y)
      length(x)==length(y)} )
  
  melted <- x
  for (i in 1:length(x))
  {
    if(is.vector(melted[[i]]))
      {
      melted[[i]] <- c(melted[[i]], y[[i]])
    }else{
      if(is.matrix(melted[[i]]))
      {
        if(ncol(melted[[i]]) != ncol(y[[i]]))
        {
         stop("Matrices ncol differ at index ", i) 
        }
        melted[[i]] <- rbind(melted[[i]], y[[i]])
      }
    }
    
  }
  
  return(melted)
}#end meltlists()