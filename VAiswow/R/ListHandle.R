#' Combine two lists
#'
#' Combine two lists of vectors or matrices
#'
#' @param x A list of vector and/or matrices
#' @param y Another list of same dimension as x, with elements matching those of x
#' 
#' @return None
#'
#' @examples
#' x <- list(matrix(0,10,10), matrix(1,10,10), rnorm(100))
#' y <- list(matrix(2,10,10), matrix(3,10,10), rnorm(100))
#' meltlists(x=x, y=y)
#'
#' @export
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
      }else{
        stop("Element ",i," was not a matrix nor a vector")
      }
    }
  }
  
  return(melted)
}#end meltlists()