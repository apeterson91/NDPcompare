
#' Generate Adjacency Matrix
#'
generate_adjacency_matrix <- function(J){

  A <- matrix(0,ncol = sum(J),nrow = sum(J))
  for(k in 1:length(J)){
    if(k==1)
      start <- 1
    else
      start <- sum(J[1:k-1])+1
    stop <- sum(J[1:k])
    for(i in start:stop){
      for(j in start:i)
        A[i,j] <- 1
    }
  }
  diag(A) <- 0

  return(A)
}
