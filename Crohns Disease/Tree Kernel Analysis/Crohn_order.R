
order = function(my.matrix, n){
  result = matrix(rep(0,n*n), nrow = n) 
  for (i in 1:n){
    for (j in 1:n){
      result[i,j] = my.matrix[as.character(i), as.character(j)]
    }
  }  
  return(result)
}