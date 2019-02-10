#Solution Function.
#This function truncates the (200x200) haplotype similarity matrix to a 
#(100x100) individual based similarity matrix.

#Find the maximum distance between the pairwise distances of each 2 haplotypes between all pairs of individuals.
Solution.function = function(tree.kernel.matrix){
  n=100
  result = matrix(rep(0,n*n), nrow = n)
  for (i in 1:(n -1)) {
    for (j in (i +1): n) {
      pair.1 = tree.kernel.matrix[2*i -1 ,2*j -1] + tree.kernel.matrix[2*i ,2*j]
      pair.2 = tree.kernel.matrix[2*i -1 ,2*j] + tree.kernel.matrix[2*i ,2*j -1]
      result[i,j]  = max(pair.1, pair.2)
      result[j,i] = max(pair.1, pair.2)
      
      diag.pair.1 =(tree.kernel.matrix[2*i-1,2*i-1] + tree.kernel.matrix[2*i,2*i])
      #diag.pair.2 = (tree.kernel.matrix[2*i-1,2*i] + tree.kernel.matrix[2*i,2*i-1])
      #result[i,i]= max(diag.pair.1, diag.pair.2)
      result[i,i] = diag.pair.1
      
      
    }
  }
  last.diag.pair1 = (tree.kernel.matrix[199,199] + tree.kernel.matrix[200,200])
  #last.diag.pair2 = (tree.kernel.matrix[199,200] + tree.kernel.matrix[200,199])
  #result[100,100] = min(last.diag.pair1, last.diag.pair2)
  result[100,100] = last.diag.pair1
  
  return(result)
  
}

