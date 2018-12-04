#Today is July 25th, 2018
#The order() function just takes the (200x200) haplotype distance matrix and reorders the elements
#such that the columns and rows are ordered from haplotype 1 to haplotype 200.

#Example: Individual 1 has haplotype 1 and haplotype 2. Individual 4 has haplotype 8 and haplotype 9.


order = function(my.matrix){
  result = matrix(rep(0,200*200), nrow = 200) 
  for (i in 1:200){
    for (j in 1:200){
      result[i,j] = my.matrix[as.character(i), as.character(j)]
    }
  }  
  return(result)
}
