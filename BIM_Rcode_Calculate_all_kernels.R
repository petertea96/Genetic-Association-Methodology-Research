#This R script calculate all kernel matrices...
#Today is July 26th, 2018


.libPaths("/global/home/hpc4300/RPackages")
#Set directory where all R packages are in CAC cluster
library("SKAT")
require("ape")
require("MDMR")

#Source some R codes needed to calculate the kernels...


source("BIM_Rcode_treesimilarityMODIFIED.R") # Calculate distance matrices
source("BIM_Rcode_SKAT_Linear.R") #SKAT package
source("BIM_Rcode_Function.R") # From SKAT package
source("BIM_Rcode_ORDER.R") # Order the matrix from 1 - 200.
source("BIM_Rcode_Solution_function.R") #Get the max. distance for each individual.


get.kernels = function(G, P1, P2, n, K, treename){
  #Returns list of all 16 kernels for the ith dataset

  
  
  ##---- Similarity Based on Identity -by - State Allele Sharing ----##
  ## Denote S_ IBS_k be the n by n similarity matrix at variant k
  S_IBS_k <-list ()
  for (k in 1:K) {
    S_IBS_k[[k]] <- matrix (rep (0,n*n), nrow =n)
    for (i in 1:n) {
      for (r in i:n) {
        if (G[i,k] == G[r,k]) {
          S_IBS_k[[k]][i,r] = S_IBS_k[[k]][r,i] = 2
          
          
        } else if (G[i,k] == 0 && G[r,k] == 1) {
          S_IBS_k[[k]][i,r] = S_IBS_k[[k]][r,i] = 1
        } else if (G[i,k] == 1 && G[r,k] == 0) {
          S_IBS_k[[k]][i,r] = S_IBS_k[[k]][r,i] = 1
        } else if (G[i,k] == 1 && G[r,k] == 2) {
          S_IBS_k[[k]][i,r] = S_IBS_k[[k]][r,i] = 1
        } else if (G[i,k] == 2 && G[r,k] == 1) {
          S_IBS_k[[k]][i,r] = S_IBS_k[[k]][r,i] = 1
        } else {
          S_IBS_k[[k]][i,r] = S_IBS_k[[k]][r,i] = 0
        }
      }
    }
  }
  
  ## Denote S_ IBS be the (n by n) similarity matrix over K variants
  sim <- matrix (c(rep (0,n*n)), nrow =n)
  for (t in 1:K) {
    sim = sim +S_IBS_k[[t]]
  }
  S_IBS <- sim/(2*K)
  
  
  
  
  
  ##################################################################################################
  ##--- Allele Match Kernel -- -##
  S_AM_k <-list ()
  for (k in 1:K) {
    S_AM_k[[k]] <- matrix (rep (0,n*n), nrow =n)
    for (i in 1:n) {
      for (r in i:n) {
        if (G[i,k] == G[r,k]) {
          S_AM_k[[k]][i,r] = S_AM_k[[k]][r,i] = 4
        } else if (G[i,k] == 0 && G[r,k] == 1) {
          S_AM_k[[k]][i,r] = S_AM_k[[k]][r,i] = 2
        } else if (G[i,k] == 1 && G[r,k] == 0) {
          S_AM_k[[k]][i,r] = S_AM_k[[k]][r,i] = 2
        } else if (G[i,k] == 1 && G[r,k] == 2) {
          S_AM_k[[k]][i,r] = S_AM_k[[k]][r,i] = 2
        } else if (G[i,k] == 2 && G[r,k] == 1) {
          S_AM_k[[k]][i,r] = S_AM_k[[k]][r,i] = 2
        } else {
          S_AM_k[[k]][i,r] = S_AM_k[[k]][r,i] = 0
        }
      }
    }
  }
  
  ## Denote S_AM be the n by n similarity matrix over K variants
  sim_AM <- matrix (c(rep (0,n*n)), nrow =n)
  for (t in 1:K) {
    sim_AM = sim_AM+S_AM_k[[t]]
  }
  S_AM <- sim_AM/(2*K) ##S_AM = 2*S_ IBS
  
  
  
  #################################################################################################
  
  ##--- Allele Share Kernel -- -##
  S_AS_k <-list ()
  for (k in 1:K) {
    S_AS_k[[k]] <- matrix (rep (0,n*n), nrow =n)
    for (i in 1:n) {
      for (r in i:n) {
        if (G[i,k] == 1 && G[r,k] == 1) {
          S_AS_k[[k]][i,r] = S_AS_k[[k]][r,i] = 1
        } else if (G[i,k] == 1 && G[r,k] == 2) {
          S_AS_k[[k]][i,r] = S_AS_k[[k]][r,i] = 1
        } else if (G[i,k] == 2 && G[r,k] == 1) {
          S_AS_k[[k]][i,r] = S_AS_k[[k]][r,i] = 1
        } else if (G[i,k] == 2 && G[r,k] == 2) {
          S_AS_k[[k]][i,r] = S_AS_k[[k]][r,i] = 2
        } else {
          S_AS_k[[k]][i,r] = S_AS_k[[k]][r,i] = 0
        }
      }
    }
  }
  ## Denote S_AS be the n by n similarity matrix over K variants
  sim_AS <- matrix (c(rep (0,n*n)), nrow =n)
  for (t in 1:K) {
    sim_AS = sim_AS+S_AS_k[[t]]
  }
  S_AS <- sim_AS/(2*K)
  
  #####################################################################################################
  
  ##--- Additive kernel (LIN )---##
  S_LIN_k <-list ()
  for (k in 1:K) {
    S_LIN_k[[k]] <- matrix (rep (0,n*n), nrow =n)
    for (i in 1:n) {
      for (r in i:n) {
        S_LIN_k[[k]][i,r] = S_LIN_k[[k]][r,i] = G[i,k]+G[r,k]
      }
    }
  }
  ## Denote S_ LIN be the n by n similarity matrix over K variants
  sim_LIN <- matrix (c(rep (0,n*n)), nrow =n)
  for (t in 1:K) {
    sim_LIN = sim_LIN +S_LIN_k[[t]]
  }
  S_LIN <- sim_LIN/(2*K)
  
  
  
  
  
  #####################################################################################################
  ##--- Additive kernel (REC )---##
  S_REC_k <-list ()
  for (k in 1:K) {
    S_REC_k[[k]] <- matrix (rep (0,n*n), nrow =n)
    for (i in 1:n) {
      for (r in i:n) {
        if (G[i,k] == 1 && G[r,k] == 2) {
          S_REC_k[[k]][i,r] = S_REC_k[[k]][r,i] = 1
        } else if (G[i,k] == 2 && G[r,k] == 1) {
          S_REC_k[[k]][i,r] = S_REC_k[[k]][r,i] = 1
        } else if (G[i,k] == 2 && G[r,k] == 0) {
          S_REC_k[[k]][i,r] = S_REC_k[[k]][r,i] = 1
        } else if (G[i,k] == 0 && G[r,k] == 2) {
          S_REC_k[[k]][i,r] = S_REC_k[[k]][r,i] = 1
        } else if (G[i,k] == 2 && G[r,k] == 2) {
          S_REC_k[[k]][i,r] = S_REC_k[[k]][r,i] = 2
        } else {
          S_REC_k[[k]][i,r] = S_REC_k[[k]][r,i] = 0
        }
      }
    }
  }
  ## Denote S_ REC be the n by n similarity matrix over K variants
  sim_REC <- matrix (c(rep (0,n*n)), nrow =n)
  for (t in 1:K) {
    sim_REC = sim_REC +S_REC_k[[t]]
  }
  S_REC <- sim_REC/(2*K)
  
  
  
  ####################################################################################################
  ##--- Additive kernel ( QUAD )---##
  S_QUAD_k <-list ()
  for (k in 1:K) {
    S_QUAD_k[[k]] <- matrix (rep (0,n*n), nrow =n)
    for (i in 1:n) {
      for (r in i:n) {
        if (G[i,k] == 0 &&G[r,k] == 0) {
          S_QUAD_k[[k]][i,r] = S_QUAD_k[[k]][r,i] = 2
        } else if (G[i,k] == 1 && G[r,k] == 1) {
          S_QUAD_k[[k]][i,r] = S_QUAD_k[[k]][r,i] = 4
        } else if (G[i,k] == 2 && G[r,k] == 2) {
          S_QUAD_k[[k]][i,r] = S_QUAD_k[[k]][r,i] = 8
        } else if (G[i,k] == 1 && G[r,k] == 0) {
          S_QUAD_k[[k]][i,r] = S_QUAD_k[[k]][r,i] = 3
        } else if (G[i,k] == 0 && G[r,k] == 1) {
          S_QUAD_k[[k]][i,r] = S_QUAD_k[[k]][r,i] = 3
        } else if (G[i,k] == 1 && G[r,k] == 2) {
          S_QUAD_k[[k]][i,r] = S_QUAD_k[[k]][r,i] = 6
        } else if (G[i,k] == 2 && G[r,k] == 1) {
          S_QUAD_k[[k]][i,r] = S_QUAD_k[[k]][r,i] = 6
        } else {
          S_QUAD_k[[k]][i,r] = S_QUAD_k[[k]][r,i] = 5
        }
      }
    }
  }
  ## Denote S_ QUAD be the n by n similarity matrix over K variants
  sim_QUAD <- matrix (c(rep (0,n*n)), nrow =n)
  for (t in 1:K) {
    sim_QUAD = sim_QUAD +S_QUAD_k[[t]]
  }
  S_QUAD <- sim_QUAD /(2*K)
  
  
  
  ###################################################################################################
  ##--- Product kernel (0.1.2) - - -##
  S_012_k <-list ()
  for (k in 1:K) {
    S_012_k[[k]] <- matrix (rep (0,n*n), nrow =n)
    for (i in 1:n) {
      for (r in i:n) {
        S_012_k[[k]][i,r] = S_012_k[[k]][r,i] = G[i,k]*G[r,k]
      }
    }
  }
  ## Denote S_ 012 be the n by n similarity matrix over K variants
  sim_012 <- matrix (c(rep (0,n*n)), nrow =n)
  for (t in 1:K) {
    sim_012 = sim_012+ S_012_k[[t]]
  }
  S_012 <- sim_012/(2*K)
  
  
  #####################################################################################################
  ##--- Product kernel (1.2.3) - - -##
  S_123_k <-list ()
  for (k in 1:K) {
    S_123_k[[k]] <- matrix (rep (0,n*n), nrow =n)
    for (i in 1:n) {
      for (r in i:n) {
        if (G[i,k] == 0 && G[r,k] == 0) {
          S_123_k[[k]][i,r] = S_123_k[[k]][r,i] = 1
        } else if (G[i,k] == 1 && G[r,k] == 1) {
          S_123_k[[k]][i,r] = S_123_k[[k]][r,i] = 4
        } else if (G[i,k] == 2 && G[r,k] == 2) {
          S_123_k[[k]][i,r] = S_123_k[[k]][r,i] = 9
        } else if (G[i,k] == 1 && G[r,k] == 0) {
          S_123_k[[k]][i,r] = S_123_k[[k]][r,i] = 2
        } else if (G[i,k] == 0 && G[r,k] == 1) {
          S_123_k[[k]][i,r] = S_123_k[[k]][r,i] = 2
        } else if (G[i,k] == 1 && G[r,k] == 2) {
          S_123_k[[k]][i,r] = S_123_k[[k]][r,i] = 6
        } else if (G[i,k] == 2 && G[r,k] == 1) {
          S_123_k[[k]][i,r] = S_123_k[[k]][r,i] = 6
        } else {
          S_123_k[[k]][i,r] = S_123_k[[k]][r,i] = 3
        }
      }
    }
  }
  ## Denote S_ 123 be the n by n similarity matrix over K variants
  sim_123 <- matrix (c(rep (0,n*n)), nrow =n)
  for (t in 1:K) {
    sim_123 = sim_123+ S_123_k[[t]]
  }
  S_123 <- sim_123/(2*K)
  
  ####################################################################################################
  
  ##--- Product kernel (1.2.4) - - -##
  S_124_k <-list ()
  for (k in 1:K) {
    S_124_k[[k]] <- matrix (rep (0,n*n), nrow =n)
    for (i in 1:n) {
      for (r in i:n) {
        if (G[i,k] == 0 && G[r,k] == 0) {
          S_124_k[[k]][i,r] = S_124_k[[k]][r,i] = 1
        } else if (G[i,k] == 1 && G[r,k] == 1) {
          S_124_k[[k]][i,r] = S_124_k[[k]][r,i] = 4
        } else if (G[i,k] == 2 && G[r,k] == 2) {
          S_124_k[[k]][i,r] = S_124_k[[k]][r,i] = 16
        } else if (G[i,k] == 1 && G[r,k] == 0) {
          S_124_k[[k]][i,r] = S_124_k[[k]][r,i] = 2
        } else if (G[i,k] == 0 && G[r,k] == 1) {
          S_124_k[[k]][i,r] = S_124_k[[k]][r,i] = 2
        } else if (G[i,k] == 1 && G[r,k] == 2) {
          S_124_k[[k]][i,r] = S_124_k[[k]][r,i] = 8
        } else if (G[i,k] == 2 && G[r,k] == 1) {
          S_124_k[[k]][i,r] = S_124_k[[k]][r,i] = 8
        } else {
          S_124_k[[k]][i,r] = S_124_k[[k]][r,i] = 4
        }
      }
    }
  }
  ## Denote S_ 124 be the n by n similarity matrix over K variants
  sim_124 <- matrix (c(rep (0,n*n)), nrow =n)
  for (t in 1:K) {
    sim_124 = sim_124+ S_124_k[[t]]
  }
  S_124 <- sim_124/(2*K)
  
  
  ##--- Similarity Using Average Allelic Sharing -- -##
  S_h1_k <-list ()
  for (k in 1:K) {
    S_h1_k[[k]] <- matrix (rep (0,n*n), nrow =n)
    for (i in 1:n) {
      for (r in i:n) {
        if (G[i,k] == 0 && G[r,k] == 0) {
          S_h1_k[[k]][i,r] = S_h1_k[[k]][r,i] = 4
        } else if (G[i,k] ==1 && G[r,k] == 1) {
          S_h1_k[[k]][i,r] = S_h1_k[[k]][r,i] = 2
        } else if (G[i,k] == 2 && G[r,k] == 2) {
          S_h1_k[[k]][i,r] = S_h1_k[[k]][r,i] = 4
        } else if (G[i,k] == 0 && G[r,k] == 1) {
          S_h1_k[[k]][i,r] = S_h1_k[[k]][r,i] = 2
        } else if (G[i,k] == 1 && G[r,k] == 0) {
          S_h1_k[[k]][i,r] = S_h1_k[[k]][r,i] = 2
        } else if (G[i,k] == 1 && G[r,k] == 2) {
          S_h1_k[[k]][i,r] = S_h1_k[[k]][r,i] = 2
        } else if (G[i,k] == 2 && G[r,k] == 1) {
          S_h1_k[[k]][i,r] = S_h1_k[[k]][r,i] = 2
        } else {
          S_h1_k[[k]][i,r] = S_h1_k[[k]][r,i] = 0
        }
      }
    }
  }
  ## Denote S_h1 be the n by n similarity matrix over K variants
  sim_h1 <- matrix (c(rep (0,n*n)), nrow =n)
  for (t in 1:K) {
    sim_h1 = sim_h1+S_h1_k[[t]]
  }
  S_h1 <- sim_h1/K
  
  
  #SKAT kernel:
  MAF = Get_MAF(G) #Must source Function.R
  weights.beta=c(1,25)
  weights<-Beta.Weights(MAF,weights.beta)
  
  skat.kernel = t(t(G)*weights)
  skat.kernel = skat.kernel%*%t(skat.kernel)
  
  
  require(ape)

  mytree=read.tree(treename)
  #mytree$tip.label
  
  #Save the tree as an image to see the tip labels...
  #png("plot.png", width=25, height=25, units="in", res=500)
  #plot(mytree) 
  #dev.off()
  
  # Various (tree) similarity matrices 
  sim1=treeSimilarity(mytree) # TMRCA - time to present of pairs mrca
  sim2=treeSimilarity(mytree,propn=TRUE) # (TMRCA - time to present of pairs mrca)/TMRCA
  sim3=treeSimilarity(mytree, method="ranks") # Rank of coalescent event for each pairs mrca (oldest ranked 1)
  sim4=treeSimilarity(mytree, method="sdscaled") # Like 1, but intercoalescence times divided by mean/sd 
  sim5=treeSimilarity(mytree, method="sdscaled",propn=T) # Like 2, but intercoalescence times divided by mean/sd 
  
  #Sorted issue:
  a1=order(sim1)
  a2=order(sim2)
  a3=order(sim3)
  a4=order(sim4)
  a5=order(sim5)
  
  #Find the max distance for all pair of individuals
  tree1 = Solution.function(a1)
  tree2 = Solution.function(a2)
  tree3 = Solution.function(a3)
  tree4 = Solution.function(a4)
  tree5 = Solution.function(a5)
  
  
  
  
  #Kernel Matrices
  kernel.list = list(S_IBS, S_AM, S_AS, S_LIN, S_REC, S_QUAD, S_012, S_123, S_124, S_h1, skat.kernel,
                     tree1, tree2, tree3, tree4, tree5)
  names(kernel.list) = c("S_IBS", "S_AM", "S_AS", "S_LIN", "S_REC", "S_QUAD", "S_012", "S_123", "S_124",
                         "S_h1", "skat.kernel","tree1", "tree2", "tree3", "tree4", "tree5")
  return(kernel.list)
}


