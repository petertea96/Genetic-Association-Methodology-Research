###############################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                   treeSimilarity.R MODIFIED VERSION....          
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################################################
#Edited: July 18th, 2018.


## This Function gets the tree-based similarity
## Assumes that the tree is in format from APE package
#Same file that Kelly sent me, except that I added/changed 2 things.
#.libPaths("/global/home/hpc4300/RPackages")
treeSimilarity=function(tree,method="branchlength",propn=FALSE){
  
  require("ape")
  
  tmrca=coalescent.intervals(tree)$total.depth
  
  if (method=="branchlength")  {
    
    # Similarity is the time back to their first parent
    # scaled by the time back to the mrca of the sample
    
    #The argument used to read "mytree", I changed it to "tree"...
    dist.raw=cophenetic.phylo(tree)/2
    x=tmrca-dist.raw
    
    if (propn==TRUE){
      x=x/tmrca
    }
  } else if (method=="ranks"){
    
    numtips=length(tree$tip.label)
    nodes=as.numeric(names(sort(branching.times(tree), decreasing=T)))
    
    # This gives us the current internal node labels. We need
    # to change these using the ordering of the sorted nodes
    x=mrca(tree)
    newx=x
    
    # Recode the matrix x so that instead of giving the
    # node id of the common ancestor, it gives the 
    # number corresponding to which coalescent event this is
    for (i in 1:length(nodes)){
      newx[x==nodes[i]]=i
    }
    
    diag(newx)=numtips
    x=newx
    
  } else  if (method=="sdscaled"){
    
    #***I added in this code***:
    ##############################
    numtips=length(tree$tip.label)
    ##############################
    
    # Get node labels for internal nodes
    nodes=as.numeric(names(sort(branching.times(tree))))
    
    # Get new branching times
    times=coalescent.intervals(tree)$interval.length
    sds=1/(choose(numtips:2,2))
    newtimes=times/sds
    newtime.topresent=newtimes
    
    for (i in 2:length(newtimes)){
      newtime.topresent[i]=sum(newtimes[1:i])
    }
    
    names(newtime.topresent)=nodes
    
    # This gives us the current internal node labels. We need
    # to change these using the ordering of the sorted nodes
    x=mrca(tree)
    newx=x
    
    for (i in 1:length(newtime.topresent)){
      index=names(newtime.topresent)[i]
      newx[x==index]=newtime.topresent[i]
    }
    
    tmrca=newtime.topresent[length(newtime.topresent)]
    x=tmrca-newx
    diag(x)=tmrca
    
    if (propn==TRUE){
      x=x/tmrca
    }
    
  }
  
  return(x)
  
}

## Example using the function above
#mytree=read.tree(text="((5:0.5330782044,9:0.5330782044):2.193516845,((1:0.07665280018,3:0.07665280018):1.092632651,(((6:0.2006021859,2:0.2006021859):0.2398099728,((7:0.001004649848,8:0.001004649848):0.06105826045,10:0.0620629103):0.3783492484):0.08891256889,4:0.5293247276):0.6399607239):1.557309598);")
#plot(mytree)
#nodelabels()


# Various similarity matrices 
#sim1=round(treeSimilarity(mytree),3) # TMRCA - time to present of pairs mrca
#sim2=round(treeSimilarity(mytree,propn=T), 3) # (TMRCA - time to present of pairs mrca)/TMRCA
#sim3=round(treeSimilarity(mytree, method="ranks"),3) # Rank of coalescent event for each pairs mrca (oldest ranked 1)
#sim4=round(treeSimilarity(mytree, method="sdscaled"),3) # Like 1, but intercoalescence times divided by mean/sd 
#sim5=round(treeSimilarity(mytree, method="sdscaled",propn=T),3) # Like 2, but intercoalescence times divided by mean/sd 


#ntips=100
#mytree=rcoal(ntips)



