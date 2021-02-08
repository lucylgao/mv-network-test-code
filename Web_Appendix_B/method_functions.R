## @knitr methods

make_p2lrt <- function(algK) { 
  new_method(paste("p2lrt-K", algK, sep=""), paste("p2lrt-K", algK, sep=""), 
             method=function(model, draw) { 
               out <- test_indep_com(list(as.matrix(draw$adj[[1]])*1, 
                                          as.matrix(draw$adj[[2]])*1), 
                                     K1=algK, K2=algK) 
               list(stat=out$P2LRstat, pval=out$pval)
             })
}

p2lrt_vect <- list()

for(j in 2:50) { 
  p2lrt_vect[[j - 1]] <- make_p2lrt(j)  
}

p2lrt_estK <- new_method("p2lrt-estK", "p2lrt-estK", 
             method=function(model, draw) { 
               out <- test_indep_com(list(as.matrix(draw$adj[[1]])*1, 
                                          as.matrix(draw$adj[[2]])*1)) 
               Kused <- dim(out$Pi.est)[1]
               list(stat=out$P2LRstat, pval=out$pval)
})

