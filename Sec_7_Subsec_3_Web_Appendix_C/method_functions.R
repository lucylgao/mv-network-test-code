## @knitr methods

gtest <- function(adj, node, K1, K2, B) { 
  n <- nrow(adj)
  
  fit1 <- tryCatch({pl_est_com(as.matrix(adj)*1, K1)}, error=function(err) return(NULL))
  fit2 <- tryCatch({Mclust(node, G=K2, modelNames=c("EII"),
                           initialization=list(hcPairs=hc(node, modelName="EII")))},  
                   err=function(err) return(NULL)) 
  
  if(is.null(fit1)) { 
    return(list(stat=NA, pval=NA))
  }
  
  if(is.null(fit2)) { 
    return(list(stat=NA, pval=NA))
  }
  
  clustTable <- table(fit1$class, fit2$class)
  piEst <- clustTable/n 
  piEst1 <- rowSums(clustTable)/n 
  piEst2 <- colSums(clustTable)/n
  
  g2Stat <- sum(clustTable*log(t(t(piEst/piEst1)/piEst2)), na.rm=TRUE)
  
  g2Statp <- rep(0, B)
  for(b in 1:B) { 
    clustTablep <- table(fit1$class, fit2$class[sample(1:n, replace=F)])
    piEstp <- clustTablep/n 
    
    g2Statp[b] <- sum(clustTablep*log(t(t(piEstp/piEst1)/piEst2)), na.rm=TRUE)
  }
  
  pvalp <- mean(ifelse(g2Statp >= g2Stat, 1, 0))
  
  return(list(stat=g2Stat, pval=pvalp))
}



peel <- function(adj, node, K, B) { 
  N <- nrow(node)
  
  fit_FMM <- tryCatch({Mclust(node, G=K, modelNames=c("EII"),
                           initialization=list(hcPairs=hc(node, modelName="EII")))},  
                   err=function(err) return(NULL)) 
  
  if(is.null(fit_FMM)) { 
    return(list(stat=NA, pval=NA))
  }
  
  class <- fit_FMM$class
  n <- table(class) 
  
  bernEnt <- function(x){ 
    if(x == 1 || x == 0) { 
      return(0)  
    }  else { 
      return(-x*log2(x) - (1-x)*log2(1-x))  
    }
  }
  
  w <- matrix(0, K, K)
  ent <- matrix(0, K, K)
  for(k1 in 1:K) { 
    for(k2 in 1:K) { 
      if(k1 != k2) { 
        w[k1, k2] <- sum(adj[class == k1, class == k2])/(n[k1]*n[k2])
        ent[k1, k2] <- (n[k1]*n[k2])*bernEnt(w[k1, k2])
      } else { 
        w[k1, k1] <- sum(adj[class == k1, class == k1])/(n[k1]*(n[k1]-1))
        ent[k1,k1] <- n[k1]*(n[k1]-1)*bernEnt(w[k1, k1])
      }
    }  
  }
  
  entropy_stat <- 0.5*sum(ent)
  
  entropy_perm <- rep(0,B)
  for(b in 1:B) { 
      class_perm <- class[sample(1:N, replace=F)]
      
      w <- matrix(0, K, K)
      ent <- matrix(0, K, K)
      for(k1 in 1:K) { 
        for(k2 in 1:K) { 
          if(k1 != k2) { 
            w[k1, k2] <- sum(adj[class_perm == k1, class_perm == k2])/(n[k1]*n[k2])
            ent[k1, k2] <- (n[k1]*n[k2])*bernEnt(w[k1, k2])
          } else { 
            w[k1, k1] <- sum(adj[class_perm == k1, class_perm == k1])/(n[k1]*(n[k1]-1))
            ent[k1,k1] <- n[k1]*(n[k1]-1)*bernEnt(w[k1, k1])
          }
        }  
      }
      
      entropy_perm[b] <- 0.5*sum(ent)
  }
  
  pvalp <- mean(ifelse(entropy_stat >= entropy_perm, 1, 0))
  
  return(list(stat=entropy_stat, pval=pvalp))
  
}

pplr_test <- new_method(name="p2lrt", 
             label="Pseudo pseudolikelihood ratio test",
             settings=list(), 
             method = function(model, draw) {
               K <- model$K
               testResult <- test_indep_com_clust(
                 list(as.matrix(draw$dat[[1]])*1, draw$dat[[2]]), 
                  K1=K, K2=K, model2="EII", init2="EII", nperm=200) 
               
               return(list(stat=testResult$P2LRstat, 
                           pval=testResult$pval))
              })

g_test <- new_method("gtest", "G-test",
                     method = function(model, draw) {
                       K <- model$K
                       testResult <- gtest(draw$dat[[1]], draw$dat[[2]], 
                                           K, K, 200)
                       stat <- testResult$stat
                       pval <- testResult$pval
                       return(list(stat=stat, pval=pval))
                     })

peel_test <- new_method("peeltest", "BeST",
                     method = function(model, draw) {
                       K <- model$K
                       testResult <- peel(draw$dat[[1]], draw$dat[[2]], 
                                           K, 200)
                       stat <- testResult$stat
                       pval <- testResult$pval
                       return(list(stat=stat, pval=pval))
                     })

peel_test_spectral <- new_method("peeltest-spectral", "BeST",
                        method = function(model, draw) {
                          estK2 <- Mclust(draw$dat[[2]], G=2:9, modelNames=c("EII"),
                                          initialization=list(hcPairs=hc(draw$dat[[2]],
                                                                         modelName="EII")))$G
                          testResult <- peel(draw$dat[[1]], draw$dat[[2]], 
                                             estK2, 200)
                          stat <- testResult$stat
                          pval <- testResult$pval
                          return(list(stat=stat, pval=pval))
                        })

g_test_spectral <- new_method("gtest-spectral", "G-test-spectral",
                     method = function(model, draw) {
                       estK1 <- BHMC.estimate(draw$dat[[1]], 15)$K
                       estK2 <- Mclust(draw$dat[[2]], G=2:9, modelNames=c("EII"),
                                       initialization=list(hcPairs=hc(draw$dat[[2]],
                                                                      modelName="EII")))$G
                       testResult <- gtest(draw$dat[[1]], draw$dat[[2]], 
                                           max(2, estK1), estK2, 200)
                       stat <- testResult$stat
                       pval <- testResult$pval
                       return(list(stat=stat, pval=pval))
                     })

pplr_test_spectral <- new_method("p2lrt-spectral", "p2lrt-spectral",
                       method = function(model, draw) {
                         n <- nrow(draw$dat[[1]])
                         testResult <- test_indep_com_clust(
                           list(as.matrix(draw$dat[[1]])*1, draw$dat[[2]]), 
                           model2="EII", init2="EII", nperm=200)
                         stat <- testResult$P2LRstat
                         pval <- testResult$pval
                         return(list(stat=stat, pval=pval))
                       })
