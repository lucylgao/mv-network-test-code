## @knitr methods
gtest <- function(adj1, adj2, K1, K2, B) { 
  n <- nrow(adj1)
  
  fit1 <- tryCatch({pl_est_com(as.matrix(adj1)*1, K1)}, error=function(err) return(NULL))
  fit2 <- tryCatch({pl_est_com(as.matrix(adj2)*1, K2)}, error=function(err) return(NULL))
  
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

pplr_test <- new_method(name="p2lrt", 
             label="Pseudo pseudolikelihood ratio test",
             settings=list(), 
             method = function(model, draw) {
               K <- model$K
               testResult <- test_indep_com(X=list(as.matrix(draw$adj[[1]])*1, 
                                                 as.matrix(draw$adj[[2]])*1), 
                                            K1=K, K2=K, nperm=200)
               
               return(list(stat=testResult$P2LRstat, 
                           pval=testResult$pval))
             })

g_test <- new_method("gtest", "G-test",
                     method = function(model, draw) {
                       K <- model$K
                       testResult <- gtest(draw$adj[[1]], draw$adj[[2]], 
                                           K, K, 200)
                       stat <- testResult$stat
                       pval <- testResult$pval
                       return(list(stat=stat, pval=pval))
                     })

g_test_spectral <- new_method("gtest-spectral", "G-test-spectral",
                     method = function(model, draw) {
                       estK1 <- BHMC.estimate(draw$adj[[1]], 100)$K
                       estK2 <- BHMC.estimate(draw$adj[[2]], 100)$K
                       testResult <- gtest(draw$adj[[1]], draw$adj[[2]], 
                                           max(2, estK1), max(estK2, 2), 200)
                       stat <- testResult$stat
                       pval <- testResult$pval
                       return(list(stat=stat, pval=pval))
                     })

pplr_test_spectral <- new_method("p2lrt-spectral", "p2lrt-spectral",
                       method = function(model, draw) {
                         testResult <- test_indep_com(X=list(as.matrix(draw$adj[[1]])*1, 
                                                             as.matrix(draw$adj[[2]])*1), 
                                                      nperm=200)
                         stat <- testResult$P2LRstat
                         pval <- testResult$pval
                         return(list(stat=stat, pval=pval))
                       })
