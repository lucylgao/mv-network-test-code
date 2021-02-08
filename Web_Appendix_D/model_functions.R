## @knitr models

# generates cluster memberships given a Pi matrix
generate_memberships <- function(n, Pi) {
  stopifnot(Pi >= 0, abs(sum(Pi) - 1) < 1e-15)
  ii <- sample(length(Pi), n, replace = TRUE, prob = Pi)
  K1 <- nrow(Pi)
  K2 <- ncol(Pi)
  cpairs <- cbind(rep(1:K1, K2), rep(1:K2, each = K1))
  cpairs[ii, ]
}

# output: random adjacency matrix X with E[X] = pmat 
generate_random_graph <- function(pmat) { 
  n <- nrow(pmat)
  upper.ind <- which(upper.tri(pmat), arr.ind=T)
  upper.unif <- runif(n = nrow(upper.ind))
  edge.ind <- upper.ind[which(upper.unif < pmat[upper.ind]), ]
  sparseMatrix(edge.ind[, 1], edge.ind[, 2], dims=c(n, n), symmetric=T)
}

make_dcsbm_prob_mat <- function(B, pop, cl) { 
  n <- length(cl) 
  K <- nrow(B)
  
  # expected value of adjacency matrix: 
  # n x n matrix with ijth entry = pop[i] pop[j] B[cl[i], cl[j]]
  outer <- matrix(0, n, K)
  outer[cbind(1:n, cl)] <- 1
  outer <- pop*outer
  
  outer%*%tcrossprod(B, outer)
}

make_Pi <- function(delta) {
  probs <- c(0.05, 0.05, 0.15, 0.15, 0.3, 0.3)
  
  Pi1 <- tcrossprod(probs, probs)
  Pi2 <- rbind(c(0.02, 0, 0, 0, 0.015, 0.015), c(0, 0.02, 0, 0, 0.015, 0.015), 
               c(0, 0, 0.06, 0, 0.045, 0.045), c(0, 0, 0, 0.06, 0.045, 0.045), 
               c(0.015, 0.015, 0.045, 0.045, 0.09, 0.09), 
               c(0.015, 0.015, 0.045, 0.045, 0.09, 0.09))
  
  return((1-delta)*Pi1 + delta*Pi2)
} 

# B has one number on the diagonal and one number on the off-diagonal
# lambda = expected node degree 
# r = ratio of expected number of within-community edges to 
# expected number of between-community edges 
# s = sparsity
# K = number of subgroups 
make_B <- function(r, s, K) { 
  gamma <- (K/(K-1+2*r))*s
  B <- matrix(gamma, K, K)
  diag(B) <- 2*gamma*r
  
  return(B)
}

generate_mv_dcsbm <- function(n, K, r, s, delta, hub) { 
  B1 <- make_B(r, s, K) #assortative 
  B2 <- matrix(B1[1,1], K, K) # not assortative
  diag(B2) <- B1[1, 2]
  Pi <- make_Pi(delta)
  cl <- generate_memberships(n, Pi)
  
  if(hub) { 
    pop1 <- sample(c(0.625, 0.625*4), n, replace=T, prob=c(0.8, 0.2))
    pop2 <- sample(c(0.625, 0.625*4), n, replace=T, prob=c(0.8, 0.2))
  } else { 
    pop1 <- rep(1, n)
    pop2 <- pop1
  }
  pmat1 <- make_dcsbm_prob_mat(B1, pop1, cl[, 1])
  pmat2 <- make_dcsbm_prob_mat(B2, pop2, cl[, 2])
  adj1 <- generate_random_graph(pmat1)
  adj2 <- generate_random_graph(pmat2)
  
  return(list(adj=list(adj1, adj2), cl=cl, pop=list(pop1, pop2)))
}

# generates multi-view DCSBM with 10% s nodes for simulation
make_mv_dcsbm_mod <- function(n, K, r, s, delta, hub) { 
  new_model(name = "mv-dcsbm-ind",
            label = sprintf("mv-dcsbm-ind (n = %s, K = %s, r = %s, s  = %s, delta = %s, hub=%s)", 
                            n, K, r, s, delta, ifelse(hub==T, "T", "F")),
            params = list(n = n, K=K, r=r, s=s, delta=delta, hub=hub), 
            simulate = function(nsim, n, K, r, s, delta) { 
              lapply(1:nsim, function(x) 
                generate_mv_dcsbm(n, K, r, s, delta, hub))
            })
}
