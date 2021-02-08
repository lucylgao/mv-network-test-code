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

make_Pi <- function(delta, smallcl=1, K) { 
  probs <- rep(0, K) 
  probs[1] <- smallcl/K 
  probs[2:K] <- (1 - smallcl/K)/(K-1)
  return(delta * diag(probs) + (1 - delta) * tcrossprod(probs, probs) )  
}

gen_gauss_cop <- function(r, n){
  rho <- 2 * sin(r * pi/6)        # Pearson correlation
  P <- toeplitz(c(1, rho))        # Correlation matrix
  d <- nrow(P)                    # Dimension
  ## Generate sample
  U <- pnorm(matrix(rnorm(n*d), ncol = d) %*% chol(P))
  return(U)
}


# offdiag = off-diagonal elements 
# sparsecl = B[1, 1]
# K = number of subgroups 
make_B <- function(offdiag, sparsecl = 1, K) { 
  B <- diag(c(sparsecl, rep(1, (K-1))) - offdiag) + matrix(offdiag, K, K)
  
  return(B)
}

generate_mv_dcsbm <- function(n, K, offdiag, sparsecl, delta, smallcl, rho) { 
  B <- make_B(offdiag, sparsecl, K)
  Pi <- make_Pi(delta, smallcl, K)
  cl <- generate_memberships(n, Pi)
  pops <- gen_gauss_cop(rho, n)
  pops <- 0.7*(pops + 0.2)
  pmat1 <- make_dcsbm_prob_mat(B, pops[, 1], cl[, 1])
  pmat2 <- make_dcsbm_prob_mat(B, pops[, 2], cl[, 2])
  adj1 <- generate_random_graph(pmat1)
  adj2 <- generate_random_graph(pmat2)
  
  while(any(rowSums(adj1) == 0) | any(rowSums(adj2) == 0)) { 
    adj1 <- generate_random_graph(pmat1)
    adj2 <- generate_random_graph(pmat2)
  }
  
  return(list(adj=list(adj1, adj2), cl=cl, pop=pops))
}

make_mv_dcsbm_mod <- function(n, K, offdiag, sparsecl, delta, smallcl, rho, batch) { 
  
  new_model(name = "mv-dcsbm",
            label = sprintf("n = %s, K = %s, offdiag=%s, sparsecl=%s, delta=%s, smallcl=%s, rho=%s, batch=%s", 
                            n, K, offdiag, sparsecl, delta, smallcl, rho, batch), 
            params = list(n=n, K=K, offdiag=offdiag, sparsecl=sparsecl, delta=delta, smallcl=smallcl, rho=rho, batch=batch), 
            simulate = function(nsim, n, K, offdiag, sparsecl, delta, smallcl, rho) { 
             lapply(1:nsim, function(x) 
                generate_mv_dcsbm(n, K, offdiag, sparsecl, delta, smallcl, rho))
            })
}

