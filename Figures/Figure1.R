library(igraph)

generate_memberships <- function(n, Pi) {
  stopifnot(Pi >= 0, abs(sum(Pi) - 1) < 1e-15)
  ii <- sample(length(Pi), n, replace = TRUE, prob = Pi)
  K1 <- nrow(Pi)
  K2 <- ncol(Pi)
  cpairs <- cbind(rep(1:K1, K2), rep(1:K2, each = K1))
  cpairs[ii, ]
}

# makes variable Pi for power curve according to Gao (2018)
make_Pi <- function(delta, K) { 
  return(delta * diag(K) / K + (1 - delta) * matrix(1, K, K) / K^2 )  
}

# makes B according to Amini (2013)
make_B <- function(beta, lambda, pi.vect, pop) { 
  K <- length(pi.vect)
  n <- length(pop)
  pop.mean <- mean(pop)
  
  # 1/beta times more likely to stay within community 
  B0 <- diag(1/beta - 1, K) + matrix(1, K, K)
  
  # Rescale to get overall expected network degree lambda
  B <- lambda*B0/((n-1)*c(t(pi.vect)%*%B0%*%pi.vect)*pop.mean^2)
  return(B)  
}

make_adjacency_matrix_dcsbm <- function(B, cl, pop) { 
  n <- length(cl)
  K <- nrow(B)  
  
  # n x n matrix with ijth entry = theta_i theta_j B_{Zi, Zj}
  Zmat <- matrix(0, n, K)
  Zmat[cbind(1:n, cl)] <- 1
  outer <- pop*Zmat
  probmat <- outer%*%B%*%t(outer)
  diag(probmat) <- 0
  
  # Generate Bernoulli upper triangle entries of adjacency matrix 
  upper.index <- which(upper.tri(probmat))
  upper.probs <- probmat[upper.index]
  upper.unif <- runif(n = length(upper.probs))
  upper.adj <- rep(0, length(upper.probs))
  upper.adj[upper.unif < upper.probs] <- 1
  
  # Generate adjacency matrix
  adj <- matrix(0, n, n)
  adj[upper.index] <- upper.adj
  adj <- adj + t(adj)
  diag(adj) <- 0
  
  return(adj)
}

set.seed(1)
n <- 10
K <- 2
lambda <- 2
beta <- 1/10
pi.vect <- rep(1/K, K)
pop1 <- rep(1, n)
pop2 <- rep(1, n)
B1 <- make_B(beta, lambda, pi.vect, pop1)
B2 <- make_B(beta, lambda, pi.vect, pop2)
Pi <- make_Pi(0.5, K)
cl <- generate_memberships(n, Pi)
adj1 <- make_adjacency_matrix_dcsbm(B1, cl[, 1], pop1)
adj2 <- make_adjacency_matrix_dcsbm(B2, cl[, 2], pop2)
net1 <- graph_from_adjacency_matrix(adj1, mode="undirected")
l1 <- layout_in_circle(net1)

pdf("Figure-1i.pdf", height=4, width=8)
par(mfrow=c(1,2)) 
par(mar=c(0,0,1.5,0))
plot.igraph(net1, vertex.color = "#737373", 
                  vertex.shape = "circle",
                  vertex.label.color="white", layout=l1,  asp=0, 
                  vertex.size=30, vertex.label.cex=2.5)
title("View 1", cex.main=2)
net2 <- graph_from_adjacency_matrix(adj2, mode="undirected")
l2 <- layout_in_circle(net2)
p2 <- plot.igraph(net2, vertex.color = "#C4C4C4", 
                  vertex.shape = "circle",
                  vertex.label.color="white", 
                  label.cex=3, layout=l2,  asp=0, 
                  vertex.size=30, vertex.label.cex=2.5)
title("View 2", cex.main=2)
dev.off()

pdf("Figure-1ii.pdf", height=4, width=8)
par(mfrow=c(1,2)) 
par(mar=c(0,0,1.5,0))
plot.igraph(net1, vertex.color = "#737373", 
                  vertex.shape = "circle",
                  vertex.label.color="white", layout=l1,  asp=0, 
                  vertex.size=30, vertex.label.cex=2.5)
title("View 1", cex.main=2)
plot(c(100, 200), c(300, 450), type= "n", xlab = "", ylab = "", axes=F)
rect(110, 325, 200, 425, density = NA, border = "#C4C4C4", col="#C4C4C4")
text(104, 418, c("1"), cex=2)
text(104, 403, c("2"), cex=2)
text(104, 333, c("10"), cex=2)
text(115, 433, c("1"), cex=2)
text(123, 433, c("2"), cex=2)
text(197, 433, c("p"), cex=2)
title("View 2", cex.main=2)
dev.off()