## @knitr init
library(methods) 
library(simulator) 
library(Matrix)
library(multiviewtest)

n <- 500 
K <- 3 
p <- 10
s <- 0.015 # s is the expected edge density of the network

# All possible settings for Section 7.3 and Web Appendix C
# for hub = 0, X1 follows SBM 
# for hub = 1, X1 follows DCSBM 
# r is the strength of the communities 
# sig is the variance of the clusters in the multivariate view
# (r, sig, hub)
all.settings <- rbind(c(1.5, 8.4, 0), c(2.0, 6.0, 0), c(1.5, 8.4, 1), c(2.0, 6.0, 1))
colnames(all.settings) <- c("r",  "sig", "hub")

# For submitting as a batch job on the cluster
if(Sys.getenv("SGE_TASK_ID") != "") { 
  this.sim.id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
} else { 
  this.sim.id <- 1
}

this.setting <- all.settings[this.sim.id, ]
r <- as.numeric(this.setting[1])
sig <- as.numeric(this.setting[2])
hub <- as.numeric(this.setting[3])

source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")

name_of_simulation <- paste("mv-node-ind-n", n, "-K", K, 
                            "-inout", r, "-sp", s, 
                            "-p", p, "-sig", sig, "-hub", hub, sep="")

sim <- new_simulation(name = name_of_simulation, label = name_of_simulation)

sim <- sim %>% generate_model(make_mv_node_mod,
                              n = n, K=K, r=r, s=s, p=p, 
                              sig=sig, hub=hub,
                              delta = as.list(seq(0, 1, length=9)),
                              vary_along = "delta")

sim <- sim  %>%  simulate_from_model(nsim = 167, index=1:12) 

## @knitr main
sim <- sim  %>%
  run_method(list(g_test, pplr_test, 
                  g_test_spectral, pplr_test_spectral, 
                  peel_test, peel_test_spectral),
             parallel=list(socket_names = 12, 
                           libraries=c("randnet", "matrixStats", "mclust", "multiviewtest"))) %>% 
  evaluate(list(pval, rejects, stat))

ev <- sim %>% evals %>% as.data.frame
