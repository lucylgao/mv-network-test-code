## @knitr init
library(methods) 
library(simulator) 
library(multiviewtest)
library(Matrix)

n <- 1000
K <- 6
hub <- 0


# r is the strength of the communities 
# s is the expected edge density of the networks
all.settings <- rbind(c(2, 0.025), c(2.25, 0.025), c(2.5, 0.025))
colnames(all.settings) <- c("r", "s")

# For submitting as a batch job on the cluster
if(Sys.getenv("SGE_TASK_ID") != "") { 
  this.sim.id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
} else { 
  this.sim.id <- 1
}

this.setting <- all.settings[this.sim.id, ]
r <- as.numeric(this.setting[1])
s <- as.numeric(this.setting[2])

source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")

name_of_simulation <- paste("mv-dcsbm-ind-n", n, "-K", K, 
                            "-inout", r, "-sp", s, "-hub", hub, sep="")

sim <- new_simulation(name = name_of_simulation,
                      label = paste("mv-dcsbm-ind-n", n, "-K", K, 
                                    "-inout", r, "-sp", s, "-hub", hub, sep=""))

sim <- sim %>% generate_model(make_mv_dcsbm_mod,
                              n=n, K=K, r=r, s=s, hub=hub,
                              delta = as.list(seq(0, 1, length=9)),
                              vary_along = "delta")
sim <- sim  %>% simulate_from_model(nsim = 167, index=1:12) 

## @knitr main
sim <- sim  %>%
  run_method(list(g_test, pplr_test, g_test_spectral, pplr_test_spectral),
             parallel=list(socket_names = 12, 
                           libraries=c("randnet", "multiviewtest", "matrixStats"))) %>% 
  evaluate(list(pval, rejects, stat))

# compile results
ev <- sim %>% evals %>% as.data.frame