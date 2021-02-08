library(simulator) 
library(Matrix)
library(multiviewtest)

n <- 50
K <- 2
delta <- 0 # independent communities
# settings for varying pi and theta
offdiag <- 0.25 
smallcl <- 1
sparsecl <- 0.5
rho <- 1 # correlation between popularities

source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")

all.settings <- 1:10

# For running as a batch job
if(Sys.getenv("SGE_TASK_ID") != "") { 
  this.setting <- as.numeric(Sys.getenv("SGE_TASK_ID"))
} else { 
  this.setting <- 1
}

batch <- all.settings[this.setting]

name_of_simulation <- paste("mv-dcsbm-n", n, "-K", K, 
                            "-popmod-1-offdiag", offdiag, "-sparsecl", sparsecl, 
                            "-smallcl", smallcl,  "-rho", rho, 
                            "-delta", delta, "-batch", batch, sep="")

sim <- new_simulation(name = name_of_simulation,
                      label = name_of_simulation)

sim <- sim %>% generate_model(make_mv_dcsbm_mod,
                              n=n, K=K, offdiag=offdiag, 
                              sparsecl=sparsecl, smallcl= smallcl, rho = rho,  
                              delta = delta, batch=batch) 

sim <- sim %>% simulate_from_model(nsim=2, index=((batch - 1)*10 +1):(batch*10))