library(simulator) 
library(Matrix)
library(multiviewtest)

n <- 50
K <- 2
# settings for varying pi and theta
offdiag <- 0.25 
smallcl <- 1
sparsecl <- 0.5
rho <- 1 # correlation between popularities
delta <- 0

all.settings <- 1:10

# For running as a batch job
if(Sys.getenv("SGE_TASK_ID") != "") { 
  this.setting <- as.numeric(Sys.getenv("SGE_TASK_ID"))
} else { 
  this.setting <- 1
}

batch <- all.settings[this.setting]

source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")

runMethods <- list(p2lrt_estK) 
  
for(i in 1:length(c(2:50))) { 
  runMethods[[i+1]] <- p2lrt_vect[[(2:50)[i] - 1]]
} 

name_of_simulation <- paste("mv-dcsbm-n", n, "-K", K, 
                            "-popmod-1-offdiag", offdiag, "-sparsecl", sparsecl, 
                            "-smallcl", smallcl,  "-rho", rho, 
                            "-delta", delta, "-batch", batch, sep="")

sim <- load_simulation(name = name_of_simulation)

sim <- sim %>% run_method(runMethods, parallel = list(socket_names = 10, 
                             libraries = c("randnet", "matrixStats", "multiviewtest"))) 

sim <- sim %>% evaluate(list(pval, stat, rejects))

df <- as.data.frame(evals(sim))