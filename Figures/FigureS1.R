library(ggplot2)
library(simulator)

offdiag <- 0.25
smallcl <- 1
sparsecl <- 0.5
n <- 50
K <- 2
rho <- 1
delta <- 0

batch <- 1
name_of_simulation <- paste("mv-dcsbm-n", n, "-K", K, 
                            "-popmod-1-offdiag", offdiag, "-sparsecl", sparsecl, 
                            "-smallcl", smallcl,  "-rho", rho, 
                            "-delta", delta, "-batch", batch, sep="")
sim <- load_simulation(name_of_simulation)
df <- sim %>% evals %>% as.data.frame

for(batch in 2:10) { 
  name_of_simulation <- paste("mv-dcsbm-n", n, "-K", K, 
                            "-popmod-1-offdiag", offdiag, "-sparsecl", sparsecl, 
                            "-smallcl", smallcl,  "-rho", rho, 
                            "-delta", delta, "-batch", batch, sep="")
  sim <- load_simulation(name_of_simulation)
  df.append <- sim %>% evals %>% as.data.frame
  df <- rbind(df, df.append)
}

maxalgK <- length(unique(df$Method)) 
power <- rep(NA, 49) 
nSim <- rep(NA, 49)

for(algK in 2:maxalgK) { 
  method <- paste("p2lrt-K", algK, sep="")
  pvals <- df[df$Method == method, ]$pval
  nSim[algK-1] <- length(pvals)
  power[algK-1] <- mean(pvals <= 0.05)
}

power.frame <- data.frame(K=2:50, power= power, nSim=200) 
power.frame$se <- sqrt(power.frame$power*(1 - power.frame$power)/power.frame$nSim)

title <- expression("Type I error rate of the"~P^2~"LRT")
g <- ggplot(power.frame, aes(x=K, y=power)) + 
  geom_point()+
  geom_errorbar(aes(ymin=pmax(0, power-1.96*se), ymax=pmin(pmax(power+1.96*se, 0), 1)), width=.2,
                position=position_dodge(0.05)) + geom_hline(yintercept=0.05, lty="dotted") +
  ylim(0, 1) + ylab("Power at nominal 5% significance level") + 
  ggtitle(title) + xlab("Number of communities used") + 
  theme_gray(base_size = 18) + coord_cartesian(ylim = c(0, 1)) + theme(legend.position="none") 
ggsave("FigureS1.pdf", g, "pdf", width=12, height=6)

pvalSpectral <- df[df$Method == "p2lrt-estK", ]$pval
type1Spectral <- mean(pvalSpectral <= 0.05)
varType1Spectral <- type1Spectral*(1 - type1Spectral)/200 
type1Spectral
round(type1Spectral - qnorm(0.975)*sqrt(varType1Spectral), 4)
round(type1Spectral + qnorm(0.975)*sqrt(varType1Spectral), 4)