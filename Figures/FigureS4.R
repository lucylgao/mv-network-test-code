library(ggplot2)
library(patchwork)
library(simulator)
library(dplyr)

n <- 250
all.settings <- cbind(expand.grid(rev(c(1.5)), 
                                  rev(c(0.05)), c(0), c(2, 3, 4)))
colnames(all.settings) <- c("r", "s", "hub", "K")

results <- data.frame()

for(i in 1:3) { 
  this.setting <- all.settings[i, ]
  r <- as.numeric(this.setting[1])
  s <- as.numeric(this.setting[2])
  hub <- as.numeric(this.setting[3])
  K <- as.numeric(this.setting[4])
  name_of_simulation <- paste("mv-dcsbm-ind-n", n, "-K", K, 
                              "-inout", r, "-sp", s, "-hub", hub, sep="")
  sim <- load_simulation(name_of_simulation)
  ev <- sim %>% evals %>% as.data.frame 
  m <- sim %>% model %>% as.data.frame
  m$Model <- factor(m$name, levels(ev$Model)) # preparing for "join" of these two data frames
  df <- ev %>% group_by(Model, Method) %>% 
    summarise(power = mean(rejects, na.rm=T), 
              se = sd(rejects, na.rm=T)/sqrt(n() - sum(is.na(rejects))), 
              num.na = sum(is.na(rejects))) %>%
    left_join(m, by = "Model") # matches up the models and corresponding evals
  
  results <- rbind(results, as.data.frame(df))
}

options(ggplot2.discrete.colour=list(RColorBrewer::brewer.pal(4, "Set2")))
theme_set(theme_bw(base_size=18) + theme(legend.position="bottom"))

p1 <- ggplot(aes(x=delta, y = power, ymin = power - qnorm(0.975)*se, 
                 ymax = power + qnorm(0.975)*se, 
                 lty=factor(K), colour=factor(Method), 
                 group=interaction(factor(K), factor(Method))), 
             data=results[results$Method %in% c("gtest", "p2lrt"), ]) +  geom_linerange(size=0.75,  linetype="solid", show.legend=F) +  geom_line(size=0.75) + 
  geom_hline(yintercept=0.05, linetype="solid", size=0.25) + 
  scale_colour_discrete(name="Method", labels=c("G-test", expression(P^2~"LRT"))) + 
  scale_linetype_manual(name="K", values=c("twodash", "solid", "22"), 
                        guide=guide_legend(keywidth=3)) + 
  annotate(geom="text", label=as.character(expression(alpha==~"0.05")),parse=T,  
           x=0.08, y=0.02, vjust=0.3) + 
  xlab(expression("Dependence between views ("~Delta~")")) + 
  ylab("Power at nominal 5% significance level")

ggsave(p1, file="FigureS4.pdf", height=6, width=12)
