library(simulator)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

n <- 1000
K <- 6
hub <- 0


all.settings <- rbind(c(2, 0.025), c(2.25, 0.025), c(2.5, 0.025))
colnames(all.settings) <- c("r", "s")

results <- data.frame()

for(i in 1:3) { 
  this.setting <- all.settings[i, ]
  r <- as.numeric(this.setting[1])
  s <- as.numeric(this.setting[2])
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

theme_set(theme_bw(base_size=18) + theme(legend.position="bottom"))

p <- ggplot(aes(x=delta, y = power, ymin = power - qnorm(0.975)*se, 
                ymax = power + qnorm(0.975)*se, 
                lty=factor(r), colour=factor(Method), 
                group=interaction(factor(r), factor(Method))), 
            data=results) + 
  geom_linerange(size=0.75,  linetype="solid", show.legend=F) +  geom_line(size=0.75) + 
  geom_hline(yintercept=0.05, linetype="solid", size=0.25) + 
  annotate(geom="text", label=as.character(expression(alpha==~"0.05")),parse=T,  
           x=0.08, y=0.02, vjust=0.3) + 
  scale_colour_manual(name="Method", guide=guide_legend(nrow=2, byrow=T), 
                        labels=c(expression("G-test, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression("G-test, estimated "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, estimated "~K^{(1)}~"and"~K^{(2)})), 
                      values=c("#f16913", "#807dba","#084594", "#238b45")) + 
  scale_linetype_manual(name="r", values=c("22", "solid", "twodash"), 
                        guide=guide_legend(nrow=2, byrow=T, keywidth=3, keyheight=1)) + 
  xlab(expression("Dependence between views ("~Delta~")")) + 
  ylab("Power at nominal 5% significance level")

ggsave(p, file="FigureS3.pdf", height=8, width=12)