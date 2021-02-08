library(ggplot2)
library(patchwork)
library(dplyr)
library(simulator)

n <- 1000
K <- 6

all.settings <- cbind(expand.grid(rev(c(1.5, 1.75, 2)), 
                                  rev(c(0.01, 0.02)), rev(c(0, 1))))
colnames(all.settings) <- c("r", "s", "hub")

results <- data.frame()

for(i in 1:12) { 
  this.setting <- all.settings[i, ]
  r <- as.numeric(this.setting[1])
  s <- as.numeric(this.setting[2])
  hub <- as.numeric(this.setting[3])
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
                 lty=factor(r), colour=factor(Method), 
                 group=interaction(factor(r), factor(Method))), 
             data=results[results$s == 0.01 & results$hub == 0 & 
                            results$r %in% c(1.5, 1.75, 2), ]) + 
  geom_linerange(size=0.75,  linetype="solid", show.legend=F) +  geom_line(size=0.75) + 
  geom_hline(yintercept=0.05, linetype="solid", size=0.25) + 
  annotate(geom="text", label=as.character(expression(alpha==~"0.05")),parse=T,  
           x=0.08, y=0.02, vjust=0.3) + 
  scale_colour_discrete(name="Method", guide=guide_legend(nrow=2, byrow=T), 
                        labels=c(expression("G-test, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression("G-test, estimated "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, estimated "~K^{(1)}~"and"~K^{(2)}))) + 
  scale_linetype_manual(name="r", values=c("22", "solid", "twodash"), 
                        guide=guide_legend(nrow=2, byrow=T, keywidth=3, keyheight=1)) + 
  ggtitle("s = 0.01") + xlab(expression("Dependence between views ("~Delta~")")) + 
  ylab("Power at nominal 5% significance level")

p2 <- ggplot(aes(x=delta, y = power, ymin = power - qnorm(0.975)*se, 
                 ymax = power + qnorm(0.975)*se, 
                 lty=factor(r), colour=factor(Method), 
                 group=interaction(factor(r), factor(Method))), 
             data=results[results$s == 0.02 & results$hub == 0 & 
                            results$r %in% c(1.5, 1.75, 2), ]) + 
  geom_linerange(size=0.75,  linetype="solid", show.legend=F) +  geom_line(size=0.75) + 
  geom_hline(yintercept=0.05, linetype="solid", size=0.25) + 
  annotate(geom="text", label=as.character(expression(alpha==~"0.05")),parse=T,  
           x=0.08, y=0.02, vjust=0.3) + 
  scale_colour_discrete(name="Method", guide=guide_legend(nrow=2, byrow=T), 
                        labels=c(expression("G-test, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression("G-test, estimated "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, estimated "~K^{(1)}~"and"~K^{(2)}))) + 
  scale_linetype_manual(name="r", values=c( "22", "solid", "twodash"), 
                        guide=guide_legend(nrow=2, byrow=T, keywidth=3, keyheight=1)) + 
  ggtitle("s = 0.02") + xlab(expression("Dependence between views ("~Delta~")")) + 
  ylab("Power at nominal 5% significance level")

ggsave(p1 + p2 + plot_layout(guides="collect"), file="Figure2.pdf", height=6, width=12)

g1 <- ggplot(aes(x=delta, y = power, ymin = power - qnorm(0.975)*se, 
                 ymax = power + qnorm(0.975)*se, 
                 lty=factor(r), colour=factor(Method), 
                 group=interaction(factor(r), factor(Method))), 
             data=results[results$s == 0.01 & results$hub == 1 & 
                            results$r %in% c(1.5, 1.75, 2), ]) + 
  geom_linerange(size=0.75,  linetype="solid", show.legend=F) +  geom_line(size=0.75) + 
  geom_hline(yintercept=0.05, linetype="solid", size=0.25) + 
  annotate(geom="text", label=as.character(expression(alpha==~"0.05")),parse=T,  
           x=0.08, y=0.02, vjust=0.3) + 
  scale_colour_discrete(name="Method", guide=guide_legend(nrow=2, byrow=T), 
                        labels=c(expression("G-test, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression("G-test, estimated "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, estimated "~K^{(1)}~"and"~K^{(2)}))) + 
  scale_linetype_manual(name="r", values=c("22", "solid", "twodash"), 
                        guide=guide_legend(nrow=2, byrow=T, keywidth=3, keyheight=1)) + 
  ggtitle("s = 0.01") + xlab(expression("Dependence between views ("~Delta~")")) + 
  ylab("Power at nominal 5% significance level")

g2 <- ggplot(aes(x=delta, y = power, ymin = power - qnorm(0.975)*se, 
                 ymax = power + qnorm(0.975)*se, 
                 lty=factor(r), colour=factor(Method), 
                 group=interaction(factor(r), factor(Method))), 
             data=results[results$s == 0.02 & results$hub == 1 & 
                            results$r %in% c(1.5, 1.75, 2), ]) + 
  geom_linerange(size=0.75,  linetype="solid", show.legend=F) +  geom_line(size=0.75) + 
  geom_hline(yintercept=0.05, linetype="solid", size=0.25) + 
  annotate(geom="text", label=as.character(expression(alpha==~"0.05")),parse=T,  
           x=0.08, y=0.02, vjust=0.3) + 
  scale_colour_discrete(name="Method", guide=guide_legend(nrow=2, byrow=T), 
                        labels=c(expression("G-test, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression("G-test, estimated "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, estimated "~K^{(1)}~"and"~K^{(2)}))) + 
  scale_linetype_manual(name="r", values=c( "22", "solid", "twodash"), 
                        guide=guide_legend(nrow=2, byrow=T, keywidth=3, keyheight=1)) + 
  ggtitle("s = 0.02") + xlab(expression("Dependence between views ("~Delta~")")) + 
  ylab("Power at nominal 5% significance level")

ggsave(g1 + g2 + plot_layout(guides="collect"), file="Figure3.pdf", height=6, width=12)
