library(simulator)
library(dplyr)
library(ggplot2)

p <- 10
n <- 500 
s <- 0.015 
K <- 3

all.settings <- rbind(as.matrix(cbind(c(1.5, 2), 
                                      c(8.4, 6.0), 0)), 
                      as.matrix(cbind(expand.grid(c(1.25, 1.5, 1.75, 2)), 
                                      expand.grid(seq(9.6, 6, length=4)), 1)))
colnames(all.settings) <- c("r", "sig", "hub")

results <- data.frame()

for(i in 1:6) { 
  this.setting <- all.settings[i, ]
  r <- as.numeric(this.setting[1])
  sig <- as.numeric(this.setting[2])
  hub <- as.numeric(this.setting[3])
  
  name_of_simulation <- paste("mv-node-ind-n", n, "-K", K, 
                              "-inout", r, "-sp", s, 
                              "-p", p, "-sig", sig, "-hub", hub, sep="")
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

options(ggplot2.discrete.colour=list(RColorBrewer::brewer.pal(10, "Paired")[c(3:4, 7:10)]))
theme_set(theme_bw(base_size=18) + theme(legend.position="bottom"))

p1 <- ggplot(aes(x=delta, y = power, ymin = power - qnorm(0.975)*se, 
                 ymax = power + qnorm(0.975)*se, 
                 colour=factor(Method), lty=factor(r), 
                 group=interaction(factor(r), factor(Method))), 
             data=results[results$hub == 0 & results$r %in% c(1.5, 2), ]) + 
  geom_hline(yintercept=0.05, linetype="solid", size=0.25) + 
  geom_linerange(size=0.75,  linetype="solid", show.legend=F) +  geom_line(size=0.75) + 
  xlab("Dependence between views") +
  ylab("Power at nominal 5% significance level") +
  scale_colour_discrete(name = "Methods", 
                        labels=c(expression("G-test, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression("G-test, estimated "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, estimated "~K^{(1)}~"and"~K^{(2)}), 
                                 expression("BESTest, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression("BESTest, estimated "~K^{(1)}~"and"~K^{(2)})), 
                        guide=guide_legend(nrow=3, byrow=T)) + 
  scale_linetype_manual(name = expression("(r, "~sigma~")"), 
                        values = c("solid", "22"), 
                        labels=c("(1.5, 8.4)", "(2, 6.0)"), 
                        guide=guide_legend(nrow=2, byrow=T)) + 
  annotate(geom="text", label=as.character(expression(alpha==~"0.05")),parse=T,  
           x=0.08, y=0.02, vjust=0.3) + 
  ggtitle(expression(X^{(1)}~"follows a SBM"))  

ggsave("Figure4.pdf", p1, "pdf", width=12, height=6)

options(ggplot2.discrete.colour=list(RColorBrewer::brewer.pal(4, "Set2")))

p2 <- ggplot(aes(x=delta, y = power, ymin = power - qnorm(0.975)*se, 
                 ymax = power + qnorm(0.975)*se, 
                 colour=factor(Method), lty=factor(r), 
                 group=interaction(factor(r), factor(Method))), 
             data=results[results$hub == 1 & results$r %in% c(1.5, 2) & 
                            results$Method %in% c("gtest", "p2lrt", "p2lrt-spectral", 
                                                  "gtest-spectral"), ]) + 
  geom_hline(yintercept=0.05, linetype="solid", size=0.25) + 
  geom_linerange(size=0.75,  linetype="solid", show.legend=F) +  geom_line(size=0.75) + 
  xlab("Dependence between views") +
  ylab("Power at nominal 5% significance level") +
  scale_colour_discrete(name = "Methods", 
                        labels=c(expression("G-test, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression("G-test, estimated "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, true "~K^{(1)}~"and"~K^{(2)}),
                                 expression(P^2~"LRT, estimated "~K^{(1)}~"and"~K^{(2)})), 
                        guide=guide_legend(nrow=2, byrow=T)) + 
  scale_linetype_manual(name = expression("(r, "~sigma~")"), 
                        values = c("solid", "22"), 
                        labels=c("(1.5, 8.4)", "(2, 6.0)"), 
                        guide=guide_legend(nrow=2, byrow=T)) + 
  annotate(geom="text", label=as.character(expression(alpha==~"0.05")),parse=T,  
           x=0.08, y=0.02, vjust=0.3) + 
  ggtitle(expression(X^{(1)}~"follows a DCSBM"))  


ggsave("FigureS2.pdf", p2, "pdf", width=12, height=6)
