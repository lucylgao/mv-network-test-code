library(ggplot2)
library(reshape2)
library(scales)

load("final-results-spectral.Rdata")
pi1 <- p2lrt.res$modelfit1$pi
pi2 <- p2lrt.res$modelfit2$pi
Pi <- p2lrt.res$Pi
C <- Pi/tcrossprod(pi1, pi2)

order1 <- order(pi1, decreasing=F)
order2 <- order(pi2, decreasing=F)
pi1.order <- pi1[order1]
pi2.order <- pi2[order2]
C.order <- C[order1, order2]

C.frame <- melt(C.order)
names(C.frame) <- c("v1", "v2", "value")
C.frame$log.value <- log2(C.frame$value)

Cplot <- ggplot(C.frame, aes(v2, rev(v1))) + geom_tile(aes(fill = value)) + 
  scale_x_discrete(breaks=1:14, labels=as.character(1:14),  limits=1:14, position="top") + 
  scale_y_discrete(breaks=1:14, labels=as.character(14:1),  limits=1:14) + 
  scale_fill_gradient2(low="steelblue", mid="white", high="red", 
                       midpoint=0, limits=range(C.frame$value), 
                       trans="log2", breaks=2^seq(-14, 4, by=2), 
                       labels=parse(text=paste("2^", seq(-14, 4, by=2))), 
                       name="") + 
  theme_grey(base_size = 20) +
  theme(legend.position="bottom", legend.key.width=unit(2,"cm")) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        legend.background=element_rect(fill="grey95")) + 
  ggtitle(expression("Heatmap of"~hat(C))) + 
  guides(fill = guide_colourbar(ticks.colour = "black",
                                ticks.linewidth = 1))

ggsave("Figure5-C.pdf", Cplot, width=5, height=6)

pi1.frame <- melt(pi1.order)
pi1.frame$v1 <- 1
pi1.frame$v2 <- rev(1:14)

pi1plot <- ggplot(pi1.frame, aes(v1, v2)) +  geom_tile(aes(fill = value)) + 
  scale_x_discrete(breaks=1, labels="", limits=1, name="") + 
  scale_y_discrete(breaks=1:14, labels=as.character(14:1),  limits=1:14, name="") + 
  scale_fill_gradient(low="steelblue", high="white", 
                      trans="log2", breaks=2^seq(-10, 0, by=1), 
                      labels=parse(text=paste("2^", seq(-10, 0, by=1))),
                      name="") + theme_grey(base_size = 20) +
  theme(legend.position="right", legend.key.height=unit(2,"cm")) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        legend.background=element_rect(fill="grey95")) + 
  ggtitle(expression("Heatmap of"~hat(pi)^"(1)")) + 
  guides(fill = guide_colourbar(ticks.colour = "black",
                                ticks.linewidth = 1))

ggsave("Figure5-pi1.pdf", pi1plot, width=3, height=6)

pi2.frame <- melt(pi2.order)
pi2.frame$v1 <- 1
pi2.frame$v2 <- rev(1:14)

pi2plot <- ggplot(pi2.frame, aes(v1, v2)) +  geom_tile(aes(fill = value)) + 
  scale_x_discrete(breaks=1, labels="", limits=1, name="") + 
  scale_y_discrete(breaks=1:14, labels=as.character(14:1),  limits=1:14, name="") + 
  scale_fill_gradient(low="steelblue", high="white", 
                      trans="log2", breaks=2^seq(-8, 0, by=1), 
                      labels=parse(text=paste("2^", seq(-8, 0, by=1))),
                      name="") + theme_grey(base_size = 20) +
  theme(legend.position="right", legend.key.height=unit(2,"cm")) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), 
        legend.background=element_rect(fill="grey95")) + 
  ggtitle(expression("Heatmap of"~hat(pi)^"(2)")) + 
  guides(fill = guide_colourbar(ticks.colour = "black",
                                ticks.linewidth = 1))

ggsave("Figure5-pi2.pdf", pi2plot, width=3, height=6)
