## @knitr metrics

stat <- new_metric(name = "stat", label = "Test Statistic",
                   metric = function(model, out) out$stat)

pval <- new_metric(name = "pval", label = "P value",
                   metric = function(model, out) out$pval)

rejects <- new_metric(name="rejects", label="Test rejects at 0.05", 
                      metric = function(model, out)
                        out$pval <= 0.05)


getK1 <- new_metric(name="getK1", label="getK1", 
                      metric = function(model, out)
                        out$K1)

getK2 <- new_metric(name="getK2", label="getK2", 
                     metric = function(model, out)
                       out$K2)

