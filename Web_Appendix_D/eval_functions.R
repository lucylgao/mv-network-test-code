## @knitr metrics

stat <- new_metric(name = "stat", label = "Test Statistic",
                   metric = function(model, out) out$stat)

pval <- new_metric(name = "pval", label = "P value",
                   metric = function(model, out) out$pval)

rejects <- new_metric(name="rejects", label="Test rejects at 0.05", 
                      metric = function(model, out)
                        out$pval <= 0.05)

