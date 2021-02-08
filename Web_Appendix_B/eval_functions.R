## @knitr metrics

pval <- new_metric(name = "pval", label = "pval",
                   metric = function(model, out) out$pval)
stat <- new_metric(name = "stat", label = "stat",
                   metric = function(model, out) out$stat)
rejects <- new_metric(name = "rejects", label = "rejects",
                   metric = function(model, out) out$pval <= 0.05)

