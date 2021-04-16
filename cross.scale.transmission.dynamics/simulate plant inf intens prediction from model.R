# to be insertd into code similar to that in 'predict plant inf intens change funcs.R'
bootMer(n.pustules.model, FUN=function(x)predict(x, pred.data, exclude='s(site)'),nsim=1)$t
