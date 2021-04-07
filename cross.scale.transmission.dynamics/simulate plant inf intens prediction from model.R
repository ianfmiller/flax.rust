# to be insertd into code similar to that in 'predict plant inf intens change funcs.R'
bootMer(plant.model, FUN=function(x)predict(x, pred.data, re.form=NA),nsim=1)$t