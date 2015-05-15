library(RColorBrewer)
RFun<-function(n) {
  return(mean(rexp(n))) #maybe rlnorm, rnorm, runif, rexp, or even a combination of these
}

nrep<-10000
nsim.per.rep.vector<-c(1, 5, 10, 100, 1000)
results<-list()
densities<-list()
mypalette<-brewer.pal(length(nsim.per.rep.vector),"Dark2")
for (nsim.index in sequence(length(nsim.per.rep.vector))) {
  per.rep.results<-replicate(n=nrep, RFun(nsim.per.rep.vector[nsim.index]))
  results[[length(results)+1]]<-per.rep.results
  hist.output<-hist(per.rep.results, plot=FALSE)
  hist.density<-list(x= hist.output$mids, y=hist.output$density)
  #  densities[[length(densities)+1]]<-density(per.rep.results, kernel="rectangular", cut=0)
  densities[[length(densities)+1]]<-hist.density
  
}
x.min<-NA
x.max<-NA
y.min<-NA
y.max<-NA
for (i in sequence(length(densities))) {
  x.min<-min(c(densities[[i]]$x,x.min),na.rm=TRUE)
  x.max<-quantile(c(densities[[i]]$x,x.max), probs=0.999,na.rm=TRUE)
  y.min<-min(c(densities[[i]]$y,y.min),na.rm=TRUE)
  y.max<-max(c(densities[[i]]$y,y.max),na.rm=TRUE)
}
par(mfcol=c(1,3))
hist(results[[1]], col=mypalette[1], xlab="value", main="")
plot(x=c(x.min, 1.1*x.max), y=c(y.min, 1.1*y.max), type="n", bty="n", ylab="density", xlab="value")
for (i in sequence(length(densities))) {
  lines(densities[[i]], col=mypalette[i])
  text(x=densities[[i]]$x[which.max(densities[[i]]$y)], y=max(densities[[i]]$y), labels=nsim.per.rep.vector[i], pos=3, col=mypalette[i])
}


last.index<-length(densities)
plot(x=c(min(densities[[last.index]]$x), 1.1*max(densities[[last.index]]$x)), y=c(y.min, 1.1*y.max), type="n", bty="n", ylab="density", xlab="value")
lines(densities[[last.index]], col=mypalette[last.index])
text(x=densities[[last.index]]$x[which.max(densities[[last.index]]$y)], y=max(densities[[last.index]]$y), labels=nsim.per.rep.vector[last.index], pos=2, col=mypalette[last.index])

sample.mean <- mean(results[[last.index]])
sample.sd <- sd(results[[last.index]])
y.val <- dnorm(densities[[last.index]]$x, mean=sample.mean, sd=sample.sd)
y.val <- max(densities[[last.index]]$y)*y.val/max(y.val)
lines(x=densities[[i]]$x, y=y.val, lty="dotted", col="black")
text(x=densities[[last.index]]$x[which.max(densities[[last.index]]$y)], y=max(densities[[last.index]]$y), labels="normal", pos=4)

