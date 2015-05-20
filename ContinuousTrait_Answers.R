#Learning outcomes:

#Skills
#	Ability to simulate a tree
#	Ability to simulate data
#	How to plot trees and data
#	Running independent contrasts
#	Estimating rates or other parameters
#	Pulling in data and trees

#Understanding:
#	Why to check data and results
#	Why to look at all program options
#	How to use different measures of fit
#	Uncertainty in estimates



library('ape')
library('phytools')
library('TreeSim')
library('geiger')
library('OUwie')
library('treebase')
library('taxize')
library('RColorBrewer')
install.packages('rfishbase') #b/c this isn't in the task view
library('rfishbase')

#first, let's simulate a tree
#because TreeSim returns a list of trees, even with just one tree simulated, take first element in list.
phy <- sim.bd.taxa(n=20,numbsim=1,lambda=1,mu=0.9)[[1]] 

#let's look at this:
plot(phy)

#huh, that's weird. Why are the tips at different times?
?sim.bd.taxa

#ah, we forgot to repress the extinct lineages. We could do this in TreeSim or in geiger (drop.extinct). Let's do it in TreeSim
phy <- sim.bd.taxa(n=30,numbsim=1,lambda=1,mu=0.9, complete=FALSE)[[1]] 


#look at it again
plot(phy)

#ah, that looks right (except for those studying taxa at different time points, like viruses or fossils).

#let's simulate some correlated data
vcv.data <- matrix(data=c(.5,0.3, 0.3, 0.45), nrow=2, ncol=2, byrow=TRUE)
print(vcv.data)
data <- sim.corrs(phy, vcv.data, anc=c(17,3))

#now let's look at the results
plot(x=data[,1], y=data[,2], pch=16, bty="n")

#that's the raw data; how does it look on a tree?
par(mfcol=c(1,2))
contMap(phy, data[,1])
contMap(phy, data[,2], direction="leftwards")

#or let's look at it in space:
par(mfcol=c(1,1))
phylomorphospace(phy, data)


#Use independent contrasts

pic.X <- pic(data[,1], phy)
pic.Y <- pic(data[,2], phy)

par(mfcol=c(1,2))
plot(data[,1], data[,2], pch=16, bty="n", col="gray")
plot(pic.X, pic.Y, pch=16, bty="n", col="gray")

#for contrasts, you should positivize them, since the order doesn't matter. This is NOT taking absolute value.
PositivizeContrasts <- function(x, y) {
	x.positivized <- x * sign(x)
	y.positivized <- y * sign(x)
	return(cbind(x.positivized, y.positivized))
}

positivized.results <- PositivizeContrasts(pic.X, pic.Y)
points(positivized.results[,1], positivized.results[,2], pch="x", col="red")

print(cor.test(positivized.results[,1], positivized.results[,2]))

#There are many other multivariate models in R for phylogenetic data, but let's go to simpler models for univariate traits

help(package="geiger") #to get help on all the functions

#note all the deprecated functions. These are functions that currently work, but are slated for eventual deletion. You'll want to use the non-deprecated ones. 

data.x <- data[,1] #just to save on typing

#ok, let's try brownian motion
geiger.brownian <-  fitContinuous(phy, data.x, SE=NA, ncores=1)
print(geiger.brownian)
#what does all this stuff mean?
#convergence, log likelihood, etc.

#did we do well estimating the rate?
print(paste("true rate was", vcv.data[1,1], "; estimated rate was", geiger.brownian$opt$sigsq))

#looks good (probably). But let's dig into it some more:
?fitContinuous

#there are different kinds of models, but there are other options, too: bounds, niter, etc. What are those for?

#Let's try running this multiple times, changing the number of iterations to just 2, and try that a few times

num.reps <- 25
optimal.sigsq.vector <- rep(NA, num.reps)

for (rep.index in sequence(num.reps)) {
	optimal.sigsq.vector[rep.index] <- fitContinuous(phy, data.x, SE=NA, ncores=1, control=list(niter=2))$opt$sigsq
}

plot(optimal.sigsq.vector, log="y", xlab="replicate", ylab="sigma squared", pch=16, bty="n")
abline(h=vcv.data[1,1], col="red")

#why are some of the estimates off the red line showing the truth?

#geiger now does many starting points and different optimizers (earlier versions weren't as thorough). You can see this in your first geiger run:

print(geiger.brownian$res)

#this shows for each replicate which optimizer was used, what it found as the best estimate for sigma squared and SE, the log likelihood, and the convergence diagnostic. A value of 0 for convergence indicates it worked fine, according to R, but note that many of those that converged actually didn't get at the right answer. This isn't a flaw in geiger but rather something about how optimization functions, even for simple problems like this. 



#GETTING DATA IN

fish.trees <- search_treebase("menadoensis", by="taxon", branch_lengths=TRUE, max_trees=30) #only one fish tree in treebase matches this currently

fish.phy <- fish.trees[[1]]

fish.chronogram <- chronos(fish.phy) #chronos isn't an ideal dating function, but it may be the best in R right now. r8s, treepl, or Beast may be better
class(fish.chronogram) <- "phylo" #we don't need this to also be class chronos

fish.sizes <- getSize(value="length")

fish.chronogram.name.resolve <- gnr_resolve(names = gsub("_", " ",fish.chronogram$tip.label), best_match_only=TRUE)

print(fish.chronogram.name.resolve)

fish.chronogram.name.resolve.ordered <- fish.chronogram.name.resolve$results[order(as.numeric(row.names(fish.chronogram.name.resolve$results))),] #puts it back in input order

fish.chronogram.renamed <- fish.chronogram
fish.chronogram.renamed$tip.label <- fish.chronogram.name.resolve.ordered$matched_name

#now, check again!
print(cbind(fish.chronogram$tip.label, fish.chronogram.renamed$tip.label))

#we would ideally resolve the size data to the same names. It will take too long. So we'll instead just change
#underscores to spaces. We'll miss some matches.
fish.sizes.renamed <- fish.sizes
names(fish.sizes.renamed) <- gsub("_", " ", names(fish.sizes))

#and check!
print(cbind(names(fish.sizes.renamed), names(fish.sizes)))

fish.sizes.renamed <- fish.sizes.renamed[!is.na(fish.sizes.renamed)] #get rid of NAs

fish.resolved.all <- treedata(fish.chronogram.renamed, fish.sizes.renamed, sort=TRUE)
fish.phy.final <- fish.resolved.all$phy
fish.data.final <- fish.resolved.all$data

#now you should look at these


#DOING ANALYSES

data(tworegime) #a sample dataset included with OUwie. It's made up data but runs well for examples. It is not identifiable, though (Ho and Ane, 2014) and will be pulled from OUwie soon

ls() #what does that dataset have in it?

print(tree) #prints a summary of the tree object

#Notice it has node labels. This is how regimes are mapped onto the tree. 

print(trait) #prints what the trait object is.

#But is it a matrix, a list, a data.frame, what? Two ways to look at this:
class(trait) #just tells what the class is
dput(trait) #prints the internal structure of the object. Can be handy with debugging








#Plot the tree and the internal nodes to highlight the selective regimes:
select.reg<-character(length(tree$node.label)) 
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"
plot(tree) 
nodelabels(pch=21, bg=select.reg)








#Now to run it!
#This takes longer than you may be used to. 
#We're a bit obsessive about doing multiple starts and in general
#performing a thorough numerical search. It took you 3+ years
#to get the data, may as well take an extra five minutes to 
#get an accurate answer
nodeBased.OUMV <- OUwie(tree,trait,model="OUMV", simmap.tree=FALSE, diagn=FALSE)
print(nodeBased.OUMV)
#What do the numbers mean?









#Let's try other models, too
#The models are "BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA" (and you've already done "OUMV")

#Start with BM1 and OU1 (YOU figure out how to do this). What do these models mean? What are their results?













#Now let's go over all the models. What does BMS mean?

#This is an easier way to write what you were doing. A faster way to analyze it would be to use mclapply rather than lapply (see the parallels package)
OUwie.model<-function(model, phy, data) {
	print(paste("Now starting model",model))
	return(OUwie(phy, data, model, simmap.tree=FALSE, diagn=FALSE))	
}
models <- c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA")
results <- lapply(models, OUwie.model, phy=tree, data=trait)










#Get the AICc values, put in a vector
AICc.values<-sapply(results, "[[", "AICc")
names(AICc.values)<-models
AICc.values<-AICc.values-min(AICc.values)


print(AICc.values) #The best model is the one with smallest AICc score

best<-results[[which.min(AICc.values)]] #store for later

print(best) #prints info on best model

#Ho and Ane (2014) point out that the Zhang & Siegmund (2007) modified BIC works better than AIC for determining the best number of regimes.
#This will take some work to implement (but see the phylolm package for their model, but which currently allows shifts in theta regime only). 
#However, it may be of interest how robust the best model selection is to different penalties for numbers of parameters.

#Get the AIC values, put in a vector
AIC.values<-sapply(results, "[[", "AIC")
names(AIC.values)<-models

#Get the loglik values, put in a vector
loglik.values<-sapply(results, "[[", "loglik")
names(loglik.values)<-models

#AIC = -2lnL + 2K
#AIC + 2 lnL = 2K
#K = (AIC + 2 lnl)/2

K.values = (AIC.values + 2 * loglik.values)/2

#Check
print(K.values)

my.palette<-brewer.pal(length(K.values),"Dark2")

multiplier.range <- seq(from=0, to=10, by=1)

plot(x=c(-0.5+min(multiplier.range), 1.05*max(multiplier.range)), y=range(-2*max(loglik.values), -2*min(loglik.values)+max(multiplier.range)*max(K.values)), xlab="Multiplier", ylab="*IC", bty="n", type="n")

for (model.index in sequence(length(K.values))) {
	lines(x=multiplier.range, y=-2*loglik.values[model.index]+multiplier.range * K.values[model.index], col=my.palette[model.index])	
	text(max(multiplier.range), y=-2*loglik.values[model.index]+max(multiplier.range) * K.values[model.index], names(K.values[model.index]), col=my.palette[model.index], pos=4, cex=0.5)
	text(min(multiplier.range), y=-2*loglik.values[model.index]+min(multiplier.range) * K.values[model.index], names(K.values[model.index]), col=my.palette[model.index], pos=2, cex=0.5)
}

abline(v=2, lty="dotted") #line for regular AIC


#We get SE for the optima (see nodeBased.OUMV$theta) but not for the other parameters. Let's see how hard they are to estimate. 
#First, look at ?OUwie.fixed to see how to calculate likelihood at a single point.
?OUwie.fixed

#Next, keep all parameters but alpha at their maximum likelihood estimates (better would be to fix just alpha and let the others optimize given this constraint, but this is harder to program for this class). Try a range of alpha values and plot the likelihood against this.
alpha.values<-seq(from=0.00001, to=2, length.out=50)

#keep it simple (and slow) and do a for loop:
likelihood.values <- rep(NA, length(alpha.values))
for (iteration in sequence(length(alpha.values))) {
	likelihood.values[iteration] <- OUwie.fixed(tree, trait, model="OUMV", alpha=rep(alpha.values[iteration],2), sigma.sq=best$solution[2,], theta=best$theta[,1])$loglik
}

plot(x=alpha.values, y=likelihood.values, xlab="Alpha", ylab="Neg lnL", type="l", bty="n")

points(x=best$solution[1,1], y=best$loglik, pch=16, col="red")
text(x=best$solution[1,1], y=best$loglik, "unconstrained best", pos=4, col="red")

abline(h=best$loglik-2, lty="dotted") #Two log-likelihood 

#Now, let's try looking at both theta parameters at once, keeping the other parameters at their MLEs
library("akima")

nreps<-400
theta1.points<-c(best$theta[1,1], rnorm(nreps-1, best$theta[1,1], 5*best$theta[1,2])) #center on optimal value, have extra variance
theta2.points<-c(best$theta[2,1], rnorm(nreps-1, best$theta[2,1], 5*best$theta[2,2])) #center on optimal value, have extra variance
likelihood.values<-rep(NA,nreps)

for (iteration in sequence(nreps)) {
	likelihood.values[iteration] <- OUwie.fixed(tree, trait, model="OUMV", alpha=best$solution[1,], sigma.sq=best$solution[2,], theta=c(theta1.points[iteration], theta2.points[iteration]))$loglik
}
#think of how long that took to do 400 iterations. Now remember how long the search took (longer).

likelihood.differences<-(-(likelihood.values-max(likelihood.values)))

#We are interpolating here: contour wants a nice grid. But by centering our simulations on the MLE values, we made sure to sample most thoroughly there
interpolated.points<-interp(x=theta1.points, y=theta2.points, z= likelihood.differences, linear=FALSE, extrap=TRUE, xo=seq(min(theta1.points), max(theta1.points), length = 400), yo=seq(min(theta2.points), max(theta2.points), length = 400))
	
contour(interpolated.points, xlim=range(c(theta1.points, theta2.points)),ylim=range(c(theta1.points, theta2.points)), xlab="Theta 1", ylab="Theta 2", levels=c(2,5,10),add=FALSE,lwd=1, bty="n", asp=1)

points(x=best$theta[1,1], y=best$theta[2,1], col="red", pch=16)

points(x=trait$X[which(trait$Reg==1)],y=rep(min(c(theta1.points, theta2.points)), length(which(trait$Reg==1))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 1, plotted along x axis
points(y=trait$X[which(trait$Reg==2)],x=rep(min(c(theta1.points, theta2.points)), length(which(trait$Reg==2))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 2, plotted along y axis


#The below only works if the discrete trait rate is low, so you have a good chance of estimating where the state is.
#If it evolves quickly, hard to estimate where the regimes are, so some in regime 1 are incorrectly mapped in
#regime 2 vice versa. This makes the models more similar than they should be.
#See Revell 2013, DOI:10.1093/sysbio/sys084 for an exploration of this effect.
library(phytools)
trait.ordered<-data.frame(trait[,2], trait[,2],row.names=trait[,1])
trait.ordered<- trait.ordered[tree$tip.label,]
z<-trait.ordered[,1]
names(z)<-rownames(trait.ordered)
tree.mapped<-make.simmap(tree,z,model="ER",nsim=1)
leg<-c("black","red")
names(leg)<-c(1,2)
plotSimmap(tree.mapped,leg,pts=FALSE,ftype="off", lwd=1)

simmapBased<-OUwie(tree.mapped,trait,model="OUMV", simmap.tree=TRUE, diagn=FALSE)
print(simmapBased)
#How does this compare to our best model from above? Should they be directly comparable?
print(best)


#Bonus: How could you do a model that has BM on some branches, OU on others?