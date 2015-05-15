rm(list=ls()) #let's tidy up our workspace
library(corHMM)

#example dataset
data(primates)

#what's in this?
ls()

#ah, it's a tree (ape phylo) and a trait data.frame

#let's consider the tree:
print(primates$tree)

#and the data.frame:
print(primates$trait)

#note that the taxon names are stored in the first column, not the rownames. The first trait column is mating system

#Usually trees have N-1 nodes. This has N-2. Is it fully resolved?
print(ifelse(is.binary.tree(primates$tree), "Yes, this is fully resolved", "Aargh, it has one or more polytomies"))


#In many comparative methods implementations (though often not the methods themselves) this is a problem. It's fine for corHMM

one.cat <- corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=1,node.states="marginal", ip=1)

two.cat <- corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=2,node.states="marginal", ip=1) #this will take a while
#if this takes too long, stop and load the saved file corHMM.sample.RData instead: load("corHMM.sample.RData")

#How do they compare? What does the rate matrix mean? Are all the zeros estimated or forced?

#Let's plot the tree and reconstruction:
plotRECON(two.cat$phy, two.cat$states)

#Let's specify a root state
one.cat.root.0 <- corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=1,node.states="marginal", ip=1, root.p=c(1,0))


#Now do the same for state 1 (out of 0,1) at the root
one.cat.root.1 <- ???
stop("replace the ??? and delete this stop")

#And equilibrium:
one.cat.root.maddfitz <- ???
stop("replace the ??? and delete this stop")


#How does this affect rates?
print(one.cat)
print(???)
print(???)
print(???)
stop("replace the ??? and delete this stop")

#how about ancestral state reconstructions?
par(mfrow=c(2,2))
plotRECON(one.cat$phy, one.cat$states, title="Flat freq", pie.cex=.8)
plotRECON(one.cat.root.maddfitz$phy, one.cat.root.maddfitz$states, title="Equilibrium freq", pie.cex=.8)
plotRECON(one.cat.root.0$phy, one.cat.root.0$states, title="Root 0", pie.cex=.8)
plotRECON(one.cat.root.1$phy, one.cat.root.1$states, title="Root 1", pie.cex=.8)