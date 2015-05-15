#Brian O'Meara
#Released under CC0 license, though attribution is welcome

#First, let's get relevant packages installed.
install.packages("ctv") #this is for task views
library(ctv)
install.views("Phylogenetics") #Installs all packages in this area

#One of the major utility packages in phylogenetics is ape (Paradis et al)
library(ape)

#Generate a tree from a simulation
phy <- rbdtree(birth=1, death=0, Tmax=4)

#Get a short summary of the tree:
print(phy)

#Get a bit longer summary:
summary(phy)

#Look at the tree:
plot(phy)

#Maybe we should get rid of the labels? But how to find out how to do this?  
stop("Use your R skills, then erase this line")

#Now plot without tip labels:
plot(phy, stop("Replace this with argument you learned above"))

#Maybe save the tree as a pdf?
pdf(file="SimulatedTree.pdf", width=5, height=5)
plot(phy)
dev.off()

#Now open the PDF. If on a mac, you can do
system("open SimulatedTree.pdf")
#If not, open the PDF the way you normally open files

#Write your tree to a file:
write.tree(phy, file="SimulatedTree.phy")

#To see the tree format, can look at it here or the file:
write.tree(phy, file="")

#This nested format is named Newick, after Newick's Lobster House (http://www.newicks.com), where the format was developed. Think of it as the Hennig set diagram, but with the top and bottom of the ovals erased, leaving only the sides, as parentheses.

#Read the tree from the file
phy2 <- read.tree(file="SimulatedTree.phy")

#The beauty of R is that it makes even complex analyses pretty simple. For example, imagine we wanted to know how many tips a given birth death rate is expected to have given time. We could do math, but that's hard (ascertainment bias, for one thing). So let's just simulate the heck out of it.

GetNtipFromSim <- function(Tmax=4, birth=1, death=0) {
  phy <- rbdtree(birth=birth, death=death, Tmax=Tmax)
  return(Ntip(phy))
}

nreps <- 100
Tmax <- 4
distribution.of.tips <- replicate(nreps, GetNtipFromSim(Tmax=Tmax))

plot(density(distribution.of.tips, from=0), xlab="Number of tips", bty="n", main=paste("Number of tips at Tmax =", Tmax))

#Now a line to show the median
abline(v=median(distribution.of.tips))

#How would you do this for Tmax=2?

#Could you plot them on the same plot (hint: ?lines)

