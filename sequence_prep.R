# Workflow for sequence data
setwd("/Users/christineewers/Desktop")

# import data as DNAbin object
library(ape)
dna.align <- read.DNA()
data(Laurasiatherian)
dna.align <- Laurasiatherian

# estimate genealogy
# UPGMA, neighbour joining, bio-nj and fast ME methods of phylogenetic reconstruction
library(ape)
dna.dist <- dist.ml(dna.align)
dna.tree <- NJ(dna.dist)
plot(dna.tree)
# distance, parsimony, and likelihood methods of phylogenetic reconstruction
library(phangorn)
#modelTest(dna.align, dna.tree)
treeRA <- random.addition(dna.align)
treeNNI <- optim.parsimony(treeRA, dna.align)
treeRatchet <- pratchet(dna.align, start=treeNNI)
dna.tree2 <- treeRatchet
# wrapper for Mr. Bayes, BEAST, RAxML
library(ips)

# convert DNAbin object to genind object
library(mmod)
dna.genind <- as.genind.DNAbin(dna.align)
library(adegenet)
dna.genind2 <- DNAbin2genind(dna.align)

# calculate coalescent Ne
library(ape)
dna.tree <- read.tree()
# input: phylo object

# calculate other Ne estimates
# input: genind object dna.genind