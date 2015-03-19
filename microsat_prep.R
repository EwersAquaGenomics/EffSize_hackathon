# Workflow for microsatellite data
setwd("/Users/christineewers/Desktop")

# import microsatellite data
library(adegenet)
micro.df <- read.table()

# convert dataframe into genind object
micro.genind <- df2genind(micro.df, sep=NULL, ncode=NULL, ind.names=NULL, loc.names=NULL,
                          pop=NULL, missing=NA, ploidy=2, type="codom")

# claculate simple distance matrix and estimate genealogy
# alternatively, the user can input trees into te coalescent Ne function
library(mmod)
micro.dist <- dist.codom(micro.genind)
library(phangorn)
micro.tree <- NJ(micro.dist)
 
# calculate coalescent Ne
# input: phylo object

# calculate other Ne estimates
# input: genind object micro.genind