### Simulating microsatellite data with forward time simulations:

library(rmetasim)
library(adegenet)
library(pegas)
source("rmetasim2adegenet.R")

ne <- function(t, mu) {
  x <- t / (4*mu)
return(x)
} 

exampleS <- matrix(c(0.5, 0, 0.5, 0.3), byrow=T, nrow = 2)
exampleR <- matrix(c(0, 1.1, 0, 0), byrow=T, nrow = 2)
exampleM <- matrix(c(0, 0, 0, 1), byrow=T, nrow = 2)

landscapeS <- matrix(0,2,2)
landscapeR <- landscapeS
landscapeM <- landscapeS

mu <- 0.01

dummy <- landscape.new.empty()
dummy <- landscape.new.intparam(dummy,s=2,cg=0,ce=0,totgen=1000)
dummy <- landscape.new.floatparam(dummy, s=0)
dummy <- landscape.new.switchparam(dummy, mp=0)
dummy <- landscape.new.local.demo(dummy, S=exampleS, R=exampleR, M=exampleM)
dummy <- landscape.new.epoch(dummy,S=landscapeS,R=landscapeR,M=landscapeM,extinct=c(0),carry=(1000)) #carry limits final popsize
dummy <- landscape.new.locus(dummy,type=1,ploidy=2,transmission=0,mutationrate=mu,numalleles=1) # the number of alleles will not increase if mutationrate=0
dummy <- landscape.new.individuals(dummy,PopulationSizes=c(50,50)) # need a popsize for each stage in each population, two stages here.

dummy.gen0 <- landscape.simulate(dummy,500) #simulate 500 time-clicks -> equilibrium?
dummy.gen1 <- landscape.simulate(dummy.gen0, 1)

#make the output into an adegenet object
g0 <- landscape.make.genind(dummy.gen0)
g1 <- landscape.make.genind(dummy.gen1)
g1@pop.names <- "2"
g01 <- repool(g0, g1)
g01@pop.names
