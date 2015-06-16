### Simulating microsatellite data with forward time simulations:

setwd("/Users/christineewers/Desktop/EffSize_package/EffSize_hackathon")

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
dummy <- landscape.new.epoch(dummy,S=landscapeS,R=landscapeR,M=landscapeM,extinct=c(0),carry=(1000)) #carry limits final popsize, the larger carry, the more alleles
dummy <- landscape.new.locus(dummy,type=1,ploidy=2,transmission=0,mutationrate=mu,numalleles=1) # the number of alleles will not increase if mutationrate=0
dummy <- landscape.new.locus(dummy,type=1,ploidy=2,transmission=0,mutationrate=mu,numalleles=1) # the number of alleles will not increase if mutationrate=0
dummy <- landscape.new.individuals(dummy,PopulationSizes=c(50,50)) # need a popsize for each stage in each population, two stages here.

dummy.gen0 <- landscape.simulate(dummy, 500) #simulate 500 time-clicks -> equilibrium? if I simulate 1000 generations, there are no differences in allele freq -> why?
dummy.gen1 <- landscape.simulate(dummy.gen0, 1) # advance to the next generation (+1)
dummy.gen2 <- landscape.simulate(dummy.gen1, 1) # advance to the next generation (+1)


#make the output into an adegenet object
g0 <- landscape.make.genind(dummy.gen0)
g1 <- landscape.make.genind(dummy.gen1)
g2 <- landscape.make.genind(dummy.gen2)
g0@pop.names <- "0"
g1@pop.names <- "1"
g2@pop.names <- "2"
g <- repool(g0, g1, g2)
g@loc.nall # 2 loci: 6 and 5 alleles



# subsample 50 individuals 
g0S <- g0[sample(1:999, size=50, replace=F)]
g1S <- g1[sample(1:999, size=50, replace=F)]
g2S <- g2[sample(1:999, size=50, replace=F)]
simG3 <-repool(g0S, g1S, g2S)
simG <- g0S 
