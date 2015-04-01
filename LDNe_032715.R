### LDNe ###

# implemented: accounting for missing data by recalculating S and dropping low frequency alleles (R script: )
# missing: confidence intervals, bootstrapping or jackknifing (implemented in NeEstimator)

require(pegas)
require(adegenet)
#source("dropLowFreqAlleles.R")

# generate Ne values for each subset of data with varying critical values of allele frequenceis
for (l in 1:length(crit)) {
  genind.by.locus <- genind.crit[[l]]
  lociObj.ls <- vector("list", length=length(genind.by.locus))
  for (i in 1:length(genind.by.locus)) {
    lociObj.ls[[i]] <- genind2loci(genind.by.locus[[i]]) # list of loci objects, one per locus
    if (i == 1) (lociObj <- lociObj.ls[[1]]) else (lociObj <- cbind(lociObj, lociObj.ls[[i]]))
  }
  
  
  K <- length(genObj@all.names)
  no.comparisons <- K*(K-1)/2


### calculating numbers and frequencies of alleles per locus, of heterozygotes and homozygotes of alleles per locus
# the following is also needed in HENe, so maybe this should be more general calculations?

p <- num.allele <- het.per.allele <- hom.per.allele <- p.het <- p.hom <- vector("list", length=K)
num.alleles.per.locus <- 0

for (j in 1:K) {
  for (i in 1:ncol(genind.by.locus[[j]])) {
    num.alleles.per.locus[j] <- sum(genind.by.locus[[j]], na.rm=T) # total number of genotyped alleles for locus 1, not the number of unique alleles
    num.allele[[j]][i] <- sum(genind.by.locus[[j]][,i], na.rm=T) # number of occurrences of allele 1 of locus 1
    p[[j]][i] <- num.allele[[j]][i] / num.alleles.per.locus[j] # allele frequency
    het.per.allele[[j]][i] <- length(genind.by.locus[[j]][,i][genind.by.locus[[j]][,i] == 1])
    hom.per.allele[[j]][i] <- length(genind.by.locus[[j]][,i][genind.by.locus[[j]][,i] == 2])
    p.het[[j]][i] <- het.per.allele[[j]][i] / num.alleles.per.locus[j]
    p.hom[[j]][i] <- hom.per.allele[[j]][i] / num.alleles.per.locus[j]
  }
}


### calculating delta for all locus combinations

delta.ls <- vector("list", length=paste(K,K-1,sep=""))
counter <- S.loci <- n.loci <- 0

#i <- 1
for (i in 1:(K-1)) {
  j <- i
  repeat {  
    j <- j+1
    delta.ls[[paste(i,j,sep="")]] <- LD2(lociObj, locus=c(i,j))$Delta
    counter[paste(i,j,sep="")] <- paste(i,j,sep="")
    S.loci[paste(i,j,sep="")] <- length(which(!is.na(lociObj[i]) & !is.na(lociObj[j]))) # number of genotyped individuals for each locus pair
    n.loci[paste(i,j,sep="")] <- (ncol(genind.by.locus[[i]])-1)*(ncol(genind.by.locus[[j]])-1)
    if (j==K) break
  }
}
# takes a second or so...

w.loci <- n.loci*S.loci^2 # huge numbers

delta <- delta.ls[lapply(delta.ls,length)>0]
#length(delta) #36
delta$"36" # this means information on which combination of loci is where in the ls is retained


### calculating delta hat

delta.hat <- vector("list", length=no.comparisons)

for (i in 1:no.comparisons) {
  delta.hat[[i]] <- delta[[i]]*S.loci[i]/(S.loci[i]-1)
}
names(delta.hat) <- names(delta)


### calculating r hat

#r.hat <- r.hat.denom <- vector("list", length=paste(K,K-1,sep=""))
r.hat.denom <- r.hat2 <- delta.hat # is this problematic for NA?

for (i in 1:(K-1)) { # locus i
  j <- i 
  repeat {  
    j <- j+1 # locus j
    for (k in 1:ncol(genind.by.locus[[i]])) { # different alleles in locus i
      for (l in 1:ncol(genind.by.locus[[j]])) { # alleles of locus j
        r.hat.denom[[paste(i,j,sep="")]][k,l] <- sqrt((p[[i]][k]*(1-p[[i]][k])+(p.hom[[i]][k]-p[[i]][k]^2))*
                                             (p[[j]][l]*(1-p[[j]][l])+(p.hom[[j]][l]-p[[j]][l]^2)))
      }
    }
    if (j==K) break
  }
}
# There were 50 or more warnings (use warnings() to see the first 50)
# the warnings are all "Na introduced". Which is expected given missing data
# but the results look fine...

r.hat2.loci <- 0

for (i in 1:no.comparisons) {
   r.hat[[i]] <-  delta.hat[[i]]/r.hat.denom[[i]]
   r.hat2[[i]] <- r.hat[[i]]*r.hat[[i]]
   r.hat2.loci[i] <- mean(r.hat2[[i]], na.rm=T) # mean for each pair of loci
}

# taken from NeCalc.pdf
r2.mean <- sum(r.hat2.loci*w.loci[2:(no.comparisons+1)])/sum(w.loci) #-0.003870057

# weighted harmonic mean <- sum(x*weights)/sum(weight)) for each x
# x is S.loci
# weights are proportional to n.loci
S <- sum(S.loci*n.loci)/sum(n.loci) #225.76


# calculate Ne under random mating:
within.sqrt <- 0
mating <- "random"

if (mating == "random") {
if (S>=30) {
  expect.r2 <- 1/S + 3.19/S^2
  r2.drift <- r2.mean - expect.r2
  within.sqrt <- 0.1111111 - 2.76*r2.drift
  #if (within.sqrt < 0) (within.sqrt <- 0) # set the value within the sqrt to zero if < 0
  Ne <- 1/(2*r2.drift) * (0.3333333 + sqrt(within.sqrt))
} else {
  expect.r2 <- 0.0018 + 0.907/S + 4.44/S^2
  r2.drift <- r2.mean - expect.r2
  within.sqrt <- 0.094864 - 2.08*r2.drift # set the value within the sqrt to zero if < 0
  #if (within.sqrt < 0) (within.sqrt <- 0)
  Ne <- 1/(2*r2.drift) * (0.308 + sqrt(within.sqrt))
}      
}

mating <- "mono"

if (mating == "mono") {
if (S>=30) {
  expect.r2 <- 1/S + 3.19/S^2
  r2.drift <- r2.mean - expect.r2
  within.sqrt <- 0.4444444 - 7.2*r2.drift
  #if (within.sqrt < 0) (within.sqrt <- 0) # set the value within the sqrt to zero if < 0
  Ne <- 1/(2*r2.drift) * (0.6666666 + sqrt(within.sqrt))
} else {
  expect.r2 <- 0.0018 + 0.907/S + 4.44/S^2
  r2.drift <- r2.mean - expect.r2
  within.sqrt <- 0.381924 - 5.24*r2.drift
  #if (within.sqrt < 0) (within.sqrt <- 0) # set the value within the sqrt to zero if < 0
  Ne <- 1/(2*r2.drift) * (0.618 + sqrt(within.sqrt))
}      
}

Ne.crit[l] <- Ne
}

