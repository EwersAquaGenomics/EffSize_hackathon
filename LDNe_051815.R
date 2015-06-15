### LDNe ###

require(pegas)
require(adegenet)

#' @title Allele frequencies per locus
#' @description Takes a genind object and calculates the frequency of each allele per locus
#' 
#' @author Christine Ewers-Saucedo <ewers.christine@@gmail.com>
#' 
#' @param genind.by.locus list of genind objects, one list element per locus
#' 
#' @examples library(adegenet)
#' @examples data(nancycats)
#' @examples makeP(nancycats)
# testing functions
#genind.by.locus <- seploc(simG)

makeP <- function(genind.by.locus) {
  
  #genind.by.locus <- seploc(genObj)
  K <- length(genind.by.locus) # number of loci
  p <- occ.allele <- vector("list", length=K)
  alleles.per.locus <- 0
  sample.loc <- 0
  
  for (j in 1:K) {
    
    sample.loc[j] <- sum(genind.by.locus[[j]]@tab, na.rm=T)  # total number of genotyped alleles for locus, not the number of unique alleles
    
    for (i in 1:ncol(genind.by.locus[[j]]@tab)) {
      occ.allele[[j]][i] <- sum(genind.by.locus[[j]]@tab[,i], na.rm=T) # number of occurrences of allele
      p[[j]][i] <- occ.allele[[j]][i] / sample.loc[j] # allele frequency
    }
  }
  return(p)
}


#' @title Removing alleles below critical frequency
#' @description Takes a genind object, removes alleles below a certain allele frequency, and returns a list of genind objects, one object per locus
#' 
#' @author Christine Ewers-Saucedo <ewers.christine@@gmail.com>
#' 
#' @param genObj genind object
#' @param crit the lowest allele frequency allowed to occur. Alleles at a lower frequency are excluded.
#' 
#' @examples library(adegenet)
#' @examples data(nancycats)
#' @examples makeCrit(nancycats)
#' 

# testing the function
#genObj <- simG
#j <- 1
#crit <- 0.02

makeCrit <- function(genObj, crit=0.01) {
  
  genind.by.locus <- seploc(genObj) # separates the genind object by locus
  genCrit <- genind.by.locus
  p <- makeP(genind.by.locus) # makes list with each list element containing the allele frequencies for one locus
  K <- length(genind.by.locus) # number of loci
  pp <- 0
  
  for (j in 1:K) {
      pp <- p[[j]] >= crit  # boolean vector (TRUE FALSE TRUE etc)
      genCrit[[j]]@tab <- genind.by.locus[[j]]@tab[,pp] # subsets each locus to include only alleles above critical allele frequency
      rs <- rowSums(genCrit[[j]]@tab, na.rm = T)
      for (k in 1:nrow(genCrit[[j]]@tab)) {
        if ((rs[k]) < 2) {
          genCrit[[j]]@tab[k,] <- NA      
        }
      }
      genCrit[[j]]@loc.nall <- ncol(genCrit[[j]]@tab)
      genCrit[[j]]@loc.fac <- genind.by.locus[[j]]@loc.fac[pp]
      genCrit[[j]]@all.names[[1]] <- genind.by.locus[[j]]@all.names[[1]][pp]
  }
return(genCrit)
}


#' @title Calculating point estimate of linkage disequilibrium effective population size
#' @description Takes a genind object and calculates effective population size based on the magnitude of linkage disequilibrium found in a sample
#' 
#' @author Christine Ewers-Saucedo <ewers.christine@@gmail.com>
#' 
#' @param genObj genind object
#' @param crit the lowest allele frequency used in the calculation. Alleles with lower allele frequency are excluded
#' @param mating mating system information
#' 
#' @examples library(adegenet)
#' @examples data(nancycats)
#' @examples LDNe(nancycats, crit=0.05, mating="random")
#' 

# to test the function, I executed the lines below and ran through each part of the code without making the function (omitted first and last line)
#genObj <- simG
#crit <- 0.01
#mating <- "random

LDNe_point <- function(genObj, crit=0.01, mating="random") {
  
  genCrit <- makeCrit(genObj, crit) # a list of genind objects, one locus per list element
  
  K <- length(genCrit) # number of loci
  no.comparisons <- K*(K-1)/2 
  p <- makeP(genCrit)


  ### calculating numbers and frequencies of alleles per locus, of heterozygotes and homozygotes of alleles per locus

  p.hom <- vector("list", length=K)
  num.alleles.per.locus <- hom.per.allele <- 0

  for (j in 1:K) {
    num.alleles.per.locus[j] <- sum(genCrit[[j]]@tab, na.rm=T)
    for (i in 1:ncol(genCrit[[j]]@tab)) {
      hom.per.allele <- sum(genCrit[[j]]@tab[,i] == 2, na.rm=T)
      p.hom[[j]][i] <- hom.per.allele / num.alleles.per.locus[j]
    }
  }


  ### calculating delta for all locus combinations
  
  # turns each genind list element into a loci object
  ploidy <- 0
  bb <- vector("list", length=length(genCrit))
  for (i in 1:length(genCrit)) {
    bb[[i]] <- as.loci(genCrit[[i]])
    #ploidy[i] <- getPloidy(bb[[i]])
  }
  # merge loci list into a single loci object
  lociObj <- bb[[1]] # generates the first element of lociObj - starting point for cbind
  for (i in 2:length(bb)) {
    lociObj <- cbind(lociObj, bb[[i]])
  }
  #getPloidy(lociObj)

  delta.ls <- vector("list", length=paste(K,K-1,sep=""))
  counter <- S.loci <- n.loci <- 0
  
  for (i in 1:(K-1)) {
    j <- i
    repeat {  
      j <- j+1
      delta.ls[[paste(i,j,sep="")]] <- LD2(lociObj, locus=c(i,j))$Delta # matrix
      counter[paste(i,j,sep="")] <- paste(i,j,sep="")
      S.loci[paste(i,j,sep="")] <- length(which(!is.na(genCrit[[i]]@tab[,1]) & (!is.na(genCrit[[j]]@tab[,1])))) # number of genotyped individuals for each locus pair
      n.loci[paste(i,j,sep="")] <- (ncol(genCrit[[i]]@tab)-1)*(ncol(genCrit[[j]]@tab)-1) # number of allele combinations for each locus pair
      if (j==K) break
    }
  }
  
  S.loci <- S.loci[-1] # vector, number of samples (individuals)
  n.loci <- n.loci[-1] # vector, number of alleles per comparison
  w.loci <- n.loci*S.loci^2 # vectors, huge numbers
  delta <- delta.ls[lapply(delta.ls,length)>0] # list, one matrix per comparison, many NA, rows and columns are alleles
  
  ### calculating delta hat
  delta.hat <- vector("list", length=no.comparisons) # one matrix per comparison, rows and columns are alleles
  for (i in 1:no.comparisons) {
    delta.hat[[i]] <- delta[[i]]*S.loci[i]/(S.loci[i]-1) # weighing delta, numbers similar to delta
  }
  names(delta.hat) <- names(delta)

  ### calculating r hat
  r.hat.denom <- vector("list", length=paste(K,K-1,sep=""))
  
  for (i in 1:(K-1)) { # locus i
    j <- i 
    repeat {  
      j <- j+1 # locus j
      r.hat.denom[[paste(i,j,sep="")]] <- matrix(0, nrow=ncol(genCrit[[i]]@tab), ncol=ncol(genCrit[[j]]@tab))
      for (k in 1:ncol(genCrit[[i]]@tab)) { # different alleles in locus i
        for (l in 1:ncol(genCrit[[j]]@tab)) { # alleles of locus j
          r.hat.denom[[paste(i,j,sep="")]][k,l] <- sqrt((p[[i]][k]*(1-p[[i]][k])+(p.hom[[i]][k]-p[[i]][k]^2))*
                                             (p[[j]][l]*(1-p[[j]][l])+(p.hom[[j]][l]-p[[j]][l]^2)))
        }
      }
      if (j==K) break
    }
  }
  r.hat.denom.s <- r.hat.denom[lapply(r.hat.denom,length)>0] # list, one matrix per comarison, rows and columns are alleles
  
  r.hat <- r.hat2 <- vector("list", length=no.comparisons)
  r.hat2.loci <- 0
  for (i in 1:no.comparisons) {
    r.hat[[i]] <-  delta.hat[[i]]/r.hat.denom.s[[i]]
    r.hat2[[i]] <- r.hat[[i]]*r.hat[[i]] # do I need to take the squareroot of this?
    r.hat2.loci[i] <- mean(r.hat2[[i]], na.rm=T) # vector, mean for each pair of loci
  }

  # taken from NeCalc.pdf
  r2.mean <- sum(r.hat2.loci*w.loci, na.rm=T)/sum(w.loci, na.rm=T) # 0.01934703 ; some NaN -> problem?

  # weighted harmonic mean <- sum(x*weights)/sum(weight)) for each x
  # x is S.loci
  # weights are proportional to n.loci
  S <- sum(S.loci*n.loci)/sum(n.loci) #49


  # calculate Ne under random mating:
  within.sqrt <- 0

  if (mating == "random") {
    if (S>=30) {
      expect.r2 <- 1/S + 3.19/S^2
      r2.drift <- r2.mean - expect.r2
      within.sqrt <- 0.1111111 - 2.76*r2.drift
      if (within.sqrt < 0) (within.sqrt <- 0)# set the value within the sqrt to zero if < 0
      Ne <- 1/(2*r2.drift) * (0.3333333 + sqrt(within.sqrt))
      } else {
        expect.r2 <- 0.0018 + 0.907/S + 4.44/S^2
        r2.drift <- r2.mean - expect.r2
        within.sqrt <- 0.094864 - 2.08*r2.drift # set the value within the sqrt to zero if < 0
        if (within.sqrt < 0) (within.sqrt <- 0)
        Ne <- 1/(2*r2.drift) * (0.308 + sqrt(within.sqrt))
      }      
  }

  if (mating == "mono") {
    if (S>=30) {
      expect.r2 <- 1/S + 3.19/S^2
      r2.drift <- r2.mean - expect.r2
      within.sqrt <- 0.4444444 - 7.2*r2.drift
      if (within.sqrt < 0) (within.sqrt <- 0) # set the value within the sqrt to zero if < 0
      Ne <- 1/(2*r2.drift) * (0.6666666 + sqrt(within.sqrt))
      } else {
        expect.r2 <- 0.0018 + 0.907/S + 4.44/S^2
        r2.drift <- r2.mean - expect.r2
        within.sqrt <- 0.381924 - 5.24*r2.drift
        if (within.sqrt < 0) (within.sqrt <- 0) # set the value within the sqrt to zero if < 0
        Ne <- 1/(2*r2.drift) * (0.618 + sqrt(within.sqrt))
      }      
  }
  return(Ne)
}


# jackknifing
# modified function summarise_bootstrap (package mmod)

jack_samples <- function(genObj) {
  nsam <- nrow(genObj@tab)
  x <- vector("list", length=nsam)
  for (i in 1:nsam) {
    x[[i]] <- genObj[-i]
  }
  return(x)
}

# if output is a vector (LDNe, HENe)
sum_jack <- function (jn, statistic, crit, mating) {
  stats <- sapply(jn, statistic, crit, mating)
  return(c(median = median(stats), quantile(stats, c(0.025, 0.975), na.rm = TRUE), variance = var(stats)))
}

#genObj <- simG
#crit = 0.02
#mating = "random"
LDNe <- function(genObj, crit=0.01, mating="random") {
  Ne_point <- LDNe_point(genObj, crit, mating)
  jn <- jack_samples(genObj)
  Ne_jack <- sum_jack(jn, LDNe_point, crit, mating) # Error in delta.hat[[i]]/r.hat.denom.s[[i]] : non-conformable arrays
  return(c(point_estimate=Ne_point, Ne_jack))
}


# troubleshooting
#LDNe(simG)

#genObj <- nancycats

#for (i in 1:length(jn)) {
#  x[i] <- LDNe_point(jn[[i]])
#}
# no errors


