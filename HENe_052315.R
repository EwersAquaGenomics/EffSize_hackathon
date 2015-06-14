###### METHOD 2: Heterozygote excess Ne ######

# missing: removing samples with allele frequency < X
# missing: confidence intervals / bootstrapping

#' @title getSampleSize_byLocus 
#' @description Takes a genind object and returns the number of samples occurring at each locus
#' 
#' @author George Shirreff <georgeshirreff@@gmail.com>
#' 
#' @examples library(adegenet)
#' @examples data(nancycats)
#' @examples getSampleSize_byLocus(nancycats)
#' 
getSampleSize_byLocus <- function(genind_obj)
{
  K<-length(genind_obj@all.names)
  
  ends<-cumsum(genind_obj@loc.nall)
  starts<-c(0,ends[-length(ends)])+1
  
  sample_sizes<-numeric(K)
  #i=4
  for (i in 1:K){
    sample_sizes[i]<-sum(na.rm=T,genind_obj@tab[,starts[i]:ends[i]])
  }
  return(sample_sizes)
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
makeCrit <- function(genObj, crit) {
  
  # testing the function
  #genObj <- nancycats
  #j <- 1
  #crit <- 0.05
  
  genind.by.locus <- seploc(genObj) # separates the genind object by locus
  p <- makeP(genind.by.locus) # makes list with each list element containing the allele frequencies for one locus
  K <- length(genind.by.locus) # number of loci
  pp <- 0
  
  for (j in 1:K) {
    pp <- p[[j]] >= crit  # boolean vector (TRUE FALSE TRUE etc)
    genCrit <- genind.by.locus
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


#' @title HENe 
#' @description Calculates population size from a genind object by heterozygote excess
#' 
#' @author Christine Ewers-Saucedo <ewers.christine@@gmail.com>, George Shirreff <georgeshirreff@@gmail.com>
#' 
#' @examples library(adegenet)
#' @examples data(nancycats)
#' @examples HENe(nancycats)
#' 

library(adegenet)
data(nancycats)
genind_obj <- nancycats

HENe <- function(genind_obj, crit) {
  
  genind.by.locus <- makeCrit(genind_obj, crit)
  K <- length(genind.by.locus) # number of loci
  
  sample_sizes <- numeric(K)
  for (i in 1:K) {
    sample_sizes[i] <- sum(genind.by.locus[[i]]@tab, na.rm=T) # alleles sampled per locus 
  }
  
  
  # generating empty files
  num.alleles.per.locus <- numeric(K)
  w <- numeric(K)
  d.per.locus <- numeric(K)
  dw <- numeric(K)
  
  p <- vector("list", length=K)
  num.allele <- vector("list", length=K)
  Hexp <- vector("list", length=K)
  no.het.per.allele <- vector("list", length=K)
  Hobs <- vector("list", length=K)
  d <- vector("list", length=K)
  
  for (j in 1:K) {
    L <- ncol(genind.by.locus[[j]]@tab) # number of alleles per locus
    num.alleles.per.locus[j] <- sum(genind.by.locus[[j]]@tab, na.rm=T) # total number of alleles for locus 1
    num.allele[[j]]<-p[[j]]<-Hexp[[j]]<-Hobs[[j]]<-no.het.per.allele[[j]]<-d[[j]]<-numeric(L)
    for (i in 1:L) {    
      num.allele[[j]][i] <- sum(genind.by.locus[[j]]@tab[,i], na.rm=T) # number of occurrences of allele i of locus j
      p[[j]][i] <- num.allele[[j]][i] / num.alleles.per.locus[j] # allele frequency
      Hexp[[j]][i] <- 2*p[[j]][i]*(1-p[[j]][i])*sample_sizes[j] / (sample_sizes[j]-1) # not totally sure why this is normalized by sample size -> weighing?    
      no.het.per.allele[[j]][i] <- sum(genind.by.locus[[j]]@tab[,i] == 1, na.rm=T)
      Hobs[[j]][i] <- no.het.per.allele[[j]][i] / sample_sizes[j]
      d[[j]][i] <- (Hobs[[j]][i]-Hexp[[j]][i]) / Hexp[[j]][i]
    }
    
    w[j] <- sqrt(sample_sizes[j]/2)*(genind_obj@loc.nall[j]-1) / genind_obj@loc.nall[j] # weights are the number of alleles per locus
    d.per.locus[j] <- sample_sizes[j]/2*sum(d[[j]])
    dw[j] <- d.per.locus[j]*w[j]
  }
  
  D <- sum(dw)/sum(w)
  Nb <- 1/2*D + 1/(2*(D+1))
  return(Nb)
} 
