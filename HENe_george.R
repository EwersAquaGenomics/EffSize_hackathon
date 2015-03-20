###### METHOD 2: Heterozygote excess Ne ######

## Christine ##

# problems: Ne estimate is too small - where is the error in my formula? Is it low frequency alleles?
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
  i=4
  for (i in 1:K){
    sample_sizes[i]<-sum(na.rm=T,genind_obj@tab[,starts[i]:ends[i]])
  }
  return(sample_sizes)
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

HENe <- function(genind_obj) {
  
  K <- length(genind_obj@all.names)
  genind.by.locus <- seploc(genind_obj, res.type="matrix", truenames=F)
  sample_sizes <- getSampleSize_byLocus(genind_obj)
  
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
  
  
  #i <- 1
  #j <- 1
  for (j in 1:K) {
    L<-ncol(genind.by.locus[[j]])
    num.alleles.per.locus[j] <- sum(genind.by.locus[[j]], na.rm=T)*2 # total number of alleles for locus 1
    num.allele[[j]]<-p[[j]]<-Hexp[[j]]<-Hobs[[j]]<-no.het.per.allele[[j]]<-d[[j]]<-numeric(L)
    for (i in 1:L) {    
      num.allele[[j]][i] <- sum(genind.by.locus[[j]][,i], na.rm=T)*2 # number of occurrences of allele 1 of locus 1
      p[[j]][i] <- num.allele[[j]][i] / num.alleles.per.locus[j]
      Hexp[[j]][i] <- 2*p[[j]][i]*(1-p[[j]][i])*2*sample_sizes[j] / (2*sample_sizes[j]-1)    
      no.het.per.allele[[j]][i] <- sum(!is.na(temp_genind.locus) & temp_genind.locus == 0.5)
      Hobs[[j]][i] <- no.het.per.allele[[j]][i] / num.allele[[j]][i]
      d[[j]][i] <- (Hobs[[j]][i]-Hexp[[j]][i]) / Hexp[[j]][i]
    }
    
    w[j] <- sqrt(sample_sizes[j])*(genind_obj@loc.nall[j]-1) / genind_obj@loc.nall[j]
    d.per.locus[j] <- sample_sizes[j]*sum(d[[j]])
    dw[j] <- d.per.locus[j]*w[j]
  }
  
  D <- sum(dw)/sum(w)
  Nb <- 1/2*D + 1/(2*(D+1))
  Nb
} 
