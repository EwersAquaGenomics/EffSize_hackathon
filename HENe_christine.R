###### METHOD 2: Heterozygote excess Ne ######

## Christine ##

# missing: confidence intervals / bootstrapping / jackknifing


getSampleSize_byLocus <- function(theObject)
{
  K<-length(genindObj@all.names)
  
  ends<-cumsum(theObject@loc.nall)
  starts<-c(0,ends[-length(ends)])+1
  
  sample_sizes<-numeric(K)
  i=4
  for (i in 1:K){
    sample_sizes[i]<-sum(na.rm=T,theObject@tab[,starts[i]:ends[i]])
  }
  return(sample_sizes)
}


require(adegenet)

for (l in 1:length(genind.crit)) {
genind.by.locus <- genind.crit[[l]]
K <- length(genind.by.locus)
sample_sizes <- getSampleSize_byLocus(genObj)

p <- num.allele <- Hexp <- no.het.per.allele <- Hobs <- d <- vector("list", length=K)
w <- d.per.locus <- dw <- num.alleles.per.locus <- 0

p <- vector("list", length=K)
num.allele <- vector("list", length=K)
num.alleles.per.locus <- 0
Hexp <- vector("list", length=K)
no.het.per.allele <- vector("list", length=K)
Hobs <- vector("list", length=K)
d <- vector("list", length=K)
w <- 0
d.per.locus <- 0
dw <- 0


#j <- 1
#i <- 1

for (j in 1:K) {
 for (i in 1:ncol(genind.by.locus[[j]]@tab)) {
    num.alleles.per.locus[j] <- sum(genind.by.locus[[j]]@tab, na.rm=T) # total number of alleles for locus 1
    num.allele[[j]][i] <- sum(genind.by.locus[[j]]@tab[,i], na.rm=T) # number of occurrences of allele 1 of locus 1
    p[[j]][i] <- num.allele[[j]][i] / num.alleles.per.locus[j] # allele frequency
    Hexp[[j]][i] <- 2*p[[j]][i]*(1-p[[j]][i])*2*sample_sizes[j] / (2*sample_sizes[j]-1) # expected heterozygosity
    no.het.per.allele[[j]][i] <- length(genind.by.locus[[j]]@tab[,i][!is.na(genind.by.locus[[j]]@tab[,i]) & genind.by.locus[[j]]@tab[,i] == 1])


    Hobs[[j]][i] <- no.het.per.allele[[j]][i] / num.allele[[j]][i]
    d[[j]][i] <- (Hobs[[j]][i]-Hexp[[j]][i]) / Hexp[[j]][i]
  }
}


D[l] <- sum(dw)/sum(w)
Nb[l] <- 1/2*D + 1/(2*(D+1))

}

