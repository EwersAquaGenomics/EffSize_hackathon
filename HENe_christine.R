###### METHOD 2: Heterozygote excess Ne ######

## Christine ##

# missing: removing samples with allele frequency < X
# missing: confidence intervals / bootstrapping / jackknifing

#l=1
for (l in 1:length(genind.crit)) {
genind.by.locus <- genind.crit[[l]]
K <- length(genind.by.locus)
sample_sizes <- getSampleSize_byLocus(genObj)

p <- num.allele <- Hexp <- no.het.per.allele <- Hobs <- d <- vector("list", length=K)
w <- d.per.locus <- dw <- num.alleles.per.locus <- 0

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
    w[j] <- sqrt(sample_sizes[j])*(nancycats@loc.nall[j]-1) / nancycats@loc.nall[j]
    d.per.locus[j] <- 1/sample_sizes[j]*sum(d[[j]])
    dw[j] <- d.per.locus[j]*w[j]
  }
}

D[l] <- sum(dw)/sum(w)
Nb[l] <- 1/2*D + 1/(2*(D+1))

}