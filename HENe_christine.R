###### METHOD 2: Heterozygote excess Ne ######

## Christine ##

# problems: Ne estimate is too small - where is the error in my formula? Is it low frequency alleles?
# missing: removing samples with allele frequency < X
# missing: confidence intervals / bootstrapping

require(adegenet)
data(nancycats)

getSampleSize_byLocus <- function(theObject)
{
  K<-length(theObject@all.names)
  
  ends<-cumsum(theObject@loc.nall)
  starts<-c(0,ends[-length(ends)])+1
  
  sample_sizes<-numeric(K)
  i=4
  for (i in 1:K){
    sample_sizes[i]<-sum(na.rm=T,theObject@tab[,starts[i]:ends[i]])
  }
  return(sample_sizes)
}

K <- length(nancycats@all.names)
genind.by.locus <- seploc(nancycats, res.type="matrix", truenames=F)
sample_sizes <- getSampleSize_byLocus(nancycats)

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
  num.alleles.per.locus[j] <- sum(genind.by.locus[[j]], na.rm=T)*2 # total number of alleles for locus 1
  w[j] <- sqrt(sample_sizes[j])*(nancycats@loc.nall[j]-1) / nancycats@loc.nall[j]
  d.per.locus[j] <- 1/sample_sizes[j]*sum(d[[j]])
  dw[j] <- d.per.locus[j]*w[j]
  
  for (i in 1:ncol(genind.by.locus[[j]])) {    
    num.allele[[j]][i] <- sum(genind.by.locus[[j]][,i], na.rm=T)*2 # number of occurrences of allele 1 of locus 1
    p[[j]][i] <- num.allele[[j]][i] / num.alleles.per.locus[j]
    Hexp[[j]][i] <- 2*p[[j]][i]*(1-p[[j]][i])*2*sample_sizes[j] / (2*sample_sizes[j]-1)
    no.het.per.allele[[j]][i] <- length(genind.by.locus[[j]][,i][!is.na(genind.by.locus[[j]][,i]) & genind.by.locus[[j]][,i] == 0.5])
    Hobs[[j]][i] <- no.het.per.allele[[j]][i] / num.allele[[j]][i]
    d[[j]][i] <- (Hobs[[j]][i]-Hexp[[j]][i]) / Hexp[[j]][i]
  }
}

D <- sum(dw)/sum(w)
Nb <- 1/2*D + 1/(2*(D+1))

