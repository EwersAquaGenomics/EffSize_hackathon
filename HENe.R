###### METHOD 2: Heterozygote excess Ne ######

K<-length(nancycats@all.names)
sample_sizes <- getSampleSize_byLocus(nancycats)
sample_sizes <- getSampleSize_byLocus(nancycats)

genind_sum<-summary(nancycats)
genind_sum$Hobs
genind_sum$Hexp

Dindex=Hexp=list()
Dlocus<-numeric(K)
Hobs=getHeterozygoteNumbers(nancycats)

j=1
for(j in 1:K){
  p<-colSums(nancycats@tab[,starts[j]:ends[j]],na.rm=T)/sample_sizes[j]
  Hexp[[j]] <- 2*p*(1-p)*(1+1/(2*sample_sizes[j]))
  #Hexp[[j]] <- 2*p*(1-p)*(1+1/(2*sample_sizes[j]-1))  #alternative version which corrects for bias
  Dindex[[j]]=(Hobs[[j]]-Hexp[[j]])/Hexp[[j]]
  Dlocus[j]<-mean(Dindex[[j]])
}

## Method 2b ##


## Christine ##

K <- length(nancycats@all.names)
genind.by.locus <- seploc(nancycats, res.type="matrix", truenames=F)

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
sample_sizes <- getSampleSize_byLocus(nancycats)

for (j in 1:K) {
  for (i in 1:ncol(genind.by.locus[[j]])) {
    num.alleles.per.locus[j] <- sum(genind.by.locus[[j]], na.rm=T)*2 # total number of alleles for locus 1
    num.allele[[j]][i] <- sum(genind.by.locus[[j]][,i], na.rm=T)*2 # number of occurrences of allele 1 of locus 1
    p[[j]][i] <- num.allele[[j]][i] / num.alleles.per.locus[j]
    Hexp[[j]][i] <- 2*p[[j]][i]*(1-p[[j]][i])*2*sample_sizes[j] / (2*sample_sizes[j]-1)
    no.het.per.allele[[j]][i] <- length(genind.by.locus[[j]][,i][!is.na(genind.by.locus[[j]][,i]) & genind.by.locus[[j]][,i] == 0.5])
    Hobs[[j]][i] <- no.het.per.allele[[j]][i] / num.allele[[j]][i]
    d[[j]][i] <- (Hobs[[j]][i]-Hexp[[j]][i]) / Hexp[[j]][i]
    w[j] <- sqrt(sample_sizes[j])*(nancycats@loc.nall[j]-1) / nancycats@loc.nall[j]
    d.per.locus[j] <- 1/sample_sizes[j]*sum(d[[j]])
    dw[j] <- d.per.locus[j]*w[j]
  }
}

D <- sum(dw)/sum(w)
Nb <- 1/2*D + 1/(2*(D+1))


