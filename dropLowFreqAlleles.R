### remove low frequency alleles ###

require(adegenet)
data(nancycats)

#genind.by.locus <- seploc(nancycats, res.type="matrix", truenames=T)
genind.by.locus <- seploc(nancycats, truenames=T)

K <- length(nancycats@all.names)
crit <- c(0, 0.01, 0.02, 0.05)

### calculating numbers and frequencies of alleles per locus

p <- num.allele <- vector("list", length=K)
num.alleles.per.locus <- 0

#j=1
#i=1
for (j in 1:K) {
  for (i in 1:ncol(genind.by.locus[[j]]@tab)) {
    num.alleles.per.locus[j] <- sum(genind.by.locus[[j]]@tab, na.rm=T) # total number of genotyped alleles for locus 1, not the number of unique alleles
    num.allele[[j]][i] <- sum(genind.by.locus[[j]]@tab[,i], na.rm=T) # number of occurrences of allele 1 of locus 1
    p[[j]][i] <- num.allele[[j]][i] / num.alleles.per.locus[j] # allele frequency
  }
}


### removing alleles below critical frequency (provided by vector "crit")
# storing filtered genind objects in list "genind.crit"

genind.crit <- vector("list", length=length(crit))
pp <- vector("list", length=length(crit))

for (i in 1:length(crit)) {
  for (j in 1:K) {
    pp <- p[[j]] >= crit[i]
    genind.crit[[i]] <- genind.by.locus
    genind.crit[[i]][[j]]@tab <- genind.by.locus[[j]]@tab[,pp]
  }
}

# genind.crit is a list of lists of genind objects. The top level list separates the genind objects by critical values
# the second level contains genind objects, one for each locus