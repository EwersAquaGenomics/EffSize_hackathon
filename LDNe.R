# Linkage disequilibrium Ne

# number of polymorphic loci
k <- length(nancycats@all.names)
# pair of loci
Er2 <- 
  
library(adegenet)
data(nancycats)
nancycats@tab
x <- summary(nancycats)
x$loc.nall[1]
table(nancycats@pop)
nancycats@tab

sample.size.per.locus <- 0
for (i in 2:(k-1)) {
sample.size.per.locus[1] <- sum(nancycats@tab[,1:nancycats@loc.nall[2]], na.rm=T) #217
sample.size.per.locus[i] <- sum(cumsum(nancycats@tab[,((nancycats@loc.nall[1]):nancycats@tab(nancycats@loc.nall[8])]:(nancycats@loc.nall[8]+nancycats@loc.nall[8+1]))], na.rm=T)
    }

#nancycats@loc.nall
#L1 L2 L3 L4 L5 L6 L7 L8 L9 
#16 11 10  9 12  8 12 12 18 

i from 2 to k
k <- 2
# starting point
sum(nancycats@tab[,1:nancycats@loc.nall[1]], na.rm=T) #217
nancycats@tab[,nancycats@loc.nall[1]]:sum(nancycats@loc.nall[1:k])]

colnames(nancycats@tab[,(nancycats@loc.nall[1]+1):nancycats@loc.nall[2]])

