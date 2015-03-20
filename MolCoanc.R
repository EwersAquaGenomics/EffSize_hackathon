## Molecular coancestry ##

library(adegenet)

data(nancycats)

SbyL<-getSampleSize_byLocus(nancycats)

n<-sum(SbyL)
L<-length(SbyL)
nP<-(n-1)*n/2


