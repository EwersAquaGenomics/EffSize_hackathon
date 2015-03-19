### METHOD 1: Linkage disequilibrium Ne ###

# problems: we dont have a formula to calculate r hat, only its expectation
# missing: implementation that accounts for missing data
# missing: confidence intervals

require(adegenet)
harmonic<-function(vec)  1/mean(1/vec)

?adegenet
ls(package:adegenet)

data(nancycats)
nancycats
head(nancycats@tab)

nancycats@ind.names
K<-length(nancycats@all.names)
nancycats@loc.nall
nancycats@loc.fac

# sample_sizes<-lapply(nancycats@all.names,function(x) sum(as.numeric(x)))

setGeneric(name="getSampleSize_byLocus",
#   signature="genind",
   def=function(theObject)
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
)


setGeneric(name="getHeterozygoteNumbers",
   #   signature="genind",
   def=function(theObject)
   {
     K<-length(theObject@all.names)
     
     ends<-cumsum(theObject@loc.nall)
     starts<-c(0,ends[-length(ends)])+1
     
     heterozygote_numbers<-list()
     i=1
     for (i in 1:K){
       heterozygote_numbers[[i]]<-colSums(na.rm=T,theObject@tab[,starts[i]:ends[i]]==0.5)/colSums(!is.na(theObject@tab[,starts[i]:ends[i]]))
     }
     return(heterozygote_numbers)
   }
)

sample_sizes <- getSampleSize_byLocus(nancycats)



s<-num.alleles<-Expect_r2<-w<-matrix(NA,nrow=K,ncol=K)
N=0
NonS=0
for(i in 1:(K-1)){
  for(j in (i+1):K){
    s[i,j] <- harmonic(c(sample_sizes[i],sample_sizes[j]))
    num.alleles[i,j] <- (nancycats@loc.nall[i]-1)*(nancycats@loc.nall[j]-1)
 
    if(s[i,j]>=30)
    {
      Expect_r2[i,j] <- 1/s[i,j] + 3.19/s[i,j]^2
    } else {
      Expect_r2[i,j] <- 0.0018 + 0.907/s[i,j] + 4.44/s[i,j]^2
    }      
  }
}

w=num.alleles*s^2

N=sum(num.alleles,na.rm=T)
NonS=sum(num.alleles/s,na.rm=T)
S=1/(NonS)*N
W=sum(w,na.rm=T)

R2_hat=sum(w*Expect_r2,na.rm=T)/W



if(S>=30)
{
  Expect_R2 <- 1/S + 3.19/S^2
  R2_drift <- R2_hat - Expect_R2
  Ne <- 1/(2*R2_drift)*(1/3 + sqrt(1/9 - 2.76*R2_drift))
} else {
  Expect_R2 <- 0.0018 + 0.907/S + 4.44/S^2
  R2_drift <- R2_hat - Expect_R2
  Ne <- 1/(2*R2_drift)*(1/3 + sqrt(0.94864 - 2.08*R2_drift))  
}      

NOT WORKING!!!
  
###### METHOD 2 ######

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



