### METHOD 4: Temporal method ###

require(adegenet)

#' @title varNe 
#' @description Calculates population size from a genind object by calculating the variance of allele frequencies between generations. 
#' @description The variance F is inversely proportional to the effective population size Ne.
#' @description This function calculates F based on Jorde and Ryman 2007. They state hat low frequency alleles are not problematic / do not bias Ne estimates, thus they are not removed here.
#'  
#' @param genObj A genind object with years as populations.
#' 
#' @return a vector with 
#'
#' @author Christine Ewers-Saucedo <ewers.christine@@gmail.com>
#' 
#' @examples library(adegenet)
#' @examples data(nancycats)
#' @examples varNe(nancycats)
#' 
varNe <- function(genObj) {
  t <- seq(1:length(genObj@pop.names)) # provided by the user, the generations at which the samples were taken
  G <- length(genObj@pop.names) # number of generations sampled
  no.com <-  G*(G-1)/2 # number of possible comparisons between two generations
  genPopObj <- genind2genpop(genObj)
  f <- makefreq(genPopObj) # columns are populations/generations. Allele freq per locus per pop are calculated

  K <- genObj@loc.nall # number of alleles per loci overall
  sum.K <- sum(K)
  L <- length(K) # number of loci
  genind.by.loc <- seploc(genObj) # one genind object per locus as list elements

  loc.pop <- vector("list", length=L)
  loc.pop.s <- matrix(0, nrow=G, ncol=L)
  for (i in 1:L) {
    loc.pop[[i]] <- genind2genpop(genind.by.loc[[i]])
    loc.pop.s[,i] <- rowSums(loc.pop[[i]]@tab)/2 # loci are columns, generations are rows
  }

  dif <- ff <- z <- matrix(0, nrow=as.numeric(paste(G,G-1,sep="")), ncol=ncol(f))
  s <- matrix(0, nrow=as.numeric(paste(G,G-1,sep="")), ncol=ncol(loc.pop.s))
  counter <- counter.j <- counter2 <- time.gen <- 0

  for (i in 1:(G-1)) {
    j <- i
    repeat {  
      j <- j+1
      dif[as.numeric(paste(i,j,sep="")),] <- (f[i,]-f[j,])^2
      z[as.numeric(paste(i,j,sep="")),] <- (f[i,]+f[j,])/2 * (1-(f[i,]+f[j,])/2)
      ff[as.numeric(paste(i,j,sep="")),] <- f[i,] != 0 & f[j,] != 0
      s[as.numeric(paste(i,j,sep="")),] <- 1/(0.5*(1/loc.pop.s[i,]+1/loc.pop.s[j,]))
      time.gen[as.numeric(paste(i,j,sep=""))] <- t[j]-t[i]
      counter[paste(i,j,sep="")] <- paste(i,j,sep="")
      counter2[paste(i,j,sep="")] <- paste(i,j,sep="-")
      counter.j[paste(i,j,sep="")] <- j
      if (j == G) break}
  }

# remove empty columns which do not contain comparisons
counter <- as.numeric(counter[2:length(counter)])
counter2 <- counter2[2:length(counter2)]
counter.j <- counter.j[-1]
dif <- dif[counter,]
ff <- ff[counter,]
z <- z[counter,]
s <- s[counter,] # harmonic mean sample sizes: columns are loci, rows are comparisons
time.gen <- time.gen[counter]

# there are some NA in ff2. I do not know where they come from (likely missing data), but they need to be replaced
z[is.na(z)] <- 0

for (i in 1:nrow(dif)) {
     for (j in 1:ncol(dif)) {
       if (z[i,j] == 0) (dif[i,j] <- NA)
     }
}

# sum allele frequency differences per comparison between generations (rows)
numerator <- rowSums(dif, na.rm=T)
denominator <- rowSums(z, na.rm=T)
Fs <- numerator/denominator # one Fs value per generation comparison. Normally 2-10 values.

# generate weighted harmonic means of sample sizes
S <- s.per.gen <- 0

for (i in 1:no.com) {
  S[i] <- 1/(1/sum.K*sum(K/s[i,], na.rm=T)) # harmonic mean sample size per comparison between two generations
}

for (i in 1:G) {
  s.per.gen[i] <- 1 / (1/sum.K * sum(K/loc.pop.s[i,], na.rm=T)) # mean of sample size per generation, weigthed by number of alleles per locus   
}

# calculate F "slash", one value per comparison between generations
Fslash <- 0
for (i in 1:no.com) {
  j <- counter.j[i]
  Fslash[i] <- (1/((1+Fs[i]/4)*(1-1/(2*s.per.gen[j])))) * (Fs[i]*(1-(1/(4*S[i])))-1/S[i])
}
Fslash

# Ne
Ne <- data.frame(gen=counter2, Ne=t/(2*Fslash))
            
return(Ne)
}
