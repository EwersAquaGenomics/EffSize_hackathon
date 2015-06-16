### METHOD 3: Temporal method ###

require(adegenet)

#' @title varNe_point
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
genObj <- simG3
t <- c(0,1,2) # provided by the user, the generations at which the samples were taken

varNe_point <- function(genObj, t = seq(1:length(genObj@pop.names))) {
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
    loc.pop[[i]] <- genind2genpop(genind.by.loc[[i]]) # one list element per locus
    loc.pop.s[,i] <- rowSums(loc.pop[[i]]@tab)/2 # number of ind. sampled per generation. Matrix with loci in columns, generations in rows
  }

  dif <- ff <- z <- matrix(0, nrow=as.numeric(paste(G,G-1,sep="")), ncol=ncol(f))
  s <- matrix(0, nrow=as.numeric(paste(G,G-1,sep="")), ncol=ncol(loc.pop.s))
  counter <- counter.j <- counter2 <- time.gen <- 0

  for (i in 1:(G-1)) {
    j <- i
    repeat {  
      j <- j+1
      dif[as.numeric(paste(i,j,sep="")),] <- sqrt((f[i,]-f[j,])^2) # allele frequency differences between generations
      z[as.numeric(paste(i,j,sep="")),] <- (f[i,]+f[j,])/2 * (1-(f[i,]+f[j,])/2) # average allele frequency?
      #ff[as.numeric(paste(i,j,sep="")),] <- (f[i,] != 0) & f[j,] != 0
      s[as.numeric(paste(i,j,sep="")),] <- 1/(0.5*(1/loc.pop.s[i,]+1/loc.pop.s[j,]))
      time.gen[as.numeric(paste(i,j,sep=""))] <- t[j]-t[i]
      counter[paste(i,j,sep="")] <- paste(i,j,sep="")
      counter2[paste(i,j,sep="")] <- paste(i,j,sep="-")
      counter.j[paste(i,j,sep="")] <- j
      if (j == G) break}
  }

  # remove empty columns which do not contain comparisons
  counter <- as.numeric(counter[-1])
  counter2 <- counter2[-1]
  counter.j <- counter.j[-1]
  dif <- dif[counter,]
  #ff <- ff[counter,]
  z <- z[counter,]
  s <- s[counter,] # harmonic mean sample sizes: columns are loci, rows are comparisons
  time.gen <- time.gen[counter]

  # there are some NA in z. I do not know where they come from (likely missing data), but they need to be replaced
  z[is.na(z)] <- 0

  if (no.com == 1) {
    for (j in 1:length(dif)) { 
      if (z[j] == 0) (dif[j] <- NA)
    }
    numerator <- sum(dif, na.rm=T) # sum allele frequency differences per comparison between generations (rows)
    denominator <- sum(z, na.rm=T)
    Fs <- numerator/denominator # single Fs value
    
    S <- 1/(1/sum.K*sum(K/s, na.rm=T)) # harmonic mean sample size per comparison between two generations
    s.per.gen <- 0

    for (i in 1:G) {
      s.per.gen[i] <- 1 / (1/sum.K * sum(K/loc.pop.s[i,], na.rm=T)) # mean of sample size per generation, weigthed by number of alleles per locus   
    }
  }
  
  if (no.com > 1 & L > 1) {
    for (i in 1:nrow(dif)) {
      for (j in 1:ncol(dif)) {
        if (z[i,j] == 0) (dif[i,j] <- NA)
      }
    }
    numerator <- rowSums(dif, na.rm=T) # sum allele frequency differences per comparison between generations (rows)
    denominator <- rowSums(z, na.rm=T)
    Fs <- numerator/denominator # one Fs value per generation comparison. Normally 2-10 values.
    
    # generate weighted harmonic means of sample sizes
    S <- s.per.gen <- 0

    for (i in 1:no.com) {
      S[i] <- 1/(1/sum.K*sum(K/s[i,], na.rm=T)) # harmonic mean sample size per comparison between two generations
      }

    for (i in 1:G) {
      s.per.gen[i] <- 1/(1/sum.K * sum(K/loc.pop.s[i,], na.rm=T)) # mean of sample size per generation, weigthed by number of alleles per locus   
    }
  }
    
  if (no.com > 1 & L == 1) {
    for (i in 1:nrow(dif)) {
      for (j in 1:ncol(dif)) {
        if (z[i,j] == 0) (dif[i,j] <- NA)
      }
    }
    numerator <- rowSums(dif, na.rm=T) # sum allele frequency differences per comparison between generations (rows)
    denominator <- rowSums(z, na.rm=T)
    Fs <- numerator/denominator # one Fs value per generation comparison. Normally 2-10 values.
    
    # generate weighted harmonic means of sample sizes
    S <- s.per.gen <- 0
    
    for (i in 1:no.com) {
      S[i] <- 1/(1/sum.K*sum(K/s[i], na.rm=T)) # harmonic mean sample size per comparison between two generations
    }
    
    for (i in 1:G) {
      s.per.gen[i] <- 1/(1/sum.K * sum(K/loc.pop.s[i,], na.rm=T)) # mean of sample size per generation, weigthed by number of alleles per locus   
    }
  }
  
  
  # calculate F "slash", one value per comparison between generations
  Fslash <- 0
  for (i in 1:no.com) {
    j <- as.numeric(counter.j[i])
    Fslash[i] <- (1/((1+Fs[i]/4)*(1-1/(2*s.per.gen[j])))) * (Fs[i]*(1-(1/(4*S[i])))-1/S[i])
  }
  #Fslash

  Ne <- data.frame(between.generations=counter2, Ne=(time.gen/(2*Fslash)))
  #Ne
  
  return(Ne)
}

# jackknifing function

#' @title jackknifing samples
#' @description jackknifes samples
#'  
#' @param bs jackknifed populations (list of genind objects)
#' @param Ne_point output of function varNe_point. Matrix containing point estimates of Ne, one estimate per comparison
#' @param statistic function to be employed, here varNe_point
#' 
#' @return a list. Each element is a jackknifed genind object
#'
#' @author Christine Ewers-Saucedo <ewers.christine@@gmail.com>
#' 
sum_jack2 <- function (Ne_point, bs, statistic) {
  nreps <- length(bs)
  stats <- lapply(bs, statistic)
  ncom <- nrow(stats[[1]])
  
  s <- vector("list", length=ncom)
  for (i in 1:ncom) {
    for (j in 1:nreps) {
      s[[i]][j] <- stats[[j]][i,2]
    }
  }
  
  o <- matrix(0, nrow=ncom, ncol=4)
  for (i in 1:ncom) {
    o[i,] <- c(median(s[[i]]), quantile(s[[i]], c(0.025, 0.975), na.rm = TRUE), var(s[[i]]))
  }
  
  oo <- data.frame(Ne_point, o)
  colnames(oo) <- c("between.generations","point.estimate", "median", "0.025", "0.975", "variance")
  rownames(oo) <- rownames(stats[[1]])
  
  return(oo)
}

# UNDER DEVELOPMENT ###
jacknife_populations <- function(genObj) {
  nsam <- as.vector(table(genObj@pop))
  for (i in 1:length(genObj@pop.names)) {
    sub.pop <- genObj[genObj@pop == genObj@pop.names[i]]
    x <- vector("list", length=nsam)
    for (i in 1:nsam) {
      x[[i]] <- sub.pop[-i]
    }
  }
  return(x) 
}
####


#' @title jackknifing samples
#' @description jackknifes samples
#'  
#' @param genObj A genind object with years as populations.
#' 
#' @return a list. Each element is a jackknifed genind object
#'
#' @author Christine Ewers-Saucedo <ewers.christine@@gmail.com>
#' 
jack_samples <- function(genObj) {
  nsam <- nrow(genObj@tab)
  x <- vector("list", length=nsam)
  for (i in 1:nsam) {
    x[[i]] <- genObj[-i]
  }
  return(x)
}


# combine point estimate and jackknifing summaries

#' @title varNe
#' @description Calculates effective population size from a genind object by calculating the variance of allele frequencies between generations, and jackknifes the samples. 
#' @description The variance F is inversely proportional to the effective population size Ne.
#' @description This function calculates F based on Jorde and Ryman 2007. They state hat low frequency alleles are not problematic / do not bias Ne estimates, thus they are not removed here.
#'  
#' @param genObj A genind object with years as populations.
#' 
#' @return a data.frame. Each column reports the result for a comparison of two generations, 
#' @return and the rows correspond to the generations compared, the point estimate of variance Ne, and the summary statistics of the jackknifing (median, 0.025 quantile, 0.975 quantile, variance)
#'
#' @author Christine Ewers-Saucedo <ewers.christine@@gmail.com>
#' 
#' @examples data(simG3)
#' @examples varNe(simG3)
#' 
varNe <- function(genObj) {
  Ne_point <- varNe_point(genObj)
  jn <- jack_samples(genObj)
  Ne_jk <- sum_jack2(Ne_point, jn, varNe_point)
  
  return(Ne_jk)
}
