### makeP and makeCrit: remove low frequency alleles ###

library(adegnet)

#' @title makeP 
#' @description Calculates the frequency of each allele per locus
#' 
#' @author Christine Ewers-Saucedo <ewers.christine@@gmail.com>, George Shirreff <georgeshirreff@@gmail.com>
#' 
#' @examples library(adegenet)
#' @examples data(nancycats)
#' @examples makeP(nancycats)
#' 
makeP <- function(genObj) {
  genind.by.locus <- seploc(genObj)
  K <- length(genind.by.locus) # number of loci
  sample.loc <- getSampleSize_byLocus(genObj)  # total number of genotyped alleles for locus 1, not the number of unique alleles
  
  p <- occ.allele <- vector("list", length=K)
  alleles.per.locus <- 0
  
  for (j in 1:K) {
    for (i in 1:ncol(genind.by.locus[[j]]@tab)) {
      occ.allele[[j]][i] <- sum(genind.by.locus[[j]]@tab[,i], na.rm=T) # number of occurrences of allele 1 of locus 1
      p[[j]][i] <- occ.allele[[j]][i] / sample.loc[j] # allele frequency
  }
}
return(p)
}


#' @title makeCrit
#' @description Removes alleles below critical frequency (provided by vector "crit"). The output is a list of lists of genind objects. 
#' @description The top level list separates the genind objects by critical values. The second level contains genind objects, one for each locus
#' 
#' @author Christine Ewers-Saucedo <ewers.christine@@gmail.com>
#' 
#' @param genObj genind object
#' @param crit vector with critical values. All alleles with frequencies below each critical value are removed
#' 
#' @return list of genind object, one object per critical value
#' 
#' @examples library(adegenet)
#' @examples data(nancycats)
#' @examples makeCrit(nancycats, crit=c(0.05, 0.01))
#' 
makeCrit <- function(genObj, crit) {
  
  genind.by.locus <- seploc(genObj)
  p <- makeP(genObj)
  K <- length(genind.by.locus) # number of loci
  
  genind.crit <- vector("list", length=length(crit))
  pp <- vector("list", length=length(crit))

  for (i in 1:length(crit)) {
    for (j in 1:K) {
      pp <- p[[j]] >= crit[i]   
      genind.crit[[i]] <- genind.by.locus
      genind.crit[[i]][[j]]@tab <- genind.by.locus[[j]]@tab[,pp]
    }
  }
  return(genind.crit)
}