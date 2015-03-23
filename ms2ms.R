#' @title ms2ms
#' @description Simulates microsatellite data by implementing the code of Pidugu and Schlotterer (2006) in R.
#' @description Depends on a input simulated in ms. 
#' 
#' @author Rebecca Harris <rbharris@@uw.edu>
#' 
#' @return A genind object of simulated microsatellite data.
#' 
#' @param ms.output Output from ms (library phyclust)
#' @examples library(phyclust)
#' @examples ms.out <- ms(nsam = 20, nreps = 10, args = "-t 3 -I 2 10 10 -ej 0.2 2 1")
#' @examples ms.out <- ms(nsam = 20, nreps = 10, opts = "-t 3 -I 2 10 10 -ej 0.2 2 1")
#' @examples ms2ms(ms.out)

#ms.output<-sim_ms
ms2ms <- function(ms.output){
	cmd <- strsplit(ms.output[1], " ")[[1]]
	nsam <- as.numeric(cmd[2])
	nreps <- as.numeric(cmd[3])
	find.I <- which(cmd == "-I")
	# Did you specify more than one population?
	
  if (length(find.I) != 0 && cmd[find.I + 1] > 1) {
		n.pops <- as.numeric(cmd[find.I + 1])
		n.ind.p <- as.numeric(cmd[(find.I + 2):(find.I + n.pops + 1)])
		if (!any(n.ind.p %% 2 == 0)) {
			print("The number of individuals in each population must be even.")
			return(NULL)
		} else {
			ind.nam <- NULL
			for (np in 1:n.pops){
				ind.nam <- c(ind.nam, paste("Ind", c(1:(n.ind.p[np]/2)), ".p", np, sep = ""))
			}
		}
	} else {
		ind.nam <- paste("Ind", c(1:(nsam/2)), ".p1", sep = "")
	}
	pos <- grep("positions", ms.output)
	full.mat <- NULL
	for (i in 1: nreps){
		loci <- ms.output[(pos[i]+1):c(pos[i]+ nsam)]
		n.pos <- nchar(loci[i])
		mut <- sample(c("+", "-"), n.pos, replace = TRUE)
		val.vec <- rep(0, nsam)
		for (j in 1: nsam){
			locus <- strsplit(loci[j], "")[[1]]
			val <- 100
			for (k in 1:n.pos){
				if (mut[k] == "+" & locus[k] == 1) val <- val + 1
				if (mut[k] == "-" & locus[k] == 1) val <- val - 1
			}
			val.vec[j] <- val
		}	
		alleles <- unique(val.vec)
		al.mat <- matrix(0, nrow = nsam/2, ncol = length(alleles))
		colnames(al.mat) <- paste("L", i, ".", c(1:length(alleles)), sep = "")
		for (ii in 1:(nsam/2)){
			if (ii == 1) jj <- ii
			else jj <- jj + 2
			al.mat[ii, which(alleles == val.vec[jj])] <- sum(0.5, al.mat[ii, which(alleles == val.vec[jj])])
			al.mat[ii, which(alleles == val.vec[jj+1])] <- sum(0.5, al.mat[ii, which(alleles == val.vec[jj+1])])
		}
		full.mat <- cbind(full.mat, al.mat)
		rownames(full.mat) <- ind.nam
	}
	full.mat <- genind(tab = full.mat, pop = as.numeric(unlist(strsplit(ind.nam, ".p"))[seq(2, nsam, 2)]))
	return(full.mat)
}
