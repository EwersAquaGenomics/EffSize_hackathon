#' @title sequence.wrapper
#' @description A wrapper function for the sequence simulator programs ms, seqgen, simSeq, and PhyloSim.
#' 
#' @author Rebecca Harris <rbharris@@uw.edu>
#' 
#' @return Writes out sequences in phylip format to your working directory.
#' 
#' @param n number of sequences
#' @param L length of sequence
#' @param method choice of sequence simulation
#' @param ms.args arguments passed to ms
#' @param tree tree of format phylo
#' @param outfile name of the file that writes out
#' @examples sequence.wrapper(n = 10, L = 10, method = "gen.seq.HKY", tree = phylo)


#do we want multiple trees or to run though

sequence.wrapper <- function(n = 10, L = 10, method = NULL, ms.args = NULL, tree = NULL, pi.n = NULL, kappa = "", outfile = "out", ancseq = NULL, rate.in = 1){
	if (method == "ms") {
		if (is.null(ms.args))
			print("Error: No ms arguments specified")
		else if (length(grep("-T", ms.args, fixed = TRUE)) == 1)
			print("Error: This is a sequence generator. Do not use -T.")
		else {
			simseq <- ms(nsam = n, nreps = 1, opts = ms.args)
			ms2nuc(simseq, file = outfile)
		}
	} 
	if (method == "gen.seq.HKY"){
		if (is.null(tree)) print("Error: No tree given.")
		pi.n <- eq.prob(pi.n, 4)
		simseq <- gen.seq.HKY(tree, pi.n, kappa = 1, L, anc.seq = ancseq, rate.scale = rate.in)
	}
	if (method == "gen.seq.SNP"){
		if (is.null(tree)) print("Error: No tree given.")
		pi.n <- eq.prob(pi.n, 2)
		simseq <- gen.seq.SNP(tree, pi.n, L, anc.seq = ancseq)
	}
	if (method == "simSeq"){
		if (is.null(tree)) print("Error: No tree given.")
		simseq <- simSeq(tree, l = L, rootseq = ancseq, bf = pi.n, rate = rate.in)
		phang2nuc(simseq, file = outfile)
	}
	if (method == "PhyloSim"){
		if  (is.null(tree)) print("Error: No tree given.")
		else {
			if (!is.null(ancseq))
				ancseq <- NucleotideSequence(string = ancseq, )
			else {
				ancseq <- NucleotideSequence(length = L, processes = list(list(JC69()))) 
				ancseq$states <- c("A", "G", "C", "T")
			}
			simseq <- PhyloSim(phylo = tree, root.seq = ancseq)
			Simulate(simseq)
			phylosim2nuc(simseq, file = outfile)
		}
	}
	return(read.dna(outfile))
}

#' @title ms2nuc
#' @description Converts output of ms2ms to nucleotides.
#' 
#' @author Rebecca Harris <rbharris@@uw.edu>
#' 
#' @return 
#' 
#' @param
#' @examples

ms2nuc <- function(ms.res, fileout){
	if (length(grep("positions", ms.res)) == 0)
		print("No segregating sites. Reexecute loop or increase theta.")
	else {
		ms.res <- ms.res[(grep("positions", ms.res)+1):length(ms.res)]
		nn <- length(ms.res)
		out <- cbind(paste("s", 1:nn, sep = ""), ms.res)
		write.table(out, fileout, quote = FALSE, row.names = FALSE, col.names = c(nn, nchar(ms.res[1])))
	}	
}

#' @title eq.prob
#' @description Checks that appropriate equilibrium probabilities are specified, otherwise specifies equal probabilities. 
#' 
#' @author Rebecca Harris <rbharris@@uw.edu>
#' 
#' @return 
#' 
#' @param
#' @examples

eq.prob <- function(pi = NULL, n.prob = NULL){
	if (is.null(pi)) {
		pi <- rep(1/n.prob, n.prob)
		print("Equilibrium probabilies assumed to be equal. Otherwise specify pi")
	} else if (length(pi) != n.prob){
		print("Pi wrong length. Making equibrium probabilites equal.")
		pi <- rep(1/n.prob, n.prob)
	}
	return(pi)
}

#' @title phang2nuc
#' @description Converts output of phangorn to nucleotides.
#' 
#' @author Rebecca Harris <rbharris@@uw.edu>
#' 
#' @return 
#' 
#' @param
#' @examples

phang2nuc <- function(sq = phang.seq, file = "test.out"){
	out <- NULL
	for (i in 1:length(sq)){
		if (any(sq[[i]] == 1))
			sq[[i]][which(sq[[i]] == 1)] <- "A"
		if (any(sq[[i]] == 2))
			sq[[i]][which(sq[[i]] == 2)] <- "C"
		if (any(sq[[i]] == 3))
			sq[[i]][which(sq[[i]] == 3)] <- "G"
		if (any(sq[[i]] == 4))
			sq[[i]][which(sq[[i]] == 4)] <- "T"
		seqline <- paste("s", i, "        ", paste(sq[[i]], collapse = ""), sep = "")
		out <- rbind(out, seqline)
	}
	out2 <- rbind(paste(length(sq), length(sq[[i]])), out)
	write(out2, file)
}

#' @title phylosim2nuc
#' @description Converts PhyloSim output to phylip format.
#' 
#' @author Rebecca Harris <rbharris@@uw.edu>
#' 
#' @return 
#' 
#' @param
#' @examples

phylosim2nuc <- function(align = sim, file = "test.out"){
	tip.seqs <- grep("s.*", rownames(align$alignment))
	out <- apply(align $alignment[tip.seqs,], 1, function(x) paste(x, collapse = ""))
	out <- c(paste(length(out), nchar(out)[1], collapse = " "), out)
	write.table(out, file, col.names = FALSE, quote = FALSE, sep = "  ")
}
