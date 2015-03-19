#' @title .George
#' @description Thibaut's desire
#' 
#' @examples .George(1)

.George <- function(x)
{
  print("This function doesn't do anything")
}

#' @title coalescent.intervals.datedPhylo
#' @description Functions to find the skyline of a binary, rooted tree with heterochronous
#' @description tips, i.e. sequences were sampled sequentially in time.
#' @description The likelihood function is based on Equations 3 and 4 in
#' @description Drummond et al (2005) Mol. Biol. Evol.
#' 
#' @author Lucy Mengqi Li <mengqi.li09@@imperial.ac.uk>
#' @examples library(ape)
#' @examples trees <- rmtree(N=5,n=20)
#' @examples coalescent.intervals.datedPhylo(trees)

coalescent.intervals.datedPhylo <- function (tr) {
  #
  # Find the coalescent intervals of a tree in units of substitutions
  #
  # Input:
  #   tr - a phylogeny of class 'datedPhylo'
  #
  # Returns a datedCI object with coalescent intervals
  #
#   if (!inherits(tr, "datedPhylo")) 
#     stop("object \"tr\" is not of class \"datedPhylo\"")
  n <- length(tr$tip.label)
  e1 <- tr$edge[, 1]
  e2 <- tr$edge[, 2]
  EL <- tr$edge.length
  i=1
  depths <- lapply(1:n, function (i) { #for each tip, create a matrix with the indexes of the nodes that are travelled on the way to the root, with the branch length of the node that ends there (the first is zero because the tip has zero length)
    start <- match(i, e2) #find the row number of the edge that this tip starts (going backwards)
    d <- c(0, EL[start]) #the edge.length leading from this tip
    ends <- c() #vector to store nodes that end branches
    while(!is.na(start)) { #provided you're still in the tree...
      end <- e1[start] #the node that ends this branch (going backwards)
      ends <- append(ends, end) #add this end to the vector of ends
      start <- match(end, e2) #start then becomes the row number of the edge
      if (!is.na(start)) d <- append(d, d[length(d)] + EL[start]) #provided you're still in the tree, add the branch length to the vector d
    }
    cbind(d, c(i, ends))
  })
  
  positions <- do.call(rbind, depths)[, 1]
  max.depth <- max(positions) #longest path from root to tip in the tree
  adjusted.times <- unlist(lapply(depths, function (x) { 
    max.depth - max(x[, 1]) + x[, 1]
  }))
  is.coalescent <- unlist(lapply(depths, function (x) {
    c(FALSE, rep(TRUE, nrow(x)-1))
    }))
  nodes <- do.call(rbind, depths)[, 2]
  times <- adjusted.times[!duplicated(nodes)]
  coal <- is.coalescent[!duplicated(nodes)]
  ord <- order(times)
  ltt <- cumsum(-coal[ord]+!coal[ord])
  obj <- list(lineages=ltt[-length(ltt)],
              interval.length=diff(times[ord]),
              is.coalescent=coal[ord][-1],
              interval.count=tr$Nnode + n - 1,
              coalescent.count=tr$Nnode,
              total.depth=max.depth)
  class(obj) <- "datedCI"
  obj
}


#' @title skyline.with.sampling
#' @description For a given datedCI object return the skyline for each coalescent interval
#' 
#' @author Lucy Mengqi Li <mengqi.li09@@imperial.ac.uk>
#' @examples library(ape)
#' @examples trees <- rmtree(N=5,n=20)
#' @examples ci <- coalescent.intervals.datedPhylo(trees)
#' @examples skyline.with.sampling(ci)

skyline.with.sampling <- function (ci) {
  #
  # For a given datedCI object return the skyline for each coalescent interval
  #
  # Inputs:
  #   ci - datedCI object
  # Returns an object of class 'datedSkyline'
  #
  coal.which <- which(ci$is.coalescent)
  diffs <- diff(c(0, coal.which))
  coal.interval.lengths <- ci$interval.length
  skylines <- sapply(1:ci$coalescent.count, function (i) {
    sub <- rev(seq(coal.which[i], by=-1, length.out=diffs[i]))
    ltt <- ci$lineages[sub]
    interval.lengths <- ci$interval.length[sub]
    sk <- sum(choose(ltt, 2) * interval.lengths)
    L <- log(choose(ltt[length(ltt)], 2)/sk) - 
      sum(choose(ltt, 2)/sk*interval.lengths) # likelihood
    c(sk, L, sum(interval.lengths))
  })
  logL <- sum(skylines[2, ])
  obj <- list(time=cumsum(skylines[3, ]),
              interval.length=skylines[3, ],
              population.size=skylines[1, ],
              parameter.count=ci$coalescent.count - 1,
              logL=logL
              )
  class(obj) <- "skyline"
  obj
}

#' @title skyline.datedPhylo
#' @description Like skyline() from ape but takes the class "datedPhylo"
#' 
#' @author Lucy Mengqi Li <mengqi.li09@@imperial.ac.uk>
#' @examples library(ape)
#' @examples trees <- rmtree(N=5,n=20)
#' @examples skyline.datedPhylo(trees)

# ----- Main functions ------#
skyline.datedPhylo <- function (tr) {
  # tr should be of class 'datedPhylo'
  require(ape)
  ci <- coalescent.intervals(tr)
  skyline.with.sampling(ci)
}

# #' @title Phylos2Skylines
# #' @description Converts a (multi)Phylo object to a series of skylines
# #' @description For a set of unrooted trees, remove the burnin trees and root the rest at 
# #' @description the root.node. Limit the total number of output trees to be less than max.trees
# #' @description Requires input from posterior sample of trees and parameters e.g. BEAST
# #' 
# #' @author Lucy Mengqi Li <mengqi.li09@@imperial.ac.uk>
# #' @examples library(ape)
# #' @examples trees <- rmtree(N=5,n=20)
# #' @examples Phylos2Skylines(trees)
# #' 
# 
# Phylos2Skylines <- function (trees.file.name, nex=TRUE, root.node=NULL,
#                              param.file.name, burninfrac=0.5, max.trees=1000) {
#   #
#   # For a set of unrooted trees, remove the burnin trees and root the rest at 
#   # the root.node. Limit the total number of output trees to be less than max.trees
#   #
#   if (nex)  trs <- read.nexus(trees.file.name)
#   else trs <- read.tree(trees.file.name)
#   
#   index.range <- c(floor(burninfrac*length(trs)+1), length(trs))
#   if (diff(index.range) >= max.trees) {
#     index.seq <- seq(floor(burninfrac*length(trs)+1), length(trs),
#                      length.out=max.trees)
#   } else {
#     index.seq <- floor(burninfrac*length(trs)+1):length(trs)
#   }
#   if(is.null(root.node)) {
#     rooted.trs <- trs
#   } else {
#     rooted.trs <- lapply(index.seq, function (i){
#       rooted <- root(trs[[i]], root.node, resolve.root=TRUE)  # MrBayes trees are already rooted
#       drop.tip(rooted, root.node, FALSE)
#       #drop.tip(trs[[i]], root.node, FALSE)
#     })
#   }
#   
#   pars <- read.table(param.file.name, header=TRUE, skip=1)
#   clock.rates <- pars$Clockrate[index.seq]
#   lapply(1:length(rooted.trs), function (i) {
#     x <- rooted.trs[[i]]
#     class(x) <- "datedPhylo"
#     sk <- skyline(x)
#     list(Skyline=sk, MolClockRate=clock.rates[i])
#   })
# }