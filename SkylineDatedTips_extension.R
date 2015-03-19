#' @title getNodeAges
#' @description Returns a vector with the depths, or heights of the nodes in the tree
#' @description Depends on a hidden function nodeHeight, from phangorn. 
#' @description Values by default are backwards from the youngest tip which is zero. They can be forward relative to the root if from_past=TRUE.
#' 
#' @author George Shirreff <georgeshirreff@@gmail.com>
#' 
#' @return For a tree with n tips, the vector is <tip1, tip2...tipn, root, internal_node2, internal_node3>
#' 
#' @examples trees <- rmtree(N=5,n=20)
#' @examples Phylos2Skylines_anchor(trees)

getNodeAges<-function(x,from_past=F)
{
  require(phangorn)
  nodeAges<-phangorn:::nodeHeight(x)
  names(nodeAges)<-c(x$tip.label,"root",rep(NA,length(nodeAges)-length(x$tip.label)-1))
  
  if(from_past)
  {
    return(max(nodeAges)-nodeAges)
  }
  return(nodeAges)
}

#' @title Phylos2Skylines_anchor
#' @description Converts a (multi)Phylo object to a series of skylines
#' @description For a set of unrooted trees, remove the burnin trees and root the rest at 
#' @description the root.node. Limit the total number of output trees to be less than max.trees
#' @description Robust to trees with tips sampled at different times
#' @description Can be directly applied to BEAST or MrBayes output, or files can be directly parsed by argument
#' 
#' @param trees
#' @param output_type
#' @param file_type
#' @param package
#' @param root_node
#' @param param.file.name
#' @param burninfrac
#' @param max.trees
#' @param youngestTip
#' @param timeForward
#' @param scaling
#' 
#' @author Lucy Mengqi Li <mengqi.li09@@imperial.ac.uk>, adapted by George Shirreff <georgeshirreff@@gmail.com>
#' @examples trees <- rmtree(N=5,n=20)
#' @examples Phylos2Skylines_anchor(trees)
#' 

Phylos2Skylines_anchor <- function (trees, output_type=c("list","matrices","master","conf.int","conf.int.plot") ,file_type="nex", package=NULL, root.node=NULL,
                             param.file.name=NULL, burninfrac=0, max.trees=1000,youngestTip=0,timeForwards=F,scaling=1) {

  require(ape)
  if(is.character(trees))
  {
    if (file_type=="nex")  
      trs <- read.nexus(trees)
    else if (file_type=="nwk")
      trs <- read.tree(trees)
  } else {
    trs<-trees
  }

  if(class(trs)=="phylo")
  {
    tmp<-list()
    tmp[[1]]=trs
    class(tmp)<-"multiPhylo"
    trs<-tmp
  }
  
  index.range <- c(floor(burninfrac*length(trs)+1), length(trs))
  if (diff(index.range) >= max.trees) {
    index.seq <- seq(floor(burninfrac*length(trs)+1), length(trs),
                     length.out=max.trees)
  } else {
    index.seq <- floor(burninfrac*length(trs)+1):length(trs)
  }
  if(is.null(root.node)) {
    rooted.trs <- trs
  } else {
    rooted.trs <- lapply(index.seq, function (i){
      rooted <- root(trs[[i]], root.node, resolve.root=TRUE)  # MrBayes trees are already rooted
      drop.tip(rooted, root.node, FALSE)
      #drop.tip(trs[[i]], root.node, FALSE)
    })
  }
  
  if(!is.null(param.file.name))
  {
    pars <- read.table(param.file.name, header=TRUE, skip=1)
    if(package=="BEAST")
    {
      clock.rates <- pars$clock.rate[index.seq]
    } else if (package=="MrBayes"){
      clock.rates <- pars$Clockrate[index.seq]  		
    } else {
      stop("parameter file package should be specified")
    }
  } else {
    clock.rates=rep(NA,length(rooted.trs))
  }

  #i=1
  skyline_list<-lapply(1:length(rooted.trs), function (i) {
    x <- rooted.trs[[i]]
    if(class(x)!="phylo") stop(paste("object",i,"is not a phylo object"))    
    if(!is.binary.tree(x)) stop(paste("tree",i,"is not binary, with unknown consequences"))    
    class(x) <- "datedPhylo"
    sk <- skyline(x)
    
    if(timeForwards)
    {
      sk$time = youngestTip + scaling*(sk$time - max(sk$time))
    } else {
      sk$time = scaling*(sk$time) + youngestTip
    }
    
    
    list(Skyline=sk, MolClockRate=clock.rates[i])
  })
  
  length(sk$time)
  length(sk$interval)
  
  if(output_type=="list")
  {
    return(skyline_list)
  } 
  
  time_mat<-sapply(1:length(skyline_list),function(x) skyline_list[[x]]$Skyline$time)
  pop_mat<-sapply(1:length(skyline_list),function(x) skyline_list[[x]]$Skyline$population.size)  
  
  if (output_type=="matrices") {
    return(list(time_mat=time_mat,pop_mat=pop_mat))
  }
  
  treenum_mat<-matrix(1:length(skyline_list),nrow=nrow(time_mat),ncol=length(skyline_list),byrow=T)
  posnum_mat<-matrix(1:nrow(time_mat),nrow=nrow(time_mat),ncol=length(skyline_list),byrow=F)
  
  treenum_vec<-treenum_mat[order(time_mat)]
  posnum_vec<-posnum_mat[order(time_mat)]
  
  master_time_vec=sort(time_mat)
  master_time_vec=c(master_time_vec[1]-1e-300,master_time_vec)
  master_pop_mat<-matrix(NA,nrow=length(master_time_vec),ncol=length(skyline_list))
  master_pop_mat[1,]<-0
  
  s=1
  for(s in 1:length(time_mat))
  {
    master_pop_mat[s+1,]<-master_pop_mat[s,]
    master_pop_mat[s+1,treenum_vec[s]]<-pop_mat[posnum_vec[s],treenum_vec[s]]    
  }
  
  rownames(master_pop_mat)<-master_time_vec
  if(output_type=="master")
  {
    return(master_pop_mat)
  }
  
  conf.int_mat<-t(apply(master_pop_mat,1, function(x) quantile(x,probs=c(0.5,0.025,0.975))))
  
  if(output_type=="conf.int")
  {
    return(conf.int_mat)
  }
  
  
  if(output_type=="conf.int.plot")
  {
    conf.int.skyline(conf.int_mat,epsilon=epsilon,...)
    return(conf.int_mat)
  }
  
  
  stop("no output type specified")
}

#' @title conf.int.skyline
#' @description Makes a skyline plot from a matrix with medians and confidence intervals as produced by Phylos2Skylines_anchor
#' 
#' @author George Shirreff <georgeshirreff@@gmail.com>
#' 
#' @return A plot with median in black and confidence intervals in red
#' 
#' @examples trees <- rmtree(N=5,n=20)
#' @examples conf.int<-Phylos2Skylines_anchor(trees,output_type="conf.int")

conf.int.skyline<-function(conf.int,epsilon=0,...)
{
  end_times<-as.numeric(rownames(conf.int))
  plot(type="n",x=range(end_times),y=range(conf.int),ylab=expression(N[e]*tau),cex.lab=1,cex.axis=1,...)
  
  draw_skyline_smooth(end_times,pop_size=conf.int[,1],col=1,epsilon=epsilon,...)
  draw_skyline_smooth(end_times,pop_size=conf.int[,2],col=2,epsilon=epsilon,...)  
  draw_skyline_smooth(end_times,pop_size=conf.int[,3],col=2,epsilon=epsilon,...)       
}

#' @title draw_skyline
#' @description Draws a classic skyline onto a plot
#' 
#' @author George Shirreff <georgeshirreff@@gmail.com>
#' 
#' @return Draws the lines of a skyline
#' 
#' @examples trees <- rmtree(N=5,n=20)
#' @examples conf.int<-Phylos2Skylines_anchor(trees,output_type="conf.int")
#' @examples end_times<-as.numeric(rownames(conf.int))
#' @examples plot(type="n",x=range(end_times),y=range(conf.int),ylab=expression(N[e]*tau))
#' @examples draw_skyline(end_times,pop_size=conf.int[,1],col=1)

draw_skyline<-function(end_times,pop_size,...)
{
  len<-length(pop_size)
  if(length(end_times)!=len) stop("vectors are different lengths")
  
  arrows(x0=end_times[-len],x1=end_times[-1],y0=pop_size[-len],code=0,...)
  arrows(x0=end_times[-1],y0=pop_size[-len],y1=pop_size[-1],code=0,...)		
}

#' @title draw_skyline
#' @description Draws a generalized skyline onto a plot
#' 
#' @author George Shirreff <georgeshirreff@@gmail.com>
#' 
#' @return Draws the lines of a skyline, but with a minimum time interval specified by epsilon
#' 
#' @examples trees <- rmtree(N=5,n=20)
#' @examples conf.int<-Phylos2Skylines_anchor(trees,output_type="conf.int")
#' @examples end_times<-as.numeric(rownames(conf.int))
#' @examples plot(type="n",x=range(end_times),y=range(conf.int),ylab=expression(N[e]*tau))
#' @examples draw_skyline_smooth(end_times,pop_size=conf.int[,1],col=1,epsilon=0.1)

draw_skyline_smooth<-function(end_times,pop_size,epsilon=0,...)
{
  if(epsilon==0)
  {
    draw_skyline(end_times,pop_size,...)
    return()
  }
  harmonic<-function(vec)  1/mean(1/vec)
  
  len<-length(pop_size)
  
  if(length(end_times)!=len) stop("vectors are different lengths")
  
  end_times_out=numeric(0)
  pop_size_out=numeric(0)
  intervals=c(0,diff(end_times)) #interval sizes
  
  cbind(1:length(end_times), end_times,intervals,pop_size)

  i=len
  while(i>=1)
  {
    combine_vec=i #the start of the intervals to be combined, even if only a single one
    while(sum(intervals[combine_vec])<epsilon) #if the combined interval is less than epsilon
    #shouldn't need to include an i-limit here, because the combine_interval should always be greater than one (provided epsilon isn't bigger than the total time!!)      
    {
      if(i>1) 
      {
        i=i-1
        combine_vec=c(combine_vec,i) #add another interval to the combined interval
      } else {
        combine_vec=c(last_combine_vec,combine_vec) #if there aren't any left, add the previous run of intervals
        end_times_out<-end_times_out[-1] #if this is the case, you have to take out the last entry of each vector so a combined value can be put in 
        pop_size_out<-pop_size_out[-1]
      }
    }
    
    end_times_out<-c(end_times[max(combine_vec)],end_times_out)
    pop_size_out<-c(harmonic(pop_size[combine_vec]),pop_size_out)
    
    last_combine_vec=combine_vec #in case the last run of intervals don't add up to >epsilon and you need to combine with the previous run of intervals
  
    i=i-1
  }
  
#   end_times_out<-c(end_times_out,end_times[1])
#   pop_size_out<-c(pop_size_out,pop_size[1])

  len_out<-length(end_times_out)
  arrows(x0=end_times_out[-len_out],x1=end_times_out[-1],y0=pop_size_out[-len_out],code=0,...)
  arrows(x0=end_times_out[-1],y0=pop_size_out[-len_out],y1=pop_size_out[-1],code=0,...)  	
}




