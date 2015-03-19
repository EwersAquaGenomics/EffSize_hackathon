Simulation of Genealogies with the Coalescent
========================================================

Ape, phyclust (ms) and phylodyn have functions that simulates genealogies under the coalescent model for specific demographic scenarios.

A Genealogy (tree) consist of Topology and coalescent times. 

* Topology can be *isochronous* (all tips are "sampled" at time 0) or *heterochronous* (tips have different "sampling times"). Lineages merge (coalesce) at random between the existant lineages at the coalescent times. For isochronous coalescent we can generate a topology independently of the coalescent times but that is not the case for heterochronous topologies.   

* Coalescent times are exponentially distributed random variables with a rate that depends on the number of lineages and the population size trajectory. Again, for isochronous we can generate those times independently of the topology but not for the heterochronous case. 


The demographic models are:

1. Constant Population Size 

2. Exponential Growth 

3. Step Model with one change point

4. Exponential double growth model 

5. Any R function for defining $$N_{e}$$ 

We will use the function simulate.tree to generate a simulation from any type of topology (isochronous or heterochronous) and from any demographic model.


Examples
==================================
The following example generate 2 isochronous trees from a constant population size (Ne=1) with 10 tips.

```r
source("simulate_trees.R")
library(ape)

my.out1<-simulate.tree(n=10,N=2)
my.out1
```

```
## $out
## 2 phylogenetic trees
## 
## $description
## [1] "Simulation from a constant population size Ne=1"
```

```r
par(mfrow=c(1,2))
plot(my.out1$out[[1]])
axisPhylo()
plot(my.out1$out[[2]])
axisPhylo()
mtext("Standard Constant",side=3,line=-1.5,outer=TRUE)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png) 

The following example uses the popular "ms" program for simulation 2 isochronous trees from an exponential growth model

```r
my.out2<-simulate.tree(n=10,N=2,simulator="ms",args="-T -G 0.1")
my.out2
```

```
## $out
## 2 phylogenetic trees
## 
## $description
## [1] "Simulation using ms with args: -T -G 0.1"
```

```r
par(mfrow=c(1,2))
plot(my.out2$out[[1]])
axisPhylo()
plot(my.out2$out[[2]])
axisPhylo()
mtext("Exponential Growth with ms",side=3,line=-1.5,outer=TRUE)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png) 

Here, we consider a more general demographic model such as bottleneck. We specify our trajectory through a function.

```r
bottleneck_traj<-function(t){
  result=rep(0,length(t))
  result[t<=0.5]<-1
  result[t>0.5 & t<1]<-.1
  result[t>=1]<-1
  return(result)
}

my.out3<-simulate.tree(n=10,N=4,simulator="thinning",Ne=bottleneck_traj,max=11)
my.out3
```

```
## $out
## 4 phylogenetic trees
## 
## $description
## [1] "isochronous simulation using thinning with trajectory provided"
```

```r
par(mfrow=c(2,2))
plot(my.out3$out[[1]])
axisPhylo()
plot(my.out3$out[[2]])
axisPhylo()
plot(my.out3$out[[3]])
axisPhylo()
plot(my.out3$out[[4]])
axisPhylo()
mtext("Isochronous Bottleneck",side=3,line=-1.5,outer=TRUE)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png) 
And finally, we simulate heterochronous genealogies from a bottleneck. Here, we provide the sampling times

```r
set.seed(123)
samp_times = c(0, sort(runif(40, 0, .8)))
n_sampled = c(10, rep(1, 40))
sample<-cbind(n_sampled,samp_times)
my.out4<-simulate.tree(n=10,N=2,simulator="thinning",Ne=bottleneck_traj,max=10,sampling="hetero",sample=sample)
```

```
## Error in coalgen_thinning_hetero(sample, fun_inv, max): could not find function "generate_newick"
```

```r
my.out4
```

```
## Error in eval(expr, envir, enclos): object 'my.out4' not found
```

```r
par(mfrow=c(1,2))
plot(my.out4$out[[1]],show.tip.label=FALSE)
```

```
## Error in plot(my.out4$out[[1]], show.tip.label = FALSE): object 'my.out4' not found
```

```r
axisPhylo()
```

```
## Error in axis(side = side, at = x, labels = lab, ...): plot.new has not been called yet
```

```r
plot(my.out4$out[[2]],show.tip.label=FALSE)
```

```
## Error in plot(my.out4$out[[2]], show.tip.label = FALSE): object 'my.out4' not found
```

```r
axisPhylo()
```

```
## Error in axis(side = side, at = x, labels = lab, ...): plot.new has not been called yet
```

```r
mtext("Heterochronous Bottleneck",side=3,line=-1.5,outer=TRUE)
```

```
## Error in mtext("Heterochronous Bottleneck", side = 3, line = -1.5, outer = TRUE): plot.new has not been called yet
```



