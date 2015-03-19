coalgen_thinning_iso<-function(sample, trajectory, upper = 25, ...) {
  s = sample[2]
  n <- sample[1]
  out <- rep(0, n - 1)
  time <- 0
  j <- n
  while (j > 1) {
    time_p <- time + rexp(1, upper * j * (j - 1) * 0.5)
    if (upper< trajectory(time_p, ...) ){
      upper<-trajectory(time_p, ...)*1.2
      time_p <- time + rexp(1, upper * j * (j - 1) * 0.5)
    }
    time<-time_p
    if (runif(1) <= trajectory(time, ...)/upper) {
      out[n - j + 1] <- time
      j <- j - 1
    }
  }
  return(list(intercoal_times = c(out[1], diff(out)), lineages = seq(n, 2, -1),upper=upper))
}

coalgen_thinning_hetero<-function (sample, trajectory, upper = 25, ...) {
  n_sampled<-sample[,1]
  s_times<-sample[,2]
  s = sample[1, 2]
  b <- sample[1, 1]
  n <- sum(sample[, 1]) - 1
  m <- n
  nsample <- nrow(sample)
  sample <- rbind(sample, c(0, 10 * max(sample, 2)))
  out <- rep(0, n)
  branches <- rep(0, n)
  i <- 1
  while (i < (nsample)) {
    if (b < 2) {
      b <- b + sample[i + 1, 1]
      s <- sample[i + 1, 2]
      i <- i + 1
    }
    E <- rexp(1, upper * b * (b - 1) * 0.5)
    if (upper< trajectory(E+s, ...) ){
      upper<-trajectory(E+s, ...)*1.2
      E <- rexp(1, upper * b * (b - 1) * 0.5)
    }
    
    if (runif(1) <= trajectory(E + s, ...)/upper) {
      if ((s + E) > sample[i + 1, 2]) {
        b <- b + sample[i + 1, 1]
        s <- sample[i + 1, 2]
        i <- i + 1
      }
      else {
        s <- s + E
        out[m - n + 1] <- s
        branches[m - n + 1] <- b
        n <- n - 1
        b <- b - 1
      }
    }
    else {
      s <- s + E
    }
  }
  while (b > 1) {
    E <- rexp(1, upper * b * (b - 1) * 0.5)
    if (runif(1) <= trajectory(E + s, ...)/upper) {
      s <- s + E
      out[m - n + 1] <- s
      branches[m - n + 1] <- b
      n <- n - 1
      b <- b - 1
    }
    else {
      s <- s + E
    }
  }
  intercoal_times = c(out[1], diff(out))
  lineages = branches
  coal_times<-cumsum(intercoal_times)
  args<-gen_INLA_args(coal_times,s_times,n_sampled)
  out<-generate_newick(args,cbind(n_sampled,s_times))$newick
  return(out)
}

gen_INLA_args<-function (coal_times, s_times, n_sampled) {
  n = length(coal_times) + 1
  data = matrix(0, nrow = n - 1, ncol = 2)
  data[, 1] = coal_times
  s_times = c(s_times, max(data[, 1]) + 1)
  data[1, 2] = sum(n_sampled[s_times <= data[1, 1]])
  tt = length(s_times[s_times <= data[1, 1]]) + 1
  for (j in 2:nrow(data)) {
    if (data[j, 1] < s_times[tt]) {
      data[j, 2] = data[j - 1, 2] - 1
    }
    else {
      data[j, 2] = data[j - 1, 2] - 1 + sum(n_sampled[s_times > 
                                                        data[j - 1, 1] & s_times <= data[j, 1]])
      tt = length(s_times[s_times <= data[j, 1]]) + 1
    }
  }
  s = unique(sort(c(data[, 1], s_times[1:length(s_times) - 1])))
  event1 = sort(c(data[, 1], s_times[1:length(s_times) - 1]), index.return = TRUE)$ix
  n = nrow(data) + 1
  l = length(s)
  event = rep(0, l)
  event[event1 < n] = 1
  y = diff(s)
  coal.factor = rep(0, l - 1)
  t = rep(0, l - 1)
  indicator = cumsum(n_sampled[s_times < data[1, 1]])
  indicator = c(indicator, indicator[length(indicator)] - 1)
  ini = length(indicator) + 1
  for (k in ini:(l - 1)) {
    j = data[data[, 1] < s[k + 1] & data[, 1] >= s[k], 2]
    if (length(j) == 0) {
      indicator[k] = indicator[k - 1] + sum(n_sampled[s_times < 
                                                        s[k + 1] & s_times >= s[k]])
    }
    if (length(j) > 0) {
      indicator[k] = j - 1 + sum(n_sampled[s_times < s[k +1] & s_times >= s[k]])
    }
  }
  coal_factor = indicator * (indicator - 1)/2
  return(list(coal_factor = coal_factor, s = s, event = event, indicator = indicator))
}

##This is the main function

simulate.tree<-function(n=10,N=1,sampling="iso",args="-T -G 0.1",Ne=1,max=1,simulator=NULL,sample=NULL, ...){
  if (is.null(simulator)) {
    #Generates samples using rcoal with Ne=1. If this is a mistake, specify your simulator (ms,thinning,standard)
    out<-replicate(N,rcoal(n),simplify=FALSE)
    class(out)<-"multiPhylo"
    return(list(out=out,description="Simulation from a constant population size Ne=1"))
  }
  if (simulator=="ms"){
    library("phyclust")
    if (is.null(args)){args<-"-T"}
    out<-read.tree(text=paste(rep(ms(nsam = n, opts = args)[3],N),sep="\n"))  
    return(list(out=out,description=paste("Simulation using ms with args:",args,sep=" ")))
  }
  if (simulator=="thinning"){
    if (is.function(Ne)){
      ##A function is needed for thinning
      fun_inv<-function(t,...){
        return(1/Ne(t,...))
      }
    }
    if (is.null(Ne)){
      fun_inv<-function(t){
        return(1)
      }
    }
     
    #upper is an upper bound on fun_inv
    if (sampling=="iso"){
      
      if (is.null(sample)){
        sample<-c(n,0)
      }
      out<-replicate(N,rcoal(n,br=coalgen_thinning_iso(sample,fun_inv,max)$intercoal_times),simplify=FALSE)
      class(out)<-"multiPhylo"
      
      return(list(out=out,description="isochronous simulation using thinning with trajectory provided"))
      
      
    }
    
    else{
      #        #sampling is heterochronous
      if (is.null(sample)){ #it is actually isochronous
        sample<-c(n,0)
        out<-replicate(N,rcoal(n,br=coalgen_thinning_iso(sample,fun_inv,max)$intercoal_times),simplify=FALSE)
        class(out)<-"multiPhylo"
        
        return(list(out=out,description="isochronous simulation using thinning with trajectory provided"))
      }else{
        #It is indeed heterochronous
        out<-replicate(N,coalgen_thinning_hetero(sample, fun_inv,max),simplify=FALSE)     
        class(out)<-"multiPhylo" 
        return(list(out=out,description="Simulation with thinning for heterochronous coalescent with a specific function"))
        #   
      }
    }
  }
  
}