#'@title .calculate.moller.hetero
#'@description Approximates the posterior distribution of Ne from a single genealogy at a regular grid of points using INLA package
#'@param coal.factor is a vector with coalescent times in increasing order 
#'@param s sampling times and coalescent times in increasing order
#'@param event and indicator vector with 1 for coalescent times and 0 for sampling times that correspond to the s vector
#'@param lengthout number of grid points
#'@param prec_alpha alpha gamma hyperparameter of precision parameter for GP prior
#'@param prec_beta  beta gamma hyperparameter of precision parameter for GP prior
#'@param E.log.zero internal
#'@param alpha TODO
#'@param beta TODO
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}

.calculate.moller.hetero<-function (coal.factor, s, event, lengthout, prec_alpha = 0.01, 
          prec_beta = 0.01, E.log.zero = -100, alpha = NULL, beta = NULL) 
{
  
  if (prec_alpha == 0.01 & prec_beta == 0.01 & !is.null(alpha) & 
        !is.null(beta)) {
    prec_alpha = alpha
    prec_beta = beta
  }
  grid <- seq(0, max(s), length.out = lengthout + 1)
  u <- diff(grid)
  field <- grid[-1] - u/2
  sgrid <- grid
  event_new <- 0
  time <- 0
  where <- 1
  E.factor <- 0
  for (j in 1:lengthout) {
    count <- sum(s > sgrid[j] & s <= sgrid[j + 1])
    if (count > 1) {
      points <- s[s > sgrid[j] & s <= sgrid[j + 1]]
      u <- diff(c(sgrid[j], points))
      event_new <- c(event_new, event[(where):(where + 
                                                 count - 1)])
      time <- c(time, rep(field[j], count))
      E.factor <- c(E.factor, coal.factor[where:(where + 
                                                   count - 1)] * u)
      where <- where + count
      if (max(points) < sgrid[j + 1]) {
        event_new <- c(event_new, 0)
        time <- c(time, field[j])
        E.factor <- c(E.factor, coal.factor[where] * 
                        (sgrid[j + 1] - max(points)))
      }
    }
    if (count == 1) {
      event_new <- c(event_new, event[where])
      points <- s[s > sgrid[j] & s <= sgrid[j + 1]]
      if (points == sgrid[j + 1]) {
        E.factor <- c(E.factor, coal.factor[where] * 
                        (sgrid[j + 1] - sgrid[j]))
        time <- c(time, field[j])
        where <- where + 1
      }
      else {
        event_new <- c(event_new, 0)
        E.factor <- c(E.factor, coal.factor[where] * 
                        (points - sgrid[j]))
        E.factor <- c(E.factor, coal.factor[where + 1] * 
                        (sgrid[j + 1] - points))
        time <- c(time, rep(field[j], 2))
        where <- where + 1
      }
    }
    if (count == 0) {
      event_new <- c(event_new, 0)
      E.factor <- c(E.factor, coal.factor[where] * (sgrid[j + 
                                                            1] - sgrid[j]))
      time <- c(time, field[j])
    }
  }
  time2 <- time
  event_new2 <- event_new
  E.factor2 <- E.factor
  for (j in 1:lengthout) {
    count <- sum(time2 == field[j])
    if (count > 1) {
      indic <- seq(1:length(event_new2))[time2 == field[j]]
      if (sum(event_new2[indic]) == 0) {
        event_new2 <- event_new2[-indic[-1]]
        time2 <- time2[-indic[-1]]
        temp <- sum(E.factor2[indic])
        E.factor2[indic[1]] <- temp
        E.factor2 <- E.factor2[-indic[-1]]
      }
    }
  }
  E.factor2.log = log(E.factor2)
  E.factor2.log[E.factor2 == 0] = E.log.zero
  data <- list(y = event_new2[-1], event = event_new2[-1], 
               time = time2[-1], E = E.factor2.log[-1])
  formula <- y ~ -1 + f(time, model = "rw1", hyper = list(prec = list(param = c(prec_alpha, 
                                                                                prec_beta))), constr = FALSE)
  mod4 <- inla(formula, family = "poisson", data = data, offset = E, 
               control.predictor = list(compute = TRUE))
  return(list(result = mod4, grid = grid, data = data, E = E.factor2.log))
}

#'@title .calculate.moller
#'@description Approximates the posterior distribution of Ne from a single genealogy at a regular grid of points using INLA package
#'@param data1 is a dataframe with two columns. The first column has the intercoalescent times and the second column the number of lineages
#'@param lengthout number of grid points
#'@param L the length for the definition of the grid
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}

.calculate.moller<-function(data1,lengthout,L){
  s<-cumsum(data1[,1])
  coal.factor<-data1[,2]*(data1[,2]-1)*.5
  length.min<-min(data1[,1])
  #maxval<-sum(data1[,1])
  grid<-seq(0,L,length.out=lengthout+1)
  #grid<-c(grid[grid<=maxval],maxval)
  u<-diff(grid)
  field<-grid[-1]-u/2
  sgrid<-grid
  event<-0
  time<-0
  E.factor<-0
  where<-1
  for (j in 1:lengthout){
    count<-sum(s>sgrid[j] & s<=sgrid[j+1]) 
    if (count>1){
      points<-s[s>sgrid[j] & s<=sgrid[j+1]]
      u<-diff(c(sgrid[j],points))
      event<-c(event,rep(1,count))
      time<-c(time,rep(field[j],count))
      E.factor<-c(E.factor,coal.factor[where:(where+count-1)]*u)
      where<-where+count
      if (max(points)<sgrid[j+1]){
        event<-c(event,0)
        time<-c(time,field[j])
        E.factor<-c(E.factor,coal.factor[where]*(sgrid[j+1]-max(points)))
      }
    }
    if (count==1){
      event<-c(event,1)
      points<-s[s>sgrid[j] & s<=sgrid[j+1]]
      if (points==sgrid[j+1]){
        E.factor<-c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
        time<-c(time,field[j])
        where<-where+1
      }else {
        event<-c(event,0)
        E.factor<-c(E.factor,coal.factor[where]*(points-sgrid[j]))
        E.factor<-c(E.factor,coal.factor[where+1]*(sgrid[j+1]-points))
        time<-c(time,rep(field[j],2))
        where<-where+1}
    }
    if (count==0){
      event<-c(event,0)
      E.factor<-c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
      time<-c(time,field[j])
    }
    
  }
  ##Fixing for when there are no observations
  combine<-unique(sort(c(grid[-1],s)))
  find2<-max(seq(1,length(combine))[combine<=max(s)])  
  event<-event[-1]
  time<-time[-1]
  E.factor<-E.factor[-1]
  grid<-c(grid[grid<=max(s)],max(s))
  data<-list(y=event[1:find2],event=event[1:find2],time=time[1:find2],E=log(E.factor[1:find2]))
  
  #  data<-list(y=event[-1],event=event[-1],time=time[-1],E=log(E.factor[-1]))
  formula<-y~-1+f(time,model="rw1",hyper=list(prec = list(param = c(.001, .001))),constr=FALSE)
  mod.moller.constant<-inla(formula,family="poisson",data=data,offset=E,control.predictor=list(compute=TRUE))
  
  return(stan_output(list(result=mod.moller.constant,grid=grid)))
    
}


#'@title plot_INLA
#'@description Plots the output from the inla functions for Ne
#'@param INLA_out Otput from the inla functions
#'@param traj the true trajectory
#'@param xlim
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}

plot_INLA = function(INLA_out, traj=NULL, xlim=NULL, ...)
{
  mod = INLA_out$result$summary.random$time
  
  grid = mod$"ID"
  if (is.null(xlim))
  {
    xlim=c(max(grid),0)
  }
  plot(grid,exp(-mod$"0.5quant"),type="l",lwd=2.5,col="blue",log="y",
       xlab="Time (past to present)",ylab="Scaled Effective Pop. Size",
       xlim=xlim,ylim=c(min(exp(-mod$"0.975quant"[grid > min(xlim) & grid < max(xlim)])),
                        max(exp(-mod$"0.025quant"[grid > min(xlim) & grid < max(xlim)]))), ...)
  lines(grid,exp(-mod$"0.975quant"),lwd=2.5,col="blue",lty=2)
  lines(grid,exp(-mod$"0.025quant"),lwd=2.5,col="blue",lty=2)
  if (!is.null(traj))
    lines(grid, traj(grid))
}

stan_output<-function(INLA_out){
  mod = INLA_out$result$summary.random$time
  grid = mod$"ID"
  
  results<-matrix(NA,nrow=length(grid),ncol=4)
  results[,1]<-grid
  results[,3]<-exp(-mod$"0.975quant")
  results[,2]<-exp(-mod$"0.025quant")
  results[,4]<-exp(-mod$"0.5quant")
  return(results)
}
