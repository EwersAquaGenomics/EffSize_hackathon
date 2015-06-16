#### testing functions ###

library(adegenet)
data(nancycats)
source("sim.R")


### LDNe_051815 ####

LDNe(nancycats) # takes some time
#point_estimate         median           2.5%          97.5%       variance 
#105.66359      105.20562       93.99704      114.29843       20.54779 

LDNe(nancycats, crit=0.03) # takes a minute or so
#point_estimate         median           2.5%          97.5%       variance 
#-67.9677954    -67.6958719    -68.9073450    -66.1378856      0.9308135 
#There were 50 or more warnings (use warnings() to see the first 50)

LDNe_point(nancycats, crit=0.03)
#[1] -67.9678
#There were 50 or more warnings (use warnings() to see the first 50)
warnings()
#50: In sqrt((p[[i]][k] * (1 - p[[i]][k]) + (p.hom[[i]][k] -  ... : NaNs produced

LDNe_point(nancycats, crit=0.01)
#[1] 105.6636

LDNe_point(nancycats, crit=0.02)
#[1] -131.367

genind.by.locus <- seploc(nancycats)
ps <- makeP(genind.by.locus)
pps <- unlist(ps) 
hist(pps, breaks = 30)
sapply(ps, mean)

LDNe(simG, crit=0.02, mating="random")
#point_estimate         median           2.5%          97.5%       variance 
#40.75978       37.69495       15.13091      112.37229    15364.28135 


### HENe_052315 ####

HENe(g0, crit=0.01)
# -1726.6 -> there is no Heterozygote excess?
# insensitive to low frequency alleles
HENe(g1, crit=0.01)
HENe(g2, crit=0)


### variance Ne ####

varNe(simG3)
#   between.generations point.estimate   median    0.025    0.975    variance
#12                 1-2       1.727976 1.727976 1.629624 1.817130 0.002412062
#13                 1-3       3.302347 3.302347 3.056548 3.433667 0.008607956
#23                 2-3       5.174282 5.174282 3.976703 6.861758 0.454806603

varNe_point(simG3)
#between.generations       Ne
#12                 1-2 1.727976
#13                 1-3 3.302347
#23                 2-3 5.174282
