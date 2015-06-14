#### testing functions ###

library(adegenet)
data(nancycats)


### LDNe_051815 ####

LDNe(nancycats, crit=0.03, mating="mono")
# 31.98824


### HENe_052315 ####

HENe(g0, crit=0.01)
# -1726.6 -> there is no Heterozygote excess?
# insensitive to low frequency alleles
HENe(g1, crit=0.01)
HENe(g2, crit=0)


### variance Ne ####

varNe(g)
#   between.generations       Ne
#12                 1-2 15.64455
#13                 1-3 16.44186
#23                 2-3 14.03297


