#### testing functions ###

library(adegenet)
data(nancycats)

getSampleSize_byLocus(nancycats)
# works

makeP(nancycats)
# works

makeCrit(nancycats, crit=c(0, 0.05))
# works

LDNe(nancycats, crit=c(0.05, 0.01), mating="mono")
#invalid length for loc.fac
#Error in validObject(x) : invalid class “genind” object: FALSE
#Called from: (function () 
#{
#  .rs.breakOnError(TRUE)
#})()

HENe(nancycats)
#Error: (list) object cannot be coerced to type 'double'