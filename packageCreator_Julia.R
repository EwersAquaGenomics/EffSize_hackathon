require(devtools)

use_package("ape")
use_package("phangorn")
use_package("phyclust")
use_package("adegenet")
use_package("INLA")

# create("multiNe")

setwd("~/Documents/Hackathon/multiNe")

Rmd_files=c("SkylinePlottingExample.Rmd",
            "Simulation_Genealogies.Rmd",
            "INLAExample.Rmd")


sapply(Rmd_files,function(x) system(paste("cp ~/Documents/Hackathon/EffSize_hackathon/",x," R/",sep="")))

document()
check()
install()



