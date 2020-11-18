## library(CompQuadForm)
## library(SKAT)
library(MASS)


## source("sim.R")
## source("gen.R")
## source("mtd.R")
## source("imp.R")
## source("err.R")
## source("cor.R")
source("tmp.R")

for(. in dir("../R", "[.]R$", ful=TRUE))
{
    source(.)
}

for(. in grep("ini.R", dir(".", "[.]R$"), value = TRUE, invert = TRUE))
{
    source(.)
}
