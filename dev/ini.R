library(CompQuadForm)
library(MASS)
## library(dglm)
## library(bestNormalize)
load("c17")
for(. in dir("R", "[.]R$", ful=TRUE))
{
    source(.)
}

for(. in grep("ini.R", dir(".", "[.]R$"), value = TRUE, invert = TRUE))
{
    source(.)
}
