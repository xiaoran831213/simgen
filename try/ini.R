library(CompQuadForm)
## library(SKAT)
library(MASS)
library(Matrix)


for(. in dir("R", "[.]R$", ful=TRUE))
{
    source(.)
}

for(. in grep("ini.R", dir(".", "[.]R$"), value = TRUE, invert = TRUE))
{
    source(.)
}
