library(CompQuadForm)
library(MASS)
library(bestNormalize)

maf <- function(g) {a <- colMeans(g, na.rm=TRUE) / 2; pmin(a, 1 - a)}
if(!exists("c17"))
{
    c17 <- readRDS("17q12.rds")
    n17 <- nrow(c17)
    m17 <- ncol(c17)
    f17 <- maf(c17)
}

for(. in dir("R", "[.]R$", ful=TRUE))
{
    source(.)
}

for(. in grep("ini.R", dir(".", "[.]R$"), value = TRUE, invert = TRUE))
{
    source(.)
}
