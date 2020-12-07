library(ggplot2)

.th <- theme(
    ## axis.title.x=element_blank(), axis.title.y=element_blank(), 
    strip.text.x = element_text(size=12, face="bold"),
    strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="red", fill="#CCCCFF"),
    legend.title=element_blank(), legend.position='bottom')

## cap the values
.cp <- function(dat, grp, val='val', cap=0.01, mtd=c('both', 'upper', 'lower'))
{
    grp <- split(dat, dat[, grp])
    mtd <- match.arg(mtd, c('both', 'upper', 'lower'))
    grp <- lapply(grp, function(g)
    {
        v <- g[, val]
        if(mtd == 'upper')
            v <- pmin(v, quantile(v, 1-cap, na.rm=TRUE))
        else if(mtd == 'lower')
            v <- pmax(v, quantile(v, 0+cap, na.rm=TRUE))
        else
        {
            v <- pmin(v, quantile(v, 1-cap/2, na.rm=TRUE))
            v <- pmax(v, quantile(v, 0+cap/2, na.rm=TRUE))
        }
        g[, val] <- v
        g
    })
    dat <- do.call(rbind, grp)
    dat
}

get.rpt <- function(sim, cache=TRUE)
{
    rds=paste0(sim, '.rds')
    if(file.exists(rds) && cache)
        agg <- readRDS(rds)
    else
    {
        agg <- lapply(dir(sim, '^[0-9]+.rds$', full=TRUE), function(f)
        {
            print(f)
            readRDS(f)
        })
        agg <- do.call(rbind, agg)
        saveRDS(agg, rds)
    }
    invisible(agg)
}

get.pow <- function(sim, cache=TRUE)
{
    rds=paste0(sim, '.pow')
    if(file.exists(rds) && cache)
        pow <- readRDS(rds)
    else
    {
        rpt <- get.rpt(sim)
        grp <- subset(rpt, se=-c(pow, egv, mcr))
        pow <- by(rpt, grp, function(g)
        {
            cfg <- subset(g, se=-c(pow, egv, rep))[1, ]
            pow <- with(g, sum(pow * rep) / sum(rep))
            egv <- with(g, sum(egv * rep) / sum(rep))
            mcr <- with(g, sum(mcr * rep) / sum(rep))
            rep <- with(g, sum(rep))
            cbind(cfg, pow=pow, egv=egv, rep=rep)
        })
        pow <- do.call(rbind, pow)
        saveRDS(pow, rds)
    }
    invisible(pow)
}

get.cfg <- function(rpt)
{
    dct <- c(N="N", a="Genetics", b="Covariates", d="GxU", e="Noise",
             evt="Min-Eigen", psd="PD-Threshold", rep="Repeats")
    unq <- sapply(lapply(rpt, unique), length) < 2
    cfg <- rpt[1, unq]
    str <- dct[colnames(cfg)]
    msg <- paste(paste(str, cfg, sep="="), collapse = ", ")
    msg
}

plt.pow <- function(sim, sub, out=paste0(sim, '.pdf'), ttl=NULL)
{
    rpt <- get.pow(sim)
    if(!missing(sub))
    {
        sub <- substitute(sub)
        sub <- eval(sub, rpt, parent.frame())
        rpt <- rpt[sub, ]
    }
    rpt <- within(rpt,
    {
        ## key <- ifelse(d > 0, "H1: GxU > 0", "H0: GxU = 0")
        key <- ifelse(a > 0, "H1: G > 0", "H0: G = 0")
        lhs <- base::sub("^LHS = ", "", lhs)
        tag <- mapply(base::sub, list("^LHS"), lhs, mdl)
    })
    
    g <- ggplot(rpt, aes(x=M, y=pow))
    g <- g + geom_line(aes(color=mtd), alpha=.5, size=1)
    g <- g + facet_grid(key ~ tag)
    g <- g + labs(title=ttl, subtitle=get.cfg(rpt))
    g <- g + xlab("Number of Traits")
    g <- g + ylab("Ratio of P < 0.05")
    g <- g + .th

    nfy <- length(unique(rpt$key))
    ufy <- 5
    nfx <- length(unique(rpt$tag))
    ufx <- 5
    ## if(ufx / ufy < 19 / 10)
    ##     ufy <- ufx / 19 * 10
    ## else
    ##     ufx <- ufy / 10 * 19
    options(bitmapType = 'cairo', device = 'pdf')
    ggsave(out, g, width=min(19, ufx * nfx), height=min(10, ufy * nfy), dpi=400, scale=.7)
    invisible(g)
}

plt.main <- function()
{
    plt.pow("run/mx6", out="run/mx6.png")
    plt.pow("run/mx6", out="run/mx6.pdf")
    plt.pow("run/mx7", out="run/mx7.png")
    plt.pow("run/mx7", out="run/mx7.pdf")
}
