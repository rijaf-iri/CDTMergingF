
rain_no_rain.mask <- function(locations.stn, newgrid, nmax)
{
    glm.binom <- tryCatch(
            glm(rnr.stn ~ grd, data = locations.stn, family = binomial(link = "logit")),
            error=function(e) e, warning=function(w) w
        )

    if(inherits(glm.binom, "warning") | inherits(glm.binom, "error")) return(NULL)

    rnr <- NULL
    if(!is.na(glm.binom$coef[2])){
        locations.stn$rnr.res <- residuals(glm.binom)
        rnr.trend <- predict(glm.binom, newdata = newgrid, type = 'link')

        rnr.res.grd <- gstat::krige(rnr.res~1, locations = locations.stn, newdata = newgrid, nmax = nmax, debug.level = 0)
        rnr <- rnr.trend + rnr.res.grd$var1.pred

        rnr <- exp(rnr) / (1 + exp(rnr))
        ### decision boundary 0.5
        rnr[rnr >= 0.5] <- 1
        rnr[rnr < 0.5] <- 0
        rnr[is.na(rnr)] <- 1
    }

    return(rnr)
}
