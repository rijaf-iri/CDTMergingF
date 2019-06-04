
rain_no_rain.mask <- function(locations.stn, newgrid, pars.RnoR, pass.nmax)
{
    wet.day <- pars.RnoR$wet.day
    if(wet.day <= 0) wet.day <- wet.day + 1e-13
    ij <- over(locations.stn, as(newgrid, "SpatialPixels"))
    locations.stn$rnr.stn <- ifelse(locations.stn$stn < wet.day, 0, 1)
    locations.stn$grd <- newgrid@data$grd[ij]
    locations.stn <- locations.stn[!is.na(locations.stn$grd), ]

    glm.binom <- glm(rnr.stn ~ grd, data = locations.stn, family = binomial(link = "logit"))

    nlon <- newgrid@grid@cells.dim[1]
    nlat <- newgrid@grid@cells.dim[2]
    rnr <- matrix(1, ncol = nlat, nrow = nlon)
    if(!is.na(glm.binom$coef[2])){
        nmax <- pass.nmax[length(pass.nmax)]

        locations.stn$rnr.res <- residuals(glm.binom)
        rnr.trend <- predict(glm.binom, newdata = newgrid, type = 'link')

        rnr.res.grd <- gstat::krige(rnr.res~1, locations = locations.stn, newdata = newgrid, nmax = nmax, debug.level = 0)
        rnr <- rnr.trend + rnr.res.grd$var1.pred

        rnr <- exp(rnr) / (1 + exp(rnr))
        ### decision boundary 0.5
        rnr[rnr >= 0.5] <- 1
        rnr[rnr < 0.5] <- 0
        rnr <- matrix(rnr, ncol = nlat, nrow = nlon)
        rnr[is.na(rnr)] <- 1
        if(pars.RnoR$smooth){
            rnr[ij] <- locations.stn$rnr.stn
            rnr <- (2 * rnr + smooth.matrix(rnr, 2))/3
        }
    }

    return(rnr)
}
