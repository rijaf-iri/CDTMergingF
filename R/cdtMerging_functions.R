
merging.functions <- function(locations.stn, newgrid, 
                            merging.method, interp.method,
                            maxdist, pass.ratio, pass.nmin, pass.nmax,
                            vgm.model, spheric, nc.date, neg.value,
                            use.RnoR, pars.RnoR, log.file)
{
    interp.method <- switch(merging.method,
                            "CSc" = "cressman",
                            "BSc" = "barnes",
                            interp.method)
    if(use.RnoR){
        wet.day <- pars.RnoR$wet.day
        if(wet.day <= 0) wet.day <- wet.day + 1e-13
        locations.stn$rnr.stn <- ifelse(locations.stn$stn < wet.day, 0, 1)
    }

    nlon <- newgrid@grid@cells.dim[1]
    nlat <- newgrid@grid@cells.dim[2]

    nrun <- if(length(pass.ratio) > 1) length(pass.ratio) else 2
    fac <- c(4, 3, rep(2, nrun - 2))

    for(pass in seq_along(pass.ratio)){
        rmax <- maxdist * pass.ratio[pass]
        nmin <- pass.nmin[pass]
        nmax <- pass.nmax[pass]

        xy.grid <- create_grid_buffer(locations.stn, newgrid, rmax, fac[pass], spheric)
        # xy.grid <- create_grid_buffer(locations.stn, newgrid, rmax, spheric)
        newdata0 <- xy.grid$grid.buff

        igrid <- xy.grid$ij
        coarsegrid <- xy.grid$coarse
        coarsegrid$res <- rep(0, length(coarsegrid))

        locations.stn$grd <- over(locations.stn, newdata0)$grd
        locations.stn <- locations.stn[!is.na(locations.stn$grd), ]

        if(length(locations.stn) < nmin){
            cat(paste(nc.date, ":", paste("not enough station data pass#", pass), "|",
                "Output: gridded data", "\n"), file = log.file, append = TRUE)
            out.mrg <- matrix(newgrid@data$grd,
                            ncol = newgrid@grid@cells.dim[2],
                            nrow = newgrid@grid@cells.dim[1])
            return(out.mrg)
        }

        ######### sp.trend & residuals

        ## CSc, BSc, SBA
        sp.trend <- newdata0@data$grd
        xres <- locations.stn$stn - locations.stn$grd

        ## RK
        if(merging.method == "RK"){
            if(var(locations.stn$stn) < 1e-07 | var(locations.stn$grd, na.rm = TRUE) < 1e-07){
                cat(paste(nc.date, ":", paste("Zero variance @ GLM pass#", pass), "|",
                    "Simple Bias Adjustment", "\n"), file = log.file, append = TRUE)
            }else{
                glm.stn <- glm(stn ~ grd, data = locations.stn, family = gaussian)
                if(is.na(glm.stn$coefficients[2]) | glm.stn$coefficients[2] < 0){
                    cat(paste(nc.date, ":", paste("Invalid GLM coeffs pass#", pass), "|",
                        "Simple Bias Adjustment", "\n"), file = log.file, append = TRUE)
                }else{
                    sp.trend <- predict(glm.stn, newdata = newdata0)
                    ina.out <- is.na(sp.trend)
                    sp.trend[ina.out] <- newdata0@data$grd[ina.out]
                    xres <- rep(NA, length(locations.stn))
                    if(length(glm.stn$na.action) > 0)
                        xres[-glm.stn$na.action] <- glm.stn$residuals
                    else
                        xres <- glm.stn$residuals
                }
            }
        }

        #########

        loc.stn <- locations.stn
        loc.stn$res <- xres
        loc.stn <- loc.stn['res']

        #########

        if(interp.method == 'kriging'){
            calc.vgm <- if(length(loc.stn$res) > 7 & var(loc.stn$res) > 1e-15) TRUE else FALSE
            if(calc.vgm){
                ## spherical distance
                vgm <- try(automap::autofitVariogram(res~1, input_data = loc.stn, model = vgm.model, cressie = TRUE), silent = TRUE)
                if(!inherits(vgm, "try-error")){
                    vgm <- vgm$var_model
                    if(spheric){
                        ctr <- rowSums(loc.stn@bbox)/2
                        pts <- c(ctr[1] + vgm$range[2] * cos(pi/4), ctr[2] + vgm$range[2] * sin(pi/4))
                        vgm$range[2] <- distance.vector(pts, matrix(ctr, nrow = 1), spheric)
                    }
                }else{
                    cat(paste(nc.date, ":", paste("Unable to compute variogram pass#", pass), "|",
                        "Interpolation using IDW", "\n"), file = log.file, append = TRUE)
                    interp.method <- "idw"
                }
            }else{
                cat(paste(nc.date, ":", paste("Unable to compute variogram pass#", pass), "|",
                    "Interpolation using IDW", "\n"), file = log.file, append = TRUE)
                interp.method <- "idw"
            }
        }

        #########

        row.names(loc.stn) <- 1:length(loc.stn)
        row.names(coarsegrid) <- length(loc.stn) + (1:length(coarsegrid))
        loc.stn <- maptools::spRbind(loc.stn, coarsegrid)

        #########
        interp.res <- switch(interp.method,
            "barnes" = barnes.interp(loc.stn@coords, loc.stn$res, newdata0@coords, nmin, nmax, spheric, p = 2),
            "cressman" = cressman.interp(loc.stn@coords, loc.stn$res, newdata0@coords, nmin, nmax, spheric),
            "idw" = idw.interp(loc.stn@coords, loc.stn$res, newdata0@coords, nmin, nmax, spheric, p = 2),
            "shepard" = shepard.interp(loc.stn@coords, loc.stn$res, newdata0@coords, nmin, nmax, spheric, p = 2),
            "spheremap" = spheremap.interp(loc.stn@coords, loc.stn$res, newdata0@coords, nmin, nmax, spheric),
            "kriging" = kriging.interp(loc.stn@coords, loc.stn$res, newdata0@coords, vgm, nmin, nmax, spheric)
        )

        #########

        out.mrg <- newgrid@data$grd
        out.mrg[igrid] <- sp.trend + interp.res[, 3]
        if(!neg.value) out.mrg[out.mrg < 0] <- 0
        ina <- is.na(out.mrg)
        out.mrg[ina] <- newgrid@data$grd[ina]

        #########

        if(use.RnoR){
            newdata0$grd <- out.mrg[igrid]
            rnr0 <- rain_no_rain.mask(locations.stn, newdata0, nmax)
            if(!is.null(rnr0)){
                rnr <- matrix(1, ncol = nlat, nrow = nlon)
                rnr[igrid] <- rnr0

               if(pars.RnoR$smooth){
                    if(pass == length(pass.ratio)){
                        ij <- over(locations.stn, as(newgrid, "SpatialPixels"))
                        rnr[ij] <- locations.stn$rnr.stn
                    }
                    rnr <- (2 * rnr + smooth.matrix(rnr, 2))/3
                }
                out.mrg <- out.mrg * c(rnr)
            }
        }

        #########

        newgrid@data$grd <- out.mrg
    }

    out.mrg <- matrix(out.mrg, ncol = nlat, nrow = nlon)

    return(out.mrg)
}
