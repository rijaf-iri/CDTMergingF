#' CDT Merging Leave-One-Out Cross-Validation
#' 
#' This function perform a Leave-One-Out Cross-Validation.
#' 
#' @param time.step Time step of the data. Should be \strong{"daily"}, \strong{"pentad"}, \strong{"dekadal"} or \strong{"monthly"}.
#' @param start.date A vector of the start date to merge. 
#'  \itemize{
#'      \item \emph{Daily data}: \code{c(year, month, day)}
#'      \item \emph{Pentad data}: \code{c(year, month, pentad)},  \code{pentad} must be 1, 2, 3, 4, 5, or 6
#'      \item \emph{Dekadal data}: \code{c(year, month, dekad)},  \code{dekad} must be 1, 2 or 3
#'      \item \emph{Monthly data}: \code{c(year, month)}
#'  }
#' @param end.date A vector of the start date to merge. Same format as \code{start.date}.
#' @param station A list containing the station information in the form \code{list(file, sep, miss)}.
#' \itemize{
#'   \item \code{file}: full path to the station file (CDT station data format)
#'   \item \code{sep}: separator of the column
#'   \item \code{miss}: missing value flag
#' }
#' @param netcdf A list containing the netcdf data information in the form \code{list(dir, order.lon, order.lat, format)}.
#' \itemize{
#'  \item \code{dir}: full path to the directory containing the ne netcdf data
#'  \item \code{varid}: name of the variable to be used
#'  \item \code{order.lon}: order of the longitude (netcdf dimension order)
#'  \item \code{order.lat}: order of the latitude
#'  \item \code{format}: filename format of the netcdf files
#' }
#' @param merging.method Merging method to be used. Valid options: \strong{"CSc"}, \strong{"BSc"}, \strong{"SBA"} or \strong{"RK"}.
#' \itemize{
#'  \item \strong{"CSc"}: Cressman Scheme
#'  \item \strong{"BSc"}: Barnes Scheme
#'  \item \strong{"SBA"}: Simple Bias Adjustment
#'  \item \strong{"RK"}: Regression Kriging
#' }
#' @param interp.method Interpolation method to be used for \strong{"SBA"} and \strong{"RK"}.
#'  Valid options: \strong{"idw"}, \strong{"shepard"}, \strong{"spheremap"} or \strong{"kriging"}.
#' \itemize{
#'  \item \strong{"idw"}: Inverse distance weighted
#'  \item \strong{"shepard"}: Modified Shepard interpolation
#'  \item \strong{"spheremap"}: Spheremap interpolation method
#'  \item \strong{"kriging"}: Ordiranry kriging
#' }
#' @param spheric If \code{FALSE} (default), then a cartesian distance will be computed. If set to \code{TRUE}, a spherical distance (using a standard great circle method) will be computed.
#' @param maxdist Maximum radius of influence in decimal degree.
#' @param pass.ratio A vector giving the fraction of \code{maxdist} to be used for each pass.
#' @param pass.nmin A vector giving the minimum number of stations to be used to interpolate a grid point for each pass. Must be the same length as \code{pass.ratio}.
#' @param pass.nmax A vector giving the maximum number of stations to be used to interpolate a grid point for each pass. Must be the same length as \code{pass.ratio}.
#' @param neg.value If \code{TRUE}, negative values will be kept. If \code{FALSE} (default), negative values will be set to zero.
#' @param output.file Full path name to save the result.
#' @param use.RnoR If \code{TRUE}, apply rain-no-rain mask for rainfall data.
#' @param pars.RnoR A list of parameters to be used if \code{use.RnoR} is \code{TRUE}.
#' \itemize{
#'  \item \code{wet.day}: wet day definition
#'  \item \code{smooth}: if \code{TRUE}, smooth the mask
#' }
#' @param vgm.model A vector of variogram model to be used if \code{interp.method} is \strong{"kriging"}. Default is \code{c("Exp", "Gau", "Sph", "Pen")}.
#' @param parallel A list of parameters for parallel processing.
#' \itemize{
#' \item \code{dopar}: if \code{TRUE} (default), use parallel computing
#' \item \code{detect.cores}: if \code{TRUE} (default), detect the number of CPU cores.
#' \item \code{nb.cores}: if \code{dopar} is \code{TRUE} and \code{detect.cores} is \code{FALSE}, the number of CPU cores to be used 
#' }
#'
#' @examples
#' 
#' \dontrun{
#' cdtMergingLOOCV(
#'      time.step = "dekadal",
#'      start.date = c(2018, 1, 1),
#'      end.date = c(2018, 12, 3),
#'      station = list(file = "~/DONNEES/STN_DATA/CDT_RR_dekad_1981-2018_CLM.csv",
#'                  sep = ",", miss = "-99"),
#'      netcdf = list(dir = "~/DONNEES/GRIDDED_DATA/tamsat_v3/dekadal", varid = "rfe",
#'                  order.lon = 1, order.lat = 2, format = "rfe_%s%s%s.nc"),
#'      merging.method = "SBA", 
#'      interp.method = "idw",
#'      spheric = TRUE,
#'      maxdist = 1.0,
#'      pass.ratio = c(1, 0.75, 0.5),
#'      pass.nmin = c(5, 4, 3),
#'      pass.nmax = c(15, 10, 7),
#'      neg.value = FALSE,
#'      output.file = "~/DONNEES/OUT_Merging/LOOCV_SBA-idw.csv",
#'      use.RnoR = TRUE,
#'      pars.RnoR = list(wet.day = 1.0, smooth = FALSE),
#'      vgm.model = c("Exp", "Gau", "Sph", "Pen")
#'  )
#' }
#' 
#' @export

cdtMergingLOOCV <- function(
                        time.step = "dekadal",
                        start.date = c(1981, 1, 1),
                        end.date = c(2018, 12, 3),
                        station = list(file = NULL, sep = ",", miss = "-99"),
                        netcdf = list(dir = NULL, varid = NA, order.lon = 1, order.lat = 2, format = "rfe_%s%s%s.nc"),
                        merging.method = "SBA", 
                        interp.method = "idw",
                        spheric = FALSE,
                        maxdist = 2.5,
                        pass.ratio = c(1, 0.5, 0.2),
                        pass.nmin = c(5, 4, 3),
                        pass.nmax = c(15, 10, 7),
                        neg.value = FALSE,
                        output.file = NULL,
                        use.RnoR = FALSE,
                        pars.RnoR = list(wet.day = 1.0, smooth = FALSE),
                        vgm.model = c("Exp", "Gau", "Sph", "Pen"),
                        parallel = list(dopar = TRUE, detect.cores = TRUE, nb.cores = 2),
                        ...
                    )
{
    # test missing start and end date
    if(time.step == "monthly"){
        if(length(start.date) < 2) stop('The length of "start.date" must be 2')
        if(length(end.date) < 2) stop('The length of "en.date" must be 2')
        start.date <- c(start.date[1:2], 1)
        end.date <- c(end.date[1:2], 1)
    }else{
        if(length(start.date) < 3) stop('The length of "start.date" must be 2')
        if(length(end.date) < 3) stop('The length of "en.date" must be 2')
        start.date <- start.date[1:3]
        end.date <- end.date[1:3]
    }

    ##############

    xdeb <- as.Date(paste(start.date, collapse = '-'))
    xfin <- as.Date(paste(end.date, collapse = '-'))
    ###
    if(time.step == 'daily') daty <- seq(xdeb, xfin, 'day')
    if(time.step == 'pentad'){
        daty <- seq(xdeb, xfin, 'day')
        daty <- daty[as.numeric(format(daty, '%d')) <= 6]
    }
    if(time.step == 'dekadal'){
        daty <- seq(xdeb, xfin, 'day')
        daty <- daty[as.numeric(format(daty, '%d')) <= 3]
    }
    if(time.step == 'monthly') daty <- seq(xdeb, xfin, 'month')

    ###

    if(time.step == 'daily'){
        xdeb <- format(daty[1], '%Y%m%d')
        xfin <- format(daty[length(daty)], '%Y%m%d')
    }
    if(time.step %in% c('pentad', 'dekadal')){
        xdeb <- paste0(format(daty[1], '%Y%m'), as.numeric(format(daty[1], '%d')))
        xfin <- paste0(format(daty[length(daty)], '%Y%m'), as.numeric(format(daty[length(daty)], '%d')))
    }
    if(time.step == 'monthly'){
        xdeb <- format(daty[1], '%Y%m')
        xfin <- format(daty[length(daty)], '%Y%m')
    }

    ##############

    # Test missing output 
    if(is.null(output.file)) output.file <- file.path(getwd(), "LOOCV_output_file.csv")
    log.file <- file.path(dirname(output.file), "LOOCV_log_file.txt")

    ##############

    # Test missing station 
    stnData <- read.table(station$file, sep = station$sep, na.strings = station$miss, colClasses = "character", stringsAsFactors = FALSE)
    seph <- rle(grepl('[[:digit:]]', as.character(stnData[, 1])))
    ipos <- which(!seph$values & seph$lengths >= 3 & seph$lengths <= 4)
    if(length(ipos) == 0 | ipos[1] != 1) stop('Station data is not in a standard unambiguous CDT format')
    pos <- seph$lengths[ipos[1]]
    heads <- as.matrix(stnData[1:pos, ])

    stnData <- list(id = as.character(stnData[1, -1]),
                    lon = as.numeric(stnData[2, -1]),
                    lat = as.numeric(stnData[3, -1]),
                    elv = if(pos == 4) as.numeric(stnData[4, -1]) else NULL,
                    dates = as.character(stnData[-(1:pos), 1]),
                    data = local({
                                tmp <- stnData[-(1:pos), -1, drop = FALSE]
                                ntmp <- dim(tmp)
                                tmp <- as.numeric(unlist(tmp))
                                dim(tmp) <- ntmp
                                tmp
                        })
                    )

    ##############

    ## Test netcdf missing
    xdeb <- as.Date(paste(start.date, collapse = '-'))
    xfin <- as.Date(paste(end.date, collapse = '-'))
    errmsg <- "NetCDF data not found"
    ncInfo <- ncFilesInfo(time.step, xdeb, xfin, 1:12, netcdf$dir, netcdf$format, errmsg)
    if(is.null(ncInfo)) stop(paste(errmsg, "\n"))

    nc <- nc_open(ncInfo$nc.files[ncInfo$exist][1])
    grd.lon <- nc$var[[netcdf$varid]]$dim[[netcdf$order.lon]]$vals
    grd.lat <- nc$var[[netcdf$varid]]$dim[[netcdf$order.lat]]$vals
    nc_close(nc)
    xo <- order(grd.lon)
    yo <- order(grd.lat)
    grd.lon <- grd.lon[xo]
    grd.lat <- grd.lat[yo]
    nc.order <- list(ilon = netcdf$order.lon, ilat = netcdf$order.lat, olon = xo, olat = yo)

    ##############

    newgrid <- defSpatialPixels(list(lon = grd.lon, lat = grd.lat))
    ijGrd <- grid2pointINDEX(list(lon = stnData$lon, lat = stnData$lat),
                            list(lon = grd.lon, lat = grd.lat))
    out.stn <- rep(NA, length(ijGrd))

    ##############

    parsL = c(condition = length(which(ncInfo$exist)) >= 20, parallel)

    exports <- list(...)$export

    out.cv <- cdt.foreach(seq_along(ncInfo$nc.files), parsL = parsL,
                    .packages = c('sp', 'ncdf4'),
                    .export = exports,
                    FUN = function(jj)
    {
        if(ncInfo$exist[jj]){
            nc <- nc_open(ncInfo$nc.files[jj])
            nc.val <- ncvar_get(nc, varid = netcdf$varid)
            nc_close(nc)
            nc.val <- transposeNCDFData(nc.val, nc.order)
        }else{
            cat(paste(ncInfo$dates[jj], ":", "no netcdf data", "|",
                "no file generated", "\n"), file = log.file, append = TRUE)
            return(out.stn)
        }

        if(all(is.na(nc.val))){
            cat(paste(ncInfo$dates[jj], ":", "all netcdf data are missing", "|",
                "no file generated", "\n"), file = log.file, append = TRUE)
            return(out.stn)
        }

        ######
        newgrid$grd <- c(nc.val)

        ######
        donne.stn <- stnData$data[stnData$date == ncInfo$dates[jj], , drop = FALSE]
        if(nrow(donne.stn) == 0 | ncol(donne.stn) < pass.nmin[length(pass.nmin)]){
            cat(paste(ncInfo$dates[jj], ":", "no station data", "|",
                "no file generated", "\n"), file = log.file, append = TRUE)
            return(out.stn)
        }

        locations.stn <- data.frame(lon = stnData$lon, lat = stnData$lat, stn = c(donne.stn[1, ]))
        coordinates(locations.stn) <- ~lon+lat
        noNA <- !is.na(locations.stn$stn)
        locations.stn <- locations.stn[noNA, ]
        ijNA <- ijGrd[noNA]

        if(length(locations.stn) < pass.nmin[length(pass.nmin)]){
            cat(paste(ncInfo$dates[jj], ":", "no station data", "|",
                "no file generated", "\n"), file = log.file, append = TRUE)
            return(out.stn)
        }

        xstn <- lapply(seq_along(locations.stn), function(ii){
            loc.stn <- locations.stn[-ii, ]
            out.mrg <- merging.functions(loc.stn, newgrid, 
                                        merging.method, interp.method,
                                        maxdist, pass.ratio, pass.nmin, pass.nmax, 
                                        vgm.model, spheric, ncInfo$dates[jj],
                                        neg.value, log.file)
            if(use.RnoR){
                rnr <- rain_no_rain.mask(loc.stn, newgrid, pars.RnoR, pass.nmax)
                out.mrg <- out.mrg * rnr
            }

            out.mrg[ijNA[ii]]
        })
        xstn <- do.call(c, xstn)
        out.stn[noNA] <- round(xstn, 1)

        return(out.stn)
    })

    out.cv <- do.call(rbind, out.cv)
    out.cv <- cbind(ncInfo$dates, out.cv)
    out.cv <- rbind(heads, out.cv)

    write.table(out.cv, output.file, sep = station$sep, na = station$miss,
                col.names = FALSE, row.names = FALSE, quote = FALSE)

    invisible()
}
