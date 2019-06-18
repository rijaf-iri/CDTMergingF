#' Distance matrix
#' 
#' Distance matrix between two sets of points
#'
#' @param x a 2 column matrix or data.frame (first column is longitude, second is latitude); or a SpatialPoints*, SpatialPixels* or SpatialGrid* object.
#' @param y Same as \code{x}. If \code{x} is a Spatial object, \code{y} can take a different Spatial object other than \code{x}.
#' @param spheric If \code{FALSE} (default), then a Cartesian distance will be computed. If set to \code{TRUE}, a spherical distance (using a standard great circle method) will be computed.
#' @return Matrix of distances. The row represent \code{y} and the column \code{x}.
#' 
#' @export

distance.Matrix <- function(x, y, spheric = FALSE){
    if(inherits(x, "matrix") && inherits(y, "matrix")){
        x <- x[, 1:2, drop = FALSE]
        y <- y[, 1:2, drop = FALSE]
    }else if(inherits(x, "data.frame") && inherits(y, "data.frame")){
        x <- as.matrix(x[, 1:2, drop = FALSE])
        y <- as.matrix(y[, 1:2, drop = FALSE])
    }else{
        class.x <- substr(class(x), 1, 7)
        class.y <- substr(class(y), 1, 7)
        if(class.x == "Spatial" && class.y == "Spatial"){
            if(!identical(proj4string(x), proj4string(y)))
                stop("x and y have different coordinate reference systems")
            x <- as.matrix(coordinates(x))
            y <- as.matrix(coordinates(y))
        }else{
            stop("x and y must be a matrix, data.frame or sp Spatial object")
        }
    }

    distance.matrix(x, y, spheric)
}

distance.vector <- function(x, y, spheric){
    # x: vector c(x, y)
    # y: matrix X Y
    x <- as.numeric(x)
    y <- as.matrix(y)
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    nr <- nrow(y)
    out <- vector("numeric", length = nr)
    out <- .Fortran("distance_vector", x, y, as.integer(nr), as.integer(spheric), dst = out)
    out$dst
}

distance.matrix <- function(x, y, spheric){
    # x: matrix X Y
    # y: matrix X Y
    x <- as.matrix(x)
    y <- as.matrix(y)
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    n1 <- nrow(x)
    n2 <- nrow(y)
    dst <- matrix(double(1), nrow = n2, ncol = n1)
    out <- .Fortran("distance_matrix", x, y, as.integer(n1), as.integer(n2),
                    as.integer(spheric), dst = dst)
    out$dst
}

########################################

defSpatialPixels <- function(grd_Coords){
    # grd_Coords: named list(lon, lat)
    newgrid <- expand.grid(lon = grd_Coords$lon, lat = grd_Coords$lat)
    coordinates(newgrid) <- ~lon+lat
    newgrid <- try(SpatialPixels(points = newgrid,
                            tolerance = sqrt(sqrt(.Machine$double.eps)),
                            proj4string = CRS(as.character(NA))), silent = TRUE)
    if(inherits(newgrid, "try-error")){
        newgrid <- expand.grid(lon = grd_Coords$lon, lat = grd_Coords$lat)
        coordinates(newgrid) <- ~lon+lat
        newgrid <- SpatialPixels(points = newgrid, tolerance = 0.001,
                                proj4string = CRS(as.character(NA)))
    }

    return(newgrid)
}

grid2pointINDEX <- function(pts_Coords, grd_Coords){
    # grd_Coords: named list(lon, lat)
    # pts_Coords: named list(lon, lat)
    newgrid <- expand.grid(lon = grd_Coords$lon, lat = grd_Coords$lat)
    coordinates(newgrid) <- ~lon+lat
    newgrid <- try(SpatialPixels(points = newgrid,
                            tolerance = sqrt(sqrt(.Machine$double.eps)),
                            proj4string = CRS(as.character(NA))), silent = TRUE)
    if(inherits(newgrid, "try-error")){
        newgrid <- expand.grid(lon = grd_Coords$lon, lat = grd_Coords$lat)
        coordinates(newgrid) <- ~lon+lat
        newgrid <- SpatialPixels(points = newgrid, tolerance = 0.001,
                                proj4string = CRS(as.character(NA)))
    }

    pts.loc <- data.frame(lon = pts_Coords$lon, lat = pts_Coords$lat)
    pts.loc <- SpatialPoints(pts.loc)
    ijGrd <- unname(over(pts.loc, geometry(newgrid)))
    return(ijGrd)
}

########################################

# smooth.matrix <- function(mat, ns){
#     rm <- nrow(mat) + 2 * ns
#     cm <- ncol(mat) + 2 * ns
#     M <- matrix(NA, rm, cm)
#     sqC <- (ns + 1):(cm - ns)
#     sqR <- (ns + 1):(rm - ns)
#     M[sqR, sqC] <- mat
#     sqN <- -ns:ns
#     for(j in sqC)
#         for(i in sqR)
#             mat[i - ns, j - ns] <- mean(M[i + sqN, j + sqN], na.rm = TRUE)
#     mat[is.nan(mat)] <- NA
#     return(mat)
# }

smooth.matrix <- function(mat, ns){
    res <- cdt.matrix.mw(mat, ns, ns, mean, na.rm = TRUE)
    res[is.nan(res)] <- NA
    return(res)
}

# Calculate moving window values for the neighborhood of a center grid
cdt.matrix.mw <- function(x, sr, sc, fun, ...){
    fun <- match.fun(fun)
    nr <- nrow(x)
    nc <- ncol(x)
    res <- x * NA
    for(j in 1:nc){
        for(i in 1:nr){
            ir <- i + (-sr:sr)
            ir <- ir[ir > 0 & ir <= nr]
            ic <- j + (-sc:sc)
            ic <- ic[ic > 0 & ic <= nc]
            res[i, j] <- fun(c(x[ir, ic]), ...)
        }
    }
    return(res)
}

cdt.2matrices.mv <- function(x, y, sr, sc, fun, ...){
    stopifnot(dim(x) == dim(y))
    fun <- match.fun(fun)
    nr <- nrow(x)
    nc <- ncol(x)
    res <- x * NA
    for(j in 1:nc){
        for(i in 1:nr){
            ir <- i + (-sr:sr)
            ir <- ir[ir > 0 & ir <= nr]
            ic <- j + (-sc:sc)
            ic <- ic[ic > 0 & ic <= nc]
            res[i, j] <- fun(c(x[ir, ic]), c(y[ir, ic]), ...)
        }
    }
    return(res)
}

########################################

create_grid_buffer <- function(locations.stn, newgrid, radius, spheric)
{
    nx <- newgrid@grid@cells.dim[1]
    ny <- newgrid@grid@cells.dim[2]
    rx <- as.integer(radius/newgrid@grid@cellsize[1])
    ry <- as.integer(radius/newgrid@grid@cellsize[2])
    ix <- seq(1, newgrid@grid@cells.dim[1], rx)
    iy <- seq(1, newgrid@grid@cells.dim[2], ry)
    if(nx - ix[length(ix)] > rx/3) ix <- c(ix, nx)
    if(ny - iy[length(iy)] > ry/3) iy <- c(iy, ny)
    ixy <- expand.grid(ix, iy)
    ixy <- ixy[, 1] + ((ixy[, 2] - 1) * nx)
    coarsegrid <- as(newgrid[ixy, ], "SpatialPixels") 
    if(spheric){
        ctr <- rowSums(coarsegrid@bbox)/2
        pts <- c(ctr[1] + radius * cos(pi/4), ctr[2] + radius * sin(pi/4))
        radS <- distance.vector(pts, matrix(ctr, nrow = 1), spheric)
    }else radS <- radius

    dst <- distance.matrix(locations.stn@coords, coarsegrid@coords, spheric)
    dst <- rowSums(dst < 0.5 * radS) == 0
    coarsegrid <- coarsegrid[dst, ]

    buffer.out <- rgeos::gBuffer(locations.stn, width = 2 * radius)
    icoarse.out <- as.logical(over(coarsegrid, buffer.out))
    icoarse.out[is.na(icoarse.out)] <- FALSE
    coarsegrid <- coarsegrid[icoarse.out, ]

    buffer.grid <- rgeos::gBuffer(locations.stn, width = radius)
    igrid <- as.logical(over(newgrid, buffer.grid))
    igrid[is.na(igrid)] <- FALSE
    newdata0 <- newgrid[igrid, ]
    list(grid.buff = newdata0, ij = igrid, coarse = coarsegrid)
}

########################################

## Create parallel loop
cdt.doparallel <- function(condition, dopar = TRUE, detect.cores = TRUE, nb.cores = 2)
{
    okpar <- FALSE
    if(dopar){
        cores <- detectCores()
        if(detect.cores){
            nb.cores <- cores - 1
            okpar <- if(nb.cores >= 2) TRUE else FALSE
        }else{
            okpar <- if(cores >= 2 && nb.cores >= 2) TRUE else FALSE
        }
    }

    if(okpar & condition){
        klust <- makeCluster(nb.cores)
        registerDoParallel(klust)
        `%dofun%` <- `%dopar%`
        closeklust <- TRUE
    }else{
        klust <- NULL
        `%dofun%` <- `%do%`
        closeklust <- FALSE
    }

    list(dofun = `%dofun%`, cluster = klust, parLL = closeklust)
}

## foreach, use lapply if not parallel
utils::globalVariables(c('jloop'))
cdt.foreach <- function(loopL, parsL, ..., FUN)
{
    FUN <- match.fun(FUN)
    if(missing(parsL)) parsL <- list(condition = FALSE)
    is.parallel <- do.call(cdt.doparallel, parsL)

    if(is.parallel$parLL){
        `%parLoop%` <- is.parallel$dofun
        ret.loop <- foreach(jloop = loopL, ...) %parLoop% FUN(jloop)
        stopCluster(is.parallel$cluster)
    }else{
        .lapply <- function(X, FUN) lapply(X, FUN)
        ret.loop <- .lapply(loopL, FUN)
    }

    return(ret.loop)
}

########################################

#' Read CDT station data format
#' 
#' This function reads a CDT station data format.
#' 
#' @param stn.file full path name to the station file (CDT station data format)
#' @param sep the field separator character of the column, default is comma.
#' @param na missing value flag, the default is "-99".
#' @return A list of the id, longitude, latitude, elevation (if available) of the stations as a vector, date and the station data as a matrix of dimensions nrow= number of stations and ncol= length of date.
#' 
#' @export

read.CDTstation <- function(stn.file, sep = ",", na = "-99"){
    stnData <- read.table(stn.file, sep = sep, na.strings = na,
                          colClasses = "character",
                          stringsAsFactors = FALSE)
    seph <- rle(grepl('[[:digit:]]', as.character(stnData[, 1])))
    ipos <- which(!seph$values & seph$lengths >= 3 & seph$lengths <= 4)
    if(length(ipos) == 0 | ipos[1] != 1) stop('Station data is not in a standard unambiguous CDT format')
    pos <- seph$lengths[ipos[1]]

    list(id = as.character(stnData[1, -1]),
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
}

