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

## Define spatialPixels
defSpatialPixels <- function(grd_Coords, projCRS = CRS(as.character(NA)), regrid = FALSE)
{
    if(regrid){
        x <- grd_Coords$lon
        xrg <- diff(range(diff(x)))
        if(xrg > 0.0001){
            xr <- range(x)
            x <- seq(xr[1], xr[2], length.out = length(x))
        }
        y <- grd_Coords$lat
        yrg <- diff(range(diff(y)))
        if(yrg > 0.0001){
            yr <- range(y)
            y <- seq(yr[1], yr[2], length.out = length(y))
        }

        grd0 <- expand.grid(lon = x, lat = y)
        coordinates(grd0) <- ~lon+lat
        grd <- SpatialPixels(points = grd0, tolerance = 0.0002, proj4string = projCRS)
    }else{
        grd0 <- expand.grid(lon = grd_Coords$lon, lat = grd_Coords$lat)
        coordinates(grd0) <- ~lon+lat

        foo <- function(tol) SpatialPixels(points = grd0, tolerance = tol, proj4string = projCRS)
        grd <- try(foo(sqrt(sqrt(.Machine$double.eps))), silent = TRUE)
        if(inherits(grd, "try-error")) grd <- foo(0.005)
    }

    return(grd)
}

##############################################

## Get index of points at grid
grid2pointINDEX <- function(pts_Coords, grd_Coords, projCRS = CRS(as.character(NA)), regrid = FALSE)
{
    newgrid <- defSpatialPixels(grd_Coords, projCRS, regrid)
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

##############################################

## nx and ny for as.image
# x: diff(range( lon or lat ))
nx_ny_as.image <- function(x) round(x / (0.0167323 * x^0.9602))

## copy of fields::as.image
cdt.as.image <- function(pts.val, pts.xy, grid = NULL, nx = 64, ny = 64, weighted = FALSE, regrid = FALSE)
{
    if(is.null(grid)){
        xlim <- range(pts.xy[, 1], na.rm = TRUE)
        ylim <- range(pts.xy[, 2], na.rm = TRUE)
        xlim <- xlim + diff(xlim) * c(-1, 1) * 0.01
        ylim <- ylim + diff(ylim) * c(-1, 1) * 0.01
        grid <- list(lon = seq(xlim[1], xlim[2], length.out = nx),
                     lat = seq(ylim[1], ylim[2], length.out = ny))
    }
    xy <- do.call(expand.grid, grid)
    ijGrd <- grid2pointINDEX(list(lon = pts.xy[, 1], lat = pts.xy[, 2]), grid, regrid = regrid)
    out <- list(x = grid$lon, y = grid$lat, z = matrix(NA, length(grid$lon), length(grid$lat)))

    ij <- !is.na(pts.val)
    pts.val <- pts.val[ij]
    if(length(pts.val) == 0) return(out)
    pts.xy <- pts.xy[ij, , drop = FALSE]
    ijGrd <- ijGrd[ij]
    idx <- split(seq_along(ijGrd), ijGrd)

    if(any(sapply(idx, length) > 1)){
        w <- rep(1, length(ijGrd))
        if(weighted){
            idup <- duplicated(ijGrd) | duplicated(ijGrd, fromLast = TRUE)
            stn.grd <- xy[ijGrd, ]
            dist <- 1 / ((stn.grd[idup, 1] - pts.xy[idup, 1])^2 + (stn.grd[idup, 2] - pts.xy[idup, 2])^2)
            dist[is.infinite(dist)] <- 2 * max(dist[!is.infinite(dist)])
            w[idup] <- dist
        }
        val <- sapply(idx, function(j) sum(w[j] * pts.val[j]) / sum(w[j]))
    }else val <- pts.val[unlist(idx)]
    ij <- as.numeric(names(idx))
    out$z[ij] <- val
    return(out)
}

########################################

# create_grid_buffer <- function(locations.stn, newgrid, radius, spheric)
# {
#     nx <- newgrid@grid@cells.dim[1]
#     ny <- newgrid@grid@cells.dim[2]
#     rx <- as.integer(radius/newgrid@grid@cellsize[1])
#     ry <- as.integer(radius/newgrid@grid@cellsize[2])
#     ix <- seq(1, newgrid@grid@cells.dim[1], rx)
#     iy <- seq(1, newgrid@grid@cells.dim[2], ry)
#     if(nx - ix[length(ix)] > rx/3) ix <- c(ix, nx)
#     if(ny - iy[length(iy)] > ry/3) iy <- c(iy, ny)
#     ixy <- expand.grid(ix, iy)
#     ixy <- ixy[, 1] + ((ixy[, 2] - 1) * nx)
#     coarsegrid <- as(newgrid[ixy, ], "SpatialPixels") 
#     if(spheric){
#         ctr <- rowSums(coarsegrid@bbox)/2
#         pts <- c(ctr[1] + radius * cos(pi/4), ctr[2] + radius * sin(pi/4))
#         radS <- distance.vector(pts, matrix(ctr, nrow = 1), spheric)
#     }else radS <- radius

#     dst <- distance.matrix(locations.stn@coords, coarsegrid@coords, spheric)
#     dst <- rowSums(dst < 0.5 * radS) == 0
#     coarsegrid <- coarsegrid[dst, ]

#     buffer.out <- rgeos::gBuffer(locations.stn, width = 2 * radius)
#     icoarse.out <- as.logical(over(coarsegrid, buffer.out))
#     icoarse.out[is.na(icoarse.out)] <- FALSE
#     coarsegrid <- coarsegrid[icoarse.out, ]

#     buffer.grid <- rgeos::gBuffer(locations.stn, width = radius)
#     igrid <- as.logical(over(newgrid, buffer.grid))
#     igrid[is.na(igrid)] <- FALSE
#     newdata0 <- newgrid[igrid, ]
#     list(grid.buff = newdata0, ij = igrid, coarse = coarsegrid)
# }

create_grid_buffer <- function(locations.stn, newgrid, radius, fac = 4, spheric = FALSE)
{
    nx <- newgrid@grid@cells.dim[1]
    ny <- newgrid@grid@cells.dim[2]

    rx <- as.integer(radius / (fac * newgrid@grid@cellsize[1]))
    ry <- as.integer(radius / (fac * newgrid@grid@cellsize[2]))

    ix <- seq(1, nx, rx)
    iy <- seq(1, ny, ry)
    if(nx - ix[length(ix)] > rx/3) ix <- c(ix, nx)
    if(ny - iy[length(iy)] > ry/3) iy <- c(iy, ny)
    ixy <- expand.grid(ix, iy)
    ixy <- ixy[, 1] + ((ixy[, 2] - 1) * nx)
    coarsegrid <- as(newgrid[ixy, ], "SpatialPixels")

    if(spheric){
        ctr <- rowSums(coarsegrid@bbox)/2
        pts <- c(ctr[1] + radius * cos(pi/4), ctr[2] + radius * sin(pi/4))
        radius <- distance.vector(pts, matrix(ctr, nrow = 1), spheric)
    }

    xgrd <- apply(coarsegrid@coords, 2, unique)
    loc.stn <- cdt.as.image(locations.stn$stn, locations.stn@coords, xgrd, regrid = TRUE)
    loc.stn <- cbind(do.call(expand.grid, loc.stn[c('x', 'y')]), z = c(loc.stn$z))
    loc.stn <- loc.stn[!is.na(loc.stn$z), , drop = FALSE]
    coordinates(loc.stn) <- c('x', 'y')

    dst <- fields::rdist(locations.stn@coords, coarsegrid@coords)
    dst <- colSums(dst < 0.5 * radius) == 0
    coarsegrid <- coarsegrid[dst, ]

    buffer.out <- rgeos::gBuffer(loc.stn, width = 0.75 * radius)
    icoarse.out <- as.logical(over(coarsegrid, buffer.out))
    icoarse.out[is.na(icoarse.out)] <- FALSE
    coarsegrid <- coarsegrid[icoarse.out, ]

    buffer.grid <- rgeos::gBuffer(loc.stn, width = 0.75 * radius)
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

