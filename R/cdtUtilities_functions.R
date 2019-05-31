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
	newgrid <- SpatialPixels(points = newgrid,
							tolerance = sqrt(sqrt(.Machine$double.eps)),
							proj4string = CRS(as.character(NA)))
	pts.loc <- data.frame(lon = pts_Coords$lon, lat = pts_Coords$lat)
	pts.loc <- SpatialPoints(pts.loc)
	ijGrd <- unname(over(pts.loc, geometry(newgrid)))
	return(ijGrd)
}

########################################

smooth.matrix <- function(mat, ns){
	mat0 <- mat
	M <- matrix(NA, nrow(mat) + 2 * ns, ncol(mat) + 2 * ns)
	sqC <- (ns + 1):(ncol(M) - ns)
	sqR <- (ns + 1):(nrow(M) - ns)
	M[sqR, sqC] <- mat
	sqN <- -ns:ns
	for(j in sqC)
		for(i in sqR)
			mat[i - ns, j - ns] <- mean(M[i + sqN, j + sqN], na.rm = TRUE)
	mat[is.nan(mat)] <- NA
	mat <- (2 * mat0 + mat) / 3
	return(mat)
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
cdt.foreach <- function(loopL, parsL = NULL, ..., FUN)
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




