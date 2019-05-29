
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

distance.vector <- function(xy.pt, XY, spheric){
	# xy.pt: vector c(x, y)
	# XY: matrix X Y
	xy.pt <- as.numeric(xy.pt)
	XY <- as.matrix(XY)
	nr <- nrow(XY)
	out <- vector("numeric", length = nr)
	out <- .Fortran("distance_vector", xy.pt, XY, as.integer(nr), as.integer(spheric), dst = out)
	out$dst
}

distance.matrix <- function(XY1, XY2, spheric){
	# XY1: matrix X Y
	# XY2: matrix X Y
	XY1 <- as.matrix(XY1)
	XY2 <- as.matrix(XY2)
	n1 <- nrow(XY1)
	n2 <- nrow(XY2)
	dst <- matrix(double(1), nrow = n2, ncol = n1)
	out <- .Fortran("distance_matrix", XY1, XY2, as.integer(n1), as.integer(n2),
					as.integer(spheric), dst = dst)
	out$dst
}

