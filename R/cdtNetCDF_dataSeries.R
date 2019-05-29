
ncFilesInfo <- function(Tstep, start.date, end.date, months,
						ncDir, ncFileFormat, error.msg)
{
	if(Tstep == 'daily'){
		dates <- format(seq(start.date, end.date, 'day'), '%Y%m%d')
		ncDataFiles <- file.path(ncDir, sprintf(ncFileFormat, substr(dates, 1, 4),
										substr(dates, 5, 6), substr(dates, 7, 8)))
	}
	if(Tstep == 'pentad'){
		dates <- seq(start.date,  end.date, 'day')
		dates <- paste0(format(dates[which(as.numeric(format(dates, '%d')) <= 6)], '%Y%m'),
					as.numeric(format(dates[which(as.numeric(format(dates, '%d')) <= 6)], '%d')))
		ncDataFiles <- file.path(ncDir, sprintf(ncFileFormat, substr(dates, 1, 4),
										substr(dates, 5, 6), substr(dates, 7, 7)))
	}
	if(Tstep == 'dekadal'){
		dates <- seq(start.date,  end.date, 'day')
		dates <- paste0(format(dates[which(as.numeric(format(dates, '%d')) <= 3)], '%Y%m'),
					as.numeric(format(dates[which(as.numeric(format(dates, '%d')) <= 3)], '%d')))
		ncDataFiles <- file.path(ncDir, sprintf(ncFileFormat, substr(dates, 1, 4),
										substr(dates, 5, 6), substr(dates, 7, 7)))
	}
	if(Tstep == 'monthly'){
		dates <- format(seq(start.date, end.date, 'month'), '%Y%m')
		ncDataFiles <- file.path(ncDir, sprintf(ncFileFormat, substr(dates, 1, 4),
												substr(dates, 5, 6)))
	}

	months.dates <- as(substr(dates, 5, 6), 'numeric')
	imo <- months.dates %in% months
	dates <- dates[imo]
	ncDataFiles <- ncDataFiles[imo]

	existFl <- unlist(lapply(ncDataFiles, file.exists))
	if(!any(existFl)){
		cat(error.msg, "\n")
		return(NULL)
	}
	return(list(dates = dates, nc.files = ncDataFiles, exist = existFl))
}

transposeNCDFData <- function(x, ixy){
	if(ixy$ilon < ixy$ilat)
		x[ixy$olon, ixy$olat]
	else
		t(x)[ixy$olon, ixy$olat]
}
