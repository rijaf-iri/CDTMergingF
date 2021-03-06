#' Parse the NetCDF files
#'
#' Parse all NetCDF data inside a directory according the provided filename format
#' 
#' @param Tstep The time step of the NetCDF files. Valid options: 'daily', 'pentad', 'dekadal', 'monthly'.
#' @param start.date An object of class \code{Date} indicating the start date to be parsed.
#' @param end.date An object of class \code{Date} indicating the end date to be parsed.
#' @param months An integer vector indicating the months to be parsed.
#' @param ncDir The full path to the directory containing the NetCDF files.
#' @param ncFileFormat Filename format of the NetCDF files.
#' @param error.msg A message to be displayed when no files were found.
#' @return A list of a character vector indicating the date, a character vector containing the full path name of the NetCDF files and a logical vector indicating if the file exists.
#' 
#' @export

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
