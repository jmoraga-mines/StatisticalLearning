#### Utility functions, used for DOE project #
####

if (!require("pacman")) {
  install.packages("pacman")
  require("pacman")
}
pacman::p_load(raster, rgdal)

extent_hymap <-
  raster(
    xmn = 325209,
    xmx = 336531,
    ymn = 4395438,
    ymx = 4412103,
    res = c(3, 3),
    crs = "+init=epsg:32611"
  )

#' Save a raster object to a standard format (GRD) and creates
#' an ENVI-style header files to preserve raster layer(s) name(s).
#' This is a wrapper for the raster::writeRaster function
#'
#' @param x Raster. A raster or stack object
#' @param filename Character. A file name (with no extension or .grd extension)
#' @param format Character. a format from the raster::writeFormats list, default
#' is "raster"
#' @param overwrite Boolean. either TRUE or FALSE
#' @param bandorder Character. 'BIL', 'BIP', or 'BSQ'. For 'native' file formats only.
#'
#' @return A file handle to the written file
#'
#' @examples
#' doe_write_raster(raster01, "my_doe_file")
#' doe_write_raster(raster01, "my_doe_file", overwrite = FALSE)
#' 
doe_write_raster <-
  function(x,
           filename,
           format = "raster",
           overwrite = TRUE,
           bandorder = "BSQ") {
    if (tools::file_ext(filename) != "grd") {
      filename <- tools::file_path_sans_ext(filename)
      filename <- paste(filename, ".grd", sep = "")
    }
    f1 <- raster::writeRaster(
      x = x,
      filename = filename,
#      bandorder = bandorder,
      format = format,
      overwrite = overwrite
    )
    hdr(f1, "ENVI")
    return(f1)
  }
