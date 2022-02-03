#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
.onAttach <- function(libname, pkgname){ 
  msg <- paste0("Please cite ", "ASURAT", " using the following:\n", 
                "  Iida, K., Kondo, J., Wibisana, J.N., Inoue M., Okada M.,
  bioRxiv (2021), https://doi.org/10.1101/2021.06.09.447731")

  packageStartupMessage(msg)
}
