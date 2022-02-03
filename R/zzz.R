.onAttach <- function(libname, pkgname){ 
  version_msg <- paste0(pkgname, " v", version, "\n\n")

  msg <- paste0("Thank you for using ", pkgname, "\n\n", 
                "Please cite ", pkgname, " using the following:\n", 
                "Iida, K., Kondo, J., Wibisana, J.N., Inoue M., Okada M., 2021. BiorXiv.
  ASURAT: functional annotation-driven unsupervised clustering of single-cell transcriptomes.\nhttps://doi.org/10.1101/2021.06.09.447731
  ")

  msg2 <- "\nPlease refer to https://keita-iida.github.io/ASURAT/ for package usage"

  packageStartupMessage(paste0(version_msg, msg, msg2))
}
