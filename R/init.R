.First.lib <- function(libname, pkgname) {

  library.dynam("preprocessCore",pkgname,libname,now=FALSE)

  .C("R_init_preprocessCore",PACKAGE="preprocessCore")
  
}
