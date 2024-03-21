.onLoad <- function(libname, pkgname)
{
}

.onUnload <- function(libpath)
{
    library.dynam.unload("pwalign", libpath)
}

.test <- function() BiocGenerics:::testPackage("pwalign")

