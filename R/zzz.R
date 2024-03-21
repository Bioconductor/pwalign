.onLoad <- function(libname, pkgname)
{
}

.onUnload <- function(libpath)
{
    library.dynam.unload("pairwiseAlignment", libpath)
}

.test <- function() BiocGenerics:::testPackage("pairwiseAlignment")

