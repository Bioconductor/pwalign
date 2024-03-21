### =========================================================================
### Percent Sequence Identity
### -------------------------------------------------------------------------
###


setGeneric("pid", signature="x", 
    function(x, type="PID1") standardGeneric("pid")
)

setMethod("pid", "PairwiseAlignments",
    function(x, type="PID1") {
        type <- match.arg(type, c("PID1", "PID2", "PID3", "PID4"))
        denom <- switch(type,
            PID1 = nchar(x),
            PID2 = nmatch(x) + nmismatch(x),
            PID3 = pmin(nchar(unaligned(pattern(x))),
                        nchar(unaligned(subject(x)))),
            PID4 = (nchar(unaligned(pattern(x))) +
                    nchar(unaligned(subject(x)))) / 2
        )
        100 * nmatch(x)/denom
})

