### =========================================================================
### Utilities to generate substitution matrices
### -------------------------------------------------------------------------
###


nucleotideSubstitutionMatrix <-
function(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA", symmetric = TRUE) {
  "%safemult%" <- function(x, y) ifelse(is.infinite(x) & y == 0, 0, x * y)
  type <- match.arg(type, c("DNA", "RNA"))
  if (!isSingleNumber(match) || !isSingleNumber(mismatch)) {
    stop("'match' and 'mismatch' must be non-missing numbers")
  }
  if (baseOnly) {
    letters <- IUPAC_CODE_MAP[DNA_BASES]
  }
  else {
    letters <- IUPAC_CODE_MAP
  }
  if (type == "RNA") {
    names(letters) <- chartr("T", "U", names(letters))
  }
  nLetters <- length(letters)
  splitLetters <- strsplit(letters,split="")
  submat <- matrix(0, nrow = nLetters, ncol = nLetters, dimnames = list(names(letters), names(letters)))
  if (symmetric) {
    for (i in seq_len(nLetters)) {
      for (j in i:nLetters) {
        submat[i,j] <- submat[j,i] <- mean(outer(splitLetters[[i]], splitLetters[[j]], "=="))
      }
    }
  }
  else {
    for (i in seq_len(nLetters)) {
      for (j in i:nLetters) {
        submat[i,j] <- mean(outer(splitLetters[[i]], splitLetters[[j]], "%in%"))
        submat[j,i] <- mean(outer(splitLetters[[j]], splitLetters[[i]], "%in%"))
      }
    }
  }
  abs(match) * submat - abs(mismatch) %safemult% (1 - submat)
}

errorSubstitutionMatrices <-
function(errorProbability, fuzzyMatch = c(0, 1), alphabetLength = 4L, bitScale = 1) {
  if (!is.numeric(errorProbability) || !all(!is.na(errorProbability) & errorProbability >= 0 & errorProbability <= 1))
    stop("'errorProbability' must be a numeric vector with values between 0 and 1 inclusive")
  if (!is.numeric(fuzzyMatch) || !all(!is.na(fuzzyMatch) & fuzzyMatch >= 0 & fuzzyMatch <= 1))
    stop("'fuzzyMatch' must be a numeric vector with values between 0 and 1 inclusive")
  errorMatrix <-
    outer(errorProbability, errorProbability,
          function(e1,e2,n) e1 + e2 - (n/(n - 1)) * e1 * e2,
          n = alphabetLength)
  adjMatchProbs <-
    lapply(list(match = (1 - errorMatrix) * alphabetLength,
                mismatch = errorMatrix * (alphabetLength / (alphabetLength - 1))),
                function(x) {dimnames(x) <- list(names(errorProbability), names(errorProbability)); x})
  output <-
    array(NA_real_, dim = c(length(errorProbability), length(errorProbability), length(fuzzyMatch)),
          dimnames = list(names(errorProbability), names(errorProbability), as.character(fuzzyMatch)))
  for (i in seq_len(length(fuzzyMatch))) {
    output[,,i] <-
      bitScale *
        log2(fuzzyMatch[i] * adjMatchProbs[["match"]] + (1 - fuzzyMatch[i]) * adjMatchProbs[["mismatch"]])
  }
  output
}

qualitySubstitutionMatrices <-
function(fuzzyMatch = c(0, 1), alphabetLength = 4L, qualityClass = "PhredQuality", bitScale = 1) {
  if (!is.numeric(fuzzyMatch) || !all(!is.na(fuzzyMatch) & fuzzyMatch >= 0 & fuzzyMatch <= 1))
    stop("'fuzzyMatch' must be a numeric vector with values between 0 and 1 inclusive")
  if (!is(new(qualityClass), "XStringQuality"))
    stop("'qualityClass' must be one of the 'XStringQuality' classes")
  qualityIntegers <- minQuality(new(qualityClass)):maxQuality(new(qualityClass))
  errorProbability <- qualityConverter(qualityIntegers, qualityClass, "numeric")
  names(errorProbability) <- as.character(qualityIntegers)
  errorSubstitutionMatrices(errorProbability, fuzzyMatch, alphabetLength = alphabetLength, bitScale = bitScale)
}

