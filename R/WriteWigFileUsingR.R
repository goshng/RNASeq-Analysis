##
## Writes a wiggle file using R.
##
require(GROseq)
load("../SOAP.0-160.combined.noRNAreps.RData")
nr0m <- NROW(E20m)
nr10m <- NROW(E210m)
nr40m <- NROW(E240m)
nr160m <- NROW(E2160m)

expCounts <- mean(nr0m, nr10m, nr40m, nr160m)

WriteWiggle(p= E20m[,c(1:3,6)], file="E20m_Plus", size=5, reverse=FALSE, str="+", debug=FALSE)
WriteWiggle(p= E20m[,c(1:3,6)], file="E20m_Minus", size=5, reverse=TRUE, str="-", debug=FALSE)

WriteWiggle(p= E210m[,c(1:3,6)], file="E210m_Plus", size=5, reverse=FALSE, str="+", debug=FALSE)
WriteWiggle(p= E210m[,c(1:3,6)], file="E210m_Minus", size=5, reverse=TRUE, str="-", debug=FALSE)

WriteWiggle(p= E240m[,c(1:3,6)], file="E240m_Plus", size=5, reverse=FALSE, str="+", debug=FALSE)
WriteWiggle(p= E240m[,c(1:3,6)], file="E240m_Minus", size=5, reverse=TRUE, str="-", debug=FALSE)

WriteWiggle(p= E2160m[,c(1:3,6)], file="E2160m_Plus", size=5, reverse=FALSE, str="+", debug=FALSE)
WriteWiggle(p= E2160m[,c(1:3,6)], file="E2160m_Minus", size=5, reverse=TRUE, str="-", debug=FALSE)

