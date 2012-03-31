y <- rep("",times=length(genesrle)) 
for (i in unique(genesrle)) {
  x <- c()
  for (j in 1:sum(genesrle == i)) {
    x <- c(x, sprintf("E%03d", j))
  }
  y[which(genesrle == i)] <- x
}
