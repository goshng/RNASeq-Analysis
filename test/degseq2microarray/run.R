
pdf("wt2tw.pdf")
x <-read.table("wt2tw.txt")
plot(x$V2~x$V3,xlab="Log2 fold changes in microarray",ylab="Log2 fold change in RNA-Seq",main="U159 versus TW1")
abline(0,1)
dev.off()

pdf("tw.pdf")
y <-read.table("tw.txt")
plot(y$V2~y$V3,xlab="Log2 fold changes in microarray",ylab="Log2 fold change in RNA-Seq",main="TW1-glucose versus TW1-galactose")
abline(0,1)
dev.off()
