args <- commandArgs(trailingOnly = TRUE)                                           
if (length(args) != 2)
{
  cat ("Rscript R/find-il-position.R IL1-1_SNP\n")
  quit("no")
}
#args<- c("email/from/sara-palmer/090712/ua159csp-f-chris.csv", "email/to/sara-palmer/071412/smu86csp.csv")

x <- read.csv(args[1],head=T)
d <- read.csv(args[2],head=T)

# SMU_1056's fake annotation is added for the sake of running the script without
# error. Sara should be able to let me know the annotation.
stopifnot(sum(d$id %in% x[,1]) == nrow(d))

# Select 2nd column (OralGen), 5th column (gene name), 6th column function .
d <- data.frame(d[,-1], oralgen.annotation=paste("=HYPERLINK(\"",x[match(d$id,x[,1]),4],"\",\"",x[match(d$id,x[,1]),2],"\")",sep=""),
                        genename=x[match(d$id,x[,1]),5],
                        some.annotation=x[match(d$id,x[,1]),6])
                   
write.csv(d,file=paste(sub(".csv","",args[2]),"-annotation.csv",sep=""),row.names=FALSE)

