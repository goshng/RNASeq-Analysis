library(rtracklayer)

#args<- c("output/sara/gff/44.gff", "email/to/sara-palmer/080812/UA159-11VS1-SMU44.csv")
args <- commandArgs(trailingOnly = TRUE)                                           
if (length(args) != 2)
{
  cat ("Rscript R/sara-create-chris-gff.R output/sara/gff/44.gff email/to/sara-palmer/080812/UA159-11VS1-SMU44.csv\n")
  quit("no")
}

gff.range <- import.gff3(args[1])
gff.range.cds <- gff.range[gff.range$type=="CDS",]
x <- gff.range.cds

# x <- read.csv(args[1],head=T)
d <- read.csv(args[2],head=T)

# SMU_1056's fake annotation is added for the sake of running the script without
# error. Sara should be able to let me know the annotation.
stopifnot(sum(d$id %in% gff.range.cds$locus_tag) == nrow(d))

# Select 2nd column (OralGen), 5th column (gene name), 6th column function .

x.note <- unlist(lapply(x$Note,function(x) if(length(x)==0) NA else x))
d <- data.frame(d[,-1],
                note.annotation=x.note[match(d$id,x$locus_tag)],
                product.annotation=x$product[match(d$id,x$locus_tag)])

write.csv(d,file=paste(sub(".csv","",args[2]),"-annotation.csv",sep=""),row.names=FALSE)
