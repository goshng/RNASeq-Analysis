sa <- data.frame(condition=factor(c(rep("treated",3),rep("untreated",3))),replicate=c(1,2,3,1,2,3),type=c("single-read","paired-end","paired-end","single-read","single-read","paired-end")) 
rownames(sa) <- c("treated1fb","treated2fb","treated3fb","untreated1fb","untreated2fb","untreated3fb")
countfiles = file.path(inDir, paste(rownames(sa), "txt", sep="."))
annotationfile = file.path(inDir, "Dmel.BDGP5.25.62.DEXSeq.chr.gff")
flattenedfile = annotationfile 

lf <- lapply( countfiles, function(x) read.table( x, header=FALSE,stringsAsFactors=FALSE ) )
dcounts <- sapply( lf, `[[`, "V2" )
rownames(dcounts) <- lf[[1]][,1]
dcounts <- dcounts[ substr(rownames(dcounts),1,1)!="_", ]
genesrle <- sapply( strsplit( rownames(dcounts), ":"), "[[", 1)
index <- order(genesrle)
dcounts <- dcounts[index,]
genesrle <- genesrle[index]
exons <- paste("E", sapply(strsplit(rownames(dcounts), ":"), function(x){return(x[2])}), sep="")
colnames(dcounts) <- countfiles
   get_attribute <- function(x, attribute){
      x <- gsub("\"", "", x)
      x <- gsub("=", " ", x)
      x <- strsplit(x, "; ")
      sapply(1:length(x), function(i){
        sapply(strsplit(x[[i]][grepl(attribute, x[[i]])], " "), "[[", 2)
      })  
   }  

aggregates<-read.delim(flattenedfile, stringsAsFactors=FALSE, header=FALSE)
colnames(aggregates)<-c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
aggregates$gene_id <- get_attribute(aggregates$attr, attribute="gene_id")
aggregates<-aggregates[order(aggregates$gene_id),]
transcripts <- get_attribute(aggregates$attr, attribute="transcripts")
transcripts <- gsub("\\+", ";", transcripts)
exoninfo<-data.frame(chr=aggregates$chr, start=aggregates$start, end=aggregates$end, strand=aggregates$strand)
if (!all(aggregates$gene_id == genesrle)){
  stop("Count files do not correspond to the flattened annotation file")
}
ecs<-newExonCountSet(countData=dcounts, design=sa, geneIDs=genesrle, exonIDs=exons, exonIntervals=exoninfo, transcripts=transcripts)
ecs@annotationFile <- flattenedfile
pData(ecs)$countfiles <- countfiles
# FINISHED
return(ecs)
