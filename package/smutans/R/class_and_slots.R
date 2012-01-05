setClass( "Smutans", 
  representation = representation( 
    countDataSet = "CountDataSet", 
    cds = "ANY", 
    res = "list", 
    padj = "numeric",
    pval = "numeric",
    title = "character"
  )
)

newSmutans <- function( countData, title = "default" )
{
  cds <- new( "Smutans", 
              countDataSet = countData,
              cds = countData,
              title = title )
  cds 
}

