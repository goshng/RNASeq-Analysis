setClass( "Smutans", 
  representation = representation( 
    countDataSet = "CountDataSet", 
    cds = "ANY", 
    res = "list", 
    padj = "numeric",
    pval = "numeric"
    )
)

newSmutans <- function( countData )
{
  cds <- new( "Smutans", 
              countDataSet = countData,
              cds = countData )
  cds 
}

