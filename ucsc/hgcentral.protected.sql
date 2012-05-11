INSERT INTO `dbDb` VALUES 
("SmuUA159v2", "March 2012", "/gbdb_cornell/SmuUA159v2", "S. mutans UA159", "chr1:1-10000", 1, 202, "S. mutans UA159", "Streptococcus mutans UA159", "/gbdb_cornell/SmuUA159v2/html/description.html", 0, 0, "genbank", 0),
("AEXT01", "Aug 24 2011", "/gbdb_cornell/AEXT01", "S. agalactiae FSL S3-026", "chr1:1-10000", 1, 1, "S. agalactiae FSL S3-026", "Streptococcus agalactiae FSL S3-026", "/gbdb_cornell/AEXT01/html/description.html", 0, 0, "genbank"),
("NC_004368", "June 2010", "/gbdb_cornell/NC_004368", "S. agalactiae NEM316", "chr1:1-10000", 1, 1, "S. agalactiae NEM316", "Streptococcus agalactiae NEM316", "/gbdb_cornell/NC_004368/html/description.html", 0, 0, "genbank"),
("NC_007432", "March 2010", "/gbdb_cornell/NC_007432", "S. agalactiae A909", "chr1:1-10000", 1, 1, "S. agalactiae A909", "Streptococcus agalactiae A909", "/gbdb_cornell/NC_007432/html/description.html", 0, 0, "genbank");


INSERT INTO `defaultDb` VALUES 
        ("S. agalactiae FSL S3-026", "AEXT01"),
        ("S. agalactiae NEM316", "NC_004368"),
        ("S. agalactiae A909", "NC_007432");

INSERT INTO `genomeClade` VALUES 
        ("S. agalactiae FSL S3-026", "streptococcus", 22),
        ("S. agalactiae NEM316", "streptococcus", 23),
        ("S. agalactiae A909", "streptococcus", 24);
