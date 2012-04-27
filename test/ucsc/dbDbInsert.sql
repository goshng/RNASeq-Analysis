INSERT INTO dbDb
    (name, description, nibPath, organism,
     defaultPos, active, orderKey, genome, scientificName,
     htmlPath, hgNearOk, hgPbOk, sourceName)
VALUES
    ("Seeq4047", "Sep 23 2011", "/gbdb/Seeq4047", "S. equi subsp. equi 4047",
     "chr1:1-100000", 1, 1, "S. equi subsp. equi 4047", "Streptococcus equi equi",
     "/gbdb/Seeq4047/html/description.html", 0, 0, "Streptococcus recombination study version 1.0");
INSERT INTO defaultDb (genome, name) VALUES ("S. equi subsp. equi 4047", "Seeq4047");
INSERT INTO genomeClade (genome, clade, priority) VALUES ("S. equi subsp. equi 4047", "streptococcus", 9);
/* INSERT INTO clade (name, label, priority) VALUES ("streptococcus", "Streptococcus", 45); */



/* SdeqATCC12394 S. dys. equi. ATCC12394 
   Streptococcus dysgalactiae equisimilis
 */
/* SdeqGGS124 S. dys. equi. GGS124
   Streptococcus dysgalactiae equisimilis
 */

/* sdd? */

/* SpyMGAS315 S. pyogenes MGAS315
   Streptococcus pyogenes
 */
/* SpyMGAS10750 S. pyogenes MGAS10750
   Streptococcus pyogenes
 */
/* Seeq4047 S. equi subsp. equi 4047
   Streptococcus equi equi
 */
