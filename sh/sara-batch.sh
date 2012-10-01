###############################################################################
# Copyright (C) 2012 Sang Chul Choi
#
# This file is part of RNASeq Analysis.
# 
# RNASeq Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RNASeq Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RNASeq Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

function sara-batch {
#  rm -f email/to/sara-palmer/080812/*-annotation.csv;
#  for i in email/to/sara-palmer/080812/*.csv; do
#    Rscript R/sara-add-annotation.R \
#      email/from/sara-palmer/090712/ua159csp-f-chris.csv $i
#  done

#  rm -rf output/sara/gff; mkdir output/sara/gff
#  for i in 44 20 21 56 69 57 86; do
#    cat data/denovo/SMU${i}_allcontigs.gbf \
#      | perl pl/bp_genbank2gff3.pl -in stdin -out stdout \
#      | perl -e "while (<>) { if (/^>/) { last; } print; }" \
#      > output/sara/gff/$i.gff
#  done
  Rscript R/sara-create-chris-gff.R \
    output/sara/gff/44.gff email/to/sara-palmer/080812/UA159-11VS1-SMU44.csv
  Rscript R/sara-create-chris-gff.R \
    output/sara/gff/20.gff email/to/sara-palmer/080812/UA159-15JP3-SMU20.csv
  Rscript R/sara-create-chris-gff.R \
    output/sara/gff/21.gff email/to/sara-palmer/080812/UA159-1SM1-SMU21.csv
  Rscript R/sara-create-chris-gff.R \
    output/sara/gff/56.gff email/to/sara-palmer/080812/UA159-N29-SMU56.csv
  Rscript R/sara-create-chris-gff.R \
    output/sara/gff/69.gff email/to/sara-palmer/080812/UA159-NLML4-SMU69.csv
  Rscript R/sara-create-chris-gff.R \
    output/sara/gff/57.gff email/to/sara-palmer/080812/UA159-NMT4863-SMU57.csv
  Rscript R/sara-create-chris-gff.R \
    output/sara/gff/21.gff email/to/sara-palmer/080812/smu21phSMU21.csv
  Rscript R/sara-create-chris-gff.R \
    output/sara/gff/86.gff email/to/sara-palmer/080812/smu86cspSMU86.csv

  Rscript R/sara-nine-pattern.R \
    email/to/sara-palmer/092412/ua159ph-annotation.csv \
    email/to/sara-palmer/092412/smu21ph-annotation.csv
  Rscript R/sara-nine-pattern.R \
    email/to/sara-palmer/092412/ua159csp-annotation.csv \
    email/to/sara-palmer/092412/smu86csp-annotation.csv

}
