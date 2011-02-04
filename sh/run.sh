#!/bin/bash
# Author: Sang Chul Choi

# READID=SRR031222
# READID=SRR031212
# READID=SRR031130
# READID=test-illumina

PRINSEQ=/Users/goshng/bin/prinseq-lite.pl
GENOMEDATADIR=/Volumes/Elements/Documents/Projects/mauve/bacteria
READDATADIR=/Volumes/Elements/Documents/Projects/rnaseq/data

# HpyloriGENOME=$GENOMEDATADIR/Helicobacter_pylori_26695_uid57787/NC_000915.fna
# HyploriREADDIR=/Volumes/Elements/Documents/Projects/rnaseq/data/SRP001481

HpyloriQUERY=output/$READID.fasta
HpyloriMAP=output/$READID.map

STUDYID=SRP001481
# Directories
function prepare-filesystem {
  GENOMEACCESSION=NC_000915
  READACCESSION=SRP001481
  RNASEQANALYSISDIR=`pwd`
  BASEDIR=`pwd`/output/$STUDYID

  RUNMAPDIR=$BASEDIR/run-map

  CACBASEDIR=/Volumes/sc2265/Documents/Projects/rnaseq/output/$STUDYID
  CACDATADIR=$CACBASEDIR/data
  CACRUNMAPDIR=$CACBASEDIR/run-map
  CACRUNGENOMEDIR=$CACBASEDIR/run-genome

  BATCH_SH_RUN_MAP=$RUNMAPDIR/batch.sh
  TMPDIR=/tmp/$JOBID.scheduler.v4linux
  TMPINPUTDIR=$TMPDIR/input
  SWIFTGENDIR=choi@swiftgen:Documents/Projects/map/output/$STUDYID
  SWIFTGENRUNMAP=$SWIFTGENDIR/run-map
  RUNLOG=$BASEDIR/run.log


  RUNGENOMEDIR=$BASEDIR/run-genome
  RUNGENOMEOUTPUTDIR=$RUNGENOMEDIR/output
  REFGENOMEINDEX=$RUNGENOMEOUTPUTDIR/$GENOMEACCESSION.idx

  # This must be changed when adding more analyses.
  # Reference genome and its query reads must be chosen by users.
  REFGENOMEFASTA=output/$GENOMEACCESSION.fna # must be from a directory
  READDATA=$READDATADIR/$READACCESSION # must be from some selectoin
}

function mkdir-STUDYID {
  mkdir $BASEDIR
  mkdir $CACBASEDIR
  mkdir $CACDATADIR
  mkdir $CACRUNMAPDIR
  mkdir $RUNMAPDIR
  mkdir $RUNGENOMEDIR
  mkdir $RUNGENOMEOUTPUTDIR
}

function copy-batch-sh-run-map {
cat>$BATCH_SH_RUN_MAP<<EOF
#!/bin/bash
#PBS -l walltime=59:00,nodes=1
#PBS -A acs4_0001
#PBS -j oe
#PBS -N Strep-${STUDYID}-Map
#PBS -q v4dev
##PBS -m e
##PBS -M schoi@cornell.edu
WORKDIR=\$PBS_O_WORKDIR
DATADIR=\$WORKDIR/../data
FASTQDUMP=\$HOME/usr/bin/fastq-dump
PRINSEQ=\$HOME/usr/bin/prinseq-lite.pl
SEGEMEHL=\$HOME/usr/bin/segemehl.x

OUTPUTDIR=\$TMPDIR/output
INPUTDIR=\$TMPDIR/input
mkdir \$INPUTDIR
mkdir \$OUTPUTDIR
cp \$FASTQDUMP \$TMPDIR/
cp \$PRINSEQ \$TMPDIR/
cp \$SEGEMEHL \$TMPDIR/
cp -r \$DATADIR \$INPUTDIR/
cd \$TMPDIR

#    convert-sra-to-fastq $s
#    prinseq-a-read $s
#    match-a-read $s
# remove-garbage $s

SRP001481SRA=( SRR031126 SRR031127 SRR031128 SRR031129 SRR031130 SRR031131 SRR031142 SRR031143 SRR031144 SRR031145 SRR031146 SRR031178 SRR031179 SRR031180 SRR031181 SRR031182 SRR031183 SRR031184 SRR031185 SRR031186 SRR031187 SRR031188 SRR031189 SRR031190 SRR031191 SRR031192 SRR031193 SRR031194 SRR031195 SRR031196 SRR031197 SRR031198 SRR031199 SRR031200 SRR031201 SRR031202 SRR031203 SRR031204 SRR031205 SRR031206 SRR031207 SRR031208 SRR031209 SRR031210 SRR031211 SRR031212 SRR031213 SRR031214 SRR031215 SRR031216 SRR031217 SRR031218 SRR031219 SRR031220 SRR031221 SRR031222 SRR031223 SRR031224 SRR031225 SRR031226 SRR031276 SRR031277 SRR031278 SRR031279 SRR031280 SRR031281 SRR031282 SRR031283 SRR031284 SRR031285 SRR031286 SRR031287 SRR031288 SRR031289 SRR031290 SRR031291 SRR031292 SRR031293 SRR031294 SRR031295 SRR034162 SRR034163 SRR034165 SRR034166 SRR034167 SRR034168 SRR034169 SRR034170 SRR034171 SRR034172 SRR034173 SRR034174 SRR034175 SRR034176 SRR034177 SRR034178 SRR034179 SRR034180 SRR034181 SRR034182 SRR034183 SRR034185 SRR034186 SRR034187 SRR034193 SRR034195 SRR034197 SRR034198 SRR034199 SRR034201 SRR034204 SRR034207 SRR034208 SRR034209) 

# SRP001481SRA=( SRR031222 SRR031223 )

for s in \${SRP001481SRA[@]}; do
  #./fastq-dump -DB "@\\\$ac.\\\$si" -DQ "+\\\$ac.\\\$si" \\
    #-A \$s -D input/data/$READACCESSION/\$s/\$s.sra -O output

  #perl prinseq-lite.pl -fastq output/\$s.fastq \\
    #-out_format 1 \\
    #-out_good output/\$s \\
    #-out_bad null \\
    #-min_len 12 \\
    #-min_qual_score 15 \\
    #-ns_max_p 50 \\
    #-noniupac \\
    #-trim_tail_left 5 \\
    #-trim_tail_right 5

  ./segemehl.x -i input/data/$GENOMEACCESSION.idx \\
    -d input/data/$GENOMEACCESSION.fna \\
    -q output/\$s.fasta \\
    --threads 8 \\
    --silent \\
    > output/\$s.map

  rm output/\$s.fasta
  rm output/\$s.fastq
done

cp -r \$OUTPUTDIR \$WORKDIR/
cd
rm -rf \$TMPDIR
EOF
  chmod a+x $BATCH_SH_RUN_MAP
  cp $BATCH_SH_RUN_MAP $CACRUNMAPDIR/
}

function copy-genome-to-cac {
  cp $REFGENOMEFASTA $CACDATADIR
  cp $REFGENOMEINDEX $CACDATADIR
}

function copy-reads-to-cac {
  cp $READDATA/*.fasta $CACDATADIR
}

function copy-data {
  PS3="Choose the study to analyze with rnaseq analysis: "
  select STUDYID in srp001481; do 
    if [ "$STUDYID" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -e "Wait for map-analysis file system preparation...\n"
      prepare-filesystem
      copy-genome-to-cac
      copy-reads-to-cac
      echo -e "Do the next step of run-map.\n"
      break
    fi
  done
}

function run-map {
  PS3="Choose the study to analyze with rnaseq analysis: "
  select STUDYID in srp001481; do 
    if [ "$STUDYID" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -e "Wait for map-analysis file system preparation...\n"
      prepare-filesystem
      copy-batch-sh-run-map
      echo -e "Go to CAC's $STUDYID run-map, and execute nsub batch.sh\n"
      break
    fi
  done
}


# fastx is not used. Instead, I use prinseq.
function filter-reads {
  fastq_quality_filter -Q 33 \
    -i output/$READID.fastq \
    -o output/$READID-filtered.fastq \
    -q 20 -p 90
  fastq_to_fasta -Q 33 \
    -n \
    -i output/$READID.fastq \
    -o output/$READID.fasta
}

# SRA files are converted to a single FASTA file.
# 1. I found a tool
# fastq_to_fasta
# at
# http://hannonlab.cshl.edu/fastx_toolkit/commandline.html
# It was not working with the FASTQ files converted from SRA.
# The fastx_toolkit seems to include some potentially useful tools that I have
# not used yet.
# 
# I read a thread of discussion about fastx -Q option.
# http://seqanswers.com/forums/archive/index.php/t-7399.html
# fastx_quality_stats -Q 33 -i output/SRR031130.fastq -o output/SRR031130.stats
# fastx_nucleotide_distribution_graph.sh -i output/SRR031130.stats -o output/SRR031130-nuc.png -t "My Library"
# fastq_quality_boxplot_graph.sh -i output/SRR031130.stats -o output/SRR031130.png -t "My Library"
function convert-sra-to-fastq {
  fastq-dump -DB "@\$ac.\$si" -DQ "+\$ac.\$si" -A $1 -D $READDATA/$1/$1.sra -O $READDATA
  #fastq_to_fasta -Q 33 -i output/$READID.fastq -o output/$READID.fasta
  #awk 'NR%4==1 || NR%4==2' output/$READID.fastq | sed 's/@/>/' > output/$READID.fasta
}

function prinseq-a-read {
  perl $PRINSEQ -fastq $READDATA/$1.fastq \
    -out_format 1 \
    -out_good $READDATA/$1 \
    -out_bad null \
    -min_len 12 \
    -min_qual_score 15 \
    -ns_max_p 50 \
    -noniupac \
    -trim_tail_left 5 \
    -trim_tail_right 5
}

function match-a-read {
  HpyloriQUERY=output/$1.fasta
  HpyloriMAP=output/$1.map
  segemehl.x -i $HpyloriINDEX -d $HpyloriGENOME -q $HpyloriQUERY > $HpyloriMAP
}

function remove-garbage {
  HpyloriQUERYFASTQ=output/$1.fastq
  HpyloriQUERY=output/$1.fasta
  rm $HpyloriQUERY
  rm $HpyloriQUERYFASTQ
}

# Full read data set.
# SRP001481SRA=( SRR031126 SRR031127 SRR031128 SRR031129 SRR031130 SRR031131 SRR031142 SRR031143 SRR031144 SRR031145 SRR031146 SRR031178 SRR031179 SRR031180 SRR031181 SRR031182 SRR031183 SRR031184 SRR031185 SRR031186 SRR031187 SRR031188 SRR031189 SRR031190 SRR031191 SRR031192 SRR031193 SRR031194 SRR031195 SRR031196 SRR031197 SRR031198 SRR031199 SRR031200 SRR031201 SRR031202 SRR031203 SRR031204 SRR031205 SRR031206 SRR031207 SRR031208 SRR031209 SRR031210 SRR031211 SRR031212 SRR031213 SRR031214 SRR031215 SRR031216 SRR031217 SRR031218 SRR031219 SRR031220 SRR031221 SRR031222 SRR031223 SRR031224 SRR031225 SRR031226 SRR031276 SRR031277 SRR031278 SRR031279 SRR031280 SRR031281 SRR031282 SRR031283 SRR031284 SRR031285 SRR031286 SRR031287 SRR031288 SRR031289 SRR031290 SRR031291 SRR031292 SRR031293 SRR031294 SRR031295 SRR034162 SRR034163 SRR034165 SRR034166 SRR034167 SRR034168 SRR034169 SRR034170 SRR034171 SRR034172 SRR034173 SRR034174 SRR034175 SRR034176 SRR034177 SRR034178 SRR034179 SRR034180 SRR034181 SRR034182 SRR034183 SRR034185 SRR034186 SRR034187 SRR034193 SRR034195 SRR034197 SRR034198 SRR034199 SRR034201 SRR034204 SRR034207 SRR034208 SRR034209) 

# Some of them
SRP001481SRA=( SRR031188 SRR031189 SRR031190 SRR031191 SRR031192 SRR031193 SRR031194 SRR031195 SRR031196 SRR031197 SRR031198 SRR031199 SRR031200 SRR031201 SRR031202 SRR031203 SRR031204 SRR031205 SRR031206 SRR031207 SRR031208 SRR031209 SRR031210 SRR031211 SRR031212 SRR031213 SRR031214 SRR031215 SRR031216 SRR031217 SRR031218 SRR031219 SRR031220 SRR031221 SRR031222 SRR031223 SRR031224 SRR031225 SRR031226 SRR031276 SRR031277 SRR031278 SRR031279 SRR031280 SRR031281 SRR031282 SRR031283 SRR031284 SRR031285 SRR031286 SRR031287 SRR031288 SRR031289 SRR031290 SRR031291 SRR031292 SRR031293 SRR031294 SRR031295 SRR034162 SRR034163 SRR034165 SRR034166 SRR034167 SRR034168 SRR034169 SRR034170 SRR034171 SRR034172 SRR034173 SRR034174 SRR034175 SRR034176 SRR034177 SRR034178 SRR034179 SRR034180 SRR034181 SRR034182 SRR034183 SRR034185 SRR034186 SRR034187 SRR034193 SRR034195 SRR034197 SRR034198 SRR034199 SRR034201 SRR034204 SRR034207 SRR034208 SRR034209) 

# For test ...
# SRP001481SRA=( SRR031222 SRR031223 )

function run-filter {

  PS3="Choose the study to analyze with rnaseq analysis: "
  select STUDYID in srp001481; do 
    if [ "$STUDYID" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -e "Wait for map-analysis file system preparation...\n"
      prepare-filesystem
      for s in ${SRP001481SRA[@]}; do
        # convert-sra-to-fastq $s
        prinseq-a-read $s
        # match-a-read $s
        # remove-garbage $s
      done
      echo -e "Do the next step of copy-data!\n"
      break
    fi
  done


}

# A bacterial genome is used to build an index file.
# See http://www.bioinf.uni-leipzig.de/Software/segemehl/ for the manual of the
# program.
function index-genome {
  PS3="Choose the study to analyze with rnaseq analysis: "
  select STUDYID in srp001481; do 
    if [ "$STUDYID" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -e "Wait for map-analysis file system preparation...\n"
      prepare-filesystem
      mkdir-STUDYID
      segemehl.x -x $REFGENOMEINDEX -d $REFGENOMEFASTA
      echo -e "Do the next step of run-map!\n"
      break
    fi
  done
}

#####################################################################
# Main part of the script.
#####################################################################
PS3="Select what you want to do with rnaseq-analysis: "
CHOICES=( index-genome run-filter copy-data run-map srp001481 build-reads prinseq-reads filter-reads build-index match-reads )
select CHOICE in ${CHOICES[@]}; do 
 
  if [ "$CHOICE" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  elif [ "$CHOICE" == "index-genome" ];  then
    index-genome
    break
  elif [ "$CHOICE" == "run-filter" ];  then
    run-filter
    break
  elif [ "$CHOICE" == "copy-data" ];  then
    copy-data
    break
  elif [ "$CHOICE" == "run-map" ];  then
    run-map
    break

  elif [ "$CHOICE" == "srp001481" ];  then
    srp001481
    break
  elif [ "$CHOICE" == "prinseq-reads" ];  then
    prinseq-reads
    break
  elif [ "$CHOICE" == "build-reads" ];  then
    build-reads
    break
  elif [ "$CHOICE" == "match-reads" ];  then
    match-reads
    break
  elif [ "$CHOICE" == "receive-run-mauve" ];  then
    break
  elif [ "$CHOICE" == "prepare-run-clonalframe" ];  then
    break
  elif [ "$CHOICE" == "compute-watterson-estimate-for-clonalframe" ];  then
    break
  elif [ "$CHOICE" == "receive-run-clonalframe" ];  then
    break
  elif [ "$CHOICE" == "prepare-run-clonalorigin" ];  then
    break
  elif [ "$CHOICE" == "receive-run-clonalorigin" ];  then
    break
  elif [ "$CHOICE" == "receive-run-2nd-clonalorigin" ];  then
    break
  else
    echo -e "You need to enter something\n"
    continue
  fi
done

