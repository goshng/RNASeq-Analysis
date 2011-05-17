#!/bin/bash
# Author: Sang Chul Choi

# READID=SRR031222
# READID=SRR031212
# READID=SRR031130
# READID=test-illumina
#####################################################################################
# Sample: PL TEX plus library from H. pylori strain 26695(SRS008287)
SRX014058SRA=( SRR031222 SRR031223 SRR031224 SRR031225 SRR031226 SRR031294 SRR031295 SRR034198 SRR034199 SRR034201 )
# Sample: PL TEX minus library from H. pylori strain 26695(SRS008286)
SRX014057SRA=( SRR031217 SRR031218 SRR031219 SRR031220 SRR031221 SRR031292 SRR031293 SRR034193 SRR034195 SRR034197 )
# Sample: AGS TEX plus library from H. pylori strain 26695(SRS008283)
SRX014056SRA=( SRR031212 SRR031213 SRR031214 SRR031215 SRR031216 SRR031290 SRR031291 SRR034166 SRR034167 SRR034168 )
# Sample: AGS TEX minus library from H. pylori strain 26695(SRS008282)
SRX014055SRA=( SRR031207 SRR031208 SRR031209 SRR031210 SRR031211 SRR031288 SRR031289 SRR034162 SRR034163 SRR034165 )
# Sample: Huh7 TEX  plus library from H. pylori strain 26695(SRS008281)
SRX014054SRA=( SRR031202 SRR031203 SRR031204 SRR031205 SRR031206 SRR031286 SRR031287 SRR034178 SRR034179 SRR034180 )
# Sample: Huh7 TEX minus library from H. pylori strain 26695(SRS008280)
SRX014053SRA=( SRR031197 SRR031198 SRR031199 SRR031200 SRR031201 SRR031284 SRR031285 SRR034175 SRR034176 SRR034177 )
# Sample: Acid-stress TEX minus library from H. pylori 26695(SRS008274)
SRX014052SRA=( SRR031192 SRR031193 SRR031194 SRR031195 SRR031196 SRR031282 SRR031283 SRR034169 SRR034170 SRR034171 SRR034204 )
# Sample: Acid-stress TEX  plus library from H. pylori 26695(SRS008273)
SRX014051SRA=( SRR031187 SRR031188 SRR031189 SRR031190 SRR031191 SRR031280 SRR031281 SRR034172 SRR034173 SRR034174 SRR034207 )
# Sample: Mid-log growthTEX minus library from H. pylori 26695(SRS008272)
SRX014050SRA=( SRR031182 SRR031183 SRR031184 SRR031185 SRR031186 SRR031278 SRR031279 SRR034181 SRR034182 SRR034183 SRR034208 )
# Sample: Mid-log growthTEX plus library from H. pylori 26695(SRS008271)
SRX014018SRA=( SRR031131 SRR031178 SRR031179 SRR031180 SRR031181 SRR031276 SRR031277 SRR034185 SRR034186 SRR034187 SRR034209 )
# Sample: Solexa RNA-seq library from H. pylori 26695 grown in the presence of Huh7 liver cells(SRS008279)
SRX014017SRA=( SRR031130 )
# Sample: Solexa RNA-seq library from H. pylori 26695 grown in the presence of AGS cells(SRS008278)
SRX014016SRA=( SRR031129 )
# Sample: Solexa RNA-seq library of H. pylori 26695 grown in cell culture medium(SRS008277)
SRX014015SRA=( SRR031128 )
# Sample: Solexa RNA-seq library of  total RNA from H. pylori grown under acid stress(SRS008276)
SRX014014SRA=( SRR031127 )
# Sample: Solexa RNA-seq library of  total RNA from H. pylori grown to mid-log phase(SRS008275)
SRX014013SRA=( SRR031126 )

# Full read data set.
GA454SRA=( SRX014058 SRX014057 SRX014056 SRX014055 SRX014054 SRX014053 SRX014052 SRX014051 SRX014050 SRX014018 ) 
SOLEXASRA=( SRR031130 SRR031129 SRR031128 SRR031127 SRR031126 ) 
TOTALSRA=( "${GA454SRA[@]}" "${SOLEXASRA[@]}" )
#####################################################################################


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
  TEMPDIRECTORY=/tmp/$JOBID.scheduler.v4linux
  TMPINPUTDIR=$TEMPDIRECTORY/input
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
#PBS -l walltime=48:00:00,nodes=1
#PBS -A acs4_0001
#PBS -j oe
#PBS -N Strep-${STUDYID}-Map
#PBS -q v4
##PBS -m e
##PBS -M schoi@cornell.edu
WORKDIR=\$PBS_O_WORKDIR
DATADIR=\$WORKDIR/../data
FASTQDUMP=\$HOME/usr/bin/fastq-dump
PRINSEQ=\$HOME/usr/bin/prinseq-lite.pl
SEGEMEHL=\$HOME/usr/bin/segemehl.x

OUTPUTDIR=\$TEMPDIRECTORY/output
INPUTDIR=\$TEMPDIRECTORY/input
mkdir \$WORKDIR/output
mkdir \$INPUTDIR
mkdir \$OUTPUTDIR
cp \$FASTQDUMP \$TEMPDIRECTORY/
cp \$PRINSEQ \$TEMPDIRECTORY/
cp \$SEGEMEHL \$TEMPDIRECTORY/
cp \$DATADIR/$GENOMEACCESSION.* \$INPUTDIR/

cd \$TEMPDIRECTORY
for s in ${TOTALSRA[@]}; do
  cp \$DATADIR/filtered-\$s.fasta \$INPUTDIR

  ./segemehl.x -i input/$GENOMEACCESSION.idx \\
    -d input/$GENOMEACCESSION.fna \\
    -q input/filtered-\$s.fasta \\
    --threads 8 \\
    --silent \\
    > output/\$s.map
  mv output/\$s.map \$WORKDIR/output/
  rm input/filtered-\$s.fasta
done

cd
rm -rf \$TEMPDIRECTORY
EOF
  chmod a+x $BATCH_SH_RUN_MAP
  cp $BATCH_SH_RUN_MAP $CACRUNMAPDIR/
}

function copy-genome-to-cac {
  cp $REFGENOMEFASTA $CACDATADIR/
  cp $REFGENOMEINDEX $CACDATADIR/
}

function copy-reads-to-cac {
  cp $READDATA/filtered*.fasta $CACDATADIR
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

# Manual preprocessing may be necessary.
# For 454 data sets, 5' 18 and 3' 120 are trimmed.
# -trim_left 18
# -trim_right 120
function prinseq-a-read {
  if [ "$2" == "GA454" ]; then
    perl $PRINSEQ -fastq $READDATA/$1.fastq \
      -out_format 3 \
      -out_good $READDATA/filtered-$1 \
      -out_bad null \
      -min_len 12 \
      -min_qual_mean 25 \
      -ns_max_p 20 \
      -max_a 70 \
      -noniupac \
      -trim_left 18 \
      -trim_right 120 \
      -trim_tail_left 5 \
      -trim_tail_right 5
  else
    perl $PRINSEQ -fastq $READDATA/$1.fastq \
      -out_format 3 \
      -out_good $READDATA/filtered-$1 \
      -out_bad null \
      -min_len 12 \
      -min_qual_mean 25 \
      -ns_max_p 20 \
      -max_a 70 \
      -noniupac \
      -trim_tail_left 5 \
      -trim_tail_right 5
  fi
    #-stats_info -stats_len -stats_dinuc -stats_tag -stats_dupl -stats_ns \
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

SRP001481SRA=( SRR031126 SRR031127 SRR031128 SRR031129 SRR031130 SRR031131 SRR031142 SRR031143 SRR031144 SRR031145 SRR031146 SRR031178 SRR031179 SRR031180 SRR031181 SRR031182 SRR031183 SRR031184 SRR031185 SRR031186 SRR031187 SRR031188 SRR031189 SRR031190 SRR031191 SRR031192 SRR031193 SRR031194 SRR031195 SRR031196 SRR031197 SRR031198 SRR031199 SRR031200 SRR031201 SRR031202 SRR031203 SRR031204 SRR031205 SRR031206 SRR031207 SRR031208 SRR031209 SRR031210 SRR031211 SRR031212 SRR031213 SRR031214 SRR031215 SRR031216 SRR031217 SRR031218 SRR031219 SRR031220 SRR031221 SRR031222 SRR031223 SRR031224 SRR031225 SRR031226 SRR031276 SRR031277 SRR031278 SRR031279 SRR031280 SRR031281 SRR031282 SRR031283 SRR031284 SRR031285 SRR031286 SRR031287 SRR031288 SRR031289 SRR031290 SRR031291 SRR031292 SRR031293 SRR031294 SRR031295 SRR034162 SRR034163 SRR034165 SRR034166 SRR034167 SRR034168 SRR034169 SRR034170 SRR034171 SRR034172 SRR034173 SRR034174 SRR034175 SRR034176 SRR034177 SRR034178 SRR034179 SRR034180 SRR034181 SRR034182 SRR034183 SRR034185 SRR034186 SRR034187 SRR034193 SRR034195 SRR034197 SRR034198 SRR034199 SRR034201 SRR034204 SRR034207 SRR034208 SRR034209) 

function run-fastq-dump {

  PS3="Choose the study to analyze with rnaseq analysis: "
  select STUDYID in srp001481; do 
    if [ "$STUDYID" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -e "Wait for map-analysis file system preparation...\n"
      prepare-filesystem
      for s in ${SRP001481SRA[@]}; do
        convert-sra-to-fastq $s
      done
      echo -e "Do the next step of combine-runs!\n"
      break
    fi
  done
}

function combine-runs {

  PS3="Choose the study to analyze with rnaseq analysis: "
  select STUDYID in srp001481; do 
    if [ "$STUDYID" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -e "Wait for map-analysis file system preparation...\n"
      prepare-filesystem

      rm -f $READDATA/SRX014058.fastq 
      for s in ${SRX014058SRA[@]}; do
        cat $READDATA/$s.fastq >> $READDATA/SRX014058.fastq 
      done
      rm -f $READDATA/SRX014057.fastq 
      for s in ${SRX014057SRA[@]}; do
        cat $READDATA/$s.fastq >> $READDATA/SRX014057.fastq 
      done
      rm -f $READDATA/SRX014056.fastq 
      for s in ${SRX014056SRA[@]}; do
        cat $READDATA/$s.fastq >> $READDATA/SRX014056.fastq 
      done
      rm -f $READDATA/SRX014055.fastq 
      for s in ${SRX014055SRA[@]}; do
        cat $READDATA/$s.fastq >> $READDATA/SRX014055.fastq 
      done
      rm -f $READDATA/SRX014054.fastq 
      for s in ${SRX014054SRA[@]}; do
        cat $READDATA/$s.fastq >> $READDATA/SRX014054.fastq 
      done
      rm -f $READDATA/SRX014053.fastq 
      for s in ${SRX014053SRA[@]}; do
        cat $READDATA/$s.fastq >> $READDATA/SRX014053.fastq 
      done
      rm -f $READDATA/SRX014052.fastq 
      for s in ${SRX014052SRA[@]}; do
        cat $READDATA/$s.fastq >> $READDATA/SRX014052.fastq 
      done
      rm -f $READDATA/SRX014051.fastq 
      for s in ${SRX014051SRA[@]}; do
        cat $READDATA/$s.fastq >> $READDATA/SRX014051.fastq 
      done
      rm -f $READDATA/SRX014050.fastq 
      for s in ${SRX014050SRA[@]}; do
        cat $READDATA/$s.fastq >> $READDATA/SRX014050.fastq 
      done
      rm -f $READDATA/SRX014018.fastq 
      for s in ${SRX014018SRA[@]}; do
        cat $READDATA/$s.fastq >> $READDATA/SRX014018.fastq 
      done

      echo -e "Do the next step of run-filter!\n"
      break
    fi
  done
}

function run-filter {
  PS3="Choose the study to analyze with rnaseq analysis: "
  select STUDYID in srp001481; do 
    if [ "$STUDYID" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -e "Wait for map-analysis file system preparation...\n"
      prepare-filesystem
      for s in ${GA454SRA[@]}; do
        prinseq-a-read $s GA454
      done
      for s in ${SOLEXASRA[@]}; do
        prinseq-a-read $s
      done
      echo -e "Do the next step of run-fastq-to-fasta!\n"
      break
    fi
  done
}


function run-fastq-to-fasta {
  PS3="Choose the study to analyze with rnaseq analysis: "
  select STUDYID in srp001481; do 
    if [ "$STUDYID" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -e "Wait for map-analysis file system preparation...\n"
      prepare-filesystem
      for s in ${GA454SRA[@]}; do
        fastq_to_fasta -Q 33 -n -i $READDATA/filtered-$s.fastq -o $READDATA/filtered-$s.fasta 
      done
      for s in ${SOLEXASRA[@]}; do
        fastq_to_fasta -Q 33 -n -i $READDATA/filtered-$s.fastq -o $READDATA/filtered-$s.fasta 
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

function receive-run-map {
  PS3="Choose the study to analyze with rnaseq analysis: "
  select STUDYID in srp001481; do 
    if [ "$STUDYID" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -e "Receiving map ...\n"
      prepare-filesystem
      cp -r $CACRUNMAPDIR/output $RUNMAPDIR/
      echo -e "Do the next step of run-map!\n"
      break
    fi
  done
}

function run-compute-graph {
  PS3="Choose the study to analyze with rnaseq analysis: "
  select STUDYID in srp001481; do 
    if [ "$STUDYID" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -e "Computing graphs ...\n"
      prepare-filesystem
      # perl pl/map2graph.pl -map output/srp001481/run-map/output/SRR031126.map -genome_length 1667867
      #for s in ${GA454SRA[@]}; do
        #perl pl/map2graph.pl -map output/$STUDYID/run-map/output/$s.map -genome_length 1667867 
      #done
      for s in ${SOLEXASRA[@]}; do
        perl pl/map2graph.pl -map output/$STUDYID/run-map/output/$s.map -genome_length 1667867 
      done
      echo -e "Do the next step of run-???!\n"
      break
    fi
  done
}


function run-fastq-to-fasta {
  PS3="Choose the study to analyze with rnaseq analysis: "
  select STUDYID in srp001481; do 
    if [ "$STUDYID" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -e "Wait for map-analysis file system preparation...\n"
      prepare-filesystem
      for s in ${GA454SRA[@]}; do
        fastq_to_fasta -Q 33 -n -i $READDATA/filtered-$s.fastq -o $READDATA/filtered-$s.fasta 
      done
      for s in ${SOLEXASRA[@]}; do
        fastq_to_fasta -Q 33 -n -i $READDATA/filtered-$s.fastq -o $READDATA/filtered-$s.fasta 
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

function receive-run-map {
  PS3="Choose the study to analyze with rnaseq analysis: "
  select STUDYID in srp001481; do 
    if [ "$STUDYID" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -e "Receiving map ...\n"
      prepare-filesystem
      cp -r $CACRUNMAPDIR/output $RUNMAPDIR/
      echo -e "Do the next step of run-map!\n"
      break
    fi
  done
}

#####################################################################
# Main part of the script.
#####################################################################
PS3="Select what you want to do with rnaseq-analysis: "
CHOICES=( index-genome run-fastq-dump combine-runs run-filter run-fastq-to-fasta copy-data run-map receive-run-map run-compute-graph srp001481 build-reads prinseq-reads filter-reads build-index match-reads )
    
    
    

select CHOICE in ${CHOICES[@]}; do 
 
  if [ "$CHOICE" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  elif [ "$CHOICE" == "index-genome" ];  then
    index-genome
    break
  elif [ "$CHOICE" == "run-fastq-dump" ];  then
    run-fastq-dump
    break
  elif [ "$CHOICE" == "combine-runs" ];  then
    combine-runs
    break
  elif [ "$CHOICE" == "run-filter" ];  then
    run-filter
    break
  elif [ "$CHOICE" == "run-fastq-to-fasta" ];  then
    run-fastq-to-fasta
    break
  elif [ "$CHOICE" == "copy-data" ];  then
    copy-data
    break
  elif [ "$CHOICE" == "run-map" ];  then
    run-map
    break
  elif [ "$CHOICE" == "receive-run-map" ];  then
    receive-run-map
    break
  elif [ "$CHOICE" == "run-compute-graph" ];  then
    run-compute-graph
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
  else
    echo -e "You need to enter something\n"
    continue
  fi
done

