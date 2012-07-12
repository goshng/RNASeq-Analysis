function conf {
  CONFFILE=conf/README
  PROJECTNAME=$(grep ^PROJECTNAME\: $CONFFILE | cut -d":" -f2)
  RUNMODE=$(grep ^RUNMODE\: $CONFFILE | cut -d":" -f2)
  NUMBERCPU=$(grep ^NUMBERCPU\: $CONFFILE | cut -d":" -f2)
  ROOTANALYSISDIR=$(grep ^ROOTANALYSISDIR\: $CONFFILE | cut -d":" -f2)

  CAC_USERNAME=$(grep ^CAC_USERNAME\: $CONFFILE | cut -d":" -f2)
  CAC_LOGIN=$(grep ^CAC_LOGIN\: $CONFFILE | cut -d":" -f2)
  CAC_ROOT=$(grep ^CAC_ROOT\: $CONFFILE | cut -d":" -f2)
  X11_USERNAME=$(grep ^X11_USERNAME\: $CONFFILE | cut -d":" -f2)
  X11_LOGIN=$(grep ^X11_LOGIN\: $CONFFILE | cut -d":" -f2)
  X11_ROOT=$(grep ^X11_ROOT\: $CONFFILE | cut -d":" -f2)
  BATCHEMAIL=$(grep ^BATCHEMAIL\: $CONFFILE | cut -d":" -f2)
  EMAILON=$(grep ^EMAILON\: $CONFFILE | cut -d":" -f2)
  if [ $EMAILON == "TRUE" ]; then
    EMAILON=""
  else
    EMAILON="# "
  fi
  BATCHACCESS=$(grep ^BATCHACCESS\: $CONFFILE | cut -d":" -f2)
  QUEUENAME=$(grep ^QUEUENAME\: $CONFFILE | cut -d":" -f2)

  MAQ=$(grep ^MAQ\: $CONFFILE | cut -d":" -f2)
  BWA=$(grep ^BWA\: $CONFFILE | cut -d":" -f2)
  PICARD=$(grep ^PICARD\: $CONFFILE | cut -d":" -f2)
  PRINSEQ=$(grep ^PRINSEQ\: $CONFFILE | cut -d":" -f2)
  CUTADAPT=$(grep ^CUTADAPT\: $CONFFILE | cut -d":" -f2)
  PYTHON=$(grep ^PYTHON\: $CONFFILE | cut -d":" -f2)
  CACRSCRIPT=$(grep ^CACRSCRIPT\: $CONFFILE | cut -d":" -f2)
  CRAMTOOLS=$(grep ^CRAMTOOLS\: $CONFFILE | cut -d":" -f2)
  PARSERNASEQ=$(grep ^PARSERNASEQ\: $CONFFILE | cut -d":" -f2)
  BOWTIE=$(grep ^BOWTIE\: $CONFFILE | cut -d":" -f2)
  SAMTOOLS=$(grep ^SAMTOOLS\: $CONFFILE | cut -d":" -f2)
  SUBREADBUILDINDEX=$(grep ^SUBREADBUILDINDEX\: $CONFFILE | cut -d":" -f2)
  SUBREADALIGN=$(grep ^SUBREADALIGN\: $CONFFILE | cut -d":" -f2)
  SEGEMEHL=$(grep ^SEGEMEHL\: $CONFFILE | cut -d":" -f2)

  # Other global variables
  CAC_USERHOST=$CAC_USERNAME@$CAC_LOGIN
  CAC_MAUVEANALYSISDIR=$CAC_USERHOST:$CAC_ROOT
  CAC_OUTPUTDIR=$CAC_ROOT/output
  CACBASE=$CAC_ROOT/output

  # X11 linux ID setup
  X11_USERHOST=$X11_USERNAME@$X11_LOGIN
  X11_MAUVEANALYSISDIR=$X11_USERHOST:$X11_ROOT
  X11_OUTPUTDIR=$CAC_ROOT/output
  X11BASE=$X11_ROOT/output

  # The main base directory contains all the subdirectories.
  OUTPUTDIR=$ROOTANALYSISDIR/output
  ROUTPUTDIR=$(grep ^ROUTPUTDIR $CONFFILE | cut -d":" -f2)
}
