function conf {
  CONFFILE=conf/README
  PROJECTNAME=$(grep PROJECTNAME $CONFFILE | cut -d":" -f2)
  CAC_USERNAME=$(grep CAC_USERNAME $CONFFILE | cut -d":" -f2)
  CAC_LOGIN=$(grep CAC_LOGIN $CONFFILE | cut -d":" -f2)
  CAC_ROOT=$(grep CAC_ROOT $CONFFILE | cut -d":" -f2)
  X11_USERNAME=$(grep X11_USERNAME $CONFFILE | cut -d":" -f2)
  X11_LOGIN=$(grep X11_LOGIN $CONFFILE | cut -d":" -f2)
  X11_ROOT=$(grep X11_ROOT $CONFFILE | cut -d":" -f2)
  BATCHEMAIL=$(grep BATCHEMAIL $CONFFILE | cut -d":" -f2)
  BATCHACCESS=$(grep BATCHACCESS $CONFFILE | cut -d":" -f2)
  QUEUENAME=$(grep QUEUENAME $CONFFILE | cut -d":" -f2)

  BWA=$(grep BWA $CONFFILE | cut -d":" -f2)
  SAMTOOL=$(grep SAMTOOL $CONFFILE | cut -d":" -f2)

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
}
