###############################################################################
# Functions: structuring file system (or creating directories)
###############################################################################
# Create initial directories
# --------------------------
# After checking out the source code from the repository output directories need
# to be created in the local, remote cluster, remote X11 machines.
# input: nothing
# output: 3 output directories
function init-file-system {
  echo -n "Creating $BASEDIR/output ..." 
  mkdir $BASEDIR/output 
  echo -e " done"
  echo -n "Creating $CAC_ROOT/output at $CAC_USERHOST ..."
  ssh -x $CAC_USERHOST mkdir -p $CAC_ROOT/output
  echo -e " done"
  echo -n "Creating $X11_ROOT/output at $X11_USERHOST ..."
  ssh -x $X11_USERHOST mkdir -p $X11_ROOT/output
  echo -e " done"
}
