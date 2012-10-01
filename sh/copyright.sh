function short-notice {
  cat <<EOF
RNASeq Analysis Copyright (C) 2011 Sang Chul Choi
This program comes with ABSOLUTELY NO WARRANTY; for details select menu warranty.
This is free software, and you are welcome to redistribute it
under certain conditions; select menu copyright for details.
EOF
}

function warranty {
  cat <<EOF
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
EOF
}

function copyright {
  less COPYING
}

function quit {
  cat <<EOF
RNA-seq Analysis - Bye!
EOF
}

