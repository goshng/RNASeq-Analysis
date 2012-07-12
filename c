#!/bin/bash
# asciidoc -b html5 -a icons -a toc doc/Manual; open doc/Manual.html
ASCIIDOC=/opt/local/bin/asciidoc
$ASCIIDOC -a toc doc/Manual
#open doc/Manual.html
#asciidoc -a toc ChangeLog; open ChangeLog.html
 
