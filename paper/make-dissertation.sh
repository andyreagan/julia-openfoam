#!/bin/bash

FILENAME_BASE="reagan-thesis"

# set STYLE to JA for a journal article style thesis
#      (and set the \input's to \include's in the main)
STYLE="standard"

# first test the body versus our style scripts
for FILE in $FILENAME_BASE.chapter*
do
  echo "checking $FILE into warnings.txt"
  perl bin/lexical-illusion-spotter.pl $FILE
  bin/weasel-word-spotter.sh $FILE
  bin/passive-word-spotter.sh $FILE
done

pdflatex $FILENAME_BASE-main

if [ $STYLE == "JA" ]; then
  for auxfile in $FILENAME_BASE.chapter*-JA.aux
  do
      bibtex `basename $auxfile .aux`
  done
  
  for bblfile in $FILENAME_BASE.chapter*-JA.bbl
  do
      sed -i 's/thebibliography/references/' $bblfile
  done
fi

bibtex $FILENAME_BASE-main #>> make-dissertation.texout
pdflatex $FILENAME_BASE-main #>> make-dissertation.texout
pdflatex $FILENAME_BASE-main #>> make-dissertation.texout
pdflatex $FILENAME_BASE-main #>> make-dissertation.texout

echo " "
echo " "
echo "There are $(($(grep todo *.tex | wc -l | awk '{print $1;}')-2)) to-dos!!"

\rm *.toc *.bbl *.aux *.log *.lot



