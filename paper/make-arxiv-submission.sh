#!/bin/bash
#
# usage: ./make-arxiv-submission.sh thermosyphon-paper-revtex4

if [ -d arxiv-packagee ]; then
    \rm -r arxiv-package
fi
mkdir -p arxiv-package
echo "making a single latex file..."
./make-single-latex-file.pl $1.tex
echo "finding figures and linking..."
./make-links.pl $1-combined.tex arxiv-package

## this makes pdf's and such: moved to makelinks.pl
## cd package
## \ls *.pdf | xargs -n 1 pdf2ps
## \rm *.pdf
## \ls *.ps | xargs -n 1 ps2eps
## \rm *.ps
## for figure in $(\ls -1 *.png *.jpg); do
##   convert $figure{,.ps}
## done

echo "creating $1-arxiv.tgz..." 1>&2;
cd arxiv-package;
tar cvfzph ../$1-arxiv.tgz * 1>&2;
cd ..

# cleanup
\rm -r arxiv-package
\rm $1-combined.tex
\rm $1-arxiv.tgz
