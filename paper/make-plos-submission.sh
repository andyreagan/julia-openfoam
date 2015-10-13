#!/bin/bash
#
# usage: ./make-plos-submission.sh thermosyphon-paper-plos

if [ -d plos ]; then
    \rm -r plos-package
fi
mkdir -p plos-package
echo "making a single latex file..."
./make-single-latex-file.pl $1.tex
echo "finding figures and linking..."
./make-links.pl $1-combined.tex plos-package

echo 'creating $1-package.tgz...'

cd plos-package

# cover letter
# ln -sf ../correspondence/2015-10-11plos-one.pdf .

# main pdf
ln -sf ../$1.pdf

# create that tar archive
tar cvfzph ../$1-package.tgz fig[0-9][0-9]*.tiff $1-combined.tex  $1.pdf
