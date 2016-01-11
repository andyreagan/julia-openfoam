#!/bin/bash
#
# usage: ./make-plos-submission.sh thermosyphon-paper-plos

pdflatex $1
bibtex $1
pdflatex $1
pdflatex $1
# ./make-plos-pieces.pl thermosyphon-paper.tex
# ./make-plos-pieces.pl $1.tex
# don't need that...

# if [ -d plos-package ]; then
#     \rm -r plos-package
#     mkdir -p plos-package
# else
#     mkdir -p plos-package
# fi
echo "making a single latex file..."
./make-single-latex-file-plos.pl $1.tex
echo "finding figures and linking..."
./make-links-plos.pl $1-combined.tex plos-package
# # already done by the above
# # cp $1-combined.tex plos-package

cd plos-package

# # already done inside the make-links.py
# for pdffile in $(\ls -1 *.pdf)
# do
#     echo "will tiffify $pdffile"
# done

pdflatex $1-combined
bibtex $1-combined
pdflatex $1-combined
pdflatex $1-combined

# does the supplementary and the main pdf need to be split??

echo 'creating $1-package.tgz...'

# # cover letter
# # ln -sf ../correspondence/2015-10-11plos-one.pdf .
ln -sf ../thermosyphon-paper-coverletter.pdf

# # main pdf
# ln -sf ../$1.pdf

../breakup-manuscript-supp.pl thermosyphon-paper-plos-combined.pdf $(cat startsupp-plos.txt)

# create that tar archive
# tar cvfzph ../$1-package.tgz fig[0-9][0-9]*.tiff $1-combined.tex  $1.pdf
tar cvzphf ../$1-package.tgz fig*.pdf fig*.tiff thermosyphon-paper-coverletter.pdf $1-combined.tex  $1-combined-{manuscript,supplementary}.pdf
# tar cvzphf ../thermosyphon-paper-plos-package.tgz fig*.pdf fig*.tiff thermosyphon-paper-coverletter.pdf thermosyphon-paper-plos-combined.tex  thermosyphon-paper-plos-combined.pdf

# cd ..



# NOTES

# made supp sections into subsections and removed the subsections that had existed.
# replace all sections with \section{ -> \section*{
                                         
