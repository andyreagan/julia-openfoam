#!/usr/bin/perl
# make links into specified directory

$usage = "
usage: make-links-plos.pl foo-combined.tex directory
(works for any latex file)

use after make-single-latex-file.pl

copies foo-combined.tex to ./directory
and replaces figure names by fig1, fig2, etc.

see end of file for localized pieces
";

if ($#ARGV < 0)
{
    print $usage;
    exit;
}

$file = $ARGV[0];
$folder = $ARGV[1];

unless (-e $file)
{
    print $usage;
    exit;
}

# clean out the package directory
# `\\rm package/*`;

# find figures

$outfile = "$folder/$file";
open (FILE,$file) or die "can't open $file: $!\n";
open (OUTFILE,">$outfile") or die "can't open $outfile: $!\n";

$i = 0;
$j = 1;
$caption = 0;
foreach $line (<FILE>)
{
    unless ($line =~ m/\s*%/) {
	if ($line =~ m/\\includegraphics\[.*?\]\{(.*?)\}/) {
	    $figfile = $1;
	    push @fullfigures, $figfile;
	    $fullfigure = $figfile;
            print "fullfigure=$fullfigure\n";
            # remove the extension
	    $figfile =~ s/^.*\///;
            print "figfile=$figfile\n";
	    push @figures, $figfile;
            # don't 0 pad...
	    $fignum = $i;
            # prefix the new file
	    $prefix = "fig".$fignum."_";
            # copy it in
            print "cp $fullfigure $folder/$prefix$figfile\n";
	    `cp $fullfigure $folder/$prefix$figfile`;
            if ($figfile =~ m/\.ps/) {
                $figfile =~ s/\.ps//;
		`pstopdf $folder/$prefix$figfile.ps`;
		`rm $folder/$prefix$figfile.ps`;
                # `./tiffify.pl $folder/$prefix$figfile.{pdf,tiff}`;
                # `rm $folder/$prefix$figfile.pdf`;
	    }
            if ($figfile =~ m/\.eps/) {
                $figfile =~ s/\.eps//;
		`epstopdf $folder/$prefix$figfile.eps`;
		`rm $folder/$prefix$figfile.eps`;
                if ($i > 0) {
                    # `./tiffify.pl $folder/$prefix$figfile.{pdf,tiff}`;
                    `rm $folder/$prefix$figfile.pdf`;
                }
	    }
            if ($figfile =~ m/\.pdf/) {
                $figfile =~ s/\.pdf//;
                # `./tiffify.pl $folder/$prefix$figfile.{pdf,tiff}`;
                `rm $folder/$prefix$figfile.pdf`;
	    }
            if ($figfile =~ m/\.jpg/) {
                $figfile =~ s/\.jpg//;
                `convert $folder/$prefix$figfile.{jpg,pdf}`;
                # `./tiffify.pl $folder/$prefix$figfile.{pdf,tiff}`;
                `rm $folder/$prefix$figfile.pdf`;
	    }
            if ($figfile =~ m/\.png/) {
                $figfile =~ s/\.png//;
                `convert $folder/$prefix$figfile.{png,pdf}`;
                `rm $folder/$prefix$figfile.png`;
                # `./tiffify.pl $folder/$prefix$figfile.{pdf,tiff}`;
                `rm $folder/$prefix$figfile.pdf`;
	    }
            if ($i > 0) {
                print "commenting out figure\n";
                $line =~ s/(\\includegraphics.*?)\{(.*?)\}/%% $1\{$prefix$figfile.tiff\}/;
            }
            else {
                print "including first figure\n";
                $line =~ s/(\\includegraphics.*?)\{(.*?)\}/$1\{$prefix$figfile.pdf\}/;
            }
            # $line =~ s/\.ps/.tiff/;
            # $line =~ s/\.eps/.tiff/;
            # $line =~ s/\.pdf/.tiff/;
	    # $line =~ s/\.jpg/.tiff/;
	    # $line =~ s/\.png/.tiff/;
 	    $i = $i + 1;
	}
        if ($line =~ m/\\lstinputlisting\[.*?\]\{(.*?)\}/) {
            print "found line of code: $line";
	    $codefile = $1;
            print "codefile=$codefile\n";
            $fullcode = $codefile;
	    $codefile =~ s/^.*\///;
            if ($j<10) {$codenum = "0$j";} else {$codenum = $j;}
	    $prefix = "code".$codenum."_";
	    `cp $fullcode $folder/$prefix$codefile`;
	    $line =~ s/(\\lstinputlisting.*?)\{(.*?)\}/$1\{$prefix$codefile\}/;
 	    $j = $j + 1;
	}
        if ($caption) {
            $caption = 0;
            print "\\textbf\{".$line."\}\n";
            $line = "\\textbf\{".$line."\}";
        }        
        if ($line =~ m/\\caption/) {
            # $figures =~ s/caption\{(.*?\.[\s\}])/caption{\\textbf{\1} /msg;
            print "found a caption\n";
            $caption = 1;
        }
    }
    
    print OUTFILE $line;
}

close FILE;
close OUTFILE;

# foreach $i (0..$#figures)
# {
#     $figure = $figures[$i];
#     $tmp = `find . -follow -name $figure -print 2>/dev/null`;
#     @tmp = split("\n",$tmp);
#     $fullfigures[$i] = $tmp[0];
#     chomp($fullfigures[$i]);
# }

# $i=1;
# foreach $fullfigure (@fullfigures)
# {
#     ($figname = $fullfigure) =~ s/.*\///;
#     $figname =~ s/\.ps/.eps/;
#     if ($i<10) {$fignum = "0$i";} else {$fignum = $i;}
#     $prefix = "fig".$fignum."_";
#     `cp $fullfigure $folder/$prefix$figname`;
#     `ln -s ../$fullfigure $folder/.`;
#     $i = $i+1;
# }
