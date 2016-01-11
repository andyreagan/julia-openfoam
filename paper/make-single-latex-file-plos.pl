#!/usr/bin/perl
$usage = "make-single-latex-file.pl foo.tex

replaces \inputs{stuff} with stuff.tex,
inputs bbl file, etc.

improve input detection!!! (need to include file find)

(usually run on foo-revtex4.tex)

presumes that comments start with two percentage signs: %%
(which works nicely for tab indenting in emacs so you
should be doing it anyway)

creates foo-combined.tex
";
if ($#ARGV < 0)
{
    print $usage;
    exit;
}
($basefile = $ARGV[0]) =~ s/\.tex//;
($realbasefile = $basefile) =~ s/-plos//;
$file = "$basefile.tex";
($bblfile = $file) =~ s/tex$/bbl/;

$outfile = "$basefile-combined.tex";

$help = "Parses through to remove all indexing commands,
inputs all \inputs directly,
and inputs the .bbl file

does not include supplementary per PLoS madness
";

open (FILE,"$file") or die "can't open $file: $!\n";
open (OUTFILE,">$outfile") or die "can't open $outfile: $!\n";

undef $tex;
foreach $line (<FILE>)
{
   $line =~ s/%%.*$//;
    if ($line =~ m/\\currfile/) {
	$line = "";
    }
    
    $line =~ s/\\filenamebase/$realbasefile/;
    
    if ($line =~ m/\\input/)
    {
	$line =~ m/\\input\{(.*?)\}/;
#	print "$1\n";
	$inputfile = $1.".tex";
	open (INPUT,"$inputfile") or die "can't open $inputfile: $!\n";
	undef $/;
	$input = <INPUT>;
	$/ = "\n";
	close INPUT;

	$line =~ s/\\input\{(.*?)\}/$input/;

# one extra layer
	while ($line =~ m/\\input\{(.*?)\}/) {
#	    print "$1\n";
	    $inputfile = $1.".tex";
	    open (INPUT,"$inputfile") or die "can't open $inputfile: $!\n";
	    undef $/;
	    $input = <INPUT>;
	    $/ = "\n";
	    close INPUT;
	    $line =~ s/\\input\{(.*?)\}/$input/;
	}
    }
    $tex = $tex."$line";
}

# bbl file
open (BBL,"$bblfile") or die "can't open $bblfile: $!\n";
undef $/;
$bbl = <BBL>;
$/ = "\n";
close BBL;

$tex =~ s/\\bibliography\{.*?\}/$bbl/;

# convert pdf links to tiff
# $tex =~ s/\.pdf/\.tiff/g;

# $tex =~ s/\\newcommand\{\\nindex.*?$//msg;
# $tex =~ s/\\newcommand\{\\sindex.*?$//msg;
# $tex =~ s/\\nindex.*?$//msg;
# $tex =~ s/\\sindex.*?$//msg;

$tex =~ s/%%.*?\n//msg;

# one last sweep
# not working
# $tex =~ s/(^\\)%.*?\n//msg; 

$tex =~ s/\\%/PERCENTAGEARAMA/g;
$tex =~ s/%.*?\n//msg;
$tex =~ s/PERCENTAGEARAMA/\\%/g;

# clean out repeated blank lines
$tex =~ s/\n\n+/\n\n/msg;
$tex =~ s/^\s+//;
$tex =~ s/\s+$/\n/;

# remove supplementary
# $tex =~ s/\\newwrite\\tempfile.*\\end{document}/\\end{document}/ms;

print OUTFILE $tex;

close FILE;
close OUTFILE;


