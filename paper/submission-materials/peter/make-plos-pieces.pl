#!/usr/bin/perl
$usage = "make-plos-pieces.pl foo.tex

presumes that foo.tex has been latexed and bibtexed.

presumes that comments start with two percentage signs: %%
(which works nicely for tab indenting in emacs so you
should be doing it anyway)

creates numerous plos-foo*.tex pieces
";
if ($#ARGV < 0)
{
    print $usage;
    exit;
}
($basefile = $ARGV[0]) =~ s/\.tex//;
($realbasefile = $basefile) =~ s/-revtex4//;

$file = "$basefile.tex";
($bblfile = $file) =~ s/tex$/bbl/;

$outfile = "plos-$basefile.tex";
$outfile =~ s/-revtex4//; # if it's there

$help = "Parses through to remove all indexing commands,
inputs all \inputs directly, with one level of recursion
and inputs the .bbl file
";

open (FILE,"$file") or die "can't open $file: $!\n";
# open (OUTFILE,">$outfile") or die "can't open $outfile: $!\n";

undef $tex;
foreach $line (<FILE>)
{
    $line =~ s/%%.*$//;
    
    if ($line =~ m/\\currfile/) {
	$line = "";
    }

    $line =~ s/\\filenamebase/$realbasefile/;

    ## operate on files included into revtex4
    if ($line =~ m/\\input/) 
    {
	$line =~ m/\\input\{(.*?)\}/;
#	print "$1\n";
	$inputfile = $1.".tex";
	open (INPUT,"$inputfile") or die "can't open $inputfile: $!\n";
	$/ = "\n";
	$input = "";
	foreach $inputline (<INPUT>) {
	    if ($inputline =~ m/\\input\{(.*?)\}/) {
		$inputfile2 = $1.".tex";
		$inputfile2 =~ s/\\filenamebase/$realbasefile/;
		open (INPUTFILE2,"$inputfile2") or die "can't open $inputfile2: $!\n";
		undef $/;
		$inputline = <INPUTFILE2>;
		$/ = "\n";
	    }
	    $input = $input.$inputline;
	}

	close INPUT;



	$input =~ s/\\%/PERCENTAGEARAMA/g;
	$input =~ s/%.*?\n//msg;
	$input =~ s/PERCENTAGEARAMA/\\%/g;

# clean out repeated blank lines
	$input =~ s/\n\n+/\n\n/msg;
	$input =~ s/^\s+//;
	$input =~ s/\s+$/\n/;

# remove revtexonly lines
	$input =~ s/\\revtexonly.*?$//msg;
# make sure plainlatex commands end at a newline
	$input =~ s/\\plainlatexonly{(.*?)}$/\1/msg;

	## remove section numbers for plos one
	$input =~ s/section\{/section*{/g;


	# extract tables and figures
	# check if file is supplementary
	$i=1;
	if ($inputfile !~ m/\.supp/) {
	    foreach $fig ($input =~ m/\\begin{[sideways]{0,8}figure.*?\\end{.*?igure\*{0,1}}/msg) {
		$figures = $figures."\n\n\n".$fig;
		$i++;
		if ($i==10) { 
		    $figures = $figures."\n\\clearpage\n"; 
		    $i=0;
		}
#		$figures = $figures."\n\n\\clearpage\n".$fig;
	    }
	    # turn figures into .eps
	    $figures =~ s/\.pdf/\.eps/g;
	    foreach $tab ($input =~ m/\\begin{table.*?\\end{table\*{0,1}}/msg) {
		$tables = $tables."\n\n".$tab;
	    }
	} else { # supplementary
	    $i=1;
	    foreach $fig ($input =~ m/\\begin{[sideways]{0,8}figure.*?\\end{.*?igure\*{0,1}}/msg) {
		$suppfigures = $suppfigures."\n\n\n".$fig;
		$i++;
		if ($i==10) { 
		    $suppfigures = $suppfigures."\n\\clearpage\n"; 
		    $i=0;
		}
#		$suppfigures = $suppfigures."\n\n\\clearpage\n".$fig;
	    }
	    # turn figures into .eps
	    $suppfigures =~ s/\.pdf/\.eps/g;
	    foreach $tab ($input =~ m/\\begin{table.*?\\end{table\*{0,1}}/msg) {
		$supptables = $supptables."\n\n".$tab;
	    }
	}

	# now remove figures and tables
	$input =~ s/\\begin{[sideways]{0,8}figure.*?\\end{.*?igure\*{0,1}}//msg;
	$input =~ s/\\begin{tabl.*?\\end{table\*{0,1}}//msg;

	$inputfileplos = "plos-".$inputfile;
	open (INPUTFILEPLOS,">$inputfileplos") or die "can't open $inputfileplos: $!\n";
	print INPUTFILEPLOS $input;
	close INPUTFILEPLOS;
	print "$inputfileplos created...\n";

#	$line =~ s/\\input\{(.*?)\}/$input/;
    }
#    $tex = $tex."$line";
}



$figures =~ s/caption\{(.*?\.[\s\}])/caption{\\textbf{\1} /msg;

# comment out figures
$figures =~ s/\\includegraphics/%\\includegraphics/g;

# remove any direct hard coded references to figures
$figures =~ s/(\\includegraphics.*?)\{.*\/(.*?)\}/\1\{\2\}/g;

$tables =~ s/caption\{(.*?\.[\s\}])/caption{\\textbf{\1} /msg;

($figfile = $outfile) =~ s/\.tex/-figures.tex/;
open (FIGURES,">$figfile") or die "can't open $figfile: $!\n";
print FIGURES "$figures
\\clearpage";
print "$figfile created...\n";
close FIGURES;

($tabfile = $outfile) =~ s/\.tex/-tables.tex/;
open (TABLES,">$tabfile") or die "can't open $tabfile: $!\n";
print TABLES "$tables
\\clearpage";
print "$tabfile created...\n";
close TABLES;

# supplementary

$suppfigures =~ s/caption\{(.*?\.[\s\}])/caption{\\textbf{\1} /msg;

# comment out figures
$suppfigures =~ s/\\includegraphics/%\\includegraphics/g;

# remove any direct hard coded references to figures
$suppfigures =~ s/(\\includegraphics.*?)\{.*\/(.*?)\}/$1{$2}/g;

$supptables =~ s/caption\{(.*?\.[\s\}])/caption{\\textbf{\1} /msg;
# help tables resize
# beautiful
$supptables =~ s/\\begin{tabular}/\\resizebox{\\textwidth}{!}{\\begin{tabular}/msg;
$supptables =~ s/\\end{tabular}/\\end{tabular}}/msg;

# put tables on separate pages:
$supptables =~ s/\\end{table}/\\end{table}\\clearpage/msg;
$supptables =~ s/\\end{table\*}/\\end{table\*}\\clearpage/msg;

($suppfigfile = $outfile) =~ s/\.tex/-suppfigures.tex/;
open (SUPPFIGURES,">$suppfigfile") or die "can't open $suppfigfile: $!\n";
print SUPPFIGURES $suppfigures;
print "$suppfigfile created...\n";
close SUPPFIGURES;

($supptabfile = $outfile) =~ s/\.tex/-supptables.tex/;
open (SUPPTABLES,">$supptabfile") or die "can't open $supptabfile: $!\n";
print SUPPTABLES $supptables;
print "$supptabfile created...\n";
close SUPPTABLES;


# bbl file
# open (BBL,"$bblfile") or die "can't open $bblfile: $!\n";
# undef $/;
# $bbl = <BBL>;
# $/ = "\n";
# close BBL;

# $tex =~ s/\\bibliography\{.*?\}/$bbl/;

# $tex =~ s/\\newcommand\{\\nindex.*?$//msg;
# $tex =~ s/\\newcommand\{\\sindex.*?$//msg;
# $tex =~ s/\\nindex.*?$//msg;
# $tex =~ s/\\sindex.*?$//msg;

# $tex =~ s/%%.*?\n//msg;

# one last sweep
# not working
# $tex =~ s/(^\\)%.*?\n//msg; 


# print OUTFILE $tex;

close FILE;

# close OUTFILE;


