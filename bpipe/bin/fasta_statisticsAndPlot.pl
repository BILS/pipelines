#!/usr/bin/perl
# Takes a fasta-file (usually a genomic or transcriptomic assembly) and checks for
# potential problems as well as calculates a few basic statistics.
# 
# By Henrik Lantz, BILS/Uppsala University, Sweden
# Modified by Jacques Dainat to plot contig distribution and N50

use warnings;
use strict;
use Statistics::R;
use POSIX qw(strftime);
use File::Basename;
use Pod::Usage;
use Getopt::Long;
use Carp;
use IO::File;

my $header;
my %sequence=();
my $problemcount=0;
my $Ncount=0;
my $totalcount=0;
my $gccount=0;
my $total_noNs=0;
my @sequencelength=();

my $opt_infile;
my $opt_dirRes;
my $opt_help = 0;


if ( !GetOptions( 'f|infile=s' => \$opt_infile,
                  'o|out|output=s' => \$opt_dirRes,
                  'h|help!'         => \$opt_help ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 0 } );
}

if ( !( defined($opt_infile) ) ) {
    pod2usage( {
           -message => "Must specify at leat one parameter:\nfasta_statisticsAndPlot.pl -f InputFastaFile [-o Ouput_directory] ",
           -verbose => 0,
           -exitval => 2 } );
}



my $outstream = IO::File->new();
my $outstreamError = IO::File->new();

if ( defined($opt_dirRes) ) {
	if (-d $opt_dirRes) {
 		print "$opt_dirRes output directory already exits !\n";exit;
		}
	mkdir $opt_dirRes;

    $outstream->open( $opt_dirRes."/fasta_report.txt", 'w' ) or
      croak( sprintf( "Can not open '%s' for writing %s", $opt_dirRes, $! ) );
    $outstreamError->open( $opt_dirRes."/problem_sequences.txt", 'w' ) or
      croak( sprintf( "Can not open '%s' for writing %s", $opt_dirRes, $! ) );
}
else {
    $outstream->fdopen( fileno(STDOUT), 'w' ) or
      croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
    $outstreamError->fdopen( fileno(STDOUT), 'w' ) or
      croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
}

my $output=$outstream;
my $outputProblem=$outstreamError;


open FASTA, $opt_infile or die "Couldn't open fasta-file";

#Populate a hash with the fasta-data

while (<FASTA>) {
  chomp;
  if (/^>(.*)$/){
    $header=$1;
  }
  elsif (/^(\S+)$/){
    $sequence{$header} .= $1 if $header;
  }
}
close FASTA;

#Calculate the statistics from the entries in the hash
foreach my $key (keys %sequence){ 
  push @sequencelength, length $sequence{$key}; #save sequence length for N50 calculation

  #Check for Ns at the beginning or end of sequence  
  if ($sequence{$key} =~ /^N/){
    print $outputProblem "\>$key\n$sequence{$key}\n";
    $problemcount++;
  }
  if ($sequence{$key} =~ /N$/){
    print $outputProblem "\>$key\n$sequence{$key}\n";
    $problemcount++;
  }
  
  
  #Count number of NNN regions
  my $match=0;
  $match++ while $sequence{$key} =~ /[ACGT]N+[ACGT]/g;
  $Ncount += $match;
  #Count GC
  $gccount += ($sequence{$key} =~ tr/gGcC/gGcC/);
  $totalcount += length $sequence{$key};
  my $noNs=$sequence{$key};
  $noNs =~ s/N//g;
  $total_noNs += length $noNs;
}

#Calculate some statistics
my $GCpercentage = ($gccount/$totalcount*100);
my $GCnoNs = ($gccount/$total_noNs*100);
my $totalNs =$totalcount-$total_noNs;
@sequencelength = reverse sort { $a <=> $b } @sequencelength;
my $N50=$totalcount/2;
my $sum=0;
my $entry;

my @sequencelengthForN50Calcul=@sequencelength;


my $nbcontig=0;
while ($sum < $N50){
  $entry = shift @sequencelengthForN50Calcul;
  $sum += $entry;
  $nbcontig+=1;
} 

#print out the statistics 
my $date = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
my $StingToPrint;
$StingToPrint .= "\n========================================\n";
$StingToPrint .= "Fasta-statistics\ launched the $date:\n";
$StingToPrint .= sprintf("There are ".scalar keys %sequence);
$StingToPrint .= " sequences\n";
$StingToPrint .= "There are $totalcount nucleotides, of which $totalNs are Ns\n";
$StingToPrint .= "There are $Ncount N-regions (possibly links between contigs)\n";
$StingToPrint .= "There are $problemcount sequence(s) that begin or end with Ns (see problem_sequences.txt)\n";
$StingToPrint .= sprintf("The GC-content is %.1f",$GCpercentage);
$StingToPrint .= "\%";
$StingToPrint .= sprintf(" (not counting Ns %.1f", $GCnoNs);
$StingToPrint .= "\%)\n";
$StingToPrint .= "The N50 is $entry\n";
$StingToPrint .= "========================================\n";
print $outstream "$StingToPrint";

#######
#
# Plot
#
######
if($opt_dirRes){
	
	print $StingToPrint;
	print "This result was saved in the $opt_dirRes directory.\nThe plots are in <pdf> format and available in the directory.\n";

	# temporary file name
	my $tempFile1="dump.tmp";

	# write the data in temporary file
	open(FILE, ">$tempFile1") || die "Erreur E/S:$!\n";
	foreach my $size ( @sequencelength ) {
	  print FILE "$size\n";
	}
	close(FILE);



	my $ouputName=basename($opt_infile);
	my $ouputPlot=$opt_dirRes."/".$ouputName;



	# Calcul percentage of contig right and left to N50
	my $percentContigRightN50=($nbcontig*100)/($#sequencelength+1);
	my $percentContigLeftN50=100-$percentContigRightN50;
	$percentContigRightN50=sprintf ("%0.2f",$percentContigRightN50)."%";
	$percentContigLeftN50=sprintf ("%0.2f",$percentContigLeftN50)."%";
	# Name of different outputs
	my $outputPlotLog=$ouputPlot."_PlotLog.pdf";
	my $outputPlotDensity=$ouputPlot."_PlotDensity.pdf";
	my $outputPlotHist=$ouputPlot."_PlotHisto.pdf";
	# calcul right and left position to write percentContigRightN50 and percentContigLeftN50
	my $biggestValue=shift @sequencelength;
	my $positionright=(5*$biggestValue)/100+$entry;
	my $positionleft=$entry-(5*$biggestValue)/100;
	# Tab=as.matrix(read.table("$tempFile1", sep="\t", he=T)) 
	# R object Declaration
	my $R = Statistics::R->new() or die "Problem with R : $!\n";

	# R command
	 $R->send(
	     qq`
	  listValues=as.matrix(read.table("$tempFile1", sep="\t", he=F))
	  myhist<-hist(listValues)
	  legendToDisplay=paste("Number of value used : ",length(listValues))

	  pdf("$outputPlotLog")
	  plot(x = log(myhist\$mids), y = log(myhist\$counts), xlab="log(Contig size)", ylab="log(Frequency)", main="Size distribution of contigs")
	  abline(v =log($entry), col=2)
	  axisValues=par("usr")
	  ymax=(axisValues[4]*90)/100
	  ymax2=(axisValues[4]*85)/100
	  shiftFivePercent=(5*(axisValues[2]-axisValues[1]))/100
	  text(x = log($entry)+shiftFivePercent, y = ymax2, paste("$percentContigRightN50"), cex = 1, col = "red")
	  text(x = log($entry)-shiftFivePercent, y = ymax2, paste("$percentContigLeftN50"), cex = 1, col = "red")
	  text(x = log($entry), y = ymax, paste("N50"), cex = 1, col = "red")
	  legend("topright", col=(1), lty=1, c(legendToDisplay))
	  dev.off()

	  pdf("$outputPlotDensity")
	  plot(density(listValues), xlab="Contig size", main="Size distribution of contigs")
	  abline(v =$entry, col=2)
	  axisValues=par("usr")
	  ymax=(axisValues[4]*90)/100
	  ymax2=(axisValues[4]*85)/100
	  text(x = $positionright, y = ymax2, paste("$percentContigRightN50"), cex = 1, col = "red")
	  text(x = $positionleft, y = ymax2, paste("$percentContigLeftN50"), cex = 1, col = "red")
	  text(x = $entry, y = ymax, paste("N50"), cex = 1, col = "red")  
	  legend("topright", col=(1), lty=1, c(legendToDisplay))
	  dev.off()
	  
	  pdf("$outputPlotHist")
	  hist(listValues, xlab="Contig size",main="Size distribution of contigs")
	  abline(v =$entry, col=2)
	  axisValues=par("usr")
	  ymax=(axisValues[4]*90)/100
	  ymax2=(axisValues[4]*85)/100
	  text(x = $positionright, y = ymax2, paste("$percentContigRightN50"), cex = 1, col = "red")
	  text(x = $positionleft, y = ymax2, paste("$percentContigLeftN50"), cex = 1, col = "red")
	  text(x = $entry, y = ymax, paste("N50"), cex = 1, col = "red")
	  legend("topright", col=(1), lty=1, c(legendToDisplay))
	  dev.off()`
	     );

	# Close the bridge
	$R->stopR();

	# Delete temporary file
	unlink $tempFile1;
	unlink "Rplots.pdf"; #created by " myhist<-hist(listValues)". I do not know how do differently...
}

__END__

=head1 NAME

fasta_statisticsAndPlot.pl - get some basic statistics about a nucleotide fasta file. (Number of sequence, Number of nucleotide, N50, GC-content,etc). It will also create R plots about contig size distribution. 
The R output plot will be perform only if an output is given.
This script is not yet designed for AA sequences or IUPAC Nucleotides.

=head1 SYNOPSIS

    ./fasta_statisticsAndPlot.pl --f=infile [--output=Directory]
    ./fasta_statisticsAndPlot.pl --help

=head1 OPTIONS

=over 8

=item B<--f>, B<--infile> or B<-f>

Input fasta file containing DNA sequences.

=item  B<--out>, B<--output> or B<-o>

[OPTIONAL] Output directory where diffrent output files will be written. If no output is specified, the result will written to STDOUT.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
