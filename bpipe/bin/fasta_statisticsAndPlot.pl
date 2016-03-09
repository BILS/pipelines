#!/usr/bin/perl
# Takes a fasta-file (usually a genomic or transcriptomic assembly) and checks for
# potential problems as well as calculates a few basic statistics.
# 
# By Henrik Lantz, BILS/Uppsala University, Sweden
# Modified by Jacques Dainat to plot contig distribution and N50

# usage: perl fasta_statistics fastafile [option]
#Not use# option: integer for Breaklines in R (default=1000)

use warnings;
use strict;
use Statistics::R;
use POSIX qw(strftime);
use File::Basename;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
        The name of the genome assembly to read. 
  Ouput:
    [--outfile filename]
        The name of the output file(s).
};

my $infile = undef;
my $outfile="fasta_report.txt";
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

my $header;
my %sequence=();
my $problemcount=0;
my $Ncount=0;
my $totalcount=0;
my $gccount=0;
my $total_noNs=0;
my @sequencelength=();
(my $ouputPlot = $outfile) =~ s/\.[^.]+$//;
#my $ouputPlot=basename($outfile);

if (! -f $infile) {
  print "$infile file does not exit !";exit;
}

open FASTA, $infile or die "Couldn't open fasta-file";
open (OUTFILE,">$outfile");
open (PROBLEMS,">$outfile.problems");

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
    print PROBLEMS "\>$key\n$sequence{$key}\n";
    $problemcount++;
  }
  if ($sequence{$key} =~ /N$/){
    print PROBLEMS "\>$key\n$sequence{$key}\n";
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
print OUTFILE "\n========================================\n";
print OUTFILE "Fasta-statistics\ launched the $date:\n";
print OUTFILE "There are ", scalar keys %sequence, " sequences\n";
print OUTFILE "There are $totalcount nucleotides, of which $totalNs are Ns\n";
print OUTFILE "There are $Ncount N-regions (possibly links between contigs)\n";
print OUTFILE "There are $problemcount sequence(s) that begin or end with Ns (see problem_sequences.txt)\n";
print OUTFILE "The GC-content is ";
printf OUTFILE ("%.1f", $GCpercentage);
print OUTFILE "\% (not counting Ns ";
printf OUTFILE ("%.1f", $GCnoNs);
print OUTFILE "\%)\n";
print OUTFILE "The N50 is $entry\n";
print OUTFILE "========================================\n";

close OUTFILE;
close PROBLEMS;

#Display result on screen:
open RESULT, $outfile or die "Couldn't open fasta-file";
while (<RESULT>){
print "$_";
}
print "The results is also availble in the file $outfile.\nThe plots are in <pdf> format and available in the directory.\n";

# temporary file name
my $tempFile1="dump.tmp";

# write the data in temporary file
open(FILE, ">$tempFile1") || die "Erreur E/S:$!\n";
foreach my $size ( @sequencelength ) {
  print FILE "$size\n";
}
close(FILE);

#######
#
# Plot
#
######

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
