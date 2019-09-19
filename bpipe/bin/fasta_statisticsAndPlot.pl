#!/usr/bin/env perl
# Takes a fasta-file (usually a genomic or transcriptomic assembly) and checks for
# potential problems as well as calculates a few basic statistics.
# 
# NBIS 2018
# jacques.dainat@nbis.se

use warnings;
use strict;
use Statistics::R;
use POSIX qw(strftime);
use File::Basename;
use Pod::Usage;
use Getopt::Long;
use Try::Tiny;
use Carp;
use IO::File;
use Bio::SeqIO;

my $nb_seq = 0;
my $problemcount=0;
my $nbseq_withLowerCase=0;
my $total_lowerCaseCount;
my $Ncount=0;
my $pureNseq=0;
my $totalcount=0;
my $totalcountOver1000=0;
my $totalcountOver10000=0;
my $gccount=0;
my $total_noNs=0;
my @sequencelength=();
my @sequencelengthOver1000=();
my @sequencelengthOver10000=();

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



my $outstream = IO::File->new();;
my $outstreamFix;
my $outstreamError;

if ( defined($opt_dirRes) ) {
	if (-d $opt_dirRes) {
 		print "$opt_dirRes output directory already exits !\n";exit;
		}
	mkdir $opt_dirRes;

    $outstream->open( "$opt_dirRes/fasta_report.txt", 'w' ) or
        croak(sprintf( "Can not open '%s' for writing %s", "$opt_dirRes/fasta_report.txt", $! ));

    open(my $fh, '>', "$opt_dirRes/fix_sequences.fa") or die "Could not open file '$opt_dirRes/fix_sequences.txt' $!";
    $outstreamFix = Bio::SeqIO->new(-fh => $fh , -format => 'fasta');

    open(my $fhe, '>', "$opt_dirRes/problem_sequences.fa") or die "Could not open file '$opt_dirRes/problem_sequences.txt' $!";
    $outstreamError = Bio::SeqIO->new(-fh => $fhe , -format => 'fasta');
}
else {
    $outstream->fdopen( fileno(STDOUT), 'w' ) or
	        croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
    $outstreamError = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'fasta');
}

# INPUT FILE
my $inseq = Bio::SeqIO->new(-file   => "<$opt_infile", -format => 'fasta');


#Calculate the statistics from the entries in the hash
my $cp1kb=0;
my $cp10kb=0;
while( my $seqObj = $inseq->next_seq() ) {

  my $sequence = $seqObj->seq();
  $nb_seq++;
  my $problemNstart=undef;
  my $problemNend =undef;
  my $problem_lowerCase = undef;

  #save sequence length for N50 calculation
  push @sequencelength, length $sequence; 

  #Check for Ns at the beginning or end of sequence  
  if ($sequence =~ /^N/){
    print $outstreamError->write_seq( $seqObj );
    $problemNstart="yes";
    $problemcount++;
  }
  if ($sequence =~ /N$/){
    print $outstreamError->write_seq( $seqObj );
    $problemNend = "yes";
    $problemcount++;
  }
  if (length($sequence) > 1000){
  	$cp1kb++;
	#Count total size over 1000
	$totalcountOver1000 += length $sequence;
 	 #save sequence length for N50 calculation
  	push @sequencelengthOver1000, length $sequence; 

  	if (length($sequence) > 10000){
  		$cp10kb++;
		  #Count total size over 10000
  		  $totalcountOver10000 += length $sequence;
		  #save sequence length for N50 calculation
  		  push @sequencelengthOver10000, length $sequence; 
  	}
  }
  
  
  #Count number of NNN regions
  my $match=0;
  $match++ while $sequence =~ /[ACGT]N+[ACGT]/g;
  $Ncount += $match;
  #Count GC
  $gccount += ($sequence =~ tr/gGcC/gGcC/);
  #Count size total
  $totalcount += length $sequence;
  #Count size with No Ns
  my $noNs=$sequence;
  $noNs =~ s/N//g;
  $total_noNs += length $noNs;
  if(length $noNs == 0){
  	$pureNseq++;
  }
  #Count lowercase outside Ns
  my $lowerCaseCount += ($noNs =~ tr/atgc/atgc/);
  $total_lowerCaseCount += $lowerCaseCount;
  if($lowerCaseCount > 0){
	$nbseq_withLowerCase++;
	$problem_lowerCase="yes";	
  }
  # Now we can print the fixed sequnece
  if ($opt_dirRes){
	
 	#Fix the lowercase issue
	if($problem_lowerCase){
	  	$sequence =~ tr/atgcn/ATGCN/;
 		$seqObj->seq($sequence);
	}

	#Fix starting N
	if($problemNstart or $problemNend){
		$sequence =~ s/^N+|N+$//g;
		$seqObj->seq($sequence);
        }
	print $outstreamFix->write_seq( $seqObj );
  }
}

#Calculate some statistics
my $GCpercentage = ($gccount/$totalcount*100);
my $GCnoNs = ($gccount/$total_noNs*100);
my $totalNs =$totalcount-$total_noNs;

#################
# Calculate N50 #
@sequencelength = reverse sort { $a <=> $b } @sequencelength;
my $N50=$totalcount/2;
my $sum=0;
my $entry;
# copy of sequencelength to keep it intactfor R calculation purpose later
my @sequencelengthForN50Calcul=@sequencelength;

my $nbcontig=0;

while ($sum < $N50){
  $entry = shift @sequencelengthForN50Calcul;
  $sum += $entry;
  $nbcontig+=1;
}

#################
# Calculate N90 #
# copy of sequencelength to keep it intactfor R calculation purpose later
my @sequencelengthForN90Calcul=@sequencelength;
@sequencelengthForN90Calcul = reverse sort { $a <=> $b } @sequencelengthForN90Calcul;
my $NinetyPercGenomeSize=( 90*$totalcount / 100);
$sum=0;
my $N90;
while ($sum < $NinetyPercGenomeSize){
  $N90 = shift @sequencelengthForN90Calcul;
  $sum += $N90;
} 


###########################
# Calculate N50 over 1000 #
@sequencelengthOver1000 = reverse sort { $a <=> $b } @sequencelengthOver1000;
my $HalfGenomeSizeOver1000=$totalcountOver1000/2;
$sum=0;
my $N50over1000=0;
while ($sum < $HalfGenomeSizeOver1000){
  $N50over1000 = shift @sequencelengthOver1000;
  $sum += $N50over1000;
} 


############################
# Calculate N50 over 10000 #
@sequencelengthOver10000 = reverse sort { $a <=> $b } @sequencelengthOver10000;
my $HalfGenomeSizeOver10000=$totalcountOver10000/2;
$sum=0;
my $N50over10000=0;
while ($sum < $HalfGenomeSizeOver10000){
  $N50over10000 = shift @sequencelengthOver10000;
  $sum += $N50over10000;
}


###########################
#print out the statistics #
my $date = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
my $StingToPrint;
$StingToPrint .= "\n========================================\n";
$StingToPrint .= "Fasta-statistics\ launched the $date:\n";
$StingToPrint .= "There are $nb_seq sequences\n";
$StingToPrint .= "There are $cp10kb sequences > 10kb \n";
$StingToPrint .= "There are $cp1kb sequences > 1kb \n";
$StingToPrint .= "There are $totalcount nucleotides, of which $totalNs are Ns\n";
$StingToPrint .= "There are $Ncount N-regions (possibly links between contigs)\n";
$StingToPrint .= "There are $pureNseq pure (only) N sequences. Assembler doing that must be notified ! \n";
$StingToPrint .= "There are $problemcount sequences that begin or end with Ns (see problem_sequences.txt)\n";
$StingToPrint .= sprintf("The GC-content is %.1f",$GCpercentage);
$StingToPrint .= "\%";
$StingToPrint .= sprintf(" (not counting Ns %.1f", $GCnoNs);
$StingToPrint .= "\%)\n";
$StingToPrint .= "There are $nbseq_withLowerCase sequences with lowercase nucleotides (Ns not considered)\n";
$StingToPrint .= "There are $total_lowerCaseCount lowercase nucleotides (Ns not considered)\n";
$StingToPrint .= "The N50 is $entry\n";
$StingToPrint .= "The N90 is $N90\n";
$StingToPrint .= "The N50 for sequences over 1000bp is $N50over1000\n";
$StingToPrint .= "The N50 for sequeces over 10000bp is $N50over10000\n";
$StingToPrint .= "========================================\n";
if($opt_dirRes){
	$StingToPrint .= "You will find corrected sequences in fix_sequences.txt (Extremities N trimed and lowercase nucleotides).\n";
}
print $outstream "$StingToPrint";

#######
#
# Plot
#
######
if($opt_dirRes){
		print $StingToPrint;
	if($nb_seq > 1){ #If only 1 seq we get error like: density.default(listValues) : need at least 2 points to select a bandwidth automatically
        # temporary file name
        my $tempFile1="dump.tmp";

		try {
			print $StingToPrint;
			print "This result was saved in the $opt_dirRes directory.\nThe plots are in <pdf> format and available in the directory.\n";

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

		## info ##
		#myhist$breaks contient la valeur minimale de chaque intervalle
		#myhist$mids contient la valeur au milieu de chaque intervalle
		#myhist$counts contient le nombre de valeurs situées dans cet intervalle
		#myhist$density contient la proportion de valeurs situées dans cet intervalle (autrement dit, tec.hist$counts / length (tec)).

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
		catch{
			warn "caught error: $_";
			unlink $tempFile1;
		}
	}
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
