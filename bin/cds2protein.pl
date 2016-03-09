#!/usr/bin/perl

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Data::Dumper;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the CDS file to read
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
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

my $fasta_in = Bio::SeqIO->new(-file => $infile, -format => "fasta" ) or die "Couldn't load the input file";
my $fasta_out = Bio::SeqIO->new(-file => '>'.$outfile, -format => "fasta" ) or die "Couldn't create the output file";


while (my $sequence = $fasta_in->next_seq){   
    # print the sequence   
    next if ($sequence->translate(-complete => 1) =~ /\*/ );
    $fasta_out->write_seq($sequence->translate(-complete => 1));
}

$fasta_in->close;
$fasta_out->close;

# --------------



