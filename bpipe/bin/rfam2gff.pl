#!/usr/bin/perl

use strict;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;


my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the file to read. 
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my @annotations;
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

my $gff_formatter = Bio::Tools::GFF->new(-gff_version => 3);

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

while (<$IN>) {
	chomp;
        my $line = $_;
        next if ($line =~ /^#.*$/); # Skipping comment lines

        my $annotation = parse_line($line);
        push(@annotations,$annotation);

}

open (my $OUT, '>', $outfile) or die "FATAL: Can't open file: $outfile for reading.\n$!\n";

foreach my $feature (@annotations) {
        $feature->gff_format($gff_formatter);
        print $OUT $feature->gff_string, "\n";
}

close($OUT);

# --------------------

sub parse_line {
        # chomp;
        my $line = shift ;
        my ($tn,$tacc,$qn,$qacc,$mdl,$mdlf,$mdlt,$seqf,$seqt,$strand,$trunc,$pass,$gc,$bias,$score,$evalue,$inc,$desc) = split(/\s+/,$line);

        my %tags = ( 'rfam-id' => $qn,
                         'rfam-acc' => ($qacc || 'unknown'),
                         'model_start' => $mdlf,
                         'model_end' => $mdlf,
                         'gc-content' => $gc,
                         'ID' => $qacc . "_" .  $tn . "_" . $seqf,
                         'Name' => $qacc . "_" .  $tn . "_" . $seqf,
            );

                my($from,$to) = sort($seqf,$seqt); # cmsearch reports coordinates in orientation of annotation, not chromosome. Need to sort from low to high for gff

            if( $evalue =~ /[0-9]/ ) {
                $tags{'evalue'} = $evalue;
            }

            my $f = Bio::SeqFeature::Generic->new( -seq_id => $tn,
                                                   -start => $from,
                                                   -end => $to,
                                                   -strand => $strand,
                                                   -frame => 0,
                                                   -primary_tag => 'ncRNA',     # may argue over whether this is an exon feature, but anything else will be ignored by Maker
                                                   -source_tag => 'Rfam',
                                                   -score => $score,
                                                   -tag => \%tags,
                );

                return $f;
}



# --------------

sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
}

sub runcmd {
  msg("Running:", @_);
  system(@_)==0 or err("Could not run command:", @_);
}

sub err {
  msg(@_);
  exit(2);
}


