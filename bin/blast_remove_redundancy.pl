#!/usr/bin/perl

use strict;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the blast file to read (blast+, outfmt 6 required)
		
	[--cutoff score]
		Identity score above which to report hits as redundant (default: 80.0)
		
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my $cutoff = 80.0 ;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
	"cufoff=i" => \$cutoff,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}


open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

# Streaming the file, line by line
while (<$IN>) {
	chomp; 
	my $line = $_; 

	my @fields = split("\t",$line) ;
	
	die "Does not look like a tab-formatted blast+ output:\n$line\n" if (scalar @fields != 12);
	
	my ($query,$target,$ident) = @fields[0..2] ;
	
	next unless ($cutoff <= $ident) and ($target ne $query);
	
	print STDOUT $target , "\n" ;

}

close ($IN);
close(STDOUT);


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


