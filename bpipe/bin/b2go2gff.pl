#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Bio::Tools::GFF;


my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--gff filename]
		The name of the gff file to read. 
		
	[--b2go filename]
		The name of the Blast2Go annotation file to read
		
  Ouput:    
    [--outfile filename]
        The name of the output file. 
};

my $outfile = undef;
my $gff = undef;
my $b2go = undef;
my $help;

GetOptions(
    "help" => \$help,
    "gff=s" => \$gff,
	"b2go=s" => \$b2go,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}


#### Read the annotation file

my %lookup = read_annotation_file($b2go);


### Parse GFF input file and add annotations

open(GFF, "<$gff") || die("Can't open $gff.");

my $gffio = Bio::Tools::GFF->new(-file => $gff, -gff_version => 3);
my $gffout = Bio::Tools::GFF->new(-fh => $outfile, -gff_version => 3);


while( my $feature = $gffio->next_feature()) {
	
	if ($feature->primary_tag =~ /mRNA/) { 
	
		my @values = $feature->get_tag_values('ID');
		my $id = shift @values;
	
		my $annotation = $lookup{$id} ;
	
		if (defined $annotation) {
			$feature->remove_tag('Description') if ($feature->has_tag('Description'));
			$feature->add_tag_value('Description','Predicted: ' . $annotation);
			print $feature->gff_string($gffout) , "\n";	
		} else { 
			print $feature->gff_string($gffout) , "\n";
		}
	
	} else {
		print $feature->gff_string($gffout) , "\n";
	}
	
}

$gffio->close();



sub read_annotation_file(file) {
	
	my $file = shift;
	my %lookup = {};
		
	open (my $IN, '<', $file) or die "FATAL: Can't open file: $file for reading.\n$!\n";

	while (<$IN>) {
		chomp; 
		my $line = $_; 

		my ($id,$go,$name) = split("\t", $line);
		$lookup{$id} = $name ;	
	}

	close ($IN);
	
	return %lookup;
	
}



