#!/usr/bin/perl

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
		
	[--id_trunk ]
		6-character base for all IDs (e.g. HOMSAP for Homo sapiens)
	[--clean ]
		Set to anything to remove Maker quality scores from GFF file	
	[--ccount ]
		Set the cdna counter to something other than 0 
	[--gcount ]
		Set the gene counter to something other than 0
	[--exount ]
		Set the exon counter to something other than 0
	[--tcount ]
		Set the transcript counter to something other than 0

  Ouput:    
    [--outfile filename]
        The name of the output file. 
};

my $outfile = undef;
my $gff = undef;

my $gcount = undef; 
my $tcount = undef;
my $ecount = undef;
my $ccount = undef;

my $clean = undef;
my $id_trunk = undef;
my $current_gene = undef;
my %transcript_hash ;

my $help;

GetOptions(
    "help" => \$help,
    "gff=s" => \$gff,
    "gcount=i" => \$gcount,
    "tcount=i" => \$tcount,
    "ecount=i" => \$ecount,
    "ccount=i" => \$ccount,
    "clean=s" => \$clean,
    "id_trunk=s" => \$id_trunk,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

### Create ID trunk

my $gcounter = $gcount || 100000;
my $tcounter = $tcount || 100000;
my $ecounter = $ecount || 100000;
my $ccounter = $ccount || 100000;
my $ucounter = 100000;

### Parse GFF input file and add annotations

open(GFF, "<$gff") || die("Can't open $gff.");

my $gffio = Bio::Tools::GFF->new(-file => $gff, -gff_version => 3);
my $gffout = Bio::Tools::GFF->new(-gff_version => 3);

while( my $feature = $gffio->next_feature()) {
	
	my @values = $feature->get_tag_values('ID');
	my $old_id = shift @values;
	
	
	my @parents = undef;
	unless ($feature->primary_tag =~/gene/) {
		@parents = $feature->get_tag_values('Parent');
	}

	if ($feature->primary_tag =~ /gene/) {
		$gcounter += 1;
		my $new_id = $id_trunk . "G" . $gcounter ;
		$current_gene = $new_id;
		$feature->remove_tag('ID');
		$feature->add_tag_value('ID',$new_id);
		
		if ($feature->has_tag('Name')) {
                	$feature->remove_tag('Name');
		}
                
		$feature->add_tag_value('Name',$new_id);

		$feature->add_tag_value("OId",$old_id);
	
	} elsif ($feature->primary_tag =~ /mRNA/ or $feature->primary_tag =~ /RNA/) { 
		$tcounter += 1;
		my $new_id = $id_trunk . "T" . $tcounter ;

		$transcript_hash{$old_id} = $new_id;

		$feature->remove_tag('Parent');
		$feature->add_tag_value('Parent',$current_gene);

		$feature->remove_tag('ID');
		$feature->add_tag_value('ID',$new_id);

		if ($feature->has_tag('Name')) {
                        $feature->remove_tag('Name');
                }

                $feature->add_tag_value('Name',$new_id);
		$feature->add_tag_value('OId',$old_id);

		if (defined $clean) {
			$feature->remove_tag('_AED') if ($feature->has_tag('_AED'));
			$feature->remove_tag('_eAED') if ($feature->has_tag('_eAED'));
			$feature->remove_tag('_QI') if ($feature->has_tag('_QI'));
		}



	} elsif ($feature->primary_tag =~ /exon/) {
		$ecounter += 1;
                my $new_id = $id_trunk . "E" . $ecounter ;
		
		# Parents can be > 1, so we need to look up all new IDs
		my @new_parents = ();

		foreach my $parent (@parents) {
			
			my $this_mapped_id = $transcript_hash{$parent} ;
			push(@new_parents,$this_mapped_id);
		
		}

                $feature->remove_tag('Parent');
                $feature->add_tag_value('Parent',@new_parents);

                $feature->remove_tag('ID');
                $feature->add_tag_value('ID',$new_id);

                if ($feature->has_tag('Name')) {
                        $feature->remove_tag('Name');
                }

                $feature->add_tag_value('Name',$new_id);

	} elsif ($feature->primary_tag =~ /CDS/ ) {

		$ccounter += 1;
                my $new_id = $id_trunk . "C" . $ccounter ;
		
		my $parent = shift @parents;
	
		my $current_transcript = $transcript_hash{$parent};

                $feature->remove_tag('Parent');
                $feature->add_tag_value('Parent',$current_transcript);

                $feature->remove_tag('ID');
                $feature->add_tag_value('ID',$new_id);

                if ($feature->has_tag('Name')) {
                        $feature->remove_tag('Name');
                }

                $feature->add_tag_value('Name',$new_id);	
	} elsif ($feature->primary_tag =~ /utr/ or $feature->primary_tag =~ /UTR/ ) {

		$ucounter += 1;
		my $new_id = $id_trunk . "U" . $ucounter ; 

		my $parent = shift @parents;

                my $current_transcript = $transcript_hash{$parent};

		$feature->remove_tag('ID');
                $feature->add_tag_value('ID',$new_id);


                $feature->remove_tag('Parent');
                $feature->add_tag_value('Parent',$current_transcript);
	
        	if ($feature->has_tag('Name')) {
                        $feature->remove_tag('Name');
                }

		$feature->add_tag_value('OId',$old_id);

	}

	
	print $feature->gff_string($gffout) , "\n";	
}

$gffio->close();





