#!/usr/bin/env perl
# PURPOSE: for each mRNA declared in -gff, write to STDOUT its spliced CDS sequence taken from -fasta
# USAGE: gff2cdsfa.pl -fasta t/VV10.fa -gff t/t1.gff
# AUTHOR: malcolm.cook@gmail.com

use strict;
use warnings;

use Bio::SeqIO;
use Bio::DB::SeqFeature::Store;
use Bio::SeqFeatureI;
sub Bio::SeqFeatureI::spliced_cds {
    ## PURPOSE: for a (presumed) mRNA SeqFeature, create new Bio::Seq
    ## holding its spliced CDS, reverse complemented as needed, ided
    ## after the display_id, if available, else the load_id.
    my ($self,$db) = @_;
    my @CDS = $self->CDS;
    @CDS=sort {$a cmp $b} @CDS;    # don't depend upon input GFF lines being sorted 
    my $seqstr;
    $seqstr .= $_->seq->seq for @CDS;    
    my $seq = Bio::Seq->new( -id => $self->display_id || $self->load_id,
                 -seq => $seqstr);
    $self->strand eq "-" ? $seq->revcom : $seq;
}

my %ARGV=@ARGV;
$ARGV{-adaptor}||='memory';
my $db = Bio::DB::SeqFeature::Store->new(%ARGV) or die $!;

my $ss = $db->get_seq_stream(-type =>  'mRNA');
my $seqout = Bio::SeqIO->new(-format=>'fasta');
while (my $mRNA = $ss->next_seq) {
    $seqout->write_seq($mRNA->spliced_cds)
}

0;
