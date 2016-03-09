#!/usr/bin/ruby

require 'rubygems'
require 'bio'

fasta = Bio::FastaFormat.open(ARGV.shift)

fasta.each_entry do |entry|

	new_seq = entry.naseq.gsub(/\.$/, '')

	# Filter illegal sequences
	next if new_seq.include?(".")

	seq = Bio::Sequence.new(new_seq)

	puts seq.output_fasta(entry.definition,80)
	

end




	
