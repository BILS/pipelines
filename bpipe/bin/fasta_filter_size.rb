#!/usr/bin/env ruby
# == NAME
# fasta_filter_size.rb
#
# == USAGE
#  ./fasta_filter_size.rb [ -h | --help ]
#                    [ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# Parses a fasta-formatted sequence file and returns only sequences >= -s (size)
#
# == OPTIONS
#  -h,--help::                  Show help
#  -i,--infile=INFILE::         input file
#  -o,--outfile=OUTFILE::       output file

#
# == EXPERT OPTIONS
#
# == AUTHOR
#   Marc Hoeppner, mphoeppner@gmail.com

require 'rubygems'
require 'rdoc/usage'
require 'optparse'
require 'ostruct'
require 'logger'
require 'bio'


### Define modules and classes here

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.on("-h","--help","Display the usage information") {RDoc::usage}
opts.on("-i","--infile", "=INFILE","Input") {|argument| options.infile = argument }
opts.on("-s","--size", "=SIZE","Size") {|argument| options.size = argument }
opts.on("-o","--outfile", "=OUTFILE","Output") {|argument| options.outfile = argument }

opts.parse! rescue RDoc::usage('usage')

options.outfile ? output_stream = File.new(options.outfile,'w') : output_stream = $stdout

entries = Bio::FastaFormat.open(options.infile)

entries.each_entry do |entry|
	next unless entry.naseq.length > options.size.to_i
	
	output_stream.puts entry.to_seq.output(:fasta, :header => entry.definition , :width => 80)
	
end

output_stream.close




