#!/bin/bash
# Reads a cufflinks-produced protein fasta file(gffread -y) and removes stop codons

INPUT=$1

sed 's/\.$//g' $INPUT
