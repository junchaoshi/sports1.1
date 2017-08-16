#!/usr/bin/perl

use strict;
use warnings;

my $usage = <<"USAGE";

Description: This script combine reads in the fasta file to get unique sequence and its reads number. output format:

>t00000001 1234567
TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGC

't00000001' represents the abundance ranking among all the sequences. In this case, the abundance of this sequence is the highest. 
'1234567' represents the reads number of the sequence 'TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGC'.

Usage: combine_reads.pl in_file.fa > out_file.fa

USAGE

my $in_file = shift or die;

open IN_FILE, $in_file
	or die "Can't open '$in_file': $!";

my $seq;
my %reads;
my $num = 1;

while (<IN_FILE>){
	chomp;
	if (/^>/){
		next;
	}
	$seq = $_;
	if($reads{$seq}){
		$reads{$seq} += 1;
		next;
	}
	$reads{$seq} = 1;
}

foreach my $key (sort  { $reads{$b} <=> $reads{$a} } keys %reads){
	$key=~tr/[acgtun\.]/[ACGTTNN]/;
	printf ">t%08d\t$reads{$key}\n$key\n", $num;
	$num += 1;
}

