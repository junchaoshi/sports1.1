#!/usr/bin/perl

use strict;
use warnings;

my $usage = <<"USAGE";

Description: This script is used to filter reads out of selected length range.

Usage: filter_length.pl in_file.fa min_length max_length > out_file.fa

USAGE

my $in_file = shift or die;

open IN_FILE, $in_file
	or die "Can't open '$in_file': $!";

my $min_len = shift or die;
my $max_len = shift or die;

my $id;
my $seq;
my $read;

while (<IN_FILE>){
	chomp;
	if (/^>/){
		($id, $read) = (split /\s+/)[0, 1];
		next;
	}
	$seq = $_;
	if ((length $seq) >= $min_len && (length $seq) <= $max_len){
		print "$id\t$read\n$seq\n";
	}
}


