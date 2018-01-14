#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';

my $usage = <<"USAGE";
Description:	This script adds CCA end to original genomic tRNA 3’ end and add G to Histidine tRNA 5’ end.
		input_tRNA_file is in .fa format.
Usage: tRNA_db_processing.pl input_tRNA_file.fa
USAGE

my $in_file = shift or die $usage;

my @input;
my $input_name;
my $input_address = abs_path($in_file);
my $input_suffix;

@input = split(/\//, $input_address);
$input_name = pop (@input);
$input_address = join('/', @input) . '/';
@input = split(/\./, $input_name);
$input_suffix = pop (@input);
$input_name = join('.', @input);

my $out_file = $input_address . $input_name . "_CCA.fa";

open INPUT, $in_file
	or die "Can't open '$in_file': $!";
open OUTPUT, ">$out_file"
	or die "Can't open '$out_file': $!";

my @seqs = ();
my @annos;
my $end = "CCA";



while (<INPUT>){
	chomp;
	if($_ =~ /^>(.*)His(.*)/){
		if(@seqs){
			print OUTPUT join("", @seqs);
			print OUTPUT "$end\n";
			@seqs = ();
		}
		@annos = split(/ /, $_);
		print OUTPUT join(" ", @annos);
		print OUTPUT "\nG";
		next;
	}
	elsif($_ =~ /^>(.*)/ && $_ !~ /^>(.*)His(.*)/ ){
		if(@seqs){
			print OUTPUT join("", @seqs);
			print OUTPUT "$end\n";
			@seqs = ();
		}
		@annos = split(/ /, $_);
		print OUTPUT join(" ", @annos);
		print OUTPUT "\n";
		next;
	}else{
		push @seqs, $_;
	}
}
print OUTPUT join("", @seqs);
print OUTPUT "$end\n";
