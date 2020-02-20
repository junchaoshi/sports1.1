#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';

my $usage = <<"USAGE";
Description:	This script transfers mitotRNAdb tRNA sequences to standard tRNA sequences. 
		input_tRNA_file is in .fa format.
Usage: mt_tRNA_transfer.pl input_tRNA_file output_file_name
USAGE

my $in_file = shift or die $usage;

my $file_name = shift or die $usage;


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

my $out_file_1 = $input_address . $file_name . ".fa";

my $out_file_2 = $input_address . $file_name . "_CCA.fa";


open INPUT, $in_file
	or die "Can't open '$in_file': $!";
open OUTPUT_1, ">$out_file_1"
	or die "Can't open '$out_file_1': $!";
open OUTPUT_2, ">$out_file_2"
	or die "Can't open '$out_file_2': $!";

my @annos;
my $anno;
my $end = "CCA";
my $tRNA_len;
my $id;
my $seq;
my %seqs;
my @species;
my %tRNA_species;
my $tRNA_type;


while (<INPUT>){
	chomp;
	if($_ =~ /^>/){
		@annos = split(/\|/, $_);
		next;
	}
	$seq = $_;
	unless (exists($seqs{$seq})){
		$seqs{$seq} = 1;
		@species = split(/_/, $annos[1]);
		$annos[3] =~ s/\d//g;
		$tRNA_type = $annos[3] . "-" . $annos[4];
		unless (defined $tRNA_species{$tRNA_type}){
			$tRNA_species{$tRNA_type} = 0;
		}
		$tRNA_species{$tRNA_type} += 1;
		$tRNA_len = (length $seq) + 3;
		$anno = ">" . join("_", @species[0,1]) . "_mt_tRNA-" . $tRNA_type . "-" . $tRNA_species{$tRNA_type} . "-1 (mitotRNAdb ChrM) " . $annos[3] . " (" . $annos[4] . ") " . $tRNA_len . " bp mature sequence";
		print OUTPUT_2 "$anno\n$seq$end\n";
		if ($annos[3] eq "His"){
			$seq = substr($seq, 1);
		}
		$tRNA_len = length $seq;
		$anno = ">" . join("_", @species[0,1]) . "_mt_tRNA-" . $tRNA_type . "-" . $tRNA_species{$tRNA_type} . "-1 (mitotRNAdb ChrM) " . $annos[3] . " (" . $annos[4] . ") " . $tRNA_len . " bp sequence";
		print OUTPUT_1 "$anno\n$seq\n";
	}

}
