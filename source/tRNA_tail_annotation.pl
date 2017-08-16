#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';

my $usage = <<"USAGE";
Description:	This script annotate RNA sequence to tRNA 5' end, tRNA 3' end (without -CCA) and tRNA 3' end (with -CCA). 
		refer_file and in_file files are .fa format.
Usage: tRNA_tail_annotation.pl refer_file in_file
USAGE

my $refer_file = shift or die $usage;
my $in_file = shift or die $usage;

my @input;
my $input_name;
my $input_address = abs_path($in_file);;
my $input_suffix;


@input = split(/\//, $input_address);
$input_name = pop (@input);
$input_address = join('/', @input) . '/';
@input = split(/\./, $input_name);
$input_suffix = pop (@input);
$input_name = join('.', @input);

my $output_5_end = $input_address . $input_name . "_output_tRNA_5_tail";
my $output_3_end = $input_address . $input_name . "_output_tRNA_3_tail";
my $output_CCA_end = $input_address . $input_name . "_output_tRNA_CCA_tail";
my $output_not_match = $input_address . $input_name . "_not_match_tRNA_tail.fa";

open REF_FILE, $refer_file
	or die "Can't open '$refer_file': $!\n";
open IN_FILE, $in_file
	or die "Can't open '$in_file': $!\n";
open OUTPUT1, ">$output_5_end"
	or die "Can't open '$output_5_end': $!\n";
open OUTPUT2, ">$output_3_end"
	or die "Can't open '$output_3_end': $!\n";
open OUTPUT3, ">$output_CCA_end"
	or die "Can't open '$output_CCA_end': $!\n";

my %seqs_ref;
my @seq;
my $seq;
my $id;
my $read;
my %seqs_in;
my %annos_in;
my %lens_in;
my $tmp;
	
while (<REF_FILE>){
	chomp;
	if(/^>(.*)/){
		$id = $1;
        @seq = ();
        next;

	}
	push @seq, $_;
    $seq = join('', @seq);
	$seq =~ s/(a|t|c|g|n|)/\U$1/g;
	$seqs_ref{$id} = $seq;
}

while (<IN_FILE>){
	chomp;
	if(/^>(.*)/){

		($id, $read) = (split /\s+/, $1);
		next;
	}
	$seq = $_;
	$lens_in{$id} = length $seq;
	$seqs_in{$id} = $seq;
	$annos_in{$id} = $read;
}

{

foreach my $key1 (sort keys %seqs_in){
	foreach my $key2 (sort keys %seqs_ref){
 		my $subseq_5 = 0;
		my $subseq_3 = 0; 
		my $subseq_CCA = 0;
		if ($seqs_ref{$key2} =~ /(^[ATCGN]{$lens_in{$key1}})/){
			$subseq_5 = $1;
		}
		my $temp_length = $lens_in{$key1}-3;
		if ($seqs_ref{$key2} =~ /([ATCGN]{$temp_length}$)/){
			$subseq_CCA = $1 . 'CCA';
		}
		if ($seqs_ref{$key2} =~ /([ATCGN]{$lens_in{$key1}}$)/){
			$subseq_3 = $1;
		}
		if ($seqs_in{$key1} eq $subseq_5){
            $key2 =~ /\)\s+(.{3,6}?)\s+\((...)\)\s+\d+/;
            $tmp = 'tRNA-' . $1 . '-' . $2 . '_5_end';
			print OUTPUT1 "$key1\t$seqs_in{$key1}\t$lens_in{$key1}\t$tmp\n";
			last;
		}
		elsif ($seqs_in{$key1} eq $subseq_3){
            $key2 =~ /\)\s+(.{3,6}?)\s+\((...)\)\s+\d+/;
            $tmp = 'tRNA-' . $1 . '-' . $2 . '_3_end';
			print OUTPUT2 "$key1\t$seqs_in{$key1}\t$lens_in{$key1}\t$tmp\n";
			last;
		}
		elsif ($seqs_in{$key1} eq $subseq_CCA){
            $key2 =~ /\)\s+(.{3,6}?)\s+\((...)\)\s+\d+/;
            $tmp = 'tRNA-' . $1 . '-' . $2 . '_CCA_end';
			print OUTPUT3 "$key1\t$seqs_in{$key1}\t$lens_in{$key1}\t$tmp\n";
			last;
		}
	}
}
}

close OUTPUT1;
close OUTPUT2;
close OUTPUT3;
close REF_FILE;
close IN_FILE;

{
my %hash;
open OUTPUT1, "$output_5_end"
	or die "Can't open '$output_5_end': $!\n";
open OUTPUT2, "$output_3_end"
	or die "Can't open '$output_3_end': $!\n";
open OUTPUT3, "$output_CCA_end"
	or die "Can't open '$output_CCA_end': $!\n";
open OUTPUT4, ">$output_not_match"
	or die "Can't open '$output_not_match': $!\n";


while (<OUTPUT1>){
	chomp;
	$id = (split /\s+/)[0];
	$hash{$id} = 1;
}

while (<OUTPUT2>){
	chomp;
	$id = (split /\s+/)[0];
	$hash{$id} = 1;
}

while (<OUTPUT3>){
	chomp;
	$id = (split /\s+/)[0];
	$hash{$id} = 1;
}

foreach my $key (keys %seqs_in){
	if (exists $hash{$key}){
		delete $seqs_in{$key};
	}
}
foreach my $key (sort keys %seqs_in){
	print OUTPUT4 ">$key\t$annos_in{$key}\n$seqs_in{$key}\n";
}
}

