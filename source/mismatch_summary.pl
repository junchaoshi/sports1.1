#!/usr/bin/perl

use strict;
use warnings;

my $in_file = $ARGV[0] or die "Usage: \n\n $0 input_file threshold > output_file\n";
my $threshold = $ARGV[1] or die "Usage: \n\n $0 input_file threshold > output_file\n";

open INPUT, $in_file
	or die "Can't open '$in_file': $!";

my ($id, $read, $chain, $anno, $pos, $seq, $len, $muts_sum, $ref_name, $mut_pos, $ref_base, $mut_base, $uniq_name, $compare_base);
my (@muts, @mut_info, @base_info);
my (%fq, %reads, %seqs, %ref_ids, %out_ids, %out_seqs, %out_mut_pos, %out_ref_base, %out_mut_base, %out_ref_reads, %out_mut_reads, %out_adj_mut_reads, %out_ref_ids);
while (<INPUT>){
	chomp;
	($id, $read, $chain, $anno, $pos, $seq, $len) = (split /\t/)[0, 1, 2, 3, 4, 5, 6];
	$fq{$id} = 0 unless defined $fq{$id};
	$fq{$id} += 1;
	$reads{$id} = $read;
	$seqs{$id} = $seq;
	$ref_name = $chain . $anno . $pos . $len;
	push @{$ref_ids{$ref_name}}, $id;
}
close INPUT;

open INPUT, $in_file;
while (<INPUT>){
	chomp;
	($id, $read, $chain, $anno, $pos, $seq, $len, $muts_sum) = (split /\t/)[0, 1, 2, 3, 4, 5, 6, 8]; 
	if ($muts_sum){
		$ref_name = $chain . $anno . $pos . $len;
		@muts = (split /,/, $muts_sum);
		foreach(@muts){
			@mut_info = (split /:/, $_);
			$mut_pos = $mut_info[0];
			@base_info = (split />/, $mut_info[1]);
			$ref_base = $base_info[0];
			$mut_base = $base_info[1];
			$uniq_name = $id . $seq . $mut_pos . $ref_base . $mut_base;
			$out_ids{$uniq_name} = $id;
			$out_seqs{$uniq_name} = $seq;
			$out_mut_pos{$uniq_name} = $mut_pos + 1;
			$out_ref_base{$uniq_name} = $ref_base;
			$out_mut_base{$uniq_name} = $mut_base;
			push @{$out_ref_ids{$uniq_name}}, @{$ref_ids{$ref_name}};
			$out_mut_reads{$uniq_name} = $read;
			$out_adj_mut_reads{$uniq_name} = 0 unless defined $out_adj_mut_reads{$uniq_name};
			$out_adj_mut_reads{$uniq_name} += $read/$fq{$id};

		}
	}
}
my %uniq_ids;
foreach $uniq_name (sort keys %out_ids){
	$out_ref_reads{$uniq_name} = 0 unless defined $out_ref_reads{$uniq_name};
	if(@{$out_ref_ids{$uniq_name}}){
		@{$out_ref_ids{$uniq_name}} = do{ my %uniq_ids; grep {!$uniq_ids{$_}++ } @{$out_ref_ids{$uniq_name}}};
		foreach(@{$out_ref_ids{$uniq_name}}){
			$compare_base = substr $seqs{$_}, ($out_mut_pos{$uniq_name}-1), 1; 
			if ($compare_base eq $out_ref_base{$uniq_name}){
				$out_ref_reads{$uniq_name} += $reads{$_};
			}
		}
	}
	if(($out_ref_reads{$uniq_name} + $out_mut_reads{$uniq_name}) >= $threshold){
		print "$out_ids{$uniq_name}\t$out_seqs{$uniq_name}\t$out_mut_pos{$uniq_name}\t$out_ref_base{$uniq_name}\t$out_mut_base{$uniq_name}\t$out_ref_reads{$uniq_name}\t$out_mut_reads{$uniq_name}\t$out_adj_mut_reads{$uniq_name}\n";
	}
}


