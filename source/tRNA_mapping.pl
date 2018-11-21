#!/usr/bin/perl

use strict;
use warnings;


my $usage = <<"USAGE";
Description:	Perl script calculates expression level of genomic loci by bowtie 1 output information.

perl tRNA_mapping.pl bowtie_output_file summary_file
 
USAGE

my $in_file = shift or die $usage;

my $summary = shift or die $usage;


open INPUT1, $in_file
	or die "Can't open '$in_file': $!";

my $total_reads;
{
	open INPUT2, $summary
		or die "Can't open '$summary': $!";
	readline INPUT2;
	$total_reads = readline INPUT2;
	$total_reads = (split /\t/, $total_reads)[2];
	close INPUT2;
}

my ($read, $annotation, $start_site, $seq, $repeat_site, @annos, $temp_anno);
my (%ref_lens, %reads);

while (<INPUT1>){
	($read, $annotation, $start_site, $seq, $repeat_site) = (split /\t/, $_)[1, 3, 4, 5, 7];
	@annos = split(/\s+/, $annotation);
	$annos[0] =~ /(.+)-.+$/;
	$temp_anno = $1;
	$ref_lens{$temp_anno} = $annos[-4];
	for (my $i = $start_site; $i < $start_site + (length $seq); $i++ ){
		$reads{$temp_anno}[$i] += $read / ($repeat_site + 1) *1000000 / $total_reads;
	}
}

foreach my $uniq_id (sort keys %ref_lens){
	for (my $i = 0; $i < ($ref_lens{$uniq_id} +3); $i++){
		unless ($reads{$uniq_id}[$i]){
			$reads{$uniq_id}[$i] = 0;
		}
	}
	print "$uniq_id\t";
	print join( ',', @{$reads{$uniq_id}});
	print "\n";
}


