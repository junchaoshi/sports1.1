#!/usr/bin/perl

use strict;
use warnings;

my $usage = "annotation.pl filename\n";
my $name = shift or die $usage;
my $in_file = $name . '.fa';
my @out_file = ($name . '_output_match_genome', 
		$name . '_output_miRNA_match_genome',
		$name . '_output_miRNA_unmatch_genome',
		$name . '_output_rRNA_match_genome',
		$name . '_output_rRNA_unmatch_genome',
		$name . '_output_tRNA_match_genome',
		$name . '_output_tRNA_mature_match_genome',
		$name . '_output_tRNA_unmatch_genome',
		$name . '_output_tRNA_mature_unmatch_genome',
		$name . '_output_piRNA_match_genome',
		$name . '_output_piRNA_unmatch_genome',
		$name . '_output_ensembl_match_genome',
		$name . '_output_ensembl_unmatch_genome',
		$name . '_output_rfam_match_genome',
		$name . '_output_rfam_unmatch_genome'
		);
my $out_final = $name . '_output.txt';
my $summary  = $name . '_summary.txt';
my $len_dis = $name . '_length_distribution.txt';

my $file_number = @out_file;
my %file_handle;
my $fh;

open FILE, $in_file
	or die "Can't open '$in_file': $!";

for (1 .. $file_number){
	open $file_handle{$_}, $out_file[$_-1] ;
}

open OUTPUT1, ">$out_final"
	or die "Can't open '$out_final': $!";
printf OUTPUT1 "ID\tSequence\tLength\tReads\tMatch_Genome\tAnnotation\n";

open OUTPUT2, ">$summary"
	or die "can't open '$summary': $!";
print OUTPUT2 "Class\tSub_Class\tReads\n";

open OUTPUT3, ">$len_dis"
	or die "can't open '$len_dis': $!";
printf OUTPUT3 "Class\tLength\tReads\n";

my %reads;
my %match_genome;
my %annos;
my %seqs;
my %lens;
my $id;
my %unannos_match;
my %unannos_unmatch;

######summarize and annotate clean reads######
{
	my %distr;
	my $sums  = 0;
	my $len;
	my $seq;
	my $read;
	while (<FILE>){
		chomp;
		if(/^>(.*)/){
			($id, $read) = (split /\t/, $1);
			next;
		}
		$seq = $_;
		$len = length $seq;
		$lens{$id} = $len;
		$seqs{$id} = $seq;
		$reads{$id} = $read;
		$sums += $read;
		$distr{$len} += $read;
		$unannos_unmatch{$id} = 1;
	}
	print OUTPUT2 "Clean_Reads\t-\t$sums\n";
	foreach $len (sort keys %distr){
		print OUTPUT3 "Clean_Reads\t$len\t$distr{$len}\n";
	}
}

######summarize and annotate match genome reads######
{
	my %distr;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{1};
	while (<$fh>){
		chomp;
		($id, $anno, $seq) = (split /\t/)[0, 3, 5];
		$match_genome{$id} = 'Yes';
		$len = length $seq;
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
	}
	print OUTPUT2 "Match_Genome\t-\t$sums\n";
	foreach $len (sort keys %distr){
	 	print OUTPUT3 "Match_Genome\t$len\t$distr{$len}\n";
	 }
}

######summarize and annotate miRNA reads: match-genome######
if (-e $out_file[1] && !-z $out_file[1]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{2};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[1];
	while (<$fh>){
		chomp;
		($id, $anno, $seq) = (split /\t/)[0, 3, 5];
		$anno = (split /\s+/, $anno)[0];
		$sum{$anno} += ($reads{$id} / $repeat_num{$id});
		$len = length $seq;
		if ($annos{$id}){
			$anno =~ s/\?/\\\?/g;
			unless ($annos{$id} =~ /$anno/){
				$anno =~ s/\\\?/\?/g;
				$annos{$id} = $annos{$id} . ';' . $anno;
			}
		}else{
			$anno =~ s/\\\?/\?/g;
			$annos{$id} = $anno;
		}
		$sums += ($reads{$id} / $repeat_num{$id});
		$distr{$len} += ($reads{$id} / $repeat_num{$id});
		if ($unannos_unmatch{$id}){
			delete $unannos_unmatch{$id};
		}
	}
	if ($sums > 0){
		printf OUTPUT2 "miRBase-miRNA_Match_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "-\t$anno\t%.2f\n", $sum{$anno};
		}
		foreach $len (sort keys %distr){
			printf OUTPUT3 "miRBase-miRNA_Match_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
}

######summarize and annotate miRNA reads: unmatch-genome######
if (-e $out_file[2] && !-z $out_file[2]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{3};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[2];
	while (<$fh>){
		chomp;
		($id, $anno, $seq) = (split /\t/)[0, 3, 5];
		$anno = (split /\s+/, $anno)[0];
		$sum{$anno} += ($reads{$id} / $repeat_num{$id});
		$len = length $seq;
		if ($annos{$id}){
			$anno =~ s/\?/\\\?/g;
			unless ($annos{$id} =~ /$anno/){
				$anno =~ s/\\\?/\?/g;
				$annos{$id} = $annos{$id} . ';' . $anno;
			}
		}else{
			$anno =~ s/\\\?/\?/g;
			$annos{$id} = $anno;
		}
		$sums += ($reads{$id} / $repeat_num{$id});
		$distr{$len} += ($reads{$id} / $repeat_num{$id});
		if ($unannos_unmatch{$id}){
			delete $unannos_unmatch{$id};
		}
	}
	if ($sums > 0){
		printf OUTPUT2 "miRBase-miRNA_Unmatch_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "-\t%.2f\n", $sum{$anno};
		}
		foreach $len (sort keys %distr){
			printf OUTPUT3 "miRBase-miRNA_Unmatch_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
}

######summarize and annotate rRNA reads: match-genome######
if (-e $out_file[3] && !-z $out_file[3]){

	my %distr_a;
	my %distr_b;
	my %sum;
	my %repeat_num;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{4};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[3];
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*\s([0-9]+\.[0-9]+S)\s/){
			$id = $1;
			$anno = $2 . '-rRNA';
			if (/([ATCG]+?)\s[I]+/){
				$seq = $1;
			}
		}
		elsif (/^(t[0-9]+?)\s.*\s([0-9]+S)\s/){
			$id = $1;
			$anno = $2 . '-rRNA';
			if (/([ATCG]+?)\s[I]+/){
				$seq = $1;
			}		
		}
		elsif (/^(t[0-9]+?)\s/){
			$id = $1;
			$anno = "other-rRNA";
			if (/([ATCG]+?)\s[I]+/){
				$seq = $1;
			}
		}
		$sum{$anno} += ($reads{$id} / $repeat_num{$id});
		$len = length $seq;
		if ($annos{$id}){
			$anno =~ s/\?/\\\?/g;
			unless ($annos{$id} =~ /$anno/){
				$anno =~ s/\\\?/\?/g;
				$annos{$id} = $annos{$id} . ';' . $anno;
			}
		}else{
			$anno =~ s/\\\?/\?/g;
			$annos{$id} = $anno;
		}
		$sums += ($reads{$id} / $repeat_num{$id});
		$distr_a{$len} += ($reads{$id} / $repeat_num{$id});
		$distr_b{$anno}{$len} = 0 unless defined $distr_b{$anno}{$len};
		$distr_b{$anno}{$len} += ($reads{$id} / $repeat_num{$id});
		if ($unannos_unmatch{$id}){
			delete $unannos_unmatch{$id};
		}

	}
	if ($sums > 0){
		printf OUTPUT2 "rRNAdb-rRNA_Match_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "-\t${anno}_Match_Genome\t%.2f\n", $sum{$anno};
		}
		foreach $len(sort keys %distr_a){
			printf OUTPUT3 "rRNAdb-rRNA_Match_Genome\t$len\t%.2f\n", $distr_a{$len};
		}
		foreach $anno(sort keys %distr_b){
			foreach $len(sort keys %{$distr_b{$anno}}){
				printf OUTPUT3 "rRNAdb-${anno}_Match_Genome\t$len\t%0.2f\n", $distr_b{$anno}{$len};
			}
		}
	}
}

######summarize and annotate rRNA reads: unmatch-genome######
if (-e $out_file[4] && !-z $out_file[4]){

	my %distr_a;
	my %distr_b;
	my %sum;
	my %repeat_num;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{5};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[4];
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*\s([0-9]+\.[0-9]+S)\s/){
			$id = $1;
			$anno = $2 . '-rRNA';
			if (/([ATCG]+?)\s[I]+/){
				$seq = $1;
			}
		}
		elsif (/^(t[0-9]+?)\s.*\s([0-9]+S)\s/){
			$id = $1;
			$anno = $2 . '-rRNA';
			if (/([ATCG]+?)\s[I]+/){
				$seq = $1;
			}		
		}
		elsif (/^(t[0-9]+?)\s/){
			$id = $1;
			$anno = "other-rRNA";
			if (/([ATCG]+?)\s[I]+/){
				$seq = $1;
			}
		}
		$sum{$anno} += ($reads{$id} / $repeat_num{$id});
		$len = length $seq;
		if ($annos{$id}){
			$anno =~ s/\?/\\\?/g;
			unless ($annos{$id} =~ /$anno/){
				$anno =~ s/\\\?/\?/g;
				$annos{$id} = $annos{$id} . ';' . $anno;
			}
		}else{
			$anno =~ s/\\\?/\?/g;
			$annos{$id} = $anno;
		}
		$sums += ($reads{$id} / $repeat_num{$id});
		$distr_a{$len} += ($reads{$id} / $repeat_num{$id});
		$distr_b{$anno}{$len} = 0 unless defined $distr_b{$anno}{$len};
		$distr_b{$anno}{$len} += ($reads{$id} / $repeat_num{$id});
		if ($unannos_unmatch{$id}){
			delete $unannos_unmatch{$id};
		}

	}
	if ($sums > 0){
		printf OUTPUT2 "rRNAdb-rRNA_Unmatch_Genome\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "-\t${anno}_Unmatch_Genome\t%.2f\n", $sum{$anno};
		}
		foreach $len(sort keys %distr_a){
			printf OUTPUT3 "rRNAdb-rRNA_Unmatch_Genome\t$len\t%.2f\n", $distr_a{$len};
		}
    		foreach $anno(sort keys %distr_b){
			foreach $len(sort keys %{$distr_b{$anno}}){
				printf OUTPUT3 "rRNAdb-${anno}_Unmatch_Genome\t$len\t%0.2f\n", $distr_b{$anno}{$len};
			}
		}
	}
}

######summarize and annotate tRNA reads: match-genome######
{
	my %distr;
	my %distr_5_end;
	my %distr_3_end;
	my %distr_CCA_end;
	my %sum;
	my %sum_5_end;
	my %sum_3_end;
	my %sum_CCA_end;
	my %repeat_num;
	my $sums  = 0;
	my $sums_5_end = 0;
	my $sums_3_end = 0;
	my $sums_CCA_end = 0;
	my $len;
	my $anno;
	my $seq;
	my $temp_anno;
	my $tRNA_len;
	my $start_site;
	my @annotation;
	my $item;
	if (-e $out_file[5] && !-z $out_file[5]){
		$fh = $file_handle{6};
		while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[5];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			@annotation = split(/ +/, $anno);
			$tRNA_len = $annotation[-4];
			$len = length $seq;
			if ($start_site == 0){
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1 . '_5_end';
				$anno = $temp_anno;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$sum_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site + 1)){
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1 . '_3_end';
				$anno = $temp_anno;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$sum_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_3_end += ($reads{$id} / $repeat_num{$id});
			}else {
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$anno = $temp_anno;
				$sum{$temp_anno} += ($reads{$id} / $repeat_num{$id});
			}
			if ($anno =~ /Undet/){
				$anno =~ s/\?/\\\?/g;
			}
			if ($annos{$id}){
				unless ($annos{$id} =~ /$anno/){
					$anno =~ s/\\\?/\?/g;
					$annos{$id} = $annos{$id} . ';' . $anno;
				}
			}else{
				$anno =~ s/\\\?/\?/g;
				$annos{$id} = $anno;
			}
			$sums += ($reads{$id} / $repeat_num{$id});
			$distr{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
	}
	if (-e $out_file[6] && !-z $out_file[6]){
		$fh = $file_handle{7};
		while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[6];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			@annotation = split(/ +/, $anno);
			$tRNA_len = $annotation[-4];
			$len = length $seq;
			if ($start_site == 0){
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1 . '_5_end';
				$anno = $temp_anno;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$sum_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site + 1)){
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1 . '_3_end';
				$anno = $temp_anno;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$sum_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_3_end += ($reads{$id} / $repeat_num{$id} );
			}elsif ($tRNA_len == ($len + $start_site - 2)){
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1 . '_CCA_end';
				$anno = $temp_anno;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$sum_CCA_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_CCA_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_CCA_end += ($reads{$id} / $repeat_num{$id});
			}else{
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$anno = $temp_anno;
				$sum{$temp_anno} += ($reads{$id} / $repeat_num{$id});
			}
			if ($anno =~ /Undet/){
				$anno =~ s/\?/\\\?/g;
			}
			if ($annos{$id}){
				unless ($annos{$id} =~ /$anno/){
					$anno =~ s/\\\?/\?/g;
					$annos{$id} = $annos{$id} . ';' . $anno;
				}
			}else{
				$anno =~ s/\\\?/\?/g;
				$annos{$id} = $anno;
			}
			$sums += ($reads{$id} / $repeat_num{$id});
			$distr{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}

	}
	if ($sums > 0){
		printf OUTPUT2 "GtRNAdb-tRNA_Match_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "-\t$anno\t%.2f\n", $sum{$anno};
		}
		foreach $len (sort keys %distr){
			printf OUTPUT3 "GtRNAdb-tRNA_Match_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
	if ($sums_5_end > 0){
		printf OUTPUT2 "GtRNAdb-tRNA_5_end_Match_Genome\t-\t%d\n", $sums_5_end;
		foreach $anno(sort keys %sum_5_end){
			printf OUTPUT2 "-\t$anno\t%.2f\n", $sum_5_end{$anno};
		}
		foreach $len (sort keys %distr_5_end){
			printf OUTPUT3 "GtRNAdb-tRNA_5_end_Match_Genome\t$len\t%.2f\n", $distr_5_end{$len};
		}
	}
	if ($sums_3_end > 0){
		printf OUTPUT2 "GtRNAdb-tRNA_3_end_Match_Genome\t-\t%d\n", $sums_3_end;
		foreach $anno(sort keys %sum_3_end){
			printf OUTPUT2 "-\t$anno\t%.2f\n", $sum_3_end{$anno};
		}
		foreach $len (sort keys %distr_3_end){
			printf OUTPUT3 "GtRNAdb-tRNA_3_end_Match_Genome\t$len\t%.2f\n", $distr_3_end{$len};
		}
	}
	if ($sums_CCA_end > 0){
		printf OUTPUT2 "GtRNAdb-tRNA_CCA_end_Match_Genome\t-\t%d\n", $sums_CCA_end;
		foreach $anno(sort keys %sum_CCA_end){
			printf OUTPUT2 "-\t$anno\t%.2f\n", $sum_CCA_end{$anno};
		}
		foreach $len (sort keys %distr_CCA_end){
			printf OUTPUT3 "GtRNAdb-tRNA_CCA_end_Match_Genome\t$len\t%.2f\n", $distr_CCA_end{$len};
		}
	}
}

######summarize and annotate tRNA reads: unmatch-genome######
{
	my %distr;
	my %distr_5_end;
	my %distr_3_end;
	my %distr_CCA_end;
	my %sum;
	my %sum_5_end;
	my %sum_3_end;
	my %sum_CCA_end;
	my %repeat_num;
	my $sums  = 0;
	my $sums_5_end = 0;
	my $sums_3_end = 0;
	my $sums_CCA_end = 0;
	my $len;
	my $anno;
	my $seq;
	my $temp_anno;
	my $tRNA_len;
	my $start_site;
	my @annotation;
	my $item;
	if (-e $out_file[7] && !-z $out_file[7]){
		$fh = $file_handle{8};
		while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[7];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			@annotation = split(/ +/, $anno);
			$tRNA_len = $annotation[-4];
			$len = length $seq;
			if ($start_site == 0){
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1 . '_5_end';
				$anno = $temp_anno;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$sum_5_end{$anno} += ($reads{$id} / ($repeat_num{$id} + 1));
				$sum{$temp_anno} += ($reads{$id} / ($repeat_num{$id} + 1));
				$distr_5_end{$len} += ($reads{$id} / ($repeat_num{$id} + 1));
				$sums_5_end += ($reads{$id} / ($repeat_num{$id} + 1));
			}elsif ($tRNA_len == ($len + $start_site + 1)){
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1 . '_3_end';
				$anno = $temp_anno;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$sum_3_end{$anno} += ($reads{$id} / ($repeat_num{$id} + 1));
				$sum{$temp_anno} += ($reads{$id} / ($repeat_num{$id} + 1));
				$distr_3_end{$len} += ($reads{$id} / ($repeat_num{$id} + 1));
				$sums_3_end += ($reads{$id} / ($repeat_num{$id} + 1));
			}else {
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$anno = $temp_anno;
				$sum{$temp_anno} += ($reads{$id} / ($repeat_num{$id} + 1));
			}
			if ($anno =~ /Undet/){
				$anno =~ s/\?/\\\?/g;
			}
			if ($annos{$id}){
				unless ($annos{$id} =~ /$anno/){
					$anno =~ s/\\\?/\?/g;
					$annos{$id} = $annos{$id} . ';' . $anno;
				}
			}else{
				$anno =~ s/\\\?/\?/g;
				$annos{$id} = $anno;
			}
			$sums += ($reads{$id} / ($repeat_num{$id} + 1));
			$distr{$len} += ($reads{$id} / ($repeat_num{$id} + 1));
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
	}
	if (-e $out_file[8] && !-z $out_file[8]){
		$fh = $file_handle{9};
		while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[8];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			@annotation = split(/ +/, $anno);
			$tRNA_len = $annotation[-4];
			$len = length $seq;
			if ($start_site == 0){
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1 . '_5_end';
				$anno = $temp_anno;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$sum_5_end{$anno} += ($reads{$id} / ($repeat_num{$id} + 1));
				$sum{$temp_anno} += ($reads{$id} / ($repeat_num{$id} + 1));
				$distr_5_end{$len} += ($reads{$id} / ($repeat_num{$id} + 1));
				$sums_5_end += ($reads{$id} / ($repeat_num{$id} + 1));
			}elsif ($tRNA_len == ($len + $start_site + 1)){
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1 . '_3_end';
				$anno = $temp_anno;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$sum_3_end{$anno} += ($reads{$id} / ($repeat_num{$id} + 1));
				$sum{$temp_anno} += ($reads{$id} / ($repeat_num{$id} + 1));
				$distr_3_end{$len} += ($reads{$id} / ($repeat_num{$id} + 1));
				$sums_3_end += ($reads{$id} / ($repeat_num{$id} + 1));
			}elsif ($tRNA_len == ($len + $start_site - 2)){
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1 . '_CCA_end';
				$anno = $temp_anno;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$sum_CCA_end{$anno} += ($reads{$id} / ($repeat_num{$id} + 1));
				$sum{$temp_anno} += ($reads{$id} / ($repeat_num{$id} + 1));
				$distr_CCA_end{$len} += ($reads{$id} / ($repeat_num{$id} + 1));
				$sums_CCA_end += ($reads{$id} / ($repeat_num{$id} + 1));
			}else{
				$annotation[-5] =~ /\((.+)\)/;
				$temp_anno = 'tRNA-' . $annotation[-6] . '-' . $1;
				$anno = $temp_anno;
				$sum{$temp_anno} += ($reads{$id} / ($repeat_num{$id} + 1));
			}
			if ($anno =~ /Undet/){
				$anno =~ s/\?/\\\?/g;
			}
			if ($annos{$id}){
				unless ($annos{$id} =~ /$anno/){
					$anno =~ s/\\\?/\?/g;
					$annos{$id} = $annos{$id} . ';' . $anno;
				}
			}else{
				$anno =~ s/\\\?/\?/g;
				$annos{$id} = $anno;
			}
			$sums += ($reads{$id} / ($repeat_num{$id} + 1));
			$distr{$len} += ($reads{$id} / ($repeat_num{$id} + 1));
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}

	}
	if ($sums > 0){
		printf OUTPUT2 "GtRNAdb-tRNA_Unmatch_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "-\t$anno\t%.2f\n", $sum{$anno};
		}
		foreach $len (sort keys %distr){
			printf OUTPUT3 "GtRNAdb-tRNA_Unmatch_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
	if ($sums_5_end > 0){
		printf OUTPUT2 "GtRNAdb-tRNA_5_end_Unmatch_Genome\t-\t%d\n", $sums_5_end;
		foreach $anno(sort keys %sum_5_end){
			printf OUTPUT2 "-\t$anno\t%.2f\n", $sum_5_end{$anno};
		}
		foreach $len (sort keys %distr_5_end){
			printf OUTPUT3 "GtRNAdb-tRNA_5_end_Unmatch_Genome\t$len\t%.2f\n", $distr_5_end{$len};
		}
	}
	if ($sums_3_end > 0){
		printf OUTPUT2 "GtRNAdb-tRNA_3_end_Unmatch_Genome\t-\t%d\n", $sums_3_end;
		foreach $anno(sort keys %sum_3_end){
			printf OUTPUT2 "-\t$anno\t%.2f\n", $sum_3_end{$anno};
		}
		foreach $len (sort keys %distr_3_end){
			printf OUTPUT3 "GtRNAdb-tRNA_3_end_Unmatch_Genome\t$len\t%.2f\n", $distr_3_end{$len};
		}
	}
	if ($sums_CCA_end > 0){
		printf OUTPUT2 "GtRNAdb-tRNA_CCA_end_Unmatch_Genome\t-\t%d\n", $sums_CCA_end;
		foreach $anno(sort keys %sum_CCA_end){
			printf OUTPUT2 "-\t$anno\t%.2f\n", $sum_CCA_end{$anno};
		}
		foreach $len (sort keys %distr_CCA_end){
			printf OUTPUT3 "GtRNAdb-tRNA_CCA_end_Unmatch_Genome\t$len\t%.2f\n", $distr_CCA_end{$len};
		}
	}
}

######summarize and annotate using piRNA database: match-genome######
if (-e $out_file[9] && !-z $out_file[9]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{10};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[9];
	while (<$fh>){
		chomp;
		($id, $seq) = (split /\s+/)[0, 5];
		$anno = 'piRNA';
		$sum{$anno} += ($reads{$id} / $repeat_num{$id});
		$len = length $seq;
		$annos{$id} = $anno;
		$sums += ($reads{$id} / $repeat_num{$id});
		$distr{$len} += ($reads{$id} / $repeat_num{$id});
		if ($unannos_unmatch{$id}){
			delete $unannos_unmatch{$id};
		}
	}
	if ($sums > 0){
		printf OUTPUT2 "piRNAdb-piRNA_Match_Genome\t-\t%d\n", $sums;
		foreach $len (sort keys %distr){
			printf OUTPUT3 "piRNAdb-piRNA_Match_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
}

######summarize and annotate using piRNA database: unmatch-genome######
if (-e $out_file[10] && !-z $out_file[10]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{11};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[10];
	while (<$fh>){
		chomp;
		($id, $seq) = (split /\s+/)[0, 5];
		$anno = 'piRNA';
		$sum{$anno} += ($reads{$id} / $repeat_num{$id});
		$len = length $seq;
		$annos{$id} = $anno;
		$sums += ($reads{$id} / $repeat_num{$id});
		$distr{$len} += ($reads{$id} / $repeat_num{$id});
		if ($unannos_unmatch{$id}){
			delete $unannos_unmatch{$id};
		}
	}
	if ($sums > 0){
		printf OUTPUT2 "piRNAdb-piRNA_Unmatch_Genome\t-\t%d\n", $sums;
		foreach $len (sort keys %distr){
			printf OUTPUT3 "piRNAdb-piRNA_Unmatch_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
}

######summarize and annotate using ensembl database: match-genome######
if (-e $out_file[11] && !-z $out_file[11]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums  = 0;
	my $len;
	my $anno;
	$fh = $file_handle{12};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[11];
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*gene_biotype:(\S+?)\s+/){
			$id = $1;
			$anno = $2;
			if ($annos{$id}){
				$anno =~ s/\?/\\\?/g;
				unless ($annos{$id} =~ /$anno/){
					$anno =~ s/\\\?/\?/g;
					$annos{$id} = $annos{$id} . ';' . $anno;
				}
			}else{
				$anno =~ s/\\\?/\?/g;
				$annos{$id} = $anno;
			}
			$sum{$anno} += ($reads{$id} / $repeat_num{$id});
			$len = $lens{$id};
			$distr{$anno}{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
	}
	if ($sums > 0){
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "ensembl-${anno}_Match_Genome\t-\t%d\n", $sum{$anno};
		}
		foreach	my $key1(sort keys %distr){
			foreach my $key2 (sort keys %{$distr{$key1}}){
				printf OUTPUT3 "ensembl-${key1}_Match_Genome\t$key2\t%.2f\n", $distr{$key1}{$key2};
			}
		}
	}
}

######summarize and annotate using ensembl database: unmatch-genome######
if (-e $out_file[12] && !-z $out_file[12]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums  = 0;
	my $len;
	my $anno;
	$fh = $file_handle{13};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[12];
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*gene_biotype:(\S+?)\s+/){
			$id = $1;
			$anno = $2;
			if ($annos{$id}){
				$anno =~ s/\?/\\\?/g;
				unless ($annos{$id} =~ /$anno/){
					$anno =~ s/\\\?/\?/g;
					$annos{$id} = $annos{$id} . ';' . $anno;
				}
			}else{
				$anno =~ s/\\\?/\?/g;
				$annos{$id} = $anno;
			}
			$sum{$anno} += ($reads{$id} / $repeat_num{$id});
			$len = $lens{$id};
			$distr{$anno}{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
	}


	if ($sums > 0){
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "ensembl-${anno}_Unmatch_Genome\t-\t%d\n", $sum{$anno};
		}
		foreach	my $key1(sort keys %distr){
			foreach my $key2 (sort keys %{$distr{$key1}}){
				printf OUTPUT3 "ensembl-${key1}_Unmatch_Genome\t$key2\t%.2f\n", $distr{$key1}{$key2};
			}
		}
	}
}

######summarize and annotate using rfam database: match-genome######
if (-e $out_file[13] && !-z $out_file[13]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums  = 0;
	my $len;
	my $anno;
	$fh = $file_handle{14};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[13];
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*(mRNA|cds|protein)/){
			$id = $1;
			$anno = 'mRNA';
		}
		elsif (/^(t[0-9]+?)\s.*(miRNA|microRNA)/){
			$id = $1;
			$anno = 'miRNA';
		}
		elsif (/^(t[0-9]+?)\s.*(tRNA-like)/){
			$id = $1;
			$anno = 'tRNA-like';			
		}
		elsif (/^(t[0-9]+?)\s.*(tRNA|transfer RNA)/){
			$id = $1;
			$anno = 'tRNA';		
		}
		elsif (/^(t[0-9]+?)\s.*(rRNA|ribosomal RNA|ribosomal DNA)/){
			$id = $1;
			$anno = 'rRNA';			
		}
		elsif (/^(t[0-9]+?)\s.*(piRNA)/){
			$id = $1;
			$anno = 'piRNA';		
		}
		elsif (/^(t[0-9]+?)\s.*(snRNA|small nuclear RNA|small-nuclear RNA)/){
			$id = $1;
			$anno = 'snRNA';	
		}
		elsif (/^(t[0-9]+?)\s.*(snoRNA)/){
			$id = $1;
			$anno = 'snoRNA';		
		}
		elsif (/^(t[0-9]+?)\s.*(scRNA)/){
			$id = $1;
			$anno = 'scRNA';		
		}
		elsif (/^(t[0-9]+?)\s.*(lncRNA|lincRNA)/){
			$id = $1;
			$anno = 'lncRNA';
		}
		elsif (/^(t[0-9]+?)\s.*(antisense)/){
			$id = $1;
			$anno = 'antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(vault RNA)/){
			$id = $1;
			$anno = 'vault_RNA';
		}
		elsif (/^(t[0-9]+?)\s.*(telomerase RNA)/){
			$id = $1;
			$anno = 'telomerase_RNA';
		}
		elsif (/^(t[0-9]+?)\s.*(noncoding)/){
			$id = $1;
			$anno = 'noncoding_RNA';
		}
		elsif (/^(t[0-9]+?)\s.*(RNA)/){
			$id = $1;
			$anno = 'other_RNA';			
		}
		elsif (/^(t[0-9]+?)\s.*(pseudogene)/){
			$id = $1;
			$anno = 'pseudogene';		
		}
		elsif (/^(t[0-9]+?)\s.*(gene)/){
			$id = $1;
			$anno = 'gene_region';		
		}
		elsif (/^(t[0-9]+?)\s.*(mitochondrial|mitochondrion)/){
			$id = $1;
			$anno = 'mt_DNA';			
		}
		elsif (/^(t[0-9]+?)\s.*(genome|chromosome|BAC|genomic|X-inactivation center)/){
			$id = $1;
			$anno = 'DNA';		
		}
		elsif (/^(t[0-9]+?)\s/){
			$id = $1;
			$anno = 'other';	
		}
		if ($annos{$id}){
			unless ($annos{$id} =~ /$anno/){
				$annos{$id} = $annos{$id} . ';' . $anno;
			}
		}else{
			$annos{$id} = $anno;
		}
		$sum{$anno} += ($reads{$id} / $repeat_num{$id});
		$len = $lens{$id};
		$distr{$anno}{$len} += ($reads{$id} / $repeat_num{$id});
		if ($unannos_unmatch{$id}){
			delete $unannos_unmatch{$id};
		}
	}
	foreach $anno(sort keys %sum){
		printf OUTPUT2 "Rfam-${anno}_Match_Genome\t-\t%d\n", $sum{$anno};
	}
	foreach	my $key1(sort keys %distr){
		foreach my $key2 (sort keys %{$distr{$key1}}){
			printf OUTPUT3 "Rfam-${key1}_Match_Genome\t$key2\t%.2f\n", $distr{$key1}{$key2};
		}
	}	
}

######summarize and annotate using rfam database: unmatch-genome######
if (-e $out_file[14] && !-z $out_file[14]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums  = 0;
	my $len;
	my $anno;
	$fh = $file_handle{15};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[14];
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*(mRNA|cds|protein)/){
			$id = $1;
			$anno = 'mRNA';
		}
		elsif (/^(t[0-9]+?)\s.*(miRNA|microRNA)/){
			$id = $1;
			$anno = 'miRNA';
		}
		elsif (/^(t[0-9]+?)\s.*(tRNA-like)/){
			$id = $1;
			$anno = 'tRNA-like';		
		}
		elsif (/^(t[0-9]+?)\s.*(tRNA|transfer RNA)/){
			$id = $1;
			$anno = 'tRNA';			
		}
		elsif (/^(t[0-9]+?)\s.*(rRNA|ribosomal RNA|ribosomal DNA)/){
			$id = $1;
			$anno = 'rRNA';			
		}
		elsif (/^(t[0-9]+?)\s.*(piRNA)/){
			$id = $1;
			$anno = 'piRNA';		
		}
		elsif (/^(t[0-9]+?)\s.*(snRNA|small nuclear RNA|small-nuclear RNA)/){
			$id = $1;
			$anno = 'snRNA';	
		}
		elsif (/^(t[0-9]+?)\s.*(snoRNA)/){
			$id = $1;
			$anno = 'snoRNA';		
		}
		elsif (/^(t[0-9]+?)\s.*(scRNA)/){
			$id = $1;
			$anno = 'scRNA';		
		}
		elsif (/^(t[0-9]+?)\s.*(lncRNA|lincRNA)/){
			$id = $1;
			$anno = 'lncRNA';
		}
		elsif (/^(t[0-9]+?)\s.*(antisense)/){
			$id = $1;
			$anno = 'antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(vault RNA)/){
			$id = $1;
			$anno = 'vault_RNA';
		}
		elsif (/^(t[0-9]+?)\s.*(telomerase RNA)/){
			$id = $1;
			$anno = 'telomerase_RNA';
		}
		elsif (/^(t[0-9]+?)\s.*(noncoding)/){
			$id = $1;
			$anno = 'noncoding_RNA';
		}
		elsif (/^(t[0-9]+?)\s.*(RNA)/){
			$id = $1;
			$anno = 'other_RNA';			
		}
		elsif (/^(t[0-9]+?)\s.*(pseudogene)/){
			$id = $1;
			$anno = 'pseudogene';	
		}
		elsif (/^(t[0-9]+?)\s.*(gene)/){
			$id = $1;
			$anno = 'gene_region';		
		}
		elsif (/^(t[0-9]+?)\s.*(mitochondrial|mitochondrion)/){
			$id = $1;
			$anno = 'mt_DNA';		
		}
		elsif (/^(t[0-9]+?)\s.*(genome|chromosome|BAC|genomic|X-inactivation center)/){
			$id = $1;
			$anno = 'DNA';	
		}
		elsif (/^(t[0-9]+?)\s/){
			$id = $1;
			$anno = 'other';	
		}
		if ($annos{$id}){
			unless ($annos{$id} =~ /$anno/){
				$annos{$id} = $annos{$id} . ';' . $anno;
			}
		}else{
			$annos{$id} = $anno;
		}
		$sum{$anno} += ($reads{$id} / $repeat_num{$id});
		$len = $lens{$id};
		$distr{$anno}{$len} += ($reads{$id} / $repeat_num{$id});
		if ($unannos_unmatch{$id}){
			delete $unannos_unmatch{$id};
		}
	}
	foreach $anno(sort keys %sum){
		printf OUTPUT2 "Rfam-${anno}_Unmatch_Genome\t-\t%d\n", $sum{$anno};
	}
	foreach	my $key1(sort keys %distr){
		foreach my $key2 (sort keys %{$distr{$key1}}){
			printf OUTPUT3 "Rfam-${key1}_Unmatch_Genome\t$key2\t%.2f\n", $distr{$key1}{$key2};
		}
	}	
}

######summarize unannotated match genome reads######
{
	my %distr;
	my $sums  = 0;
	my $len;
	my $seq;
	my $read;
	foreach $id(sort keys %unannos_unmatch){
		if($match_genome{$id}){
			$unannos_match{$id} = 1;
			delete $unannos_unmatch{$id};
		}

	}
	foreach $id(sort keys %unannos_match){
		$len = length $seqs{$id};
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
	}
	print OUTPUT2 "Unannotated_Match_Genome\t-\t$sums\n";
	foreach $len (sort keys %distr){
		print OUTPUT3 "Unannotated_Match_Genome\t$len\t$distr{$len}\n";
	}
}

######summarize unannotated unmatch genome reads######
{
	my %distr;
	my $sums  = 0;
	my $len;
	my $seq;
	my $read;
	foreach $id(sort keys %unannos_unmatch){
		$len = length $seqs{$id};
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
	}
	print OUTPUT2 "Unannotated_Unmatch_Genome\t-\t$sums\n";
	foreach $len (sort keys %distr){
		print OUTPUT3 "Unannotated_Unmatch_Genome\t$len\t$distr{$len}\n";
	}
}


######final summarize and annotate######

foreach $id(sort keys %seqs){
    if($match_genome{$id} && $annos{$id}){
        print OUTPUT1 "$id\t$seqs{$id}\t$lens{$id}\t$reads{$id}\t$match_genome{$id}\t$annos{$id}\n";
    }elsif($match_genome{$id}){
        print OUTPUT1 "$id\t$seqs{$id}\t$lens{$id}\t$reads{$id}\t$match_genome{$id}\tNO_Annotation\n";
    }elsif($annos{$id}){
        print OUTPUT1 "$id\t$seqs{$id}\t$lens{$id}\t$reads{$id}\tNo\t$annos{$id}\n";
    }else{
        print OUTPUT1 "$id\t$seqs{$id}\t$lens{$id}\t$reads{$id}\tNo\tNO_Annotation\n";
    }
}

