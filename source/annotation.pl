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
		$name . '_output_tRNA_pre_match_genome',
		$name . '_output_tRNA_mature_match_genome',
		$name . '_output_tRNA_pre_unmatch_genome',
		$name . '_output_tRNA_mature_unmatch_genome',
		$name . '_output_mt_tRNA_pre_match_genome',
		$name . '_output_mt_tRNA_mature_match_genome',
		$name . '_output_mt_tRNA_pre_unmatch_genome',
		$name . '_output_mt_tRNA_mature_unmatch_genome',
		$name . '_output_piRNA_match_genome',
		$name . '_output_piRNA_unmatch_genome',
		$name . '_output_ensembl_match_genome',
		$name . '_output_ensembl_unmatch_genome',
		$name . '_output_rfam_match_genome',
		$name . '_output_rfam_unmatch_genome',
		$name . '_output_miRNA-antisense_match_genome',
		$name . '_output_miRNA-antisense_unmatch_genome',
		$name . '_output_rRNA-antisense_match_genome',
		$name . '_output_rRNA-antisense_unmatch_genome',
		$name . '_output_tRNA-antisense_pre_match_genome',
		$name . '_output_tRNA-antisense_mature_match_genome',
		$name . '_output_tRNA-antisense_pre_unmatch_genome',
		$name . '_output_tRNA-antisense_mature_unmatch_genome',
		$name . '_output_mt_tRNA-antisense_pre_match_genome',
		$name . '_output_mt_tRNA-antisense_mature_match_genome',
		$name . '_output_mt_tRNA-antisense_pre_unmatch_genome',
		$name . '_output_mt_tRNA-antisense_mature_unmatch_genome',
		$name . '_output_piRNA-antisense_match_genome',
		$name . '_output_piRNA-antisense_unmatch_genome',
		$name . '_output_ensembl-antisense_match_genome',
		$name . '_output_ensembl-antisense_unmatch_genome',
		$name . '_output_rfam-antisense_match_genome',
		$name . '_output_rfam-antisense_unmatch_genome'
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
my $i = 1;

######summarize and annotate clean reads######
{
	my %distr;
	my $sums = 0;
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
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{$i};
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

######summarize and annotate miRNA seqs: match-genome######
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
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
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "miRBase-miRNA_Match_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "miRBase-miRNA_Match_Genome\t$anno\t%.2f\n", $sum{$anno};
		}
		foreach $len (sort keys %distr){
			printf OUTPUT3 "miRBase-miRNA_Match_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
}

######summarize and annotate miRNA seqs: unmatch-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
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
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "miRBase-miRNA_Unmatch_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "miRBase-miRNA_Unmatch_Genome\t$anno\t%.2f\n", $sum{$anno};
		}
		foreach $len (sort keys %distr){
			printf OUTPUT3 "miRBase-miRNA_Unmatch_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
}

######summarize and annotate rRNA seqs: match-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr_a;
	my %distr_b;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		if(!/^(t[0-9]+?)\s.*\s(RNY[0-9])\s/){
			if (/^(t[0-9]+?)\s.*\s([0-9]+\.[0-9]+S)\s/){
				$id = $1;
				$anno = $2 . '-rRNA';
				if (/([ATCGN]+?)\s[I]+/){
					$seq = $1;
				}
			}
			elsif (/^(t[0-9]+?)\s.*\s([0-9]+S)\s/){
				$id = $1;
				$anno = $2 . '-rRNA';
				if (/([ATCGN]+?)\s[I]+/){
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
	}
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "rRNAdb-rRNA_Match_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "rRNAdb-rRNA_Match_Genome\t${anno}\t%.2f\n", $sum{$anno};
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

if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr_a;
	my %distr_b;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		if(/^(t[0-9]+?)\s.*\s(RNY[0-9])\s/){
			$id = $1;
			$anno = $2 . '-YRNA';
			if (/([ATCGN]+?)\s[I]+/){
				$seq = $1;
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
	}
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "YRNAdb-YRNA_Match_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "YRNAdb-YRNA_Match_Genome\t${anno}\t%.2f\n", $sum{$anno};
		}
		foreach $len(sort keys %distr_a){
			printf OUTPUT3 "YRNAdb-YRNA_Match_Genome\t$len\t%.2f\n", $distr_a{$len};
		}
		foreach $anno(sort keys %distr_b){
			foreach $len(sort keys %{$distr_b{$anno}}){
				printf OUTPUT3 "YRNAdb-${anno}_Match_Genome\t$len\t%0.2f\n", $distr_b{$anno}{$len};
			}
		}
	}
}

######summarize and annotate rRNA seqs: unmatch-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr_a;
	my %distr_b;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		if(!/^(t[0-9]+?)\s.*\s(RNY[0-9])\s/){
			if (/^(t[0-9]+?)\s.*\s([0-9]+\.[0-9]+S)\s/){
				$id = $1;
				$anno = $2 . '-rRNA';
				if (/([ATCGN]+?)\s[I]+/){
					$seq = $1;
				}
			}
			elsif (/^(t[0-9]+?)\s.*\s([0-9]+S)\s/){
				$id = $1;
				$anno = $2 . '-rRNA';
				if (/([ATCGN]+?)\s[I]+/){
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
	}
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "rRNAdb-rRNA_Unmatch_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "rRNAdb-rRNA_Unmatch_Genome\t${anno}\t%.2f\n", $sum{$anno};
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

if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr_a;
	my %distr_b;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		if(/^(t[0-9]+?)\s.*\s(RNY[0-9])\s/){
			$id = $1;
			$anno = $2 . '-YRNA';
			if (/([ATCGN]+?)\s[I]+/){
				$seq = $1;
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
	}
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "YRNAdb-YRNA_Unmatch_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "YRNAdb-YRNA_Unmatch_Genome\t${anno}\t%.2f\n", $sum{$anno};
		}
		foreach $len(sort keys %distr_a){
			printf OUTPUT3 "YRNAdb-YRNA_Unmatch_Genome\t$len\t%.2f\n", $distr_a{$len};
		}
		foreach $anno(sort keys %distr_b){
			foreach $len(sort keys %{$distr_b{$anno}}){
				printf OUTPUT3 "YRNAdb-${anno}_Unmatch_Genome\t$len\t%0.2f\n", $distr_b{$anno}{$len};
			}
		}
	}
}

######summarize and annotate tRNA seqs: match-genome######
$i += 1;
{
	my %distr_pre;
	my %distr_pre_5_end;
	my %distr_pre_3_end;
	my %distr_mature;
	my %distr_mature_5_end;
	my %distr_mature_3_end;
	my %distr_mature_CCA_end;
	my %sum_pre;
	my %sum_pre_5_end;
	my %sum_pre_3_end;
	my %sum_mature;
	my %sum_mature_5_end;
	my %sum_mature_3_end;
	my %sum_mature_CCA_end;
	my %repeat_num;
	my $sums_pre = 0;
	my $sums_pre_5_end = 0;
	my $sums_pre_3_end = 0;
	my $sums_mature = 0;
	my $sums_mature_5_end = 0;
	my $sums_mature_3_end = 0;	
	my $sums_mature_CCA_end = 0;
	my $len;
	my $anno;
	my $seq;
	my $temp_anno;
	my $tRNA_len;
	my $tRNA_name;
	my $tRNA_codon;
	my $start_site;
	my $temp_end_str;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end';
				$anno = $temp_anno;
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_pre_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end';
				$anno = $temp_anno;
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_pre_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_3_end += ($reads{$id} / $repeat_num{$id});
			}else {
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$anno = $temp_anno;
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_pre += ($reads{$id} / $repeat_num{$id});
			$distr_pre{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	$i += 1;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5]; 
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end';
				$anno = $temp_anno;
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_mature_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site + 3)){
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end';
				$anno = $temp_anno;
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_mature_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_3_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_end_str = substr $seq, -3, 3;
				if ($temp_end_str eq 'CCA'){
					$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_CCA_end';
					$anno = $temp_anno;
					$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
					$sum_mature_CCA_end{$anno} += ($reads{$id} / $repeat_num{$id});
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
					$distr_mature_CCA_end{$len} += ($reads{$id} / $repeat_num{$id});
					$sums_mature_CCA_end += ($reads{$id} / $repeat_num{$id});
				}else {
					$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
					$anno = $temp_anno;
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				}
			}else {
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$anno = $temp_anno;
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_mature += ($reads{$id} / $repeat_num{$id});
			$distr_mature{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	if ($sums_pre > 0){
		printf OUTPUT2 "GtRNAdb-pre-tRNA_Match_Genome\t-\t%d\n", $sums_pre;
		foreach $anno(sort keys %sum_pre){
			printf OUTPUT2 "GtRNAdb-pre-tRNA_Match_Genome\t$anno\t%.2f\n", $sum_pre{$anno};
		}
		foreach $len (sort keys %distr_pre){
			printf OUTPUT3 "GtRNAdb-pre-tRNA_Match_Genome\t$len\t%.2f\n", $distr_pre{$len};
		}
	}
	if ($sums_pre_5_end > 0){
		printf OUTPUT2 "GtRNAdb-pre-tRNA_5_end_Match_Genome\t-\t%d\n", $sums_pre_5_end;
		foreach $anno(sort keys %sum_pre_5_end){
			printf OUTPUT2 "GtRNAdb-pre-tRNA_5_end_Match_Genome\t$anno\t%.2f\n", $sum_pre_5_end{$anno};
		}
		foreach $len (sort keys %distr_pre_5_end){
			printf OUTPUT3 "GtRNAdb-pre-tRNA_5_end_Match_Genome\t$len\t%.2f\n", $distr_pre_5_end{$len};
		}
	}
	if ($sums_pre_3_end > 0){
		printf OUTPUT2 "GtRNAdb-pre-tRNA_3_end_Match_Genome\t-\t%d\n", $sums_pre_3_end;
		foreach $anno(sort keys %sum_pre_3_end){
			printf OUTPUT2 "GtRNAdb-pre-tRNA_3_end_Match_Genome\t$anno\t%.2f\n", $sum_pre_3_end{$anno};
		}
		foreach $len (sort keys %distr_pre_3_end){
			printf OUTPUT3 "GtRNAdb-pre-tRNA_3_end_Match_Genome\t$len\t%.2f\n", $distr_pre_3_end{$len};
		}
	}
	if ($sums_mature > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_Match_Genome\t-\t%d\n", $sums_mature;
		foreach $anno(sort keys %sum_mature){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_Match_Genome\t$anno\t%.2f\n", $sum_mature{$anno};
		}
		foreach $len (sort keys %distr_mature){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_Match_Genome\t$len\t%.2f\n", $distr_mature{$len};
		}
	}
	if ($sums_mature_5_end > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_5_end_Match_Genome\t-\t%d\n", $sums_mature_5_end;
		foreach $anno(sort keys %sum_mature_5_end){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_5_end_Match_Genome\t$anno\t%.2f\n", $sum_mature_5_end{$anno};
		}
		foreach $len (sort keys %distr_mature_5_end){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_5_end_Match_Genome\t$len\t%.2f\n", $distr_mature_5_end{$len};
		}
	}
	if ($sums_mature_3_end > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_3_end_Match_Genome\t-\t%d\n", $sums_mature_3_end;
		foreach $anno(sort keys %sum_mature_3_end){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_3_end_Match_Genome\t$anno\t%.2f\n", $sum_mature_3_end{$anno};
		}
		foreach $len (sort keys %distr_mature_3_end){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_3_end_Match_Genome\t$len\t%.2f\n", $distr_mature_3_end{$len};
		}
	}
	if ($sums_mature_CCA_end > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_CCA_end_Match_Genome\t-\t%d\n", $sums_mature_CCA_end;
		foreach $anno(sort keys %sum_mature_CCA_end){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_CCA_end_Match_Genome\t$anno\t%.2f\n", $sum_mature_CCA_end{$anno};
		}
		foreach $len (sort keys %distr_mature_CCA_end){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_CCA_end_Match_Genome\t$len\t%.2f\n", $distr_mature_CCA_end{$len};
		}
	}
}

######summarize and annotate tRNA seqs: unmatch-genome######
$i += 1;
{
	my %distr_pre;
	my %distr_pre_5_end;
	my %distr_pre_3_end;
	my %distr_mature;
	my %distr_mature_5_end;
	my %distr_mature_3_end;
	my %distr_mature_CCA_end;
	my %sum_pre;
	my %sum_pre_5_end;
	my %sum_pre_3_end;
	my %sum_mature;
	my %sum_mature_5_end;
	my %sum_mature_3_end;
	my %sum_mature_CCA_end;
	my %repeat_num;
	my $sums_pre = 0;
	my $sums_pre_5_end = 0;
	my $sums_pre_3_end = 0;
	my $sums_mature = 0;
	my $sums_mature_5_end = 0;
	my $sums_mature_3_end = 0;	
	my $sums_mature_CCA_end = 0;
	my $len;
	my $anno;
	my $seq;
	my $temp_anno;
	my $tRNA_len;
	my $tRNA_name;
	my $tRNA_codon;
	my $start_site;
	my $temp_end_str;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end';
				$anno = $temp_anno;
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_pre_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end';
				$anno = $temp_anno;
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_pre_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_3_end += ($reads{$id} / $repeat_num{$id});
			}else {
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$anno = $temp_anno;
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_pre += ($reads{$id} / $repeat_num{$id});
			$distr_pre{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	$i += 1;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end';
				$anno = $temp_anno;
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_mature_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site + 3)){
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end';
				$anno = $temp_anno;
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_mature_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_3_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_end_str = substr $seq, -3, 3;
				if ($temp_end_str eq 'CCA'){
					$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_CCA_end';
					$anno = $temp_anno;
					$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
					$sum_mature_CCA_end{$anno} += ($reads{$id} / $repeat_num{$id});
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
					$distr_mature_CCA_end{$len} += ($reads{$id} / $repeat_num{$id});
					$sums_mature_CCA_end += ($reads{$id} / $repeat_num{$id});
				}else {
					$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
					$anno = $temp_anno;
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				}
			}else {
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$anno = $temp_anno;
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_mature += ($reads{$id} / $repeat_num{$id});
			$distr_mature{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	if ($sums_pre > 0){
		printf OUTPUT2 "GtRNAdb-pre-tRNA_Unmatch_Genome\t-\t%d\n", $sums_pre;
		foreach $anno(sort keys %sum_pre){
			printf OUTPUT2 "GtRNAdb-pre-tRNA_Unmatch_Genome\t$anno\t%.2f\n", $sum_pre{$anno};
		}
		foreach $len (sort keys %distr_pre){
			printf OUTPUT3 "GtRNAdb-pre-tRNA_Unmatch_Genome\t$len\t%.2f\n", $distr_pre{$len};
		}
	}
	if ($sums_pre_5_end > 0){
		printf OUTPUT2 "GtRNAdb-pre-tRNA_5_end_Unmatch_Genome\t-\t%d\n", $sums_pre_5_end;
		foreach $anno(sort keys %sum_pre_5_end){
			printf OUTPUT2 "GtRNAdb-pre-tRNA_5_end_Unmatch_Genome\t$anno\t%.2f\n", $sum_pre_5_end{$anno};
		}
		foreach $len (sort keys %distr_pre_5_end){
			printf OUTPUT3 "GtRNAdb-pre-tRNA_5_end_Unmatch_Genome\t$len\t%.2f\n", $distr_pre_5_end{$len};
		}
	}
	if ($sums_pre_3_end > 0){
		printf OUTPUT2 "GtRNAdb-pre-tRNA_3_end_Unmatch_Genome\t-\t%d\n", $sums_pre_3_end;
		foreach $anno(sort keys %sum_pre_3_end){
			printf OUTPUT2 "GtRNAdb-pre-tRNA_3_end_Unmatch_Genome\t$anno\t%.2f\n", $sum_pre_3_end{$anno};
		}
		foreach $len (sort keys %distr_pre_3_end){
			printf OUTPUT3 "GtRNAdb-pre-tRNA_3_end_Unmatch_Genome\t$len\t%.2f\n", $distr_pre_3_end{$len};
		}
	}
	if ($sums_mature > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_Unmatch_Genome\t-\t%d\n", $sums_mature;
		foreach $anno(sort keys %sum_mature){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature{$anno};
		}
		foreach $len (sort keys %distr_mature){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_Unmatch_Genome\t$len\t%.2f\n", $distr_mature{$len};
		}
	}
	if ($sums_mature_5_end > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_5_end_Unmatch_Genome\t-\t%d\n", $sums_mature_5_end;
		foreach $anno(sort keys %sum_mature_5_end){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_5_end_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature_5_end{$anno};
		}
		foreach $len (sort keys %distr_mature_5_end){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_5_end_Unmatch_Genome\t$len\t%.2f\n", $distr_mature_5_end{$len};
		}
	}
	if ($sums_mature_3_end > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_3_end_Unmatch_Genome\t-\t%d\n", $sums_mature_3_end;
		foreach $anno(sort keys %sum_mature_3_end){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_3_end_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature_3_end{$anno};
		}
		foreach $len (sort keys %distr_mature_3_end){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_3_end_Unmatch_Genome\t$len\t%.2f\n", $distr_mature_3_end{$len};
		}
	}
	if ($sums_mature_CCA_end > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_CCA_end_Unmatch_Genome\t-\t%d\n", $sums_mature_CCA_end;
		foreach $anno(sort keys %sum_mature_CCA_end){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_CCA_end_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature_CCA_end{$anno};
		}
		foreach $len (sort keys %distr_mature_CCA_end){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_CCA_end_Unmatch_Genome\t$len\t%.2f\n", $distr_mature_CCA_end{$len};
		}
	}
}

######summarize and annotate mt_tRNA seqs: match-genome######
$i += 1;
{
	my %distr_pre;
	my %distr_pre_5_end;
	my %distr_pre_3_end;
	my %distr_mature;
	my %distr_mature_5_end;
	my %distr_mature_3_end;
	my %distr_mature_CCA_end;
	my %sum_pre;
	my %sum_pre_5_end;
	my %sum_pre_3_end;
	my %sum_mature;
	my %sum_mature_5_end;
	my %sum_mature_3_end;
	my %sum_mature_CCA_end;
	my %repeat_num;
	my $sums_pre = 0;
	my $sums_pre_5_end = 0;
	my $sums_pre_3_end = 0;
	my $sums_mature = 0;
	my $sums_mature_5_end = 0;
	my $sums_mature_3_end = 0;	
	my $sums_mature_CCA_end = 0;
	my $len;
	my $anno;
	my $seq;
	my $temp_anno;
	my $tRNA_len;
	my $tRNA_name;
	my $tRNA_codon;
	my $start_site;
	my $temp_end_str;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end';
				$anno = $temp_anno;
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_pre_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end';
				$anno = $temp_anno;
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_pre_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_3_end += ($reads{$id} / $repeat_num{$id});
			}else {
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$anno = $temp_anno;
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_pre += ($reads{$id} / $repeat_num{$id});
			$distr_pre{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	$i += 1;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5]; 
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end';
				$anno = $temp_anno;
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_mature_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site + 3)){
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end';
				$anno = $temp_anno;
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_mature_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_3_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_end_str = substr $seq, -3, 3;
				if ($temp_end_str eq 'CCA'){
					$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_CCA_end';
					$anno = $temp_anno;
					$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
					$sum_mature_CCA_end{$anno} += ($reads{$id} / $repeat_num{$id});
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
					$distr_mature_CCA_end{$len} += ($reads{$id} / $repeat_num{$id});
					$sums_mature_CCA_end += ($reads{$id} / $repeat_num{$id});
				}else {
					$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
					$anno = $temp_anno;
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				}
			}else {
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$anno = $temp_anno;
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_mature += ($reads{$id} / $repeat_num{$id});
			$distr_mature{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	if ($sums_pre > 0){
		printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA_Match_Genome\t-\t%d\n", $sums_pre;
		foreach $anno(sort keys %sum_pre){
			printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA_Match_Genome\t$anno\t%.2f\n", $sum_pre{$anno};
		}
		foreach $len (sort keys %distr_pre){
			printf OUTPUT3 "mitotRNAdb-pre-mt_tRNA_Match_Genome\t$len\t%.2f\n", $distr_pre{$len};
		}
	}
	if ($sums_pre_5_end > 0){
		printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA_5_end_Match_Genome\t-\t%d\n", $sums_pre_5_end;
		foreach $anno(sort keys %sum_pre_5_end){
			printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA_5_end_Match_Genome\t$anno\t%.2f\n", $sum_pre_5_end{$anno};
		}
		foreach $len (sort keys %distr_pre_5_end){
			printf OUTPUT3 "mitotRNAdb-pre-mt_tRNA_5_end_Match_Genome\t$len\t%.2f\n", $distr_pre_5_end{$len};
		}
	}
	if ($sums_pre_3_end > 0){
		printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA_3_end_Match_Genome\t-\t%d\n", $sums_pre_3_end;
		foreach $anno(sort keys %sum_pre_3_end){
			printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA_3_end_Match_Genome\t$anno\t%.2f\n", $sum_pre_3_end{$anno};
		}
		foreach $len (sort keys %distr_pre_3_end){
			printf OUTPUT3 "mitotRNAdb-pre-mt_tRNA_3_end_Match_Genome\t$len\t%.2f\n", $distr_pre_3_end{$len};
		}
	}
	if ($sums_mature > 0){
		printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_Match_Genome\t-\t%d\n", $sums_mature;
		foreach $anno(sort keys %sum_mature){
			printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_Match_Genome\t$anno\t%.2f\n", $sum_mature{$anno};
		}
		foreach $len (sort keys %distr_mature){
			printf OUTPUT3 "mitotRNAdb-mature-mt_tRNA_Match_Genome\t$len\t%.2f\n", $distr_mature{$len};
		}
	}
	if ($sums_mature_5_end > 0){
		printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_5_end_Match_Genome\t-\t%d\n", $sums_mature_5_end;
		foreach $anno(sort keys %sum_mature_5_end){
			printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_5_end_Match_Genome\t$anno\t%.2f\n", $sum_mature_5_end{$anno};
		}
		foreach $len (sort keys %distr_mature_5_end){
			printf OUTPUT3 "mitotRNAdb-mature-mt_tRNA_5_end_Match_Genome\t$len\t%.2f\n", $distr_mature_5_end{$len};
		}
	}
	if ($sums_mature_3_end > 0){
		printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_3_end_Match_Genome\t-\t%d\n", $sums_mature_3_end;
		foreach $anno(sort keys %sum_mature_3_end){
			printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_3_end_Match_Genome\t$anno\t%.2f\n", $sum_mature_3_end{$anno};
		}
		foreach $len (sort keys %distr_mature_3_end){
			printf OUTPUT3 "mitotRNAdb-mature-mt_tRNA_3_end_Match_Genome\t$len\t%.2f\n", $distr_mature_3_end{$len};
		}
	}
	if ($sums_mature_CCA_end > 0){
		printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_CCA_end_Match_Genome\t-\t%d\n", $sums_mature_CCA_end;
		foreach $anno(sort keys %sum_mature_CCA_end){
			printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_CCA_end_Match_Genome\t$anno\t%.2f\n", $sum_mature_CCA_end{$anno};
		}
		foreach $len (sort keys %distr_mature_CCA_end){
			printf OUTPUT3 "mitotRNAdb-mature-mt_tRNA_CCA_end_Match_Genome\t$len\t%.2f\n", $distr_mature_CCA_end{$len};
		}
	}
}

######summarize and annotate mt_tRNA seqs: unmatch-genome######
$i += 1;
{
	my %distr_pre;
	my %distr_pre_5_end;
	my %distr_pre_3_end;
	my %distr_mature;
	my %distr_mature_5_end;
	my %distr_mature_3_end;
	my %distr_mature_CCA_end;
	my %sum_pre;
	my %sum_pre_5_end;
	my %sum_pre_3_end;
	my %sum_mature;
	my %sum_mature_5_end;
	my %sum_mature_3_end;
	my %sum_mature_CCA_end;
	my %repeat_num;
	my $sums_pre = 0;
	my $sums_pre_5_end = 0;
	my $sums_pre_3_end = 0;
	my $sums_mature = 0;
	my $sums_mature_5_end = 0;
	my $sums_mature_3_end = 0;	
	my $sums_mature_CCA_end = 0;
	my $len;
	my $anno;
	my $seq;
	my $temp_anno;
	my $tRNA_len;
	my $tRNA_name;
	my $tRNA_codon;
	my $start_site;
	my $temp_end_str;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end';
				$anno = $temp_anno;
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_pre_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end';
				$anno = $temp_anno;
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_pre_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_3_end += ($reads{$id} / $repeat_num{$id});
			}else {
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$anno = $temp_anno;
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_pre += ($reads{$id} / $repeat_num{$id});
			$distr_pre{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	$i += 1;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end';
				$anno = $temp_anno;
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_mature_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site + 3)){
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end';
				$anno = $temp_anno;
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$sum_mature_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_3_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_end_str = substr $seq, -3, 3;
				if ($temp_end_str eq 'CCA'){
					$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_CCA_end';
					$anno = $temp_anno;
					$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
					$sum_mature_CCA_end{$anno} += ($reads{$id} / $repeat_num{$id});
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
					$distr_mature_CCA_end{$len} += ($reads{$id} / $repeat_num{$id});
					$sums_mature_CCA_end += ($reads{$id} / $repeat_num{$id});
				}else {
					$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
					$anno = $temp_anno;
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				}
			}else {
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon;
				$anno = $temp_anno;
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_mature += ($reads{$id} / $repeat_num{$id});
			$distr_mature{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	if ($sums_pre > 0){
		printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA_Unmatch_Genome\t-\t%d\n", $sums_pre;
		foreach $anno(sort keys %sum_pre){
			printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA_Unmatch_Genome\t$anno\t%.2f\n", $sum_pre{$anno};
		}
		foreach $len (sort keys %distr_pre){
			printf OUTPUT3 "mitotRNAdb-pre-mt_tRNA_Unmatch_Genome\t$len\t%.2f\n", $distr_pre{$len};
		}
	}
	if ($sums_pre_5_end > 0){
		printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA_5_end_Unmatch_Genome\t-\t%d\n", $sums_pre_5_end;
		foreach $anno(sort keys %sum_pre_5_end){
			printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA_5_end_Unmatch_Genome\t$anno\t%.2f\n", $sum_pre_5_end{$anno};
		}
		foreach $len (sort keys %distr_pre_5_end){
			printf OUTPUT3 "mitotRNAdb-pre-mt_tRNA_5_end_Unmatch_Genome\t$len\t%.2f\n", $distr_pre_5_end{$len};
		}
	}
	if ($sums_pre_3_end > 0){
		printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA_3_end_Unmatch_Genome\t-\t%d\n", $sums_pre_3_end;
		foreach $anno(sort keys %sum_pre_3_end){
			printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA_3_end_Unmatch_Genome\t$anno\t%.2f\n", $sum_pre_3_end{$anno};
		}
		foreach $len (sort keys %distr_pre_3_end){
			printf OUTPUT3 "mitotRNAdb-pre-mt_tRNA_3_end_Unmatch_Genome\t$len\t%.2f\n", $distr_pre_3_end{$len};
		}
	}
	if ($sums_mature > 0){
		printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_Unmatch_Genome\t-\t%d\n", $sums_mature;
		foreach $anno(sort keys %sum_mature){
			printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature{$anno};
		}
		foreach $len (sort keys %distr_mature){
			printf OUTPUT3 "mitotRNAdb-mature-mt_tRNA_Unmatch_Genome\t$len\t%.2f\n", $distr_mature{$len};
		}
	}
	if ($sums_mature_5_end > 0){
		printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_5_end_Unmatch_Genome\t-\t%d\n", $sums_mature_5_end;
		foreach $anno(sort keys %sum_mature_5_end){
			printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_5_end_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature_5_end{$anno};
		}
		foreach $len (sort keys %distr_mature_5_end){
			printf OUTPUT3 "mitotRNAdb-mature-mt_tRNA_5_end_Unmatch_Genome\t$len\t%.2f\n", $distr_mature_5_end{$len};
		}
	}
	if ($sums_mature_3_end > 0){
		printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_3_end_Unmatch_Genome\t-\t%d\n", $sums_mature_3_end;
		foreach $anno(sort keys %sum_mature_3_end){
			printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_3_end_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature_3_end{$anno};
		}
		foreach $len (sort keys %distr_mature_3_end){
			printf OUTPUT3 "mitotRNAdb-mature-mt_tRNA_3_end_Unmatch_Genome\t$len\t%.2f\n", $distr_mature_3_end{$len};
		}
	}
	if ($sums_mature_CCA_end > 0){
		printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_CCA_end_Unmatch_Genome\t-\t%d\n", $sums_mature_CCA_end;
		foreach $anno(sort keys %sum_mature_CCA_end){
			printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA_CCA_end_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature_CCA_end{$anno};
		}
		foreach $len (sort keys %distr_mature_CCA_end){
			printf OUTPUT3 "mitotRNAdb-mature-mt_tRNA_CCA_end_Unmatch_Genome\t$len\t%.2f\n", $distr_mature_CCA_end{$len};
		}
	}
}

######summarize and annotate using piRNA database: match-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
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
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "piRNAdb-piRNA_Match_Genome\t-\t%d\n", $sums;
		foreach $len (sort keys %distr){
			printf OUTPUT3 "piRNAdb-piRNA_Match_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
}

######summarize and annotate using piRNA database: unmatch-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
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
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "piRNAdb-piRNA_Unmatch_Genome\t-\t%d\n", $sums;
		foreach $len (sort keys %distr){
			printf OUTPUT3 "piRNAdb-piRNA_Unmatch_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
}

######summarize and annotate using ensembl database: match-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
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
			$sums += ($reads{$id} / $repeat_num{$id});
			$sum{$anno} += ($reads{$id} / $repeat_num{$id});
			$len = $lens{$id};
			$distr{$anno}{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
	}
	close $fh;
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
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
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
			$sums += ($reads{$id} / $repeat_num{$id});
			$sum{$anno} += ($reads{$id} / $repeat_num{$id});
			$len = $lens{$id};
			$distr{$anno}{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
	}
	close $fh;
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
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
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
	close $fh;
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
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
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
	close $fh;
	foreach $anno(sort keys %sum){
		printf OUTPUT2 "Rfam-${anno}_Unmatch_Genome\t-\t%d\n", $sum{$anno};
	}
	foreach	my $key1(sort keys %distr){
		foreach my $key2 (sort keys %{$distr{$key1}}){
			printf OUTPUT3 "Rfam-${key1}_Unmatch_Genome\t$key2\t%.2f\n", $distr{$key1}{$key2};
		}
	}	
}

######summarize and annotate miRNA-antisense seqs: match-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		($id, $anno, $seq) = (split /\t/)[0, 3, 5];
		$anno = (split /\s+/, $anno)[0] . '-antisense';
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
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "miRBase-miRNA-antisense_Match_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "miRBase-miRNA-antisense_Match_Genome\t$anno\t%.2f\n", $sum{$anno};
		}
		foreach $len (sort keys %distr){
			printf OUTPUT3 "miRBase-miRNA-antisense_Match_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
}

######summarize and annotate miRNA-antisense seqs: unmatch-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		($id, $anno, $seq) = (split /\t/)[0, 3, 5];
		$anno = (split /\s+/, $anno)[0] . '-antisense';
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
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "miRBase-miRNA-antisense_Unmatch_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "miRBase-miRNA-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum{$anno};
		}
		foreach $len (sort keys %distr){
			printf OUTPUT3 "miRBase-miRNA-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
}

######summarize and annotate rRNA-antisense seqs: match-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr_a;
	my %distr_b;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		if(!/^(t[0-9]+?)\s.*\s(RNY[0-9])\s/){
			if (/^(t[0-9]+?)\s.*\s([0-9]+\.[0-9]+S)\s/){
				$id = $1;
				$anno = $2 . '-rRNA-antisense';
				if (/([ATCGN]+?)\s[I]+/){
					$seq = $1;
				}
			}
			elsif (/^(t[0-9]+?)\s.*\s([0-9]+S)\s/){
				$id = $1;
				$anno = $2 . '-rRNA-antisense';
				if (/([ATCGN]+?)\s[I]+/){
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
	}
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "rRNAdb-rRNA-antisense_Match_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "rRNAdb-rRNA-antisense_Match_Genome\t${anno}\t%.2f\n", $sum{$anno};
		}
		foreach $len(sort keys %distr_a){
			printf OUTPUT3 "rRNAdb-rRNA-antisense_Match_Genome\t$len\t%.2f\n", $distr_a{$len};
		}
		foreach $anno(sort keys %distr_b){
			foreach $len(sort keys %{$distr_b{$anno}}){
				printf OUTPUT3 "rRNAdb-${anno}_Match_Genome\t$len\t%0.2f\n", $distr_b{$anno}{$len};
			}
		}
	}
}

if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr_a;
	my %distr_b;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		if(/^(t[0-9]+?)\s.*\s(RNY[0-9])\s/){
			$id = $1;
			$anno = $2 . '-YRNA-antisense';
			if (/([ATCGN]+?)\s[I]+/){
				$seq = $1;
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
	}
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "YRNAdb-YRNA-antisense_Match_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "YRNAdb-YRNA-antisense_Match_Genome\t${anno}\t%.2f\n", $sum{$anno};
		}
		foreach $len(sort keys %distr_a){
			printf OUTPUT3 "YRNAdb-YRNA-antisense_Match_Genome\t$len\t%.2f\n", $distr_a{$len};
		}
		foreach $anno(sort keys %distr_b){
			foreach $len(sort keys %{$distr_b{$anno}}){
				printf OUTPUT3 "YRNAdb-${anno}_Match_Genome\t$len\t%0.2f\n", $distr_b{$anno}{$len};
			}
		}
	}
}

######summarize and annotate rRNA-antisense seqs: unmatch-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr_a;
	my %distr_b;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		if(!/^(t[0-9]+?)\s.*\s(RNY[0-9])\s/){
			if (/^(t[0-9]+?)\s.*\s([0-9]+\.[0-9]+S)\s/){
				$id = $1;
				$anno = $2 . '-rRNA-antisense';
				if (/([ATCGN]+?)\s[I]+/){
					$seq = $1;
				}
			}
			elsif (/^(t[0-9]+?)\s.*\s([0-9]+S)\s/){
				$id = $1;
				$anno = $2 . '-rRNA-antisense';
				if (/([ATCGN]+?)\s[I]+/){
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
	}
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "rRNAdb-rRNA-antisense_Unmatch_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "rRNAdb-rRNA-antisense_Unmatch_Genome\t${anno}\t%.2f\n", $sum{$anno};
		}
		foreach $len(sort keys %distr_a){
			printf OUTPUT3 "rRNAdb-rRNA-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_a{$len};
		}
		foreach $anno(sort keys %distr_b){
			foreach $len(sort keys %{$distr_b{$anno}}){
				printf OUTPUT3 "rRNAdb-${anno}_Unmatch_Genome\t$len\t%0.2f\n", $distr_b{$anno}{$len};
			}
		}
	}
}

if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr_a;
	my %distr_b;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		if(/^(t[0-9]+?)\s.*\s(RNY[0-9])\s/){
			$id = $1;
			$anno = $2 . '-YRNA-antisense';
			if (/([ATCGN]+?)\s[I]+/){
				$seq = $1;
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
	}
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "YRNAdb-YRNA-antisense_Unmatch_Genome\t-\t%d\n", $sums;
		foreach $anno(sort keys %sum){
			printf OUTPUT2 "YRNAdb-YRNA-antisense_Unmatch_Genome\t${anno}\t%.2f\n", $sum{$anno};
		}
		foreach $len(sort keys %distr_a){
			printf OUTPUT3 "YRNAdb-YRNA-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_a{$len};
		}
		foreach $anno(sort keys %distr_b){
			foreach $len(sort keys %{$distr_b{$anno}}){
				printf OUTPUT3 "YRNAdb-${anno}_Unmatch_Genome\t$len\t%0.2f\n", $distr_b{$anno}{$len};
			}
		}
	}
}

######summarize and annotate tRNA-antisense seqs: match-genome######
$i += 1;
{
	my %distr_pre;
	my %distr_pre_5_end;
	my %distr_pre_3_end;
	my %distr_mature;
	my %distr_mature_5_end;
	my %distr_mature_3_end;
	my %distr_mature_CCA_end;
	my %sum_pre;
	my %sum_pre_5_end;
	my %sum_pre_3_end;
	my %sum_mature;
	my %sum_mature_5_end;
	my %sum_mature_3_end;
	my %sum_mature_CCA_end;
	my %repeat_num;
	my $sums_pre = 0;
	my $sums_pre_5_end = 0;
	my $sums_pre_3_end = 0;
	my $sums_mature = 0;
	my $sums_mature_5_end = 0;
	my $sums_mature_3_end = 0;	
	my $sums_mature_CCA_end = 0;
	my $len;
	my $anno;
	my $seq;
	my $temp_anno;
	my $tRNA_len;
	my $tRNA_name;
	my $tRNA_codon;
	my $start_site;
	my $temp_end_str;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_pre_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_pre_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_3_end += ($reads{$id} / $repeat_num{$id});
			}else {
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$anno = $temp_anno;
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_pre += ($reads{$id} / $repeat_num{$id});
			$distr_pre{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	$i += 1;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_mature_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site + 3)){
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_mature_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_3_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_end_str = substr $seq, -3, 3;
				if ($temp_end_str eq 'CCA'){
					$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_CCA_end-antisense';
					$anno = $temp_anno;
					$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
					$sum_mature_CCA_end{$anno} += ($reads{$id} / $repeat_num{$id});
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
					$distr_mature_CCA_end{$len} += ($reads{$id} / $repeat_num{$id});
					$sums_mature_CCA_end += ($reads{$id} / $repeat_num{$id});
				}else {
					$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
					$anno = $temp_anno;
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				}
			}else {
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$anno = $temp_anno;
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_mature += ($reads{$id} / $repeat_num{$id});
			$distr_mature{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	if ($sums_pre > 0){
		printf OUTPUT2 "GtRNAdb-pre-tRNA-antisense_Match_Genome\t-\t%d\n", $sums_pre;
		foreach $anno(sort keys %sum_pre){
			printf OUTPUT2 "GtRNAdb-pre-tRNA-antisense_Match_Genome\t$anno\t%.2f\n", $sum_pre{$anno};
		}
		foreach $len (sort keys %distr_pre){
			printf OUTPUT3 "GtRNAdb-pre-tRNA-antisense_Match_Genome\t$len\t%.2f\n", $distr_pre{$len};
		}
	}
	if ($sums_pre_5_end > 0){
		printf OUTPUT2 "GtRNAdb-pre-tRNA_5_end-antisense_Match_Genome\t-\t%d\n", $sums_pre_5_end;
		foreach $anno(sort keys %sum_pre_5_end){
			printf OUTPUT2 "GtRNAdb-pre-tRNA_5_end-antisense_Match_Genome\t$anno\t%.2f\n", $sum_pre_5_end{$anno};
		}
		foreach $len (sort keys %distr_pre_5_end){
			printf OUTPUT3 "GtRNAdb-pre-tRNA_5_end-antisense_Match_Genome\t$len\t%.2f\n", $distr_pre_5_end{$len};
		}
	}
	if ($sums_pre_3_end > 0){
		printf OUTPUT2 "GtRNAdb-pre-tRNA_3_end-antisense_Match_Genome\t-\t%d\n", $sums_pre_3_end;
		foreach $anno(sort keys %sum_pre_3_end){
			printf OUTPUT2 "GtRNAdb-pre-tRNA_3_end-antisense_Match_Genome\t$anno\t%.2f\n", $sum_pre_3_end{$anno};
		}
		foreach $len (sort keys %distr_pre_3_end){
			printf OUTPUT3 "GtRNAdb-pre-tRNA_3_end-antisense_Match_Genome\t$len\t%.2f\n", $distr_pre_3_end{$len};
		}
	}
	if ($sums_mature > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA-antisense_Match_Genome\t-\t%d\n", $sums_mature;
		foreach $anno(sort keys %sum_mature){
			printf OUTPUT2 "GtRNAdb-mature-tRNA-antisense_Match_Genome\t$anno\t%.2f\n", $sum_mature{$anno};
		}
		foreach $len (sort keys %distr_mature){
			printf OUTPUT3 "GtRNAdb-mature-tRNA-antisense_Match_Genome\t$len\t%.2f\n", $distr_mature{$len};
		}
	}
	if ($sums_mature_5_end > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_5_end-antisense_Match_Genome\t-\t%d\n", $sums_mature_5_end;
		foreach $anno(sort keys %sum_mature_5_end){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_5_end-antisense_Match_Genome\t$anno\t%.2f\n", $sum_mature_5_end{$anno};
		}
		foreach $len (sort keys %distr_mature_5_end){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_5_end-antisense_Match_Genome\t$len\t%.2f\n", $distr_mature_5_end{$len};
		}
	}
	if ($sums_mature_3_end > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_3_end-antisense_Match_Genome\t-\t%d\n", $sums_mature_3_end;
		foreach $anno(sort keys %sum_mature_3_end){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_3_end-antisense_Match_Genome\t$anno\t%.2f\n", $sum_mature_3_end{$anno};
		}
		foreach $len (sort keys %distr_mature_3_end){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_3_end-antisense_Match_Genome\t$len\t%.2f\n", $distr_mature_3_end{$len};
		}
	}
	if ($sums_mature_CCA_end > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_CCA_end-antisense_Match_Genome\t-\t%d\n", $sums_mature_CCA_end;
		foreach $anno(sort keys %sum_mature_CCA_end){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_CCA_end-antisense_Match_Genome\t$anno\t%.2f\n", $sum_mature_CCA_end{$anno};
		}
		foreach $len (sort keys %distr_mature_CCA_end){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_CCA_end-antisense_Match_Genome\t$len\t%.2f\n", $distr_mature_CCA_end{$len};
		}
	}
}

######summarize and annotate tRNA-antisense seqs: unmatch-genome######
$i += 1;
{
	my %distr_pre;
	my %distr_pre_5_end;
	my %distr_pre_3_end;
	my %distr_mature;
	my %distr_mature_5_end;
	my %distr_mature_3_end;
	my %distr_mature_CCA_end;
	my %sum_pre;
	my %sum_pre_5_end;
	my %sum_pre_3_end;
	my %sum_mature;
	my %sum_mature_5_end;
	my %sum_mature_3_end;
	my %sum_mature_CCA_end;
	my %repeat_num;
	my $sums_pre = 0;
	my $sums_pre_5_end = 0;
	my $sums_pre_3_end = 0;
	my $sums_mature = 0;
	my $sums_mature_5_end = 0;
	my $sums_mature_3_end = 0;	
	my $sums_mature_CCA_end = 0;
	my $len;
	my $anno;
	my $seq;
	my $temp_anno;
	my $tRNA_len;
	my $tRNA_name;
	my $tRNA_codon;
	my $start_site;
	my $temp_end_str;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_pre_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_pre_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_3_end += ($reads{$id} / $repeat_num{$id});
			}else {
				$temp_anno = 'pre-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$anno = $temp_anno;
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_pre += ($reads{$id} / $repeat_num{$id});
			$distr_pre{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	$i += 1;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_mature';
				$sum_mature_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site + 3)){
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_mature_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_3_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_end_str = substr $seq, -3, 3;
				if ($temp_end_str eq 'CCA'){
					$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_CCA_end-antisense';
					$anno = $temp_anno;
					$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
					$sum_mature_CCA_end{$anno} += ($reads{$id} / $repeat_num{$id});
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
					$distr_mature_CCA_end{$len} += ($reads{$id} / $repeat_num{$id});
					$sums_mature_CCA_end += ($reads{$id} / $repeat_num{$id});
				}else {
					$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
					$anno = $temp_anno;
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				}
			}else {
				$temp_anno = 'mature-tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$anno = $temp_anno;
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_mature += ($reads{$id} / $repeat_num{$id});
			$distr_mature{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	if ($sums_pre > 0){
		printf OUTPUT2 "GtRNAdb-pre-tRNA-antisense_Unmatch_Genome\t-\t%d\n", $sums_pre;
		foreach $anno(sort keys %sum_pre){
			printf OUTPUT2 "GtRNAdb-pre-tRNA-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_pre{$anno};
		}
		foreach $len (sort keys %distr_pre){
			printf OUTPUT3 "GtRNAdb-pre-tRNA-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_pre{$len};
		}
	}
	if ($sums_pre_5_end > 0){
		printf OUTPUT2 "GtRNAdb-pre-tRNA_5_end-antisense_Unmatch_Genome\t-\t%d\n", $sums_pre_5_end;
		foreach $anno(sort keys %sum_pre_5_end){
			printf OUTPUT2 "GtRNAdb-pre-tRNA_5_end-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_pre_5_end{$anno};
		}
		foreach $len (sort keys %distr_pre_5_end){
			printf OUTPUT3 "GtRNAdb-pre-tRNA_5_end-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_pre_5_end{$len};
		}
	}
	if ($sums_pre_3_end > 0){
		printf OUTPUT2 "GtRNAdb-pre-tRNA_3_end-antisense_Unmatch_Genome\t-\t%d\n", $sums_pre_3_end;
		foreach $anno(sort keys %sum_pre_3_end){
			printf OUTPUT2 "GtRNAdb-pre-tRNA_3_end-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_pre_3_end{$anno};
		}
		foreach $len (sort keys %distr_pre_3_end){
			printf OUTPUT3 "GtRNAdb-pre-tRNA_3_end-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_pre_3_end{$len};
		}
	}
	if ($sums_mature > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA-antisense_Unmatch_Genome\t-\t%d\n", $sums_mature;
		foreach $anno(sort keys %sum_mature){
			printf OUTPUT2 "GtRNAdb-mature-tRNA-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature{$anno};
		}
		foreach $len (sort keys %distr_mature){
			printf OUTPUT3 "GtRNAdb-mature-tRNA-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_mature{$len};
		}
	}
	if ($sums_mature_5_end > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_5_end-antisense_Unmatch_Genome\t-\t%d\n", $sums_mature_5_end;
		foreach $anno(sort keys %sum_mature_5_end){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_5_end-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature_5_end{$anno};
		}
		foreach $len (sort keys %distr_mature_5_end){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_5_end-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_mature_5_end{$len};
		}
	}
	if ($sums_mature_3_end > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_3_end-antisense_Unmatch_Genome\t-\t%d\n", $sums_mature_3_end;
		foreach $anno(sort keys %sum_mature_3_end){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_3_end-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature_3_end{$anno};
		}
		foreach $len (sort keys %distr_mature_3_end){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_3_end-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_mature_3_end{$len};
		}
	}
	if ($sums_mature_CCA_end > 0){
		printf OUTPUT2 "GtRNAdb-mature-tRNA_CCA_end-antisense_Unmatch_Genome\t-\t%d\n", $sums_mature_CCA_end;
		foreach $anno(sort keys %sum_mature_CCA_end){
			printf OUTPUT2 "GtRNAdb-mature-tRNA_CCA_end-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature_CCA_end{$anno};
		}
		foreach $len (sort keys %distr_mature_CCA_end){
			printf OUTPUT3 "GtRNAdb-mature-tRNA_CCA_end-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_mature_CCA_end{$len};
		}
	}
}

######summarize and annotate mt_tRNA-antisense seqs: match-genome######
$i += 1;
{
	my %distr_pre;
	my %distr_pre_5_end;
	my %distr_pre_3_end;
	my %distr_mature;
	my %distr_mature_5_end;
	my %distr_mature_3_end;
	my %distr_mature_CCA_end;
	my %sum_pre;
	my %sum_pre_5_end;
	my %sum_pre_3_end;
	my %sum_mature;
	my %sum_mature_5_end;
	my %sum_mature_3_end;
	my %sum_mature_CCA_end;
	my %repeat_num;
	my $sums_pre = 0;
	my $sums_pre_5_end = 0;
	my $sums_pre_3_end = 0;
	my $sums_mature = 0;
	my $sums_mature_5_end = 0;
	my $sums_mature_3_end = 0;	
	my $sums_mature_CCA_end = 0;
	my $len;
	my $anno;
	my $seq;
	my $temp_anno;
	my $tRNA_len;
	my $tRNA_name;
	my $tRNA_codon;
	my $start_site;
	my $temp_end_str;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_pre_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_pre_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_3_end += ($reads{$id} / $repeat_num{$id});
			}else {
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$anno = $temp_anno;
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_pre += ($reads{$id} / $repeat_num{$id});
			$distr_pre{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	$i += 1;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_mature_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site + 3)){
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_mature_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_3_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_end_str = substr $seq, -3, 3;
				if ($temp_end_str eq 'CCA'){
					$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_CCA_end-antisense';
					$anno = $temp_anno;
					$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
					$sum_mature_CCA_end{$anno} += ($reads{$id} / $repeat_num{$id});
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
					$distr_mature_CCA_end{$len} += ($reads{$id} / $repeat_num{$id});
					$sums_mature_CCA_end += ($reads{$id} / $repeat_num{$id});
				}else {
					$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
					$anno = $temp_anno;
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				}
			}else {
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$anno = $temp_anno;
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_mature += ($reads{$id} / $repeat_num{$id});
			$distr_mature{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	if ($sums_pre > 0){
		printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA-antisense_Match_Genome\t-\t%d\n", $sums_pre;
		foreach $anno(sort keys %sum_pre){
			printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA-antisense_Match_Genome\t$anno\t%.2f\n", $sum_pre{$anno};
		}
		foreach $len (sort keys %distr_pre){
			printf OUTPUT3 "mitotRNAdb-pre-mt_tRNA-antisense_Match_Genome\t$len\t%.2f\n", $distr_pre{$len};
		}
	}
	if ($sums_pre_5_end > 0){
		printf OUTPUT2 "mitotRNAdb-pre-tRNA_5_end-antisense_Match_Genome\t-\t%d\n", $sums_pre_5_end;
		foreach $anno(sort keys %sum_pre_5_end){
			printf OUTPUT2 "mitotRNAdb-pre-tRNA_5_end-antisense_Match_Genome\t$anno\t%.2f\n", $sum_pre_5_end{$anno};
		}
		foreach $len (sort keys %distr_pre_5_end){
			printf OUTPUT3 "mitotRNAdb-pre-tRNA_5_end-antisense_Match_Genome\t$len\t%.2f\n", $distr_pre_5_end{$len};
		}
	}
	if ($sums_pre_3_end > 0){
		printf OUTPUT2 "mitotRNAdb-pre-tRNA_3_end-antisense_Match_Genome\t-\t%d\n", $sums_pre_3_end;
		foreach $anno(sort keys %sum_pre_3_end){
			printf OUTPUT2 "mitotRNAdb-pre-tRNA_3_end-antisense_Match_Genome\t$anno\t%.2f\n", $sum_pre_3_end{$anno};
		}
		foreach $len (sort keys %distr_pre_3_end){
			printf OUTPUT3 "mitotRNAdb-pre-tRNA_3_end-antisense_Match_Genome\t$len\t%.2f\n", $distr_pre_3_end{$len};
		}
	}
	if ($sums_mature > 0){
		printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA-antisense_Match_Genome\t-\t%d\n", $sums_mature;
		foreach $anno(sort keys %sum_mature){
			printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA-antisense_Match_Genome\t$anno\t%.2f\n", $sum_mature{$anno};
		}
		foreach $len (sort keys %distr_mature){
			printf OUTPUT3 "mitotRNAdb-mature-mt_tRNA-antisense_Match_Genome\t$len\t%.2f\n", $distr_mature{$len};
		}
	}
	if ($sums_mature_5_end > 0){
		printf OUTPUT2 "mitotRNAdb-mature-tRNA_5_end-antisense_Match_Genome\t-\t%d\n", $sums_mature_5_end;
		foreach $anno(sort keys %sum_mature_5_end){
			printf OUTPUT2 "mitotRNAdb-mature-tRNA_5_end-antisense_Match_Genome\t$anno\t%.2f\n", $sum_mature_5_end{$anno};
		}
		foreach $len (sort keys %distr_mature_5_end){
			printf OUTPUT3 "mitotRNAdb-mature-tRNA_5_end-antisense_Match_Genome\t$len\t%.2f\n", $distr_mature_5_end{$len};
		}
	}
	if ($sums_mature_3_end > 0){
		printf OUTPUT2 "mitotRNAdb-mature-tRNA_3_end-antisense_Match_Genome\t-\t%d\n", $sums_mature_3_end;
		foreach $anno(sort keys %sum_mature_3_end){
			printf OUTPUT2 "mitotRNAdb-mature-tRNA_3_end-antisense_Match_Genome\t$anno\t%.2f\n", $sum_mature_3_end{$anno};
		}
		foreach $len (sort keys %distr_mature_3_end){
			printf OUTPUT3 "mitotRNAdb-mature-tRNA_3_end-antisense_Match_Genome\t$len\t%.2f\n", $distr_mature_3_end{$len};
		}
	}
	if ($sums_mature_CCA_end > 0){
		printf OUTPUT2 "mitotRNAdb-mature-tRNA_CCA_end-antisense_Match_Genome\t-\t%d\n", $sums_mature_CCA_end;
		foreach $anno(sort keys %sum_mature_CCA_end){
			printf OUTPUT2 "mitotRNAdb-mature-tRNA_CCA_end-antisense_Match_Genome\t$anno\t%.2f\n", $sum_mature_CCA_end{$anno};
		}
		foreach $len (sort keys %distr_mature_CCA_end){
			printf OUTPUT3 "mitotRNAdb-mature-tRNA_CCA_end-antisense_Match_Genome\t$len\t%.2f\n", $distr_mature_CCA_end{$len};
		}
	}
}

######summarize and annotate mt_tRNA-antisense seqs: unmatch-genome######
$i += 1;
{
	my %distr_pre;
	my %distr_pre_5_end;
	my %distr_pre_3_end;
	my %distr_mature;
	my %distr_mature_5_end;
	my %distr_mature_3_end;
	my %distr_mature_CCA_end;
	my %sum_pre;
	my %sum_pre_5_end;
	my %sum_pre_3_end;
	my %sum_mature;
	my %sum_mature_5_end;
	my %sum_mature_3_end;
	my %sum_mature_CCA_end;
	my %repeat_num;
	my $sums_pre = 0;
	my $sums_pre_5_end = 0;
	my $sums_pre_3_end = 0;
	my $sums_mature = 0;
	my $sums_mature_5_end = 0;
	my $sums_mature_3_end = 0;	
	my $sums_mature_CCA_end = 0;
	my $len;
	my $anno;
	my $seq;
	my $temp_anno;
	my $tRNA_len;
	my $tRNA_name;
	my $tRNA_codon;
	my $start_site;
	my $temp_end_str;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_pre_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_pre_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_pre_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_pre_3_end += ($reads{$id} / $repeat_num{$id});
			}else {
				$temp_anno = 'pre-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$anno = $temp_anno;
				$sum_pre{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_pre += ($reads{$id} / $repeat_num{$id});
			$distr_pre{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	$i += 1;
	if (-e $out_file[$i] && !-z $out_file[$i]){
		$fh = $file_handle{$i+1};
		while (<$fh>){
			chomp;
			$id = (split /\t/)[0];
			$repeat_num{$id} = 0 unless defined $repeat_num{$id};
			$repeat_num{$id} += 1;
		}
		close $fh;
		open $fh, $out_file[$i];
		while (<$fh>){
			chomp;
			($id, $anno, $start_site, $seq) = (split /\t/)[0, 3, 4, 5];
			$anno =~ / ([0-9]+) bp/;
			$tRNA_len = $1;
			$anno =~ /([A-Za-z]+) \([A-Za-z]+\)/;
			$tRNA_name = $1;
			$anno =~ /[A-Za-z]+ \(([A-Za-z]+)\)/;
			$tRNA_codon = $1;
			$len = length $seq;
			if ($start_site == 0){
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_5_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_mature';
				$sum_mature_5_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_5_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_5_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site + 3)){
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_3_end-antisense';
				$anno = $temp_anno;
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$sum_mature_3_end{$anno} += ($reads{$id} / $repeat_num{$id});
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				$distr_mature_3_end{$len} += ($reads{$id} / $repeat_num{$id});
				$sums_mature_3_end += ($reads{$id} / $repeat_num{$id});
			}elsif ($tRNA_len == ($len + $start_site)){
				$temp_end_str = substr $seq, -3, 3;
				if ($temp_end_str eq 'CCA'){
					$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '_CCA_end-antisense';
					$anno = $temp_anno;
					$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
					$sum_mature_CCA_end{$anno} += ($reads{$id} / $repeat_num{$id});
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
					$distr_mature_CCA_end{$len} += ($reads{$id} / $repeat_num{$id});
					$sums_mature_CCA_end += ($reads{$id} / $repeat_num{$id});
				}else {
					$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
					$anno = $temp_anno;
					$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
				}
			}else {
				$temp_anno = 'mature-mt_tRNA-' . $tRNA_name . '-' . $tRNA_codon . '-antisense';
				$anno = $temp_anno;
				$sum_mature{$temp_anno} += ($reads{$id} / $repeat_num{$id});
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
			$sums_mature += ($reads{$id} / $repeat_num{$id});
			$distr_mature{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
		close $fh;
	}
	if ($sums_pre > 0){
		printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA-antisense_Unmatch_Genome\t-\t%d\n", $sums_pre;
		foreach $anno(sort keys %sum_pre){
			printf OUTPUT2 "mitotRNAdb-pre-mt_tRNA-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_pre{$anno};
		}
		foreach $len (sort keys %distr_pre){
			printf OUTPUT3 "mitotRNAdb-pre-mt_tRNA-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_pre{$len};
		}
	}
	if ($sums_pre_5_end > 0){
		printf OUTPUT2 "mitotRNAdb-pre-tRNA_5_end-antisense_Unmatch_Genome\t-\t%d\n", $sums_pre_5_end;
		foreach $anno(sort keys %sum_pre_5_end){
			printf OUTPUT2 "mitotRNAdb-pre-tRNA_5_end-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_pre_5_end{$anno};
		}
		foreach $len (sort keys %distr_pre_5_end){
			printf OUTPUT3 "mitotRNAdb-pre-tRNA_5_end-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_pre_5_end{$len};
		}
	}
	if ($sums_pre_3_end > 0){
		printf OUTPUT2 "mitotRNAdb-pre-tRNA_3_end-antisense_Unmatch_Genome\t-\t%d\n", $sums_pre_3_end;
		foreach $anno(sort keys %sum_pre_3_end){
			printf OUTPUT2 "mitotRNAdb-pre-tRNA_3_end-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_pre_3_end{$anno};
		}
		foreach $len (sort keys %distr_pre_3_end){
			printf OUTPUT3 "mitotRNAdb-pre-tRNA_3_end-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_pre_3_end{$len};
		}
	}
	if ($sums_mature > 0){
		printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA-antisense_Unmatch_Genome\t-\t%d\n", $sums_mature;
		foreach $anno(sort keys %sum_mature){
			printf OUTPUT2 "mitotRNAdb-mature-mt_tRNA-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature{$anno};
		}
		foreach $len (sort keys %distr_mature){
			printf OUTPUT3 "mitotRNAdb-mature-mt_tRNA-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_mature{$len};
		}
	}
	if ($sums_mature_5_end > 0){
		printf OUTPUT2 "mitotRNAdb-mature-tRNA_5_end-antisense_Unmatch_Genome\t-\t%d\n", $sums_mature_5_end;
		foreach $anno(sort keys %sum_mature_5_end){
			printf OUTPUT2 "mitotRNAdb-mature-tRNA_5_end-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature_5_end{$anno};
		}
		foreach $len (sort keys %distr_mature_5_end){
			printf OUTPUT3 "mitotRNAdb-mature-tRNA_5_end-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_mature_5_end{$len};
		}
	}
	if ($sums_mature_3_end > 0){
		printf OUTPUT2 "mitotRNAdb-mature-tRNA_3_end-antisense_Unmatch_Genome\t-\t%d\n", $sums_mature_3_end;
		foreach $anno(sort keys %sum_mature_3_end){
			printf OUTPUT2 "mitotRNAdb-mature-tRNA_3_end-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature_3_end{$anno};
		}
		foreach $len (sort keys %distr_mature_3_end){
			printf OUTPUT3 "mitotRNAdb-mature-tRNA_3_end-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_mature_3_end{$len};
		}
	}
	if ($sums_mature_CCA_end > 0){
		printf OUTPUT2 "mitotRNAdb-mature-tRNA_CCA_end-antisense_Unmatch_Genome\t-\t%d\n", $sums_mature_CCA_end;
		foreach $anno(sort keys %sum_mature_CCA_end){
			printf OUTPUT2 "mitotRNAdb-mature-tRNA_CCA_end-antisense_Unmatch_Genome\t$anno\t%.2f\n", $sum_mature_CCA_end{$anno};
		}
		foreach $len (sort keys %distr_mature_CCA_end){
			printf OUTPUT3 "mitotRNAdb-mature-tRNA_CCA_end-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr_mature_CCA_end{$len};
		}
	}
}

######summarize and annotate using piRNA-antisense database: match-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		($id, $seq) = (split /\s+/)[0, 5];
		$anno = 'piRNA-antisense';
		$sum{$anno} += ($reads{$id} / $repeat_num{$id});
		$len = length $seq;
		$annos{$id} = $anno;
		$sums += ($reads{$id} / $repeat_num{$id});
		$distr{$len} += ($reads{$id} / $repeat_num{$id});
		if ($unannos_unmatch{$id}){
			delete $unannos_unmatch{$id};
		}
	}
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "piRNAdb-piRNA-antisense_Match_Genome\t-\t%d\n", $sums;
		foreach $len (sort keys %distr){
			printf OUTPUT3 "piRNAdb-piRNA-antisense_Match_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
}

######summarize and annotate using piRNA-antisense database: unmatch-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		($id, $seq) = (split /\s+/)[0, 5];
		$anno = 'piRNA-antisense';
		$sum{$anno} += ($reads{$id} / $repeat_num{$id});
		$len = length $seq;
		$annos{$id} = $anno;
		$sums += ($reads{$id} / $repeat_num{$id});
		$distr{$len} += ($reads{$id} / $repeat_num{$id});
		if ($unannos_unmatch{$id}){
			delete $unannos_unmatch{$id};
		}
	}
	close $fh;
	if ($sums > 0){
		printf OUTPUT2 "piRNAdb-piRNA-antisense_Unmatch_Genome\t-\t%d\n", $sums;
		foreach $len (sort keys %distr){
			printf OUTPUT3 "piRNAdb-piRNA-antisense_Unmatch_Genome\t$len\t%.2f\n", $distr{$len};
		}
	}
}

######summarize and annotate using ensembl-antisense database: match-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*gene_biotype:(\S+?)\s+/){
			$id = $1;
			$anno = $2 . '-antisense';
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
			$sum{$anno} += ($reads{$id} / $repeat_num{$id});
			$len = $lens{$id};
			$distr{$anno}{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
	}
	close $fh;
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

######summarize and annotate using ensembl-antisense database: unmatch-genome######
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*gene_biotype:(\S+?)\s+/){
			$id = $1;
			$anno = $2 . '-antisense';
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
			$sum{$anno} += ($reads{$id} / $repeat_num{$id});
			$len = $lens{$id};
			$distr{$anno}{$len} += ($reads{$id} / $repeat_num{$id});
			if ($unannos_unmatch{$id}){
				delete $unannos_unmatch{$id};
			}
		}
	}
	close $fh;
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
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*(mRNA|cds|protein)/){
			$id = $1;
			$anno = 'mRNA-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(miRNA|microRNA)/){
			$id = $1;
			$anno = 'miRNA-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(tRNA-like)/){
			$id = $1;
			$anno = 'tRNA-like-antisense';			
		}
		elsif (/^(t[0-9]+?)\s.*(tRNA|transfer RNA)/){
			$id = $1;
			$anno = 'tRNA-antisense';		
		}
		elsif (/^(t[0-9]+?)\s.*(rRNA|ribosomal RNA|ribosomal DNA)/){
			$id = $1;
			$anno = 'rRNA-antisense';			
		}
		elsif (/^(t[0-9]+?)\s.*(piRNA)/){
			$id = $1;
			$anno = 'piRNA-antisense';		
		}
		elsif (/^(t[0-9]+?)\s.*(snRNA|small nuclear RNA|small-nuclear RNA)/){
			$id = $1;
			$anno = 'snRNA-antisense';	
		}
		elsif (/^(t[0-9]+?)\s.*(snoRNA)/){
			$id = $1;
			$anno = 'snoRNA-antisense';		
		}
		elsif (/^(t[0-9]+?)\s.*(scRNA)/){
			$id = $1;
			$anno = 'scRNA-antisense';		
		}
		elsif (/^(t[0-9]+?)\s.*(lncRNA|lincRNA)/){
			$id = $1;
			$anno = 'lncRNA-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(antisense)/){
			$id = $1;
			$anno = 'antisense-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(vault RNA)/){
			$id = $1;
			$anno = 'vault_RNA-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(telomerase RNA)/){
			$id = $1;
			$anno = 'telomerase_RNA-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(noncoding)/){
			$id = $1;
			$anno = 'noncoding_RNA-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(RNA)/){
			$id = $1;
			$anno = 'other_RNA-antisense';			
		}
		elsif (/^(t[0-9]+?)\s.*(pseudogene)/){
			$id = $1;
			$anno = 'pseudogene-antisense';		
		}
		elsif (/^(t[0-9]+?)\s.*(gene)/){
			$id = $1;
			$anno = 'gene_region-antisense';		
		}
		elsif (/^(t[0-9]+?)\s.*(mitochondrial|mitochondrion)/){
			$id = $1;
			$anno = 'mt_DNA-antisense';			
		}
		elsif (/^(t[0-9]+?)\s.*(genome|chromosome|BAC|genomic|X-inactivation center)/){
			$id = $1;
			$anno = 'DNA-antisense';		
		}
		elsif (/^(t[0-9]+?)\s/){
			$id = $1;
			$anno = 'other-antisense';	
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
	close $fh;
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
$i += 1;
if (-e $out_file[$i] && !-z $out_file[$i]){
	my %distr;
	my %sum;
	my %repeat_num;
	my $sums = 0;
	my $len;
	my $anno;
	$fh = $file_handle{$i+1};
	while (<$fh>){
		chomp;
		$id = (split /\t/)[0];
		$repeat_num{$id} = 0 unless defined $repeat_num{$id};
		$repeat_num{$id} += 1;
	}
	close $fh;
	open $fh, $out_file[$i];
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*(mRNA|cds|protein)/){
			$id = $1;
			$anno = 'mRNA-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(miRNA|microRNA)/){
			$id = $1;
			$anno = 'miRNA-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(tRNA-like)/){
			$id = $1;
			$anno = 'tRNA-like-antisense';		
		}
		elsif (/^(t[0-9]+?)\s.*(tRNA|transfer RNA)/){
			$id = $1;
			$anno = 'tRNA-antisense';			
		}
		elsif (/^(t[0-9]+?)\s.*(rRNA|ribosomal RNA|ribosomal DNA)/){
			$id = $1;
			$anno = 'rRNA-antisense';			
		}
		elsif (/^(t[0-9]+?)\s.*(piRNA)/){
			$id = $1;
			$anno = 'piRNA-antisense';		
		}
		elsif (/^(t[0-9]+?)\s.*(snRNA|small nuclear RNA|small-nuclear RNA)/){
			$id = $1;
			$anno = 'snRNA-antisense';	
		}
		elsif (/^(t[0-9]+?)\s.*(snoRNA)/){
			$id = $1;
			$anno = 'snoRNA-antisense';		
		}
		elsif (/^(t[0-9]+?)\s.*(scRNA)/){
			$id = $1;
			$anno = 'scRNA-antisense';		
		}
		elsif (/^(t[0-9]+?)\s.*(lncRNA|lincRNA)/){
			$id = $1;
			$anno = 'lncRNA-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(antisense)/){
			$id = $1;
			$anno = 'antisense-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(vault RNA)/){
			$id = $1;
			$anno = 'vault_RNA-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(telomerase RNA)/){
			$id = $1;
			$anno = 'telomerase_RNA-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(noncoding)/){
			$id = $1;
			$anno = 'noncoding_RNA-antisense';
		}
		elsif (/^(t[0-9]+?)\s.*(RNA)/){
			$id = $1;
			$anno = 'other_RNA-antisense';			
		}
		elsif (/^(t[0-9]+?)\s.*(pseudogene)/){
			$id = $1;
			$anno = 'pseudogene-antisense';	
		}
		elsif (/^(t[0-9]+?)\s.*(gene)/){
			$id = $1;
			$anno = 'gene_region-antisense';		
		}
		elsif (/^(t[0-9]+?)\s.*(mitochondrial|mitochondrion)/){
			$id = $1;
			$anno = 'mt_DNA-antisense';		
		}
		elsif (/^(t[0-9]+?)\s.*(genome|chromosome|BAC|genomic|X-inactivation center)/){
			$id = $1;
			$anno = 'DNA-antisense';	
		}
		elsif (/^(t[0-9]+?)\s/){
			$id = $1;
			$anno = 'other-antisense';	
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
	close $fh;
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
	my $sums = 0;
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
	my $sums = 0;
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
        print OUTPUT1 "$id\t$seqs{$id}\t$lens{$id}\t$reads{$id}\tNO\t$annos{$id}\n";
    }else{
        print OUTPUT1 "$id\t$seqs{$id}\t$lens{$id}\t$reads{$id}\tNO\tNO_Annotation\n";
    }
}
