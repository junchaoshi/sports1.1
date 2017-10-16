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
		$name . '_output_tRNA_5_tail_match_genome',		
		$name . '_output_tRNA_3_tail_match_genome',		
		$name . '_output_tRNA_CCA_tail_match_genome',		
		$name . '_output_tRNA_unmatch_genome',
		$name . '_output_tRNA_5_tail_unmatch_genome',
		$name . '_output_tRNA_3_tail_unmatch_genome',
		$name . '_output_tRNA_CCA_tail_unmatch_genome',
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
	open $file_handle{$_}, $out_file[$_-1];
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
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{2};
	while (<$fh>){
		chomp;
		($id, $anno, $seq) = (split /\t/)[0, 3, 5];
		$anno = (split /\s+/, $anno)[0];
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$len = length $seq;
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
		delete $unannos_unmatch{$id};


	}
	print OUTPUT2 "miRBase-miRNA_Match_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t$anno\t$sum{$anno}\n";
	}
	foreach $len (sort keys %distr){
		print OUTPUT3 "miRBase-miRNA_Match_Genome\t$len\t$distr{$len}\n";
	}
}

######summarize and annotate miRNA reads: unmatch-genome######
if (-e $out_file[2] && !-z $out_file[2]){
	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{3};
	while (<$fh>){
		chomp;
		($id, $anno, $seq) = (split /\t/)[0, 3, 5];
		$anno = (split /\s+/, $anno)[0];
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$len = length $seq;
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
		delete $unannos_unmatch{$id};
	}
	print OUTPUT2 "miRBase-miRNA_Unmatch_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t$anno\t$sum{$anno}\n";
	}
	foreach $len (sort keys %distr){
		print OUTPUT3 "miRBase-miRNA_Unmatch_Genome\t$len\t$distr{$len}\n";
	}
}

######summarize and annotate rRNA reads: match-genome######
if (-e $out_file[3] && !-z $out_file[3]){

	my %distr_a;
	my %distr_b;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{4};
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
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$len = length $seq;
		$sums += $reads{$id};
		$distr_a{$len} += $reads{$id};
		$distr_b{$anno}{$len} = 0 unless defined $distr_b{$anno}{$len};
		$distr_b{$anno}{$len} += $reads{$id};
		delete $unannos_unmatch{$id};

	}
	print OUTPUT2 "rRNAdb-rRNA_Match_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t${anno}_Match_Genome\t$sum{$anno}\n";
	}
	foreach $len(sort keys %distr_a){
		print OUTPUT3 "rRNAdb-rRNA_Match_Genome\t$len\t$distr_a{$len}\n";
    }
    foreach $anno(sort keys %distr_b){
	foreach $len(sort keys %{ $distr_b{$anno} }){
	    print OUTPUT3 "rRNAdb-${anno}_Match_Genome\t$len\t$distr_b{$anno}{$len}\n";
	}
    }
}

######summarize and annotate rRNA reads: unmatch-genome######
if (-e $out_file[4] && !-z $out_file[4]){

	my %distr_a;
	my %distr_b;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{5};
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
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$len = length $seq;
		$sums += $reads{$id};
		$distr_a{$len} += $reads{$id};
		$distr_b{$anno}{$len} = 0 unless defined $distr_b{$anno}{$len};
		$distr_b{$anno}{$len} += $reads{$id};
		delete $unannos_unmatch{$id};

	}
	print OUTPUT2 "rRNAdb-rRNA_Unmatch_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t${anno}_Unmatch_Genome\t$sum{$anno}\n";
	}
	foreach $len(sort keys %distr_a){
		print OUTPUT3 "rRNAdb-rRNA_Unmatch_Genome\t$len\t$distr_a{$len}\n";
    }
    foreach $anno(sort keys %distr_b){
	foreach $len(sort keys %{ $distr_b{$anno} }){
	    print OUTPUT3 "rRNAdb-${anno}_Unmatch_Genome\t$len\t$distr_b{$anno}{$len}\n";
	}
    }
	
}

######summarize and annotate tRNA reads: match-genome######
if (-e $out_file[5] && !-z $out_file[5]){

	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	my $temp_anno;
	$fh = $file_handle{6};
	while (<$fh>){
		chomp;
		($id, $anno, $seq) = (split /\t/)[0, 3, 5];
		$anno =~ /\)\s+(.{3,6}?)\s+\((...)\)\s+\d+/;
		$temp_anno = 'tRNA-' . $1 . '-' . $2;
		$anno = $temp_anno;
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$len = length $seq;
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
		delete $unannos_unmatch{$id};
	}
	print OUTPUT2 "tRNAdb-tRNA_Match_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t$anno\t$sum{$anno}\n";
	}
	foreach $len (sort keys %distr){
		print OUTPUT3 "tRNAdb-tRNA_Match_Genome\t$len\t$distr{$len}\n";
	}
}

######summarize and annotate tRNA-5-end reads: match-genome######
if (-e $out_file[6] && !-z $out_file[6]){
	
	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{7};
	while (<$fh>){
		chomp;
		($id, $seq, $len, $anno) = split /\t/;
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
		delete $unannos_unmatch{$id};
	}
	print OUTPUT2 "tRNAdb-tRNA_5_end_Match_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t$anno\t$sum{$anno}\n";
	}
	foreach $len (sort keys %distr){
		print OUTPUT3 "tRNAdb-tRNA_5_end_Match_Genome\t$len\t$distr{$len}\n";
	}
}

######summarize and annotate tRNA-3-end reads: match-genome######
if (-e $out_file[7] && !-z $out_file[7]){
	
	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{8};
	while (<$fh>){
		chomp;
		($id, $seq, $len, $anno) = split /\t/;
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
		delete $unannos_unmatch{$id};
	}
	print OUTPUT2 "tRNAdb-tRNA_3_end_Match_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t$anno\t$sum{$anno}\n";
	}
	foreach $len (sort keys %distr){
		print OUTPUT3 "tRNAdb-tRNA_3_end_Match_Genome\t$len\t$distr{$len}\n";
	}
}

######summarize and annotate tRNA-CCA-end reads: match-genome######
if (-e $out_file[8] && !-z $out_file[8]){
	
	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{9};
	while (<$fh>){
		chomp;
		($id, $seq, $len, $anno) = split /\t/;
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
		delete $unannos_unmatch{$id};
	}
	print OUTPUT2 "tRNAdb-tRNA_CCA_end_Match_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t$anno\t$sum{$anno}\n";
	}
	foreach $len (sort keys %distr){
		print OUTPUT3 "tRNAdb-tRNA_CCA_end_Match_Genome\t$len\t$distr{$len}\n";
	}
}

######summarize and annotate tRNA reads: unmatch-genome######
if (-e $out_file[9] && !-z $out_file[9]){

	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	my $temp_anno;
	$fh = $file_handle{10};
	while (<$fh>){
		chomp;
		($id, $anno, $seq) = (split /\t/)[0, 3, 5];
		$anno =~ /\)\s+(.{3,6}?)\s+\((...)\)\s+\d+/;
		$temp_anno = 'tRNA-' . $1 . '-' . $2;
		$anno = $temp_anno;
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$len = length $seq;
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
		delete $unannos_unmatch{$id};
	}
	print OUTPUT2 "GtRNAdb-tRNA_Unmatch_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t$anno\t$sum{$anno}\n";
	}
	foreach $len (sort keys %distr){
		print OUTPUT3 "GtRNAdb-tRNA_Unmatch_Genome\t$len\t$distr{$len}\n";
	}
}

######summarize and annotate tRNA-5-end reads unmatch-genome######
if (-e $out_file[10] && !-z $out_file[10]){
	
	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{11};
	while (<$fh>){
		chomp;
		($id, $seq, $len, $anno) = split /\t/;
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
		delete $unannos_unmatch{$id};
	}
	print OUTPUT2 "tRNAdb-tRNA_5_end_Unmatch_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t$anno\t$sum{$anno}\n";
	}
	foreach $len (sort keys %distr){
		print OUTPUT3 "GtRNAdb-tRNA_5_end_Unmatch_Genome\t$len\t$distr{$len}\n";
	}
}

######summarize and annotate tRNA-3-end reads unmatch-genome######
if (-e $out_file[11] && !-z $out_file[11]){
	
	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{12};
	while (<$fh>){
		chomp;
		($id, $seq, $len, $anno) = split /\t/;
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
		delete $unannos_unmatch{$id};

	}
	print OUTPUT2 "tRNAdb-tRNA_3_end_Unmatch_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t$anno\t$sum{$anno}\n";
	}
	foreach $len (sort keys %distr){
		print OUTPUT3 "GtRNAdb-tRNA_3_end_Unmatch_Genome\t$len\t$distr{$len}\n";
	}
}

######summarize and annotate tRNA-CCA-end reads unmatch-genome######
if (-e $out_file[12] && !-z $out_file[12]){
	
	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{13};
	while (<$fh>){
		chomp;
		($id, $seq, $len, $anno) = split /\t/;
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
		delete $unannos_unmatch{$id};
	}
	print OUTPUT2 "tRNAdb-tRNA_CCA_end_Unmatch_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t$anno\t$sum{$anno}\n";
	}
	foreach $len (sort keys %distr){
		print OUTPUT3 "tRNAdb-tRNA_CCA_end_Unmatch_Genome\t$len\t$distr{$len}\n";
	}
}

######summarize and annotate using piRNA database: match-genome######
if (-e $out_file[13] && !-z $out_file[13]){
	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{14};
	while (<$fh>){
		chomp;
		($id, $anno, $seq) = (split /\s+/)[0, 3, 5];
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$len = length $seq;
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
		delete $unannos_unmatch{$id};
	}
	print OUTPUT2 "piRNAdb-piRNA_Match_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t$anno\t$sum{$anno}\n";
	}
	foreach $len (sort keys %distr){
		print OUTPUT3 "piRNAdb-piRNA_Match_Genome\t$len\t$distr{$len}\n";
	}
}

######summarize and annotate using piRNA database: match-genome######
if (-e $out_file[14] && !-z $out_file[14]){
	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	my $seq;
	$fh = $file_handle{15};
	while (<$fh>){
		chomp;
		($id, $anno, $seq) = (split /\s+/)[0, 3, 5];
		$annos{$id} = $anno;
		$sum{$anno} += $reads{$id};
		$len = length $seq;
		$sums += $reads{$id};
		$distr{$len} += $reads{$id};
		delete $unannos_unmatch{$id};
	}
	print OUTPUT2 "piRNAdb-piRNA_Unmatch_Genome\t-\t$sums\n";
	foreach $anno(sort keys %sum){
		print OUTPUT2 "-\t$anno\t$sum{$anno}\n";
	}
	foreach $len (sort keys %distr){
		print OUTPUT3 "piRNAdb-piRNA_Unmatch_Genome\t$len\t$distr{$len}\n";
	}
}

######summarize and annotate using ensembl database: match-genome######
if (-e $out_file[15] && !-z $out_file[15]){

	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	$fh = $file_handle{16};
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*gene_biotype:(\S+?)\s+/){
			$id = $1;
			$anno = $2;
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
			delete $unannos_unmatch{$id};
		}
	}
	foreach $anno(sort keys %sum){
		print OUTPUT2 "ensembl-${anno}_Match_Genome\t-\t$sum{$anno}\n";
	}
	foreach	my $key1(sort keys %distr){
		foreach my $key2 (sort keys %{$distr{$key1}}){
			print OUTPUT3 "ensembl-${key1}_Match_Genome\t$key2\t$distr{$key1}{$key2}\n";
		}
	}
}

######summarize and annotate using ensembl database: unmatch-genome######
if (-e $out_file[16] && !-z $out_file[16]){

	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	$fh = $file_handle{17};
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*gene_biotype:(\S+?)\s+/){
			$id = $1;
			$anno = $2;
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
			delete $unannos_unmatch{$id};
		}
	}
	foreach $anno(sort keys %sum){
		print OUTPUT2 "ensembl-${anno}_Unmatch_Genome\t-\t$sum{$anno}\n";
	}
	foreach	my $key1(sort keys %distr){
		foreach my $key2 (sort keys %{$distr{$key1}}){
			print OUTPUT3 "ensembl-${key1}_Unmatch_Genome\t$key2\t$distr{$key1}{$key2}\n";
		}
	}
}

######summarize and annotate using rfam database: match-genome######
if (-e $out_file[17] && !-z $out_file[17]){
	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	$fh = $file_handle{18};
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*(mRNA|cds|protein)/){
			$id = $1;
			$anno = 'mRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(miRNA|microRNA)/){
			$id = $1;
			$anno = 'miRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(tRNA-like)/){
			$id = $1;
			$anno = 'tRNA-like';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(tRNA|transfer RNA)/){
			$id = $1;
			$anno = 'tRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(rRNA|ribosomal RNA|ribosomal DNA)/){
			$id = $1;
			$anno = 'rRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(piRNA)/){
			$id = $1;
			$anno = 'piRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(snRNA|small nuclear RNA|small-nuclear RNA)/){
			$id = $1;
			$anno = 'snRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};		
		}
		elsif (/^(t[0-9]+?)\s.*(snoRNA)/){
			$id = $1;
			$anno = 'snoRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(scRNA)/){
			$id = $1;
			$anno = 'scRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};		
		}
		elsif (/^(t[0-9]+?)\s.*(lncRNA|lincRNA)/){
			$id = $1;
			$anno = 'lncRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(antisense)/){
			$id = $1;
			$anno = 'antisense';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(vault RNA)/){
			$id = $1;
			$anno = 'vault_RNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(telomerase RNA)/){
			$id = $1;
			$anno = 'telomerase_RNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(noncoding)/){
			$id = $1;
			$anno = 'noncoding_RNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(RNA)/){
			$id = $1;
			$anno = 'other_RNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(pseudogene)/){
			$id = $1;
			$anno = 'pseudogene';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};		
		}
		elsif (/^(t[0-9]+?)\s.*(gene)/){
			$id = $1;
			$anno = 'gene_region';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(mitochondrial|mitochondrion)/){
			$id = $1;
			$anno = 'mt_DNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(genome|chromosome|BAC|genomic|X-inactivation center)/){
			$id = $1;
			$anno = 'DNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};		
		}
		elsif (/^(t[0-9]+?)\s/){
			$id = $1;
			$anno = 'other';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};		
		}
		delete $unannos_unmatch{$id};
	}
	foreach $anno(sort keys %sum){
		print OUTPUT2 "Rfam-${anno}_Match_Genome\t-\t$sum{$anno}\n";
	}
	foreach	my $key1(sort keys %distr){
		foreach my $key2 (sort keys %{$distr{$key1}}){
			print OUTPUT3 "Rfam-${key1}_Match_Genome\t$key2\t$distr{$key1}{$key2}\n";
		}
	}	
}


######summarize and annotate using rfam database: unmatch-genome######
if (-e $out_file[18] && !-z $out_file[18]){
	my %distr;
	my %sum;
	my $sums  = 0;
	my $len;
	my $anno;
	$fh = $file_handle{19};
	while (<$fh>){
		chomp;
		if (/^(t[0-9]+?)\s.*(mRNA|cds|protein)/){
			$id = $1;
			$anno = 'mRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(miRNA|microRNA)/){
			$id = $1;
			$anno = 'miRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(tRNA-like)/){
			$id = $1;
			$anno = 'tRNA-like';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(tRNA|transfer RNA)/){
			$id = $1;
			$anno = 'tRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(rRNA|ribosomal RNA|ribosomal DNA)/){
			$id = $1;
			$anno = 'rRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(piRNA)/){
			$id = $1;
			$anno = 'piRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(snRNA|small nuclear RNA|small-nuclear RNA)/){
			$id = $1;
			$anno = 'snRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};		
		}
		elsif (/^(t[0-9]+?)\s.*(snoRNA)/){
			$id = $1;
			$anno = 'snoRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(scRNA)/){
			$id = $1;
			$anno = 'scRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};		
		}
		elsif (/^(t[0-9]+?)\s.*(lncRNA|lincRNA)/){
			$id = $1;
			$anno = 'lncRNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(antisense)/){
			$id = $1;
			$anno = 'antisense';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(vault RNA)/){
			$id = $1;
			$anno = 'vault_RNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(telomerase RNA)/){
			$id = $1;
			$anno = 'telomerase_RNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(noncoding)/){
			$id = $1;
			$anno = 'noncoding_RNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};
		}
		elsif (/^(t[0-9]+?)\s.*(RNA)/){
			$id = $1;
			$anno = 'other_RNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(pseudogene)/){
			$id = $1;
			$anno = 'pseudogene';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};		
		}
		elsif (/^(t[0-9]+?)\s.*(gene)/){
			$id = $1;
			$anno = 'gene_region';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(mitochondrial|mitochondrion)/){
			$id = $1;
			$anno = 'mt_DNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};			
		}
		elsif (/^(t[0-9]+?)\s.*(genome|chromosome|BAC|genomic|X-inactivation center)/){
			$id = $1;
			$anno = 'DNA';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};		
		}
		elsif (/^(t[0-9]+?)\s/){
			$id = $1;
			$anno = 'other';
			$annos{$id} = $anno;
			$sum{$anno} += $reads{$id};
			$len = $lens{$id};
			$distr{$anno}{$len} += $reads{$id};		
		}
		delete $unannos_unmatch{$id};
	}
	foreach $anno(sort keys %sum){
		print OUTPUT2 "Rfam-${anno}_Unmatch_Genome\t-\t$sum{$anno}\n";
	}
	foreach	my $key1(sort keys %distr){
		foreach my $key2 (sort keys %{$distr{$key1}}){
			print OUTPUT3 "Rfam-${key1}_Unmatch_Genome\t$key2\t$distr{$key1}{$key2}\n";
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
        print OUTPUT1 "$id\t$seqs{$id}\t$lens{$id}\t$reads{$id}\tNO\t$annos{$id}\n";
    }else{
        print OUTPUT1 "$id\t$seqs{$id}\t$lens{$id}\t$reads{$id}\tNO\tNO_Annotation\n";
    }
}

