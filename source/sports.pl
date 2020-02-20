#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Std;
use vars qw($opt_i $opt_a $opt_b $opt_x $opt_y $opt_l $opt_L $opt_M $opt_s $opt_p $opt_g $opt_m $opt_r $opt_t $opt_e $opt_f $opt_w $opt_o $opt_k $opt_z $opt_v $opt_h);
getopts('i:abx:y:l:L:M:sp:g:m:r:t:e:f:w:o:kzvh');
my $input_file		= $opt_i;
my $opt_adapter		= $opt_a ? 1 : 0;
my $opt_format		= $opt_b ? 1 : 0;
my $adapter_5		= $opt_x ? $opt_x : "default";
my $adapter_3		= $opt_y ? $opt_y : "default";
my $min_length		= $opt_l ? $opt_l : 15;
my $max_length		= $opt_L ? $opt_L : 45;
my $mismatch		= $opt_M ? $opt_M : 0;
my $sense			= $opt_s ? 1 : 0;
my $thread		= $opt_p ? $opt_p : 1;
my $genome_address	= $opt_g;
my $miRNA_db_address	= $opt_m ? $opt_m : "NULL";
my $rRNA_db_address	= $opt_r ? $opt_r : "NULL";
my $tRNA_db_address	= $opt_t ? $opt_t : "NULL";
my $ensembl_nc_address	= $opt_e ? $opt_e : "NULL";
my $rfam_address	= $opt_f ? $opt_f : "NULL";
my $piRNA_db_address	= $opt_w ? $opt_w : "NULL";
my $output_address	= $opt_o;
my $keep_all		= $opt_k;
my $mismatch_stat	= $opt_z;
my $version		= $opt_v ? 1 : 0;
my $help		= $opt_h ? 1 : 0;

my $threshold = 10;
my $seq_err = 0.01;

my $version_info = "1.1.0";

my $usage = <<"USAGE";
Description:	Perl script used to annotate small RNA sequences in batch.
		Query input files support Solexa / Illumina sequence formats: .sra, .fastq/.fq, .fasta/.fa.
		Attention: compressed files need to be unpacked before input! 

Usage Examples:	sports.pl -i reads.fa -g /foo/bar/Genome/Mouse/UCSC/mm10/Sequence/BowtieIndex/Genome
		or: 
		sports.pl -i seq_address.txt -p 4 -g /foo/bar/Genome/Human/UCSC/hg38/Sequence/BowtieIndex/Genome -m /foo/bar/Database/Human/miRBase_21/miRbase_21 -r /foo/bar/Database/Human/rRNA_db/human_rRNA -t /foo/bar/Database/Human/GtRNAdb/hg19-tRNAs -w /foo/bar/Database/Human/piRBase/piR_human_v1.0 -o /foo/bar/output/
		or:
		sports.pl -i /foo/bar/download_seq/ -p 4 -g /foo/bar/Genome/rat/UCSC/rn6/Sequence/BowtieIndex/Genome -m /foo/bar/Database/Rat/miRBase_21/miRbase_21 -r /foo/bar/Database/Rat/rRNA_db/rat_rRNA -t /foo/bar/Database/Rat/GtRNAdb/rn5-tRNAs -w /foo/bar/Database/Rat/piRBase/piR_rat_v1.0 -e /foo/bar/Database/Rat/Ensembl/release-89/Rattus_norvegicus.Rnor_6.0.ncrna -f /foo/bar/Database/Rat/Rfam_12.3/Rfam-12.3-rat -o /foo/bar/output/ -k

Input:
  -i <file>	Input could be: 
		a directory (will run all qualified files in the directory recursively); 
		a .txt file (for batch processing data, which should contain absolute path of input files or directories); 
		a .fastq/.fq or .fasta/.fa file. 
  -p <int>	number of threads to launch (default = 1)


Index:
  -g <str>	reference genome bowtie index
  -m <str>	miRNA database bowtie index (optional)
  -r <str>	rRNA database bowtie index (optional)
  -t <str>	tRNA database bowtie index (optional)
  -e <str>	ensembl noncoding RNA database bowtie index (optional)
  -f <str>	rfam database bowtie index (optional)
  -w <str>	piRNA database bowtie index (optional)

Output:
  -o <str>	output address of annotation results (default: input address)
  -k		keep all the intermediate files generated during the running progress
  -z		output mismatch statistics (optional, -M should larger than 0)
  
Alignment:
  -l <int>	the minimal length of the output sequences (default = 15)
  -L <int>	the maximal length of the output sequences (default = 45)
  -M <int>	the total number of mismatches in the entire alignment (default = 0)
  -s 		align to forward/reverse-complement reference strand of annotation database
  -a		Remove 5\'/3\' adapters
  -x <str>	(When -a is selected) Your 5\' adapter sequence or use "default" adapter. Default = "GTTCAGAGTTCTACAGTCCGACGATC"
  -y <str>	(When -a is selected) Your 3\' adapter sequence or use "default" adapter. Default = "TGGAATTCTCGGGTGCCAAGG"

Others:
  -v		print version information
  -h		print this usage message
USAGE


if ($version) {
	print "sports version : $version_info\n";
	exit;
}elsif($help) {
	print $usage;
	exit;
}

##determine if input file and genome file are defined
unless (defined $input_file && 
	defined $genome_address
	){
	print "Input file genome bowtie index should be specified!\n\n";
	print $usage;
	exit;
}

unless (-e $input_file){
	print "Input file is not exist!\n\n";
	print $usage;
	exit;
}

##determine if adapter sequences are valid
if($opt_adapter){
	unless ($adapter_5 =~ /^[atgcu]+$/i || $adapter_5 =~ /^default$/i ){
	print "Invalid 5\' end adapter input!\n\n";
	print $usage;
	exit;
	}
	unless ($adapter_3 =~ /^[atgcu]+$/i || $adapter_3 =~ /^default$/i ){
	print "Invalid 3\' end adapter input!\n\n";
	print $usage;
	exit;
	}
}

##define parameters
my @input;
my $input_address = abs_path($input_file);

my @dirlist;
my @filelist;
my @tmp_filelist;

my $input_query_address;
my $input_query_name;
my $input_query_suffix;
my $strand;

my $script_address = `which sports.pl`;
   @input = split(/\//, $script_address);
   pop @input;
   $script_address = join('/', @input) . '/';

##generate genome bowtie index
my $genome_bowtie_file = $genome_address . ".1.ebwt";
unless (-e $genome_bowtie_file){
	print "\n\nGenerating genome bowtie index...\n\n";
	system ("bowtie-build -q ${genome_address}.fa ${genome_address}");
}

##generate miRNA database bowtie index
my $miRNA_db_bowtie_file = $miRNA_db_address . ".1.ebwt";
unless($miRNA_db_address eq "NULL"){
	unless (-e $miRNA_db_bowtie_file){
		print "\n\nGenerating miRNA database bowtie index...\n\n";
		system ("bowtie-build -q ${miRNA_db_address}.fa ${miRNA_db_address}");
	}
}

##generate ensembl noncoding RNA database bowtie index
my $ensembl_nc_bowtie_file = $ensembl_nc_address . ".1.ebwt";
unless($ensembl_nc_address eq "NULL"){
	unless (-e $ensembl_nc_bowtie_file){
		print "\n\nGenerating ensembl database bowtie index...\n\n";
		system ("bowtie-build -q ${ensembl_nc_address}.fa ${ensembl_nc_address}");
	}
}

##generate rfam database bowtie index
my $rfam_bowtie_file = $rfam_address . ".1.ebwt";
unless($rfam_address eq "NULL"){
	unless (-e $rfam_bowtie_file){
		print "\n\nGenerating rfam database bowtie index...\n\n";
		system ("bowtie-build -q ${rfam_address}.fa ${rfam_address}");
	}
}

##generate piRNA database bowtie index
my $piRNA_db_bowtie_file = $piRNA_db_address . ".1.ebwt";
unless($piRNA_db_address eq "NULL"){
	unless (-e $piRNA_db_bowtie_file){
		print "\n\nGenerating piRNA database bowtie index...\n\n";
		system ("bowtie-build -q ${piRNA_db_address}.fa ${piRNA_db_address}");
	}
}

##generate mature tRNA database and bowtie index
my $tRNA_db_mature_tRNA_file = $tRNA_db_address;
$tRNA_db_mature_tRNA_file =~ s/-tRNAs$//;
$tRNA_db_mature_tRNA_file = $tRNA_db_mature_tRNA_file . "-mature-tRNAs.fa";

my $tRNA_db_mito_tRNA_file = $tRNA_db_address;
$tRNA_db_mito_tRNA_file =~ s/-tRNAs$//;
$tRNA_db_mito_tRNA_file = $tRNA_db_mito_tRNA_file . "-mt_tRNAs";

unless($tRNA_db_address eq "NULL"){
	##genomic-tRNA-bowtie-build
	my $tRNA_db_tRNA_pre_file = $tRNA_db_address. ".1.ebwt";
	unless (-e $tRNA_db_tRNA_pre_file){
		print "\n\nGenerating pre-tRNA database bowtie index...\n\n";
		system ("bowtie-build -q ${tRNA_db_address}.fa ${tRNA_db_address}");
	}
	my $tRNA_db_tRNA_CCA_file = $tRNA_db_address;
	$tRNA_db_tRNA_CCA_file = $tRNA_db_address . "_CCA.1.ebwt";
	unless (-e $tRNA_db_tRNA_CCA_file){
		if (-e $tRNA_db_mature_tRNA_file){
			system ("perl ${script_address}tRNA_db_processing.pl ${tRNA_db_mature_tRNA_file}");
			$tRNA_db_mature_tRNA_file =~ s/\.fa$//;
			system ("mv ${tRNA_db_mature_tRNA_file}_CCA.fa ${tRNA_db_address}_CCA.fa");
			$tRNA_db_mature_tRNA_file = $tRNA_db_mature_tRNA_file . ".fa";
		}else{
			system ("perl ${script_address}tRNA_db_processing.pl ${tRNA_db_address}.fa");
		}
		print "\n\nGenerating mature-tRNA database bowtie index...\n\n";
		system ("bowtie-build -q ${tRNA_db_address}_CCA.fa ${tRNA_db_address}_CCA");
	}
	##mito-tRNA-bowtie-build
	my $tRNA_db_mito_tRNA_bowtie_file = $tRNA_db_mito_tRNA_file . ".1.ebwt";
	unless (-e $tRNA_db_mito_tRNA_bowtie_file){
		print "\n\nGenerating mito-tRNA database bowtie index...\n\n";
		system ("bowtie-build -q ${tRNA_db_mito_tRNA_file}.fa ${tRNA_db_mito_tRNA_file}");
		system ("bowtie-build -q ${tRNA_db_mito_tRNA_file}_CCA.fa ${tRNA_db_mito_tRNA_file}_CCA");
	}
}

##push all the query into an array: @filelist
if (-d $input_address){
	@input = split(/\//, $input_address);
	$input_address = join('/', @input) . '/';
	@dirlist = ($input_address);
	print "\n\nSearching input files...\n\n";
	while (@dirlist){
		my $tmp_d = $dirlist[0];
		opendir DIR, $tmp_d || die "Cannot open directory: $tmp_d !\n";
		my @files = readdir DIR;
		closedir DIR;
		my $tmp_f;
		foreach (@files){
			$tmp_f = $tmp_d . $_;
			if ($_ eq "." || $_ eq ".."){
				next;
			}
			if (-d $tmp_f){
				$tmp_f = $tmp_d . $_ . '/';
				push (@dirlist, $tmp_f);
			}
			if (-f $tmp_f){
				@input = split(/\//, $tmp_f);
				$input_query_name = pop (@input);
				$input_query_address = join('/', @input) . '/';
				@input = split(/\./, $input_query_name);
				$input_query_suffix = pop (@input);
				$input_query_name = join('.', @input);
				if ($input_query_suffix eq 'fastq' || $input_query_suffix eq 'fq' || 
				    $input_query_suffix eq 'fasta' || $input_query_suffix eq 'fa' ||
		                    $input_query_suffix eq 'sra'){
					push (@filelist, $tmp_f);
					print "$input_query_address$input_query_name.$input_query_suffix\n";
				}
			}
		}
		shift @dirlist;
	}
	@tmp_filelist = sort { lc($a) cmp lc($b) } @filelist;
	@filelist = @tmp_filelist;
}
elsif (-f $input_address){
	@input = split(/\//, $input_address);
	$input_query_name = pop (@input);
	$input_query_address = join('/', @input) . '/';
	@input = split(/\./, $input_query_name);
	$input_query_suffix = pop (@input);
	$input_query_name = join('.', @input);
	my $tmp_f = $input_query_address . $input_query_name . '.' . $input_query_suffix;
	if ($input_query_suffix eq 'fastq' || $input_query_suffix eq 'fq' || 
	    $input_query_suffix eq 'fasta' || $input_query_suffix eq 'fa' ||
	    $input_query_suffix eq 'sra'){
		push (@filelist, $tmp_f);
	}
	elsif ($input_query_suffix eq 'txt'){
		open FILE_IN, $input_address || die "Cannot open file ${input_query_name}.${input_query_suffix} !\n";
		while (<FILE_IN>){
			if ($_ =~ /^\s+$/){
				next;
			}
			chomp;
			my $tmp = abs_path($_);
			if (-d $tmp){
				@input = split(/\//, $tmp);
				$tmp = join('/', @input) . '/';
				@dirlist = ($tmp);
				while (@dirlist){
					my $tmp_d = $dirlist[0];
					opendir DIR, $tmp_d || die "Cannot open directory: $tmp_d !\n";
					my @files = readdir DIR;
					closedir DIR;
					my $tmp_f;
					foreach (@files) {
						$tmp_f = $tmp_d . $_;
						if ($_ eq "." || $_ eq ".."){
							next;
						}
						if (-d $tmp_f){
							push (@dirlist, $tmp_f);
						}
						if (-f $tmp_f){
							@input = split(/\//, $tmp_f);
							$input_query_name = pop (@input);
							$input_query_address = join('/', @input) . '/';
							@input = split(/\./, $input_query_name);
							$input_query_suffix = pop (@input);
							$input_query_name = join('.', @input);
							if ($input_query_suffix eq 'fastq' || $input_query_suffix eq 'fq' || 
				  			    $input_query_suffix eq 'fasta' || $input_query_suffix eq 'fa' ||
		                   			    $input_query_suffix eq 'sra'){
								push (@filelist, $tmp_f);
							}
						}
					}
				}
				shift @dirlist;
			}
			elsif (-f $tmp){
				@input = split(/\//, $tmp);
				$input_query_name = pop (@input);
				$input_query_address = join('/', @input) . '/';
				@input = split(/\./, $input_query_name);
				$input_query_suffix = pop (@input);
				$input_query_name = join('.', @input);
				if ($input_query_suffix eq 'fastq' || $input_query_suffix eq 'fq' || 
				    $input_query_suffix eq 'fasta' || $input_query_suffix eq 'fa' ||
		                    $input_query_suffix eq 'sra'){
					push (@filelist, $tmp);
				}
			}
	
		}
	}
	@tmp_filelist = @filelist;
}
else{
	print "Input is not valid!\n";
	print $usage;
	exit;
}


if (-f $input_address){
	@input = split(/\//, $input_address);
	pop (@input);
	$input_address = join('/', @input) . '/';
}


unless (defined $output_address){
	$output_address = $input_address;
}

unless (-e $output_address){
	`mkdir -p $output_address`;
}
$output_address = abs_path($output_address);

{
@input = split(/\//, $output_address);
$output_address = join('/', @input) . '/';
}

my $output_file = $output_address . 'run_this.sh';

open FILE_OUT, ">$output_file" || die $!;

my $count = 0;
my @rRNA_length;

unless ($rRNA_db_address eq "NULL"){
	my @rRNA = ('2S', '4.5S', '5S', '5.3S', '5.8S', '12S', '16S', '17S', '18S', '25S', '26S', '28S', '45S', 'RNY1', 'RNY3', 'RNY4', 'RNY5', 'other');
	my $bowtie_fa;
	my $bowtie_index;
	my $temp_length;
	foreach (@rRNA){
		$bowtie_fa = "${rRNA_db_address}_${_}.fa";
		$bowtie_index = "${rRNA_db_address}_${_}.1.ebwt";
		if (-e $bowtie_fa){
			unless (-e $bowtie_index){
				print "\n\nGenerating ${_} rRNA database bowtie index...\n\n";
				system ("bowtie-build -q ${rRNA_db_address}_${_}.fa ${rRNA_db_address}_${_}");
			}
			$temp_length = `grep -v '>' $bowtie_fa | wc -m` - `grep -v '>' $bowtie_fa | wc -w`;
			push (@rRNA_length, "${_}=${temp_length}");
		}
	}
}

while (@tmp_filelist){
	my $tmp_f = $tmp_filelist[0];
	$count += 1;
	@input = split(/\//, $tmp_f);
	$input_query_name = pop (@input);
	$input_query_address = join('/', @input) . '/';
	@input = split(/\./, $input_query_name);
	$input_query_suffix = pop (@input);
	$input_query_name = join('.', @input);
	
	my $tmp_sh = $output_address . "sh_file";
	unless (-e $tmp_sh){
	`mkdir -p $tmp_sh`;
	}
	my $tmp_repo = $output_address . "processing_report/";
	unless (-e $tmp_repo){
	`mkdir -p $tmp_repo`;
	}
	print FILE_OUT "echo \"$input_query_address$input_query_name.$input_query_suffix\"\n";
	print FILE_OUT "sh \"${tmp_sh}/${count}_$input_query_name.sh\" > \"$tmp_repo${count}_${input_query_name}.txt\" 2>&1\n";
	open FILE, '>',"${tmp_sh}/${count}_$input_query_name.sh"
		or die "can not open '${tmp_sh}/${count}_$input_query_name.sh'";
	print FILE '#!/bin/bash	
##This script uses to annotate small RNA
date
echo ""
thread="' . $thread . '"
adapter5="' . $adapter_5 . '"
adapter3="' . $adapter_3 . '"
min_length="' . $min_length . '"
max_length="' . $max_length . '"
mismatch="' . $mismatch . '"
input_query_address="' . $input_query_address . '"
input_query_name="' . $input_query_name . '"
input_query_suffix="' . $input_query_suffix . '"
output_address="' . $output_address . ${count} . '_${input_query_name}/"
script_address="' . $script_address . '"


if [ ! -d "${output_address}" ]; then
	mkdir -p ${output_address}
fi';

	if($opt_format){
		print FILE '
cp ${input_query_address}${input_query_name}.${input_query_suffix} ${output_address}
cd ${output_address}';
	}else{
		print FILE '
ln -s ${input_query_address}${input_query_name}.${input_query_suffix} ${output_address}
cd ${output_address}';
	}

###input query format is .sra
	if ($input_query_suffix eq "sra"){
		print FILE '
fastq-dump ${output_address}${input_query_name}.sra -O ${output_address}
input_query_suffix=fastq
rm -rf ${output_address}${input_query_name}.sra
';
		$input_query_suffix = 'fastq';

	}

###remove adapter
	if ($opt_adapter){
		if ($adapter_5 =~ /^default$/i && 
		    $adapter_3 =~ /^default$/i){
my $adapter_5		= "GTTCAGAGTTCTACAGTCCGACGATC";
my $adapter_3		= "TGGAATTCTCGGGTGCCAAGG";
		print FILE '
echo ""
echo "remove default 5\' and 3\' adapters"
cutadapt -j ${thread} -g ' . $adapter_5 . ' -a ' . $adapter_3 . ' -o ${output_address}${input_query_name}_trim_1.${input_query_suffix} --max-n 0 ${output_address}${input_query_name}.${input_query_suffix}
rm -rf ${output_address}${input_query_name}.${input_query_suffix}
		';
		}
		elsif ($adapter_5 !~ /^default$/i && 
		       $adapter_3 =~ /^default$/i){
		print FILE '
echo ""
echo "remove 5\' adapter"
cutadapt -j ${thread} -g ' . $adapter_5 . ' -o ${output_address}${input_query_name}_trim_1.${input_query_suffix} --max-n 0 ${output_address}${input_query_name}.${input_query_suffix}
rm -rf ${output_address}${input_query_name}.${input_query_suffix}
		';
		}
		elsif ($adapter_5 =~ /^default$/i && 
		       $adapter_3 !~ /^default$/i){
		print FILE '
echo ""
echo "remove 3\' adapter"
cutadapt -j ${thread} -a ' . $adapter_3 . ' -o ${output_address}${input_query_name}_trim_1.${input_query_suffix} --max-n 0 ${output_address}${input_query_name}.${input_query_suffix}
rm -rf ${output_address}${input_query_name}.${input_query_suffix}
		';
		}
		elsif ($adapter_5 !~ /^default$/i && 
		       $adapter_3 !~ /^default$/i){
		print FILE '
echo ""
echo "remove 5-end and 3-end adapters"
cutadapt -j ${thread} -g ' . $adapter_5 . ' -a ' . $adapter_3 . ' -o ${output_address}${input_query_name}_trim_1.${input_query_suffix} --max-n 0 ${output_address}${input_query_name}.${input_query_suffix}
		';
		}
	}
	else{
		print FILE '
mv ${output_address}${input_query_name}.${input_query_suffix} ${output_address}${input_query_name}_trim_1.${input_query_suffix}
		';
	}

###organize sequence
	unless($opt_format){
		if($input_query_suffix eq 'fastq' || $input_query_suffix eq 'fq' ){
			print FILE '
perl ${script_address}fastq2fasta.pl ${output_address}${input_query_name}_trim_1.${input_query_suffix} > ${output_address}${input_query_name}_trim_1.fa
rm ${output_address}${input_query_name}_trim_1.${input_query_suffix}
input_query_suffix=fa';
			$input_query_suffix = 'fa';
		}

		print FILE '

perl ${script_address}fastaparse.pl ${output_address}${input_query_name}_trim_1.fa -b > ${output_address}${input_query_name}_trim_2.fa 2>${output_address}${input_query_name}_discarded_reads.fa
rm ${output_address}${input_query_name}_trim_1.fa

perl ${script_address}fastaparse.pl ${output_address}${input_query_name}_trim_2.fa -a ${min_length} > ${output_address}${input_query_name}_trim_3.fa 2>${output_address}${input_query_name}_too_short_reads.fa
rm ${output_address}${input_query_name}_trim_2.fa

perl ${script_address}fastaparse.pl ${output_address}${input_query_name}_trim_3.fa -c ${max_length} > ${output_address}${input_query_name}_trim_4.fa 2>${output_address}${input_query_name}_too_long_reads.fa
rm ${output_address}${input_query_name}_trim_3.fa

perl ${script_address}combine_reads.pl ${output_address}${input_query_name}_trim_4.fa > ${output_address}${input_query_name}.fa
rm ${output_address}${input_query_name}_trim_4.fa
		';
	}
	if($opt_format){
		print FILE '
tr -d "\15\32" < ${output_address}${input_query_name}_trim_1.${input_query_suffix} > ${output_address}${input_query_name}_trim_2.${input_query_suffix}
rm ${output_address}${input_query_name}_trim_1.${input_query_suffix}
perl ${script_address}filter_length.pl ${output_address}${input_query_name}_trim_2.${input_query_suffix} ${min_length} ${max_length} > ${output_address}${input_query_name}.${input_query_suffix}
rm ${output_address}${input_query_name}_trim_2.${input_query_suffix}
		';
	}

###annotation process - genome
	my $step_number = 1;
	print FILE '

input=${output_address}${input_query_name}.fa	

###step' . $step_number . ': match to genome
echo "match to genome"
bowtie_address=' . $genome_address . '
output_match=${output_address}${input_query_name}_match_genome.fa
output_unmatch=${output_address}${input_query_name}_unmatch_genome.fa
output_detail=${output_address}${input_query_name}_output_match_genome

touch ${output_match}
touch ${output_unmatch}
bowtie ${bowtie_address} -f ${input} -v ${mismatch} -k 1 -p ${thread} --al ${output_match} --un ${output_unmatch} > ${output_detail}

input_match=${output_address}${input_query_name}_match_genome.fa
input_unmatch=${output_address}${input_query_name}_unmatch_genome.fa';

###annotation process - rRNA
	unless ($rRNA_db_address eq "NULL"){
		$step_number += 1;
		print FILE '
###step' . $step_number . ': match to rRNA database
echo ""
echo "match to rRNA database"
output_detail_match_genome=${output_address}${input_query_name}_output_rRNA_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_rRNA_unmatch_genome	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
		my @rRNA = ('2S', '4.5S', '5S', '5.3S', '5.8S', '12S', '16S', '17S', '18S', '25S', '26S', '28S', '45S', 'other', 'RNY1', 'RNY3', 'RNY4', 'RNY5');
		my $bowtie_fa;
		$strand = 0;
		foreach (@rRNA){
			$bowtie_fa = "${rRNA_db_address}_${_}.fa";
			if (-e $bowtie_fa){
				BowtiePrint("rRNA_${_}", "${rRNA_db_address}_${_}", $strand);
			}
		}
	}

###annotation process - tRNA
	unless ($tRNA_db_address eq "NULL"){
		$step_number += 1;
		print FILE '
###step' . $step_number . ': match to tRNA database
echo ""
echo "match to tRNA database"

######match genome part - tRNA-mature
name=tRNA_mature
bowtie_address=' . $tRNA_db_address . '_CCA
echo ""
echo "match to tRNA_mature-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
output_detail_match_genome=${output_address}${input_query_name}_output_${name}_match_genome

touch ${output_detail_match_genome}
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} > ${output_detail_match_genome}

input_match=${output_unmatch_match_genome}

######match genome part - tRNA
name=tRNA_pre
bowtie_address=' . $tRNA_db_address . '
echo ""
echo "match to tRNA-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
output_detail_match_genome=${output_address}${input_query_name}_output_${name}_match_genome

touch ${output_detail_match_genome}
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} > ${output_detail_match_genome}

######unmatch genome part - tRNA-mature
name=tRNA_mature
bowtie_address=' . $tRNA_db_address . '_CCA
echo ""
echo "match to tRNA_mature-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
output_detail_unmatch_genome=${output_address}${input_query_name}_output_${name}_unmatch_genome

touch ${output_detail_unmatch_genome}
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} > ${output_detail_unmatch_genome}

input_unmatch=${output_unmatch_unmatch_genome}

######unmatch genome part - tRNA
name=tRNA_pre
bowtie_address=' . $tRNA_db_address . '
echo ""
echo "match to tRNA-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
output_detail_unmatch_genome=${output_address}${input_query_name}_output_${name}_unmatch_genome

touch ${output_detail_unmatch_genome}
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} > ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}';
	}

###annotation process - mito_tRNA
	if (-e "${tRNA_db_mito_tRNA_file}.fa"){
		$step_number += 1;
		print FILE '
###step' . $step_number . ': match to mito_tRNA database
echo ""
echo "match to mito_tRNA database"
output_detail_match_genome=${output_address}${input_query_name}_output_mt_tRNA_mature_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_mt_tRNA_mature_unmatch_genome
	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
		$strand = 0;
  BowtiePrint('mt_tRNA_mature', "${tRNA_db_mito_tRNA_file}_CCA", $strand);
		print FILE '
echo ""
output_detail_match_genome=${output_address}${input_query_name}_output_mt_tRNA_pre_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_mt_tRNA_pre_unmatch_genome

if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
		$strand = 0;
  BowtiePrint('mt_tRNA_pre', ${tRNA_db_mito_tRNA_file}, $strand);
	}

###annotation process - miRNA
	unless ($miRNA_db_address eq "NULL"){
		$step_number += 1;
		print FILE '
###step' . $step_number . ': match to microRNA database
echo ""
echo "match to microRNA database"
output_detail_match_genome=${output_address}${input_query_name}_output_miRNA_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_miRNA_unmatch_genome
	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
		$strand = 0;
  BowtiePrint('miRNA', $miRNA_db_address, $strand);
	}

###annotation process - ensembl
	unless ($ensembl_nc_address eq "NULL"){
		$step_number += 1;
		print FILE '
###step' . $step_number . ': match to ensembl database
echo ""
echo "match to ensembl database"
output_detail_match_genome=${output_address}${input_query_name}_output_ensembl_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_ensembl_unmatch_genome	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
		$strand = 0;
		BowtiePrint('ensembl', $ensembl_nc_address, $strand);
	}

###annotation process - rfam
	unless ($rfam_address eq "NULL"){
		$step_number += 1;
		print FILE '
###step' . $step_number . ': match to rfam database
echo ""
echo "match to rfam database"
output_detail_match_genome=${output_address}${input_query_name}_output_rfam_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_rfam_unmatch_genome	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
		$strand = 0;
		BowtiePrint('rfam', $rfam_address, $strand);
	}

###annotation process - piRNA
	unless ($piRNA_db_address eq "NULL"){
		$step_number += 1;
		print FILE '
###step' . $step_number . ': match to piRNA database
echo ""
echo "match to piRNA database"
output_detail_match_genome=${output_address}${input_query_name}_output_piRNA_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_piRNA_unmatch_genome	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
		$strand = 0;
		BowtiePrint('piRNA', $piRNA_db_address, $strand);
	}

###annotation process - antisense
	if ($sense){
		unless ($rRNA_db_address eq "NULL"){
			$step_number += 1;
			print FILE '
###step' . $step_number . ': match to rRNA database - antisense
echo ""
echo "match to rRNA database - antisense"
output_detail_match_genome=${output_address}${input_query_name}_output_rRNA-antisense_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_rRNA-antisense_unmatch_genome	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
			my @rRNA = ('2S', '4.5S', '5S', '5.3S', '5.8S', '12S', '16S', '17S', '18S', '25S', '26S', '28S', '45S', 'other', 'RNY1', 'RNY3', 'RNY4', 'RNY5');
			my $bowtie_fa;
			$strand = 1;
			foreach (@rRNA){
				$bowtie_fa = "${rRNA_db_address}_${_}.fa";
				if (-e $bowtie_fa){
					BowtiePrint("rRNA_${_}", "${rRNA_db_address}_${_}", $strand);
				}
			}
		}
		unless ($tRNA_db_address eq "NULL"){
			$step_number += 1;
			print FILE '
###step' . $step_number . ': match to tRNA database - antisense
echo ""
echo "match to tRNA database - antisense"

######match genome part - tRNA-mature - antisense
name=tRNA-antisense_mature
bowtie_address=' . $tRNA_db_address . '_CCA
echo ""
echo "match to tRNA_mature-match_genome - antisense"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
output_detail_match_genome=${output_address}${input_query_name}_output_${name}_match_genome

touch ${output_detail_match_genome}
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --al ${output_match_match_genome} --un ${output_unmatch_match_genome} > ${output_detail_match_genome}

input_match=${output_unmatch_match_genome}

######match genome part - tRNA_antisense
name=tRNA-antisense_pre
bowtie_address=' . $tRNA_db_address . '
echo ""
echo "match to tRNA_antisense-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
output_detail_match_genome=${output_address}${input_query_name}_output_${name}_match_genome

touch ${output_detail_match_genome}
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --al ${output_match_match_genome} --un ${output_unmatch_match_genome} > ${output_detail_match_genome}

######unmatch genome part - tRNA-mature - antisense
name=tRNA-antisense_mature
bowtie_address=' . $tRNA_db_address . '_CCA
echo ""
echo "match to tRNA_mature-unmatch_genome - antisense"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
output_detail_unmatch_genome=${output_address}${input_query_name}_output_${name}_unmatch_genome

touch ${output_detail_unmatch_genome}
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} > ${output_detail_unmatch_genome}

input_unmatch=${output_unmatch_unmatch_genome}

######unmatch genome part - tRNA - antisense
name=tRNA-antisense_pre
bowtie_address=' . $tRNA_db_address . '
echo ""
echo "match to tRNA-unmatch_genome - antisense"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
output_detail_unmatch_genome=${output_address}${input_query_name}_output_${name}_unmatch_genome

touch ${output_detail_unmatch_genome}
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} > ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}';
		}
	if (-e "${tRNA_db_mito_tRNA_file}.fa"){
		$step_number += 1;
		print FILE '
###step' . $step_number . ': match to mito_tRNA database - antisense
echo ""
echo "match to mito_tRNA database - antisense"
output_detail_match_genome=${output_address}${input_query_name}_output_mt_tRNA_mature-antisense_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_mt_tRNA_mature-antisense_unmatch_genome
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
			$strand = 1;
			BowtiePrint('mt_tRNA_mature', "${tRNA_db_mito_tRNA_file}_CCA", $strand);
			print FILE '
output_detail_match_genome=${output_address}${input_query_name}_output_mt_tRNA_pre-antisense_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_mt_tRNA_pre-antisense_unmatch_genome
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
			$strand = 1;
			BowtiePrint('mt_tRNA_pre', ${tRNA_db_mito_tRNA_file}, $strand);
		}
		unless ($miRNA_db_address eq "NULL"){
			$step_number += 1;
			print FILE '
###step' . $step_number . ': match to microRNA database - antisense
echo ""
echo "match to microRNA database - antisense"
output_detail_match_genome=${output_address}${input_query_name}_output_miRNA-antisense_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_miRNA-antisense_unmatch_genome
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
			$strand = 1;
			BowtiePrint('miRNA', $miRNA_db_address, $strand);
			}
		unless ($ensembl_nc_address eq "NULL"){
			$step_number += 1;
			print FILE '
###step' . $step_number . ': match to ensembl database - antisense
echo ""
echo "match to ensembl database - antisense"
output_detail_match_genome=${output_address}${input_query_name}_output_ensembl-antisense_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_ensembl-antisense_unmatch_genome	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
			$strand = 1;
			BowtiePrint('ensembl', $ensembl_nc_address, $strand);
		}
		unless ($rfam_address eq "NULL"){
			$step_number += 1;
			print FILE '
###step' . $step_number . ': match to rfam database - antisense
echo ""
echo "match to rfam database - antisense"
output_detail_match_genome=${output_address}${input_query_name}_output_rfam-antisense_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_rfam-antisense_unmatch_genome	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
			$strand = 1;
			BowtiePrint('rfam', $rfam_address, $strand);
		}
		unless ($piRNA_db_address eq "NULL"){
			$step_number += 1;
			print FILE '
###step' . $step_number . ': match to piRNA database - antisense
echo ""
echo "match to piRNA database - antisense"
output_detail_match_genome=${output_address}${input_query_name}_output_piRNA-antisense_match_genome
output_detail_unmatch_genome=${output_address}${input_query_name}_output_piRNA-antisense_unmatch_genome	
if [ -e "${output_detail_match_genome}" ]; then
rm ${output_detail_match_genome}
fi
if [ -e "${output_detail_unmatch_genome}" ]; then
rm ${output_detail_unmatch_genome}
fi
touch ${output_detail_match_genome}
touch ${output_detail_unmatch_genome}';
			$strand = 1;
			BowtiePrint('piRNA', $piRNA_db_address, $strand);
		}
	}
	
###annotation process - annotation summary
	print FILE '
perl ${script_address}annotation.pl ${output_address}${input_query_name}

if [ ! -d "${output_address}${input_query_name}_result" ]; then
	mkdir -p ${output_address}${input_query_name}_result
fi

if [ ! -d "${output_address}${input_query_name}_processed" ]; then
	mkdir -p ${output_address}${input_query_name}_processed
fi

if [ ! -d "${output_address}${input_query_name}_fa" ]; then
	mkdir -p ${output_address}${input_query_name}_fa
fi

mv ${output_address}${input_query_name}_*.txt ${output_address}${input_query_name}_result
mv ${output_address}${input_query_name}_output_* ${output_address}${input_query_name}_processed
mv ${output_address}${input_query_name}*.f* ${output_address}${input_query_name}_fa';

###mismatch statistics
	if($mismatch_stat){
		if($mismatch > 0){
			if ($miRNA_db_address ne "NULL" || 
				$rRNA_db_address ne "NULL" || 
				$tRNA_db_address ne "NULL" || 
				$ensembl_nc_address ne "NULL" || 
				$rfam_address ne "NULL"  || 
				$piRNA_db_address ne "NULL"){
				print FILE '
echo ""
echo "mismatch loci statistics"
if [ -f "${output_address}${input_query_name}_result/${input_query_name}_mismatch_summary.txt" ]; then
	echo > ${output_address}${input_query_name}_result/${input_query_name}_mismatch_summary.txt
fi
for file in ${output_address}${input_query_name}_processed/*match_genome
do
	if [ -s ${file} ] && [ ${file} != "${output_address}${input_query_name}_processed/${input_query_name}_output_match_genome" ]; then
		perl ${script_address}mismatch_summary.pl ${file} ' . $threshold .' >> ${output_address}${input_query_name}_result/${input_query_name}_mismatch_summary.txt
	fi
done
if [ -s ${output_address}${input_query_name}_result/${input_query_name}_mismatch_summary.txt ]; then
	Rscript --vanilla ${script_address}mismatch_stat.R ${output_address}${input_query_name}_result/${input_query_name}_mismatch_summary.txt ' . $seq_err . '
fi';
			}
		}
	}


###generate output figures
	unless ($miRNA_db_address eq "NULL" && 
		$rRNA_db_address eq "NULL" && 
		$tRNA_db_address eq "NULL" &&
		$piRNA_db_address eq "NULL"){
		print FILE '
echo ""
echo "generating graph"
echo ""
######overall length distribution figure
Rscript --vanilla ${script_address}overall_RNA_length_distribution.R ${output_address} ${input_query_name}
echo "overall length distribution figure generated"
';
	}
	unless ($tRNA_db_address eq "NULL"){
		print FILE '
######pre-tRNA mapping figure
cat ${output_address}${input_query_name}_processed/${input_query_name}_output_*tRNA_pre_*_genome > ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA_pre
if [ -s ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA_pre ]; then
	perl ${script_address}tRNA_mapping.pl ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA_pre ${output_address}${input_query_name}_result/${input_query_name}_summary.txt > ${output_address}${input_query_name}_processed/${input_query_name}_tRNA_pre_mapping.txt
	Rscript --vanilla ${script_address}tRNA_mapping.R ${output_address}${input_query_name}_processed/${input_query_name}_tRNA_pre_mapping.txt ${output_address}${input_query_name}_result/${input_query_name}_tRNA_pre_mapping.pdf
	echo "pre-tRNA mapping figure generated"
fi

######mature-tRNA mapping figure
cat ${output_address}${input_query_name}_processed/${input_query_name}_output_*tRNA_mature_*_genome > ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA_mature
if [ -s ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA_mature ]; then
	perl ${script_address}tRNA_mapping.pl ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA_mature ${output_address}${input_query_name}_result/${input_query_name}_summary.txt > ${output_address}${input_query_name}_processed/${input_query_name}_tRNA_mature_mapping.txt
	Rscript --vanilla ${script_address}tRNA_mapping.R ${output_address}${input_query_name}_processed/${input_query_name}_tRNA_mature_mapping.txt ${output_address}${input_query_name}_result/${input_query_name}_tRNA_mature_mapping.pdf
	echo "mature-tRNA mapping figure generated"
fi';
	}
	unless ($rRNA_db_address eq "NULL"){
		print FILE '
temp_length=' . join(',', @rRNA_length) . '
######rRNA length distribution figure
Rscript --vanilla ${script_address}rRNA_length_distribution.R ${output_address} ${input_query_name} ${temp_length}
echo "rRNA length distribution figure generated"

######rRNA mapping figure
Rscript --vanilla ${script_address}rRNA_mapping.R ${output_address} ${input_query_name} ${temp_length}
echo "rRNA mapping figure generated"';
	}
	unless ($keep_all){
		print FILE '
rm ${output_address}${input_query_name}_processed/${input_query_name}_output_tRNA
rm ${output_address}${input_query_name}_processed/${input_query_name}_tRNA_mapping.txt
rm -rf ${output_address}${input_query_name}_fa
rm -rf ${output_address}${input_query_name}_processed';
	}
	print FILE '
echo ""
date';

	shift @tmp_filelist;
}

close FILE_OUT
	or warn $!;

print "\nProcessing input files...\n\n";

system("sh $output_file");	

unless ($keep_all){
	system("rm -rf ${output_address}processing_report");
	system("rm $output_file");
	system("mv ${output_address}sh_file/*.sh ${output_address}");
	system("rmdir ${output_address}sh_file");
	system("rm ${output_address}*.sh");
}

print "Done!\n";

sub BowtiePrint{
	if ($_[2] == 0){
		print FILE '
name=' . $_[0] . '
bowtie_address=' . $_[1] . '

######match genome part
echo ""
echo "match to ${name}-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}';
	}
	elsif ($_[2] == 1){
		print FILE '
name=' . $_[0] . '
bowtie_address=' . $_[1] . '

######match genome part
echo ""
echo "match to ${name}_antisense-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}-antisense_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}-antisense_match_genome.fa
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -a -p ${thread} --fullref --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}_antisense-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}-antisense_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}-antisense_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -a -p ${thread} --fullref --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}';
	}
}
