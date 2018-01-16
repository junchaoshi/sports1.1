#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Std;
use vars qw($opt_i $opt_a $opt_b $opt_x $opt_y $opt_l $opt_L $opt_M $opt_p $opt_g $opt_m $opt_r $opt_t $opt_e $opt_f $opt_w $opt_o $opt_k $opt_v $opt_h);
getopts('i:f:abx:y:l:L:M:p:s:g:m:r:t:e:f:w:o:kvh');
my $input_file		= $opt_i;
my $opt_adapter		= $opt_a ? 1 : 0;
my $opt_format		= $opt_b ? 1 : 0;
my $adapter_5		= $opt_x ? $opt_x : "default";
my $adapter_3		= $opt_y ? $opt_y : "default";
my $min_length		= $opt_l ? $opt_l : 15;
my $max_length		= $opt_L ? $opt_L : 45;
my $mismatch		= $opt_M ? $opt_M : 0;
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
my $version		= $opt_v ? 1 : 0;
my $help		= $opt_h ? 1 : 0;


my $version_info = "1.0.1";

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
  -m <str>	miRNA database bowtie index
  -r <str>	rRNA database bowtie index
  -t <str>	tRNA database bowtie index
  -e <str>	ensembl noncoding RNA database bowtie index (optional)
  -f <str>	rfam database bowtie index (optional)
  -w <str>	piRNA database bowtie index (optional)

Output:
  -o <str>	output address of annotation results (default: input address)
  -k		keep all the intermediate files generated during the running progress

Alignment:
  -l <int>	the minimal length of the output sequences (default = 15)
  -L <int>	the maximal length of the output sequences (default = 45)
  -a		Remove 5\'/3\' adapters
  -x <str>	(When -a is selected) Your 5\' adapter sequence or use "default". Default = "GTTCAGAGTTCTACAGTCCGACGATC"
  -y <str>	(When -a is selected) Your 3\' adapter sequence or use "default". Default = "TGGAATTCTCGGGTGCCAAGG"

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

my $script_address = `which sports.pl`;
   @input = split(/\//, $script_address);
   pop @input;
   $script_address = join('/', @input) . '/';

##generate tRNA unmatch genome database
unless($tRNA_db_address eq "NULL"){
	my $tRNA_db_UMG_file = $tRNA_db_address . "_CCA.1.ebwt";
	unless (-e $tRNA_db_UMG_file){
		system ("perl ${script_address}tRNA_db_processing.pl ${tRNA_db_address}.fa");
		print "\n\nGenerating tRNA (unmatch genome) database bowtie index...\n\n";
		system ("bowtie-build -q ${tRNA_db_address}_CCA.fa ${tRNA_db_address}_CCA");
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
	$output_address = abs_path($output_address);
}

{
@input = split(/\//, $output_address);
$output_address = join('/', @input) . '/';
}

my $output_file = $output_address . 'run_this.sh';

open FILE_OUT, ">$output_file" || die $!;

my $count = 0;
my @rRNA_length;

unless ($rRNA_db_address eq "NULL"){
	my @rRNA = ('2S', '4.5S', '5S', '5.3S', '5.8S', '12S', '16S', '17S', '18S', '25S', '26S', '28S', '45S');
	my $bowtie_fa;
	my $temp_length;
	foreach (@rRNA){
		$bowtie_fa = "${rRNA_db_address}_${_}.fa";
		if (-e $bowtie_fa){
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
	print FILE_OUT "echo $input_query_address$input_query_name.$input_query_suffix\n";
	print FILE_OUT "sh ${tmp_sh}/${count}_$input_query_name.sh > $tmp_repo${count}_${input_query_name}.txt 2>&1\n";
	open FILE, '>',"${tmp_sh}/${count}_$input_query_name.sh"
		or die "can not open '${tmp_sh}/${count}_$input_query_name.sh";
	print FILE '#!/bin/bash	
##This script uses to annotate small RNA
date
thread=' . $thread . '
adapter5=' . $adapter_5 . '
adapter3=' . $adapter_3 . '
min_length=' . $min_length . '
max_length=' . $max_length . '
mismatch=' . $mismatch . '
input_query_address=' . $input_query_address . '
input_query_name=' . $input_query_name . '
input_query_suffix=' . $input_query_suffix . '
output_address=' . $output_address . ${count} . '_${input_query_name}/
script_address=' . $script_address . '


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
cutadapt -g ' . $adapter_5 . ' -a ' . $adapter_3 . ' -o ${output_address}${input_query_name}_trim_1.${input_query_suffix} --max-n 0 ${output_address}${input_query_name}.${input_query_suffix}
rm -rf ${output_address}${input_query_name}.${input_query_suffix}
		';
		}
		elsif ($adapter_5 !~ /^default$/i && 
		       $adapter_3 =~ /^default$/i){
		print FILE '
echo ""
echo "remove 5\' adapter"
cutadapt -g ' . $adapter_5 . ' -o ${output_address}${input_query_name}_trim_1.${input_query_suffix} --max-n 0 ${output_address}${input_query_name}.${input_query_suffix}
rm -rf ${output_address}${input_query_name}.${input_query_suffix}
		';
		}
		elsif ($adapter_5 =~ /^default$/i && 
		       $adapter_3 !~ /^default$/i){
		print FILE '
echo ""
echo "remove 3\' adapter"
cutadapt -a ' . $adapter_3 . ' -o ${output_address}${input_query_name}_trim_1.${input_query_suffix} --max-n 0 ${output_address}${input_query_name}.${input_query_suffix}
rm -rf ${output_address}${input_query_name}.${input_query_suffix}
		';
		}
		elsif ($adapter_5 !~ /^default$/i && 
		       $adapter_3 !~ /^default$/i){
		print FILE '
echo ""
echo "remove 5-end and 3-end adapters"
cutadapt -g ' . $adapter_5 . ' -a ' . $adapter_3 . ' -o ${output_address}${input_query_name}_trim_1.${input_query_suffix} --max-n 0 ${output_address}${input_query_name}.${input_query_suffix}
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
tr -d "\15\32" < ${output_address}${input_query_name}_trim_1.${input_query_suffix} > ${output_address}${input_query_name}.${input_query_suffix}
rm ${output_address}${input_query_name}_trim_1.${input_query_suffix}
		';
	}

###annotation process
	print FILE '

input=${output_address}${input_query_name}.fa	

###step1: match to genome
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

	unless ($miRNA_db_address eq "NULL"){
		print FILE '
###step2: match to microRNA database
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
		BowtiePrint('miRNA', $miRNA_db_address);
	}
	
	unless ($rRNA_db_address eq "NULL"){
		print FILE '
###step3: match to rRNA database
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
		my @rRNA = ('2S', '4.5S', '5S', '5.3S', '5.8S', '12S', '16S', '17S', '18S', '25S', '26S', '28S', '45S', 'other');
		my $bowtie_fa;
		foreach (@rRNA){
			$bowtie_fa = "${rRNA_db_address}_${_}.fa";
			if (-e $bowtie_fa){
				BowtiePrint("rRNA_${_}", "${rRNA_db_address}_${_}");
			}
		}
	}
	unless ($tRNA_db_address eq "NULL"){
		print FILE '
###step4: match to tRNA database
echo ""
echo "match to tRNA database"

######match genome part - tRNA
name=tRNA
bowtie_address=' . $tRNA_db_address . '
echo ""
echo "match to tRNA-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
output_detail_match_genome=${output_address}${input_query_name}_output_${name}_match_genome

touch ${output_detail_match_genome}
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -k 10000 -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} > ${output_detail_match_genome}

input_match=${output_unmatch_match_genome}

######match genome part - tRNA-CCA
name=tRNA_CCA
bowtie_address=' . $tRNA_db_address . '_CCA
echo ""
echo "match to tRNA_CCA-match_genome"
output_match_match_genome=${output_address}${input_query_name}_match_${name}_match_genome.fa
output_unmatch_match_genome=${output_address}${input_query_name}_unmatch_${name}_match_genome.fa
output_detail_match_genome=${output_address}${input_query_name}_output_${name}_match_genome

touch ${output_detail_match_genome}
touch ${output_match_match_genome}
touch ${output_unmatch_match_genome}

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -k 10000 -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} > ${output_detail_match_genome}

######unmatch genome part - tRNA
name=tRNA
bowtie_address=' . $tRNA_db_address . '
echo ""
echo "match to tRNA-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
output_detail_unmatch_genome=${output_address}${input_query_name}_output_${name}_unmatch_genome

touch ${output_detail_unmatch_genome}
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -k 10000 -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} > ${output_detail_unmatch_genome}

input_unmatch=${output_unmatch_unmatch_genome}

######unmatch genome part - tRNA-CCA
name=tRNA_CCA
bowtie_address=' . $tRNA_db_address . '_CCA
echo ""
echo "match to tRNA_CCA-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
output_detail_unmatch_genome=${output_address}${input_query_name}_output_${name}_unmatch_genome

touch ${output_detail_unmatch_genome}
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -k 10000 -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} > ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}';
	}

	unless ($piRNA_db_address eq "NULL"){
		print FILE '
###step5: match to piRNA database
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
		BowtiePrint('piRNA', $piRNA_db_address);
	}

	unless ($ensembl_nc_address eq "NULL"){
		print FILE '
###step6: match to ensembl database
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
		BowtiePrint('ensembl', $ensembl_nc_address);
	}

	unless ($rfam_address eq "NULL"){
		print FILE '
###step7: match to rfam database
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
		BowtiePrint('rfam', $rfam_address);
	}

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
mv ${output_address}${input_query_name}*.fa* ${output_address}${input_query_name}_fa';

	unless ($miRNA_db_address eq "NULL" && 
		$rRNA_db_address eq "NULL" && 
		$tRNA_db_address eq "NULL" &&
		$piRNA_db_address eq "NULL"){
		print FILE '
echo ""
echo "generating graph"
Rscript --vanilla ${script_address}overall_RNA_length_distribution.R ${output_address} ${input_query_name}';
	}

	unless ($rRNA_db_address eq "NULL"){
		print FILE '
temp_length=' . join(',', @rRNA_length) . '
Rscript --vanilla ${script_address}rRNA_length_distribution.R ${output_address} ${input_query_name} ${temp_length}
Rscript --vanilla ${script_address}rRNA_mapping.R ${output_address} ${input_query_name} ${temp_length}';
	}
	unless ($keep_all){
		print FILE '
rm -rf ${output_address}${input_query_name}_fa
rm -rf ${output_address}${input_query_name}_processed
		';
	}
	print FILE '

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
	system("mv  ${output_address}sh_file/*.sh ${output_address}");
	system("rmdir ${output_address}sh_file");
	system("rm ${output_address}*.sh");
}

print "Done!\n";

sub BowtiePrint{
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

bowtie ${bowtie_address} -f ${input_match} -v ${mismatch} -k 1 -p ${thread} --fullref --norc --al ${output_match_match_genome} --un ${output_unmatch_match_genome} >> ${output_detail_match_genome}

######unmatch genome part
echo ""
echo "match to ${name}-unmatch_genome"
output_match_unmatch_genome=${output_address}${input_query_name}_match_${name}_unmatch_genome.fa
output_unmatch_unmatch_genome=${output_address}${input_query_name}_unmatch_${name}_ummatch_genome.fa
touch ${output_match_unmatch_genome}
touch ${output_unmatch_unmatch_genome}

bowtie ${bowtie_address} -f ${input_unmatch} -v ${mismatch} -k 1 -p ${thread} --fullref --norc --al ${output_match_unmatch_genome} --un ${output_unmatch_unmatch_genome} >> ${output_detail_unmatch_genome}

######define next input
input_match=${output_unmatch_match_genome}
input_unmatch=${output_unmatch_unmatch_genome}';
}
