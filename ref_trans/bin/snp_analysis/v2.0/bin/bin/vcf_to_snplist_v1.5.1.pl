#!/usr/bin/perl 
# Copyright (c) BMK 2014/3/4
# Writer:         lium <lium@biomarker.com.cn>
# Program Date:   2014/3/4
# Last Modified:  2014/8/28   shitw 1.5 to 1.5.1
my $ver="1.5.1";

use strict;
use Cwd;
use Getopt::Long;
my $BEGIN=time();


# get opts
my %opts;
GetOptions(\%opts,"i=s","o=s","ref=i","h" );
if(! defined($opts{i}) ||! defined($opts{o})||! defined($opts{'ref'}) || defined($opts{h})){
	&USAGE;
}


###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";


# open file
open (IN,"$opts{i}") || die "Can't open file $opts{i}\n";
open (OUT,">$opts{o}") || die "open or create file $opts{o} is failed\n";


# define base encode
my %type = ('AG'=>'R', 'GA'=>'R', 'CT'=>'Y', 'TC'=>'Y', 'GT'=>'K', 'TG'=>'K', 'AC'=>'M', 'CA'=>'M', 
	'CG'=>'S', 'GC'=>'S', 'AT'=>'W', 'TA'=>'W', 'A'=>'A', 'C'=>'C', 'G'=>'G', 'T'=>'T', 
	'CGT'=>'B', 'AGT'=>'D', 'ACT'=>'H', 'ACG'=>'V', 'ACGT'=>'N',
);


# get individual info
my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format) = ();
my $get_indi_ok = 0;
my @indi = ();
my $indi_num = 0;
while(<IN>){
	chomp;
	next if (/^$/);		# ignore the blank lines
	# check individual info
	if(/^\#CHROM/){
		($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=(split/\s+/,$_);
		$indi_num = @indi;
		$get_indi_ok = 1;
		last;
	}
	# check commend
	next if(/^\#/); 	# ignore the commend
	# error: get individual info is failed
	print "error1: get individual info is failed\n";
	exit;
}


# check get indi result
if($get_indi_ok==0){
	print "error2: get individual info is failed\n";
	exit;
}


# print header
print OUT "#Chr\tPos\t";
if($opts{'ref'}){print OUT "Ref\tAlt";}
for (my $i=0; $i<@indi; $i++){
	print OUT "\t$indi[$i]\tDepth\tAlleDp" ;
}
print OUT "\tEffect\tCodon_change\tGene_id\n";


# abstract snplist line by line
my @flag = ();
my @twobase = ("","",""); # for sample alleles not the same as ref
while(<IN>){
	chomp;
	next if (/^$/);		# ignore the blank lines
	# abstract info from one line
	($chr,$pos,$id,$twobase[0],$twobase[1],$qual,$filter,$info,$format,@flag)=(split/\s+/,$_);
	next if ($format !~ /AD/ || $format !~ /DP/) ;
	# ignore indel
	#next if( (length($twobase[0])!=1) || (length($twobase[1])!=1) );
	next if (length($twobase[0]) != 1);
	my $raw_alt = $twobase[1] ;
	#next if (length($twobase[1]) != 1) ;  # filter snp diff than ref
	if (length($twobase[1]) != 1){
		my @bases = split /\,/, $twobase[1] ;
		my $flag = 0 ;
		for my $base (@bases){
			if (length $base != 1){
				$flag = 1 ;
				last ;
			}
		}
		next if ($flag == 1) ;
		for (my $i=0; $i<@bases; $i++){
			$twobase[$i+1] = $bases[$i] ;
		}
	}
	# output chr and pos
	print OUT "$chr\t$pos";
	# output ref if "-ref is show"
	if( $opts{'ref'} ){
		print OUT "\t$twobase[0]\t$raw_alt";
	}
	# get DP position
		my $dp_pos="NA";
		my $gt_pos="NA";
		my $ad_pos="NA";
		my @form=split /\:/,$format;
		for (my $i=0;$i<@form ;$i++) {
			if ($form[$i] eq "DP") {
				$dp_pos=$i;
			}
			elsif ($form[$i] eq "GT") {
				$gt_pos=$i;
			}
			elsif ($form[$i] eq "AD") {
				$ad_pos=$i;
			}
		}
	# iter process each indi
	foreach my $f (@flag){
		my @sam=split /\:/,$f;
		my $gt=$sam[$gt_pos];
		my $dp=$sam[$dp_pos];
		my $ad=$sam[$ad_pos];
		if ($gt_pos eq "NA") {
			$gt=".";
		}
		if ($dp_pos eq "NA") {
			$dp=".";
		}
		if ($ad_pos eq "NA") {
			$ad=".";
		}
		my @bases = split /\//, $gt ;
		#my @ads = split /,/, $ad ;
		my %hbase = ();
		for (my $i=0; $i<@bases; $i++){
			next if ($bases[$i] eq '.');
			$hbase{$twobase[$bases[$i]]} ++ ;
		}
		my @f_base = keys %hbase ;
		# for loss
		if (@f_base == 0){
			print OUT "\tN\t.\t." ;
			next ;
		}
		# for normal base
		my $t = join '', sort @f_base ;
		if(!exists $type{$t}){
			print "key($t) not exists\n";
			exit;
		}else{
			print OUT "\t$type{$t}\t$dp\t$ad";
		}
	}
	if (m/SNPEFF_EFFECT=(.*?);/){
		print OUT "\t$1" ;
		if (m/SNPEFF_CODON_CHANGE=(\D+?);/){
			print OUT "\t$1" ;
			if (m/SNPEFF_GENE_NAME=(.*?);/){
				(my $name = $1) =~ s/\.\d+$// ;
				print OUT "\t$name" ;
			}
		}
	}
	else{
		print OUT "\tUNKNOWN" ;
	}
	# print \n
	print OUT "\n";
}

close OUT;
close IN;


###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";
&Runtime($BEGIN);


sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub Runtime{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n\n";
}




sub USAGE{
print <<"Usage End.";
Description: convert vcf format snp file to snplist file.
	Version: $ver
	v1.0: function realization.
	v1.1: bug fixed for output snps diff from ref(eg: sample A/G, ref T).
	v1.2: support haploid and multiploid(pre-version only support diploid).
	v1.5: change transcript name to gene name.
	      ignore snp diff with ref.
	v1.5.1: fixed bug of some loci do not contain AD or DP
Usage: perl vcf_to_snplist.pl -i test.vcf -o out.hete -ref 0
	-i   SNP rawdata file(vcf format)		must be given

	-o   the filename of snplist			must be given

	-ref output ref base (0:no; 1:yes)		must be given

	-h   Help document
Usage End.
exit;
}


