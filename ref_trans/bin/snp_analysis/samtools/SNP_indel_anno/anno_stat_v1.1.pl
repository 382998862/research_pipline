#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2014
# Writer:         Dengdj <dengdj@biomarker.com.cn>
# Program Date:   2014
# Modifier:       Dengdj <dengdj@biomarker.com.cn>
# Last Modified:  2014.
my $ver="1.1";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作

my %opts;
GetOptions(\%opts,"i=s","id=s","o=s","h" );

#&help()if(defined $opts{h});
if((!defined($opts{i}) && !defined($opts{id})) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		Static annotation result of SNP/INDEL.
		Version: $ver

	Usage:

		-id          indir of annotation file        <indir>            -id/-i must have one
		-i           annotation file of GATK         <infile>           -id/-i must have one
		-o           statistic result                <outfile>          must be given
		-h           Help document

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################
# get parameters
my ($infile, $indir) ;
if (defined $opts{i}){
	$infile = $opts{i} ;
}
if (defined $opts{id}){
	$indir = $opts{id} ;
	$indir = &ABSOLUTE_DIR($indir);
}
my $outfile = $opts{o} ;
my %hanno = ();
my %htype = ();
my %hinfo = (
	'NONE'=>'--', 'CHROMOSOME'=>'--', 'CUSTOM'=>'--', 'CDS'=>'CDS', 'INTERGENIC'=>'--',
	'INTERGENIC_CONSERVED'=>'--', 'UPSTREAM'=>'--', 'UTR_5_PRIME'=>'--', 'UTR_5_DELETED'=>'--',
	'START_GAINED'=>'--', 'SPLICE_SITE_ACCEPTOR'=>'--', 'SPLICE_SITE_DONOR'=>'--', 'SPLICE_SITE_REGION'=>'--',
	'INTRAGENIC'=>'--', 'START_LOST'=>'CDS', 'SYNONYMOUS_START'=>'CDS', 'NON_SYNONYMOUS_START'=>'CDS',
	'GENE'=>'--', 'TRANSCRIPT'=>'--', 'EXON'=>'CDS', 'EXON_DELETED'=>'CDS', 'NON_SYNONYMOUS_CODING'=>'CDS',
	'SYNONYMOUS_CODING'=>'CDS', 'FRAME_SHIFT'=>'CDS', 'CODON_CHANGE'=>'CDS', 'CODON_INSERTION'=>'CDS',
	'CODON_CHANGE_PLUS_CODON_INSERTION'=>'CDS', 'CODON_DELETION'=>'CDS', 'CODON_CHANGE_PLUS_CODON_DELETION'=>'CDS',
	'STOP_GAINED'=>'CDS', 'SYNONYMOUS_STOP'=>'CDS', 'STOP_LOST'=>'CDS', 'RARE_AMINO_ACID'=>'CDS', 'INTRON'=>'--',
	'INTRON_CONSERVED'=>'--', 'UTR_3_PRIME'=>'--', 'UTR_3_DELETED'=>'--', 'DOWNSTREAM'=>'--',
	'REGULATION'=>'--', 'Other'=>'--'
);
my %horder = (
	'INTERGENIC'=>'1', 'INTERGENIC_CONSERVED'=>'2', 'INTRAGENIC'=>'3', 'INTRON'=>'4', 'INTRON_CONSERVED'=>'5', 
	'UPSTREAM'=>'6', 'DOWNSTREAM'=>'7', 'UTR_5_PRIME'=>'8', 'UTR_5_DELETED'=>'9','UTR_3_PRIME'=>'10', 
	'UTR_3_DELETED'=>'11', 'SPLICE_SITE_ACCEPTOR'=>'12', 'SPLICE_SITE_DONOR'=>'13', 'SPLICE_SITE_REGION'=>'14', 
	'START_GAINED'=>'15', 'START_LOST'=>'16', 'SYNONYMOUS_START'=>'17', 'NON_SYNONYMOUS_START'=>'18',
	'SYNONYMOUS_CODING'=>'19', 'NON_SYNONYMOUS_CODING'=>'20', 'SYNONYMOUS_STOP'=>'21','STOP_GAINED'=>'22', 
	'STOP_LOST'=>'23', 'FRAME_SHIFT'=>'17', 'CODON_CHANGE'=>'18', 'CODON_INSERTION'=>'18', 'CODON_DELETION'=>'18',
	'EXON_DELETED'=>'18', 'CODON_CHANGE_PLUS_CODON_INSERTION'=>'19', 'CODON_CHANGE_PLUS_CODON_DELETION'=>'19', 'RARE_AMINO_ACID'=>'20',
	'Other'=>'24', 'GENE'=>'25', 'TRANSCRIPT'=>'26', 'EXON'=>'27', 
	'REGULATION'=>'28', 'NONE'=>'29', 'CHROMOSOME'=>'30', 'CUSTOM'=>'31', 'CDS'=>'32'
);
#my @region_type = ('--', 'CDS');

# reading vcf file
if (defined $opts{id}){
	my @vcffiles = glob("$indir/*.vcf");
	for my $vcf_file(@vcffiles){
		&stat_vcf_file($vcf_file, \%hanno, \%htype);
	}
}
if (defined $opts{i}){
	&stat_vcf_file($infile, \%hanno, \%htype);
}

# output result
&output_result($outfile, \%hanno, \%htype, \%hinfo, \%horder);

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;
	$cur_dir =~ s/\n$//;
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;
		$dir =~ s/\n$// ;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;
		$return =~ s/\n$// ;
	}
	else
	{
		warn "Warning just for file and dir\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

# &show_log("txt")
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}

#&run_or_die($cmd);
sub run_or_die()
{
	my ($cmd) = @_ ;
	&show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		&show_log("Error: command fail: $cmd");
		exit(1);
	}
	&show_log("done.");
	return ;
}

## qsub

#&stat_vcf_file($vcf_file, \%hanno, \%htype);
sub stat_vcf_file()
{
	my ($vcf_file, $ahanno, $ahtype) = @_ ;
	open (IN, $vcf_file) || die "Can't open $vcf_file, $!\n" ;
	my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format) = ();
	my $get_indi_ok = 0;
	my @indi = ();
	my $indi_num = 0;
	my $sample_name = "all" ;
	while(<IN>){
		chomp ;
		next if (m/^\s*$/);
		if (m/^\#CHROM/){
			($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=(split/\s+/,$_);
			$indi_num = @indi;
			if ($indi_num == 1){
				$sample_name = $indi[0] ;
			}
		}
		next if (m/^\#/);
		if (m/SNPEFF_EFFECT=(.*?)\;/){
			$ahanno->{$sample_name}->{$1}++ ;
			$ahtype->{$1}++ ;
		}
		else{
			$ahanno->{$sample_name}->{'Other'}++ ;    # un-annotation snp number
		}
	}
	close(IN);
	return ;
}

#&output_result($outfile, \%hanno, \%htype, \%hinfo, \%horder);
sub output_result()
{
	my ($outfile, $ahanno, $ahtype, $ahinfo, $ahorder) = @_ ;
	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	my @types = sort {$ahorder->{$a} <=> $ahorder->{$b}} keys %{$ahtype} ;
	push @types, 'Other' ;
	my @samples = sort keys %{$ahanno} ;
	print OUT "#Region\tType\t", join ("\t", @samples), "\n" ;
	for my $type (@types){
		print OUT $ahinfo->{$type}, "\t$type" ;
		for my $sample (@samples){
			my $num = defined $ahanno->{$sample}->{$type} ? $ahanno->{$sample}->{$type} : 0 ;
			print OUT "\t", $num ;
		}
		print OUT "\n" ;
	}
	close(OUT);

	return ;
}


