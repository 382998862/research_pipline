#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";

#######################################################################################
# ------------------------------------------------------------------
# GetOptions   
# ------------------------------------------------------------------
my ($align,$key,$od);
my $keywords="$Bin/keywords.txt";
GetOptions(
				"help|?" =>\&USAGE,
				"align:s"=>\$align,
				"key:s"=>\$key,
				"od:s"=>\$od,
				) or &USAGE;
&USAGE unless ($align and $key);
################################
$align=&ABSOLUTE_DIR($align);
$od=$od || "./";
` mkdir $od ` if (!-d "$od");
$od=&ABSOLUTE_DIR($od);
open(IN,"$keywords")||die $!;
my %keywords;
while(<IN>){
	chomp;
	my @line=split/\t/,$_;
	$keywords{$line[0]}=$line[2];
}
close IN;
my $rRNA=keyword(\%keywords,"ribosomal RNA;rRNA");
$rRNA=substr($rRNA,0,length($rRNA)-1);
#print "$rRNA\n";
my $tRNA=keyword(\%keywords,"transfer RNA;tRNA");
$tRNA=substr($tRNA,0,length($tRNA)-1);
#print "$tRNA\n";
my $snRNA=keyword(\%keywords,"small nuclear;snRNA");
$snRNA=substr($snRNA,0,length($snRNA)-1);
#print "$snRNA\n";
my $snoRNA=keyword(\%keywords,"Small nucleolar;snoRNA");
$snoRNA=substr($snoRNA,0,length($snoRNA)-1);
sub keyword{
	my $keyword=shift;
	my $word=shift;
	my @word=split/;/,$word;
	my $word1=$word[0];
	my $word2=$word[1];
	my $return;
	my %keyword=%$keyword;
	foreach(keys %keyword){
		if($keyword{$_}=~/ctRNA/){next;}
		if($keyword{$_}=~/$word1/){
			$return.=$_."|";
		}
		elsif($keyword{$_}=~/$word2/){
			$return.=$_."|";
		}
	}
	return $return;
}

####################################### align process
open IN,"$align" || die $!;
open OUT,">$od/$key.bowtie_Rfam.list" || die $!;
while (<IN>) {
	chomp;
	next if (/^$/);
	my @tmp=split /\s+/,$_;
	my $index;
	if ($tmp[2]=~/$rRNA/) {
		#print "$tmp[2]\n";
		$index="rRNA";
	}
	elsif ($tmp[2]=~/$tRNA/) {
		#print "tRNA\n";
		$index="tRNA";
	}
	elsif ($tmp[2]=~/$snRNA/) {
		#print "$tmp[2]\n";
		$index="snRNA";
	}
	elsif ($tmp[2]=~/$snoRNA/) {
		$index="snoRNA";
	}else{
		next;
	}
	print OUT "$tmp[0]\t$index\t$tmp[1]\t$tmp[2]\n";
}
close OUT;


######################################### subs
sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: 
Version: $version
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:
	
Usage:
	-align               Uniq fasta file       must be given
	-key                 sample name           must be given
	-od                  output dir            choice default "./";
USAGE
	print $usage;
	exit;
}
