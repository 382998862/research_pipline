#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my ($ep,$em,$ip,$im,$out,$name);
GetOptions(
				"help|?" =>\&USAGE,
				"ep:s"=>\$ep,
				"em:s"=>\$em,
				"ip:s"=>\$ip,
				"im:s"=>\$im,
				"o:s"=>\$out,
				"name:s"=>\$name,
				) or &USAGE;
&USAGE unless ($ep and $em and $ip and $im and $out);
$ep=&ABSOLUTE_DIR($ep);
$em=&ABSOLUTE_DIR($em);
$ip=&ABSOLUTE_DIR($ip);
$im=&ABSOLUTE_DIR($im);
my %hash_ei;my %hash_gene;
###ÓÅÏÈÕýÁ´--exon
my $plus="plus";my $minus="minus";
&read_hash($ep,$plus);
&read_hash($ip,$plus);
&read_hash($em,$minus);
&read_hash($im,$minus);
my $exon_p;my $exon_m;my $intron_p;my $intron_m;
foreach(keys %hash_ei){
	my $reads=1;
	if($_=~/_x(\d+)/){
		$reads=$reads+$1-1;
	}
	if($hash_ei{$_} eq "exon_plus"){
		
		$exon_p+=$reads;
	}
	elsif($hash_ei{$_} eq "exon_minus"){
		$exon_m+=$reads;
	}
	elsif($hash_ei{$_} eq "intron_plus"){
		$intron_p+=$reads;
	}
	elsif($hash_ei{$_} eq "intron_minus"){
		$intron_m+=$reads;
	}
}
open(OU,">$out");
print OU "sample\texon_same\texon_reverse\tintron_same\tintron_reverse\n";
print OU "$name\t$exon_p\t$exon_m\t$intron_p\t$intron_m\n";
close OU;

sub read_hash{
	my $file=$_[0];
	my $direc=$_[1];
	open(IN,"$file");
	while(<IN>){
		my @line=split/\s+/,$_;
		
		unless(exists $hash_ei{$line[12]}){
			my $val=join("_",$line[3],$direc);
			$hash_ei{$line[9]}=$val;
		}
		
	}
	close IN;
	
}

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}
sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Contact:	Liuxiaoshuang <liuxs\@biomarker.com.cn> 
Program Date:   2015.8.4
Usage:
  Options:
  -ep         	<file>    mapper result of exon positive direction bed file,forced 
 
  -em       	<file>    mapper result of exon negative direction bed file,forced 
 
  -ip       	<file>    mapper result of intron positive direction bed file,forced 
 
  -im  			<file>     mapper result of intron negative direction bed file,forced 
 
  -o      		<file>     output file,forced
 
  -name			<file>	sample name.forced
  -h         	Help

USAGE
	print $usage;
	exit;
}
