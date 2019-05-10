#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($od,$query);
GetOptions(
			"help|?"=>\&USAGE,
			"od:s"=>\$od,
			"fa:s"=>\$query,
			)or &USAGE;
&USAGE unless ($query and $od);
`mkdir $od ` unless (-d $od);

$od=ABSOLUTE_DIR($od);
$query=ABSOLUTE_DIR($query);
my $name=basename$query;

$/=">";
my %sequence;
my %sequence_num;
open (IN,$query)||die $!;
while (<IN>) {
	chomp;
	next if ($_=~/^\s*$/);
	my($id,$seq)=split(/\n+/,$_,2);
	$seq=~s/\s+//g;
	if ($id=~/^\d+\s*$/) {
		$sequence_num{$id}=$seq;
	}else{
		$sequence{$id}=$seq;
	}
}
$/="\n";
close(IN);

##################################################
open OUT,">$od/$name.fa" ||die $!;
if ((keys %sequence_num)!=0) {
	foreach my $id (sort {$a<=>$b} keys %sequence_num) {
		print OUT ">$id\n$sequence_num{$id}\n";
	}
}
if ((keys %sequence)!=0) {
	foreach my $id (sort keys %sequence) {
		print OUT ">$id\n$sequence{$id}\n";
	}
}
close(OUT);

################################################################################
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
ProgramName:  blast_supply
Version:	$version
Contact:	Wang Yajing <wangyj\@biomarker.com.cn> 
Program Date:   2015/5/27
Usage:
Options:
		-od     output dir
		-fa     query file
		-h      Help

USAGE
	print $usage;
	exit;
}
