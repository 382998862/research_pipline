#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($od,$data,$list,$type);
GetOptions(
			"help|?"=>\&USAGE,
			"od:s"=>\$od,
			"data:s"=>\$data,
			"list:s"=>\$list,
			"type:s"=>\$type,
			)or &USAGE;
&USAGE unless ($data and $list );
$od||="./";

`mkdir $od ` unless (-d $od);

$od=ABSOLUTE_DIR($od);
$data=ABSOLUTE_DIR($data);
my $basename=basename$list;
$type||="extract";
#################################################
my $data_type;

open (IN,$data)||die $!;
while (<IN>) {
	chomp;
	if ($_=~/^>/) {
	$data_type='fa';
	}
	else{$data_type='list';}
	last;
}

##################################################

my %list;
open (L,$list)||die $!;
while (<L>) {
	chomp;
	if ($type eq 'combine') {
		my ($name,$inform)=split(/\s+/,$_,2);
		$list{$name}=$inform;
	}
	else {
		my ($name)=split(/\s+/,$_);
		$list{$name}=1;
	}
}
close (L);
###################################################
if ( $data_type eq 'fa' and $type eq 'combine') {
	print "Sorry!You choice is wrong ,pleace Check!";die;
}


seek (IN,0,0);
if ($data_type eq 'fa') {
	open (OUT,">$od/$basename.fa") ||die $!;
	$/=">";
	while(<IN>)
	{
		chomp;
		next if ($_=~/^\s*$/);
		my ($name,$seq)=split(/\n/,$_,2);
		($name,undef)=split(/\s+/,$name);
		if (exists $list{$name}) {
			if ($type eq 'extract') {

				print OUT ">$name\n$seq";
			}
		}
		else {
			if ($type eq 'filter') {

				print OUT ">$name\n$seq";
			}

		}
=c
		print OUT ">$name\n$seq" if (exists $list{$name} and $type eq 'extract'); 
		unless (exists $list{$name}) {
			if ($type eq 'filter') {
				print OUT ">$name\n$seq";
			}
		}
=cut
	}
	$/="\n";
}
else {
	open (OUT,">$od/$basename.list") ||die $!;
	while(<IN>)
	{
		chomp;
		next if ($_=~/^\s*$/);
		print OUT "$_\n" if ($_=~/^#/);
		my ($name,$seq)=split(/\s+/,$_,2);
#		print OUT "$name\t$list{$name}\t$seq\n" if (exists $list{$name} and $type eq 'combine' );
		if (exists $list{$name}) {
			if ($type eq 'extract') {

				print OUT "$name\t$seq\n";
			}
			if ($type eq 'combine') {

				print OUT "$name\t$list{$name}\t$seq\n";
			}
		}
		else {
			if ($type eq 'filter') {

				print OUT "$_\n";
			}
			if ($type eq 'combine' ) {

				print OUT "$name\t--\t--\t$seq\n";
			}

		}
	}

}
close(IN);
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
ProgramName:  extract_list
Version:	$version
Contact:	Wang Yajing <wangyj\@biomarker.com.cn> 
Program Date:   2014/12/11
Usage:
Options:
-od     output dir
-data      <file> list or fa 
-list       <file>   list 
-type       <str>    extract or filter or combine(-data list),default extract
-h      Help

USAGE
	print $usage;
	exit;
}
