#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $version="1.0.0";

#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fa,$db,$od,$mis,$out,$sam);
GetOptions(
				"help|?" =>\&USAGE,
				"fa:s"=>\$fa,
				"db:s"=>\$db,
				"sam:s"=>\$sam,
				"out:s"	=>\$out,
				"mis:n"=>\$mis,
				) or &USAGE;
&USAGE unless ($fa and $db and $out);
################################
my %CFG=%{readconf("$Bin/../../CFG")};
my $bowtie = $CFG{bowtie};


$fa=&ABSOLUTE_DIR($fa);
$db=&ABSOLUTE_DIR($db);
unless (defined $mis) {
	$mis = 2;
}
####################################### align process
if ((!-f "$db.1.ebwt") and   (!-f "$db.1.ebwtl")) {
	print "the reference database is not index, can't run align\n";
	print "please run \"bowtie-build -f Ref.fa\"\n";
	die;
}
print "$db\n";
print "$mis\n";
print "$fa\n";
print "$out\n";
if(-f "$db.1.ebwt"){
	if ($sam) {
		runOrDie("$bowtie $db -v $mis -p 5 -S  -f $fa $out");
	}
	else{
		runOrDie("$bowtie $db -v $mis -p 5  -f $fa $out");
	}
}
elsif(-f "$db.1.ebwtl"){
	if ($sam) {
		runOrDie("$bowtie $db -v $mis -p 5 -S  -f $fa $out --large-index");
	}
	else{
		runOrDie("$bowtie $db -v $mis -p 5  -f $fa $out --large-index");
	}
}
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
	-fa               Uniq fasta file       must be given
	-db               Reference Database    must be given
	-out              output file           must be given
	-sam			  output file format samfile
	-mis              mismatch for align    option,default 2;
USAGE
	print $usage;
	exit;
}
