#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use autodie;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0";
my @Times = localtime();
#######################################################################################
my $Time_Start = sub_format_datetime(localtime(time()));
print STDOUT "$Script start at:[$Time_Start]\n";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($ftab,$anno);
GetOptions(
				"help|h|?" =>\&USAGE,
				"tab:s"=>\$ftab,
				"anno:s"=>\$anno,
				) or &USAGE;
&USAGE unless ($ftab  and $anno);
$ftab = Cwd::realpath($ftab);

my %annotation;
open I,"$anno";
while(<I>){
	chomp;
	next if /^\s*$/ or /^\#/;
	my @tmp=split/\t/;
	$annotation{$tmp[0]}=$tmp[1];
	$annotation{$tmp[0]}="" unless $tmp[1] ;
}
close I;
open TAB,"$ftab" or die "cannot open $ftab \n";
open OUTKO,">$ftab.txt" or die "cannot open $ftab.txt ,$!\n";
while (<TAB>) {
	chomp;
	next if (/^\s*$/);
	my @temp = split /\t/, $_;
	if (defined $annotation{$temp[-1]}) {
		$temp[4]=$temp[-1];
		$temp[-1]=$annotation{$temp[-1]};
		print OUTKO join("\t",@temp),"\n";
	}	

}
close TAB;
close OUTKO;
system("mv $ftab.txt $ftab");

#######################################################################################
my $Time_End = sub_format_datetime(localtime(time()));
print STDOUT "\n$Script Done at: [$Time_End]\t\tTotal elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub sub_format_datetime {   #Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub USAGE {
	my $usage=<<"USAGE";
	Program:$Script
	Version:$version	[2013/3/26]
	Contact:Sun Huaiyu <sunhy\@biomarker.com.cn>
	Description:	Extract KEGG Pathway and KOs file From Blast tab file;

	Usage:
		-tab         The blast_tab Out of seq with KEGG           must be given
		-anno        The kobas anno                               must be given
		-h          Help document
USAGE
	print $usage;
	exit;
}
