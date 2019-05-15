#!/usr/bin/env perl
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use threads;
use newPerlBase;
my $Title="DE_target";
my $version="v1.0"; 

my ($cfg,$idir,$odir,$step,$only_step,$list,$ceRNA);
GetOptions(
    "cfg:s"	=>\$cfg,
    "DEG:s"	=>\$idir,
    "od:s"	=>\$odir,
    "id_name:s"	=>\$list,
    "ceRNA:s"	=>\$ceRNA,
    "step:s"	=>\$step,
    "only_step:s"	=>\$only_step,
    "help|h"	=>\&USAGE,
    ) or &USAGE;
&USAGE unless ($cfg and $idir);
$odir||="./";
`mkdir $odir` unless (-d "$odir");
$idir= abs_path($idir);
$cfg=abs_path($cfg);
$odir=abs_path($odir);
$ceRNA=abs_path($ceRNA);
$list=abs_path($list) if(defined $list);

my %type;my %config;
open(CFG,"$cfg") or die $!;
while(<CFG>){
	chomp;
	next if(/^#|^$|^\s$|^Diff/);
	my @tmp = split /\s+/,$_;
	if(/^Sample/){
		if($tmp[1] eq "ID"){
			for(my $i=2;$i<@tmp;$i++){
				$type{$tmp[$i]}=1;
			}
		}
		next;
	}
	$config{$tmp[0]}=$tmp[1];	
}
close(CFG);

###
my ($diff_perl,$venn_perl);
if(exists $type{"lncRNA"} && exists $type{"circRNA"} && exists $type{"sRNA"}){
	$diff_perl = "get_diffRNA_lcs.pl";
	$venn_perl = "get_veen_lcs.pl";
}elsif(exists $type{"lncRNA"} && exists $type{"circRNA"} && !exists $type{"sRNA"}){
	$diff_perl = "get_diffRNA_lc.pl";
	$venn_perl = "get_veen_lc.pl";
}elsif(!exists $type{"lncRNA"} && exists $type{"circRNA"} && exists $type{"sRNA"} && exists $type{"mRNA"}){
	$diff_perl = "get_diffRNA_csm.pl";
	$venn_perl = "get_veen_csm.pl";
}elsif(!exists $type{"lncRNA"} && exists $type{"circRNA"} && !exists $type{"sRNA"} && exists $type{"mRNA"}){
	$diff_perl = "get_diffRNA_cm.pl";
	$venn_perl = "get_veen_cm.pl";
}elsif(!exists $type{"lncRNA"} && exists $type{"circRNA"} && exists $type{"sRNA"} && !exists $type{"mRNA"}){
	$diff_perl = "get_diffRNA_cs.pl";
	$venn_perl = "get_veen_cs.pl";
}elsif(!exists $type{"lncRNA"} && !exists $type{"circRNA"} && exists $type{"sRNA"} && exists $type{"mRNA"}){
	$diff_perl = "get_diffRNA_ms.pl";
	$venn_perl = "get_veen_ms.pl";
}else{
	print "Check Your data type!\n";
	die;
}


$step ||= 1;
`mkdir $odir/work_sh` unless (-d "$odir/work_sh");
my $cmd;
if($step==1) {
	open(SH,">$odir/work_sh/s1.targetRNA_relationship.sh");
	$cmd = "$diff_perl -in $idir -cfg $cfg -od $odir ";
	my $medical = (split /\//,(dirname $config{Ref_seq}))[-2];
	if($medical =~/Homo_sapiens|Mus_musculus|Rattus_norvegicus/i){
		if(defined $list){
			$cmd .= "-name2list $list";
		}
	}
	print SH "$cmd ";
	close(SH);
	system("sh $odir/work_sh/s1.targetRNA_relationship.sh > $odir/work_sh/s1.targetRNA_relationship.sh.log");
	$step++ unless(defined $only_step);
}

if($step==2) {
	open(SH,">$odir/work_sh/s2.veen.sh");
	$cmd = "$venn_perl -in $idir -cfg $cfg -od $odir ";
	print SH "$cmd\n";
        close(SH);
	system("sh $odir/work_sh/s2.veen.sh > $odir/work_sh/s2.veen.sh.log");
	$step++ unless(defined $only_step);
}

if($step==3) {
        open(SH,">$odir/work_sh/s3.ce-cytoscape.sh");
        $cmd = "get_ceRNA_cytoscape.pl -in $idir -ce $ceRNA -od $odir -cfg $cfg";
        print SH "$cmd\n";
        close(SH);
	system("sh $odir/work_sh/s3.ce-cytoscape.sh > $odir/work_sh/s3.ce-cytoscape.sh.log");
}

########################################Sub Function########################################################
sub log_current_time {
     # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}

#############################################################################################################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: niepy <niepy\@biomarker.com.cn> 
      Date: 2018-06-041
     Usage:
            -cfg	<FILE>	config, rawdata path
            -DEG		<DIR>	analysis/DEG_Analysis
            -od		<DIR>	Combine/cytoscape
	    -id_name	<file>	gene2name.list
	    -ceRNA	<file>	ceRNA analyis result
            -step	<INT>	step to start from or steps to run
                          1     targetRNA_relationship
                          2     venn
                          3     ce-cytoscape
            -only_step		step-by-step analysis, start from the step specified by -step or not, only run the steps specified by step
            -h			help documents

   Example:
            perl $Script -cfg cfg -DEG DEG_Analysis -od Combine/cytoscape -id_name gene2name -ceRNA ceRNA.xls -step 1

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
