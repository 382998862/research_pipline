#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
# 
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($binding, $fusion, $ppi);

GetOptions(
        "help|?" =>\&USAGE,
        "bind:s"=>\$binding,
        "fusion:s"=>\$fusion,
        "ppi:s"=>\$ppi,
        ) or &USAGE;
&USAGE unless ($binding and $fusion and $ppi);

&timeLog("$Script start...");
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------
my %binding;

open (BIND, $binding) or die $!;

while (<BIND>) {
    chomp;
    next if (/^\s*$/ or /^#/);
    my ($bait_or_prey, $BD_or_AD1,$BD_or_AD2) = (split /\t/)[0,4,-1];
    my $subID;
    if($BD_or_AD1 !~ /gnl/){
    	 ($subID)=$BD_or_AD1=~/^(\S+)/;
    }
    else{
    	 ($subID)=$BD_or_AD2=~/^(\S+)/;
    }
    $binding{$subID}{$bait_or_prey} = 1;
}

close BIND;

# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------

open (PPI, ">$ppi") or die;
open (FUSION, $fusion) or die;
#print PPI "#Query_id1\tType\tQuery_id2\tSubject_id1\tSubject_id2\tScore\n";   # protein links
print PPI "#Query_id1\tType\tQuery_id2\tSubject_id1\tSubject_id2\tMode\tScore\n";   # protein actions

while (<FUSION>) {
    chomp;
    next if (/^\s*$/ or /^#/);
    my ($BD, $AD, $more_infos) = (split /\s+/, $_, 3);

    for my $k1 (keys %{$binding{$BD}} ) {
        for my $k2 (keys %{$binding{$AD}} ) {
            print PPI "$k1\tpp\t$k2\t$BD\t$AD\t$more_infos\n" if ($k1 ne $k2);
        }
    }
}

close FUSION;
close PPI;


#######################################################################################
my $elapse_time = (time()-$BEGIN_TIME)."s";
&timeLog("$Script done. Total elapsed time: $elapse_time.");
#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

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

#############################################################################################################
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

################################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
 ProgramName:   $Script
     Version:	$version
     Contact:	Simon Young <yangxh\@biomarker.com.cn> 
Program Date:	2012.07.02
      Modify:	
 Description:	This program is used to ......
       Usage:
        Options:
        --bind      <file>  input file, BD-bait bindings or AD-prey bindings, forced

        --fusion    <file>  input file, BD-AD fusions, forced

        --ppi       <file>  output file, bait-prey interacts, forced

        Example:
            perl $Script --bind test3/Maize.Unigene.blast.tab.best --fusion eukaryota.4577.protein.links.txt -ppi Maize.Unigene.ppi.xls 

USAGE
	print $usage;
	exit;
}
