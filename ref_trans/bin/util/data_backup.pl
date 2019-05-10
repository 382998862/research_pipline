#!/usr/bin/perl -w
#
# Copyright (c) BMK 2012
# Writer:         mengf <mengf@biomarker.com.cn>
# Program Date:   2012.
# Modifier:       mengf <mengf@biomarker.com.cn>
# Last Modified:  2012.
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $program_name=basename($0);

my $ver="1.0";
############################################
my %opts;
GetOptions(\%opts,"id=s","od=s","h");
if (!defined($opts{id})||!defined($opts{od})||defined($opts{h})) {
	&help();
}
###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\n[$Time_Start] $Script start ... \n\n";


###############
#my $cfg=&ABSOLUTE_DIR($opts{cfg});
my $id = &ABSOLUTE_DIR($opts{id});
&MKDIR($opts{od});
my $od=&ABSOLUTE_DIR($opts{od});
my $basic="$od/Basic_Analysis" ;
&MKDIR($basic) unless (-d "$od/Basic_Analysis");

my $mapstat="$od/Basic_Analysis/Map_Stat" ;
&MKDIR($mapstat) unless (-d "$od/Basic_Analysis/Map_Stat");

my $final_track="$od/Basic_Analysis/final_track" ;
&MKDIR($final_track) unless (-d "$od/Basic_Analysis/final_track");


my $DEG_PPI="$od/DEG_PPI" ;
&MKDIR($DEG_PPI) unless (-d "$od/DEG_PPI");


#open LOG,"$od/log.txt" || die;
###############
#my %para;
#&para_load($cfg,\%para);

############### Extract Assembly dir

if (-d "$id/Structure_and_Expression") {
	system "cp $id/Structure_and_Expression/Cuffmerge/merged.gtf   $basic  ";
	system "cp -r $id/Structure_and_Expression/Cufflinks/   $basic  ";
	system "cp $id/Structure_and_Expression/Map_Stat/*.list   $mapstat  ";
	system "cp $id/Structure_and_Expression/Map_Stat/*.xls   $mapstat  ";
	system "cp $id/Structure_and_Expression/Map_Stat/*.stat   $mapstat  ";
	system "cp $id/Structure_and_Expression/geneExpression/final_track/*.list   $final_track  ";
}

if (-d "$id/work_sh") {
	system "cp -r $id/work_sh   $od  ";
}

if (-d "$id/Structure_and_Expression/DEG_Analysis") {
    system "cp -r $id/Structure_and_Expression/DEG_Analysis/DEG_PPI/used_*   $DEG_PPI  ";
}
################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $program_name Time :[$Time_End]\n\n";
###############Subs
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

sub para_load {
	my ($file,$para)=@_;
	open IN,$file||die "$!";
	while (<IN>) {
		chomp;
		s/\r+//g;
		next if(/^$/||/^\#/);
		my ($key,$info)=(split/\s+/,$_)[0,1];
		if(!$key){print "$_\n";die;}
		$para->{$key}=$info;
	}
	close IN;
}
sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub Runtime{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "\nTotal elapsed time: ${t}s\n";
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub help{
	print << "	Usage End.";
	Description: Extract Have_Ref Transcriptome Reaults for Html Process;
	version:$ver
	Usage:
		--id  <STR>   input dir, analysis output directory   force
		--od  <STR>   result output dir                      force
		--h           help
	Usage End.
		exit;
}
