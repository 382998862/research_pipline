#!/usr/bin/perl -w
use strict;
use newPerlBase;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $Title = "fusion_circos_link";
my $version="1.0.1";
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($circos_cfg,$fusion_result,$odir,  $data_analyzer, $customer_service_exe,$test,$type);

GetOptions(
    "cfg:s" =>\$circos_cfg,
    "fr:s"  =>\$fusion_result,
    "od:s"   =>\$odir,
	"test"   =>\$test,
	"type:s" =>\$type,
    "help|h" =>\&USAGE,
    ) or &USAGE;
&USAGE unless ($circos_cfg and $fusion_result and $odir );

createLog($Title,$version,$$,"$odir/log/",$test);

# ------------------------------------------------------------------
# read configs and infile
# ------------------------------------------------------------------
stepStart(1,"circos plot"); 
system "mkdir -p $odir" unless (-d $odir);
$odir=&ABSOLUTE_DIR($odir);
$circos_cfg=&ABSOLUTE_DIR($circos_cfg);
my %config = %{readconf("$circos_cfg")};
my $link = "$odir/link.txt";
my $text = "$odir/text.txt";
my $report = "$odir/FusionReport.txt";
open(IN,"$fusion_result")||die "open $fusion_result failed!\n";
open (OUT,">$link")||die "open or creat $link failed!\n";
open (OUT2,">$text")||die "open or creat $text failed!\n";
open (OUT3,">$report") ||die "open or creat $report failed!\n";
while (<IN>) {
	chomp;
	my @ll = split /\t/,$_;
	if (/^FusionID/) {
		print OUT3 "$_\n";
	}
	if ( defined $ll[14]) {
	}else{
		print OUT "chr$ll[6]\t$ll[7]\t$ll[8]\tchr$ll[11]\t$ll[12]\t$ll[13]\n";
		print OUT2 "chr$ll[6]\t$ll[7]\t$ll[8]\t$ll[4]\nchr$ll[11]\t$ll[12]\t$ll[13]\t$ll[9]\n";
		print OUT3 "$_\n";
	}
	#print OUT "$ll[2]\t$ll[3]\t$ll[3]\t$ll[5]\t$ll[6]\t$ll[6]\n";
}

close IN;
close OUT;
if ($type eq "Human.B38") {
	`perl $Bin/circos_human_v1.pl --chr $Bin/karyotype.human.hg38.txt --circle $odir/link.txt --type link --circle $odir/text.txt --type text --od $odir `;
	print "perl $Bin/circos_human_v1.pl --chr $Bin/karyotype.human.hg38.txt --circle $odir/link.txt --type link --circle $odir/text.txt --type text --od $odir\n";
}elsif ($type eq "Human.hg19" ) {
	`perl $Bin/circos_human_v1.pl --chr $Bin/karyotype.human.hg19.txt --circle $odir/link.txt --type link --circle $odir/text.txt --type text --od $odir `;
	print "perl $Bin/circos_human_v1.pl --chr $Bin/karyotype.human.hg19.txt --circle $odir/link.txt --type link --circle $odir/text.txt --type text --od $odir\n";
}elsif ($type eq "Mouse.B38") {
	`perl $Bin/circos_mouse_v1.pl --chr $Bin/karyotype.mouse.mm10.txt --circle $odir/link.txt --type link --circle $odir/text.txt --type text --od $odir `;
	print "perl $Bin/circos_mouse_v1.pl --chr $Bin/karyotype.mouse.mm10.txt --circle $odir/link.txt --type link --circle $odir/text.txt --type text --od $odir\n";
}elsif ($type eq "Rat.B5.0") {
	`perl $Bin/circos_rat_v1.pl --chr $Bin/karyotype.rat.rn5.txt --circle $odir/link.txt --type link --circle $odir/text.txt --type text --od $odir `;
	print "perl $Bin/circos_rat_v1.pl --chr $Bin/karyotype.rat.rn5.txt --circle $odir/link.txt --type link --circle $odir/text.txt --type text --od $odir\n";
}elsif ($type eq "Rat.B6.0") {
        `perl $Bin/circos_rat_v1.pl --chr $Bin/karyotype.rat.rn6.txt --circle $odir/link.txt --type link --circle $odir/text.txt --type text --od $odir `;
        print "perl $Bin/circos_rat_v1.pl --chr $Bin/karyotype.rat.rn6.txt --circle $odir/link.txt --type link --circle $odir/text.txt --type text --od $odir\n";
}

#`perl $Bin/circos_v1.pl --chr $Bin/circos.txt --circle $odir/link.txt --type link --circle $odir/text.txt --type text --od $odir `;
#print "perl $Bin/circos_v1.pl --chr $Bin/circos.txt --circle $odir/link.txt --type link --od $odir\n";

######################### Data Assesss

stepTime(1);
totalTime();
# 
#######################################################################################

#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#############################################################################################################
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
#############################################################################################################
#############################################################################################################
sub step_cmd_process {
    my ($cmd, $step_hash, $step_n, $sh_dir) = @_;
    my $sh_file = "$sh_dir/$step_n.$step_hash->{$step_n}.sh";
    my $log_file = "$sh_file.log";
    my $flag = 0;
    my $start_time = time();
    &log_current_time("step$step_n. $step_hash->{$step_n} start ...");
    &log_current_time("CMD: $cmd");

    if (-e $sh_file) {
        system "cat $sh_file >> $sh_file.bak";
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    } else {
        open (SH, ">$sh_file") or die "$!: $sh_file\n";
        print SH "$cmd\n";
        close SH;
    }

    $flag = system("$cmd > $log_file");
    if ($flag != 0){
        log_current_time("Error: command failed: $cmd");
        exit(1);
    } else {
        my $escaped_time = (time()-$start_time)."s";
        &log_current_time("step$step_n. $step_hash->{$step_n} done, escaped time: $escaped_time.\n");
    }
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



#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $0
   Version: $version

     Usage:
            --cfg       <FILE>  circos software path
            --od        <DIR>   analysis output directory
            --fr        <file>  the fusion gene result
            --PL        <STR>   abbr. of Project Leader\'s name
            --CSE       <STR>   abbr. of Customer Service Executive\'s name
            --h                 help documents

   Example:
            perl link.pl --cfg path.cfg --fr result.txt --od Analysis/ 
----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}
