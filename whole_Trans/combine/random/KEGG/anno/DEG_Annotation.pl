#!/usr/bin/perl -w
# 
#Copyright (c) BMK 2011 
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2011 
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2011 
my $ver="1.0.0";
my $BEGIN=time();

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname);

############ 流程名称，版本，及固定配置文件
use lib "/share/nas2/genome/bmksoft/tool/newPerlBase/v1.0";
use newPerlBase;

##############################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作；
my %opts;
GetOptions(\%opts,"DEG=s","gene=s","id=s","od=s","key=s","h");
if (!defined($opts{DEG}) || !defined($opts{id}) || !defined($opts{od}) || !defined($opts{key}) || defined($opts{h})) {
	&help;
	exit;
}

#######################

my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";

######################
my $anno_dir=abs_path($opts{id});
&MKDIR($opts{od});
my $odir=abs_path($opts{od});
my $deg_list=abs_path($opts{DEG});
my $key=$opts{key};

my $cog_parser="$Bin/ann/CogFunClassDrawer.pl";

&MKDIR("$odir/Cog_Anno");

my %DEG_ID;
my %DEG_info;
my $deg_head;
open DEG,"$deg_list" || die;
while (<DEG>) {
	chomp;
	next if (/^$/);
	if (/^\#/) {
		$deg_head=$_;
		next;
	}
	my $id=(split/\s+/,$_)[0];
	$DEG_ID{$id}=1;
	$DEG_info{$id}=$_;
}
close DEG;

my $cog=glob "$anno_dir/*.Cog_class.txt";
if (-e $cog) {
	open COG,"$cog" || die "error in [$cog], $!";
	open OUT,">$odir/Cog_Anno/$key.Cog_class.txt" || die "error in [$odir/Cog_Anno/$key.Cog_class.txt], $!";
	while (<COG>) {
		chomp;
		next if (/^$/) ;
		if (/^\#/) {
			print OUT "$_\n";
			next;
		}
		my ($id,$anno)=split /\s+/,$_,2;
        my $lc_id=lc($id);  #转换成小写字母
        my $uc_id=uc($id);  #转换成大写字母
		if (exists $DEG_ID{$id} || exists $DEG_ID{$lc_id} || exists $DEG_ID{$uc_id}) {
			print OUT "$id\t$anno\n";
		}
	}
	chdir "$odir/Cog_Anno/";
    print "cd $odir/Cog_Anno/\n";
    
    runOrDie("perl $cog_parser -i $key.Cog_class.txt -o $key.Cog.classfy.svg -png");
    if ( !-f "$key.Cog.classfy.plot.log") {  
        runOrDie("Rscript $Bin/ann/cog_anno_plot.r $key.Cog_class.txt.stat $key.Cog.classfy.png >$key.Cog.classfy.plot.log 2>&1");
        runOrDie("rm $key.Cog.classfy.svg");
    }
    else {
        print "skip CMD:[Rscript $Bin/ann/cog_anno_plot.r $key.Cog_class.txt.stat $key.Cog.classfy.png >$key.Cog.classfy.plot.log]\n";
    }
    
	close COG;
	close OUT;
}

my $anno_inte="$anno_dir/Integrated_Function.annotation.xls";

if (-e $anno_inte) {
	open IN,"$anno_inte" || die;
	open OUT,">$odir/$key.annotation.xls" || die;
	#print OUT "$deg_head\t";
	print OUT "gene_ID\t";
	while (<IN>) {
		chomp;
		next if (/^$/);
		if (/^\#/) {
			$_=~s/^\#\S+\s+//;
			print OUT "$_\n";
			next;
		}
		my ($gene_id,$anno)=split/\t+/,$_,2;
        my $lc_id=lc($gene_id);  #转换成小写字母
        my $uc_id=uc($gene_id);  #转换成大写字母
		if (exists $DEG_ID{$gene_id}) {
			my@tmp=split(/\s+/,$DEG_info{$gene_id});
			print OUT "$tmp[0]\t$anno\n";
		}
        elsif (exists $DEG_ID{$lc_id} ){
            my@tmp=split(/\s+/,$DEG_info{$lc_id});
			print OUT "$tmp[0]\t$anno\n";
        }
        elsif (exists $DEG_ID{$uc_id} ){
            my@tmp=split(/\s+/,$DEG_info{$uc_id});
			print OUT "$tmp[0]\t$anno\n";
        }
	}
	close OUT;
}

my @annos=glob "$anno_dir/*.anno.txt";

foreach my $anno_file (@annos) {
	next if ($anno_file=~/GO/);
	&extract_anno(\%DEG_ID,$anno_file,$odir);
}

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";
&Runtime($BEGIN);

###########subs
sub extract_anno {#
	my $DEG_list = shift @_;
	my $anno = shift @_;
	my $od = shift @_;

	my %list = %$DEG_list;
	my $database_name;
	if ($anno =~ /.*\.(.+)\.anno.txt/) {
		$database_name=$1;
	}
	open IN,"$anno" || die $!;
	open OUT,">$od/$key.$database_name.anno.txt" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r+$//;
		next if (/^$/);
		if (/^\#/) {
			print OUT "$_\n";
			next;
		}
		my $gene_id=(split /\s+/,$_)[0];
        my $lc_id=lc($gene_id);  #转换成小写字母
        my $uc_id=uc($gene_id);  #转换成大写字母
		if (exists $list{$gene_id} || exists $list{$lc_id} || exists $list{$uc_id}) {
			print OUT "$_\n";
		}
	}
	close IN;
	close OUT;
}


sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n";
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

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

sub help
{
	print <<"	Usage End.";
	Description:
		Function : Use DEG_Analysis Out to Extract Annotation file and Draw Pictures(Cog,Kegg..);
		Version  : $ver.
		Writer   : mengf <mengf\@biomarker.com.cn>
		Usage    :
		-DEG
		    DEG List , OUT put of DEG_Analysis;
		-id
		    Annotation OUT DIR (eg. /Result/);
		-od
		    DEG_Anno OUT dir (eg. /DEG_Anno);
		-key
		    The Prefix of OUT Anno files (eg. Munro_Q1_Q2_DEG ".Kegg.anno.txt");

	Usage End.

	exit;
}
