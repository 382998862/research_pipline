#!/usr/bin/perl -w
# 
#Copyright (c) BMK 2012 
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2012 
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2012 
my $version="1.0.0";
my $BEGIN=time();

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $programe_dir=basename($0);
my $path=dirname($0);
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET=1;
use newPerlBase;
my %config=%{readconf("$Bin/../db_file.cfg")}; 

#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info
# ------------------------------------------------------------------
my ($id,$Blast_cpu,$Blast_e,$Database,$odir,$HELP,$queues);
my @anno;

GetOptions(
		"id:s"=>\$id,
		"Database:s"=>\$Database,
		"od:s"=>\$odir,
		"Blast_cpu:s"=>\$Blast_cpu,
		"Blast_e:s"=>\$Blast_e,
		"queue:s"=>\$queues,
		"help"=>\$HELP
	) or &USAGE;

&USAGE if (!defined $id || !defined $odir || !defined $Database || $HELP) ;


###------------------软件路径----------------------------###
# all the program are in this dir                    
chomp $Bin;
#my $blastall = "/share/nas2/genome/biosoft/blast/2.2.26/bin/blastall";  # 2015-02-04 ~ 
my $blastx = "$config{blast}/blastx";  # 2015-09-01 ~ 
my $blastn = "$config{blast}/blastn"; # 2015-09-01 
#=================== 一些参数 ================================
&MKDIR($odir);
$id=&ABSOLUTE_DIR($id);
$odir=&ABSOLUTE_DIR($odir);
$Database=&ABSOLUTE_DIR($Database);


$Blast_cpu||=50;
$Blast_e||=1e-5;
$queues||="middle.q";
my $notename=`hostname`;chomp $notename;

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\n$programe_dir Start Time :[$Time_Start]\n\n";


my $Q_name=(glob "$id/*.fa")[0];
$Q_name=basename $Q_name;
if ($Q_name=~/(.+)\.\d+\.fa$/) {
	$Q_name=$1;
}
else {print "Your file name is Wrong!\n";die;}

&MKDIR("$odir/Result");
my $Result_dir="$odir/Result";
&MKDIR("$odir/Cog_Dir");
&MKDIR("$odir/work_sh");
&MKDIR("$odir/02.gene-annotation");
my $Tab_dir="$odir/02.gene-annotation";
&MKDIR("$odir/work_sh");
&MKDIR("$odir/work_sh/Cog_sh");
my $sh_dir="$odir/work_sh/Cog_sh";

my $blast_shell_file = "$odir/work_sh/Cog_sh/Cog.blast.sh";
my @subfiles;
@subfiles = glob("$id/*.fa");

############creat shell file
open OUT,">$blast_shell_file" || die "fail $blast_shell_file";
foreach my $subfile (@subfiles) {
	my $name=basename $subfile;
	print OUT "$blastx -task blastx-fast -num_descriptions 100 -num_alignments 100 -evalue $Blast_e -db  $Database -query  $subfile -num_threads 2 -out  $odir/Cog_Dir/$name.Cog.blast &&\n";
}
close OUT;

####################run the shell file
&qsubOrDie("$blast_shell_file",$queues,$Blast_cpu,"6G");

`perl $Bin/bin/blast_format.pl -id $odir/Cog_Dir -od $Tab_dir -shdir $sh_dir  -middir $odir/mid --cog -queue $queues`;


`perl $Bin/bin/cog_parser.pl $Tab_dir/$Q_name.Cog.blast.tab.best `;
system "cp $Tab_dir/$Q_name.Cog.class $Result_dir/$Q_name.Cog_class.txt";
#`perl $Bin/bin/CogFunClassDrawer.pl -i $Tab_dir/$Q_name.Cog.class -o $Tab_dir/$Q_name.Cog.cluster.svg -png`;


###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\n$programe_dir End Time :[$Time_End]\n\n";
&Runtime($BEGIN);

#====================================================================================================================
#  +------------------+
#  |   subprogram     |
#  +------------------+

sub LOAD_PARA
{
	my $para_file= shift;
	my $para= shift;

	my $error_status = 0;
	open IN,$para_file || die "fail open: $para_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if (/^$/ or /^\#/) ;
		my ($para_key,$para_value) = split(/\s+/,$_);
		$para->{$para_key} = $para_value;
		if (!defined $para_value) {
			warn "Non-exist: $para_key\n";
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
}

sub runcmd
{
	my ($program,$cmd)=@_;
	open (SH, ">>work.sh") ||  die "Can't write work: $!\n";
	print SH "$cmd \n";
	system($cmd) && LOGFILE(1,$program); 
	close SH;
}

sub LOGFILE
{
	my $flog=shift;
	my $program=shift;
	my $Time= sub_format_datetime(localtime(time()));
	print LOG "[$Time +0800]\ttask\t0\tstart\t$program\tStart to analysis......\n";
	if($flog==0){
		print LOG "[$Time +0800]\ttask\t0\tend\t$program\tDone.\n";
	}else{
		print LOG "[$Time +0800]\ttask\t0\terror\t$program\tAt least one $program in this section is in error.\n";
		close LOG;
		exit;
	}
}
close LOG;
sub ann
{
	@anno = split /\t/, $_;
	print OUT "$anno[0]\t$anno[4]\t$anno[13]\t$anno[8]\t$anno[12]\t$anno[15]\n";
}

sub LOAD_SEQ 
{
	my ($fa,$info) = @_;

	open IN,"$fa" || die $!;
	$/='>';
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r+$//;
		next if (/^$/ || /^\#/);
		my ($head,$seq)=split/\n+/,$_,2;
		my $id=(split/\s+/,$head)[0];
		$info->{$id}=$seq;
	}
	$/="\n";
	close IN;
}


sub CUTFA 
{
	my ($fa,$od,$cut,$name) = @_;

	&MKDIR("$od/$name.div");
	my %seq=%$fa;
	my @aa=sort(keys %seq);
	my $index=0;
	LAB: for (my $i=1;;) {
		my $num=0;
		open OUT,">$od/$name.div/$name.$i.fa" || die $!;
		for ($index..$#aa) {
			$index++;
			if ($num<$cut) {
				print OUT ">$aa[$_]\n$seq{$aa[$_]}\n";
				$num++;
			}
			if ($num>=$cut) {
				$num=0;
				$i++;
				close OUT;
				if ($index==$#aa+1) {
					last;
				}
				else {
					next LAB;
				}
			}
		}
		if ($num) {
			close OUT;
		}
		last;
	}
}


sub nr_pie_stat{#&nr_pie_stat($in,$out)
	my ($in,$out)=@_;
	my %H;
	open (IN,$in) or die $!;
	<IN>;
	while (<IN>) {
		chomp;
		my @A=split/\t/,$_;;
		my $name=$A[-1];
		$name=~s/\s*$//;
		next unless $name=~/\[([^\]]+)\]$/;
		$H{$1}=0 unless exists $H{$1};
		$H{$1}++ if exists $H{$1};
	}
	close (IN) ;

	open (OUT,">$out") or die $!;
	my $limit=keys %H;
	if ($limit<=10) {
		if ($limit == 1) {
			my ($key,$value) = each %H;
			my $str = &cut_str($key);
			print OUT "$str\t$value\t$key\n";
		}else{
			foreach my $key (sort {$H{$b}<=>$H{$a}} keys %H) {
				my $str = &cut_str($key);
				print OUT "$str\t$H{$key}\t$key\n";
			}
		}
	}
	else {
		my $n=0;
		my $other=0;
		foreach my $key (sort {$H{$b}<=>$H{$a}} keys %H) {
			$n++;
			if($n<10){
				my $str = &cut_str($key);
				print OUT "$str\t$H{$key}\t$key\n";
			}
		$other+=$H{$key} unless $n<10;
		}
	print OUT "Other\t$other\tOther species\n";
	}
	close OUT;
}
sub cut_str
{
	my $string = shift;
	my @str = split /\s+/,$string;
	if (@str > 2) {
		return "$str[0] $str[1]"
	}else{
		return $string;
	}
}

sub parse_config
{ # load config file
	my $config_file= shift;
	my $DataBase= shift;
	
	my $error_status = 0;
	open IN,$config_file || die "fail open: $config_file";
	while (<IN>) {
		chomp;
		s/\s+//;s/\s+$//;s/\r$//;
		next if(/$/ or /\#/);
		my ($software_name,$software_address) = split(/\s+/,$_);
		$DataBase->{$software_name} = $software_address;
		if (! -e $software_address){
			warn "Non-exist:  $software_name  $software_address\n"; 
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
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

sub USAGE 
{
	print <<"	Usage End.";
	Program:Blast && interproscan with COG_DataBase Annotate New Gene
	Version: $version
	Contact: zhang xuechuan <zhangxc\@biomarker.com.cn>

	Description:

      -id                    mRNA dir,forced
      -Database              Nr Database file,forced
      -od                    OUT DIR,forced
      -Blast_cpu             cpu number,default 50
      -Blast_e               e-value of cutoff,default 1e-5
      -queue                 the queue is used for qsub jobs 

      -help

	Usage End.
	exit;
}