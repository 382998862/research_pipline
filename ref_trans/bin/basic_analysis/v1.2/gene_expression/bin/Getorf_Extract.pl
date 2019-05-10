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
use File::Basename qw(basename dirname);

##############################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作；
my %opts;
GetOptions(\%opts,"fa=s","od=s","h");
if (!defined($opts{fa}) || !defined($opts{od}) || defined($opts{h})) {
	&help;
	exit;
}

#######################

my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));

######################
my $notename=`hostname`;chomp $notename;
my $seq_fa=&ABSOLUTE_DIR ($opts{fa});
my $name=basename ($seq_fa);
$name=~m/(.*)\.fa/;
my $key=$1;
`mkdir $opts{od} ` if(!-d $opts{od});
my $odir=&ABSOLUTE_DIR ($opts{od});

chdir "$odir";
print "getorf -find 1 $seq_fa $key.orf >/dev/null 2>&1 \n";
system "getorf -find 1 $seq_fa $key.orf >/dev/null 2>&1 ";

my $orf=glob "$odir/*.orf";

open IN1,"$seq_fa" || die;
$/='>';
<IN1>;
my %seq_length;
my %seq_line;
while (<IN1>) {
	chomp;
	next if (/^$/ || /^\#/);
	my ($head,$seq)=split/\n+/,$_,2;
	my $id=(split/\s+/,$head)[0];
	$seq=~s/\s+//g;
	$seq_line{$id}=$seq;
	my $length=length($seq);
	$seq_length{$id}=$length;
}

close IN1;

open IN2,"$orf" || die;
<IN2>;
my %cds_pos;
my %pep;
my %cds_strand;
while (<IN2>) {
	chomp;
	next if (/^$/ || /^\#/);
	my ($head,$seq)=split /\n+/,$_,2;
	$seq=~s/\s+//g;
	my $length=length($seq);
	my @tmp=split/\s+/,$head;
	$tmp[0]=~/(.*)_\d+/;
	my $id=$1;
	$tmp[1]=~/\[(\d+)/;
	my $start_pos=$1;
	$tmp[3]=~/(\d+)\]/;
	my $end_pos=$1;
	my $strand="+";
	if ($head=~/REVERSE\sSENSE/) {
		$strand="-";
	}
	if (!defined $pep{$id}) {
		$pep{$id}=$seq;
		$cds_pos{$id}{$start_pos}=$end_pos;
		$cds_strand{$id}=$strand;
	}
	elsif (defined $pep{$id} && $length>length ($pep{$id})) {
		$pep{$id}=$seq;
		foreach my $s (keys %{$cds_pos{$id}}) {
			delete $cds_pos{$id};
		}
		$cds_pos{$id}{$start_pos}=$end_pos;
		$cds_strand{$id}=$strand;
	}
}
$/="\n";
close IN2;

open OUT1,">$odir/$key.pep.fa" || die;
open OUT2,">$odir/$key.cds.fa" || die;
open OUT3,">$odir/$key.cds_pep.stat.xls" || die;
foreach my $id (keys %pep) {
	print OUT1 ">$id\n$pep{$id}\n";
	my $cds;
	foreach my $start (keys %{$cds_pos{$id}}) {
		if ($cds_strand{$id} eq "+") {
			print OUT3 ">$id\tlength=$seq_length{$id}\tstrand=\'$cds_strand{$id}\'\tstart=$start\tend=$cds_pos{$id}{$start}\n";
			$cds=substr($seq_line{$id},$start-1,$cds_pos{$id}{$start}-$start+1);
		}
		else {
			my $cds_seq=substr($seq_line{$id},$cds_pos{$id}{$start}-1,$start-$cds_pos{$id}{$start}+1);
			print OUT3 ">$id\tlength=$seq_length{$id}\tstrand=\'$cds_strand{$id}\'\tstart=$cds_pos{$id}{$start}\tend=$start\n";
			$cds=reverse $cds_seq;
			$cds=~tr/ATCGatcg/TAGCtagc/;
		}
	}
	print OUT2 ">$id\n$cds\n";
	print OUT3 "$cds\n$pep{$id}\n";
}

close OUT1;
close OUT2;
close OUT3;

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);

###########subs
sub Runtime
{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
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

sub help{
	print <<"	Usage End.";
	Description:
		Function : use Getorf predict transcript cDNA and the Pep sequence;
		Version  : $ver.
		Writer   : mengf <mengf\@biomarker.com.cn>
		Usage    :
		-fa
		   Unigene seq fa;
		-od
		   pep seq Out dir ;
		-h
		    Help document
		Attention The Program Cat't Run Backplat;
	Usage End.

	exit;
}

