#!/usr/bin/perl -w
#
# Copyright (c) BMK 2012
# Writer:         He hua <heh@biomarker.com.cn>
# Program Date:   2012.
# Modifier:       He hua <heh@biomarker.com.cn>
# Last Modified:  2012-4-20.
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $program_name=basename($0);

my $ver="1.0";
############################################
my %opts;
GetOptions(\%opts,"ref=s","od=s","h");
if (!defined($opts{ref})||!defined($opts{od})||defined($opts{h})) {
	&help();
}
###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";
###############
my $ref=path($opts{ref});
mkdir "$opts{od}" if(!-d $opts{od});
my $outdir=path($opts{od});
my $outkey=basename($ref);
$/=">";
open L,">$outdir/$outkey.len"||die;
open GC,">$outdir/$outkey.GC"||die;
print L "#Chrom\tAGCTN\tAGCT\tN\n";
print GC "#Chrom\tPercent\n";
my $Total_len=0;my $Total_N=0;my $Total_GC=0;
open R,$ref||die;
while (<R>) {
	chomp;
	my $GC_L=0;
	my $N_L=0;
	next if(/^$/);
	my ($chr,$seq)=split/\n+/,$_,2;
	$chr=(split/\s+/,$chr)[0];
	print "$chr\n";
	$seq=~s/\s+//g;
	$GC_L=$seq;
#	print "$GC_L\n";
	$N_L=$seq;
	$GC_L=~s/A//ig;
	$GC_L=~s/N//ig;
	$GC_L=~s/T//ig;
#	print "A\n";
	$N_L=~s/A//ig;
	$N_L=~s/T//ig;
	$N_L=~s/C//ig;
	$N_L=~s/G//ig;
#	print "B\n";
	my $GC=length$GC_L;
	my $N=length$N_L;
	my $len=length$seq;
#	print "A:$GC\tB:$N\t$len\n";die;
	print L "$chr\t$len\t",$len-$N,"\t$N\n";
	printf GC "%s\t%.2f\n",$chr,100*$GC/($len-$N); 
	$Total_len+=$len;
	$Total_N+=$N;
	$Total_GC+=$GC;
}
close R;
$/="\n";
print L "#Total_Len\t$Total_len\t",$Total_len-$Total_N,"\t$Total_N\n";
printf GC "%s\t%.2f\n","Total_GC",100*$Total_GC/($Total_len-$Total_N);
close L;
close GC;
################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $program_name Time :[$Time_End]\n\n";
###############Subs
sub path{
	my ($in)=@_;
	my $return;
	my $cur_dir=`pwd`;
	chomp($cur_dir);
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "$in:Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
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
sub help{
	print << "	Usage End.";
	Description:
		Count fasta file per chromsome length and GC percent.
	version:$ver
	Usage:
		
		-ref     reference fasta file     must be given;
		-od      outdir                   must be given;

	Usage End.
		exit;
}
