#!/usr/bin/perl -w
#
# Copyright (c) BMK 2011
# Writer:         He hua <heh@biomarker.com.cn>
# Program Date:   2011.
# Modifier:       He hua <heh@biomarker.com.cn>
# Last Modified:  2011.
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(sum);
my $program_name=basename($0);

my $ver="1.0";
############################################
my %opts;
my ($gtf,$out,$isoform_len,$isoform_num,$isoform_cov,$no_sense);
GetOptions(
				"help|?" =>\&USAGE,
				"gtf:s"=>\$gtf,
				"l:s"=> \$isoform_len,
				"n:s"=> \$isoform_num,
				"c:s"=> \$isoform_cov,
				"out:s"=>\$out,
				"no_sense:s"=>\$no_sense,
				) or &USAGE;
&USAGE unless ($gtf and $out);

$isoform_len= $isoform_len || 200;
$isoform_num=$isoform_num ||  2;
$isoform_cov||= 2;
my %lnc_type = (
	"e" => "pre-mRNA_lncRNA",
	"o" => "exonic_overlap_lncRNA",
	"u" => "lincRNA",
        "i" => "intronic RNA",
        "x" => "anti-sense RNA",
);
if (defined $no_sense){
	delete$lnc_type{"o"};
	delete$lnc_type{"e"};
}
###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";
###############
my %trans;
#print Dumper(\%trans);die;

open GTF,"$gtf"||die "$!";
while (<GTF>) {
	chomp;
	next if (/\#/||/^$/) ;
	my ($trans_id,$type);
	my $line = $_;
	my @lines = split ("\t",$line);
	if ($lines[2]=~/transcript/){
		$line =~/transcript_id\s\"([^\"]+)\";/;
		$trans_id = $1;
		 $line =~/class_code\s\"([^\"]+)\";/;
		$type =$1;
		$trans{$trans_id}{'type'} = $type;
	}
	if ($lines[2]=~/exon/){
		$line =~/transcript_id\s\"([^\"]+)\";/;
		$trans_id = $1;
		my $exon_len = abs($lines[4]-$lines[3]);
		push (@{$trans{$trans_id}{'exon'}},$exon_len);
	}
	
	
	
	
}
close GTF;
#print Dumper(\%trans);
foreach my $trans_id (keys (%trans)) {
	$trans{$trans_id}{'length'} = sum @{$trans{$trans_id}{'exon'}};
	#$trans{$trans_id}{'exon_num'} = $trans{$trans_id}{'exon'};
	my $num=$trans{$trans_id}{'exon'};
        $trans{$trans_id}{'exon_num'} =@$num;
	my $trans_type = $trans{$trans_id}{'type'};
	if ($trans{$trans_id}{'length'} < $isoform_len) {
		delete $trans{$trans_id};
	}
	elsif($trans{$trans_id}{'exon_num'} < $isoform_num) {
		delete $trans{$trans_id};
	}
	elsif (!exists $lnc_type{$trans_type}) {
		 delete $trans{$trans_id};
	}
}
my $lnc_num = keys(%trans);
print "There are $lnc_num transcription in library\n";
open GTF,"<$gtf"||die "$!";
open OUT,">$out"||die $!;
while (<GTF>) {
	chomp;
	next if (/\#/||/^$/) ;
	my $line = $_;
	$line =~/transcript_id\s\"([^\"]+)\";/;
	my $trans_id = $1;
	#my ($class) = $line =~/class_code "(\S+)";/;
	if ( exists $trans{$trans_id} ) {
		print OUT "$line\n";
	}
}
close GTF;
close OUT ;

################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $program_name Time :[$Time_End]\n\n";

###############Subs
sub ABSOLUTE_DIR
{
        my ($in,$cur_dir)=@_;
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
                warn "Warning just for file and dir\n";
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
sub Runtime
{ # &Runtime($BEGIN);
        my ($t1)=@_;
        my $t=time()-$t1;
        print "\nTotal elapsed time: ${t}s\n";
}
sub USAGE{
	print << "	Usage End.";
Description:The process get cds or Intro sequence from genome fasta file base on gff file.
version:$ver
Usage:
	-gtf	<STR>	gtf file                         must be given;
	-l	<NUM>	transcription length;
	-n	<NUM>	Exon Num;
	-c	<NUM>	transcription coverage;
	-no_sense	filter sense lncRNA
	-out	<STR>	Out file;

	Usage End.
		exit;
}
