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
GetOptions(\%opts,"cfg=s","id=s","h");
if (!defined($opts{cfg})||!defined($opts{id})||defined($opts{h})) {
	&help();
}
###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $program_name Time :[$Time_Start]\n\n";


###############
my $cfg=&ABSOLUTE_DIR($opts{cfg});
my $indir=&ABSOLUTE_DIR($opts{id});

my %para;
&para_load($cfg,\%para);

&MKDIR("$opts{id}/html");
my $html_dir="$opts{id}/html";

system "cp -r $Bin/bin/css $html_dir";
system "cp -r $Bin/bin/images $html_dir";
system "cp -r $Bin/bin/js $html_dir";

##############
# Make enrichment html
##################
my @DEG_Combinations=glob "$indir/DEG_Analysis/*";
foreach my $com (@DEG_Combinations) {
	next unless -d $com;
	next if $com=~/\/corr$/;
	next if $com=~/\/density$/;
	$com=~/.*\/(.*)/;
	my $k=$1;
	system "perl $Bin/bin/html_go.pl -i $indir/DEG_Analysis/$k/go_enrichment -o $indir/DEG_Analysis/$k/go_enrichment -k $k";
	system "perl $Bin/bin/html_pathway.pl -i $indir/DEG_Analysis/$k/pathway -o $indir/DEG_Analysis/$k/pathway -key $k";
}

##############
# Make html report
##################
system "perl $Bin/bin/rawdata_stat_html.pl -r $indir/rawdata -c $cfg -o $indir/html";
#system "perl $Bin/bin/gene_expression_html.pl -m $indir/Map_Stat -g $indir/geneExpression -c $cfg -o $indir/html";
system "perl $Bin/bin/gene_expression_html.pl -g $indir/geneExpression -c $cfg -o $indir/html";
system "perl $Bin/bin/Alt_splice_html.pl -d $indir/Alt_splice -c $cfg -o $indir/html";
system "cp $Bin/bin/Alt_splice_info.html $indir/html";

=a
if (defined $para{IDEG6}) {
	system "perl $Bin/bin/html_IDEG6_deg.pl -d $indir/DEG_Analysis -c $cfg -o $indir/html";
}
=cut

if (defined $para{DEG}) {
	system "perl $Bin/bin/html_deg.pl -d $indir/DEG_Analysis -c $cfg -o $indir/html";
}

#############
# Make index html
##################
system "perl $Bin/bin/html_index.pl -c $cfg -o $indir";


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
		-cfg  <STR>   config file   force
		-id   <STR>   outdir        force
		-h            help
	Usage End.
		exit;
}
