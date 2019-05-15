#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $path = $Bin;
$path = substr($path,0,index($path,"bin"));
my %config=%{readconf("$path/project.cfg")};
my $Title=$config{Title};												#流程的名称，必填
my $version=$config{version};

my $notename=`hostname`;chomp $notename;

#######################
if (@ARGV!=3) {
	print "\n\tFunction: tackle graphgff Use SpliceGrapher soft\n\n";
	print "\tVersion: $version\n";
	print "\tUsage: <alt_graphgff_dir> <out_dir> <Sample>\n\n";
	exit;
}

#########################
my $python = $config{python}; # 2014-12-18 ~ 
#my $view_splicegraphs = "/share/nas2/genome/biosoft/SpliceGrapher/current/view_splicegraphs.py";
my $view_splicegraphs = $config{view_splicegraphs}; # 2014-12-18 ~ 



$ARGV[0]=&AbsolutePath($ARGV[0]);
&MKDIR($ARGV[1]);
$ARGV[1]=&AbsolutePath($ARGV[1]);

my @chrs = glob "$ARGV[0]/*";

&MKDIR ("$ARGV[1]/work_sh");

my $Gff_dir="$ARGV[1]/gff";
&MKDIR ($Gff_dir);
my $pdf_dir="$ARGV[1]/pdf";
&MKDIR ($pdf_dir);

if (-f "$ARGV[1]/$ARGV[2].SpliceGrapher.gff") {
	system "rm $ARGV[1]/$ARGV[2].SpliceGrapher.gff";
}

open OUT,">$ARGV[1]/$ARGV[2].SpliceGrapher.gff" || die $!;
foreach my $chr_dir (@chrs) {
	next if (!-d $chr_dir);
	my $chr = basename ($chr_dir);
	my @gffs = glob "$chr_dir/*.gff";

	&MKDIR ("$Gff_dir/$chr");
	&MKDIR ("$pdf_dir/$chr");

	open SH,">$ARGV[1]/work_sh/$ARGV[2].$chr.SpliceGrapher.sh" || die $!;

	my $i = 0;
	my $sh;
	foreach my $gff (@gffs) {
		$gff =~ /.*\/(.*)\.gff/;
		my $gene = $1;
		my $Alt_Form = `less -S $gff | grep AltForm |wc -l `;
		if ($Alt_Form != 0) {
			system "cat $gff >>$ARGV[1]/$ARGV[2].SpliceGrapher.gff";
			system "cp $gff $Gff_dir/$chr";
			if ($i < 10) {
				$sh .= "$python $view_splicegraphs $gff -o $ARGV[1]/pdf/$chr/$gene.splicegraph.pdf && ";
				$i++;
			}
			else {
				print SH "$sh";
				print SH "$python $view_splicegraphs $gff -o $ARGV[1]/pdf/$chr/$gene.splicegraph.pdf \n";
				$sh =();
				$i = 0;
			}
		}
		else {
			next;
		}
	}
	close SH;
}
close OUT;

my @sh_files = glob "$ARGV[1]/work_sh/$ARGV[2]*.SpliceGrapher.sh";
foreach my $sh_file (@sh_files) {
	system "sh $sh_file";
}


############ subs
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	if ($line<=1000) {
		system "$config{qsub} --queue $queue --maxproc $cpu --resource vf=$vf --reqsub  $shell";
	}
	if ($line>1000) {
		my @div=glob "$shell.div*";
		foreach (@div) {
			if (-e $_) {
				system "rm $_";
			}
		}
		@div=();
		my $div_index=1;
		my $line_num=0;
		open IN,"$shell" || die;
		while (<IN>) {
			chomp;
			open OUT,">>$shell.div.$div_index" || die;
			if ($line_num<1000) {
				print OUT "$_\n";
				$line_num++;
			}
			else {
				print OUT "$_\n";
				$div_index++;
				$line_num=0;
				close OUT;
			}
		}
		if ($line_num!=0) {
			close OUT;
		}
		@div=glob "$shell.div*";
		foreach my $div_file (@div) {
			system "$config{qsub} --queue $queue --maxproc $cpu --resource vf=$vf --reqsub $div_file";
		}
	}
}

sub Check_qsub_error {#
	# Check The qsub process if error happend 
	my $sh=shift;
	my @Check_file=glob "$sh*.qsub/*.Check";
	my @sh_file=glob "$sh*.qsub/*.sh";

	if ($#sh_file!=$#Check_file) {
		print "Their Some Error Happend in $sh qsub, Please Check..\n";
		die;
	}
	else {
		print "$sh qsub is Done!\n";
	}
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub AbsolutePath {		#获取指定目录或文件的决定路径
	my $input = shift;
	my $return;
	if ($input =~ /\.\/\s*$/) {
		chomp($return = `pwd`);
		return $return;
	}
	my $dir=dirname($input);
	my $file=basename($input);
	chomp(my $pwd = `pwd`);
	chdir($dir);
	chomp($return = `pwd`);
	$return .= "\/".$file;
	chdir($pwd);
	return $return;
}
