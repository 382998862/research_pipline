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
#my $gene_model_to_splicegraph = "/share/nas2/genome/biosoft/SpliceGrapher/current/gene_model_to_splicegraph.py";
my $gene_model_to_splicegraph = $config{gene_model_to_splicegraph}; # 2014-12-18 ~ 

$ARGV[0]=&AbsolutePath($ARGV[0]);
&MKDIR($ARGV[1]);
$ARGV[1]=&AbsolutePath($ARGV[1]);

my @chrs = glob "$ARGV[0]/*";
system "mkdir $ARGV[1]/$ARGV[2]";
system "mkdir $ARGV[1]/work_sh";

open SH,">$ARGV[1]/work_sh/$ARGV[2].SpliceGraphgff.sh" || die $!;
open OUT,">$ARGV[1]/$ARGV[2].gene.list.txt" || die $!;
foreach my $chr_dir (@chrs) {
	next if (!-d $chr_dir);
	my $chr = basename ($chr_dir);
	my $chr_AS = "$ARGV[1]/$ARGV[2]/$chr";
	system "mkdir $chr_AS";
	my @gffs = glob "$chr_dir/*.gff";

	my $i = 0;
	my $sh;
	foreach my $gff (@gffs) {
		$gff =~ /.*\/(.*)\.gff/;
		my $gene = $1;
		print OUT "$chr\t$gene\n";
		if ($i < 10) {
			$sh .= "$gene_model_to_splicegraph -a -g $gene -m $gff -o $chr_AS/$gene.splicegraph.gff && ";
			$i++;
		}
		else {
			print SH "$sh";
			print SH "$gene_model_to_splicegraph -a -g $gene -m $gff -o $chr_AS/$gene.splicegraph.gff \n";
			$sh =();
			$i = 0;
		}
	}
}
close SH;
close OUT;

&Cut_shell_qsub ("$ARGV[1]/work_sh/$ARGV[2].SpliceGraphgff.sh",30,"1G",$config{queue});
&Check_qsub_error ("$ARGV[1]/work_sh/$ARGV[2].SpliceGraphgff.sh");

############ subs
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	if ($line<=1000) {
		system "$config{qsub}--queue $queue --maxproc $cpu --resource vf=$vf --reqsub  $shell";
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
