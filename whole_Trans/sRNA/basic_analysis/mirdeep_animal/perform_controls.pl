#!/usr/bin/perl

use warnings;
use Cwd;
use strict;
use Getopt::Std;
use File::Copy;
use File::Path;
use FindBin qw($Bin $Script);
use newPerlBase;

my %CFG=%{readconf("$Bin/../../CFG")};

my $usage=
"$0 file_command_line file_structure rounds_controls

-a   Output progress to screen
";

my $file_command_line=shift or die $usage;
my $file_structure=shift or die $usage;
my $rounds=shift or die $usage;

#options
my %options=();
getopts("a",\%options);


my $notename=`hostname`;chomp $notename;
my $ltime=time();
my $work_dir=`pwd`;chomp $work_dir;
my $dir="dir_perform_controls$ltime";

my $command_line=parse_file_command_line($file_command_line);

if($options{a}){print STDERR "total number of rounds controls=$rounds\n";}

perform_controls();


sub perform_controls{

    mkdir $dir;
    mkdir "$work_dir/work_sh" unless -d "$work_dir/work_sh";

    my $round=1;

	system("permute_structure.pl $file_structure > $dir/precursors_permuted.str 2> /dev/null");
	open (SH,">$work_dir/work_sh/miRDeep2_core_algorithm.sh") or die $!;
    while($round<=$rounds){

		my $ret="cd $work_dir && $command_line > work_sh/Cycle.$round.out.xls 2> /dev/null ";
		print SH "$ret \n";

		$round++;
    }
	close SH;
	qsubOrDie("$work_dir/work_sh/miRDeep2_core_algorithm.sh",$CFG{queue},$CFG{cpu},$CFG{vf});
	$round=1;
    while($round<=$rounds){
		print "permutation $round\n\n";
		open (IN,"$work_dir/work_sh/Cycle.$round.out.xls") or die $!;
		while (<IN>) {
			print $_;
		}
		close IN;
		$round++;
    }

    rmtree($dir);

    if($options{a}){print STDERR "controls performed\n\n";}
}


sub parse_file_command_line{

    my ($file) = @_;

    open (FILE, "<$file") or die "can not open $file\n";
    while (my $line=<FILE>){

	if($line=~/(\S+)/){

	    chomp $line;

	    $line=~s/$file_structure/$dir\/precursors_permuted.str/;

	    $line=~s/>.+//;

	    return $line;

	}
    }
    die "$file is empty\n";
}


