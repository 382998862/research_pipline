#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";
my $version="v1.0.2";
my %config=%{readconf("$Bin/../../db_file.cfg")};
my $BEGIN_TIME=time();
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($fa, $key, $od,);
GetOptions(
				"help|?" =>\&USAGE,
				"fa:s"  =>\$fa,
				"od:s"   =>\$od,
				"key:s"    =>\$key,
				) or &USAGE;
&USAGE unless ($fa and $key and $od) ;
################################


&MKDIR($od);
$od  = &ABSOLUTE_DIR($od);     &MKDIR("$od/work_sh");

#&MKDIR("$od/precursor");


#==================================================================
# bins 
#==================================================================

my $MAKEBLASTDB   = $config{makeblastdb};   # 2014-12-17 ~ 
my $BLASTN    = $config{blastn};    # 2014-12-17 ~ 
my $HAIRPINFA = $config{HAIRPIN};
#==================================================================
# pipeline 
#==================================================================

################select hp_mirna.fa from mibase 
`perl $Bin/select_fa.pl -i $key -fa $HAIRPINFA -o $od/osa_hp.fa`;
############# make database format
`$MAKEBLASTDB -in $fa -input_type fasta -dbtype nucl -out $od/lnc_filter_final.fa`;
###########blastn
`$BLASTN  -evalue 1e-5 -db $od/lnc_filter_final.fa -query $od/osa_hp.fa -num_threads 2 -out $od/precursor.txt -outfmt 6`;
############### result deal
#
`perl $Bin/result.pl -in $od/precursor.txt -out $od/lncRNA_precursor.txt`;






sub GetTMR {#
	#Get Total Mapped Reads from file line which contain string "Mapped Reads\s"
	my $fStat=shift;
	open (IN,"<",$fStat) or die $!;
	while (<IN>) {
		if (/^Mapped Reads\s(\d+)/) {
			close (IN) ;
			return $1;
		}
	}
	close (IN) ;
	die "Error Reads Stat file.\n";
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
#	rmdir($dir) if(-d $dir);
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

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {
	my $usage=<<"USAGE";
Description:
	This program is a Procedure to judge  if the lncRNA is the procursor of hairpin fa  


Usage:
        -fa             lncRNA fasta sequence
	-od               output dir            must be given
	-key            the specise  abbreviation forexample "osa" in the miRBase database,   must be given;
USAGE
	print $usage;
	exit;
}
