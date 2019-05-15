#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# # GetOptions, Check format, Output info.
# # ------------------------------------------------------------------
my ($in, $out,);
GetOptions(
                                "help|?" =>\&USAGE,
                                "in:s"  =>\$in,
                                "out:s"   =>\$out,
                      #          "od:s"    =>\$od,
                                ) or &USAGE;
&USAGE unless ($in and $out) ;
#$od  = &ABSOLUTE_DIR($od);
#$in = &ABSOLUTE_DIR($in);
#$out = &ABSOLUTE_DIR($out);
open (IN, $in) or die $!;
open (OUT, ">$out")or die $!;
print OUT "#LncRNA_ID\tmiRNA_hp\tLncRNA_start\tLncRNA_end\n";
while (<IN>) {
	chomp;
	next if (/^$/);
	my @tmp=split/\s+/,$_;
	print OUT "$tmp[1]\t$tmp[0]\t$tmp[8]\t$tmp[9]\n";
}
close IN;
close OUT;






sub MKDIR
{ # &MKDIR($out_dir);
        my ($dir)=@_;
#       rmdir($dir) if(-d $dir);
               mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR{
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





sub USAGE{
	print <<"	USAGE";        

Contact: niulg <niulg\@biomarker.com.cn>

Description:
        	This program is a Procedure deal with the result of blastn 

Usage:
	        -in                <file>	 blastn result    must be given contain at least 10 colunm
        	-out               <file>	output file
perl scripts.pl -in file.txt -out out.txt;
USAGE		
	
	USAGE
	exit; 
}







