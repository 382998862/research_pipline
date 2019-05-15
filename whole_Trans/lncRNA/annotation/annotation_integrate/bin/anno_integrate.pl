#!/usr/bin/perl -w
# 
#Copyright (c) BMK 2011 
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2011 
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2011 
my $version="1.1.0";
my $BEGIN=time();

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use Encode;
use Spreadsheet::WriteExcel;
use Spreadsheet::ParseExcel;  
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $programe_dir=basename($0);

##############################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作；
my %opts;
GetOptions(\%opts,"gene=s","id=s","od=s","key=s","h");
if (!defined($opts{gene}) || !defined($opts{id}) || !defined($opts{od}) || !defined($opts{key}) || defined($opts{h})) {
	&help;
	exit;
}

#######################
print "\n[".sub_format_datetime(localtime(time()))."] $Script start....\n\n";
######################
my $gene=&ABSOLUTE_DIR($opts{gene});
my $anno_dir=&ABSOLUTE_DIR($opts{id});
&MKDIR($opts{od});
my $odir=&ABSOLUTE_DIR($opts{od});
my $cds_key=$opts{key};

my %Gene_length;
open IN,"$gene" || die;
$/=">";
while (<IN>) {
	chomp;
	s/^>// if $.==1;
	next if (/^$/ || /^\#/);
	my ($title,$seq)=split/\n+/,$_,2;
	$seq=~s/\s+//g;
	my $id=(split/\s+/,$title)[0];
	my $length=length $seq;
	$Gene_length{$id}=$length;
}
$/="\n";
close IN;

my %Anno_Gene;
my %Anno_Stat;

my @Anno_files=glob "$anno_dir/*.anno.txt";
my $ko_anno=glob "$anno_dir/*.ko";
push @Anno_files,$ko_anno;
my %Anno_Base;

foreach my $file (@Anno_files) {
	my $database_key;
	if ($file=~/$cds_key\.Kegg.ko/) {
		$database_key="KEGG";
	}
	elsif ($file=~/$cds_key/ && $file=~/anno/) {
		$file=~m/$cds_key\.(.*).anno.txt/;
		$database_key=$1;
	}
	$Anno_Base{$database_key}=1;
	if ($database_key!~/GO/) {
		open IN,"$file" || die;
		while (<IN>) {
			chomp;
			next if (/^$/ || /^\#/);
			my @annotate=split/\t+/,$_;
			my $id=$annotate[0];
			my $anno=$annotate[-1];
			$Anno_Stat{All}{$id}=1;
			$Anno_Stat{$database_key}{Anno}++;
			$Anno_Gene{$id}{$database_key}=$anno;
			if ($Gene_length{$id}>=300 && $Gene_length{$id}<1000) {
				$Anno_Stat{$database_key}{300}++;
			}
			if ($Gene_length{$id}>=1000) {
				$Anno_Stat{$database_key}{1000}++;
			}
		}
		close IN;
	}
	else {
		open IN,"$file" || die;
		while (<IN>) {
			chomp;
			next if (/^$/ || /^\#/);
			my ($id,$anno)=split/\t+/,$_,2;
			$Anno_Stat{All}{$id}=1;
			$Anno_Stat{$database_key}{Anno}++;
			$anno=~s/^\d+\t+//;
			$anno=~s/\t+/\; /g;
			$Anno_Gene{$id}{$database_key}=$anno;
			if ($Gene_length{$id}>=300 && $Gene_length{$id}<1000) {
				$Anno_Stat{$database_key}{300}++;
			}
			if ($Gene_length{$id}>=1000) {
				$Anno_Stat{$database_key}{1000}++;
			}
		}
		close IN;
	}
}

my $cog=glob "$anno_dir/*.Cog_class.txt";
my $kog=glob "$anno_dir/*.Kog_class.txt";

# cog
if ($cog) {
    $Anno_Base{COG}=1;

    open COG,"$cog" || die;
    while (<COG>) {
        chomp;
        next if (/^$/ || /^\#/) ;
        my @anno=split/\t+/,$_;
        $Anno_Stat{All}{$anno[0]}=1;
        $Anno_Stat{COG}{Anno}++;
        $Anno_Gene{$anno[0]}{COG}="$anno[-3]\t$anno[-1]";
        if ($Gene_length{$anno[0]}>=300 && $Gene_length{$anno[0]}<1000) {
            $Anno_Stat{COG}{300}++;
        }
        if ($Gene_length{$anno[0]}>=1000) {
            $Anno_Stat{COG}{1000}++;
        }
    }
    close COG;
}

# kog 
if ($kog) {
    $Anno_Base{KOG}=1;

    open KOG,"$kog" || die;
    while (<KOG>) {
        chomp;
        next if (/^$/ || /^\#/) ;
        my @anno=split/\t+/,$_;
        $Anno_Stat{All}{$anno[0]}=1;
        $Anno_Stat{KOG}{Anno}++;
        $Anno_Gene{$anno[0]}{KOG}="$anno[-3]\t$anno[-1]";
        if ($Gene_length{$anno[0]}>=300 && $Gene_length{$anno[0]}<1000) {
            $Anno_Stat{KOG}{300}++;
        }
        if ($Gene_length{$anno[0]}>=1000) {
            $Anno_Stat{KOG}{1000}++;
        }
    }
    close KOG;
}

open OUT,">$odir/Integrated_Function.annotation.xls" || die;
print OUT "#GeneID";
foreach (sort keys %Anno_Base) {
	if (/COG/) {
		print OUT "\tCOG_class\tCOG_class_annotation";
	}
    elsif (/KOG/) {
        print OUT "\tKOG_class\tKOG_class_annotation";
    }
	else {
		print OUT "\t$_"."_annotation";
	}
}

print OUT "\n";

foreach (keys %{$Anno_Stat{All}}) {
	$Anno_Stat{zz}{Anno}++;
	my $id=$_;
	print OUT "$id";
	if ($Gene_length{$id}>=1000) {
		$Anno_Stat{zz}{1000}++;
	}
	if ($Gene_length{$id}>=300 && $Gene_length{$id}<1000) {
		$Anno_Stat{zz}{300}++;
	}
	foreach (sort keys %Anno_Base) {
		if (/COG/ && !defined $Anno_Gene{$id}{COG}) {
			print OUT "\t--\t--";
		}
        elsif (/KOG/ && !defined $Anno_Gene{$id}{KOG}) {
            print OUT "\t--\t--";
        }
		elsif (defined $Anno_Gene{$id}{$_}) {
			printf OUT "\t$Anno_Gene{$id}{$_}";
		}
		else {
			printf OUT "\t--";
		}
	}
	print OUT "\n";
}
close OUT;

open STAT,">$odir/Function_Annotation.stat.xls" || die;
print STAT "#Anno_Database\tAnnotated_Number\t300<=length<1000\tlength>=1000\n";
foreach (sort keys %Anno_Stat) {
	next if (/All/) ;
	if ($_!~/zz/) {
		print STAT "$_"."_"."Annotation\t$Anno_Stat{$_}{Anno}\t$Anno_Stat{$_}{300}\t$Anno_Stat{$_}{1000}\n";
	}
	else {
		print STAT "All_Annotated\t$Anno_Stat{$_}{Anno}\t$Anno_Stat{$_}{300}\t$Anno_Stat{$_}{1000}\n";
	}
}

close STAT;

################################################

my $stat="$odir/All_Database_annotation.xls";
$stat = Spreadsheet::WriteExcel->new("$stat");
my @filter;
push @filter,"$odir/Integrated_Function.annotation.xls";
push @filter,"$odir/Function_Annotation.stat.xls";
if (defined $Anno_Base{GO}) {
	push @filter,"$anno_dir/$cds_key.GO.list.txt";
	push @filter,"$anno_dir/$cds_key.GO_tree.stat.xls";
}
if (defined $Anno_Base{KEGG}) {
	push @filter,"$anno_dir/$cds_key.Kegg.pathway";
	push @filter,"$anno_dir/$cds_key.Kegg.ko";
}
foreach my $f (@filter) {
	my $name=basename($f);
	$name=~s/\.txt$//;
	$name=~s/\.xls$//;
    $name=~s/annotation/anno/i;
	#$name=~s/.*\.fa\.//;
	$name=~s/$cds_key\.//;
#	Beta($f,$stat,$name,0);
	Beta_x($f,$stat,$name,0);
}



###############Time
print "[".sub_format_datetime(localtime(time()))."] $Script done.\n\n";
###########subs
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

sub Beta { #将linux下的xls文件转换成windows下可读的excel文件
	my ($file,$xls,$sheet,$f)=@_;
	open (IN,"$file")||die "can't open file $file\n";
	my $worksheet = $xls->add_worksheet(decode('GB2312',"$sheet"));
	my $format = $xls->add_format();
	$format->set_border($f);
	$format->set_align('left');
#	$format->set_font('Times New Roman');
	$format->set_size('11') if($f ==0);
	$format->set_font('宋体') if($f ==0);
	my $num=0;

	while(<IN>) {
		chomp;
		my @line=split/\t/,$_;
		for(my $i=0;$i<@line;$i++)
		{
			$format->set_num_format('#,###') if($f==1 && $line[$i]=~/\d+/ && $line[$i]!~/[\.\%]/);
			$worksheet->write($num,$i,decode('GB2312',"$line[$i]"),$format);
		}
		$num++;
	}

	close IN;
}

sub Beta_x { #将linux下的xls文件转换成windows下可读的excel文件
	my ($file,$xls,$sheet,$f)=@_;

	open (IN,"$file")||die "can't open file $file\n";
    chomp(my $line_num= `wc -l $file|cut -f 1 -d ' '`);
    my $sheetNum = 1;
	my $num=0;
    my $worksheet;

    if ($line_num>65000) {
        $worksheet = $xls->add_worksheet(decode('GB2312',"$sheet.$sheetNum"));
    } else {
        $worksheet = $xls->add_worksheet(decode('GB2312',"$sheet"));
    }
	
	my $format = $xls->add_format();
	$format->set_border($f);
	$format->set_align('left');
#	$format->set_font('Times New Roman');
	$format->set_size('11') if($f ==0);
	$format->set_font('宋体') if($f ==0);


	while(<IN>) {
		chomp;
		my @line=split/\t/,$_;
		for(my $i=0;$i<@line;$i++) {
			$format->set_num_format('#,###') if($f==1 && $line[$i]=~/\d+/ && $line[$i]!~/[\.\%]/);
			$worksheet->write($num,$i,decode('GB2312',"$line[$i]"),$format);
		}
		$num++;

        if ($.%65000==0) {
            $num=0;
            $sheetNum++;
            $worksheet = $xls->add_worksheet(decode('GB2312',"$sheet.$sheetNum"));
        }
	}

	close IN;
}

sub ABSOLUTE_DIR { #$pavfile=&ABSOLUTE_DIR($pavfile);
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

sub help {
    my $usage =<<"_USAGE_";
#-------------------------------------------------------------------------------------------------
    Program: $Script
    Version: $version
     Writer: mengf <mengf\@biomarker.com.cn>
       Data: 2012-00-00
   Modifier: Simon Young <simonyoung8824\@gmail.com>
       Data: 2014-10-08
Description: Extract and convert annotation result to excel format for DEG Analysis.
             v1.0.0-based, and deal with excel overflow bug by spliting file with more than 65,000 lines.

      Usage:
            -gene   Unigene fa file ,to state The gene length;
            -id     Annotation OUT DIR (eg. ./Result);
            -od     Anno Integerate and Stat dir (eg. /Result);
            -key    mRNA fa name (mRNA seq file name)

#-------------------------------------------------------------------------------------------------
_USAGE_
    print $usage;
    exit;
}
