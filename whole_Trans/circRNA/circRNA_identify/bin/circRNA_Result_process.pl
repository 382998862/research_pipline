#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $path = $Bin;
$path = substr($path,0,index($path,"bin"));
my %config=%{readconf("$path/project.cfg")};
my $Title=$config{Title};												#流程的名称，必填
my $version=$config{version};
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($indir,$odir);

GetOptions(
	"indir:s" =>\$indir,
	"odir:s"   =>\$odir,
	"help|h" =>\&USAGE,
	) or &USAGE;
&USAGE unless ($indir and $odir);
my @file = glob "$indir/*.circ";
my $flag = 0;
"rm $odir/All_CircRNA.tmp" if(-e "$odir/All_CircRNA.tmp");
foreach my $file(@file)
{
	my $samp = basename $file;
	$samp=~s/.circ//;
	`head -n 1 $file > $odir/All_CircRNA.tmp`,$flag =1 unless ($flag);
	`grep -v \"#\"  $file|awk  '{{for(i=1;i<NF;i++)printf(\"\%s\\t\",\$i)} len=split(\$NF,a,\",\");for(j=1;j<len;j++) printf(\"\%s,\",\"$samp:\"a[j]);printf(\"\\n\")}' >> $odir/All_CircRNA.tmp`;
	#print "grep -v \"#\"  $file|awk  '{{for(i=1;i<NF;i++)printf(\"\%s\\t\",\$i)} len=split(\$NF,a,\",\");for(j=1;j<len;j++) printf(\"\%s,\",\"$samp:\"a[j]);printf(\"\\n\")}' >> $odir/All_CircRNA.tmp";
}
open IN,"$odir/All_CircRNA.tmp" or die $!;
open OUT,">$odir/All_CircRNA.xls" or die $!;
my %hash;
my $head=<IN>;
print OUT "$head";
while(<IN>)
{
	chomp;
	my @line = split/\t/,$_;
	my $key = join "\t",@line[0..3],$line[8],$line[10];

	if(exists $hash{$key})
	{
		my @tmp = split/\t/,$hash{$key};
		my ($junction_reads,$SM_MS_SMS,$non_junction_reads,$junction_reads_ratio,$gene_id,$junction_reads_ID)=@tmp;
		$junction_reads+=$line[4];
		my @tmp_sms1 = split /_/,$SM_MS_SMS;
		my @tmp_sms2 = split /_/,$line[5];
		$SM_MS_SMS =  join "_",($tmp_sms1[0]+$tmp_sms1[0]),($tmp_sms1[1]+$tmp_sms1[1]),($tmp_sms1[2]+$tmp_sms1[2]);
		$non_junction_reads+=$line[6];
		$junction_reads_ratio=2*$junction_reads/(2*$junction_reads+$non_junction_reads);
		my @tmp_gene_id1=split /,/,$gene_id;
		my @tmp_gene_id2=split /,/,$line[9];
		push @tmp_gene_id1, @tmp_gene_id2;
		my %tmp_hash;
		@tmp_gene_id1 = grep { ++$tmp_hash{$_} < 2 } @tmp_gene_id1;
		$gene_id =(join ",",@tmp_gene_id1).",";
		$junction_reads_ID.=$line[11];
		$hash{$key}=join"\t",$junction_reads,$SM_MS_SMS,$non_junction_reads,$junction_reads_ratio,$gene_id,$junction_reads_ID;
	}else{
		my $value =  join "\t",@line[4..7],$line[9],$line[11];
		$hash{$key}=$value;
	}
}
close IN;
foreach my $key (keys%hash)
{
	my @tmp1=split /\t/,$key;
	my @tmp2=split /\t/,$hash{$key};
	print OUT "$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp2[0]\t$tmp2[1]\t$tmp2[2]\t$tmp2[3]\t$tmp1[4]\t$tmp2[4]\t$tmp1[5]\t$tmp2[5]\n";
}
close OUT;



#############################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: $Script
   Version: $version
   Contact: songmm <songmm\@biomarker.com.cn>
      Date: 2016-08-10

     Usage:
            --indir      <FILE>  all sample circRNA result dir
            --odir      <FILE>  output dir
            --h                 help documents

   Example:
            perl $Script --indir CircRNA_identify/ --odir  CircRNA_identify/

----------------------------------------------------------------------------------------------------------------------------
USAGE
	print $usage;
	exit;
}