use Getopt::Std;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my ($gene_count,$gene_fpkm,$trans_count,$trans_fpkm,$gene_gtf,$lnc_gtf,$novel_lnc,$od,$fpkm);
GetOptions(
	"od:s"		=>\$od,
	"gene_count:s"	=>\$gene_count,
	"gene_fpkm:s"	=>\$gene_fpkm,
	"trans_count:s"	=>\$trans_count,
	"trans_fpkm:s"	=>\$trans_fpkm,
	"gene_gtf:s"	=>\$gene_gtf,
	"lnc_gtf:s"	=>\$lnc_gtf,
	"novel_lnc:s"	=>\$novel_lnc,
	"fpkm:s"	=>\$fpkm,
	"help|h"	=>\&USAGE,
    ) or &USAGE;
&USAGE unless ($od);

my %genes=();
if(defined $gene_gtf){
	open(GTF,$gene_gtf)||die $!;
	while(<GTF>){
		chomp;
		my @tmp=split(/\t/,$_);
		$tmp[8]=~/gene_id (.*?)\;/;
		my $id=$1;$id=~s/"//g;
		$genes{$id}++;
	}
	close(GTF);
	&writeOut($gene_count,"$od/All_gene_counts.list")	if(defined $gene_count);
	&writeOut($gene_fpkm,"$od/All_gene_fpkm.list")       if(defined $gene_fpkm);
}

%genes=();
if(defined $lnc_gtf){
        open(GTF,$lnc_gtf)||die $!;
        while(<GTF>){
                chomp;
                my @tmp=split(/\t/,$_);
                $tmp[8]=~/transcript_id (.*?)\;/;
                my $id=$1;$id=~s/"//g;
                $genes{$id}++;
        }
        close(GTF);
        &writeOut($trans_count,"$od/Known_lncRNA_counts.list")       if(defined $trans_count);
        &writeOut($trans_fpkm,"$od/Known_lncRNA_fpkm.list")       if(defined $trans_fpkm);
}

%genes=();
if(defined $novel_lnc){
        open(GTF,$novel_lnc)||die $!;
        while(<GTF>){
                chomp;
                my @tmp=split(/\t/,$_);
                $tmp[8]=~/transcript_id (.*?)\;/;
                my $id=$1;$id=~s/"//g;
                $genes{$id}++;
        }
        close(GTF);
	if(defined $fpkm){
		open(FPKM,$trans_fpkm)||die $!;
		while(<FPKM>){
			chomp;my @tmp=split(/\t/,$_);
			my $id=shift @tmp;
			my $length=pop @tmp;
			next if(!exists $genes{$id});
			my $num=0;
			for(my $i=0;$i<@tmp;$i++){
				$num++	if($tmp[$i]<$fpkm);
			}
			delete $genes{$id}	if(@tmp==$num);
		}
		close(FPKM);
	}
        &writeOut($trans_count,"$od/Candicate_lncRNA_counts.list")       if(defined $trans_count);
        &writeOut($trans_fpkm,"$od/Candicate_lncRNA_fpkm.list")       if(defined $trans_fpkm);
}





sub writeOut{
	my ($i,$o)=@_;
	open(IN,$i)||die $!;
	open(OUT,">$o")||die $!;
	while(<IN>){
		chomp;
		if($_=~/^#/){
			print OUT "$_\n";
		}else{
			my @tmp=split;
			print OUT "$_\n"	if(exists $genes{$tmp[0]});
		}
	}
	close(IN);
	close(OUT);
}

sub USAGE {
        my $usage=<<"USAGE";
----------------------------------------------------------------------------------------------------------------------------
   Program: 
   Version: 
   Usage:
            --gene_count	<FILE>  gene_count_matrix.list
            --gene_fpkm		<FILE>  gene_fpkm_matrix.list
            --od        	<DIR>   analysis output directory
            --trans_count       <FILE>  transcript_count_matrix.list;
            --trans_fpkm        <FILE>  transcript_fpkm_matrix.list;
	    --gene_gtf		<FILE>	GENE GTF
	    --lnc_gtf		<FILE>	Known lncRNA GTF
	    --novel_lnc		<FILE>	Novel lncRNA GTF
	    --fpkm		<float>	default null. If given this parameter, if all samples of novel lncRNA fpkm less than this value, the novel lncRNA would be ignored.
		
            --h                 help documents

   Example:
  

----------------------------------------------------------------------------------------------------------------------------
USAGE
        print $usage;
        exit;
}

