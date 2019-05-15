#!/usr/bin/perl -w
use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);
use threads;
use newPerlBase;
use List::Util qw/sum/;;

my ($fpkm, $cfg, $odir);
GetOptions(
        "h|?"           =>\&USAGE,
        "cfg:s"        =>\$cfg,
	"fpkm:s"	=>\$fpkm,
        "od:s"           =>\$odir,
)or &USAGE;
&USAGE unless ($cfg and $fpkm and $odir);


#################调用脚本
my $Rscript="/share/nas2/genome/biosoft/R/3.5.2/bin/Rscript";
my $TF_activity="$Bin/CoRegNet_ref_trans.R";

########extract info;
$cfg = abs_path($cfg);
$fpkm = abs_path($fpkm);
my $name = basename($fpkm);
###########
#
my $tf_file;
open(CFG,$cfg)||die $!;
while (<CFG>) {
	chomp;next if(/^#/);
	if(/^TFDB/) {
		$tf_file=(split/\s+/,$_)[1];$tf_file=abs_path($tf_file);
	}
}

##########
open (TF,"$tf_file")||die $!;
open (TFLIST,">$odir/TFlist.txt")||die $!;
my %reship=();
while (<TF>) {
	chomp;next if(/^Species/);
	my @tmp=split(/\s+/,$_);
	my $tf=$tmp[1]."\-\-".$tmp[3];
	$reship{$tmp[2]}=$tf;
	print TFLIST $tf,"\n";
	
}

########
open (FPKM,"$fpkm")||die $!;
my $header=<FPKM>;
chomp($header);
my @tmp=split(/\s+/,$header);
my @headers;
for(my $i=1;$i<@tmp;$i++){
	next if($tmp[$i]=~m/(\.FDR$)/ || $tmp[$i]=~m/(\.regulated$)/ || $tmp[$i]=~m/(\.log2FC$)/ || $tmp[$i]=~m/(\.PValue$)/);
	push @headers,$tmp[$i];
}

open(OUT,">$odir/$name.fpkm")||die $!;
print OUT "ID\t",join("\t",@headers),"\n";

while (<FPKM>) {
	chomp;
	my @tmp=split(/\s+/,$_);
	my $gene=$tmp[0];
	my @samples;
	for(my $i=1;$i<@headers+1;$i++){
		push @samples,$tmp[$i];
	}
	my $sum=sum @samples;
	next if($sum < 1);
	if (exists $reship{$gene}) {
		print OUT $reship{$gene},"\t",join("\t",@samples),"\n";
	}else{	
		print OUT $gene,"\t",join("\t",@samples),"\n";
	}

}

close(OUT);
close(TFLIST);
close(TF);
close(CFG);

my $lines=`wc -l $odir/$name.fpkm`;
my $Allnum=(split(/\s+/,$lines))[0]-1;
my $gene=`grep -P "^E" $odir/$name.fpkm |wc -l`;
my $gene_num=(split(/\s+/,$gene))[0];
print $Allnum,"\t",$gene_num,"\n";

if($Allnum ne $gene_num) {
	my $cmd;
	$cmd="$Rscript $TF_activity --fpkm $odir/$name.fpkm --TF $odir/TFlist.txt --out $odir/TFs";
	print "$cmd\n";
	system($cmd);
	`convert $odir/TFs_influences_heatmap.pdf $odir/TFs_influences_heatmap.png`;
	`convert $odir/TFs_network.pdf $odir/TFs_network.png`
}else{
	print "This analysis cannot be done. There are at least 2 of the provided regulators/transcription factor (TFlist) in the rownames of the gene expression matrix.\n";
}


sub USAGE       #use function
{
        my $usage=<<"USAGE";
    Description:
        Contact: lijj\@biomarker.com.cn
        version:   v1.0.0 at 2018.10.10 
        Description:  The prediction of transcriptor factor bingding sites of genes you providing;
    Usage:   
        Options:
        -cfg            Configure file, <required>;
	-fpkm		fpkm file, <required>;
        -od             the output path
        -h              Help, display this help info.
        
USAGE
        print $usage;
        exit;  #die program
}




