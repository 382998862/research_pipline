use FindBin qw($Bin $Script);
use Cwd qw(abs_path getcwd);
use File::Basename;
my($file,$deg,$od)=@ARGV;
if(@ARGV!=2&&@ARGV!=3){
	print "perl $0\n1. *_KEGG_pathway_enrich.list\n2. DEG_final.xls \n3. outdir(default the same dir with KEGG_pathway_enrich.list)\n";
	exit;
}
$od   ||=dirname $file;
`mkdir -p $od`	unless(-d $od);
$od	=abs_path($od);
$deg	=abs_path($deg);
$file	=abs_path($file);
my %degs=();
open(DEG,$deg)||die $!;
while(<DEG>){
	chomp;next if($_=~/^#/);
	my @tmp=split;
	$degs{$tmp[0]}=$tmp[-1];
}
close(DEG);
my $base=basename $deg;
$base=(split(/\./,$base))[0];
if(!-e $file){exit;}
open(KEGG,$file)||die $!;
open(OUT,">$od/$base.KEGG.list")||die $!;
#print OUT "pathway_term\trich_factor\tqvalue\tgene_number\tshape\n";

print OUT "ID\tDescription\tGeneRatio\tBgRatio\tenrich_factor\tpvalue\tqvalue\tgeneID\tgene_number\tContained\n";
my $num=0;
while(<KEGG>){
	next if($num==20);
	chomp;next if($_=~/^#|^ID/);
	my @tmp=split(/\t/,$_);
	my @gs=split(/\/|\;/,$tmp[7]);
	foreach my $g(@gs){
		push @{$tmp[1]},$degs{$g};
	}
	my $type=&change(\@{$tmp[1]});
	print OUT join("\t",@tmp),"\t",scalar(@gs),"\t$type\n";
	$num++;
}
close(KEGG);
close(OUT);

my $Rscript="/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";
print "$Rscript $Bin/KEGG_dot.r --infile $od/$base.KEGG.list --outdir $od --pathway_term.col 2 --qvalue.col 7 --rich_factor.col 5 --gene_number.col 9 --outfile $base.Phase --height 900 --width 1200\n";
`$Rscript $Bin/KEGG_dot.r --infile $od/$base.KEGG.list --outdir $od --pathway_term.col 2 --qvalue.col 7 --rich_factor.col 5 --gene_number.col 9 --outfile $base.Phase --height 900 --width 1200`;


sub change{
	my $s=shift;
	my @set=@{$s};
	my %hash=();
	foreach my $s(@set){
		$hash{$s}++;
	}
	if(exists $hash{up} && !exists $hash{down}){
		return "up";
	}elsif(exists $hash{down} && !exists $hash{up}){
		return "down";
	}elsif(exists $hash{down} && exists $hash{up}){
		return "up&down";
	}
}

