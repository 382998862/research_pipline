use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($infile, $outdir, $p);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$infile,
				"o:s"=>\$outdir,
				"p:s"=>\$p,
				) or &USAGE;
&USAGE unless ($infile and $outdir and $p);

my %CFG=%{readconf("$Bin/../../CFG")};

my $mod=1;
mkdir $outdir unless (-e $outdir);
$outdir=&ABSOLUTE_DIR($outdir);
my %hash;
$/=">";
open I, "$infile";
<I>;
while(<I>){
	chomp;
	my ($head,$seq)=split/\n/,$_;
	$seq  = uc($seq);
	$seq=~tr/Tt/Uu/;
	my $len=length($seq);
	my $fist_base=substr($seq, 0, 1);
	$hash{$len}{$fist_base}+=1;

}
close I;
my $out_file = "$outdir/$p.first_base";
# my $out_file="$outdir/first_base" if $mod==1;
# $out_file="$outdir/first_base" if $mod==0;
open O, ">$out_file.stat";
print O join("\t",'len','A','U','C','G'),"\n";
foreach my $len (sort {$a<=>$b} keys %hash){
	last if $len>30;
	print O $len, "\t";
	print O (exists $hash{$len}{'A'})?$hash{$len}{'A'}:0,"\t";
	print O (exists $hash{$len}{'U'})?$hash{$len}{'U'}:0,"\t";
	print O (exists $hash{$len}{'C'})?$hash{$len}{'C'}:0,"\t";
	print O (exists $hash{$len}{'G'})?$hash{$len}{'G'}:0,"\t";
	print O "\n";
}close O;
my $num=keys %hash;

my $R=<<SCIRIPT;
read.table("$out_file.stat",skip=1)->a
col=c("hotpink","green","steelblue1","lightgoldenrod")
for(i in 1:nrow(a)) a[i,2:5]/sum(a[i,2:5])->a[i,2:5]
png("$out_file.png",width=1400,height=800)
par(mar=c(5,4.5,8,8))
barplot(
        t(a[,2:5]),names.arg=a[,1],col=col,las=1,cex.main=2,
        xlab="Lengthl(nt)",ylab="Percent(%)",yaxt='n',
        main="miRNA first nucleotide bias",border=F
)
axis(2,at=axTicks(2),las=2,labels=(axTicks(2)*100))
legend(
        par("usr")[2],1,pch=15,col=col,
        c("A","U","C","G"),xpd=T,bty='n',cex=2
)
abline(h=axTicks(2),col="gray")
dev.off()
for(i in 1:nrow(a)) a[i,2:5]/sum(a[i,2:5])->a[i,2:5]
pdf("$out_file.pdf",width=7,height=4)
par(mar=c(5,4.5,8,8))
barplot(
        t(a[,2:5]),names.arg=a[,1],col=col,las=1,cex.main=2,
        xlab="Lengthl(nt)",ylab="Percent(%)",yaxt='n',
        main="miRNA first nucleotide bias",border=F
)
axis(2,at=axTicks(2),las=2,labels=(axTicks(2)*100))
legend(
        par("usr")[2],1,pch=15,col=col,
        c("A","U","C","G"),xpd=T,bty='n',cex=2
)
abline(h=axTicks(2),col="gray")
dev.off()
SCIRIPT

open R,">$outdir/.first_bias.R" or die;
print R $R;
runOrDie("$CFG{Rscript} $outdir/.first_bias.R && rm $outdir/.first_bias.R");

############################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
	my $usage=<<"USAGE";
Usage:
  Options:

  -i <file>  input file,fasta form,forced

  -o <dir>   output dir,forced

  -p <str>   prefix of output,forced

  -h         Help

USAGE
	print $usage;
	exit;
}
