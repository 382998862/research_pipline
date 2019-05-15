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
## out_dir
#$outdir.="/" unless $outdir =~/\/$/;
mkdir $outdir unless (-e $outdir);
$outdir=&ABSOLUTE_DIR($outdir);

## match_genome file
my %hash;
open I, "$infile";
while(<I>){
  chomp;
  if($mod==0){
  	my @d=split;my @base=split (//,$d[6]);
  	next if (length $d[6] > 30 or length $d[6] < 20);
  	foreach my $i(0..$#base){
  		$hash{$i+1}{$base[$i]}+=$d[3];
  	}
	}
	elsif($mod==1){
    if (s/^>//){
			$_=~/(.*)_(\d+)_x([\d]+)$/;
			my $cnt=$3;
	  my $seq=<I>;chomp $seq;
	  $seq= uc $seq;
	  $seq=~tr/Tt/Uu/;
      my @base=split (//,$seq);
      next if (length $seq > 30 or length $seq < 20);
      foreach my $i(0..$#base){$hash{$i+1}{$base[$i]}+=1;}
    }
	}
	else{
		die("please set input file type for match_hairpin.txt or mir-novel.fa using 0 or 1!\n");
	}
}close I;
my $out_file = "$outdir/$p.base_bias";
# my $out_file="$outdir/base_bias"if $mod==0;
# $out_file="$outdir/base_bias"if $mod==1;
open O, ">$out_file.stat";
print O join("\t",'pos','A','U','C','G'),"\n";
foreach my $pos (sort {$a<=>$b} keys %hash){
	print O $pos, "\t";
	print O $hash{$pos}{'A'}=$hash{$pos}{'A'}?$hash{$pos}{'A'}:0 ,"\t";
	print O $hash{$pos}{'U'}=$hash{$pos}{'U'}?$hash{$pos}{'U'}:0 ,"\t";
	print O $hash{$pos}{'C'}=$hash{$pos}{'C'}?$hash{$pos}{'C'}:0, "\t";
	print O $hash{$pos}{'G'}=$hash{$pos}{'G'}?$hash{$pos}{'G'}:0, "\t";
	print O "\n";
}close O;

my $R=<<SCIRIPT;
read.table("$out_file.stat",skip=1)->a
for(i in 1:nrow(a)) a[i,2:5]/sum(a[i,2:5])->a[i,2:5]
col=c("hotpink","green","steelblue1","lightgoldenrod")
png("$out_file.png",width=1820,height=900)
par(mar=c(5,4.5,8,8))
barplot(
	t(a[,2:5]),names.arg=a[,1],col=col,las=1,cex.main=2,
	xlab="Position(nt)",ylab="Percent(%)",border=F,
	main="miRNA nucleotide bias at each position",yaxt='n'
)
axis(2,at=axTicks(2),las=2,labels=(axTicks(2)*100))
legend(
	par("usr")[2],1,pch=15,col=col,
	c("A","U","C","G"),xpd=T,bty='n',cex=2
)
abline(h=axTicks(2),col="gray")
dev.off()
pdf("$out_file.pdf",width=8,height=4)
par(mar=c(5,4.5,8,8))
barplot(
	t(a[,2:5]),names.arg=a[,1],col=col,las=1,cex.main=2,
	xlab="Position(nt)",ylab="Percent(%)",border=F,
	main="miRNA nucleotide bias at each position",yaxt='n'
)
axis(2,at=axTicks(2),las=2,labels=(axTicks(2)*100))
legend(
	par("usr")[2],1,pch=15,col=col,
	c("A","U","C","G"),xpd=T,bty='n',cex=2
)
abline(h=axTicks(2),col="gray")
dev.off()
SCIRIPT

open R,">$outdir/.base_pos.R" or die;
print R $R;
close R;
runOrDie("$CFG{Rscript} $outdir/.base_pos.R && rm $outdir/.base_pos.R");


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
