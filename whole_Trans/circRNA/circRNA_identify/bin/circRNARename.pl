#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($fIn,$fOut,$symbol,$label,$gff,$anno);

GetOptions(
                "help|?"	=>\&USAGE,
                "o:s"		=>\$fOut,
                "i:s"		=>\$fIn,
		"symbol:s"	=>\$symbol,
		"label:s"	=>\$label,
		"gff:s"		=>\$gff,
		"anno:s"	=>\$anno,
                ) or &USAGE;
&USAGE unless ($fIn and $fOut and $gff and $anno);

$label||="Novel";

my %sym=();
if(defined $symbol){
	open(SYM,$symbol)||die $!;
	while(<SYM>){
		chomp;next if($_=~/^#/);
		my @tmp=split(/\t/,$_);
		$sym{$tmp[0]}=$tmp[1];
	}
	close(SYM);
}

my(%h,%anno);
open(IN,"$gff") or die $!;
while(<IN>){
	chomp;
	next if(/^#|^$|^\s+$/);
	my @tmp =split /\t/,$_;
	if($tmp[2] eq "gene"){
		my $gene = $tmp[8];
		$gene =~ s/ID=//;
		$gene =~ s/;//;
		$h{$gene}=join("\t",$tmp[3],$tmp[4]);
	}
}
close(IN);

open(IN,$anno) or die $!;
my ($id1,$head1,$line);
while(<IN>){
	chomp;
	next if(/^$/);
	if(/^#/){
		($id1,$head1)=split /\t/,$_,2;
		my @count = split /\t/,$head1;
		my $num = @count;
		$line = "--\t" x $num;
		$line =~ s/\t$//g;
		next;
	}
	my ($g,$info)=split /\t/,$_,2;
	$anno{$g}=$info;
}

close(IN);

my $od=dirname $fOut;

open (IN,$fIn) or die $!;
open(BED,">$od/Circ.bed")||die $!;
open(TAR,">$od/Circ_source_gene.xls")||die $!;
print TAR "#ID\tGene\n";
open (OUT,">$fOut") or die $!;
my %hash=();
while (<IN>) {
	chomp;
	if(/circRNA_ID/){
		my @tmp = split /\t/,$_;
		shift @tmp;
		pop @tmp;
		my $head = join("\t",@tmp);
		if(defined $symbol){
			print OUT "#New_ID\tcircRNA_ID\tcircRNA_chr\tcircRNA_start\tcircRNA_end\tSource_gene_id\tSymbol\tSource_gene_start\tSource_gene_end\t$head\t$head1\n";
		}else{
			print OUT "#New_ID\tcircRNA_ID\tcircRNA_chr\tcircRNA_start\tcircRNA_end\tSource_gene_id\tSource_gene_start\tSource_gene_end\t$head\t$head1\n";
		}
		next;
	}
	my @lines=split(/\t/,$_);
	my $cid = shift @lines;
	my $chr = (split /:/,$cid)[0];
	my $start = (split /\|/,(split /:/,$cid)[1])[0];
	my $end = (split /\|/,(split /:/,$cid)[1])[1];
	my $gene = pop @lines;
	$gene=~s/,$//;
	my @gs=split(/,/,$gene);
	my @names=();
	foreach my $g(@gs){	
		my $flag=$g;
		$flag=$sym{$g}	if(exists $sym{$g});
		$flag="intergenic"	if($flag eq "n\/a" ||$flag eq "N\/A");
		if(exists $hash{$g}){
			$hash{$g}++;	
		}else{
			$hash{$g}=1;	
		}
		my $id="$label"."_".$flag."_".$hash{$g};
		push @names,$id;
		if(exists $h{$g}){
			if(defined $symbol && exists $sym{$g}){
				if(exists $anno{$g}){
					print OUT "$id\t$cid\t$chr\t$start\t$end\t$g\t$sym{$g}\t$h{$g}\t".join("\t",@lines)."\t$anno{$g}\n";
				}else{
					print OUT "$id\t$cid\t$chr\t$start\t$end\t$g\t$sym{$g}\t$h{$g}\t".join("\t",@lines)."\t$line\n";
				}
			}elsif(defined $symbol && !exists $sym{$g}){
				if(exists $anno{$g}){	
					print OUT "$id\t$cid\t$chr\t$start\t$end\t$g\t$g\t$h{$g}\t".join("\t",@lines)."\t$anno{$g}\n";
				}else{
					print OUT "$id\t$cid\t$chr\t$start\t$end\t$g\t$g\t$h{$g}\t".join("\t",@lines)."\t$line\n";
				}
			}elsif(!defined $symbol){
				if(exists $anno{$g}){
					print OUT "$id\t$cid\t$chr\t$start\t$end\t$g\t$h{$g}\t".join("\t",@lines)."\t$anno{$g}\n";
				}else{
					print OUT "$id\t$cid\t$chr\t$start\t$end\t$g\t$h{$g}\t".join("\t",@lines)."\t$line\n";
				}
			}else{print "War1!\n";}
		}else{
			if(defined $symbol && exists $sym{$g}){
				print OUT "$id\t$cid\t$chr\t$start\t$end\t$g\t$sym{$g}\t--\t--\t".join("\t",@lines)."\t$line\n";
			}elsif(defined $symbol && !exists $sym{$g}){
				print OUT "$id\t$cid\t$chr\t$start\t$end\t$g\t$g\t--\t--\t".join("\t",@lines)."\t$line\n";
			}elsif(!defined $symbol){
				print OUT "$id\t$cid\t$chr\t$start\t$end\t$g\t--\t--\t".join("\t",@lines)."\t$line\n";
			}else{print "War2!\n";}
		}
	}
	print TAR "$cid\t",join(",",@gs),"\n"	if($gene ne "n\/a" && $gene ne "N\/A");
	print BED join("\t",split(/:|\|/,$cid)),"\t$cid\n";
}

close IN;
close OUT;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

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

################################################################################################################

sub max{#&max(lists or arry);
    #求列表中的最大值
    my $max=shift;
    my $temp;
    while (@_) {
        $temp=shift;
        $max=$max>$temp?$max:$temp;
    }
    return $max;
}

################################################################################################################

sub min{#&min(lists or arry);
    #求列表中的最小值
    my $min=shift;
    my $temp;
    while (@_) {
        $temp=shift;
        $min=$min<$temp?$min:$temp;
    }
    return $min;
}

################################################################################################################

sub revcom(){#&revcom($ref_seq);
    #获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
    my $seq=shift;
    $seq=~tr/ATCGatcg/TAGCtagc/;
    $seq=reverse $seq;
    return uc $seq;
}

################################################################################################################

sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
    my $usage=<<"USAGE";
 ProgramName:
     Version:   $version
     Contact:   Simon Young <yangxh\@biomarker.com.cn>
Program Date:   2012.07.02
      Modify:
 Description:   This program is used to ......
       Usage:
        Options:
        -i <file>   input file,xxx format,forced
        -o <file>   output file,optional
	-symbol	<file>	contained 2 column: ENSG_ID symbol, not must
	-label	<str>	newname prefix,not must
	-gff <file>	input file,gff,circRNA source gene position
	-anno <file>	input file,source gene annotation
        -h      help

USAGE
    print $usage;
    exit;
}
