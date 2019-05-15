#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use Cwd qw(realpath fast_abs_path);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($target,$anno,$fa,$count,$exp,$degdir,$out,$cfg,$gene_name,$type,$rnatype);
GetOptions(
				"help|?" =>\&USAGE,
				"target:s"=>\$target,
				"anno:s"=>\$anno,
				"fa:s"=>\$fa,
				"count:s"=>\$count,
				"exp:s"=>\$exp,
				"degdir:s"=>\$degdir,
				"out:s"=>\$out,
				"cfg:s"=>\$cfg,
				"gnfile:s"=>\$gene_name,
				"type:s"=>\$type,
				) or &USAGE;
&USAGE unless ($anno and $target and $fa and $count and $exp and $degdir and $out and $cfg and $type);

$degdir = &realpath($degdir);
$count = &realpath($count);
$exp = &realpath($exp);
$target = &realpath($target);
$anno = &realpath($anno);
$fa = &realpath($fa);
$rnatype = "miRNA";
my %detail_cfg=&detail_cfg_read($cfg);

$detail_cfg{"$rnatype.Method_RE"}||="DESeq";
$detail_cfg{"$rnatype.Method_NR"}||="EBSeq";
$detail_cfg{"$rnatype.Method_RE"}=~s/DEseq/DESeq/g;
$detail_cfg{"$rnatype.Method_NR"}=~s/EBseq/EBSeq/g;
$detail_cfg{"$rnatype.Method_RE"}=~s/DEseq2/DESeq2/g;
my ($fdr,$fold,$thre);
$fold = $detail_cfg{"$rnatype.fold"};
if(!exists $detail_cfg{"$rnatype.FDR"} && !exists $detail_cfg{"FDR"} && exists $detail_cfg{"$rnatype.PValue"}){
	$fdr = $detail_cfg{"$rnatype.PValue"};
	$thre="PValue";
}else{
	$fdr = $detail_cfg{"$rnatype.FDR"};
	$thre="FDR";
}
my %gnlist=();
if(defined $gene_name ){
	open (Name,$gene_name) or die $!;
	while(<Name>){
		chomp;
		next if(/^#/);
		my @a = split /\t/,$_;
		$gnlist{$a[0]}=$a[1];
	}
	close Name;
}

my %anno;
my $anno_title;
my $anno_num;
open ANNO,$anno;
while (<ANNO>) {
	chomp;
	next if (/^\s*$/);
	if (/^#/) {
		$anno_title = (split /\t/,$_,2)[1];
		my @anno_title = split /\t/,$anno_title;
		$anno_num = @anno_title;
	}
	my ($key,$anno) = split /\t/,$_,2;
	$anno{$key} = $anno;
}
close ANNO;

#open TAROUT,">$out_dir/Target_gene_info.xls";
my %tar;
my %tar_count;
open TAR,$target;
while (<TAR>) {
	chomp;
	next if (/^#/ or /^\s*$/);
	my @b = split /\t/i,$_;
	my @tar = split /;|,/,$b[1];
	foreach my $i (sort @tar) {
#		print TAROUT "$i\n" unless ($tar_count{$i});
		$tar{$b[0]}{$i} = 1;
		$tar_count{$i} = 1;
	}
}
close TAR;
#close TAROUT;

my %fa;
my %len;
my $fa_id;
open FA,$fa;
while (<FA>) {
	chomp;
	next if (/^#/ or /^\s*$/);
	if (/^>/) {
		$fa_id = (split /\s+/,$_)[0];
		$fa_id =~ s/^>//;
		next;
	}
	$_ =~ tr/aucgtT/AUCGUU/;
	$fa{$fa_id} = $_;
	$len{$fa_id} = length $_;
}
close FA;


#All_miRNA.count.list
my %count_sample;
my %count;
my @num1;
open COUNT,$count;
while(<COUNT>){
	chomp;
	next if (/^\s*$/);
	if(/^#/){
		my @title = split /\t/,$_;
		my $title_len = @title;
		$title_len -= 1;
		for my $i (1..$title_len) {
			next if $title[$i] eq 'geneLength';
			$count_sample{$i} = $title[$i];
			push @num1,$i;
		}
		next;
	}
	my @coun = split /\t/,$_;
	foreach my $i (@num1){
		$count{$coun[0]}{$count_sample{$i}}=$coun[$i];
	}
}
close COUNT;
#open MIR,">$out_dir/miRNA_info.xls";
my %exp_sample;
my %exp;
my @num;
open EXP,$exp;
while (<EXP>) {
	chomp;
	next if (/^\s*$/);
	if (/^#/) {
		my @title = split /\t/;
		my $title_len = @title;
		$title_len -= 1;
		for my $i (1..$title_len) {
			next if $title[$i] eq 'geneLength';
			$exp_sample{$i} = $title[$i];
			push @num,$i;
		}
		next;
	}
	my @line = split /\t/;
#	print MIR "$line[0]\n";
	foreach my $i (@num) {
		$exp{$line[0]}{$exp_sample{$i}} = $line[$i];
	}
}
close EXP;


my %sample;
my %deg;
my %count1;
my %all_deg;
my @deg_group;
my %Group;
my ($soft,$condition);
foreach my $deg_file (glob "$degdir/*_vs_*/*_vs_*.DEG_final.xls") {
	open DEG,$deg_file;
	my $file_name = basename($deg_file);
	$file_name =~ /(.*)\.DEG_final.xls/;
	my $deg_group = $1;
	my @num = (split /_vs_/,$deg_group);
	my $num1 = split /_/,$num[0];
	my $num2 = split /_/,$num[1];
	my $num = $num1 + $num2;
	if($num==2){
		$soft = $detail_cfg{"$rnatype.Method_NR"};
	}else{
		$soft = $detail_cfg{"$rnatype.Method_RE"};
	}
	push @deg_group,$deg_group;
	$_ = <DEG>;
	chomp;
	my @title = split /\t/;
	my $title_len = @title;
	$title_len -= 1;
	my $j=1;
	for my $i (1..$title_len) {
		$sample{$i} = $title[$i];
		if($sample{$i} =~m/PValue|FDR/){
			$condition = $sample{$i};
	  	 	print "Your Screening threshold is: \"$condition\"\n";
		}
		
		if ($sample{$i} =~m/PValue|FDR/ || $sample{$i} eq 'log2FC') {
			my $count_key = $deg_group.'_'.$soft.'_'.$sample{$i};
			$Group{$deg_group}{$j}=$count_key;
			print "$deg_group\t$j\t$count_key\n";
			$j++;
			$count1{$count_key} = 1;
		}elsif($sample{$i} eq 'regulated'){
			my $count_key = $deg_group.'_'.$soft.'_('.$condition.'_'.$fdr.'_FC_'.$fold.')_'.$sample{$i};
			$Group{$deg_group}{$j}=$count_key;
			print "$deg_group\t$j\t$count_key\n";
			$j++;
			$count1{$count_key} = 1;
		}
	}
	print "The threshold that you set in your detail.cfg is: \n";
        print "\"software\" - \"$soft\"\t==";
        print "\"FC\" - \"$fold\"\t==";
        print "\"$thre\" - \"$fdr\"\n";
	if($soft=~m/EBSeq/i && ($condition=~m/PValue/i || $thre =~m/PValue/i)){
		print "Attention:EBSeq can not generate Pvalue,please check you parameters\n";
	}
        if($condition ne $thre){
        	print "Attention::Your Screening threshold is \"$condition\" ,which is not the same as the value \"$thre\",that you set in your detail.cfg\n";
                die;
        }
	while (<DEG>) {
		chomp;
		next if (/^\s*$/);
		my @line = split /\t/,$_;
		foreach my $i (sort {$a<=>$b} keys %sample) {
			my $key = $sample{$i};
			if ($key =~m/PValue|FDR/ || $key eq 'log2FC') {
				$key = $deg_group.'_'.$soft.'_'.$sample{$i};
				$deg{$line[0]}{$key} = $line[$i];
#				$count{$key} = 1;
			}elsif($key eq "regulated"){
				print "$key\t$line[$i]\n";
				$all_deg{$line[0]}{$deg_group} = $line[$i] if ($line[$i] eq "up"||$line[$i] eq "down");
				$key = $deg_group.'_'.$soft.'_('.$condition.'_'.$fdr.'_FC_'.$fold.')_'.$sample{$i};
				$deg{$line[0]}{$key} = $line[$i];
			}
		}
	}
	close DEG;
}


#get normal miRNA from *_vs_*.final.xls
foreach my $deg_file (glob "$degdir/*_vs_*/*_vs_*.all") {
        open DEG,$deg_file;
        my $file_name = basename($deg_file);
        $file_name =~ /(.*)\.all/;
        my $deg_group = $1;
        while (<DEG>) {
                chomp;
                next if (/^\s*$/);
                my @line = split /\t/;
                foreach my $i (sort {$a<=>$b} keys %sample) {
                        my $key = $sample{$i};
                        if ($key =~m/$condition/ || $key eq 'log2FC' ) {
                                $key = $deg_group.'_'.$soft.'_'.$sample{$i};
                                $deg{$line[0]}{$key} ||= $line[$i];
                        }
			my $key1 = $deg_group.'_'.$soft.'_('.$condition.'_'.$fdr.'_FC_'.$fold.')_regulated';
			if (!exists $deg{$line[0]}{$key1}){
				$deg{$line[0]}{$key1}= "normal";
			}
                }
        }
        close DEG;
}


#open ALLDEG,">$out_dir/All_Group_DEG_stat.xls";
#print ALLDEG "#miRNA_ID";
#foreach my $deg_group (@deg_group) {
#	print ALLDEG "\t$deg_group","_regulated";
#}
#print ALLDEG "\n";

#foreach my $miRNA_id (keys %all_deg) {
#	print ALLDEG "$miRNA_id";
#	foreach my $deg_group (@deg_group) {
#		$all_deg{$miRNA_id}{$deg_group} ||= '--';
#		print ALLDEG "\t$all_deg{$miRNA_id}{$deg_group}";
#	}
#	print ALLDEG "\n";
#}
#close ALLDEG;

my $out_title = "#miRNA_ID\tTarget_gene\tgene_name\tmiRNA_Seq\tmiRNA_Length";
foreach my $d (sort values %count_sample) {
        $out_title .= "\t$d"."_count";
	$out_title .= "\t$d"."_$type";
}
foreach my $g (sort keys %Group){
	foreach my $gg (sort keys %{$Group{$g}}){
		$out_title .= "\t$Group{$g}{$gg}";
	}
}

#foreach my $e (sort keys %count1){
#	$out_title .= "\t$e";
#}
$out_title .= "\t$anno_title\n";


#open PP,">$out_dir/ppi_qurey.ppi.cytoscapeInput.sif";
open DOUT,">$out";
print DOUT $out_title;
foreach my $id (sort keys %exp) {
	my $out1;
	my $exp = "";
	foreach my $c (sort values %exp_sample) {
                $exp .= "\t$count{$id}{$c}";
		$exp .= "\t$exp{$id}{$c}";
	}
	foreach my $g (sort keys %Group){
		foreach my $gg (sort keys %{$Group{$g}}){
			$deg{$id}{$Group{$g}{$gg}} = '--' unless (defined $deg{$id}{$Group{$g}{$gg}});
			$exp .= "\t$deg{$id}{$Group{$g}{$gg}}";
		}
	}
#	foreach my $d (sort keys %count1) {
#		$deg{$id}{$d} = '--' unless (defined $deg{$id}{$d});
#		$exp .= "\t$deg{$id}{$d}";
#	}
	if (defined $tar{$id}) {
		foreach my $tar (sort keys %{$tar{$id}}) {
			#print PP "$id\tpp\t$tar\n";
			if(exists $gnlist{$tar}){
				$out1 = "$id\t$tar\t$gnlist{$tar}\t$fa{$id}\t$len{$id}";
			}else{
				$out1 = "$id\t$tar\t--\t$fa{$id}\t$len{$id}";
			}
			if (defined $anno{$tar}) {
				$out1 .= "$exp\t$anno{$tar}\n";
				print DOUT $out1;
			}else{
				$out1 .= "$exp";
				$out1 .= "\t--"x$anno_num;
				$out1 .= "\n";
				print DOUT $out1;
			}
		}
	}else {
		$out1 = "$id\t--\t--\t$fa{$id}\t$len{$id}$exp";
		$out1 .= "\t--"x$anno_num;
		$out1 .= "\n";
		print DOUT $out1;
	}
}
close DOUT;
#close PP;




#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub detail_cfg_read {
    my $cfg_file=shift;
    my %detail_cfg;
    $rnatype = "miRNA";
    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+$/ or /^#/);
        my ($key, $value) = (split /\s+/)[0,1];
        next unless(defined $key);
        $detail_cfg{$key} = $value;

    }
    close CFG;

    my @para = ("Method_RE","Method_NR","fold","FDR","PValue");
    foreach my $para(@para){
        if(exists $detail_cfg{$para} && !exists $detail_cfg{"$rnatype.$para"}){
                $detail_cfg{"$rnatype.$para"}= $detail_cfg{$para};
        }
    }
    return %detail_cfg;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program: $Script
Version: $version
Contact: liux\@biomarker.com.cn
Description: 
 Usage:
	perl make_full_table.sRNA.pl -target mir2target.list -anno Integrated_Function.annotation.xls -fa All_miRNA.expressed.fa -exp All_miRNA_expression.list -degdir DEG_Analysis/sRNA -cfg detail.cfg -gnfile all.id_name.list -out lnc_trans_sRNA_full_table.xls -type TPM -count All_miRNA.count.list
Options:
	-target <file>		Target_Predict/*.mir2target.list,forced
	-anno <file>		Diff_Analysis/All_Anno/Result/Integrated_Function.annotation.xls,forced
	-fa <file>		Diff_Analysis/All_miRNA.expression.fa,forced
	-count <file>		All_miRNA.count.list,forced
	-exp <file>		Diff_Analysis/All_gene_expression.list,forced
	-degdir <dir>		Diff_Analysis/,forced
	-cfg <file>		detail.cfg,forced
	-gnfile <file>		gene_name.list
	-type <str>		sRNA:TPM,default
	-out <outputfile>	output file,forced
	-h         Help

USAGE
	print $usage;
	exit;
}
