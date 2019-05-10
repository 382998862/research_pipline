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

# ----------------------------------------------------------
# GetOptions
# ----------------------------------------------------------
my ($idir,$odir,$gene_name,$cfg);

GetOptions(
        "odir:s" => \$odir,
        "idir:s" => \$idir,
        "gnfile:s"=>\$gene_name,
        "cfg:s"=>\$cfg,
        "help|?" => \&USAGE,
) or &USAGE;
&USAGE unless ($idir and $odir);

mkdir "$odir" unless (-d $odir);
$odir = &realpath($odir);
$idir = &realpath($idir);
my $out_file="$odir/ref_trans_full_table.xls";
my (%detail_cfg);
&detail_cfg_read($cfg,\%detail_cfg);

my %full_table_hash;
my @table_field;
my %id_name;
my %pathway_ko;
&log_current_time("$Script start...");



# ----------------------------------------------------------
# load gene name
# ----------------------------------------------------------
$detail_cfg{Method_RE}||="DESeq";
$detail_cfg{Method_NR}||="EBSeq";
$detail_cfg{Method_RE}=~s/DEseq/DESeq/g;
$detail_cfg{Method_NR}=~s/EBseq/EBSeq/g;
$detail_cfg{Method_RE}=~s/DEseq2/DESeq2/g;

=pod
my $fdr;
if ($detail_cfg{Method_RE} eq "edgeR" and $detail_cfg{Method_NR} eq "edgeR" and exists $detail_cfg{Pvalue}) {
    $fdr=$detail_cfg{Pvalue};
}
else{
   $fdr=$detail_cfg{FDR};
}

=cut
if (defined $gene_name) {
   open NAME, $gene_name or die $!;
   while (<NAME>) {
    chomp;
    next if (/^#/);
    my @line = split /\t/,$_;
    $id_name{$line[0]}=$line[1];
   }
   close NAME;
}


# ----------------------------------------------------------
# load expresion value
# ----------------------------------------------------------

open(FPKM,"$idir/Structure_and_Expression/DEG_Analysis/All_gene_fpkm.list")||die $!;
open(COUNT,"$idir/Structure_and_Expression/DEG_Analysis/All_gene_counts.list")||die $!;
my $h1=<FPKM>;chomp($h1);my @h1s=split(/\t/,$h1);
while(<FPKM>){
	chomp;my @tmp=split(/\t/,$_);
	for(my $i=1;$i<@tmp;$i++){
		$full_table_hash{$tmp[0]}{"$h1s[$i]_FPKM\t$h1s[$i]_count"} = $tmp[$i];
	}
}
my $h2=<COUNT>;chomp($h2);my @h2s=split(/\t/,$h2);
while(<COUNT>){
	chomp;my @tmp=split(/\t/,$_);
	for(my $i=1;$i<@tmp;$i++){
                $full_table_hash{$tmp[0]}{"$h2s[$i]_FPKM\t$h2s[$i]_count"} .= "\t$tmp[$i]";
        }
}
close(FPKM);
close(COUNT);
shift @h2s;
push @table_field, map{"$_"."_FPKM\t$_"."_count"} @h2s;

# ----------------------------------------------------------
# load DEG
# ----------------------------------------------------------
my @deg_file = glob "$idir/Structure_and_Expression/DEG_Analysis/*_vs_*";

for my $dir (@deg_file) {
    my $set_id=basename($dir);
    my $file="$dir/$set_id.DEG_final.xls";
    print "$file\n";
    basename($file) =~/(\S+_vs_\S+)\.DEG\_final\.xls/;
    my (@num)=$set_id=~/(_)/g;
    my $num=@num;
    my $soft;my $thre;my $fdr;

    if ($num==2) {
        $soft=$detail_cfg{Method_NR};
    }else{
      $soft=$detail_cfg{Method_RE};
    }

    if ($soft =~/edgeR/ and $detail_cfg{Pvalue} !=0) {
        $thre="Pvalue";
    }else{
        $thre="FDR";
    }
    $fdr=$detail_cfg{${thre}};
    push @table_field,("${set_id}_${soft}_${thre}", "${set_id}_${soft}_log2FC", "${set_id}_${soft}_(${thre}_${fdr}_FC_$detail_cfg{fold})_regulated");
	if (-f $file) {
        open (DEG, $file) or die $!;
        while (<DEG>) {
            chomp;
            next if (/^\s*$/ || /^\#/);
                next if (/^\s*$/ || $.==1 || /^\#/);
                my @lines = split /\t+/,$_;
                my $gene_id=shift@lines;
                next if (exists $full_table_hash{$gene_id}{"${set_id}_${thre}"});
                $full_table_hash{$gene_id}{"${set_id}_${soft}_${thre}"}       = $lines[-3];
                $full_table_hash{$gene_id}{"${set_id}_${soft}_log2FC"}    =  $lines[-2];
                $full_table_hash{$gene_id}{"${set_id}_${soft}_(${thre}_${fdr}_FC_$detail_cfg{fold})_regulated"} = $lines[-1];
        }
        close DEG;
    }
    if (-f "$dir/$set_id.all") {
        open (DE, dirname($file)."/$set_id.all") or die $!;
        while (<DE>) {
            chomp;
            next if (/^\s*$/ || $.==1 || /^\#/);
            my @lines = split /\t+/,$_;
			my $gene_id=shift@lines;
            next if (exists $full_table_hash{$gene_id}{"${set_id}_${soft}_${thre}"});
            $full_table_hash{$gene_id}{"${set_id}_${soft}_${thre}"}       = $lines[-2];
            $full_table_hash{$gene_id}{"${set_id}_${soft}_log2FC"}    = $lines[-1];
            $full_table_hash{$gene_id}{"${set_id}_${soft}_(${thre}_${fdr}_FC_$detail_cfg{fold})_regulated"} = "normal";
        }
        close DE;
    }
}

#print Dumper(\%full_table_hash);
# ----------------------------------------------------------
# load function annotation
# ----------------------------------------------------------
my $anno_file = "$idir/Needed_Data/Allgene_annoNseq/Allgene_Anno/Result/Integrated_Function.annotation.xls";
my @anno_field;
open (ANNO, $anno_file) or die $!;

while (<ANNO>) {
    chomp;
    next if (/^\s*$/);
    my ($gene_id, @anno_col) = (split /\t/);

    if ($.==1) {
	for my $tmp (@anno_col){
		if($tmp=~/swiss/i){
			$tmp=~s/-/_/;
			push @anno_field,$tmp;
		}
		else{
			push @anno_field,$tmp;
		}
	}
        #@anno_field = @anno_col;
        push @table_field, @anno_field;
    } else {
        for (my $i=0; $i<@anno_col; $i++) {
            $full_table_hash{$gene_id}{$anno_field[$i]} = $anno_col[$i] if (exists $full_table_hash{$gene_id});
        }
    }
}

close ANNO;

# ----------------------------------------------------------
# output 
# ----------------------------------------------------------

open (OUT, ">$out_file") or die $!;
print OUT "#ID\tgene_name\t".(join "\t",@table_field)."\n";
for my $g (sort keys %full_table_hash) {
    my $line = "$g";
    if (defined $gene_name) {
        if (exists $id_name{$g}) {
            $line.= "\t$id_name{$g}";
        }
        else
        {
            $line.= "\t--";
        }
    }
	else{
		$line.="\t--";
	}
    for my $c (@table_field) {
        if (exists $full_table_hash{$g}{$c}) {
            $line.= "\t$full_table_hash{$g}{$c}";
        } else {
            $line.= "\t--";
        }
    }
    print OUT "$line\n";
}

close OUT;

# finish log
my $elapse_time = (time()-$BEGIN_TIME)."s";
&log_current_time("$Script done. Total elapsed time: $elapse_time.");
# ----------------------------------------------------------
# sub function
# ----------------------------------------------------------
##############################################################################################################
sub log_current_time {
     # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}

##############################################################################################################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub detail_cfg_read {
    my ($cfg_file, $detail_cfg) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+$/ or /^#/);
        my ($key, $value) = (split /\s+/)[0,1];
        next unless(defined $key);
        $detail_cfg->{$key} = $value;

    }
    close CFG;

}


##############################################################################################################
sub USAGE {
    my $usage=<<"__USAGE__";
#------------------------------------------------------------------------------------
Program: $Script
Version: $version
Contact: songmm <songmm\@biomarker.com.cn>
Data: 2016-10-08
Fuction: the script is used to ...
  Usage:
         perl $Script -idir ./Analysis_Report/ -odir ./Analysis_Report/
Options:
         --idir    <string>     input dir
         --odir    <string>     output dir
	--cfg		<file>
         --gnfile  <file>       gene_name file
         --help                 shown this docment and exit

#------------------------------------------------------------------------------------
__USAGE__
    print $usage;
    exit;
}
