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
use newPerlBase;
##############################����д����֮ǰ��һ��д��ʱ�䡢������;������˵����ÿ���޸ĳ���ʱ��Ҳ������ע�͹�����
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

###read COG NOG eggNOG 每个分类编号对应的功能描述
my %Class=(
	"J" => [1,"Translation, ribosomal structure and biogenesis"],
	"A" => [2,"RNA processing and modification"],
	"K" => [3,"Transcription"],
	"L" => [4,"Replication, recombination and repair"],
	"B" => [5,"Chromatin structure and dynamics"],
	"D" => [6,"Cell cycle control, cell division, chromosome partitioning"],
	"Y" => [7,"Nuclear structure"],
	"V" => [8,"Defense mechanisms"],
	"T" => [9,"Signal transduction mechanisms"],
	"M" => [10,"Cell wall/membrane/envelope biogenesis"],
	"N" => [11,"Cell motility"],
	"Z" => [12,"Cytoskeleton"],
	"W" => [13,"Extracellular structures"],
	"U" => [14,"Intracellular trafficking, secretion, and vesicular transport"],
	"O" => [15,"Posttranslational modification, protein turnover, chaperones"],
	"C" => [16,"Energy production and conversion"],
	"G" => [17,"Carbohydrate transport and metabolism"],
	"E" => [18,"Amino acid transport and metabolism"],
	"F" => [19,"Nucleotide transport and metabolism"],
	"H" => [20,"Coenzyme transport and metabolism"],
	"I" => [21,"Lipid transport and metabolism"],
	"P" => [22,"Inorganic ion transport and metabolism"],
	"Q" => [23,"Secondary metabolites biosynthesis, transport and catabolism"],
	"R" => [24,"General function prediction only"],
	"S" => [25,"Function unknown"],
);

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


my $pathway_anno=glob "$anno_dir/*.pathway";
my %Anno_pathway;
if ($pathway_anno) {
		open IN,"$pathway_anno" || die;
		while (<IN>) {
			chomp;
			next if (/^$/ || /^\#/);
			my ($pathway,$pathway_id,$gene_id)=(split/\t+/,$_)[0,1,3];
			my @gene_id=split/;/,$gene_id;
			foreach my $gene(@gene_id){
				push @{$Anno_pathway{$gene}},"$pathway ($pathway_id)";
				#push @{$Anno_pathway{$gene}{pathway}},$pathway;
				#push @{$Anno_pathway{$gene}{pathway_id}},$pathway_id;
			}
		}
		close IN;
		#foreach my $gene(keys %Anno_pathway){
		#	
		#}
}

my $eggNOG=glob "$anno_dir/*.eggNOG_class.txt";
my $cog=glob "$anno_dir/*.Cog_class.txt";
my $kog=glob "$anno_dir/*.Kog_class.txt";
# cog
if ($eggNOG) {
    $Anno_Base{eggNOG}=1;

    open NOG,"$eggNOG" || die;
    while (<NOG>) {
        chomp;
        next if (/^$/ || /^\#/) ;
        my @anno=split/\t+/,$_;
        $Anno_Stat{All}{$anno[0]}=1;
        $Anno_Stat{eggNOG}{Anno}++;
        ###by xugl 201610-27  for multi clas

        my $classDes="";
        for my $class (split//,$anno[4]){
            next if $class !~/[A-Z]/;
            if(exists $Class{$class}){
                $classDes.=";; $Class{$class}[1]";
            }
            else{
                $classDes.=";; --";
            }
        }
        $classDes=~s/^;; //;
       $Anno_Gene{$anno[0]}{eggNOG}="$anno[4]\t$classDes";
        ##

        if ($Gene_length{$anno[0]}>=300 && $Gene_length{$anno[0]}<1000) {
            $Anno_Stat{eggNOG}{300}++;
        }
        if ($Gene_length{$anno[0]}>=1000) {
            $Anno_Stat{eggNOG}{1000}++;
        }
    }
    close NOG;
}
if ($cog) {
    $Anno_Base{COG}=1;

    open COG,"$cog" || die;
    while (<COG>) {
        chomp;
        next if (/^$/ || /^\#/) ;
        my @anno=split/\t+/,$_;
        $Anno_Stat{All}{$anno[0]}=1;
        $Anno_Stat{COG}{Anno}++;
        ###by xugl 201610-27  for multi clas

        my $classDes="";
        for my $class (split//,$anno[-3]){
            next if $class !~/[A-Z]/;
            if(exists $Class{$class}){
                $classDes.=";; $Class{$class}[1]";
            }
            else{
                $classDes.=";; --";
            }
        }
        $classDes=~s/^;; //;
       $Anno_Gene{$anno[0]}{COG}="$anno[-3]\t$classDes";
        ##

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
        ###by xugl 201610-27  for multi clas

        my $classDes="";
        for my $class (split//,$anno[-3]){
            next if $class !~/[A-Z]/;
            if(exists $Class{$class}){
                $classDes.=";; $Class{$class}[1]";
            }
            else{
                $classDes.=";; --";
            }
        }
        $classDes=~s/^;; //;
       $Anno_Gene{$anno[0]}{KOG}="$anno[-3]\t$classDes";
        ##

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
	if (/eggNOG/) {
		print OUT "\teggNOG_class\teggNOG_class_annotation";
	}
	elsif (/COG/) {
		print OUT "\tCOG_class\tCOG_class_annotation";
	}
    elsif (/KOG/) {
        print OUT "\tKOG_class\tKOG_class_annotation";
    }
	elsif(/KEGG/i){
		print OUT "\tKEGG_annotation\tKEGG_pathway_annotation";	
	}
	else {
		my $new_name=$_;
                $new_name=~s/Swissprot/Swiss\-Prot/i;
                $new_name=~s/nr/NR/i;
                $new_name=~s/nt/NT/i;
                $new_name=~s/Cog/COG/i;
                $new_name=~s/Kog/KOG/i;
		print OUT "\t$new_name"."_annotation";
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
		if (/eggNOG/ && !defined $Anno_Gene{$id}{eggNOG}) {
			print OUT "\t--\t--";
		}
		elsif (/COG/ && !defined $Anno_Gene{$id}{COG}) {
			print OUT "\t--\t--";
		}
        elsif (/KOG/ && !defined $Anno_Gene{$id}{KOG}) {
            print OUT "\t--\t--";
        }
		elsif(/KEGG/i && !defined $Anno_Gene{$id}{KEGG}){
			print OUT "\t--\t--";
		}
		elsif(/KEGG/i && defined $Anno_Gene{$id}{KEGG}){
			if (defined $Anno_pathway{$id}) {
				my $pathway=join(";; ",@{$Anno_pathway{$id}});
                print OUT "\t$Anno_Gene{$id}{$_}\t$pathway";
            }
			else{
				print OUT "\t$Anno_Gene{$id}{$_}\t--";
			}
            
		}
		elsif (defined $Anno_Gene{$id}{$_} and $_!~/KEGG/i) {
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
		my $new_name=$_;
                $new_name=~s/Swissprot/Swiss\-Prot/i;
                $new_name=~s/nr/NR/i;
                $new_name=~s/nt/NT/i;
                $new_name=~s/Cog/COG/i;
                $new_name=~s/Kog/KOG/i;
				$Anno_Stat{$_}{Anno}||=0;
				$Anno_Stat{$_}{300}||=0;
				$Anno_Stat{$_}{1000}||=0;
		print STAT "$new_name"."_"."Annotation\t$Anno_Stat{$_}{Anno}\t$Anno_Stat{$_}{300}\t$Anno_Stat{$_}{1000}\n";
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

sub Beta { #��linux�µ�xls�ļ�ת����windows�¿ɶ���excel�ļ�
	my ($file,$xls,$sheet,$f)=@_;
	open (IN,"$file")||die "can't open file $file\n";
	my $worksheet = $xls->add_worksheet(decode('GB2312',"$sheet"));
	my $format = $xls->add_format();
	$format->set_border($f);
	$format->set_align('left');
#	$format->set_font('Times New Roman');
	$format->set_size('11') if($f ==0);
	$format->set_font('����') if($f ==0);
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

sub Beta_x { #��linux�µ�xls�ļ�ת����windows�¿ɶ���excel�ļ�
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
	$format->set_font('����') if($f ==0);


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
