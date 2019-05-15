#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;

my $version="1.0";
#######################################################################################

my $notename=`hostname`;chomp $notename;

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($annot,$index,$od);
GetOptions(
				"help|?" =>\&USAGE,
				"annot:s"=>\$annot,
				"i:s"=>\$index,
                "od=s"=>\$od,,
				) or &USAGE;
&USAGE unless ($annot and $index and $od );

$annot = &ABSOLUTE_DIR ($annot);
&MKDIR ($od);
$od = &ABSOLUTE_DIR ($od);
my %file_db=%{readconf("$Bin/../../../../db_file.cfg")};


#my $obo = "/share/nas2/database/go/gene_ontology_ext.obo";
my $obo =$file_db{go_obo_file};  # 2015-03-02 ~ 

############################################################################################################¥¶¿Ì.annot∫Õ.obo
open (AN,"$annot") or die $!;
my (%Query,%GO);

while (<AN>) {
	chomp;
	next if (/^$/);
	my ($query,$go_id)=(split/\t/,$_)[0,1];
	$Query{$query}{$go_id}=1;
	$GO{$go_id}{$query}=1;
}
close AN;

#open (OBO,$obo)||die "$!";
#$/="[Term]";
#my ($go_id,$go_name,$go_class,$anno);
#my %GO_Info;
#my %GO_anno;
#while(<OBO>){
#	chomp;
#	next if(/^$/);
#	my @Term_info=split /\n+/,$_;
#	foreach (@Term_info) {
#		if($_=~/^id: (.+)/){
#			$go_id=$1;
#		}elsif($_=~/^name: (.+)/){
#			$go_name=$1;
#		}elsif($_=~/^namespace: (.+)/){
#			my $class=$1;
#			if ($class=~/biological_process/) {
#				$go_class="Biological Process";
#			}
#			if ($class=~/cellular_component/) {
#				$go_class="Cellular Component";
#			}
#			if ($class=~/molecular_function/) {
#				$go_class="Molecular Function";
#			}
#		}elsif($_=~/^def: \"(.+)\"/){
#			$anno=$1;
#			$GO_Info{$go_id}{CLASS}=$go_class;
#			$GO_Info{$go_id}{NAME}=$go_name;
#			$GO_Info{$go_id}{ANNO}=$anno;
#			$GO_anno{$go_id}="$go_class: $go_name ($go_id);";
#		}
#	}
#}
#$/="\n";
#close OBO;
#-------------------------------- by Simon Young 2015-03-02 --------------------------------------
my (%GO_Info, %GO_anno);
open (OBO,$obo)||die "$!";
$/="[Term]";
while (<OBO>) {
    chomp;
    next if(/^$/);
    my @Term_info=split /\n+/,$_;
    my (@go_ids,$go_name,$go_class,$anno);
    foreach (@Term_info) {
        if($_=~/^id: (.+)/ or /^alt_id: (.+)/){
            push @go_ids,$1;
        }elsif($_=~/^name: (.+)/){
            $go_name=$1;
        }elsif($_=~/^namespace: (.+)/){
            my $class=$1;
            if ($class=~/biological_process/) {
                $go_class="Biological Process";
            }
            if ($class=~/cellular_component/) {
                $go_class="Cellular Component";
            }
            if ($class=~/molecular_function/) {
                $go_class="Molecular Function";
            }
        }elsif($_=~/^def: \"(.+)\"/){
            $anno=$1;
            for my $go_id (@go_ids) {
                $GO_Info{$go_id}{CLASS}=$go_class;
                $GO_Info{$go_id}{NAME}=$go_name;
                $GO_Info{$go_id}{ANNO}=$anno;
                $GO_anno{$go_id}="$go_class: $go_name ($go_id);";
            }
        }
    }
}
$/="\n";
close OBO;

my %GO_Class_stat;
my %GO_stat;
open OUT1,">$od/$index.GO.list.txt"||die"$!";
open OUT2,">$od/$index.GO.anno.txt"||die "$!";
print OUT2 "#Gene\tGO_Anno\n";
foreach my $gene (sort {$a cmp $b} keys %Query) {
	print OUT1"$gene\t";
	my @go_list=(keys %{$Query{$gene}});
	my $num=@go_list;
	print OUT2 "$gene\t$num\t";
	foreach my $go_id (sort {$a cmp $b} keys %{$Query{$gene}}) {
		if (exists $GO_anno{$go_id}) {
			print OUT1"$go_id\t";
			print OUT2 "$GO_anno{$go_id}\t";
			$GO_Class_stat{$GO_Info{$go_id}{CLASS}}++;
			$GO_stat{$GO_Info{$go_id}{CLASS}}{$GO_anno{$go_id}}{$gene}=1;
		}
	}
	print OUT1 "\n";
	print OUT2 "\n";
}
close OUT1;
close OUT2;

open OUT3,">$od/$index.GO_tree.stat.xls"||die "$!";
print OUT3"#GO_Function\tUnigene_number\tUnigene_ID\n";
foreach my $go_class (sort keys %GO_Class_stat) {
	print OUT3 "$go_class\t$GO_Class_stat{$go_class}\n";
	foreach my $go_term (sort keys %{$GO_stat{$go_class}}) {
		my @genes=keys %{$GO_stat{$go_class}{$go_term}};
		my $gene_id_str=join ";",@genes;
		my $num=@genes;
		print OUT3 "$go_term\t$num\t$gene_id_str\n";
	}
}

close (OUT3) ;

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#######################################################################################
sub file_cfg_read {
    my ($cfg_file, $file_db) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s*$/ or /^#/);
        my ($key, $value) = (split /\s+/)[0,1];
            $file_db->{$key} = $value;
    }
    close CFG;
}


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

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

################################################################
sub USAGE {#
	my $usage=<<"USAGE";
ProgramName: Extract Go Anno from annot file
Version: $version
Contact: mengf <mengf\@biomarker.com.cn>
Usage:
                 Options:
                     -i     <dir>      index of out file                         forced;
  
                     -od    <dir>      output file directory                     forced;
  
                     -annot <str>      Blat2GO annot file                        forced;

                     -h         Help

USAGE
	print $usage;
	exit;
}
