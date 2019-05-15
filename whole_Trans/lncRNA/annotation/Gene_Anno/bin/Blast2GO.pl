#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
#######################################################################################
use newPerlBase;
my %config=%{readconf("$Bin/../db_file.cfg")}; 

my $notename=`hostname`;chomp $notename;

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($od,$Index,$id,$worksh,$cpu,$step,$queues);
GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$od,
				"id:s"=>\$id,
				"sh_dir=s"=>\$worksh,
				"cpu=s"=>\$cpu,
				"queue:s"=>\$queues,
				"s=s"=>\$step,
				"k:s"=>\$Index,
				) or &USAGE;
&USAGE unless ($od and $Index and $id );

$step=$step || 1;
$id=&ABSOLUTE_DIR($id);
$cpu=$cpu || 50;
&MKDIR($od);
$od=&ABSOLUTE_DIR($od);
$worksh=$worksh || "work_sh";
&MKDIR($worksh);
$worksh=&ABSOLUTE_DIR($worksh);
$queues||="middle.q";

&MKDIR("$od/Blast2go");


my $blast2go=$config{blast2go}; #2014-12-08 ~ 
my $obo=$config{go_obo_file};


my @files=sort glob("$id/*.xml");

open (SH,">$worksh/blast2go.sh")||die "$!";
my $j=1;

foreach my $file (@files) {
	print SH " java -XX:ParallelGCThreads=5 -cp $blast2go/*:$blast2go/ext/*: es.blast2go.prog.B2GAnnotPipe -in $file -out $od/Blast2go/$Index.$j -prop $blast2go/b2gPipe.properties -v -annot \n";
	$j++;
}

close SH;

if ($step==1) {
	&qsubOrDie("$worksh/blast2go.sh",$queues,$cpu,"10G");
	$step=2;
}

if ($step==2) {
	my @stat=sort glob("$od/Blast2go/*.annot");
	foreach my $file (@stat) {
		if($file=~/$Index\.(\d+)\./){
			`cat $file >> "$od/$Index.annot"`;
		}
	}
	$step=3;
}

############################################################################################################¥¶¿Ì.annot∫Õ.obo
if ($step==3) {
	open (AN,"$od/$Index.annot") or die $!;
	my (%Query,%GO);

	while (<AN>) {
		chomp;
		next if (/^$/);
		my ($query,$go_id)=(split/\t/,$_)[0,1];
		$Query{$query}{$go_id}=1;
		$GO{$go_id}{$query}=1;
	}
	close AN;

#	open (OBO,$obo)||die "$!";
#	$/="[Term]";
#	my ($go_id,$go_name,$go_class,$anno);
#	my %GO_Info;
#	my %GO_anno;
#	while(<OBO>){
#		chomp;
#		next if(/^$/);
#		my @Term_info=split /\n+/,$_;
#		foreach (@Term_info) {
#			if($_=~/^id: (.+)/){
#				$go_id=$1;
#			}elsif($_=~/^name: (.+)/){
#				$go_name=$1;
#			}elsif($_=~/^namespace: (.+)/){
#				my $class=$1;
#				if ($class=~/biological_process/) {
#					$go_class="Biological Process";
#				}
#				if ($class=~/cellular_component/) {
#					$go_class="Cellular Component";
#				}
#				if ($class=~/molecular_function/) {
#					$go_class="Molecular Function";
#				}
#			}elsif($_=~/^def: \"(.+)\"/){
#				$anno=$1;
#				$GO_Info{$go_id}{CLASS}=$go_class;
#				$GO_Info{$go_id}{NAME}=$go_name;
#				$GO_Info{$go_id}{ANNO}=$anno;
#				$GO_anno{$go_id}="$go_class: $go_name ($go_id);";
#			}
#		}
#	}
#	$/="\n";
#	close OBO;
#-------------------------------- by Simon Young 2014-10-09 --------------------------------------
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
	open OUT1,">$od/$Index.GO.list.txt"||die"$!";
	open OUT2,">$od/$Index.GO.anno.txt"||die "$!";
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

	open OUT3,">$od/$Index.GO_tree.stat.xls"||die "$!";
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
}

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
ProgramName:
Version: $version
Contact: Guo XianChao <guoxc\@biomarker.com.cn>
Program Date:   2012.8.31
Modified Date:  2012.9.7     By mengf
Description:	This program is used to tackle xml(m7 blast) to GO Annotation;
Usage:
                 Options:
                     -id <dir>      input files directory, blast result in format -m7               forced;
  
                     -od <dir>      output file directory                                           forced;
  
                     -k  <str>       out file prefix name                                            forced;

                     -sh_dir <str>  shell dir                                                       default ./work_sh;
                     -queue  <str>  job queue 
                     -s  <int>       Step Start the program
                         1       from Blast2GO, tackle blast m7 Result          default;
                         2       from cat qsub Blast2GO's Reault file;
                         3       from GO Annotation stat ;

                     -h         Help

USAGE
	print $usage;
	exit;
}
