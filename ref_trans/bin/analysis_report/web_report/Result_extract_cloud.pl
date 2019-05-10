#!/usr/bin/perl -w
#
# Copyright (c) BMK 2012
# Writer:         mengf <mengf@biomarker.com.cn>
# Program Date:   2012.
# Modifier:       mengf <mengf@biomarker.com.cn>
# Last Modified:  2012.
use	strict;
use	Getopt::Long;
use	Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $program_name=basename($0);

my $ver="1.0";
############################################
my %opts;
GetOptions(\%opts,"id=s","od=s","species=s","cloud=s","h");
if (!defined($opts{id})||!defined($opts{od})||!defined($opts{species})||defined($opts{h})) {
	&help();
}
###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\n[$Time_Start] $Script start ... \n\n";


###############
#my $cfg=&ABSOLUTE_DIR($opts{cfg});
my $id = &ABSOLUTE_DIR($opts{id});
&MKDIR($opts{od});
my $od=&ABSOLUTE_DIR($opts{od});
my $DEG_od="$od/DEG_Analysis" ;
&MKDIR($DEG_od) if (-d "$id/Structure_and_Expression/DEG_Analysis");
my $rawdata_od="$od/rawdata";
&MKDIR($rawdata_od);
&MKDIR("$od/geneExpression");
#&MKDIR("$od/Gene_Anno");
&MKDIR("$od/NewGene");
&MKDIR("$od/NewGene/NewGene_Anno");
&MKDIR("$od/AllGene");
&MKDIR("$od/AllGene/AllGene_Anno");
&MKDIR("$od/SNP_Analysis") if (-d "$id/Structure_and_Expression/SNP_Analysis") ;
&MKDIR("$od/../Needed_Data");
&MKDIR("$od/../Needed_Data/Allgene_annoNseq");

#open LOG,"$od/log.txt" || die;
###############
#my %para;
#&para_load($cfg,\%para);

############### Extract Assembly dir
if (-d "$id/Data_Assess") {
	system "cp $id/Data_Assess/AllSample_GC_Q.stat $rawdata_od";
	system "cp $id/Data_Assess/Sample_data_relation.xls $rawdata_od";
	system "cp -r $id/Data_Assess/PNG $rawdata_od";
}
if (-d "$id/Config") {
	system "cp -r $id/Config $od/../Needed_Data";
}

if (-d "$id/Structure_and_Expression/Map_Stat/Saturation") {
	system "cp $id/Structure_and_Expression/Map_Stat/Saturation/*/*.png  $od/geneExpression";
	system "cp $id/Structure_and_Expression/Map_Stat/*.png  $od/geneExpression";
}


if (-d "$id/Structure_and_Expression/Alitsplice_Analysis") {
	&MKDIR("$od/Alt_splice");
	system "cp $id/Structure_and_Expression/Alitsplice_Analysis/*.fpkm $od/Alt_splice";
	system "cp $id/Structure_and_Expression/Alitsplice_Analysis/*.stat    $od/Alt_splice";
	system "cp $id/Structure_and_Expression/Alitsplice_Analysis/*.png    $od/Alt_splice";
}
if (-d "$id/Structure_and_Expression/Anno_Integrate/Allgene_Anno/") {
    system "cp -r $id/Structure_and_Expression/Anno_Integrate/Allgene_Anno $od/../Needed_Data/Allgene_annoNseq";
    system "rm -r $od/../Needed_Data/Allgene_annoNseq/Allgene_Anno/work_sh/" if (-d "$od/../Needed_Data/Allgene_annoNseq/Allgene_Anno/work_sh/");
	my $transcript=(glob("$od/../Needed_Data/Allgene_annoNseq/Allgene_Anno/02.gene-annotation/*.fa"))[0];
	my ($index)=split(/_/,basename$transcript);
    system "cp  $transcript $od/../Needed_Data/Allgene_annoNseq/$index.longest_transcript.fa ";
}

if (-d "$id/Structure_and_Expression/DEG_Analysis") {
	my @DEG_dir=glob "$id/Structure_and_Expression/DEG_Analysis/*";
	my %Anno;
	my %Stat;
	foreach my $dir (@DEG_dir) {
		if (-f $dir) {
			next if $dir=~/\.svg$/;
			next if $dir=~/\.log$/;
			next if $dir=~/\/used/;
			next if $dir=~/kmean_group_fpkm.list$/;
			next if $dir=~/kmeans_centers.txt$/;
			system "cp $dir $DEG_od";
		}
		if (-d $dir) {
			$dir=~m/.*\/(\S+)/;
			my $nam=$1;
			&MKDIR("$DEG_od/$nam");
			if ($dir=~/\/All_DEG$/) {
				system "cp -r $dir/* $DEG_od/$nam/";
			}
			elsif ($dir=~/\/work_sh$/) {
				`rm -r $DEG_od/$nam`;
			}
			elsif ($dir=~/\/density$/) {
				system "cp -r $dir/*.png $od/geneExpression/";
				system "cp -r $dir/*.cor $od/geneExpression/";
				system "cp -r $dir/cor_plot $od/geneExpression/";
                system "rm -r $DEG_od/density/";
			}
			elsif ($dir=~/\/htseq_count$/){
				system "cp -r $dir/*.count $DEG_od/$nam/";
			}elsif ($dir=~/\/DEG_PPI$/){
                &MKDIR("$DEG_od/DEG_PPI");
				system "cp -r $dir/*.sif  $DEG_od/DEG_PPI";
				system "cp -r $dir/*.txt  $DEG_od/DEG_PPI";
            }elsif ($dir=~/\/All_DEG$/){
				system "cp -r $dir $DEG_od/";
			}elsif ($dir=~/\/kmean_/){
				system "cp -r $dir/kmean_group_fpkm.list $DEG_od/$nam/";
				system "cp -r $dir/k-means.png $DEG_od/$nam/";
				system "cp -r $dir/kmeans_cluster.txt $DEG_od/$nam/";
			}elsif ($dir=~/_vs_/){
				system "cp -r $dir/*.DEG.final.xls* $DEG_od/$nam";
				system "cp -r $dir/*.all $DEG_od/$nam";
                system "cp -r $dir/*.png $DEG_od/$nam" unless $dir=~/.*_vs_.*_vs_.*/;


				$Stat{$nam}{up}=0;
				$Stat{$nam}{down}=0;
				$Stat{$nam}{total}=0;            
            if (-e "$dir/$nam.DEG.final.xls" && !-e "$dir/$nam.DEG.final.xls.zero") {
				system "cp -r $dir/Anno_enrichment/* $DEG_od/$nam/";
				system "cp -r $dir/DEG_Cluster $DEG_od/$nam";
				system "cp -r $dir/DEXSeqReport $DEG_od/$nam" if (-d "$dir/DEXSeqReport" );
				system "rm $DEG_od/$nam/go_enrichment/*.svg.list" if (-f "$DEG_od/$nam/go_enrichment/*.svg.list");
                system "rm $DEG_od/$nam/Cog_Anno/*.plot.log" if (-f "$DEG_od/$nam/Cog_Anno/*.plot.log");
                system "rm $DEG_od/$nam/DEG_Cluster/Rplots.pdf" if (-f "$DEG_od/$nam/DEG_Cluster/Rplots.pdf");
                system "rm $DEG_od/$nam/Graph/topGO.*";
               #system "rm -r $DEG_od/$nam/DEXSeqReport/Example" if (-d "$dir/DEXSeqReport" );

				if (-d "$dir/Graphic") {
					system "cp $dir/Graphic/*.png $DEG_od/$nam";
					system "cp $dir/Graphic/*.pdf $DEG_od/$nam";
				}
				my %Site;
				my $file=(glob "$dir/Anno_enrichment/*.annotation.xls")[0];
				open (IN,$file) or die $!;
				print $file,"\n";
				while (<IN>) {
					chomp;
					if (/^\#/) {
						my @Anno=split/\s+/,$_;
						for (my $s=0;$s<@Anno ;$s++) {
							if ($Anno[$s] eq 'COG_class') {
								$Site{'COG'}=$s;
							}
							if ($Anno[$s] eq 'KOG_class') {
								$Site{'KOG'}=$s;
							}
                                                        if ($Anno[$s] eq 'eggNOG_class') {
                                                                $Site{'eggNOG'}=$s;
                                                        }
 							if ($Anno[$s]=~/^([^_]+)_annotation/) {
								 if ($Anno[$s] eq 'Swissprot_annotation') {
								 $Site{'Swiss-Prot'}=$s;
								 next;
								 }
								 if ($Anno[$s] eq 'nt_annotation') {
								 $Site{'NT'}=$s;
								 next;
								 }
								 if ($Anno[$s] eq 'nr_annotation') {
								 $Site{'NR'}=$s;
								 next;
								 }
								$Site{$1}=$s ;
							}
						}
                        $Anno{$nam}{'Total'} = 0;
					}
					else{
						my @Info=split /\t+/,$_;
						foreach my $key (keys %Site) {
							$Anno{$nam}{$key}||=0;
							$Anno{$nam}{$key}++ if (exists $Info[$Site{$key}] and $Info[$Site{$key}]=~/\w+/) ;
						}
                        $Anno{$nam}{'Total'} ++;
					}
				}
				close IN;
				$file=(glob "$dir/*.DEG.final.xls")[0];
				open (IN,$file) or die $!;
				while (<IN>) {
					chomp;
					next if /^\#/;
					my $type=(split/\s+/,$_)[-1];
					$Stat{$nam}{up}++ if $type eq 'up';
					$Stat{$nam}{down}++ if $type eq 'down';
					$Stat{$nam}{total}++;
				}
				close IN;
			} 
        }
	}
}

my $anno_num;
	open (OUT,">$DEG_od/DEG.anno.stat") or die $!;
	my $limit_anno=0;
	foreach my $key (sort keys %Anno) {#!
		if ($limit_anno==0) {
			print OUT "#DEG Set\tTotal";
			foreach my $key1 (sort keys %{$Anno{$key}}) {
				next if ($key1 eq 'Total') ;
				print OUT "\t$key1";
			}
			print OUT "\n";
			$limit_anno++;
		}
		print OUT "$key";
		print OUT "\t$Anno{$key}{'Total'}";
		my @anno_nums=keys %{$Anno{$key}};
			$anno_num=@anno_nums;

			foreach my $key1 (sort keys %{$Anno{$key}}) {

				next if ($key1 eq 'Total') ;
				print OUT "\t$Anno{$key}{$key1}";
			}
		print OUT "\n";
	}

	foreach my $key (sort keys %Stat) {
		my @anno_nums=keys %{$Anno{$key}};
		
		if ($Stat{$key}{up}==0 and $Stat{$key}{down}==0 and $Stat{$key}{total}==0) {
			print OUT $key,"\t0" x $anno_num,"\n";
		}
	}

	close OUT;

	open (OUT,">$DEG_od/DEG.stat") or die $!;
	print OUT "DEG Set\tAll DEG\tup-regulated\tdown-regulated\n";
	foreach my $key (sort keys %Stat) {
		$Stat{$key}{up}||=0;
		$Stat{$key}{down}||=0;
		print OUT "$key\t$Stat{$key}{total}\t$Stat{$key}{up}\t$Stat{$key}{down}\n";
	}
	close OUT;
}

if (-d "$id/Structure_and_Expression/geneExpression/final_track/") {

	system "cp $id/Structure_and_Expression/geneExpression/final_track/*.newGene.longest_transcript.fa $od/NewGene";
	system "cp $id/Structure_and_Expression/geneExpression/final_track/*.newGene_final.filtered.gff $od/NewGene";
}

if (-d "$id/Structure_and_Expression/Anno_Integrate/New_Anno/Result/") {
    system "cp -r $id/Structure_and_Expression/Anno_Integrate/New_Anno/Result/* $od/NewGene/NewGene_Anno";
    system "rm -r $od/NewGene/NewGene_Anno/Blast2go" if (-d "$od/NewGene/NewGene_Anno/Blast2go");
    system "rm $od/NewGene/NewGene_Anno/*.annot"  ;
}

if (-d "$id/Structure_and_Expression/Anno_Integrate/Allgene_Anno/Result/") {
	system "cp -r $id/Structure_and_Expression/Anno_Integrate/Allgene_Anno/Result/* $od/AllGene/AllGene_Anno";
	system "cp -r $id/Structure_and_Expression/Anno_Integrate/Allgene_Anno/02.gene-annotation/*.fa  $od/AllGene/All_Gene.longest_transcript.fa";
	system "rm -r $od/AllGene/AllGene_Anno/Blast2go" if (-d "$od/AllGene/AllGene_Anno/Blast2go");
}

if(-d "$id/Structure_and_Expression/SNP_Analysis") {
	my $SNP_dir="$id/Structure_and_Expression/SNP_Analysis";

	system "cp -r $SNP_dir/stat/* $od/SNP_Analysis";
	system "cp -r $SNP_dir/All_gene.fa $od/AllGene";
	system "cp -r $SNP_dir/New_gene.fa $od/NewGene";
}


if (-d "$id/Structure_and_Expression/Gene_Structure_Optimize") {
	system "cp -r $id/Structure_and_Expression/Gene_Structure_Optimize $od";
}


if (defined $opts{cloud}) {
	system "perl $Bin/bin/get_sample.pl --input $od/ref_trans_full_table.xls --list $Bin/bin/local_and_biomaRt_database.txt --output $od/tmp --species $opts{species} --cloud";
}else {
	system "perl $Bin/bin/get_sample.pl --input $od/ref_trans_full_table.xls --list $Bin/bin/local_and_biomaRt_database.txt --output $od/tmp --species $opts{species} ";
}
if (-s "$od/tmp.xls") {
	system "perl $Bin/bin/extract_list.pl -data $od/ref_trans_full_table.xls  -list $od/tmp.xls -type combine -od $od ";
	system "mv $od/tmp.xls.list $od/ref_trans_full_table.xls";
	system "rm $od/tmp.xls";
}

system "cp $Bin/readme/readme.txt $od/";


################Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $program_name Time :[$Time_End]\n\n";
###############Subs
sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
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

sub para_load {
	my ($file,$para)=@_;
	open IN,$file||die "$!";
	while (<IN>) {
		chomp;
		s/\r+//g;
		next if(/^$/||/^\#/);
		my ($key,$info)=(split/\s+/,$_)[0,1];
		if(!$key){print "$_\n";die;}
		$para->{$key}=$info;
	}
	close IN;
}
sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub Runtime{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "\nTotal elapsed time: ${t}s\n";
}

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub help{
	print << "	Usage End.";
	Description: Extract Have_Ref Transcriptome Reaults for Html Process;
	version:$ver
	Usage:
		--id   <STR>      input dir, analysis output directory   force
		--od   <STR>      result output dir                      force
		--species  <STR>  species name for get gene name 
                          some name from  $Bin/bin/local_and_biomaRt_database.txt  force
        --cloud           analysis at biocloud
    --h           help
	Usage End.
		exit;
}
