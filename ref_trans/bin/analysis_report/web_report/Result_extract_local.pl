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
use List::Util qw/sum/;
my $program_name=basename($0);
my $readmeDir="$Bin/../../../Readme";
my $ver="1.0";
############################################
my %opts;
GetOptions(\%opts,"id=s","od=s","species=s","cloud=s","h","genename=s");
if (!defined($opts{id})||!defined($opts{od})||!defined($opts{species})||defined($opts{h})) {
	&help();
}
###############Time
my $BEGIN=time();
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\n[$Time_Start] $Script start ... \n\n";


###############
my $id = &ABSOLUTE_DIR($opts{id});
&MKDIR($opts{od});
my $od=&ABSOLUTE_DIR($opts{od});
my $DEG_od="$od/BMK_5_DEG_Analysis" ;
my @deg_sam=glob("$id/Structure_and_Expression/DEG_Analysis/*_vs_*");
&MKDIR($DEG_od) if $#deg_sam>=0;
my $rawdata_od="$od/BMK_1_rawdata";
&MKDIR($rawdata_od);
&MKDIR("$od/BMK_4_geneExpression");
&MKDIR("$od/BMK_2_NewGene");
&MKDIR("$od/BMK_2_NewGene/BMK_1_NewGene_Anno");
&MKDIR("$od/BMK_3_AllGene");
&MKDIR("$od/BMK_3_AllGene/BMK_1_AllGene_Anno");
&MKDIR("$od/BMK_7_SNP_Analysis") if (-d "$id/Structure_and_Expression/SNP_Analysis") ;
&MKDIR("$od/../Needed_Data");
&MKDIR("$od/../Needed_Data/Allgene_annoNseq");
system "cp $readmeDir/Readme_total.pdf $od/Readme.pdf";
system "cp $readmeDir/Readme_geneExpression.pdf $od/BMK_4_geneExpression/Readme.pdf";
system "cp $id/Structure_and_Expression/PPI/PPI.txt $od/../Needed_Data/";
my %config=();
if(-d "$id/Config"){
        system "cp -r $id/Config $od/../Needed_Data";
        system "mv $od/../Needed_Data/Config/*detail* $od/../Needed_Data/Config/ref_trans.detail.cfg";
        system "mv $od/../Needed_Data/Config/*data* $od/../Needed_Data/Config/ref_trans.data.cfg";
        %config=&detail_cfg_read("$od/../Needed_Data/Config/ref_trans.detail.cfg");
}
`cat $config{Ref_ann} $id/Structure_and_Expression/geneExpression/final_track/*.newGene_final.filtered.gff >$od/../Needed_Data/all.gff`;
`perl $Bin/extract_gene_position_gff.pl -i $od/../Needed_Data/all.gff -o $od`;
############### Extract Assembly dir
if (-d "$id/Data_Assess") {
	`cat $id/Data_Assess/AllSample_GC_Q.stat|cut -f 1,3,4,5,6,7,9 >$rawdata_od/AllSample_GC_Q.stat`;
	open(IN,"$id/Data_Assess/Sample_data_relation.xls")||die $!;
	open(OU,">$rawdata_od/Sample_data_relation.xls")||die $!;
	while (<IN>) {
        	if (/^\s+(.+)/) {
        		print OU "\t$1\n";
        	}else{
			print OU "$_";
		}
	}
	close IN;close OU;
	system "cp $readmeDir/Readme_rawdata.pdf $rawdata_od/Readme.pdf";
	system "cp -r $id/Data_Assess/PNG $rawdata_od";
	&MKDIR("$rawdata_od/BMK_1_PNG");
	system "cp $rawdata_od/PNG/* $rawdata_od/BMK_1_PNG/";
	system "rm -rf $rawdata_od/PNG";
	system "rm $rawdata_od/BMK_1_PNG/*cycleQ.png";
	system "cp $readmeDir/Readme_rawdata_png.pdf $rawdata_od/BMK_1_PNG/Readme.pdf";
}

if (-d "$id/Structure_and_Expression/geneExpression") {
	&MKDIR("$od/BMK_4_geneExpression/BMK_1_Mapped_Statistics");
	system "cp $id/Structure_and_Expression/Map_Stat/*.coverage.png $od/BMK_4_geneExpression/BMK_1_Mapped_Statistics";
	system "cp $id/Structure_and_Expression/Map_Stat/*.type.png $od/BMK_4_geneExpression/BMK_1_Mapped_Statistics";
	system "cp $id/Structure_and_Expression/Map_Stat/*.type.stat $od/BMK_4_geneExpression/BMK_1_Mapped_Statistics";
	system "cp $id/Structure_and_Expression/Map_Stat/*.mappedStat.xls $od/BMK_4_geneExpression/BMK_1_Mapped_Statistics";
	system "cp $readmeDir/Readme_geneExpression_MappedStatistics.pdf $od/BMK_4_geneExpression/BMK_1_Mapped_Statistics/Readme.pdf";
	&MKDIR("$od/BMK_4_geneExpression/BMK_2_Library_Assessment");
	system "cp $id/Structure_and_Expression/Map_Stat/*.insertSize.png $od/BMK_4_geneExpression/BMK_2_Library_Assessment";
	system "cp $id/Structure_and_Expression/Map_Stat/*.randcheck.png $od/BMK_4_geneExpression/BMK_2_Library_Assessment";
	system "cp $id/Structure_and_Expression/Map_Stat/Total.randchecks.png $od/BMK_4_geneExpression/BMK_2_Library_Assessment";
	system "cp $id/Structure_and_Expression/Map_Stat/*Saturation.png $od/BMK_4_geneExpression/BMK_2_Library_Assessment";
	system "cp $readmeDir/Readme_geneExpression_libraryAssesment.pdf $od/BMK_4_geneExpression/BMK_2_Library_Assessment/Readme.pdf";
	&MKDIR("$od/BMK_4_geneExpression/BMK_3_Expression_Statistics");
	system "cp $id/Structure_and_Expression/DEG_Analysis/density/all.fpkm* $od/BMK_4_geneExpression/BMK_3_Expression_Statistics" if (-d "$id/Structure_and_Expression/DEG_Analysis/density");
	system "cp $id/Structure_and_Expression/DEG_Analysis/All_gene* $od/BMK_4_geneExpression/BMK_3_Expression_Statistics" if (-d "$id/Structure_and_Expression/DEG_Analysis") ;
	system "cp $id/Structure_and_Expression/DEG_Analysis/density/*fpkm_density.png $od/BMK_4_geneExpression/BMK_3_Expression_Statistics";
	
	system "cp $readmeDir/Readme_geneExpression_Expression_Statistics.pdf $od/BMK_4_geneExpression/BMK_3_Expression_Statistics/Readme.pdf";
	&MKDIR("$od/BMK_4_geneExpression/BMK_4_Sample_Correlation");
	system "cp $id/Structure_and_Expression/DEG_Analysis/density/sample* $od/BMK_4_geneExpression/BMK_4_Sample_Correlation";
	system "cp $readmeDir/Readme_geneExpression_Sample_Correlation.pdf $od/BMK_4_geneExpression/BMK_4_Sample_Correlation/Readme.pdf";
}

if(-d "$id/Structure_and_Expression/Gene_Fusion"){
        my $medicalod="$id/Structure_and_Expression";
        my $resultpath="$od/BMK_10_Gene_Fusion";
        &MKDIR($resultpath);
        system("mkdir -p $resultpath/report")       if(!-d "$resultpath/report");
        system("mkdir -p $resultpath/png")              if(!-d "$resultpath/png");
        system("cp $medicalod/Gene_Fusion/result/png/* $resultpath/png");
        opendir(DIR,"$medicalod/Gene_Fusion/result/report")||die $!;
        my @sample=grep{$_=~/txt/} readdir(DIR);
        foreach my $s(@sample){
                system("cp $medicalod/Gene_Fusion/result/report/$s $resultpath/report/$s.xls");
        }
        closedir(DIR);
        system("cp $readmeDir/Readme_fusion.txt $resultpath/Readme.txt");
        system("cp $medicalod/Gene_Fusion/result/Fusion_gene_Pfam.anno $resultpath/Fusion_gene_Pfam.anno.xls");

}


##mark
if (-d "$id/Structure_and_Expression/Alitsplice_Analysis") {
	&MKDIR("$od/BMK_6_Alt_splice");
	system "cp $readmeDir/Readme_Altsplice.pdf $od/BMK_6_Alt_splice/Readme.pdf";
	my @fpkm=glob("$id/Structure_and_Expression/Alitsplice_Analysis/*.fpkm");
	for(my $i=0;$i<=$#fpkm;$i++){
		open(IN,$fpkm[$i])||die $!;
		my $file=(split(/\./,basename $fpkm[$i] ))[0];
		open(OU,">$od/BMK_6_Alt_splice/$file.AS.list.xls")||die $!;
		while (<IN>) {
			chomp;
			my @line=split/\s+/,$_;
			pop @line;
			my $line=join("\t",@line);
			print OU "$line\n";
        }
		close IN;close OU;
	}
	system "cp $id/Structure_and_Expression/Alitsplice_Analysis/*.stat    $od/BMK_6_Alt_splice";
	system "cp $id/Structure_and_Expression/Alitsplice_Analysis/*.png    $od/BMK_6_Alt_splice";
}
if (-d "$id/Structure_and_Expression/Anno_Integrate/Allgene_Anno/") {
    system "cp -r $id/Structure_and_Expression/Anno_Integrate/Allgene_Anno $od/../Needed_Data/Allgene_annoNseq";
    system "rm -r $od/../Needed_Data/Allgene_annoNseq/Allgene_Anno/work_sh/" if (-d "$od/../Needed_Data/Allgene_annoNseq/Allgene_Anno/work_sh/");
	my $transcript=(glob("$od/../Needed_Data/Allgene_annoNseq/Allgene_Anno/02.gene-annotation/*.fa"))[0];
	my ($index)=split(/_/,basename$transcript);
    system "cp  $transcript $od/../Needed_Data/Allgene_annoNseq/$index.longest_transcript.fa ";
}

####################DEG analysis
my @all_deg_sam=glob("$id/Structure_and_Expression/DEG_Analysis/*_vs_*");
if ($#all_deg_sam>=0) {
	system "cp $readmeDir/Readme_DEGanalysis.pdf $DEG_od/Readme.pdf";
	&MKDIR("$DEG_od/BMK_1_All_DEG");
	&MKDIR("$DEG_od/BMK_2_DEG_PPI");
	system "cp $readmeDir/Readme_DEGanalysis_AllDEG.pdf $DEG_od/BMK_1_All_DEG/Readme.pdf";
	my $DEG_dir="$id/Structure_and_Expression/DEG_Analysis";
	`cp -r $DEG_dir/All_DEG/* $DEG_od/BMK_1_All_DEG`;
	`cp -r $DEG_dir/DEG_PPI/* $DEG_od/BMK_2_DEG_PPI && rm $DEG_od/BMK_2_DEG_PPI/used_gene.list`;
	foreach my $deg(@all_deg_sam){
		my $dir=dirname $deg;
		my $vs=basename $deg;	
		&MKDIR("$DEG_od/BMK_3_$vs");
		`cp $readmeDir/GSEA.pdf $DEG_od/BMK_3_$vs/Readme.pdf`;
		&MKDIR("$DEG_od/BMK_3_$vs/BMK_1_Statistics_Visualization");
		`cp $deg/$vs*.p* $DEG_od/BMK_3_$vs/BMK_1_Statistics_Visualization`;
		`cp $deg/$vs*.xls $DEG_od/BMK_3_$vs/BMK_1_Statistics_Visualization`;
		`cp $deg/$vs.all $DEG_od/BMK_3_$vs/BMK_1_Statistics_Visualization/$vs.all.xls`;
		&MKDIR("$DEG_od/BMK_3_$vs/BMK_2_Anno_enrichment");
		&MKDIR("$DEG_od/BMK_3_$vs/BMK_2_Anno_enrichment/BMK_1_Anno");
		`cp $deg/Anno_enrichment/anno/* $DEG_od/BMK_3_$vs/BMK_2_Anno_enrichment/BMK_1_Anno`;
		&MKDIR("$DEG_od/BMK_3_$vs/BMK_2_Anno_enrichment/BMK_2_Enrichment");
		`cp $deg/Anno_enrichment/enrich/* $DEG_od/BMK_3_$vs/BMK_2_Anno_enrichment/BMK_2_Enrichment`;
		`rm $DEG_od/BMK_3_$vs/BMK_2_Anno_enrichment/BMK_2_Enrichment/*_GSEA.xls`;
		`cp -r $deg/Anno_enrichment/pathway/kegg_map $DEG_od/BMK_3_$vs/BMK_2_Anno_enrichment/BMK_3_KEGG_map`;	
		&MKDIR("$DEG_od/BMK_3_$vs/BMK_3_GSEA");
		`cp -r $deg/Anno_enrichment/enrich/*_GSEA.xls $DEG_od/BMK_3_$vs/BMK_3_GSEA`;
		`cp -r $deg/Anno_enrichment/enrich/GSEA $DEG_od/BMK_3_$vs/BMK_3_GSEA/result`;
	}
	&getDEG("$id/Structure_and_Expression/DEG_Analysis","$DEG_od/DEG.stat");
	&DEG_anno_stat("$id/Structure_and_Expression/DEG_Analysis","Anno_enrichment","$DEG_od/DEG.anno.stat");
}
#########TFs analysis
#########TFBS_Analysis
if (-d "$id/Structure_and_Expression/DEG_Analysis/TFBS_Analysis") {
	my $tfdir="$DEG_od/BMK_4_TF_Analysis/BMK_1_TFBS_Analysis";
	`mkdir -p $tfdir`;
	`mkdir $tfdir/each_DEgeneRes` unless(-d "$tfdir/each_DEgeneRes");
	if(-e "$id/Structure_and_Expression/DEG_Analysis/TFBS_Analysis/newGenes_TFBS_predictRes.txt"){
		`cp $id/Structure_and_Expression/DEG_Analysis/TFBS_Analysis/newGenes_TFBS_predictRes.txt $tfdir/newGenes_TFBS_predictRes.xls`;
	}elsif(-e "$id/Structure_and_Expression/DEG_Analysis/TFBS_Analysis/allGenes_TFBS_predictRes.txt"){
		`cp $id/Structure_and_Expression/DEG_Analysis/TFBS_Analysis/allGenes_TFBS_predictRes.txt $tfdir/allGenes_TFBS_predictRes.xls`;
	}
	`mkdir $tfdir/DEG_seqLogo` unless(-d "$tfdir/DEG_seqLogo");
	`cp $id/Structure_and_Expression/DEG_Analysis/TFBS_Analysis/DEG_seqLogo/*.p* $tfdir/DEG_seqLogo`;
	`cp $id/Structure_and_Expression/DEG_Analysis/TFBS_Analysis/each_DEgeneRes/* $tfdir/each_DEgeneRes`;
	`cp $readmeDir/Readme_TFBS_analysis.pdf $tfdir/Readme.pdf`;
}
######TF_activity
if (-d "$id/Structure_and_Expression/DEG_Analysis/TF_activity") {
	my $tfdir="$DEG_od/BMK_4_TF_Analysis/BMK_2_TF_activity";
	`mkdir -p $tfdir`;
	my $dir="$id/Structure_and_Expression/DEG_Analysis/TF_activity";
	`cp $dir/TFs_network.p* $tfdir`;
	`cp $dir/TFs_influences_heatmap.p* $tfdir`;
	`cp $dir/TFs_cornet.txt $tfdir/TFs_cornet.xls`;
	`cp $dir/TFs_activity_grn.txt $tfdir/TFs_activity_grn.xls`;
	`cp $dir/TFs_influence.txt $tfdir/TFs_influence.xls`;
	`cp $readmeDir/Readme_TF_activity.pdf $tfdir/Readme.pdf`;
}
#######differential_AS_analysis
my @diff_AS=glob("$id/Structure_and_Expression/differential_AS_analysis/*_vs_*");
if(@diff_AS>0) {
	foreach my $diffas(@diff_AS){
		my $vs=basename($diffas);
		&MKDIR("$DEG_od/BMK_3_$vs/BMK_5_diff_AS_analysis");
		`sed -i 's/\"//g' $diffas/*.JC.txt `;
		`cut -f2,4-23 $diffas/A3SS.MATS.JC.txt > $DEG_od/BMK_3_$vs/BMK_5_diff_AS_analysis/A3SS.MATS.JC.xls`;
		`cut -f2,4-23 $diffas/A5SS.MATS.JC.txt > $DEG_od/BMK_3_$vs/BMK_5_diff_AS_analysis/A5SS.MATS.JC.xls`;
		`cut -f2,4-23 $diffas/MXE.MATS.JC.txt > $DEG_od/BMK_3_$vs/BMK_5_diff_AS_analysis/MXE.MATS.JC.xls`;
		`cut -f2,4-23 $diffas/RI.MATS.JC.txt > $DEG_od/BMK_3_$vs/BMK_5_diff_AS_analysis/RI.MATS.JC.xls`;
		`cut -f2,4-23 $diffas/SE.MATS.JC.txt > $DEG_od/BMK_3_$vs/BMK_5_diff_AS_analysis/SE.MATS.JC.xls`;
		`cp $readmeDir/Readme_DEGanalysis_rMATS.pdf $DEG_od/BMK_3_$vs/BMK_5_diff_AS_analysis/Readme.pdf`;
	}
	&getDiff_AS("$id/Structure_and_Expression/differential_AS_analysis","$DEG_od/Diff.AS.stat");
}
my @deus=glob("$id/Structure_and_Expression/DEU_Analysis/*_vs_*");
if(scalar(@deus)>0){
	foreach my $deu(@deus){
		my $vs=basename $deu;
		`mkdir $DEG_od/BMK_3_$vs/BMK_4_DEXSeqReport`	unless(-d "$DEG_od/BMK_3_$vs/BMK_4_DEXSeqReport");
		`cp $deu/DEU_Result_All.xls $DEG_od/BMK_3_$vs/BMK_4_DEXSeqReport/${vs}DEU.all`;
		`cp $deu/DEU_Result_Final.xls $DEG_od/BMK_3_$vs/BMK_4_DEXSeqReport/${vs}DEU.final.xls`;
		`cp $deu/DEXSeqReport/testForDEU.html $DEG_od/BMK_3_$vs/BMK_4_DEXSeqReport/testForDEU.html`;
		`cp -r $deu/DEXSeqReport/files $DEG_od/BMK_3_$vs/BMK_4_DEXSeqReport`;
		`cp $readmeDir/Readme_DEGanalysis_sample_DEX.pdf $DEG_od/BMK_3_$vs/BMK_4_DEXSeqReport/Readme.pdf`;
	}
}

sub getDEG{
        my ($dir,$o)=@_;
        opendir(DIR,$dir)||die $!;
        my @vs=grep{/_vs_/ && -d "$dir/$_"}readdir(DIR);
        closedir(DIR);
        open(STAT,">$o")||die $!;
        print STAT "DEG Set\tAll DEG\tup-regulated\tdown-regulated\n";
        for(my $i=0;$i<@vs;$i++){
                my $all=`wc -l $dir/$vs[$i]/$vs[$i].DEG_final.xls`;
                chomp($all);$all=(split(/\s+/,$all))[0]-1;
                my $up=`awk '{print \$NF}' $dir/$vs[$i]/$vs[$i].DEG_final.xls|grep up|wc -l`;
                chomp($up);$up=(split(/\s+/,$up))[0];
                my $down=$all-$up;
                print STAT "$vs[$i]\t$all\t$up\t$down\n";
        }
        close(STAT);
}
sub getDiff_AS{
	my ($dir,$o)=@_;
	opendir(DIR,$dir)||die $!;
	my @vs=grep{/_vs_/ && -d "$dir/$_"}readdir(DIR);
	closedir(DIR);
	open(STAT,">$o")||die $!;
	print STAT "DEG Set\tA3SS\tA5SS\tMXE\tRI\tSE\n";
	for(my $i=0;$i<@vs;$i++){
		my $A3SS=`wc -l $dir/$vs[$i]/A3SS.MATS.JC.txt`;
		$A3SS=(split(/\s+/,$A3SS))[0]-1;
		my $A5SS=`wc -l $dir/$vs[$i]/A5SS.MATS.JC.txt`;
		$A5SS=(split(/\s+/,$A5SS))[0]-1;
		my $MXE=`wc -l $dir/$vs[$i]/MXE.MATS.JC.txt`;
		$MXE=(split(/\s+/,$MXE))[0]-1;
		my $RI=`wc -l $dir/$vs[$i]/RI.MATS.JC.txt`;
		$RI=(split(/\s+/,$RI))[0]-1;
		my $SE=`wc -l $dir/$vs[$i]/SE.MATS.JC.txt`;
		$SE=(split(/\s+/,$SE))[0]-1;
		print STAT "$vs[$i]\t$A3SS\t$A5SS\t$MXE\t$RI\t$SE\n";
	}
	close(STAT);
}
sub DEG_anno_stat{
        my ($dir,$path,$out)=@_;        #L01_vs_L11/Anno_enrichment/anno/L01_vs_L11.annotation.xls
        my @annos=("COG_class_annotation","GO_annotation","KEGG_annotation","KOG_class_annotation","Swiss-Prot_annotation","eggNOG_class_annotation","NR_annotation");
        my @files=glob("$dir/*_vs_*/$path/anno/*.annotation.xls");
 	open(OUT,">$out")||die $!;
	print OUT "#DEG Set\tAnnotated\t",join("\t",map{(split(/_annotation/,$_))[0]}@annos),"\n";
        if(@files>0){
                foreach my $file(@files){
                        my %info=();my %genes=();
                        open(ANNO,$file)||die $!;
                        my $head=<ANNO>;chomp($head);my @header=split(/\t/,$head);
                        my %index=();
                        for(my $i=0;$i<@header;$i++){
                                foreach my $a(@annos){
                                        $index{$a}=$i if($header[$i] eq $a);
                                }
                        }
                        while(<ANNO>){
                                chomp;next if($_=~/^#/);
                                my @tmp=split(/\t/,$_);
                                for(my $i=0;$i<@annos;$i++){
                                        if($tmp[$index{$annos[$i]}] ne "--"){
                                                $info{$annos[$i]}{$tmp[0]}++;
                                                $genes{$tmp[0]}++;
                                        }
                                }
                        }
                        close(ANNO);
                        my @values=();
                        for(my $i=0;$i<@annos;$i++){
                                push @values, scalar(keys %{$info{$annos[$i]}});
                        }
                        my $base=basename $file;$base=(split(/\./,$base))[0];
                        print OUT "$base\t",scalar(keys %genes),"\t",join("\t",@values),"\n";
                }
        }
	close(OUT);
}


############mark#################
if (-d "$id/Structure_and_Expression/geneExpression/final_track/") {
	system "cp $readmeDir/Readme_newgene.pdf $od/BMK_2_NewGene/Readme.pdf";
	system "cp $id/Structure_and_Expression/geneExpression/final_track/*.newGene.longest_transcript.fa $od/BMK_2_NewGene/";
	system "cp $id/Structure_and_Expression/geneExpression/final_track/*.newGene_final.filtered.gff $od/BMK_2_NewGene/";
}

if (-d "$id/Structure_and_Expression/Anno_Integrate/New_Anno/Result/") {
	system "cp $readmeDir/Readme_newgene_anno.pdf $od/BMK_2_NewGene/BMK_1_NewGene_Anno/Readme.pdf";
    system "cp -r $id/Structure_and_Expression/Anno_Integrate/New_Anno/Result/* $od/BMK_2_NewGene/BMK_1_NewGene_Anno/";
    system "rm -r $od/BMK_2_NewGene/BMK_1_NewGene_Anno/Blast2go" if (-d "$od/NewGene/NewGene_Anno/Blast2go");
    system "rm $od/BMK_2_NewGene/BMK_1_NewGene_Anno/*.annot";
   my @temp=glob("$od/BMK_2_NewGene/BMK_1_NewGene_Anno/Go_Mid_File*");
    system "rm -r $od/BMK_2_NewGene/BMK_1_NewGene_Anno/Go_Mid_File*"  if $#temp>=0;
    my @file=glob("$od/BMK_2_NewGene/BMK_1_NewGene_Anno/*");
	foreach my $file(@file){
		if (-f $file) {
            my $newname=$file;
			#$newname=~s/_Unigene\.GO_enrichment\.stat\.xls/_Unigene\.GO\.stat\.xls/;
			$newname=~s/newGene\.longest_transcript\.fa\.GO_enrichment\.stat\.xls/newGene\.longest_transcript\.fa\.GO\.stat\.xls/;
			$newname=~s/nr/NR/i;
			$newname=~s/nt/NT/i;
			$newname=~s/Cog/COG/i;
			$newname=~s/Kog/KOG/i;
			$newname=~s/Swissprot/Swiss\-Prot/i;
			system "mv $file $newname" unless ($file eq $newname);
			if ($file=~/newGene\.longest_transcript\.fa\.GO\.anno\.txt/) {
                open(IN_GENE_A,"$file")||die $!;
				open(OUT_GENE_A, ">temp_file")||die $!;
				while (<IN_GENE_A>) {
                    chomp;
					if (/^#/) {
                        print OUT_GENE_A "#Gene\tNumber\tGO_Anno\n";
                    }
					else{
						print OUT_GENE_A "$_\n";
					}
                }
				close IN_GENE_A;close OUT_GENE_A;
				system "rm $file";
				system "mv temp_file $file";
            }
        }
	}
	#system "rm $od/BMK_2_NewGene/BMK_1_NewGene_Anno/Function_Annotation.stat.xls";Integrated_Function.annotation.xls
	system "rm $od/BMK_2_NewGene/BMK_1_NewGene_Anno/INTegrated_Function.annotation.xls";
	system "rm $od/BMK_2_NewGene/BMK_1_NewGene_Anno/*.newGene.longest_transcript.fa.GO.list.txt";
	system "rm $od/BMK_2_NewGene/BMK_1_NewGene_Anno/*.newGene.longest_transcript.fa.GO_tree.stat.xls";
	system "rm $od/BMK_2_NewGene/BMK_1_NewGene_Anno/*.newGene.longest_transcript.fa.Kegg.globe";
	system "rm $od/BMK_2_NewGene/BMK_1_NewGene_Anno/*.newGene.longest_transcript.fa.Kegg.ko";
	system "rm $od/BMK_2_NewGene/BMK_1_NewGene_Anno/*.newGene.longest_transcript.fa.Kegg.pathway";
	my @tmp=glob("$od/BMK_2_NewGene/BMK_1_NewGene_Anno/*.svg");
	system "rm $od/BMK_2_NewGene/BMK_1_NewGene_Anno/*.svg" if ($#tmp>=0);
	
	&MKDIR("$od/BMK_2_NewGene/BMK_1_NewGene_Anno/BMK_1_PNG");
	system "mv $od/BMK_2_NewGene/BMK_1_NewGene_Anno/*stat $od/BMK_2_NewGene/BMK_1_NewGene_Anno/BMK_1_PNG/";
	system "mv $od/BMK_2_NewGene/BMK_1_NewGene_Anno/*stat.xls $od/BMK_2_NewGene/BMK_1_NewGene_Anno/BMK_1_PNG/";
	system "mv $od/BMK_2_NewGene/BMK_1_NewGene_Anno/BMK_1_PNG/Function_Annotation.stat.xls $od/BMK_2_NewGene/BMK_1_NewGene_Anno/";
	system "mv $od/BMK_2_NewGene/BMK_1_NewGene_Anno/*png $od/BMK_2_NewGene/BMK_1_NewGene_Anno/BMK_1_PNG/";
	system "mv $od/BMK_2_NewGene/BMK_1_NewGene_Anno/*pdf $od/BMK_2_NewGene/BMK_1_NewGene_Anno/BMK_1_PNG/";
	system "cp $readmeDir/Readme_newgene_anno_png.pdf $od/BMK_2_NewGene/BMK_1_NewGene_Anno/BMK_1_PNG/Readme.pdf";
	system "cp $readmeDir/Readme_newgene_anno.pdf $od/BMK_2_NewGene/BMK_1_NewGene_Anno/Readme.pdf";
	&MKDIR("$od/BMK_2_NewGene/BMK_1_NewGene_Anno/BMK_2_KEGG_map");
	system "cp $od/BMK_2_NewGene/BMK_1_NewGene_Anno/Kegg_map/* $od/BMK_2_NewGene/BMK_1_NewGene_Anno/BMK_2_KEGG_map/";
	system "rm -rf $od/BMK_2_NewGene/BMK_1_NewGene_Anno/Kegg_map";
	system "rm -rf $od/BMK_2_NewGene/BMK_1_NewGene_Anno/Blast2go";
}

if (-d "$id/Structure_and_Expression/Anno_Integrate/Allgene_Anno/Result/") {
	system "cp $readmeDir/readme_Allgene.pdf $od/BMK_3_AllGene/Readme.pdf";
	system "cp $readmeDir/Readme_allgene_anno.pdf $od/BMK_3_AllGene/BMK_1_AllGene_Anno/Readme.pdf";
	system "cp -r $id/Structure_and_Expression/Anno_Integrate/Allgene_Anno/02.gene-annotation/*.fa  $od/BMK_3_AllGene/All_Gene.longest_transcript.fa";
    system "cp -r $id/Structure_and_Expression/Anno_Integrate/Allgene_Anno/Result/* $od/BMK_3_AllGene/BMK_1_AllGene_Anno/";
	system "rm -r $od/BMK_3_AllGene/BMK_1_AllGene_Anno/Blast2go" if (-d "$od/BMK_3_AllGene/BMK_1_AllGene_Anno/Blast2go");
	my @file=glob("$od/BMK_3_AllGene/BMK_1_AllGene_Anno/*");
	foreach my $file(@file){
		if (-f $file) {
            my $newname=$file;
			$newname=~s/_Unigene\.GO_enrichment\.stat\.xls/_Unigene\.Result_extract_cloudGO\.stat\.xls/;
			$newname=~s/_Unigene\.nr/_Unigene\.NR/i;
			$newname=~s/_Unigene\.Cog/_Unigene\.COG/i;
			$newname=~s/_Unigene\.Kog/_Unigene\.KOG/i;
			$newname=~s/_Unigene\.Swissprot/_Unigene\.Swiss\-Prot/i;
			$newname=~s/Unigene/Gene/i;
			system "mv $file $newname" unless ($file eq $newname);
			if ($file=~/Gene\.GO\.anno\.txt/) {
                open(IN_GENE_A,"$file")||die $!;
				open(OUT_GENE_A, ">temp_file")||die $!;
				while (<IN_GENE_A>) {
                    chomp;
					if (/^#/) {
                        print OUT_GENE_A "#Gene\tNumber\tGO_Anno\n";
                    }
					else{
						print OUT_GENE_A "$_\n";
					}
                }
				close IN_GENE_A;close OUT_GENE_A;
				system "rm $file";
				system "mv temp_file $file";
            }
			if ($file=~/Gene\.GO\.stat\.txt/) {
                open(IN_GENE_A,"$file")||die $!;
				open(OUT_GENE_A, ">temp_file")||die $!;
				while (<IN_GENE_A>) {
                    chomp;
					$_=~s/Unigene/Gene/;
					print OUT_GENE_A "$_\n";
					
                }
				close IN_GENE_A;close OUT_GENE_A;
				system "rm $file";
				system "mv temp_file $file";
            }
        }
	}
    system "rm $od/BMK_3_AllGene/BMK_1_AllGene_Anno/Function_Annotation.stat.xls";
	system "rm $od/BMK_3_AllGene/BMK_1_AllGene_Anno/Integrated_Function.annotation.xls";
	system "rm $od/BMK_3_AllGene/BMK_1_AllGene_Anno/*Gene.GO.list.txt";
	system "rm $od/BMK_3_AllGene/BMK_1_AllGene_Anno/*Gene.GO_tree.stat.xls";
	system "rm $od/BMK_3_AllGene/BMK_1_AllGene_Anno/*.Kegg.globe";
	system "rm $od/BMK_3_AllGene/BMK_1_AllGene_Anno/*.Kegg.ko";
	system "rm $od/BMK_3_AllGene/BMK_1_AllGene_Anno/*.Kegg.pathway";
	my @tmp=glob("$od/BMK_3_AllGene/BMK_1_AllGene_Anno/*.svg");
	system "rm $od/BMK_3_AllGene/BMK_1_AllGene_Anno/*.svg" if ($#tmp>=0);
	
	&MKDIR("$od/BMK_3_AllGene/BMK_1_AllGene_Anno/BMK_1_PNG");
	system "mv $od/BMK_3_AllGene/BMK_1_AllGene_Anno/*stat $od/BMK_3_AllGene/BMK_1_AllGene_Anno/BMK_1_PNG/";
	system "mv $od/BMK_3_AllGene/BMK_1_AllGene_Anno/*stat.xls $od/BMK_3_AllGene/BMK_1_AllGene_Anno/BMK_1_PNG/";
	system "mv $od/BMK_3_AllGene/BMK_1_AllGene_Anno/*png $od/BMK_3_AllGene/BMK_1_AllGene_Anno/BMK_1_PNG/";
	system "mv $od/BMK_3_AllGene/BMK_1_AllGene_Anno/*pdf $od/BMK_3_AllGene/BMK_1_AllGene_Anno/BMK_1_PNG/";
	system "cp $readmeDir/Readme_allgene_anno_png.pdf $od/BMK_3_AllGene/BMK_1_AllGene_Anno/BMK_1_PNG/Readme.pdf";
	&MKDIR("$od/BMK_3_AllGene/BMK_1_AllGene_Anno/BMK_2_KEGG_map");
	system "cp $od/BMK_3_AllGene/BMK_1_AllGene_Anno/Kegg_map/* $od/BMK_3_AllGene/BMK_1_AllGene_Anno/BMK_2_KEGG_map/";
	system "rm -rf $od/BMK_3_AllGene/BMK_1_AllGene_Anno/Kegg_map";
	system "cp $readmeDir/Readme_allgene_anno.pdf $od/BMK_3_AllGene/BMK_1_AllGene_Anno/Readme.pdf";
}

if(-d "$id/Structure_and_Expression/SNP_Analysis") {
	my $SNP_dir="$id/Structure_and_Expression/SNP_Analysis";
	system "cp -r $SNP_dir/stat/* $od/BMK_7_SNP_Analysis/";
	system "cp -r $SNP_dir/All_gene.fa $od/BMK_3_AllGene/";
	system "cp -r $SNP_dir/New_gene.fa $od/BMK_2_NewGene/";
	system "rm $od/BMK_7_SNP_Analysis/final.indel.anno.gatk.vcf" if (-f "$od/BMK_7_SNP_Analysis/final.indel.anno.gatk.vcf");
	system "rm $od/BMK_7_SNP_Analysis/final.snp.anno.gatk.vcf" if (-f "$od/BMK_7_SNP_Analysis/final.snp.anno.gatk.vcf");
	system "cp $readmeDir/Readme_SNP.pdf $od/BMK_7_SNP_Analysis/Readme.pdf";
	&MKDIR("$od/BMK_7_SNP_Analysis/BMK_2_SNP_anno");
	system "cp $od/BMK_7_SNP_Analysis/snp_anno/* $od/BMK_7_SNP_Analysis/BMK_2_SNP_anno/";
	system "rm -rf $od/BMK_7_SNP_Analysis/snp_anno";
	system "cp $readmeDir/Readme_SNPanno.pdf $od/BMK_7_SNP_Analysis/BMK_2_SNP_anno/Readme.pdf";
	my @gatk=glob("$od/BMK_7_SNP_Analysis/BMK_2_SNP_anno/final.snp.anno.gatk.*.list");
	system "rm $od/BMK_7_SNP_Analysis/BMK_2_SNP_anno/final.snp.anno.gatk.*.list" if (@gatk>0);
	&MKDIR("$od/BMK_7_SNP_Analysis/BMK_1_InDel_anno");
	system "cp $od/BMK_7_SNP_Analysis/indel_anno/* $od/BMK_7_SNP_Analysis/BMK_1_InDel_anno/";
	system "rm -rf $od/BMK_7_SNP_Analysis/indel_anno";
	system "cp $readmeDir/Readme_Indelanno.pdf $od/BMK_7_SNP_Analysis/BMK_1_InDel_anno/Readme.pdf";
	@gatk=glob("$od/BMK_7_SNP_Analysis/BMK_1_InDel_anno/final.indel.anno.gatk.*.list");
	system "rm $od/BMK_7_SNP_Analysis/BMK_1_InDel_anno/final.indel.anno.gatk.*.list" if(@gatk>0);
	&MKDIR("$od/BMK_7_SNP_Analysis/BMK_3_SNP_type");
	system "cp $readmeDir/Readme_SNPtype.pdf $od/BMK_7_SNP_Analysis/BMK_3_SNP_type/Readme.pdf";
	system("mv $od/BMK_7_SNP_Analysis/*snp.type* $od/BMK_7_SNP_Analysis/BMK_3_SNP_type/");
	system("mv $od/BMK_7_SNP_Analysis/AllSample.snp.stat $od/BMK_7_SNP_Analysis/AllSample.SNP.stat");
	system("mv $od/BMK_7_SNP_Analysis/All.snp_type.stat $od/BMK_7_SNP_Analysis/BMK_3_SNP_type/All.SNP_type.stat");
	my @file=glob("$od/BMK_7_SNP_Analysis/BMK_2_SNP_anno/*");
	my @file2=glob("$od/BMK_7_SNP_Analysis/*.list");
	push @file,@file2;
	foreach my $file(@file){
		my $newname=$file;
		$newname=~s/snp/SNP/g;
		$newname=~s/indel/InDel/g;
		system "mv $file $newname" unless ($file eq $newname);
	}
	@file=glob("$od/BMK_7_SNP_Analysis/BMK_3_SNP_type/*");
	foreach my $file(@file){
		my $newname=$file;
		$newname=~s/snp/SNP/g;
		$newname=~s/indel/InDel/g;
		system "mv $file $newname" unless ($file eq $newname);
	}
	@file=glob("$od/BMK_7_SNP_Analysis/BMK_1_InDel_anno/*");
	foreach my $file(@file){
		my $newname=$file;
		$newname=~s/snp/SNP/g;
		$newname=~s/indel/InDel/g;
		system "mv $file $newname" unless ($file eq $newname);
	}
	@file=glob("$od/BMK_7_SNP_Analysis/*");
	foreach my $file(@file){
		if (-f $file) {
        		my $newname=$file;
			$newname=~s/snp/SNP/g;
			$newname=~s/indel/InDel/g;
			system "mv $file $newname" unless ($file eq $newname);
        	}
	}
}

if (-d "$id/Structure_and_Expression/Gene_Structure_Optimize") {
	&MKDIR("$od/BMK_8_Gene_Structure_Optimize");
	system "cp $readmeDir/Readme_Gene_stru_optimize.pdf $od/BMK_8_Gene_Structure_Optimize/Readme.pdf";
	system "cp -r $id/Structure_and_Expression/Gene_Structure_Optimize/*.geneStructure.optimize.xls $od/BMK_8_Gene_Structure_Optimize/";
}
my $path=dirname $opts{genename};
my $genefile="$path/id_name.list";
if(-f $genefile){
	&cmd_call("perl $Bin/bin/make_full_table.pl --idir $id/ --odir $od/ --gnfile $genefile --cfg $od/../Needed_Data/Config/ref_trans.detail.cfg") ;
}else{
	print "no file $genefile\n";
	&cmd_call("perl $Bin/bin/make_full_table.pl --idir $id/ --odir $od/ --cfg $od/../Needed_Data/Config/ref_trans.detail.cfg") ;
}
system("cp $od/ref_trans_full_table.xls $od/ref_trans_full_table.xls.bak");
&cmd_call("perl $Bin/bin/Go_class.pl -table $od/ref_trans_full_table.xls.bak -o $od/ref_trans_full_table.xls") ;
system "rm $od/ref_trans_full_table.xls.bak $od/ref_trans_full_table.xls.temp";

my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd $program_name Time :[$Time_End]\n\n";
###############Subs
sub cmd_call {
	print "@_\n";
	system(@_) == 0 or die "system @_ failed: $?";
}
sub detail_cfg_read {
    my ($cfg_file) = @_;
    my %detail_cfg=();
    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+/ or /^#/);
        my ($key, $value) = (split /\s+/)[0,1];
	next unless(defined $key and $value);
        $detail_cfg{$key} = $value;
        if ($key eq 'Project_name' or $key eq 'Customer_info' or $key eq 'Project_id' or $key eq 'Project_key') {
            $detail_cfg{$key} = $value;
	}
        if ($key eq 'Known_unigene' or $key eq 'Known_pep') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg{$key} = $value;
        }
        if ($key eq 'Known_anno') {
            $detail_cfg{$key} = $value;
        }
        if ($key eq 'Ref_seq' or $key eq 'Ref_ann') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg{$key} = $value;
        }
        if ($key =~/^SNP_/) {
            $detail_cfg{$key} = $value;
        }
        if ($key eq 'nr' or $key eq 'Swissprot' or $key eq 'Kegg' or $key eq 'Pfam' or $key eq 'Cog' or $key eq 'Kog') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg{anno}{$key} = $value;
        }

        if ($key eq 'Queue_type1') {
            $detail_cfg{$key} = $value;
        }

        if ($key eq 'Queue_type2') {
            $detail_cfg{$key} = $value;
        }
        print "$key: $value\n" if (exists $detail_cfg{$key});
    }
    close CFG;
    return %detail_cfg;
}
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
