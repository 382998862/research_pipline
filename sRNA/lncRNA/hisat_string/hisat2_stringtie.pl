#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../lncRNA_pip.cfg")};
my $BEGIN_TIME    = time();
my @Original_ARGV = @ARGV;
my $version="2.4.0";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ( $cfg1, $cfg2, $od, $step, $log, $oneStepOnly,$test );
GetOptions(
    "help|?"        => \&USAGE,
    "cfg1:s"        => \$cfg1,
    "cfg2:s"        => \$cfg2,
    "od:s"          => \$od,
    "s:s"           => \$step,
    "oneStepOnly:s" => \$oneStepOnly,
    "test"=>\$test,
) or &USAGE;
&USAGE unless ( $cfg1 and $cfg2 and $od );
################################

my $notename = `hostname`;
chomp $notename;

$cfg1 = &ABSOLUTE_DIR($cfg1);
$cfg2 = &ABSOLUTE_DIR($cfg2);
&MKDIR($od);
$od = &ABSOLUTE_DIR($od);
&MKDIR("$od/work_sh");

$step = $step || 1;

#
# log file
#
my $startTime = GetTime();
my $user      = `whoami`;
chomp $user;
my $workDir = `pwd`;
chomp $workDir;
my $task = "perl $Bin/$Script " . join( "\t", @Original_ARGV );

open( $log, ">", "$od/Hisat2_Stringtie." . time() . ".log" ) or die $!;

print $log "######################################\n";
print $log "$startTime\n";
print $log "user:	$user\n";
print $log "work directory:	$workDir\n";
print $log "job:	$task\n";
print $log "######################################\n\n";

print $log "data config file:  $cfg1\n";
print $log "detail config file:  $cfg2\n";
print $log "output dir:  $od\n";
print $log "start from step: $step\n\n";

#==================================================================
# bins
#==================================================================


#my $GFFREAD_BIN ="/share/nas2/genome/biosoft/cufflinks/2.2.1/gffread";   
my $GFFREAD_BIN=$config{GFFREAD_BIN};


#my $TOPHAT_BIN = "/share/nas2/genome/biosoft/tophat/2.0.13_fix/tophat2";
my $HISAT=$config{HISAT};
my $HISAT_BUILD=$config{HISAT_BUILD};

#my $GTF2FA_BIN = "/share/nas2/genome/biosoft/tophat/2.0.7/gtf_to_fasta";
my $STRINGTIE="$Bin/bin/stringtie";
#my $CUFFCOMPARE_BIN = "/share/nas2/genome/biosoft/cufflinks/2.2.1/cuffcompare"; 
my $CUFFCOMPARE_BIN=$config{CUFFCOMPARE_BIN};
my $GFFCOMPARE_BIN="$Bin/bin/gffcompare";
my $samtools=$config{samtools};
my $RSeQC="/share/nas2/genome/bmksoft/pipeline/LncRNA_pipeline/v3.1.6/bin/basic_analysis/hisat_string/RSeQC/RSeQC_analysis.pl";
#==================================================================
# load config file
#==================================================================

my %total_read;
my %para;
my %sample;

open( IN, "cat $cfg1 $cfg2|" ) || die "$!\n";
while (<IN>) {
    chomp;
    s/\r$//;
    s/^\s+//;
    s/\s+$//;
    next if ( /^\#/ || /^$/ );

    my @tmp = split /\s+/, $_;
    if ( $tmp[0] eq "Sample" ) {
        my $fq1 = <IN>;
        chomp $fq1;
        my $fq2 = <IN>;
        chomp $fq2;
        my @fq_1 = split /\s+/, $fq1;
        $sample{ $tmp[1] }{FQ1} = $fq_1[1];
        my @fq_2 = split /\s+/, $fq2;
        $sample{ $tmp[1] }{FQ2} = $fq_2[1];
    }
    $para{ $tmp[0] } = $tmp[1];
}
close IN;
if ( !-e "$od/totalRead.stat.xls" ) {
    system "perl $Bin/bin/get_totalread.pl -cfg1 $cfg1 -od $od -queue $para{Queue_type}";
}
my $strand;
if ($para{'Lib_type1'} eq "fr-firststrand" ){
	$strand="RF";
}elsif ($para{'Lib_type1'} eq "fr-secondstrand"){
    $strand="FR";
}
#==================================================================
# pipeline
#==================================================================

#######################################
#
# step 1: check and build bowtie index
#
######

my $genome;
my $genome_size;
my $gtf;
my $gff;
my @chromosome;
my $idx_prefix;
my $index=$para{'Project_key'};
if ( $step != 1 ) {
    &MKDIR("$od/Ref_Genome");
    chdir "$od/Ref_Genome";
    $genome     = basename( $para{Ref_seq} );
    my $genome_dir= dirname( $para{Ref_seq} );
    $idx_prefix = basename( $para{Ref_seq} );
    $idx_prefix =~ s/.fa$//;
    my $genome_path = dirname( $para{Ref_seq} );

    if ( !-f "$od/Ref_Genome/$genome".".hdrs" ) {
        `grep '>' $para{Ref_seq} > $od/Ref_Genome/$genome.hdrs`;    #2015/08/19,modified by niulg
    }
    if ( !-f "$idx_prefix.8.ht2" and !-f "$idx_prefix.8.ht2l"  ) {
        my $genome_abs_path = $para{Ref_seq};


        if (-e "$genome_dir/$idx_prefix.1.ht2") {
            system "ln -s $genome_dir/$genome $genome_dir/$idx_prefix.*.ht2 ./";
        } elsif (-e "$genome_dir/$idx_prefix.1.ht2l") {   # for large genomes more than about 4 billion nucleotides in length
            system "ln -s $genome_dir/$genome $genome_dir/$idx_prefix.*.ht2l ./";
        } else {
            system "ln -s $genome_dir/$genome ./";
            system "$HISAT_BUILD -p 4 $genome $idx_prefix";
        }
    }
    if (!-f "$od/Ref_Genome/$genome.fai"){
        if (-e "$genome_dir/$genome.fai"){
            system "ln -s $genome_dir/$genome.fai ./";
        }else{
            system "$config{samtools} faidx $idx_prefix.fa ";
        }

    }
    $genome_size = "$od/Ref_Genome/genome_size.txt";
    if (!-f $genome_size ) {
        open OUT, ">$genome_size";
        my %genomes =&fasta($para{Ref_seq});
	#print Dumper(\%genomes);
        foreach my $chro ( keys %genomes ) {
            #next if length($chro) > 11;    ##2015-11-12 for PSI genome by linhj
   		print "$chro\n"; 
            my $genome_len = length( $genomes{$chro} );
            if ( $genome_len > 1000000 ) {
                print OUT "$chro\t$genome_len\n";
            }
        }
        close OUT;
    }
    
    
    my $gff_name = basename( $para{Ref_ann} );
    $gff = $gff_name;
	if ($gff=~m/\.gff3?$/) {
		$gff_name =~ s/\.gff3?$//i;
		$gtf = "$gff_name.gtf";
	}elsif($gff=~m/\.gtf$/){
		$gtf=$gff_name;
	}
    chdir "../";

    #$step++ unless ($oneStepOnly) ;
    print STDOUT " check and build bowtie index  Done\n";
    print $log "check and build bowtie index  Done\n";
}
if ( $step == 1 ) {
    print STDOUT "=== check and build bowtie index ===\n";
    print $log "=== check and build bowtie index ===\n";
    &MKDIR("$od/Ref_Genome");
    chdir "$od/Ref_Genome";
    $genome     = basename( $para{Ref_seq} );
    my $genome_dir= dirname( $para{Ref_seq} );
    $idx_prefix = basename( $para{Ref_seq} );
    $idx_prefix =~ s/.fa$//;
   my $genome_path = dirname( $para{Ref_seq} );

    if ( !-f "$od/Ref_Genome/$genome".".hdrs" ) {
        `grep '>' $para{Ref_seq} > $od/Ref_Genome/$genome.hdrs`;    #2015/08/19,modified by niulg
    }
    if ( !-f "$idx_prefix.8.ht2" and !-f "$idx_prefix.8.ht2l"  ) {
        my $genome_abs_path = $para{Ref_seq};

        if (-e "$genome_dir/$idx_prefix.1.ht2") {
            system "ln -s $genome_dir/$genome $genome_dir/$idx_prefix.*.ht2 ./";
        } elsif (-e "$genome_dir/$idx_prefix.1.ht2l") {   # for large genomes more than about 4 billion nucleotides in length
            system "ln -s $genome_dir/$genome $genome_dir/$idx_prefix.*.ht2l ./";
        } else {
            system "ln -s $genome_dir/$genome ./";
            system "$HISAT_BUILD -p 4 $genome $idx_prefix";
        }
    }
    if (!-f "$od/Ref_Genome/$genome.fai"){
        if (-e "$genome_dir/$genome.fai"){
            system "ln -s $genome_dir/$genome.fai ./";
        }else{
            system "$config{samtools} faidx $idx_prefix.fa ";
        }

    }

   $genome_size = "$od/Ref_Genome/genome_size.txt";
	if (!-f "$genome_size" ) {
        open OUT, ">$genome_size";
        my %genomes =&fasta($genome);
        foreach my $chro ( keys %genomes ) {
	 my $genome_len = length( $genomes{$chro} );
            if ( $genome_len > 1000000 ) {
                print OUT "$chro\t$genome_len\n";
            }
        }
        close OUT;
	}
 
    ################################## gff2gtf
    my $gff_name = basename( $para{Ref_ann} );
    $gff = $para{Ref_ann};    #$gff_name;
    system "ln -s $para{Ref_ann} ./" unless ( -f $gff );

	if ($gff=~m/\.gff3?$/) {
	    	$gff_name =~ s/\.gff3?$//i;
    		$gtf = "$gff_name.gtf";
		system "$GFFREAD_BIN $gff -T -o $gff_name.gtf";
	}elsif($gff=~m/\.gtf$/){
    		$gtf=$gff_name;
		system "ln -s $gff ./  ";
	}
    chdir "../";

    $step++ unless ($oneStepOnly);
    print STDOUT " check and build bowtie index  Done\n";
    print $log "check and build bowtie index  Done\n";

}

####################################
#
# step 2: Align the RNA-seq read to genome using Tophat2 + bowtie2
#
#########

if ( $step == 2 ) {
    print STDOUT
      "=== Align the RNA-seq read to genome using Hisat2 ===\n";
    print STDOUT "Hisat mapping shell file: $od/work_sh/Hisat.sh\n";
    print $log
      "=== Align the RNA-seq read to genome using Hisat2   ===\n";
    print $log "Hisat mapping shell file: $od/work_sh/Hisat.sh\n";

    #
    # write shell
    #
#    open READ, "<$od/totalRead.stat.xls"|| die $!;
    my %reads=();
#    while(<READ>){
#	chomp;my @tmp=split /\s+/,$_;
#	next if(/^$/);
#	print "1:$_\n";
#	print "2:$tmp[0]\t$tmp[1]\n";
#	$reads{$tmp[0]}=$tmp[1];
#	print "3:$tmp[0]\t$reads{$tmp[0]}\n";
#   }
#    close READ;
    my @tmp = glob("$od/Mapped/*.reads.txt");
    for my $k (@tmp){
	open (IN,"<$k");
	while(<IN>){
		chomp;
		my @t = split /\s+/,$_;
		$reads{$t[0]}=$t[1];
	}
	close(IN);
    }

    ###########
    #
    ###########
    open OUT1, ">$od/work_sh/Hisat.sh" || die $!;
    open OUT2, ">$od/work_sh/Samtools.sh" || die $!;
    open MAP,">$od/work_sh/Total_bam_stat.sh" ||die $!;
    &MKDIR("$od/Hisat");
    &MKDIR("$od/Mapped");
    foreach my $sam ( sort keys %sample ) {
        &MKDIR("$od/Hisat/$sam");
        print OUT1 "cd $od/Hisat/$sam &&";
        print OUT1 " $HISAT --dta -p 6";
        print OUT1 " --rna-strandness $strand" if defined $strand;
        print OUT1 " --phred64" if (exists $para{Qphred} and $para{Qphred} eq "64");
        print OUT1 " -x $od/Ref_Genome/$idx_prefix -1 $sample{$sam}{FQ1} -2 $sample{$sam}{FQ2} -S $od/Hisat/$sam/$sam.HISAT_aln.sam > $od/Hisat/$sam/$sam.log\n";
        print OUT2  "$samtools view -F 4 -Su $od/Hisat/$sam/$sam.HISAT_aln.sam | $samtools sort - $od/Hisat/$sam/$sam.HISAT_aln.sorted && cd $od/Hisat/$sam && $samtools index $sam.HISAT_aln.sorted.bam $sam.HISAT_aln.sorted.bam.bai\n";
	print MAP "perl $Bin/bin/bam2map_stat.pl -i $sam -bam $od/Hisat/$sam/$sam.HISAT_aln.sorted.bam -totalRead $reads{$sam} -od $od/Mapped\n";
	print "perl $Bin/bin/bam2map_stat.pl -i $sam -bam $od/Hisat/$sam/$sam.HISAT_aln.sorted.bam -totalRead $reads{$sam} -od $od/Mapped\n";
    }
    close OUT1;
    close OUT2;
    close MAP;
    $para{Memory}     ||= "15G";
    $para{Queue_type} ||= "general.q";
    $para{CPU} ||= "30";
    &Cut_shell_qsub( "$od/work_sh/Hisat.sh",$para{CPU} , $para{Memory},"$para{Queue_type}" );
    &Check_qsub_error("$od/work_sh/Hisat.sh");
    &Cut_shell_qsub( "$od/work_sh/Samtools.sh",$para{CPU} , $para{Memory},"$para{Queue_type}" );
    &Check_qsub_error("$od/work_sh/Samtools.sh");

    ##################mapping efficiency##########################
    &Cut_shell_qsub( "$od/work_sh/Total_bam_stat.sh",$para{CPU} , $para{Memory},"$para{Queue_type}" );
    &Check_qsub_error("$od/work_sh/Total_bam_stat.sh");
#All.mappedStat.xls
    my $cmd="perl $Bin/map_stat/mapped_stat.pl -id $od/Mapped -o $od/Mapped/All.mappedStat.xls";
    print $cmd,"\n";system $cmd;
    $step++ unless ($oneStepOnly);
    print STDOUT "\n";
    print $log "\n";
}

#####################################
#
# step 3: Scripture Assembly Analysis
#
########
`rm $od/Hisat/*/*sam`;
if ( $step == 3 ) {
    print STDOUT "=== StrintTie analysis   ===\n";
    print STDOUT "shell file: $od/work_sh/StringTie.sh\n";
    print STDOUT "shell file: $od/work_sh/StringTie.sh\n";
    print $log "=== StrintTie analysis  ===\n";
    print $log "shell file: $od/work_sh/StringTie.sh\n";
    print $log "shell file: $od/work_sh/StringTie.sh\n";
    open OUT,">$od/work_sh/StringTie.sh" || die;
    open OUT1,">$od/work_sh/Cuffcompare.sh" || die;
    &MKDIR("$od/Cuffcompare");
    &MKDIR("$od/StringTie");
    $para{Memory}     ||= "15G";
    $para{Queue_type} ||= "general.q";
    $para{CPU} ||= "20";
    foreach my $sam (sort keys %sample){
        &MKDIR("$od/StringTie/$sam");
        &MKDIR("$od/Cuffcompare/$sam");
        print OUT1 "cd $od/Cuffcompare/$sam &&";
	if($strand eq "RF"){
		print OUT " $STRINGTIE $od/Hisat/$sam/$sam.HISAT_aln.sorted.bam -G $para{'Ref_ann'} -p 2 --rf -l $sam -o $od/StringTie/$sam/StringTie_asm.gtf\n";
	}else{
		print OUT " $STRINGTIE $od/Hisat/$sam/$sam.HISAT_aln.sorted.bam -G $para{'Ref_ann'} -p 2 --fr -l $sam -o $od/StringTie/$sam/StringTie_asm.gtf\n";
	}
        print OUT1 "$CUFFCOMPARE_BIN -r $od/Ref_Genome/$gtf $od/StringTie/$sam/StringTie_asm.gtf\n";
    }
    close OUT;
    &Cut_shell_qsub( "$od/work_sh/StringTie.sh",$para{CPU}, $para{Memory},$para{Queue_type} );
    &Check_qsub_error("$od/work_sh/StringTie.sh");
	###################Cuffcompare just for Alitsplice_Analysis 
    &Cut_shell_qsub( "$od/work_sh/Cuffcompare.sh",$para{CPU}, $para{Memory},$para{Queue_type} );
    &Check_qsub_error("$od/work_sh/Cuffcompare.sh");
    
    $step++ unless ($oneStepOnly);
    print STDOUT "StrintTie analysis Finished \n";
    print $log "StrintTie analysis Finished\n";
}
if ( $step == 4 ) {
    print STDOUT "=== StrintTie Merged analysis   ===\n";
    print STDOUT "shell file: $od/work_sh/StringTie_Merged.sh\n";
    print STDOUT "shell file: $od/work_sh/StringTie_Merged.sh\n";
    print $log "=== StrintTie Merged analysis  ===\n";
    print $log "shell file: $od/work_sh/StringTie_Merged.sh\n";
    print $log "shell file: $od/work_sh/StringTie_Merged.sh\n";
    open OUT,">$od/work_sh/StringTie_Merged.sh" || die;
    &MKDIR("$od/StringTie_Merged");
    $para{Memory}     ||= "15G";
    $para{Queue_type} ||= "general.q";
    $para{CPU} ||= "20";
   
    if($strand eq "RF"){ 
	print OUT "$STRINGTIE  --merge --rf -G $para{'Ref_ann'} -F 0.1 -T 0.1 -i -o $od/StringTie_Merged/StringTie_merged.gtf ";
    }else{
	print OUT "$STRINGTIE  --merge --fr -G $para{'Ref_ann'} -F 0.1 -T 0.1 -i -o $od/StringTie_Merged/StringTie_merged.gtf ";
    }
    foreach my $sam (sort keys %sample){
        print OUT " $od/StringTie/$sam/StringTie_asm.gtf ";
    }
    print OUT " \n";
    close OUT;
    &Cut_shell_qsub( "$od/work_sh/StringTie_Merged.sh","$para{CPU}", $para{Memory}, "$para{Queue_type}" );
    &Check_qsub_error("$od/work_sh/StringTie_Merged.sh");
    #########
    &MKDIR("$od/Compare");
    chdir "$od/Compare";
    my $cmd="$GFFCOMPARE_BIN -r $para{Ref_ann} $od/StringTie_Merged/StringTie_merged.gtf";
    system "$cmd";

    $step++ unless ($oneStepOnly);
    print STDOUT "StrintTie Merged Analysis Finished \n";
    print $log "StrintTie Merged Analysis Finished\n";
}


########################################对所有基因，已知lncRNA和剩余lncRNA进行定量####
if ( $step == 5 ) {

    `mkdir -p $od/LncPredict`    unless(-d "$od/LncPredict");
    `mkdir -p $od/genePredict`    unless(-d "$od/genePredict");
    &MKDIR("$od/Ballgown");
    #########################################新基因预测###################################
    my $cmd1 = "perl $Bin/gene_expression/newgene_expression.pl -gtf $od/Compare/gffcmp.annotated.gtf -genome $para{Ref_seq} -index $para{Project_key} -od $od/genePredict --queue $para{Queue_type} && $GFFREAD_BIN $od/genePredict/final_track/$para{Project_key}.newGene_final.filtered.gff -T -o $od/genePredict/$para{Project_key}.newGene_final.filtered.gtf && cat $od/Ref_Genome/$gtf $od/genePredict/$para{Project_key}.newGene_final.filtered.gtf >$od/genePredict/All_gene_final.gtf && cat $para{Ref_ann} $od/genePredict/final_track/$para{Project_key}.newGene_final.filtered.gff >$od/genePredict/All_gene_final.gff && perl $Bin/gene_expression/bin/filter.known_transcript.pl -gff $para{Ref_ann} -kfa $para{Known_unigene} -nfa $od/genePredict/final_track/$para{Project_key}.newGene.longest_transcript.fa -out $od/genePredict/All.longest_transcript.fa ";

    #########################################lncRNA过滤###################################
    my $cmd2	= "perl $Bin/lnc_predict/lnc_merger_filter.pl -gtf $od/Compare/gffcmp.annotated.gtf -out $od/LncPredict/merged_filter.gtf  ";
    $cmd2	.= "-n $para{'exon'} "      if (exists $para{'exon'});
    $cmd2	.= "-no_sense yes "         if (exists $para{'No_sense'});
    $cmd2	.= "&& perl $Bin/lnc_predict/lnc.stat.class_code.pl -i $od/LncPredict/merged_filter.gtf -o $od/LncPredict ";

    if(exists $para{'Lnc_ann'}){
	$cmd2	.="&& cp $para{'Lnc_ann'} $od/LncPredict/Known_lncRNA.gff ";
	$cmd2	.="&& perl $Bin/bin/lncRNA.gff_to_gtf.pl $od/LncPredict/Known_lncRNA.gff $od/LncPredict/Known_lncRNA.gtf $od/LncPredict/Known_lncRNA_id_name.list ";
	$cmd2	.="&& perl $Bin/lnc_predict/lncRNA_known_filter.pl $od/LncPredict/Known_lncRNA.gtf $od/LncPredict/merged_filter.gtf $od/LncPredict && cat $od/LncPredict/Known_lncRNA.gtf $od/LncPredict/unKnown_lncRNA.gtf >$od/LncPredict/All_LncRNA.gtf ";
    }else{
	$cmd2	.="&& cp $od/LncPredict/merged_filter.gtf $od/LncPredict/All_LncRNA.gtf && cp $od/LncPredict/merged_filter.gtf $od/LncPredict/unKnown_lncRNA.gtf";
    }
    $cmd2 .="&& cat $od/genePredict/All_gene_final.gtf $od/LncPredict/All_LncRNA.gtf >$od/Ballgown/All_RNA.gtf ";

    print "$cmd1\n"; system $cmd1;
    print "$cmd2\n"; system $cmd2;

    ########################################定量#########################################
    print STDOUT "===  Ballgown analysis   ===\n";
    print STDOUT "shell file: $od/work_sh/Ballgown.sh\n";
    print STDOUT "shell file: $od/work_sh/Ballgown.sh\n";
    print $log "=== Ballgown analysis  ===\n";
    print $log "shell file: $od/work_sh/Ballgown.sh\n";
    print $log "shell file: $od/work_sh/Ballgown.sh\n";
    open OUT,">$od/work_sh/Ballgown.sh" || die;
    &MKDIR("$od/Abundances_N_Coverage");
    &MKDIR("$od/Assembly");
    $para{Memory}     ||= "15G";
    $para{Queue_type} ||= "general.q";
    $para{CPU} ||= "20";
    foreach my $sam (sort keys %sample){
        &MKDIR("$od/Ballgown/$sam");
        &MKDIR("$od/Abundances_N_Coverage/$sam");
        &MKDIR("$od/Assembly/$sam");
	if($strand eq "RF"){
		print OUT "$STRINGTIE $od/Hisat/$sam/$sam.HISAT_aln.sorted.bam -eB -G $od/Ballgown/All_RNA.gtf -p 6 --rf -o $od/Ballgown/$sam/tringTie_asm.gtf -A $od/Abundances_N_Coverage/$sam/gene_abundence.tab -C $od/Abundances_N_Coverage/$sam/known.cov_refs.gtf \n";
	}else{
        	print OUT "$STRINGTIE $od/Hisat/$sam/$sam.HISAT_aln.sorted.bam -eB -G $od/Ballgown/All_RNA.gtf -p 6 --fr -o $od/Ballgown/$sam/tringTie_asm.gtf -A $od/Abundances_N_Coverage/$sam/gene_abundence.tab -C $od/Abundances_N_Coverage/$sam/known.cov_refs.gtf \n";
	}
    }
    close OUT;
    &Cut_shell_qsub( "$od/work_sh/Ballgown.sh",$para{CPU}, $para{Memory},$para{Queue_type} );
    &Check_qsub_error("$od/work_sh/Ballgown.sh");

    $step++ unless ($oneStepOnly);
    print STDOUT "Ballgown analysis Finished \n";
    print $log "Ballgown analysis  Finished\n";

    #################################get gene/transcript fpkm/count ##############################    
    $para{'Read_length'}||=150;
    my $prepDE_dir="$od/prepDE";
    &MKDIR("$prepDE_dir");
    open EXP,">$od/work_sh/quantity_exp.sh" || die;
    my $cmd = "python $Bin/bin/prepDE.py -i $od/Ballgown -gtf $od/Ballgown/All_RNA.gtf -g $prepDE_dir/Total_gene_count.csv -t $prepDE_dir/Total_transcript_count.csv -l $para{'Read_length'} && python $Bin/bin/prepDE_fpkm.py -i $od/Ballgown -gtf $od/Ballgown/All_RNA.gtf -g $prepDE_dir/Total_gene_fpkm.csv -t $prepDE_dir/Total_transcript_fpkm.csv -l $para{'Read_length'} && dos2unix $prepDE_dir/Total_gene_count.csv $prepDE_dir/Total_transcript_count.csv $prepDE_dir/Total_gene_fpkm.csv $prepDE_dir/Total_transcript_fpkm.csv && ";
    $cmd   .= "perl $Bin/bin/get_gene_rna_fpkm.pl -gene_count $prepDE_dir/Total_gene_count.csv -gene_fpkm $prepDE_dir/Total_gene_fpkm.csv -trans_count $prepDE_dir/Total_transcript_count.csv -trans_fpkm $prepDE_dir/Total_transcript_fpkm.csv -gene_gtf $od/genePredict/All_gene_final.gtf -novel_lnc $od/LncPredict/unKnown_lncRNA.gtf -od $prepDE_dir ";
    $cmd   .= " -lnc_gtf $od/LncPredict/Known_lncRNA.gtf "	if(exists $para{'Lnc_ann'});
    $cmd   .= " -fpkm $para{'fpkm'} "				if(exists $para{'fpkm'});
    $cmd   .= "&& perl $Bin/bin/filter_fpkm.pl -known_lncrna_count $prepDE_dir/Known_lncRNA_counts.list -known_lncrna_fpkm $prepDE_dir/Known_lncRNA_fpkm.list -known_lncrna_gff $od/LncPredict/Known_lncRNA.gff -known_lncrna_gtf $od/LncPredict/Known_lncRNA.gtf " if(exists $para{'Lnc_ann'});
    print EXP "$cmd\n";
    close(EXP);
    &Cut_shell_qsub( "$od/work_sh/quantity_exp.sh",$para{CPU}, $para{Memory},$para{Queue_type} );
    &Check_qsub_error("$od/work_sh/quantity_exp.sh");

}

###################################链特异性评估
if ( $step == 6 ) {
    print STDOUT "=== RSeQC Analysis  ===\n";
    print STDOUT "shell file: $od/work_sh/RSeQC_Analysis.sh\n";
    print STDOUT "shell file: $od/work_sh/RSeQC_Analysis.sh\n";
    print $log "=== RSeQC Analysis ===\n";
    print $log "shell file: $od/work_sh/RSeQC_Analysis.sh\n";
    print $log "shell file: $od/work_sh/RSeQC_Analysis.sh\n";
    open OUT,">$od/work_sh/RSeQC_Analysis.sh" || die;
    $para{Memory}     ||= "30G";
    $para{Queue_type} ||= "general.q";
    $para{CPU} ||= "30";
    print OUT "perl $RSeQC -i $od/Hisat -genome $para{'Ref_ann'} -o $od/Lib_type_stat -queue $para{Queue_type}  -cpu  $para{CPU}  -m  $para{Memory} ";
    close OUT;
    &Cut_shell_qsub( "$od/work_sh/RSeQC_Analysis.sh","$para{CPU}", $para{Memory}, "$para{Queue_type}" );
    &Check_qsub_error("$od/work_sh/RSeQC_Analysis.sh");

    ##$step++ unless ($oneStepOnly);
    print STDOUT "RSeQC Analysis Finished \n";
    print $log "RSeQC Analysis Finished\n";
}

#######################################################################################
my $endTime = GetTime();
print $log "$endTime\n";
print STDOUT "\nDone. Total elapsed time : ", time() - $BEGIN_TIME, "s\n";
print $log "\nDone. Total elapsed time : ", time() - $BEGIN_TIME, "s\n";
#######################################################################################
close($log);

#==================================================================
# subs
#==================================================================
sub fasta {
    my $fa = shift;
    my %Fasta;
    open IN, $fa || die;
    $/ = '>';
    <IN>;
    while (<IN>) {
        chomp;
        my ( $id, $seq ) = split /\n+/, $_, 2;
        my $seq_id = ( split /\s+/, $id )[0];

        #$seq=~s/\s+//g;
        $Fasta{$seq_id} = $seq;
    }
    close IN;
    return %Fasta;

}

sub LOAD_PARA {
    my $para_file = shift;
    my $para      = shift;

    my $error_status = 0;
    open IN, $para_file || die "fail open: $para_file";
    while (<IN>) {
        chomp;
        s/^\s+//;
        s/\s+$//;
        s/\r$//;
        next if ( /^$/ or /^\#/ );
        my ( $para_key, $para_value ) = split( /\s+/, $_ );
        $para->{$para_key} = $para_value;
        if ( !defined $para_value ) {
            warn "Non-exist: $para_key\n";
            $error_status = 1;
        }
    }
    close IN;
    die "\nExit due to error Undefined parameter\n" if ($error_status);
}

sub Cut_shell_qsub {    #Cut shell for qsub 1000 line one file
                        # &Cut_shell_qsub($shell,$cpu,$vf,$queue);
    my $shell = shift;
    my $cpu   = shift;
    my $vf    = shift;
    my $queue = shift;

    my $line = `less -S $shell |wc -l `;
	chomp $line;
    if ( $line <= 1000 ) {
        if ( $notename =~ /cluster/ ) {    ####2015-09-25
                                           #if ($notename=~/login\-0\-4/) {
            system "sh $config{qsub_sh} $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
        }
        else {
            system "sh $config{qsub_sh} $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
        }
    }
    if ( $line > 1000 ) {
        my @div = glob "$shell.div*";
        foreach (@div) {
            if ( -e $_ ) {
                system "rm $_";
            }
        }
        @div = ();
        my $div_index = 1;
        my $line_num  = 1;
        open IN, "$shell" || die;
        while (<IN>) {
            chomp;
            open OUT, ">>$shell.div.$div_index.sh" || die;
            if ( $line_num < 1000 ) {
                print OUT "$_\n";
                $line_num++;
            }
            else {
                print OUT "$_\n";
                $div_index++;
                $line_num = 1;
                close OUT;
            }
        }
        if ( $line_num != 1 ) {
            close OUT;
        }
        @div = glob "$shell.div*";
        foreach my $div_file (@div) {
            if ( $notename =~ /cluster/ ) {    ####2015-09-25
                system "sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
            }
            else {
                system	"sh $config{qsub_sh} $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
            }
        }
    }
}

sub Check_qsub_error {                         #
        # Check The qsub process if error happend
    my $sh         = shift;
    my @Check_file = glob "$sh*.qsub/*.Check";
    my @sh_file    = glob "$sh*.qsub/*.sh";

    if ( $#sh_file != $#Check_file ) {
        print "Their Some Error Happend in $sh qsub, Please Check..\n";
        die;
    }
    else {
        print "$sh qsub is Done!\n";
    }
}

sub GetTMR {    #
     #Get Total Mapped Reads from file line which contain string "Mapped Reads\s"
    my $fStat = shift;
    open( IN, "<", $fStat ) or die $!;
    while (<IN>) {
        if (/^Mapped Reads\s(\d+)/) {
            close(IN);
            return $1;
        }
    }
    close(IN);
    die "Error Reads Stat file.\n";
}

sub MKDIR {    # &MKDIR($out_dir);
    my ($dir) = @_;

    #	rmdir($dir) if(-d $dir);
    mkdir($dir) if ( !-d $dir );
}

sub ABSOLUTE_DIR {    #$pavfile=&ABSOLUTE_DIR($pavfile);
    my $cur_dir = `pwd`;
    chomp($cur_dir);
    my ($in) = @_;
    my $return = "";

    if ( -f $in ) {
        my $dir  = dirname($in);
        my $file = basename($in);
        chdir $dir;
        $dir = `pwd`;
        chomp $dir;
        $return = "$dir/$file";
    }
    elsif ( -d $in ) {
        chdir $in;
        $return = `pwd`;
        chomp $return;
    }
    else {
        warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
        exit;
    }

    chdir $cur_dir;
    return $return;
}

sub GetTime {
    my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) =
      localtime( time() );
    return sprintf(
        "%4d-%02d-%02d %02d:%02d:%02d",
        $year + 1900,
        $mon + 1, $day, $hour, $min, $sec
    );
}

sub USAGE {    #
    my $usage = <<"USAGE";
Program: Tophat&Cufflinks_Analysis Procedure
Version: $version
Contact: Meng Fei <mengf\@biomarker.com.cn>

Description:
	This program is a Procedure deal with Multiple Samples RNA_Analysis
	Tophat+Cufflinks Combination��Designed for RNA Analysis with a Reference Genome

	The program will calculate the Junctions, Transcripts(RABT) Assembly, FPKM of genes & isoforms

Usage:
    -cfg1             data config, rawdata & refseq path
    -cfg2             detail config, analysis parameters
	-od               output dir            must be given
	-s                step of the program   option,default 1;
                            1  Check the index of Genome & build index for alignment
                            2  run Hisat analysis
                            3   run StringTie analysis
                            4   merged gtf
                            5   abundence analysis
                            6   RSeQC Analysis
	-oneStepOnly      perform one of the pipeline steps above mentioned
USAGE
    print $usage;
    exit;
}
