#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Config::General;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;
my $Title= "alignment.pl";
my $version="1.0";

#-------------------------------------------------
# GetOptions
#-------------------------------------------------
my($id,$od,$main_cfg,$filter,$step,$ones,$genome,$gff,$data_cfg);
GetOptions(
	"help|?"		=>\&USAGE,
	"od:s"			=>\$od,
	"id:s"			=>\$id,
	"data_cfg:s"		=>\$data_cfg,
	"main_cfg:s"		=>\$main_cfg,
	"filter:s"		=>\$filter,
	"ones"			=>\$ones,
	"step:s"		=>\$step,
	"genome:s"		=>\$genome,
	"gff:s"			=>\$gff,
) or &USAGE;
&USAGE unless ($od and $main_cfg and $data_cfg);

$id||="$od/../../Data_Assess/miRNA";

if(!-d $id){
	print "$id\n";
	print "Parameter id must be given!\n";
	die;
}

######################### prepare
$od=abs_path($od);
mkdir $od unless (-d $od);
$data_cfg=abs_path($data_cfg);
$main_cfg=abs_path($main_cfg);
my $work_sh="$od/work_sh";
mkdir $work_sh unless (-d $work_sh);
my $Ref_Database="$od/Ref_Database";
mkdir $Ref_Database unless (-d $Ref_Database);

#createLog($Title,$version,$$,"$od/log/","test");

########################## database or database
my %CFG=%{readconf("$Bin/../CFG")};
my $tRNA	= $CFG{tRNA};
my $silver	= $CFG{silver};
my $rfam	= $CFG{rfam};
my $repbase	= $CFG{repbase};

my $intersectBed = $CFG{intersectBed};
my $bowtie_build = $CFG{bowtie_build};
my $samtools	 = $CFG{samtools};
my $bedtools	 = $CFG{bedtools};



########### default values
$filter ||=1;### whether filter rRNA,tRNA,rfam,repbase
$step 	||=0;

############detail config
my %Raw_FQ;

open(IN,$data_cfg) or die $!;
while(<IN>){
	chomp;
	next if(/^#/||/^$/);
	if(/^NAME/){
		my $name=(split/\s+/,$_)[1];
		my $fq=<IN>;
		$Raw_FQ{$name}=(split/\s+/,$fq)[1];
		unless (-e $Raw_FQ{$name}) {
			print "Check Your config file:$data_cfg!Expecially $name\n\n";die;
		}
	}
}
close IN;

my %config=&readConfig($main_cfg);


################################ Filter rRNA
if($step==0){
	&log_current_time("Filter rRNA begin:");
	#stepStart(0,"Filter rRNA");
	if($filter==1){
		unless (-e $tRNA){print "there is no $tRNA lib file\n";}
		unless (-e $silver){print "there is no $silver lib file \n";}
		unless (-e $rfam){print "there is no $rfam lib file \n";}
		unless (-e $repbase){print "there is no $repbase lib file \n"}

		open (SH,">$work_sh/step0_filter_rRNA_rep.sh") or die $!;
		foreach my $name(sort keys %Raw_FQ){
			mkdir "$od/$name" unless (-d "$od/$name");
			mkdir "$od/$name/tRNA_align" unless (-d "$od/$name/tRNA_align");
			if(!-e "$od/$name/$name.clean.fa"){
				`ln -s $id/$name.clean.fa $od/$name`;
			}
			if(!-e "$od/$name/$name.collapse.fa"){
				`ln -s $id/$name.collapse.fa $od/$name`;
			}
			print SH "perl $Bin/bin/data_align.pl -fa $od/$name/$name.clean.fa -mis 0 -db $tRNA -out $od/$name/tRNA_align/$name.tRNA.bowtie.out &&";
			print SH "perl $Bin/bin/tRNA_tackle.pl -align $od/$name/tRNA_align/$name.tRNA.bowtie.out -key $name -od $od/$name/tRNA_align \n";

			mkdir "$od/$name/silver_align" unless (-d "$od/$name/silver_align");
			print SH "perl $Bin/bin/data_align.pl -fa $od/$name/$name.clean.fa -mis 0 -db $silver -out $od/$name/silver_align/$name.silver.bowtie.out &&";
			print SH "perl $Bin/bin/silver_tackle.pl -align $od/$name/silver_align/$name.silver.bowtie.out -key $name -od $od/$name/silver_align \n";

			mkdir "$od/$name/rfam_align" unless (-d "$od/$name/rfam_align");
			print SH "perl $Bin/bin/data_align.pl -fa $od/$name/$name.clean.fa -mis 0 -db $rfam -out $od/$name/rfam_align/$name.rfam.bowtie.out &&";
			print SH "perl $Bin/bin/Rfam_tackle2.pl -align $od/$name/rfam_align/$name.rfam.bowtie.out -key $name -od $od/$name/rfam_align \n";

			mkdir "$od/$name/Repbase_align" unless (-d "$od/$name/Repbase_align");
			print SH "perl $Bin/bin/data_align.pl -fa $od/$name/$name.clean.fa -mis 0 -db $repbase -out $od/$name/Repbase_align/$name.Repbase.bowtie.out &&";
			print SH "perl $Bin/bin/repeat_tackle.pl -align $od/$name/Repbase_align/$name.Repbase.bowtie.out -key $name -od $od/$name/Repbase_align \n";
		}
		close SH;
		qsubOrDie("$od/work_sh/step0_filter_rRNA_rep.sh",$CFG{queue},$CFG{cpu},$CFG{vf});

		open (SH,">$work_sh/step1_no_rRNA_rep.sh") or die $!;
		foreach my $name(sort keys %Raw_FQ){
			print SH "perl $Bin/bin/select_fa_by_file.pl -fa $od/$name/$name.clean.fa -i $od/$name/tRNA_align/$name.bowtie_tRNA.list $od/$name/silver_align/$name.bowtie_silver.list $od/$name/rfam_align/$name.bowtie_Rfam.list $od/$name/Repbase_align/$name.bowtie_Repbase.list -o $od/$name/$name.no_rRNA_rep.fa \n";
		}
		close SH;
		qsubOrDie("$od/work_sh/step1_no_rRNA_rep.sh",$CFG{queue},$CFG{cpu},$CFG{vf});
	}
	&log_current_time("Filter rRNA done!");
	#stepTime(0);
	$step++ unless ($ones);
}


############################ rRNA_no collapse
if($step==1){
	&log_current_time("rRNA_no collapse begin:");
	#stepStart(1,"rRNA_no collapse");
	open(SH ,">$work_sh/step2_no_rep_collapse.sh") or die $!;
	foreach my $name (sort keys %Raw_FQ){
		print SH "perl $Bin/bin/collapse_reads_md.pl $od/$name/$name.no_rRNA_rep.fa $name >$od/$name/$name.no_rRNA_rep.collapsed.fa\n";
	}
	close SH;
	qsubOrDie("$od/work_sh/step2_no_rep_collapse.sh",$CFG{queue},$CFG{cpu},$CFG{vf});
	&log_current_time("rRNA_no collapse done!");
	#stepTime(1);
	$step++ unless ($ones);
}


############### Creat genome file
if($step==2){
	&log_current_time("Creat genome file begin:");
	#stepStart(2,"Creat genome file");
	chdir "$Ref_Database";
	my $genome_abs_path = $genome;
	my $name=basename($genome);
	open(SH,">$work_sh/step3_genome_build.sh") or die $!;
	print SH "cd $Ref_Database &&";
	print SH "perl $Bin/bin/select_fa.pl -i All -fa $genome -o $Ref_Database/$name >>$Ref_Database/select_genome.out";
	close SH;
	qsubOrDie("$od/work_sh/step3_genome_build.sh",$CFG{queue},$CFG{cpu},$CFG{vf});




	if(!-f "$Ref_Database/$name"){
		runOrDie("ln -s $genome_abs_path ./");
	}
	if ((!-f "$Ref_Database/$name.1.ebwt") and (!-f "$Ref_Database/$name.1.ebwtl") ) {
        if (-e "$genome_abs_path.1.ebwt") {
            #runOrDie("ln -s $genome_abs_path.*.ebwt $Ref_Database");
			system "ln -s  $genome_abs_path  $genome_abs_path.*.ebwt $Ref_Database";
        } elsif (-e "$genome_abs_path.1.ebwtl") {   # for large genomes more than about 4 billion nucleotides in length
            #runOrDie("ln -s $genome_abs_path.*.ebwtl  $Ref_Database");
			system "ln -s $genome_abs_path  $genome_abs_path.*.ebwtl  $Ref_Database";
        } else {
			#runOrDie("$bowtie_build --large-index $Ref_Database/$name  $name");
			system ("$bowtie_build $Ref_Database/$name  $name");
			#system "ln -s $genome_abs_path ./";
            #system "bowtie-build --large-index $Ref_Database/$name  $name";
        }
	}

	$genome="$Ref_Database/$name";
	&log_current_time("Creat genome file done!");
	#stepTime(2);
	$step++ unless ($ones);
}


############### genome mapping
if($step==3){
	&log_current_time("Genome mapping begin:");
	#stepStart(3,"Genome mapping");
	unless (-e $genome){ print "there is not $genome lib file\n";}
	my $genome_name=basename($genome);
	$genome="$Ref_Database/$genome_name";
	open(GEN,"$Ref_Database/$genome_name") or die $!;
	$/=">";<GEN>;
	open(GL,">$Ref_Database/genome_len.txt");
	while(<GEN>){
		chomp;
		my @l=(split /\n/,$_,2);
		my $n=$l[0];
		my $l=length($l[1]);
		print GL "$n\t$l\n";
	}
	close GL;
	close GEN;
	$/="\n";

	open(SH,">$work_sh/step4_genome_mapping.sh") or die $!;
	foreach my $name (sort keys %Raw_FQ){
		my $genome_map="$od/$name/genome_map";
		mkdir "$genome_map" unless (-d "$genome_map");
		print SH "perl $Bin/bin/data_align_v2.pl -fa $od/$name/$name.no_rRNA_rep.fa -mis 0 -db $genome -out $genome_map/$name.genome_map.bowtie.out -sam 1 &&";
		print SH "$samtools view -bS $genome_map/$name.genome_map.bowtie.out >$genome_map/$name.genome_map.bowtie.bam &&";
		print SH "$samtools sort $genome_map/$name.genome_map.bowtie.bam $genome_map/$name.genome_map.sort &&";
		print SH "$samtools depth $genome_map/$name.genome_map.sort.bam >$genome_map/$name.genome_map.sort.bam.depth \n";
	}
	close SH;
	qsubOrDie("$od/work_sh/step4_genome_mapping.sh",$CFG{queue},$CFG{cpu},$CFG{vf});
	&log_current_time("Genome mapping done!");
	#stepTime(3);
	$step ++ unless ($ones);
}


###################### genome mapping2
if($step==4){
	&log_current_time("genome mapping2 begin:");
	#stepStart(4,"genome mapping2");
	open(SH,">$work_sh/step5_genome_mapping2.sh") or die $!;
	foreach my $name (sort keys %Raw_FQ){
		print SH "$bedtools genomecov -ibam $od/$name/genome_map/$name.genome_map.sort.bam -g $Ref_Database/genome_len.txt -d -strand - |awk '\$3!=0' >$od/$name/genome_map/$name.bam.minus.depth && ";
		print SH "$bedtools genomecov -ibam $od/$name/genome_map/$name.genome_map.sort.bam -g $Ref_Database/genome_len.txt -d -strand + |awk '\$3!=0' >$od/$name/genome_map/$name.bam.plus.depth && ";
		if (exists $config{medical}){
			print SH "perl $Bin/bin/plotReadDensity2.pl -i $od/$name/genome_map/$name.bam.plus.depth:$od/$name/genome_map/$name.bam.minus.depth -o $od/$name/genome_map -k $name -queue $CFG{queue} -medical $config{medical} && ";
		}else {
			print SH "perl $Bin/bin/plotReadDensity2.pl -i $od/$name/genome_map/$name.bam.plus.depth:$od/$name/genome_map/$name.bam.minus.depth -o $od/$name/genome_map -k $name -queue $CFG{queue} && ";
		}
		print SH "perl $Bin/bin/select_fa_by_file2.pl -fa $od/$name/$name.no_rRNA_rep.fa -i $od/$name/genome_map/$name.genome_map.bowtie.out -o $od/$name/genome_map/$name.genome_map.fa \n";
	}
	close SH;
	qsubOrDie("$od/work_sh/step5_genome_mapping2.sh",$CFG{queue},$CFG{cpu},$CFG{vf});
	&log_current_time("genome mapping2 done!");
	#stepTime(4);
	$step++ unless($ones);
}


############################# base,length distribution
if($step==5){
	&log_current_time("Base Length distribution done!");
	#stepStart(5,"Base Length distribution");

	my $seq=join ",",keys %Raw_FQ;
	runOrDie("perl $Bin/bin/Total_length_stat_by_detail.pl -id $od -samples $seq");
	my %Unanno;		###Statistic reads for genome map
	my %Genome_map;	###Statistic reads mapped;
	if($filter==1){
		foreach my $name (sort keys %Raw_FQ){
			#### genbank is best
			my %NC;
			$NC{scRNA}=0;
			$NC{snoRNA}=0;
			$NC{snRNA}=0;
			$NC{tRNA}=0;
			$NC{rRNA}=0;
			$NC{Repbase}=0;

			open(IN,"$od/$name/Repbase_align/$name.bowtie_Repbase.list") or die $!;
			while(<IN>){
				chomp;
				my ($s,$t)=(split /\t/,$_)[0,1];
				$NC{tag}{$s}='Repbase';
			}
			close IN;

			open(IN,"$od/$name/rfam_align/$name.bowtie_Rfam.list") or die $!;
			while(<IN>){
				chomp;
				my ($s,$t)=(split /\t/,$_)[0,1];
				$NC{tag}{$s}=$t;
			}
			close IN;

			open(IN,"$od/$name/tRNA_align/$name.bowtie_tRNA.list") or die $!;
			while(<IN>){
				chomp;
				my ($s,$t)=(split /\t/,$_)[0,1];
				$NC{tag}{$s}=$t;
			}
			close IN;

			open(IN,"$od/$name/silver_align/$name.bowtie_silver.list") or die $!;
			while(<IN>){
				chomp;
				my ($s,$t)=(split /\t/,$_)[0,1];
				$NC{tag}{$s}=$t;
			}
			close IN;

			foreach my $key (keys %{$NC{tag}}){
				$NC{$NC{tag}{$key}}++;
			}

			open (IN,"$od/$name/Len_stat/$name.Total.stat") or die $!;
			my ($total_num,$genome_num,$Unannotated)=(0,0,0);
			my ($per_rRNA,$per_scRNA,$per_snRNA,$per_snoRNA,$per_tRNA,$per_Repbase,$per_map)=(0,0,0,0,0,0,0);
			while(<IN>){
				chomp;
				if(/^Total/){
					($total_num,$genome_num)=(split/\s+/,$_)[1,2];
				}
			}
			close IN;

			$Unannotated=$total_num-$NC{scRNA}-$NC{snoRNA}-$NC{snRNA}-$NC{tRNA}-$NC{rRNA}-$NC{Repbase};
			$Unanno{$name}=$Unannotated;
			$Genome_map{$name}=$genome_num;

			open(OUT,">$od/$name/Data.stat") or die $!;
			print OUT "Types\tNumber\tPercentage\n";
			print OUT "Total\t$total_num\t100.00%\n";

			$per_rRNA=sprintf "%.2f",$NC{rRNA}/$total_num*100;
			print OUT "rRNA\t$NC{rRNA}\t$per_rRNA%\n";

			$per_scRNA=sprintf "%.2f",$NC{scRNA}/$total_num*100;
			print OUT "scRNA\t$NC{scRNA}\t$per_scRNA%\n";

			$per_snRNA=sprintf "%.2f",$NC{snRNA}/$total_num*100;
			print OUT "snRNA\t$NC{snRNA}\t$per_snRNA%\n";

			$per_snoRNA=sprintf "%.2f",$NC{snoRNA}/$total_num*100;
			print OUT "snoRNA\t$NC{snoRNA}\t$per_snoRNA%\n";

			$per_tRNA=sprintf "%.2f",$NC{tRNA}/$total_num*100;
			print OUT "tRNA\t$NC{tRNA}\t$per_tRNA%\n";

			$per_Repbase=sprintf "%.2f",$NC{Repbase}/$total_num*100;
			print OUT "Repbase\t$NC{Repbase}\t$per_Repbase%\n";

			$per_map=(100-$per_rRNA-$per_scRNA-$per_snRNA-$per_snoRNA-$per_tRNA-$per_Repbase);
			print OUT "Unannotated\t$Unannotated\t$per_map%\n";

			close OUT;
		}
	}
	
	open(MAP,">$od/All_sample_map.stat") or die $!;
	print MAP "#BMK-ID\tTotal_Reads\tMapped_Reads\tMapped_reads(+)\tMapped_reads(-)\n";

	my %genome_plus;
	my %genome_minus;
	foreach my $name(sort keys %Raw_FQ){
		my %plus;
		my %minus;
		open(IN,"$od/$name/genome_map/$name.genome_map.bowtie.out") or die $!;
		while(<IN>){
			chomp;
			next if(/^\@/||/^\s+/);
			my @line=(split /\s+/,$_);
			if($line[1]==0){
				$plus{$line[0]}=1;
			}
			elsif($line[1]==16){
				$minus{$line[0]}=1;
			}
		}
		close IN;

		my @plus=keys %plus;
		my @minus=keys %minus;
		my $plus=$#plus+1;
		my $minus=$#minus+1;
		$genome_plus{$name}=$plus;
		$genome_minus{$name}=$minus;


		print MAP "$name\t$Unanno{$name}\t$Genome_map{$name}\t$genome_plus{$name}\t$genome_minus{$name}\n";
	}
	close MAP;
	
	&log_current_time("Base Length distribution done!");
	#stepTime(5);
	$step++ unless($ones);
}


############################ common_specific
if($step==6){
	&log_current_time("Common_specific begin:");
	#stepStart(6,"Common_specific");
	my @samples=(sort keys %Raw_FQ);
	if(@samples<=1){
		print "Only one sample ,so there are not common and specific sRNA\n";
	}
	else{
		mkdir "$od/Common_specific" unless (-d "$od/Common_specific");
		for (my $i=0;$i<@samples;$i++){
			for (my $j=$i+1;$j<@samples;$j++){
				runOrDie("perl $Bin/bin/common_specific_stat.pl -1 $od/$samples[$i]/$samples[$i].collapse.fa -2 $od/$samples[$j]/$samples[$j].collapse.fa -s1 $samples[$i] -s2 $samples[$j] -od $od/Common_specific/$samples[$i]\_vs\_$samples[$j]");
			}
		}
	}
	&log_current_time("Common_specific done!");
	#stepTime(6);
	$step++ unless($ones);
}



############## Creat gff
if(defined $gff){
	&log_current_time("Creat gff begin:");
	my $gff_name=basename($gff);
	#stepStart(7,"Creat gff");
	if($gff =~ /\.gff/){
		if(-f "$Ref_Database/exon.gff"){`rm $Ref_Database/exon.gff`;}
		if(-f "$Ref_Database/intron.gff"){`rm $Ref_Database/intron.gff`;}
		runOrDie("perl $Bin/bin/gff_prepare.pl -in $gff -prefix exon,intron -out $Ref_Database");
		runOrDie("perl $Bin/bin/gff_intron.pl -in $gff -out $Ref_Database/intron.gff");
	}
	elsif($gff=~/gtf/){
		runOrDie("perl $Bin/bin/gtf2gff3.pl --cfg $Bin/bin/gtf2gff3.cfg $gff >$Ref_Database/$gff_name.gff3");
		runOrDie("perl $Bin/bin/gff_prepare.pl -in $Ref_Database/$gff_name.gff3 -prefix exon,intron -out $Ref_Database");
		if(!-f "$Ref_Database/intron.gff"){
			runOrDie("perl $Bin/bin/gff_intron.pl -in $Ref_Database/$gff_name.gff3 -out $Ref_Database/intron.gff");
		}
	}
	else{
		print "Wrong gff file! Must be gtf or gff3\n"; die;
	}
	&log_current_time("Creat gff done!");
	#stepTime(7);

	##############Exon Intron stat
	&log_current_time("Exon_Intron stat begin:");
	#stepStart(8,"Exon_Intron stat ");

	runOrDie("perl $Bin/bin/gff2bed.pl -i $Ref_Database/exon.gff -o $Ref_Database/exon.bed");
	runOrDie("perl $Bin/bin/gff2bed.pl -i $Ref_Database/intron.gff -o $Ref_Database/intron.bed");

	open (SH,">$work_sh/step6_gff.sh") or die $!;
	foreach my $name(sort keys %Raw_FQ){
		my $gff_path="$od/$name/genome_map/gff";
		mkdir $gff_path unless (-d $gff_path);
		&cmd_call("perl $Bin/bin/bowtie2bed.pl -i $od/$name/genome_map/$name.genome_map.bowtie.out -o $gff_path/$name.genome_map.bed");
		print SH "$intersectBed -a $Ref_Database/exon.bed -b $gff_path/$name.genome_map.bed -wa -wb -s >$gff_path/$name.exon_plus.bed\n";
		print SH "$intersectBed -a $Ref_Database/exon.bed -b $gff_path/$name.genome_map.bed -wa -wb -S >$gff_path/$name.exon_minus.bed\n";
		print SH "$intersectBed -a $Ref_Database/intron.bed -b $gff_path/$name.genome_map.bed -wa -wb -s >$gff_path/$name.intron_plus.bed\n";
		print SH "$intersectBed -a $Ref_Database/intron.bed -b $gff_path/$name.genome_map.bed -wa -wb -S >$gff_path/$name.intron_minus.bed\n";
	}
	close SH;
	qsubOrDie("$od/work_sh/step6_gff.sh",$CFG{queue},$CFG{cpu},$CFG{vf});

	############# statistic
	open(SH,">$work_sh/step7_exon_intron_stat.sh") or die $!;
	foreach my $name(sort keys %Raw_FQ){
		print SH "perl $Bin/bin/exon_intron_stat.pl -ep $od/$name/genome_map/gff/$name.exon_plus.bed -em $od/$name/genome_map/gff/$name.exon_minus.bed -ip $od/$name/genome_map/gff/$name.intron_plus.bed -im $od/$name/genome_map/gff/$name.intron_minus.bed -o $od/$name/genome_map/gff/exon_intron.stat -name $name\n";
	}
	close SH;
	qsubOrDie("$work_sh/step7_exon_intron_stat.sh",$CFG{queue},$CFG{cpu},$CFG{vf});
	&log_current_time("Exon_Intron stat done!");
	#stepTime(8);
}


####################
sub cmd_call {
	print "@_\n";
	system(@_) == 0 or die "system @_ failed: $?";
}
###############################
sub readConfig{
        my $configFile=shift;
        my $d=Config::General->new(-ConfigFile => "$configFile");
        my %config=$d->getall;
        return %config;
}

#############################sub
sub log_current_time{
        # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}


sub date_time_format{
        my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}





sub USAGE {#
	my $usage=<<"USAGE";
Program: $Script
Version: $version
Usage:
  Options:
	-id		input directory contained *.clean.fa, *.collapse.fa
	-od		output directory
	-data_cfg	data.cfg
	-main_cfg	detail.cfg
	-genome		genome.fa
	-gff		genome.gff
	-filter		1/0 (default:1)
	-step	<num>	default 0
		0 Filter rRNA
		1 rRNA_no collapse
		2 Creat genome file
		3 genome mapping
		4 genome mapping2
		5 base,length distribution
		6 common_specific
	-ones		only one step
	-h		Help

USAGE
	print $usage;
	exit;
}
