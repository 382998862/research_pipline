#!/usr/bin/perl -w
# 
#Copyright (c) BMK 2013
#Writer           Mengf <mengf@biomarker.com.cn>
#Program Date   2013 
#Modifier         Mengf <mengf@biomarker.com.cn>
#Last modified  2012 
my $version="1.0";
my $BEGIN=time();

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $programe_dir=basename($0);
my $path=dirname($0);
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET=1;
use newPerlBase;

#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info
# ------------------------------------------------------------------
my ($Anno_dir,$gene,$odir,$out_index,$HELP);

my @anno;

GetOptions(
		"i:s"=>\$out_index,
		"gene:s"=>\$gene,
		"anno:s"=>\$Anno_dir,
		"od:s"=>\$odir,
		"help"=>\$HELP
	) or &USAGE;

&USAGE if (!defined $Anno_dir || !defined $odir || !defined $out_index) ;

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\n$programe_dir Start Time :[$Time_Start]\n\n";


my %file_db=%{readconf("$Bin/../../../Config/db_file.cfg")} ;

###------------------ע��·��----------------------------###



my ($cog_best,$kog_best,$nr_best,$nt_best,$trembl_best,$swissprot_best,$kegg_best,$go);
my ($cog,$kog,$eggNOG,$nr,$nt,$trembl,$swissprot,$kegg,$kobas,$eggNOG_class,$cog_class);
my $anno_file=$file_db{kegg_anno_file};

my @annos = split /,/,$Anno_dir;
foreach my $anno_dir (@annos) {
	$anno_dir = &ABSOLUTE_DIR($anno_dir);
	my @tabs = glob "$anno_dir/02.gene-annotation/*.blast.tab" ;
	#push @best_tabs,"$anno_dir/02.gene-annotation/kobas.annotation" ;
	foreach my $tab (@tabs) {
		my $file = "$tab.best";
		if ($tab =~ /\.eggNOG\./i) {
			$eggNOG .= "$tab ";
			$tab=~/^(.*).blast.tab/;
			$eggNOG_class .=" $1'.class' ";
			
		}
		if ($tab =~ /\.cog\./i) {
			$cog_best .= "$file ";
			$cog .= "$tab ";
			$tab=~/^(.*).blast.tab/;
			$cog_class .=" $1.class ";
		}
		if ($tab =~ /\.kog\./i) {
			$kog_best .= "$file ";
			$kog .= "$tab ";
		}
		if ($tab =~ /\.nr\./i) {
			$nr_best .= "$file ";
			$nr .= "$tab ";
		}
		if ($tab =~ /\.nt\./i) {
			$nt_best .= "$file ";
			$nt .= "$tab ";
		}
		if ($tab =~ /\.trembl\./i) {
			$trembl_best .= "$file ";
			$trembl .= "$tab ";
		}
		if ($tab =~ /\.swissprot\./i) {
			$swissprot_best .= "$file ";
			$swissprot .= "$tab ";
		}
		if ($tab =~ /\.kegg\./i) {
			$kegg_best .= "$file ";
			$kegg .= "$tab ";
		}
	#	if ($file =~ /kobas\.annotation/i) {
	#		$kobas .= "$file ";
	#	}
	}
	my @go_annot = glob "$anno_dir/*/*.annot";
	if ($#go_annot >= 0) {
		$go .= "$go_annot[0] ";
	}
}

&MKDIR($odir);
$odir=&ABSOLUTE_DIR($odir);

&MKDIR("$odir/Result");
&MKDIR("$odir/02.gene-annotation");
my $sh_dir="$odir/work_sh";
&MKDIR ($sh_dir);
###----------------ע�����ϳ���--------------------------###
my @genes = split/,/,$gene;
my $cat_gene;
foreach (@genes) {
	$_ = &ABSOLUTE_DIR($_);
	$cat_gene .= "$_ ";
}

system "cat $cat_gene >$odir/02.gene-annotation/$out_index.fa";


if (defined $eggNOG) {
	open OUT,">$sh_dir/eggNOG_cat.sh" || die $!;
	print OUT "cat $eggNOG >$odir/02.gene-annotation/$out_index.eggNOG.blast.tab \n";
	print OUT "cat $eggNOG_class >$odir/02.gene-annotation/$out_index.eggNOG.class \n";
	close OUT;

	system "sh $sh_dir/eggNOG_cat.sh";
	system "sed -i \'2,\$s\/^Query.*\/\/g\' $odir/02.gene-annotation/$out_index.eggNOG.blast.tab";
	system "sed -i \'2,\$s\/^#.*\/\/g\' $odir/02.gene-annotation/$out_index.eggNOG.class";
	system  "sed -i \'\/^\$\/d\'   $odir/02.gene-annotation/$out_index.eggNOG.blast.tab ";
        system "sed -i \'\/^\$\/d\'   $odir/02.gene-annotation/$out_index.eggNOG.class ";

	open OUT,">$sh_dir/eggNOG_tab2class.sh" || die $!;

	print OUT "perl $Bin/bin/eggNOGClassDrawer.pl -i $odir/02.gene-annotation/$out_index.eggNOG.class -o $odir/Result/$out_index.eggNOG.cluster \n";
	close OUT;

	system "sh $sh_dir/eggNOG_tab2class.sh";

	system "cp $odir/02.gene-annotation/$out_index.eggNOG.class $odir/Result/$out_index.eggNOG_class.txt";
	#system "cp $odir/02.gene-annotation/$out_index.Cog.cluster.svg $odir/Result/$out_index.Cog.cluster.svg";
	#system "cp $odir/02.gene-annotation/$out_index.eggNOG.cluster.png $odir/Result/$out_index.eggNOG.cluster.png";
	#system "cp $odir/02.gene-annotation/$out_index.eggNOG.class.stat $odir/Result/$out_index.eggNOG_class.stat.xls";
}
if (defined $cog) {
	open OUT,">$sh_dir/Cog_cat.sh" || die $!;
	print OUT "cat $cog >$odir/02.gene-annotation/$out_index.Cog.blast.tab \n";
	print OUT "cat $cog_best >$odir/02.gene-annotation/$out_index.Cog.blast.tab.best \n";
	my @cog_file =split / /,$cog_class;
	my $cog_class_new;
	foreach my $element (@cog_file)
	{
		next if ($element=~/^\s*$/);
		open COG,$element or die "fail $element";
		my $head = <COG>;
		if ($head =~/Organism/){
			system "cut -f 1-5,7-33 $element > $odir/02.gene-annotation/$out_index.Cog.class_tmp ";
			$cog_class_new .=" $odir/02.gene-annotation/$out_index.Cog.class_tmp ";
		}else{
			$cog_class_new .=" $element ";
		};
		close COG;
	}

	#print OUT "cat $cog_class_new >$odir/02.gene-annotation/$out_index.Cog.class \n";
	close OUT;

	system "sh $sh_dir/Cog_cat.sh";
	#system"rm -rf $odir/02.gene-annotation/$out_index.Cog.class_tmp";

	open OUT,">$sh_dir/Cog_tab2class.sh" || die $!;
	print OUT "perl $Bin/bin/cog_parser.pl $odir/02.gene-annotation/$out_index.Cog.blast.tab.best \n";

	print OUT "perl $Bin/bin/CogFunClassDrawer.pl -i $odir/02.gene-annotation/$out_index.Cog.class -o $odir/02.gene-annotation/$out_index.Cog.cluster.png \n";
	close OUT;

	system "sh $sh_dir/Cog_tab2class.sh";

	system "cp $odir/02.gene-annotation/$out_index.Cog.class $odir/Result/$out_index.Cog_class.txt";
	#system "cp $odir/02.gene-annotation/$out_index.Cog.cluster.svg $odir/Result/$out_index.Cog.cluster.svg";
	system "cp $odir/02.gene-annotation/$out_index.Cog.cluster.png $odir/Result/$out_index.Cog.cluster.png";
	system "cp $odir/02.gene-annotation/$out_index.Cog.class.stat $odir/Result/$out_index.Cog_class.stat.xls";
}
if (defined $kog) {
	open OUT,">$sh_dir/Kog_cat.sh" || die $!;
	print OUT "cat $kog >$odir/02.gene-annotation/$out_index.Kog.blast.tab \n";
	print OUT "cat $kog_best >$odir/02.gene-annotation/$out_index.Kog.blast.tab.best \n";
	close OUT;

	system "sh $sh_dir/Kog_cat.sh";

	open OUT,">$sh_dir/Kog_tab2class.sh" || die $!;
	print OUT "perl $Bin/bin/kog_parser.pl $odir/02.gene-annotation/$out_index.Kog.blast.tab.best \n";
	print OUT "perl $Bin/bin/KogFunClassDrawer.pl  -i $odir/02.gene-annotation/$out_index.Kog.class -o $odir/02.gene-annotation/$out_index.Kog.cluster.png \n";
	close OUT;

	system "sh $sh_dir/Kog_tab2class.sh";

	system "cp $odir/02.gene-annotation/$out_index.Kog.class $odir/Result/$out_index.Kog_class.txt";
	#system "cp $odir/02.gene-annotation/$out_index.Kog.cluster.svg $odir/Result/$out_index.Kog.cluster.svg";
	system "cp $odir/02.gene-annotation/$out_index.Kog.cluster.png $odir/Result/$out_index.Kog.cluster.png";
	system "cp $odir/02.gene-annotation/$out_index.Kog.class.stat $odir/Result/$out_index.Kog_class.stat.xls";
}


if (defined $nr) {
	open OUT,">$sh_dir/nr_cat.sh" || die $!;
	print OUT "cat $nr >$odir/02.gene-annotation/$out_index.nr.blast.tab \n";
	print OUT "cat $nr_best >$odir/02.gene-annotation/$out_index.nr.blast.tab.best \n";
	close OUT;

	system "sh $sh_dir/nr_cat.sh";

	open NR,"$odir/02.gene-annotation/$out_index.nr.blast.tab.best" || die $!;
	open OUT,">$odir/Result/$out_index.nr.anno.txt" || die $!;
	print OUT "#NrGeneID\tDatabase_ID\tE_value\tIdentity\tScore\tAnnotation\n";
	while (<NR>) 
	{
		chomp;
		&ann;
	}
	close OUT;
	close NR;

	open OUT,">$sh_dir/nr_pis_stat.sh" || die $!;
	print OUT "perl $Bin/bin/nr_nt_piegraph.pl -i $odir/Result/$out_index.nr.anno.txt -o $odir/Result -k $out_index.nr.anno \n";
	print OUT "perl $Bin/bin/Just_Pie.pl -i $odir/Result/$out_index.nr.lib.stat -o $odir/Result/$out_index.nr.lib.svg -w 600 -css $Bin/bin/pie12.css -note \"NR Species distribution\" \n";
	close OUT;

	system "perl $Bin/bin/nr_nt_piegraph.pl -i $odir/Result/$out_index.nr.anno.txt -o $odir/Result -k $out_index.nr.anno ";
	&nr_pie_stat("$odir/Result/$out_index.nr.anno.txt", "$odir/Result/$out_index.nr.lib.stat");
	system "perl $Bin/bin/Just_Pie.pl -i $odir/Result/$out_index.nr.lib.stat -o $odir/Result/$out_index.nr.lib.svg -w 600 -css $Bin/bin/pie12.css -note \"NR Species distribution\" ";

}

if (defined $kegg) {
	open OUT,">$sh_dir/Kegg_cat.sh" || die $!;
	print OUT "cat $kegg >$odir/02.gene-annotation/$out_index.Kegg.blast.tab \n";
	print OUT "cat $kegg_best >$odir/02.gene-annotation/$out_index.Kegg.blast.tab.best \n";
	#print OUT "cat $kobas >$odir/02.gene-annotation/kobas.annotation \n";
	close OUT;

	system "sh $sh_dir/Kegg_cat.sh";

	open OUT,">$sh_dir/Kegg_tab2path.sh" || die $!;
	print OUT "perl $Bin/bin/kegg_tab2path_ko.pl -tab $odir/02.gene-annotation/$out_index.Kegg.blast.tab.best -od $odir/Result -key $out_index.Kegg     \n";
	close OUT;

	system "perl $Bin/bin/kegg_tab2path_ko.pl -tab $odir/02.gene-annotation/$out_index.Kegg.blast.tab.best -od $odir/Result -key $out_index.Kegg ";

}

if (defined $go) {
	open OUT,">$sh_dir/GO_cat.sh" || die $!;
	print OUT "cat $go >$odir/02.gene-annotation/$out_index.GO.annot \n";
	close OUT;

	system "sh $sh_dir/GO_cat.sh";

	$out_index =~ /([^\.]+)/;
	my $go_graph_mark = $1;
	open OUT,">$sh_dir/GO_extract_graph.sh" || die $!;
	print OUT "perl $Bin/bin/GO_annot.pl -annot $odir/02.gene-annotation/$out_index.GO.annot -i $out_index -od $odir/Result \n";
	print OUT "cd $odir/Result && perl $Bin/bin/gene_ontology_graph.pl -i $odir/Result/$out_index.GO.list.txt -mark $go_graph_mark -o $odir/Result -k $out_index -c $file_db{component_file} -f $file_db{function_file} -p $file_db{process_file} \n";
	close OUT;

	system "sh $sh_dir/GO_extract_graph.sh";
}

if (defined $swissprot) {
	open OUT,">$sh_dir/Swissprot_cat.sh" || die $!;
	print OUT "cat $swissprot >$odir/02.gene-annotation/$out_index.Swissprot.blast.tab \n";
	print OUT "cat $swissprot_best >$odir/02.gene-annotation/$out_index.Swissprot.blast.tab.best \n";
	close OUT;

	system "sh $sh_dir/Swissprot_cat.sh";

	open SWISS,"$odir/02.gene-annotation/$out_index.Swissprot.blast.tab.best" || die $!;
	open OUT,">$odir/Result/$out_index.Swissprot.anno.txt" || die $!;
	print OUT "#SwissprotGeneID\tDatabase_ID\tE_value\tIdentity\tScore\tAnnotation\n";
	while (<SWISS>) 
	{
		chomp;
		&ann;
	}
	close OUT;
	close SWISS;
}

if (defined $trembl) {
	open OUT,">$sh_dir/TrEMBL_cat.sh" || die $!;
	print OUT "cat $trembl >$odir/02.gene-annotation/$out_index.TrEMBL.blast.tab \n";
	print OUT "cat $trembl_best >$odir/02.gene-annotation/$out_index.TrEMBL.blast.tab.best \n";
	close OUT;

	system "sh $sh_dir/TrEMBL_cat.sh";

	open TrEMBL,"$odir/02.gene-annotation/$out_index.TrEMBL.blast.tab.best" || die $!;
	open OUT,">$odir/Result/$out_index.TrEMBL.anno.txt" || die $!;
	print OUT "#TrEMBLGeneID\tDatabase_ID\tE_value\tIdentity\tScore\tAnnotation\n";
	while (<TrEMBL>) 
	{
		chomp;
		&ann;
	}
	close OUT;
	close TrEMBL;
}

if (defined $nt) {
	open OUT,">$sh_dir/nt_cat.sh" || die $!;
	print OUT "cat $nt >$odir/02.gene-annotation/$out_index.nt.blast.tab \n";
	print OUT "cat $nt_best >$odir/02.gene-annotation/$out_index.nt.blast.tab.best \n";
	close OUT;

	system "sh $sh_dir/nt_cat.sh";

	open NT,"$odir/02.gene-annotation/$out_index.nt.blast.tab.best" || die $!;
	open OUT,">$odir/Result/$out_index.nt.anno.txt" || die $!;
	print OUT "#NtGeneID\tDatabase_ID\tE_value\tIdentity\tScore\tAnnotation\n";
	while (<NT>) 
	{
		chomp;
		&ann;
	}
	close OUT;
	close NT;
}

####### Anno Integrate
if((defined $eggNOG|| defined $cog ||defined $kog) && (defined $nr || defined $nt || defined $swissprot || defined $trembl || defined $kegg)){
	open OUT,">$sh_dir/Anno_integrated.sh" || die $!;
	print OUT "perl $Bin/bin/anno_integrate2.pl -gene $odir/02.gene-annotation/$out_index.fa -id $odir/Result -od $odir/Result -key $out_index \n";
	close OUT;

	system "sh $sh_dir/Anno_integrated.sh";
}

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\n$programe_dir End Time :[$Time_End]\n\n";
&Runtime($BEGIN);
&show_log2("step_4: New gene functional annotation finished.");
#====================================================================================================================
#  +------------------+
#  |   subprogram     |
#  +------------------+
sub show_log2()
{
    open LOG, ">>$odir/../../../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}
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

sub ann
{
	@anno = split /\t/, $_;
	print OUT "$anno[0]\t$anno[4]\t$anno[13]\t$anno[8]\t$anno[12]\t$anno[15]\n";
}

sub nr_pie_stat{#&nr_pie_stat($in,$out)
	my ($in,$out)=@_;
	my %H;
	my $total=0;
	open (IN,$in) or die $!;
	<IN>;
	while (<IN>) {
		chomp;
		my @A=split/\t/,$_;;
		my $name=$A[-1];
		$name=~s/\s*$//;
		next unless $name=~/\[([^\]]+)\]$/;
		my $str = &cut_str($1);
		$H{$str}=0 unless exists $H{$str};
		$H{$str}++ if exists $H{$str};
		$total++;
	}
	close (IN) ;

	open (OUT,">$out") or die $!;
    print OUT "#Species_Name\tHomologous_Number\tRatio\n";
	my $limit=keys %H;
	if ($limit<=10) {
		if ($limit == 1) {
			my ($str,$value) = each %H;
			my $key = sprintf("%.2f",100*$value/$total).'%';
			print OUT "$str\t$value\t$key\n";
		}else{
			foreach my $str (sort {$H{$b}<=>$H{$a}} keys %H) {
				my $key = sprintf("%.2f",100*$H{$str}/$total).'%';
				print OUT "$str\t$H{$str}\t$key\n";
			}
		}
	}
	else {
		my $n=0;
		my $other=0;
		foreach my $str (sort {$H{$b}<=>$H{$a}} keys %H) {
			$n++;
			if($n<10){
				my $key = sprintf("%.2f",100*$H{$str}/$total).'%';
				print OUT "$str\t$H{$str}\t$key\n";
			}
		$other+=$H{$str} unless $n<10;
		}
		my $key=sprintf("%.2f",100*$other/$total).'%';
	print OUT "Other\t$other\t$key\n";
	}
	close OUT;
}


sub cut_str
{
	my $string = shift;
	my @str = split /\s+/,$string;
	if (@str > 2) {
		return "$str[0] $str[1]"
	}else{
		return $string;
	}
}

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

sub USAGE 
{
	print <<"	Usage End.";
	Program: Indegrate diff GeneUnits Annotation into One Combination
	Version: $version
	Contact: Meng Fei <mengf\@biomarker.com.cn>

	Description:

                -i             Out File index                                         must be given
                -gene          Annotation integrate gene files (split by ",")         must be given
                -anno          Annotation output dir           (split by ",")         must be given
                -od            OUT files DIR                                          must be given

	Example:
	  perl  Anno_integrated_v.1.0.pl -i AAA_Unigene -gene fa1,fa2  -anno fa1_Annotation/,fa2_Annotation/  -od allgene_Annotation/

	Usage End.
	exit;
}
