#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2012
# Writer:         Zhangyh <zhangyh@biomarker.com.cn>
# Program Date:   2012/8/20 10:00
# Modifier:       Zhangyh <zhangyh@biomarker.com.cn>
# Last Modified:  2012/9/13 18:19
# Modifier:       Xugl <xugl@biomarker.com.cn>
# Last Modified:  2015/05/25 16:00


my $ver="1.1.0";

use GD;
use GD::Polyline;
use Cwd;
use strict;
use Getopt::Long;
use FindBin qw/$Bin $Script/;
use Data::Dumper;
use File::Basename qw(basename dirname);
use Text::NSP::Measures::2D::Fisher::right;
use newPerlBase;

my %opts;
GetOptions(\%opts,
	"i=s",
	"o=s",
	"k=s",
	"d=s",
	"func=s",
	"v",
	"h",
);
&help() if ( defined($opts{h}) || !defined($opts{i}) ||
			!defined($opts{o}) || !defined($opts{k}) || 
			!defined($opts{d})  );

my $BEGIN=time();
my $Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";

#��������������  program begin ������������������

my $dir=getcwd;
$opts{o} = "$dir/$opts{o}" unless $opts{o} =~ /^[\/~]/;
my $go_out = "$opts{o}/go_enrichment";
my $pathway = "$opts{o}/pathway";
my $kegg_out = "$pathway/kegg_enrichment";
my $kegg_map = "$pathway/kegg_map";
#$opts{map} ||= "/share/nas2/database/KEGG/20150511/";
$opts{func} ||= "both";
$opts{i} = ABSOLUTE_DIR ($opts{i});
$opts{d} = ABSOLUTE_DIR ($opts{d});

my %file_db=%{readconf("$Bin/../database.config")}; 

if (!-d $opts{o}) {
	mkdir $opts{o} || die "Can't creat $opts{o},$!\n";
}
unless ($opts{func} eq "go") {
	if (!-d "$pathway") {
		mkdir "$pathway" || die "Can't creat $pathway,$!\n";
	}
	if (!-d "$pathway/kegg_map") {
		mkdir "$pathway/kegg_map" || die "Can't creat $pathway/kegg_map,$!\n";
	}
	if (!-d "$pathway/kegg_enrichment") {
		mkdir "$pathway/kegg_enrichment" || die "Can't creat $pathway/kegg_enrichment,$!\n";
	}
}
unless ($opts{func} eq "kegg") {
	if (!-d "$opts{o}/go_enrichment") {
		mkdir "$opts{o}/go_enrichment" || die "Can't creat $opts{o}/go_enrichment,$!\n";
	}
}

my ($path,$ko,$tree,$go_list);
opendir DH,$opts{i} || die "cannot open the $opts{i},$!\n";
for my $file (readdir DH) {
	$path      = $file if ($file =~ /\.path/);
	$ko        = $file if ($file =~ /\.ko/);
	$tree      = $file if ($file =~ /tree/);
	$go_list   = $file if ($file =~ /GO.list/);

}
closedir DH;

unless ($opts{func} eq "kegg") {
	if ( !defined($tree) || !defined($go_list) ) {
		die "\n\tPlease check input files!\n\n";
	}
}
unless ($opts{func} eq "go") {
	if ( !defined($path) ) {
		die "\n\tPlease check input files!\n\n";
	}
}

#################################################
#################### enrichment #################
#################################################

################# build reference ###############
my %deg; my $Marker;
open (UNI,"$opts{d}") || die "Can't open $opts{d},$!\n";
while (<UNI>) {
	chomp;
	next if ($. == 1);
	my $gene = (split /\t/,$_)[0];
	$deg{$gene} = 1;
	$Marker++;
}
if ($Marker<=3) {
    print "$opts{d} only have $Marker gene. Exit\n";
	system "touch $kegg_out/$opts{k}.KEGG.zero";
	system "touch $go_out/$opts{k}.GO.zero";
	exit;
}

unless ($opts{func} eq "kegg") 
{

	################### extract txt #################
	&extract ("$go_list","$opts{k}.GO.list.txt");
	
	################### GO distribution png #########
	chdir $go_out;
	print "perl $Bin/gene_ontology_graph_v1.2.pl -i $opts{i}/$go_list -i $opts{k}.GO.list.txt -mark All_transcript -mark DE_transcript -o $go_out -k $opts{k} -c $file_db{component_file} -f $file_db{function_file} -p $file_db{process_file}\n"	;
	system "perl $Bin/gene_ontology_graph_v1.2.pl -i $opts{i}/$go_list -i $opts{k}.GO.list.txt -mark All_transcript -mark DE_transcript -o $go_out -k $opts{k} -c $file_db{component_file} -f $file_db{function_file} -p $file_db{process_file} " ;                               #modified by sunhy  change to new script


	#################### GO #########################
	my %go;
	sub go_sub{
		my ($flag,$type)=@_;
		if ($flag eq "$type") {
			if (! keys %{$go{$flag}}){
				$go{$flag}{name}=$_;
			}else{
				my @t = split /\t/,$_;
				my @p = split /;/,($t[2]);
				my @tmp;
				my $num=0;
				for (@p) {
					$go{$flag}{all}{$_}=1;
					if (exists $deg{$_}) {
						push @tmp,$_;
						$num+=1;
						$go{$flag}{match}{$_}=1;
					}	
				}
				if (@tmp) {
					my $g = join (';',@tmp);
					$g = join ("\t",$t[0],$num,$t[1],$g);
					$go{$flag}{gene}{$g}='';
				}
			}
		}
	}
	open (GO,"$opts{i}/$tree");
	my $flag;
	while (<GO>) {
		next if /^\#/;
		$flag="Biological" if /Biological/;
		$flag="Molecular" if /Molecular/;
		$flag="Cellular" if /Cellular/;
		chomp;
		&go_sub($flag,"Biological");
		&go_sub($flag,"Molecular");
		&go_sub($flag,"Cellular");
	}
	$\= "\n";
	$, = "\t";
	foreach (sort keys %go) {
		open (OUT,">$go_out/$opts{k}.GO.$_.xls") || die "Can't creat $go_out/$opts{k}.GO.$_.xls,$!\n";
		my $gomatch=keys %{$go{$_}{match}};
		my $goall=keys %{$go{$_}{all}};

#		print OUT "\#$_";
#		print OUT "\#GO_n",$gomatch,"\#GO_N",$goall;
#		print OUT "\#Go_Term",'GO_ID','Unigene','gene','Unigene_all','gene_all';
#---------------------------------- by Simon Young 2015-03-04 --------------------------------------
        my $category = ($_ eq "Biological") ? "Biological Process" : 
                       ($_ eq "Molecular" ) ? "Molecular Function" : 
                                              "Cellular Component" ;
        my $header =<<"_HEADER_";
#! Input file to measure the differentially expressed gene(DEG)-enrichment in annotation $category GO terms by the right sided Fisher exact test.
## Term: GO term name and id, format like '$category: TERM NAME (TERM ID);'.
## mNum: DEG number in specified annotation $category GO term.
## nNum: Unigene number in specified annotation $category GO term.
## MNum: DEG number in all annotation $category GO term.
## NNum: Unigene number in all annotation $category GO term.
## GIDs: DEG IDs in specified annotation $category GO term, separated by ';'.
#Term	mNum	nNum	MNum	NNum	GIDs
_HEADER_
        $header =~s/\n$//; 
        print OUT $header;

		foreach (keys %{$go{$_}{gene}}) {
			my @tmp = (split /\t/,$_);
			print OUT @tmp[0..2],$gomatch,$goall,$tmp[3];
		}
		close OUT;
	}
	%go = ();
	$, = undef;
	$\ = undef;

	################ Fisher's exact test #############

	foreach my $stat_go (glob "$go_out/$opts{k}.GO.*.xls") {
		open (IN,"$stat_go") || die "Can't open $stat_go,$!\n";
		$stat_go =~s/.xls$//;
		open (OUT,">$stat_go.stat") || die "Can't creat $stat_go.stat,$!\n";
		$,="\t";
		$\= "\n";
		my (%hash,%new,$i);
		print OUT "\#".'Gene_Ontology_term','Cluter_frequency','Genome_frequency','P-value','Corrected_P-value';
		while (<IN>) {
			#next unless $. > 3;
            next if (/^#/);
			chomp;
			my ($term,$m,$M,$n,$N) = (split /\t/,$_)[0..4];
			next if $m <1;   
			my $p_value = &hyper($m,$M,$n,$N);
			my $m_ratio = "$m out of $n ".100 * $m / $n.'%';
			my $M_ratio = "$M out of $N ".100 * $M / $N.'%';
			$hash{$term}{name} = $term;
			$hash{$term}{p_value} = $$p_value;
			$hash{$term}{p_corr} = $$p_value;
			$hash{$term}{ratio}="$m_ratio\t$M_ratio";
			$i+=1;
		}
		for (sort {$hash{$a}->{p_value} <=> $hash{$b}->{p_value}}  keys %hash) {	#tag �ö�����ֵ����,�����½���ϣ
			$hash{$_}{p_corr} = ($hash{$_}{p_value} * $i >1) ? 1 : $hash{$_}{p_value} * $i;
			print OUT $hash{$_}->{name},$hash{$_}->{ratio},$hash{$_}{p_value},$hash{$_}->{p_corr};
		}		undef $,;
		undef $\;
		close IN;
		close OUT;
		%hash = ();
		%new = ();
	#	`rm $stat_go.xls`;
	}
}

######################## KEGG ####################

unless ($opts{func} eq "go")
{
	&extract ("$ko","$opts{k}.Kegg.ko");

	my %kegg;
	open (KEGG,"$opts{i}/$path") or die $!;
	while (<KEGG>) {
		next if /^\#/;
		chomp;
		my @t = split /\t/,$_;
		my @p = split /;/,($t[3]);
		my @q = split /\+/,($t[4]);
		my $num = 0;
		my (@p_new,@q_new);
		for my $i (0..$#p) {
			$kegg{all}{$p[$i]}=1;
			if (exists $deg{$p[$i]}) {
				push @p_new,$p[$i];
				push @q_new,$q[$i];
				$num += 1;
				$kegg{match}{$p[$i]}=1;
			}
		}
		next if (!@p_new);
		my $p = join (';',@p_new);
		my $q = join ('+',@q_new);
		my $g = join ("\t",$t[0],$t[1],$num,$t[2],$p,$q);
		$kegg{gene}{$g}='';
	}
	$\= "\n";
	$, = "\t";
	$kegg{num}{match}=keys %{$kegg{match}};
	$kegg{num}{all}=keys %{$kegg{all}};
	open (OUT,">$kegg_out/$opts{k}.KEGG.xls") || die "Can't creat $kegg_out/$opts{k}.KEGG.xls,$!\n";
	print OUT "\#KEGG_n",$kegg{num}{match},"KEGG_N",$kegg{num}{all};
	print OUT "\#Pathway",'ko_ID','Unigene','gene','Unigene_all','gene_all';
#---------------------------------- by Simon Young 2015-03-04 --------------------------------------
        my $header =<<"_HEADER_";
#! Input file to measure the differentially expressed gene(DEG)-enrichment in annotation KEGG pathways by the right sided Fisher exact test.
## MapName: KEGG pathway name.
## MapID: KEGG pathway ID.
## mNum: DEG number in specified annotation KEGG pathways.
## nNum: Unigene number in specified annotation KEGG pathways.
## MNum: DEG number in all annotation KEGG pathways.
## NNum: Unigene number in all annotation KEGG pathways.
## GIDs: DEG IDs in specified annotation KEGG pathways, separated by ';'.
## KIDs: DEG IDs in specified annotation KEGG pathways, separated by ';'.
#kMap	KPid	mNum	nNum	MNum	NNum	GIDs
_HEADER_
#        $header =~s/\n$//; 
#        print OUT $header;




	foreach (keys %{$kegg{gene}}) {
		my @tmp=(split/\t/,$_);
		print OUT @tmp[0..3],$kegg{num}{match},$kegg{num}{all},$tmp[4],$tmp[5];
	}
	$, = undef;
	$\= undef;
	%kegg = ();
	close OUT;

	################ Fisher's exact test #############

	open (IN,"$kegg_out/$opts{k}.KEGG.xls") || die "Can't open $kegg_out/$opts{k}.KEGG.xls,$!\n";
	open (OUT,">$kegg_out/$opts{k}.KEGG.stat") || die "Can't creat $kegg_out/$opts{k}.KEGG.stat,$!\n";
	open (OUT_T,">$kegg_out/$opts{k}.KEGG.stat.temp") || die "Can't creat $kegg_out/$opts{k}.KEGG.stat.temp,$!\n";
	my $Marker=0;
	$, = "\t";
	$\ = "\n";
	{
		my (%hash,$i);
		print OUT "\#Kegg_pathway",'ko_id','Cluter_frequency','Genome_frequency','P-value','Corrected_P-value';
		print OUT_T "\#Kegg_pathway",'ko_id','Cluter_frequency','Genome_frequency','P-value','Corrected_P-value','rich_factor';
		while (<IN>) {
				next if /^\#/;
				chomp;
				$Marker++;
				my ($term,$ko_id,$m,$M,$n,$N) = (split /\t/,$_)[0..5];
				next if $m <1;                                         ##reduce pathway which deg number is smaller than2;
				my $p_value = &hyper($m,$M,$n,$N);
				my $m_ratio = "$m out of $n ".100 * $m / $n.'%';
				my $M_ratio = "$M out of $N ".100 * $M / $N.'%';
				my $rich_factor= ($m / $n)/($M / $N);
				$hash{$term}{name} = $term;
				$hash{$term}{ko_id} = $ko_id;
				$hash{$term}{p_value} = $$p_value;
				$hash{$term}{p_corr} = $$p_value;
				$hash{$term}{ratio}="$m_ratio\t$M_ratio";
				$hash{$term}{rich_factor}=$rich_factor;
				$i+=1;
			}
		for (sort {$hash{$a}->{p_value} <=> $hash{$b}->{p_value}}  keys %hash) {	#tag �ö�����ֵ����,�����½���ϣ
			$hash{$_}{p_corr} = ($hash{$_}{p_value} * $i >1) ? 1 : $hash{$_}{p_value} * $i;
			print OUT $hash{$_}->{name},$hash{$_}->{ko_id},$hash{$_}->{ratio},$hash{$_}->{p_value},$hash{$_}->{p_corr};
			print OUT_T $hash{$_}->{name},$hash{$_}->{ko_id},$hash{$_}->{ratio},$hash{$_}->{p_value},$hash{$_}->{p_corr},$hash{$_}->{rich_factor};
		}
		undef $,;
		undef $\;
	%hash = ();
	}
	#`rm $opts{key}.KEGG.xls`;
	close IN;
	close OUT;
	close OUT_T;
       if($Marker==0){
	  system("touch $kegg_out/$opts{k}.KEGG.zero");
       }
	#################################################
	#################### kegg map  ##################
	#################################################

	#����������  hash-of-level form DEG.list ��������
	my %ko;
	open (UD,"$opts{d}") || die "Can't open $opts{d},$!\n";
	while (<UD>) {
		next if $.==1;
		chomp;
		my @tmp = split/\t/;
		$ko{$tmp[0]}{color} = "up";
	}
	close (UD);

	#��������������  relate-to-gene  ����������������
	my %png_new;
	open (LI,"$kegg_out/$opts{k}.KEGG.xls") || die "Can't open $kegg_out/$opts{k}.KEGG.xls,$!\n";
	while (<LI>) {
		next if $.<3;
		chomp;
		my @tmp_1 = split /\t/;
		my @tmp_2 = split /;/,$tmp_1[6];
		my @tmp_3 = split /\+/,$tmp_1[7];
		for my $i (0..$#tmp_2) {
			$png_new{$tmp_1[1]}{$tmp_3[$i]}{color}.="$ko{$tmp_2[$i]}{color};";	#tag��������Ƭ
			$png_new{$tmp_1[1]}{$tmp_3[$i]}{value}.="1;";
			$png_new{$tmp_1[1]}{$tmp_3[$i]}{gname}.="$tmp_2[$i];";
		}
		for my $i (0..$#tmp_2) {
			$png_new{$tmp_1[1]}{$tmp_3[$i]}{color} =~ s/;$//;	#tag��map����grep�޸� ��ɢ����Ƭ
			$png_new{$tmp_1[1]}{$tmp_3[$i]}{value} =~ s/;$//;
			$png_new{$tmp_1[1]}{$tmp_3[$i]}{gname} =~ s/;$//;
			$png_new{$tmp_1[1]}{$tmp_3[$i]}{COLOR} = "up" if ($png_new{$tmp_1[1]}{$tmp_3[$i]}{color} !~ /down/);
			$png_new{$tmp_1[1]}{$tmp_3[$i]}{COLOR} = "down" if ($png_new{$tmp_1[1]}{$tmp_3[$i]}{color} !~ /up/);
			$png_new{$tmp_1[1]}{$tmp_3[$i]}{COLOR} = "mix" if ( $png_new{$tmp_1[1]}{$tmp_3[$i]}{color} =~/up/ && $png_new{$tmp_1[1]}{$tmp_3[$i]}{color} =~/down/);
			$png_new{$tmp_1[1]}{$tmp_3[$i]}{COLOR} = "none" if ($png_new{$tmp_1[1]}{$tmp_3[$i]}{color} !~ /up/ && $png_new{$tmp_1[1]}{$tmp_3[$i]}{color} !~ /down/);
		}
	}
	%ko=();
	close (LI);

	#����������������  modify png  ������������������
	#����������������& creat html  ������������������
	my @lost;
	for my $koxx (sort keys %png_new) {
		if (-e "$file_db{kegg_png_file}/$koxx.png") {
			my @K_ID = keys %{$png_new{$koxx}};
#			print "@K_ID\n",'=' x 10,"\n"; #tag
			my %html;
			open (PNG,"$file_db{kegg_png_file}/$koxx.png") || die "Can' open $file_db{kegg_png_file}/$koxx.png,$!\n";
			open (CONF,"$file_db{kegg_png_file}/$koxx.conf") || die "Can't open $file_db{kegg_png_file}/$koxx.conf,$!\n";
			open (RESULT,">$kegg_map/$koxx.png") || die "Can't creat $kegg_map/$koxx.png,$!\n";
			open (HTML,">$kegg_map/$koxx.html") || die "Can't creat $kegg_map/$koxx.html,$!\n";
	#		my $im = GD::Image->new(*PNG);				#��truecolorģʽ�򿪣������ko00190��ɫȱʧ
			my $im = newFromPng GD::Image(\*PNG,1) ||die;		#truecolorģʽ
			my $red = $im->colorAllocate(255,0,0);
			my $green = $im->colorAllocate(0,255,0);
			my $blue = $im->colorAllocate(0,0,255);
			my $black = $im->colorAllocate(0,0,0);
			while (my $rect=<CONF>) {
				my ($type,$p1,$p2,$p3,$p4) = $rect =~ /^(\w+)\t\((\d+),(\d+),(\d+),(\d+)\)/;
				my ($up,$down,$mix)=(0,0,0);
				foreach my $K (@K_ID) {
					if ($rect =~ /$K/){
						$up += 1 if ($png_new{$koxx}{$K}{COLOR} eq 'up');
						$down += 1 if ($png_new{$koxx}{$K}{COLOR} eq 'down');
						$mix += 1 if ($png_new{$koxx}{$K}{COLOR} eq 'mix');
						
						if($type=~/rect/){
							$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{info}.="<ul><li>$K $png_new{$koxx}{$K}{gname} ($png_new{$koxx}{$K}{value})</li></ul>";
							$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{KO}.="+$K";
							$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{html_coord}="$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2";
						}
						elsif($type=~/line/){
							my @temp=split/\t/,$rect;
							for(my $i=1;$i<=$#temp-1;$i++){
								my ($p1,$p2,$p3,$p4) = $temp[$i] =~ /\((\d+),(\d+),(\d+),(\d+)\)/;
								my ($x1,$x2,$x3,$x4,$y1,$y2,$y3,$y4);
								if($p2==$p4){
									$y1=$p2-1;$y2=$p2+1;$y3=$p4-1;$y4=$p4+1;
									$html{$koxx}{"$p1,$p2,$p3,$p4"}{png_coord}="$p1,$p2,$p3,$p4,$p3,$y1,$p1,$y1,$p1,$p2";  ##png_coord for bold the png line 
									$html{$koxx}{"$p1,$p2,$p3,$p4"}{html_coord}="$p1,$y1,$p1,$y2,$p3,$y4,$p3,$y3,$p1,$y1";   ###html_coord for view the href of html
								}
								else{
									$x1=$p1-1;$x2=$p1+1;$x3=$p3-1;$x4=$p3+1;
									$html{$koxx}{"$p1,$p2,$p3,$p4"}{png_coord}="$p1,$p2,$p3,$p4,$x3,$p4,$x1,$p2,$p1,$p2";
									$html{$koxx}{"$p1,$p2,$p3,$p4"}{html_coord}="$x1,$p2,$x2,$p2,$x4,$p4,$x3,$p4,$x1,$p2";

								}
								$html{$koxx}{"$p1,$p2,$p3,$p4"}{info}.="<ul><li>$K $png_new{$koxx}{$K}{gname} ($png_new{$koxx}{$K}{value})</li></ul>";
								$html{$koxx}{"$p1,$p2,$p3,$p4"}{KO}.="+$K";
							}
						}
					}
				}
				if ($down == 0 && $mix == 0 && $up != 0) {
					#print "$up\t$down\t$mix\tred\t$p1\t$p2\t$p3\t$p4\n",$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{info},"\n" if (defined $opts{v});
					if($type=~/rect/){
						my $poly=new GD::Polygon;
						$poly->addPt($p1,$p2);$poly->addPt($p1,$p4);$poly->addPt($p3,$p4);$poly->addPt($p3,$p2);
						$im->openPolygon($poly,$black);
						#$im->filledPolygon($poly,$red);
						$im->fillToBorder($p1+1,$p2+1,$black,$red);  
						for my $x($p1+1..$p3-1){                    ##look for the close area and fill color
							for my $y ($p2+1..$p4-1){
								my $colorIndex=$im->getPixel($x,$y);
								my ($r,$g,$b)=$im->rgb($colorIndex);
								if($r==255 and $g ==255 and $b==255){
									$im->fillToBorder($x,$y,$black,$red);
								}
							}
						}
						$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{updown}="Up";
						$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{color}="f00";
					}
					elsif($type=~/line/){
						my @temp=split/\t/,$rect;
						for(my $i=1;$i<=$#temp-1;$i++){
							my ($p1,$p2,$p3,$p4) = $temp[$i] =~ /\((\d+),(\d+),(\d+),(\d+)\)/;
							#print "$p1,$p2,$p3,$p4\n$html{$koxx}{\"$p1,$p2,$p3,$p4\"}{png_coord}\n";
							my $poly=new GD::Polygon;
							my @new_coord=split/,/,$html{$koxx}{"$p1,$p2,$p3,$p4"}{png_coord};
							for (my $i=0;$i<=$#new_coord ;$i+=2) {
								$poly->addPt($new_coord[$i],$new_coord[$i+1]);
							}
							$im->openPolygon($poly,$red);
							$html{$koxx}{"$p1,$p2,$p3,$p4"}{updown}="Up";
							$html{$koxx}{"$p1,$p2,$p3,$p4"}{color}="f00";
						}
						
					}
				}
				if ($up == 0 && $mix == 0 && $down != 0) {
					#print "$up\t$down\t$mix\tgreen\t$p1\t$p2\t$p3\t$p4\n",$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{info},"\n" if (defined $opts{v});
					if($type=~/rect/){
						my $poly=new GD::Polygon;
						$poly->addPt($p1,$p2);$poly->addPt($p1,$p4);$poly->addPt($p3,$p4);$poly->addPt($p3,$p2);
						$im->openPolygon($poly,$black);
						$im->fillToBorder($p1+1,$p2+1,$black,$green);  
						for my $x($p1+1..$p3-1){                    ##look for the close area and fill color
							for my $y ($p2+1..$p4-1){
								my $colorIndex=$im->getPixel($x,$y);
								my ($r,$g,$b)=$im->rgb($colorIndex);
								if($r==255 and $g ==255 and $b==255){
									$im->fillToBorder($x,$y,$black,$green);
								}
							}
						}
						$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{updown}="Down";
						$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{color}="0f0";
					}
					elsif($type=~/line/){
						my @temp=split/\t/,$rect;
						for(my $i=1;$i<=$#temp-1;$i++){
							my ($p1,$p2,$p3,$p4) = $temp[$i] =~ /\((\d+),(\d+),(\d+),(\d+)\)/;
							my $poly=new GD::Polygon;
							my @new_coord=split/,/,$html{$koxx}{"$p1,$p2,$p3,$p4"}{png_coord};
							for (my $i=0;$i<=$#new_coord ;$i+=2) {
								$poly->addPt($new_coord[$i],$new_coord[$i+1]);
							}
							$im->openPolygon($poly,$green);
							$html{$koxx}{"$p1,$p2,$p3,$p4"}{updown}="Down";
							$html{$koxx}{"$p1,$p2,$p3,$p4"}{color}="0f0";
						}
						
					}
				}
				if ($mix != 0 || ($up*$down) != 0) {
					#print "$up\t$down\t$mix\tblue\t$p1\t$p2\t$p3\t$p4\n",$html{$koxx}{"$p1,$p2,$p3,$p4"}{info},"\n" if (defined $opts{v});
					if($type=~/rect/){
						my $poly=new GD::Polygon;
						$poly->addPt($p1,$p2);$poly->addPt($p1,$p4);$poly->addPt($p3,$p4);$poly->addPt($p3,$p2);
						$im->openPolygon($poly,$black);
						$im->fillToBorder($p1+1,$p2+1,$black,$blue);  
						for my $x($p1+1..$p3-1){                    ##look for the close area and fill color
							for my $y ($p2+1..$p4-1){
								my $colorIndex=$im->getPixel($x,$y);
								my ($r,$g,$b)=$im->rgb($colorIndex);
								if($r==255 and $g ==255 and $b==255){
									$im->fillToBorder($x,$y,$black,$blue);
								}
							}
						}
						$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{updown}="Mix";
						$html{$koxx}{"$p1,$p2,$p1,$p4,$p3,$p4,$p3,$p2"}{color}="00f";
					}
					elsif($type=~/line/){
						my @temp=split/\t/,$rect;
						for(my $i=1;$i<=$#temp-1;$i++){
							($p1,$p2,$p3,$p4) = $temp[$i] =~ /\((\d+),(\d+),(\d+),(\d+)\)/;
							my $poly=new GD::Polygon;
							my @new_coord=split/,/,$html{$koxx}{"$p1,$p2,$p3,$p4"}{png_coord};
							for (my $i=0;$i<=$#new_coord ;$i+=2) {
								$poly->addPt($new_coord[$i],$new_coord[$i+1]);
							}
							$im->openPolygon($poly,$blue);
							#$im->filledRectangle($p1,$p2,$p3,$p4,$blue);
							$html{$koxx}{"$p1,$p2,$p3,$p4"}{updown}="Mix";
							$html{$koxx}{"$p1,$p2,$p3,$p4"}{color}="00f";
						}
						
					}
				}
			}
			##########################
			######	HTML BEGIN	######
			##########################
			print HTML <<____________HTML;
<html>
<meta http-equiv="content-type" content="text/html; charset=utf-8">
<head>
<title>
$koxx
</title>
<style type="text/css">
<!--

area {cursor: pointer;}

-->
</style>
<link rel="stylesheet" href="/css/kegg.css" type="text/css" />
<script language="JavaScript" src="/js/dhtml.js"></script>
<script type="text/javascript">
<!---

function showInfo(info) {
	obj = document.getElementById("result");
	obj.innerHTML = "<div style='cursor: pointer; position: absolute; right: 5px; color: #000;' onclick='javascript: document.getElementById(\\"result\\").style.display = \\"none\\";' title='close'>X</div>" + info;
	obj.style.top = document.body.scrollTop;
	obj.style.left = document.body.scrollLeft;
	obj.style.display = "";
}

//--->
</script>

</head>
<body>
<img src="$koxx.png" usemap="#mapdata" border="0" />
<map name="mapdata">
____________HTML
			foreach my $coords ( keys %{$html{$koxx}} ) {
				my $color = $html{$koxx}{$coords}{color};
				my $updown = $html{$koxx}{$coords}{updown};
				my $info = $html{$koxx}{$coords}{info};
				print "$koxx  $coords\n" and exit if (!defined $info or !defined $color or !defined $updown);
				if (!($updown eq 'Mix')) {
					print HTML <<____________HTML;
<area shape='poly' coords='$html{$koxx}{$coords}{html_coord}' href='http://www.kegg.jp/dbget-bin/www_bget?$html{$koxx}{$coords}{KO}' target="_blank" onmouseover='javascript: showInfo("<ul><li style=\\"color: #$color;\\">$updown regulated$info</li></ul>");' />
____________HTML
				}else{
					my ($up,$down) = (&mix_info($info)->[0],&mix_info($info)->[1]);
					print HTML <<____________HTML;
<area shape='poly' coords='$html{$koxx}{$coords}{html_coord}' href='http://www.kegg.jp/dbget-bin/www_bget?$html{$koxx}{$coords}{KO}' target="_blank" onmouseover='javascript: showInfo("<ul><li style=\\"color: #$color;\\">$updown regulated$up$down</li></ul>");' />
____________HTML

				}
			}

			print HTML <<____________HTML;
</map>
<div id='result' style='position: absolute; width: 50%; border: 1px solid #000; background-color: #fff; filter: alpha(opacity=95); opacity: 0.95; font-size: 12px; padding-right: 20px; display: none;' onmouseover="javascript: this.style.filter = 'alpha(opacity=100)'; this.style.opacity = 1;" onmouseout="javascript: this.style.filter = 'alpha(opacity=95)'; this.style.opacity = 0.95;"></div>
</body>
</html>
____________HTML
			##########################
			######	HTML END	######
			##########################
			close (HTML);
			close (PNG);
			close (CONF);
			binmode RESULT;
			print RESULT $im->png;
			close (RESULT);
			print $koxx,"\n" if (defined $opts{v});
			print "job done!\n" if (defined $opts{v});
			print '=' x 40,"\n" if (defined $opts{v});
		} else {
			print "$koxx map not found in kegg datebase!\n" if (defined $opts{v});
			print '=' x 40,"\n" if (defined $opts{v});
			push @lost,"$koxx\n";
		}
		
	}
	print "[notice]:These maps are not in datebase\n",@lost;
}
#��������������  program  end  ������������������

my $Time_End = sub_format_datetime(localtime(time()));
&Runtime($BEGIN);
print "\nEnd Time :[$Time_End]\n\n";

#������������������������������������������������
#                  Sub Routines                ��
#������������������������������������������������

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


sub extract{
	my $var = shift;
	my $pre_key = shift;
	open (IN,"$opts{i}/$var") || die "Can't open $opts{i}/$var,$!\n";
	if ($var=~/ko/) {
		open (OUT,">$kegg_out/$pre_key") || die "Can't creat $go_out/$var,$!\n";
	}
	else {
		open (OUT,">$go_out/$pre_key") || die "Can't creat $go_out/$var,$!\n";
	}
	while (<IN>) {
		if (/^\#/) {
			print OUT $_;
			next;
		}
		my $gene = (split /\t/,$_)[0];
		if (exists $deg{$gene}) {
			print OUT $_;
		}
	}
	close (IN);
	close (OUT);
}
sub hyper{ #calculate  Fisher's exact test
	my ($n11,$n1p,$np1,$npp)=@_;
	my ($out,$errorCode,$right_value);
	$right_value = calculateStatistic(
		n11=>$n11,
		n1p=>$n1p,
		np1=>$np1,
		npp=>$npp
	);
	if( ( $errorCode = getErrorCode() ) ) {
		$out = $errorCode." - ".getErrorMessage();
	}
	else {
		$out = $right_value;
	}
	return \$out;
}
sub mix_info{	#For HTML mix color info
	my $text = shift;
	my %info;
	my @messages = $text =~ /<ul><li>(.+?)<\/li><\/ul>/g;
	for my $mess (@messages) {
		my @s = split / /,$mess;
		my @g = split /;/,$s[1];
		my @v = (split /[\(;\)]/,$s[2],-1);
		shift @v,pop @v;
		for my $i (0..$#v) {
			if ($v[$i] > 0) {
				$info{up}{$s[0]}{ID}.= "$g[$i];";
				$info{up}{$s[0]}{VA}.= "$v[$i];";
			}else{
				$info{down}{$s[0]}{ID}.= "$g[$i];";
				$info{down}{$s[0]}{VA}.= "$v[$i];";
			}
		}
	}
	my ($message_up,$message_down);
	for my $key_1 ( keys %info ) {
		if ($key_1 eq 'up') {
			for my $key_2 (sort keys %{$info{$key_1}}) {
				$info{$key_1}{$key_2}{ID} =~ s/;$//;
				$info{$key_1}{$key_2}{VA} =~ s/(.*);$/\($1\)/;
				my $inner_message = "<ul><li>$key_2 $info{$key_1}{$key_2}{ID} $info{$key_1}{$key_2}{VA}</li></ul>";
				$message_up.=$inner_message;
			}
		}
		if ($key_1 eq 'down') {
			for my $key_2 (sort keys %{$info{$key_1}}) {
				$info{$key_1}{$key_2}{ID} =~ s/;$//;
				$info{$key_1}{$key_2}{VA} =~ s/(.*);$/\($1\)/;
				my $inner_message = "<ul><li>$key_2 $info{$key_1}{$key_2}{ID} $info{$key_1}{$key_2}{VA}</li></ul>";
				$message_down.=$inner_message;
			}
		}
	}
	$message_up = defined $message_up ? "<ul><li>Up $message_up</li></ul>" : "<ul><li>Up </li></ul>" ;
	$message_down = defined $message_down ? "<ul><li>Down $message_down</li></ul>" : "<ul><li>Down </li></ul>" ;
	return [($message_up,$message_down)];
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
sub help {#usage
	print <<"	Usage End.";
	Description:
		KEGG GO enrichment & maps & web pages
		Version: $ver

	Usage:
		
		-d    <str>      DEG file
		-k    <str>      output files prefix
		-i    <str>      input dir
		-o    <str>      output dir
		-func            [option]
		      kegg       only kegg enrichment maps web;
		      go         only go enrichment
		      both       [default]
 		-v               [option]   view process on screen
	
	Example��

		perl KeggGo_enrich_map_web.pl -d xxx.DEG.final.xls -k demo -i input -o output -func kegg 

	[notice]input dir should at least contain:
			
		xxx.fa.Kegg.path                         [for kegg]
		xxx.fa.Kegg.ko
		or
		xxx.pep.fa.GO_tree*              [for  go ]
		xxx.pep.fa.GO.list*

	Usage End.

	exit;
}
