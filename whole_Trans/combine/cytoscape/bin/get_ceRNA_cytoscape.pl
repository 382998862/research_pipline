#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($idir,$odir,$data_cfg,$detail_cfg,$ce);

GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$odir,
				"in:s"=>\$idir,
				"ce:s"=>\$ce,
				"data_cfg:s"=>\$data_cfg,
				"detail_cfg:s"=>\$detail_cfg,
				) or &USAGE;
&USAGE unless ($idir and $odir and $data_cfg and $detail_cfg and $ce);
# ------------------------------------------------------------------

$data_cfg = &ABSOLUTE_DIR($data_cfg);
$detail_cfg = &ABSOLUTE_DIR($detail_cfg);
$idir = &ABSOLUTE_DIR($idir);
`mkdir $odir` unless (-d "$odir");
$odir = &ABSOLUTE_DIR($odir);
`mkdir $odir/work_sh` unless (-d "$odir/work_sh");

my $work_sh = "$odir/work_sh";
$ce = &ABSOLUTE_DIR($ce);
print("$ce");

my %sample = &relate($data_cfg);
my %config = ();
&readcfg($detail_cfg);
my @diff;my %vs;
open (CFG,"$detail_cfg") or die $!;
print ("read cfgi and diff group is:\n");
while (<CFG>){
	chomp;
	next if(/^#|^$|^\s$/);
	my @tmp = split /\s+/,$_;
	if($tmp[0] eq "Sep"){
		$tmp[1]=~s/\;/_vs_/;
		$tmp[1]=~s/,/_/g;
		push @diff,$tmp[1];
		print "$tmp[1]\n";
	}
	if($tmp[0] eq "Com"){
		$tmp[1]=~s/,/_vs_/;
		push @diff,$tmp[1];
		print "$tmp[1]\n";
	}
}
close(CFG);

my @RNA = ("gene","lncRNA","circRNA","sRNA");
my (%h,%LSM,%CSM);
#$vs{$group}{$rna}=$Group	$vs{W01_vs_W02}{$gene}=L01_vs_L02
foreach my $k1 (@diff){
	&switch($k1);
	foreach my $k2 (@RNA){
		$h{$k1}{$k2}= "$idir/$k2/$vs{$k1}{$k2}/$vs{$k1}{$k2}.DEG_final.xls";
		print "$idir/$k2/$vs{$k1}{$k2}/$vs{$k1}{$k2}.DEG_final.xls\n";
	}
}

open(SH,">$work_sh/ce.cytoscape.sh");
foreach my $vs (@diff){
	system ("cat $h{$vs}{gene} $h{$vs}{lncRNA} $h{$vs}{sRNA} > $odir/$vs.lncRNA-miRNA-mRNA.all.DEG.xls");
	print "cat $h{$vs}{gene} $h{$vs}{lncRNA} $h{$vs}{sRNA} > $odir/$vs.lncRNA-miRNA-mRNA.all.DEG.xls\n";
	system ("cat $h{$vs}{gene} $h{$vs}{circRNA} $h{$vs}{sRNA} > $odir/$vs.circRNA-miRNA-mRNA.all.DEG.xls");
	print "cat $h{$vs}{gene} $h{$vs}{circRNA} $h{$vs}{sRNA} > $odir/$vs.circRNA-miRNA-mRNA.all.DEG.xls\n";
	`mkdir $odir/lncRNA-miRNA-mRNA` unless (-d "$odir/lncRNA-miRNA-mRNA");
	`mkdir $odir/circRNA-miRNA-mRNA` unless (-d "$odir/circRNA-miRNA-mRNA");
	print SH "perl $Bin/abstract.pl -deg $odir/$vs.lncRNA-miRNA-mRNA.all.DEG.xls -ce $ce -o $odir/lncRNA-miRNA-mRNA/$vs.lncRNA-miRNA-mRNA.ceRNA_pair_adjust_p_Sig_diff.xls\n";
	print SH "perl $Bin/abstract.pl -deg $odir/$vs.circRNA-miRNA-mRNA.all.DEG.xls -ce $ce -o $odir/circRNA-miRNA-mRNA/$vs.circRNA-miRNA-mRNA.ceRNA_pair_adjust_p_Sig_diff.xls\n";
}
close(SH);
&qsub("$work_sh/ce.cytoscape.sh");
######################################################################################l
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub qsub(){
        my $shfile= shift;
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $shfile --queue $config{Queue_type}";
        &run_or_die($cmd);
}
sub run_or_die(){
        my ($cmd) = @_ ;
        &show_log($cmd);
        my $flag = system($cmd) ;
        if ($flag != 0){
                &show_log("Error: command fail: $cmd");
                exit(1);
        }
        &show_log("done.");
}
sub show_log(){
        my ($txt) = @_ ;
        my $time = time();
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
        $wday = $yday = $isdst = 0;
        my $Time=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
        print "$Time:\t$txt\n" ;
}
########################
sub readcfg{
        my $cfg=shift;
        open(CFG,$cfg)||die $!;
        while(<CFG>){
                next if($_=~/^#|^\s$/);
                chomp;my @tmp=split(/\s+/,$_);
                if($tmp[0] eq "Com"||$tmp[0] eq "Sep"){
                        push @{$config{$tmp[0]}},$tmp[1];
                }else{
                        $config{$tmp[0]}=$tmp[1];
                }
        }
        close(CFG);
}
#######################
sub relate{	
        my $file=shift;
        open(REL,$file)||die $!;
        my %sample=();
        while(<REL>){
                chomp;
                next if($_ !~/^SampleID/);
                my @tmp=split(/\s+/,$_);
                $sample{$tmp[-1]}{sRNA}=$tmp[2];
                $sample{$tmp[-1]}{lncRNA}=$tmp[3];
                $sample{$tmp[-1]}{circRNA}=$tmp[3];
                $sample{$tmp[-1]}{gene}=$tmp[3];
        }
        close(REL);
        return %sample;
}

sub switch {
	my $group = shift;
	my ($v1,$v2,@V1,@V2);
	
		($v1,$v2) = split /_vs_/,$group,2;
		@V1 = split /_/,$v1;
		@V2 = split /_/,$v2;
		foreach my $key ("sRNA","lncRNA","circRNA","gene"){
			my (@s1,@s2);
			for(my $i=0;$i<@V1;$i++){
				push @s1,$sample{$V1[$i]}{$key};
			}
			for(my $j=0;$j<@V2;$j++){
				push @s2,$sample{$V2[$j]}{$key};
			}
			my $S1 = join ("_",@s1);
			my $S2 = join ("_",@s2);
			my $Group = join("_vs_",$S1,$S2);
			$vs{$group}{$key}=$Group;
		}
}

sub getCyto {
	my $inter = shift;
	my $att = shift;
	my $line = shift;
	my $Cyto = shift;
	open (O1,">$inter");
	open (O2,">$att");
	my %ind;
	my ($center,$type)=split //,$line;
	foreach my $rna(sort keys %$Cyto){
		foreach my $stat (sort keys %{$$Cyto{$rna}}){
			if (!exists $ind{$rna}){
				print O2 "$rna\t$center\t$stat\n";
				$ind{$rna}=1;
			}
			my @p = split /\;/,$$Cyto{$rna}{$stat};
			for (my $i=0;$i<@p;$i++){
				my ($id,$ud,$k) = split /\(|\)/,$p[$i];
				print O1 "$rna\t$id\t$line\n";
				if (!exists $ind{$id}){
					print O2 "$id\t$type\t$ud\n";
					$ind{$id}=1;
				}
			}
		}
	}
	close(O1);close(O2);
	return ($inter,$att);
}

sub Veen1{
        my ($List,$name,$id)=@_;
	my $dir = dirname $name;
        my $prefix = basename $name;
        my (%info,%venny, %com);
        open (SET, ">$name.degset.xls") or die $!;
        print SET "#DEG_Set\tDEG_Num\tDEG_IDs\n";
        my @label =@$id;
        my @list =@$List;
        for my $i (1..@list){
                my @ids;

                for my $deg ($list[$i-1]){
                        @ids = split /\s/,$deg;
                        for (my $j=0;$j<@ids;$j++){
                                my $deg_id = $ids[$j];
                                $info{$deg_id}{$label[$i-1]} = 1;
                        }
                }
                my ($id_num,$ids) = ($#ids+1, (join ";",@ids));
                print SET "$label[$i-1]\t$id_num\t$ids\n";
        }

        close SET;
#	open (OUT,">$name.veen.r");
	
	
}

sub Veen{
	my ($List,$name,$id)=@_;
	my $dir = dirname $name;
	my $prefix = basename $name;
	my (%info,%venny, %com);
	open (SET, ">$name.degset.xls") or die $!;
	print SET "#DEG_Set\tDEG_Num\tDEG_IDs\n";
	my @label =@$id;
	my @list =@$List;
	for my $i (1..@list){
		my @ids;
		
		for my $deg ($list[$i-1]){
			@ids = split /\s/,$deg;
			for (my $j=0;$j<@ids;$j++){
				my $deg_id = $ids[$j];
				$info{$deg_id}{$label[$i-1]} = 1;
			}
		}
		my ($id_num,$ids) = ($#ids+1, (join ";",@ids));
		print SET "$label[$i-1]\t$id_num\t$ids\n";
	}

	close SET;
	for my $e (sort keys %info) {
		my $com = join ",",(sort keys %{$info{$e}});
		$com{$com}++;
		$venny{$com}{$e} = 1;
	}
	open (VENN, ">$name.vennset.xls");
	print VENN "#Venn_Set\tElement_Num\tElement_IDs\n";

	for my $s (sort keys %venny) {
		my $elements = join ";",(sort keys %{$venny{$s}});
		print VENN "$s\t$com{$s}\t$elements\n";
	}
	my @color=("'cornflowerblue'","'green'","'yellow'","'darkorchid1'","'red'");
	my %DEG;
	my $list_content;
	my $label_content;
	my $color_content;
	for my $i (1..@list) {
		my @id;
		for my $deg ($list[$i-1]){
			@id = split /\s/,$deg;
			for(my $j=0;$j<@id;$j++){
				my $deg_id = "'".$id[$j]."'";
			push @{$DEG{$label[$i-1]}}, $deg_id;
			}
		}
	}
	my @name=qw(A B C D E F G M);
	my $name_content;
	for my $i (0..@label-1) {
		$list_content.= "$name[$i] <- c(".(join ", ",@{$DEG{$label[$i]}}).")\n";
		$label_content.= "$name[$i] = $name[$i], ";
		$name_content.="\"$label[$i]\",";
		$color_content.= "$color[$i], ";
	}
	$list_content =~ s/\n$//;
	$label_content =~ s/, $//;
	$color_content =~ s/, $//;
	$name_content=~s/,$//;
	my $more_opts = "";
	if (@label == 5) {
		$more_opts.= "    cex = 0.5,\n";
		$more_opts.= "    cat.cex = 0.6,\n";
		$more_opts.= "    margin = 0.1,\n";
		$more_opts.= "    cat.dist = c(0.20, 0.25, 0.20, 0.20, 0.25),\n";
		$more_opts.= "    scaled = FALSE,\n";
	}
	elsif (@label == 4) {
		$more_opts.= "    cex = 0.6,\n";
		$more_opts.= "    cat.cex = 0.7,\n";
		$more_opts.= "    margin = 0.08,\n";
		$more_opts.= "    scaled = FALSE,\n";
	}
	elsif (@label == 3) {
		$more_opts.= "    cex = 0.7,\n";
		$more_opts.= "    cat.cex = 0.7,\n";
		$more_opts.= "    margin = 0.06,\n";
		$more_opts.= "    euler.d = FALSE,\n";
		$more_opts.= "    cat.pos = c(0,0,180),\n";
		$more_opts.= "    scaled = FALSE,\n";
	}
	elsif (@label == 2) {
		$more_opts.= "    cex = 0.8,\n";
		$more_opts.= "    cat.cex = 0.5,\n";
		$more_opts.= "    margin = 0.05,\n";
		$more_opts.= "    cat.pos = 180,\n";
		$more_opts.= "    euler.d = TRUE,\n";
		$more_opts.= "    scaled = FALSE,\n";
	}
my $R_script = << "EOF";
#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript
#import library
library(grid)
library(VennDiagram)
#init deg lists
$list_content
lst=list($label_content)
names(lst)=c($name_content)
#plot venn
venn.diagram(
	x = lst,
	filename = "$prefix.venn.png",
	imagetype = "png",
	fill = c($color_content),
	height=1300,
	width=1300,
	resolution=200,
	units='px',
	wd=1,
	$more_opts
);
	
EOF
	open (RS, ">$name.venn.r") or die;
		print RS $R_script;
	close RS;

	system "cd $dir && /share/nas2/genome/biosoft/R/3.1.1/bin/Rscript $name.venn.r";
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

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
 ProgramName:
     Version:	$version
     Contact:	Niepy <niepy\@biomarker.com.cn> 
Program Date:	2018.05.23
      Modify:	
 Description:	This program is used to abstract diff.lncRNA-miRNA-mRNA and diff.circRNA-miRNA-mRNA form ceRNA result
       Usage:
		Options:
		-in <dir>	input dir ,DEG_Analysis,force
		-od <dir>	output dir , Combine/Cytoscape,force
		-ce <file>	ceRNA result,Coexpression_ceRNA_pair_adjust_p_Sig.txt(sample>=5) or ceRNA_pair_adjust_p_Sig.txt
		-data_cfg <file>	data_cfg,force
		-detail_cfg <file>	detail_cfg,force
		-h		help

USAGE
	print $usage;
	exit;
}
