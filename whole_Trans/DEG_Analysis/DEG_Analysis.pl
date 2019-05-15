use Getopt::Long;
use Getopt::Std;
use Config::General;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
use Cwd qw(abs_path getcwd);
my ($od,$conf,$count,$fpkm,$target,$cis_target,$trans_target,$anno,$type,$ppi,$human,$newgene_gff);
GetOptions(
        "h|?"           =>\&USAGE,
	"od:s"     	=>\$od,
	"count:s"	=>\$count,
	"fpkm:s"	=>\$fpkm,
	"conf:s"	=>\$conf,
	"ppi:s"		=>\$ppi,
	"target:s"	=>\$target,
	"cis:s"		=>\$cis_target,
	"trans:s"	=>\$trans_target,
	"anno:s"	=>\$anno,
	"type:s"	=>\$type,
	"gff:s"		=>\$newgene_gff,
	"human:s"	=>\$human,
)or &USAGE;
&USAGE unless ($od and $conf and $human);

`mkdir -p $od`	unless(-d $od);
`mkdir -p $od/work_sh`  unless(-d "$od/work_sh");

$od	=abs_path($od);
$count	=abs_path($count);
$fpkm	=abs_path($fpkm);
$conf	=abs_path($conf);
$anno	=abs_path($anno)	if(defined $anno);
$target	=abs_path($target)	if(defined $target);
$cis_target	=abs_path($cis_target)		if(defined $cis_target);
$trans_target	=abs_path($trans_target)	if(defined $trans_target);

my %config=&readConfig($conf);

$type ||="FPKM";
my $Rscript="/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";
my $rscript="/share/nas2/genome/biosoft/R/3.3.2/bin/Rscript";
#/share/nas2/genome/biosoft/R/3.3.2/bin/Rscript

`mkdir -p $od/density`	unless(-d "$od/density");
#################fpkm density 
&run_or_die("$Rscript  $Bin/cor/fpkm_density_plot_func.r $fpkm $count $od/density $type");

my @groups=&getvs($conf);
if(@groups==0){
	print "There is not deg group!\n";
	exit;
}
##################sample corelation
&run_or_die("$Rscript $Bin/cor/fpkm_cor.r --infile $fpkm --outdir $od/density --show_rownames --show_colnames --cluster_rows --cluster_cols --legend");

#################DEG Analysis
&run_or_die("perl $Bin/bin/gene_counts_to_FDR_FC.pl -count $count -conf $conf -fpkm $fpkm -type $type -human $human -od $od");

#################All DEG
`mkdir $od/All_DEG`	unless(-d "$od/All_DEG");
my @degs=glob("$od/*_vs_*/*_vs_*DEG_final.xls");
&getAllDEG(\@degs,$fpkm,"$od/All_DEG/All.DEG_final.xls");
&run_or_die("$Rscript $Bin/bin/draw/pheatmap.r --infile $od/All_DEG/All.DEG_final.xls --outfile $od/All_DEG/All.DEG_final --scale none --color.type 1 --legend  --is.log  --show_colnames --cluster_rows ");

my @deg;
foreach my $zero (@degs){
	my $line = `wc -l $zero`;chomp $line;
	if($line >1){
		push @deg,$zero;
	}
}
if(@deg>1 && @deg<=9){
	my $cmd="$Rscript $Bin/vn.r ".join(",",@deg)." $od/All_DEG  Veen";
	&run_or_die($cmd);		
}

#################Anno enrichment
if(!defined $anno){
	print "No anno dir!\n";
	exit;
}

my $go_anno=(glob("$anno/*.GO.anno.txt"))[0];
my $kegg_anno=(glob("$anno/*.Kegg.pathway"))[0];
my $inter_anno="$anno/Integrated_Function.annotation.xls";

&run_or_die("perl $Bin/enrich/bin/process_Kegg.pathway.pl $kegg_anno $od/KEGG.info")	unless(-e "$od/KEGG.info");
&run_or_die("perl $Bin/enrich/bin/process_GO.pl $go_anno $od/GO.info")			unless(-e "$od/GO.info");
if(!defined $target && !defined $cis_target && !defined $trans_target){
	&run_or_die("perl $Bin/enrich/anno/get_info_from_backgroud.pl $od/All_DEG/All.DEG_final.xls $inter_anno $od/All_DEG/All.DEG_final_anno.xls");

####################TF_analysis
#
	if(exists $config{medical}){
		######TF_activity
		if(exists $config{TFDB}){
			my $lines=`wc -l $od/All_DEG/All.DEG_final.xls`;
        	        my $num=(split(/\s+/,$lines))[0];
                	if($num <= 700) {
				open(TF,">$od/work_sh/tfs_activity.sh")||die $!;
	                        `mkdir -p $od/TF_activity` unless (-d "$od/TF_activity");
        	                my $tf_cmd = "perl $Bin/TF_analysis/TF_activity/bin/TF_activity.pl -fpkm $od/All_DEG/All.DEG_final.xls -cfg $conf -od $od/TF_activity ";
                        	print TF "$tf_cmd\n";
	                        close(TF);
				my $qsubcmd="sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $od/work_sh/tfs_activity.sh --independent --queue medical.q ";
				system($qsubcmd);
			}
		}
		######TFBS
		if(exists $config{spe_id}){
			open(TFSH,">$od/work_sh/tfbs_analysis.sh")||die $!;
			my $tf_cmd;
			if(exists $config{score}){
				my $score=$config{score};
	        	        $tf_cmd="perl $Bin/TF_analysis/TFBStools/bin/TFBS_predict.pl -cfg $conf -score $score -od $od/TFBS_Analysis -gff2 $newgene_gff ";
			}else{
				$tf_cmd="perl $Bin/TF_analysis/TFBStools/bin/TFBS_predict.pl -cfg $conf -od $od/TFBS_Analysis -gff2 $newgene_gff ";
			}
                	print TFSH "$tf_cmd\n";
        	        close(TFSH);
			my $qsubcmd="sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $od/work_sh/tfbs_analysis.sh --independent --queue medical.q ";
			&run_or_die($qsubcmd);
		
			#open(TFSH,">$od/work_sh/tfs_analysis.sh")||die $!;
			#`mkdir -p $od/TF_activity` unless (-d "$od/TF_activity");
				#my $tf_cmd="perl $Bin/TF_analysis/TF_activity/bin/TF_activity.pl -fpkm $od/All_DEG/All.DEG_final.xls -cfg $conf -od $od/TF_activity ";
			#$tf_cmd.="&& perl $Bin/TF_analysis/TFBStools/bin/TFBS_predict.pl -cfg $conf -od $od/TFBS_Analysis -gff2 $newgene_gff ";
			#print TFSH "$tf_cmd\n";
			#close(TFSH);
			#my $qsubcmd="sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $od/work_sh/tfs_analysis.sh --independent --queue medical.q ";
			#&run_or_die($qsubcmd);
			##########
			#DEG_TFBS
			`cut -f1 $od/All_DEG/All.DEG_final.xls > $od/TFBS_Analysis/DEG.list`;
			`mkdir -p $od/TFBS_Analysis/each_DEgeneRes` unless (-d "$od/TFBS_Analysis/each_DEgeneRes");
			`mkdir -p $od/TFBS_Analysis/DEG_seqLogo` unless (-d "$od/TFBS_Analysis/DEG_seqLogo");
	
			my $TFBSdb="/share/nas1/lijj/develop/TFBStools/$config{Project_key}"."_TFBS_"."$config{medical}";
			my @files=();
			my @pdfs=();
                        my @pngs=();
			if (-d $TFBSdb) {
				push @files,glob("$TFBSdb/each_geneRes/*predictRes.txt");
	                	push @files,glob("$od/TFBS_Analysis/each_geneRes/*predictRes.txt");
				push @pdfs,glob("$TFBSdb/seqLogo/*.pdf");
                                push @pngs,glob("$TFBSdb/seqLogo/*.png");
                                push @pdfs,glob("$od/TFBS_Analysis/seqLogo/*.pdf");
                                push @pngs,glob("$od/TFBS_Analysis/seqLogo/*.png");
			}else{
        	        	push @files,glob("$od/TFBS_Analysis/each_geneRes/*predictRes.txt");
				push @pdfs,glob("$od/TFBS_Analysis/seqLogo/*.pdf");
                                push @pngs,glob("$od/TFBS_Analysis/seqLogo/*.png");
			}
		
			open(LIST,"$od/TFBS_Analysis/DEG.list")||die $!;
			my %res=();
			my %pdfres=();
                        my %pngres=();
			foreach(my $i=0;$i<@files;$i++) {
			        my $f=$files[$i];$f=abs_path($f);
		        	my $gene=basename($f);
       			 	my $geneID=(split /\_TFBS/,$gene)[0];
				$res{$geneID}=$f;
			}
			foreach(my $i=0;$i<@pdfs;$i++) {
                                my $f=$pdfs[$i];$f=abs_path($f);
                                my $gene=basename($f);
                                my $geneID=(split /\_TFBS/,$gene)[0];
                                $pdfres{$geneID}=$f;
                        }
			foreach(my $i=0;$i<@pngs;$i++) {
                                my $f=$pngs[$i];$f=abs_path($f);
                                my $gene=basename($f);
                                my $geneID=(split /\_TFBS/,$gene)[0];
                                $pngres{$geneID}=$f;
                        }
			while(<LIST>){
				chomp;next if(/^#/);
				my $gene=(split /\s+/,$_)[0];
				if(exists $res{$gene}) {
					`cp $res{$gene} $od/TFBS_Analysis/each_DEgeneRes`;
				}
				if(exists $pdfres{$gene}) {
                                        `cp $pdfres{$gene} $od/TFBS_Analysis/DEG_seqLogo`;
                                }
                                if(exists $pngres{$gene}) {
                                        `cp $pngres{$gene} $od/TFBS_Analysis/DEG_seqLogo`;
                                }
			}	
			close(LIST);
			my @DEcheck=glob("$od/TFBS_Analysis/each_DEgeneRes/EN*");
                        my @DEplots=glob("$od/TFBS_Analysis/DEG_seqLogo/EN*.p*");
                        if((@DEcheck < 0) || (@DEplots < 0)){
                                die "There are something wrong may exits";
                        }
		}
	}
#########################	
	open(SH,">$od/work_sh/s3.enrich_anno.sh")||die $!;
	foreach my $deg(@degs){
		my $deg_num = `wc -l $deg`;
		if($deg_num > 1){
			my $base=basename $deg;
			my $vs=(split(/\./,$base))[0];
			my $odir="$od/$vs/Anno_enrichment"; `mkdir $odir`   unless(-d $odir);
		
			my $cmd = "export LD_LIBRARY_PATH=/share/nas2/genome/biosoft/gcc/5.4.0/lib64:/share/nas2/genome/biosoft/gcc/5.4.0/lib:\$LD_LIBRARY_PATH && $rscript $Bin/enrich/enrich.r --all $od/$vs/$vs.all --deg $od/$vs/$vs.DEG_final.xls --go $od/GO.info --kegg $od/KEGG.info --prefix $vs --od $odir/enrich ";
			$cmd .="&& perl $Bin/enrich/bin/plot_deg_KEGG.pl $odir/enrich/${vs}_KEGG_pathway_enrich.list $od/$vs/$vs.DEG_final.xls $odir/enrich ";
			$cmd .="&& perl $Bin/enrich/anno.pl -input $od/$vs/$vs.DEG_final.xls -od $odir/anno -anno $anno ";
			$cmd .="&& perl $Bin/enrich/anno/Kegg_map_web.pl -d $od/$vs/$vs.DEG_final.xls -k $vs -i $anno -o $odir ";
			if(exists $config{medical}){
				$cmd .="&& perl $Bin/enrich/anno/cosmic_TF_anno/medical_anno.pl -deg $od/$vs/$vs.DEG_final.xls -o $od/$vs -db $config{medical} ";
			}
			print SH "$cmd\n";
		}
	}
	close(SH);
	&qsub("$od/work_sh/s3.enrich_anno.sh");

	if(defined $ppi){
		open (SH,">$od/work_sh/s4.deg_ppi.sh")||die $!;
		my $cmd = "perl $Bin/Protein_to_protein_deg.pl -od $od/DEG_PPI -ppi $ppi -idir $od/";
		print SH "$cmd\n";
		close (SH);
		&qsub("$od/work_sh/s4.deg_ppi.sh");
	}

}

if(defined $target){
	&target_find_All($target,"$od/All_DEG/All.DEG_final.xls","$od/All_DEG/All.DEG_final_target.xls");
	&run_or_die("perl $Bin/enrich/anno/get_info_from_backgroud.pl $od/All_DEG/All.DEG_final_target.xls $inter_anno $od/All_DEG/All.DEG_final_anno.xls");

	my @tars=&target_find($target,\@degs,"Target");
	print join("\n",@tars),"\n";
	&anno(\@tars,"Anno_enrichment","$od/work_sh/s3.enrich_anno.sh");
	&qsub("$od/work_sh/s3.enrich_anno.sh");
}

if(defined $cis_target){
	&target_find_All($cis_target,"$od/All_DEG/All.DEG_final.xls","$od/All_DEG/All.DEG_final_Cis_target.xls");
        &run_or_die("perl $Bin/enrich/anno/get_info_from_backgroud.pl $od/All_DEG/All.DEG_final_Cis_target.xls $inter_anno $od/All_DEG/All.DEG_final_Cis_anno.xls");
	my @tars=&target_find($cis_target,\@degs,"Cis_Target");
	&anno(\@tars,"Cis_Anno_enrichment","$od/work_sh/s3.1.cis_enrich_anno.sh");
	&qsub("$od/work_sh/s3.1.cis_enrich_anno.sh");
}

if(defined $trans_target){
	&target_find_All($trans_target,"$od/All_DEG/All.DEG_final.xls","$od/All_DEG/All.DEG_final_Trans_target.xls");
        &run_or_die("perl $Bin/enrich/anno/get_info_from_backgroud.pl $od/All_DEG/All.DEG_final_Trans_target.xls $inter_anno $od/All_DEG/All.DEG_final_Trans_anno.xls");
	my @tars=&target_find($trans_target,\@degs,"Trans_Target");
	&anno(\@tars,"Trans_Anno_enrichment","$od/work_sh/s3.2.trans_enrich_anno.sh");
	&qsub("$od/work_sh/s3.2.trans_enrich_anno.sh");
}

&run_or_die("perl $Bin/kmeans/kmean_analysis.pl -fpkm $fpkm -od $od/kmeans");

##############Functional sub fucntion
sub getAllDEG{
        my($deg,$fpkm,$out)=@_;
        my @degs=@{$deg};
        my %all=();
        my @vs=();
        my $out_head;
        foreach my $d(@degs){
                my $base=basename $d;
                $base=(split(/\./,$base))[0];
                push @vs,$base;
                open(DEG,$d)||die $!;
                my $head=<DEG>;chomp($head);
                my @header=split(/\t/,$head);
                $out_head .= "\t$base.$header[-3]\t$base.$header[-2]\t$base.$header[-1]";
                while(<DEG>){
                        chomp;next if($_=~/^#/);
                        my @tmp=split(/\t/,$_);
                        $all{$tmp[0]}{$base}=join("\t",($tmp[-3],$tmp[-2],$tmp[-1]));
                }
                close(DEG);
        }
        open(FPKM,$fpkm)||die $!;
        open(OUT,">$out")||die $!;
        my $head=<FPKM>;chomp($head);
        print OUT "$head"."$out_head\n";
        while(<FPKM>){
                chomp;my @tmp=split(/\t/,$_);
                next if(!exists $all{$tmp[0]});
                print OUT join("\t",@tmp);
                foreach my $v(@vs){
                        $all{$tmp[0]}{$v}="--\t--\t--"  if(!exists $all{$tmp[0]}{$v});
                        print OUT "\t$all{$tmp[0]}{$v}";
                }
                print OUT "\n";
        }
        close(OUT);
        close(FPKM);

}

sub getvs{
        my $cfg=shift;
        my @vs=();
        open(CFG,$conf)||die $!;
        while(<CFG>){
                chomp;next if($_=~/^#|^$|^\s+/);
                my @tmp=split(/\s+/,$_);
                push @vs,$tmp[1]        if($tmp[0] eq "Com" || $tmp[0] eq "Sep");
        }
        close(CFG);
        return @vs;
}

sub anno{
	my ($degs,$key,$sh)=@_;
	open(SH,">$sh")||die $!;
	foreach my $deg(@{$degs}){
		my $base=basename $deg; my $vs=(split(/\./,$base))[0];
		my $dir=dirname $deg;
		`mkdir $dir/$key`	unless(-d "$dir/$key");
		print SH "$rscript $Bin/enrich/enrich.r --deg $deg --go $od/GO.info --kegg $od/KEGG.info --prefix $vs --od $dir/$key/enrich && perl $Bin/enrich/anno.pl -input $deg -od $dir/$key/anno -anno $anno && perl $Bin/enrich/anno/Kegg_map_web.pl -d $deg -k $base -i $anno -o $dir/$key -target 1\n";
	}
	close(SH);

}

sub target_find{
	my ($target,$deg,$type)=@_;#taget file; deg dir(contained multi *_vs_*);Target,Cis_Target;Trans_Target
	open(TARGET,$target)||die $!;
	my %tar=();
	while(<TARGET>){
		chomp;next if($_=~/^#/);
		my ($id,$gene)=split(/\t/,$_);
		$tar{$id}=$gene;
	}
	close(TARGET);
	my @degs=@{$deg};
	my @outs=();
	foreach my $deg(@degs){
		my $base=basename $deg;
		my $d=dirname $deg;
		my $vs=(split(/\./,$base))[0];
		open(DEG,$deg)||die $!;
		my $header=<DEG>;	$header=~s/^#//;
		push @outs,"$d/$vs.DEG_final.$type.xls";
		open(OUT,">$d/$vs.DEG_final.$type.xls")||die $!;
#		print "$d/$vs.DEG_final.$type.xls\n";
		print OUT "#Gene\t$header";
		while(<DEG>){
			chomp;my @tmp=split(/\t/,$_);
			if(exists $tar{$tmp[0]}){
				my @genes=split(/\;|,/,$tar{$tmp[0]});
				foreach my $g(@genes){
					print OUT "$g\t",join("\t",@tmp),"\n";
				}
			}
		}
		close(OUT);
		close(DEG);
	}

	return @outs;
}



sub target_find_All{
        my ($target,$deg,$out)=@_;#taget file; deg file;Target,Cis_Target;Trans_Target
        open(TARGET,$target)||die $!;
        my %tar=();
        while(<TARGET>){
                chomp;next if($_=~/^#/);
                my ($id,$gene)=split(/\t/,$_);
                $tar{$id}=$gene;
        }
        close(TARGET);
        open(DEG,$deg)||die $!;
        my $header=<DEG>;       $header=~s/^#//;
        open(OUT,">$out")||die $!;
        print OUT "#Gene\t$header";
        while(<DEG>){
        	chomp;my @tmp=split(/\t/,$_);
                if(exists $tar{$tmp[0]}){
                	my @genes=split(/\;|,/,$tar{$tmp[0]});
                	foreach my $g(@genes){
                		print OUT "$g\t",join("\t",@tmp),"\n";
                        }
                }
        }
        close(OUT);
        close(DEG);

}

sub readPPI{
	my $file=shift;
	open(PPI,$file)||die $!;
	while(<PPI>){
		chomp;next if($_=~/^#/);
		my @tmp=split(/\t/,$_);
	}
	close(PPI);

}
#################
sub readConfig{
	my $configFile=shift;
	my $d=Config::General->new(-ConfigFile => "$configFile");
	my %config=$d->getall;	
	return %config;
}
sub qsub(){
        my $shfile= shift;
	my $queue="medical.q";
	$queue=$config{Queue_type}	if(defined $config{Queue_type});
        my $cmd = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub $shfile --independent --queue $queue ";
        &run_or_die($cmd);              
        return ;
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
        return ;
}
sub show_log(){
        my ($txt) = @_ ;
        my $time = time();
        my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime($time);
        $wday = $yday = $isdst = 0;
        my $Time=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
        print "$Time:\t$txt\n" ;
        return ($time) ;
}

sub USAGE{
	my $usage=<<"USAGE";
Program: $0
Version: 1.0.0
Contact: wenyh\@biomarker.com.cn
Usage:
	
	Options:
	-h	Help

	-conf   file    forced          DEG parameter
	-count	file	forced		All_gene_counts.list
	-fpkm	file	forced		All_gene_fpkm.list
	-type	str	optional	gene_exp type, default FPKM 
	-anno	file	forced		Gene Anno Dir, Only Given this parameter enrichment analysis would be done.
	-gff	file	forced		newgene gff3 file for TFBS analysis.
	-target	file	forced		For lncRNA/miRNA/circRNA, the relation of them with gene should be given
					If exists multiple target gene files, seperated by comma
					File format(2 column):	lncRNA	targetgene 
	-cis	file	forced		lncRNA cis_target
	-trans	file	forced		lncRNA trans_target

	-ppi	file	forced		ppi extract, only perform when no cis/trans/target file exists!

        -od	path	forced		output path
	
	-human forced                   if "yes"，human edger norep dispersion is 0.4^2=0.16，if not edger norep dispersion is 0.1^2=0.01
Example:
	perl $0 -conf deg.cfg -count All_gene_counts.list -fpkm All_gene_fpkm.list -type FPKM -anno Unigene_Annotation/Result -ppi PPI.txt -od DEG_Analysis/gene
 
USAGE
	print $usage;
	exit;
}


