#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
#my ($fIn,$fInn,$fOut);
my ($idir,$list);

GetOptions(
				"help|?" =>\&USAGE,
				"idir:s"=>\$idir,
				"name2list:s"=>\$list,
				) or &USAGE;
&USAGE unless ($idir and $list);
# ------------------------------------------------------------------

$idir = &ABSOLUTE_DIR($idir);
my $odir = &ABSOLUTE_DIR($idir);
open (IN,$list);
my %hash;
while (<IN>) {
       chomp;
       my @lines = split /\s+/,$_;
       $hash{$lines[0]} = $lines[1];
       $hash{"gene:$lines[0]"}=$lines[1];
       $hash{"lncRNA:$lines[0]"}=$lines[1];
}
close (IN);
my (%File_L,%File_G,%File_M,%File_M1,%File_D,%File_S,%File_Z);
@{$File_L{1}}=glob "$idir/BMK_*_LncRNA/BMK_4_LncRNA_Target/*Cis_target_gene.xls";	#more than 2 column, need 2 
@{$File_L{2}}=glob "$idir/BMK_*_LncRNA/BMK_4_LncRNA_Target/*Trans_target_gene.xls";	#2 column,need 2
@{$File_L{6}}=glob "$idir/BMK_7_Combine/BMK_5_coexp/Sig.lncRNA_mRNA.xls";	#more than 2 column,need 2
@{$File_L{7}}=glob "$idir/BMK_*_LncRNA/BMK_5_DEG_Analysis/BMK_1_All_DEG/All.DEG_final_*_anno.xls";#more than 2 column, need 2
@{$File_L{8}}=glob "$idir/BMK_*_LncRNA/BMK_5_DEG_Analysis/BMK_1_All_DEG/All.DEG_final_*_target.xls";#more than 2 column, need 2
@{$File_L{9}}=glob "$idir/BMK_*_LncRNA/BMK_5_DEG_Analysis/*_vs_*/BMK_1_Statistics_Visualization/*.DEG_final.*_Target.xls";
@{$File_L{10}}=glob "$idir/BMK_*_LncRNA/BMK_5_DEG_Analysis/*_vs_*/*Anno_enrichment/BMK_1_Annotation/*.annotation.xls";


@{$File_M{3}}=glob "$idir/BMK_7_Combine/BMK_3_ceRNA*/*ceRNA_pair_adjust_p_Sig.xls";	#more than 2 column,need 2,:
@{$File_M{1}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/*RNA-miRNA-mRNA/*_vs_*.*RNA-miRNA-mRNA.ceRNA_pair_adjust_p_Sig_diff.xls";	#more than 2 column,need 2,:

@{$File_M1{6}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/circRNA/*.cm.Interaction.list";	#more than 2 column,need 1
@{$File_M1{7}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/*RNA-miRNA-mRNA/*_vs_*.*RNA-miRNA-mRNA.Interaction.list";	#more than 2 column,need 1
@{$File_M1{8}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/gene/*_vs_*.m*.Interaction.list";	#more than 2 column,need 1
@{$File_M1{9}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/lncRNA/*_vs_*.lm*Interaction.list";	#more than 2 column,need 1
@{$File_M1{12}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/lncRNA/*_vs_*.ls.Interaction.list";	#more than 2 column,need 1
@{$File_M1{10}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/sRNA/*_vs_*.sl.Interaction.list";	#more than 2 column,need 1
@{$File_M1{11}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/sRNA/*_vs_*.sm.Interaction.list";	#more than 2 column,need 1

@{$File_G{3}}=glob "$idir/BMK_7_Combine/BMK_3_ceRNA*/random/*.txt";	#2 column,need 1
@{$File_G{4}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/circRNA/*_vs_*.cm.Attribution.list";	#3 column,need 1
@{$File_G{5}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/*RNA-miRNA-mRNA/*_vs_*.*RNA-miRNA-mRNA.Attribution.list";	#3 column,need 1
@{$File_G{7}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/gene/*_vs_*.mc.Attribution.list";	#3 column,need 1
@{$File_G{8}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/gene/*_vs_*.ml*Attribution.list";	#3 column,need 1
@{$File_G{9}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/gene/*_vs_*.ms.Attribution.list";	#3 column,need 1
@{$File_G{10}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/lncRNA/*_vs_*.lm*Attribution.list";	#3 column,need 1
@{$File_G{11}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/lncRNA/*_vs_*.ls.Attribution.list";	#3 column,need 1
@{$File_G{12}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/sRNA/*_vs_*.sl.Attribution.list";	#3 column,need 1
@{$File_G{13}}=glob "$idir/BMK_7_Combine/BMK_4_cytoscape/sRNA/*_vs_*.sm.Attribution.list";	#3 column,need 1

@{$File_D{1}}=glob "$idir/BMK_3_mRNA/BMK_3_DEG_Analysis/BMK_PPI/*_vs_*.DEG.detail.xls";
@{$File_S{1}}=glob "$idir/BMK_3_mRNA/BMK_3_DEG_Analysis/BMK_PPI/*_vs_*.ppi.cytoscapeInput.sif";

#gene;ceRNA1 up/down lncRNA:ceRNA2 up/down other_info1 ...
@{$File_Z{1}}=glob "$idir/BMK_7_Combine/BMK_3_ceRNA*/*/sub_*_network.xls";
@{$File_Z{3}}=glob "$idir/BMK_7_Combine/BMK_*_coexp/Diff.Sig.lncRNA_mRNA.xls";
#target_gene_id type gene_id/lncRNA_id
@{$File_Z{2}}=glob "$idir/BMK_7_Combine/BMK_3_ceRNA*/*/*Key_RNA_gene.xls";

foreach my $z (keys %File_Z){
        if(@{$File_Z{$z}}>=1){
        foreach my $zz (@{$File_Z{$z}}){
                open(IN,$zz)||die $!;
                open(OUT,">${zz}_Symbol.xls");
                while(<IN>){
                        chomp;
                        my @tmp = split /\s+/,$_;
                        if(@tmp>3){
                                if(/^#|^ceRNA|^RNA/){
                                        my ($ce1,$up,$ce2,$down,$info)=split /\s+/,$_,5;
                                        print OUT "$ce1\tSymbol\t$up\t$ce2\tSymbol\t$down\t$info\n";
                                        next;
                                }
                                my ($ceRNA1,$regulated1,$ceRNA2,$regulated2,$detail)=split /\s+/,$_,5;
                                if(exists $hash{$ceRNA1}){
                                        print OUT "$ceRNA1\t$hash{$ceRNA1}\t$regulated1\t";
                                }else{
					(my $s1 = $ceRNA1) =~ s/^gene:|^lncRNA:|^circRNA://;
                                        print OUT "$ceRNA1\t$s1\t$regulated1\t";
                                }
                                if(exists $hash{$ceRNA2}){
                                        print OUT "$ceRNA2\t$hash{$ceRNA2}\t$regulated2\t";
                                }else{
					(my $s2 = $ceRNA2) =~ s/^gene:|^lncRNA:|^circRNA://;
                                        print OUT "$ceRNA2\t$s2\t$regulated2\t";
                                }
                                print OUT "$detail\n";
                        }else{
                                if(/^#|^ceRNA|^RNA/){
                                        my ($Id1,$type,$Id2)=split /\s+/,$_,3;
                                        print OUT "$Id1\tSymbol\t$type\t$Id2\tSymbol\n";
                                        next;
                                }
                                my ($gene1,$t,$gene2)=split /\s+/,$_,3;
                                if(exists $hash{$gene1}){
                                        print OUT "$gene1\t$hash{$gene1}\t$t\t";
                                }else{
                                        print OUT "$gene1\t$gene1\t$t\t";
                                }
                                if(exists $hash{$gene2}){
                                        print OUT "$gene2\t$hash{$gene2}\n";
                                }else{
                                        print OUT "$gene2\t$gene2\n";
                                }
                        }
                }
                close(IN);
                close(OUT);
                `mv ${zz}_Symbol.xls $zz`;
        }
        }
}

foreach my $key (keys %File_G){
    if(@{$File_G{$key}}>=1){
     foreach my $file (@{$File_G{$key}}){
        open (IN,$file)||die $!;
        open (OUT,">${file}_new.xls")||die $!;
        while(<IN>){
            chomp;
	    next if(/^#/);        
	    my ($gene,$info)=split /\s+/,$_,2;
            if(exists $hash{$gene}){print OUT "$gene\t$hash{$gene}\t$info\n";}
	    else {print OUT "$gene\t$gene\t$info\n";}
        }
        close IN;
        close OUT;
	`mv ${file}_new.xls ${file}`;
	 }
     }
}

foreach my $k2 (keys %File_L){
	if(@{$File_L{$k2}}>0){
        foreach my $f2 (@{$File_L{$k2}}){
                open (IN,$f2)||die $!;
                open (OUT,">${f2}_new.xls")||die $!;
                while(<IN>){
                        chomp;
			my @tmp = split /\t/,$_;
			if(@tmp>2){
	                        if(/^#/){my ($id1,$id2,$info)=split /\t/,$_,3;print OUT "$id1\tSymbol\t$id2\tSymbol\t$info\n";next;}
        	                my ($lnc,$gene,$info)=split /\t/,$_,3;
                	        if(exists $hash{$lnc}){
                        	        print OUT "$lnc\t$hash{$lnc}\t";
	                        }else{
        	                        print OUT "$lnc\t$lnc\t";
                	        }
                        	my @tmp = split /\;|,/,$gene;
	                        my $len = @tmp;
        	                my $new;
                	        if($len==1){
                        	        if(exists $hash{$tmp[0]}){
                                	        print OUT "$tmp[0]\t$hash{$tmp[0]}\t";
	                                }else{
        	                                print OUT "$tmp[0]\t$tmp[0]\t";
                	                }
                        	}else{
                                	if(exists $hash{$tmp[0]}){
                                        	$new = $hash{$tmp[0]};
	                                }else{
        	                                $new = $tmp[0];
                	                }
                        	        for (my $i=1;$i<@tmp;$i++){
                                	        if(exists $hash{$tmp[$i]}){
                                        	        $new .= ";".$hash{$tmp[$i]};
	                                        }else{
        	                                        $new .= ";".$tmp[$i];
                	                        }
                        	        }
                                	print OUT "$gene\t$new\t";
	                        }
        	                print OUT "$info\n";
			}else{
				if(/^#/){my ($id1,$id2)=split /\t/,$_,2;print OUT "$id1\tSymbol\t$id2\tSymbol\n";next;}
                                my ($lnc,$gene)=split /\t/,$_,2;
                                if(exists $hash{$lnc}){
                                        print OUT "$lnc\t$hash{$lnc}\t";
                                }else{
                                        print OUT "$lnc\t$lnc\t";
                                }
                                my @tmp = split /\;|,/,$gene;
                                my $new;
                                if(@tmp==1){
                                        if(exists $hash{$tmp[0]}){
                                                print OUT "$tmp[0]\t$hash{$tmp[0]}\t";
                                        }else{
                                                print OUT "$tmp[0]\t$tmp[0]\t";
                                        }
                                }else{
                                        if(exists $hash{$tmp[0]}){
                                                $new = $hash{$tmp[0]};
                                        }else{
                                                $new = $tmp[0];
                                        }
                                        for (my $i=1;$i<@tmp;$i++){
                                                if(exists $hash{$tmp[$i]}){
                                                        $new .= ";".$hash{$tmp[$i]};
                                                }else{
                                                        $new .= ";".$tmp[$i];
                                                }
                                        }
                                        print OUT "$gene\t$new\n";
                                }
			}
                }
                close(IN);
                close(OUT);
                `mv ${f2}_new.xls ${f2}`;
        }
   }
}

foreach my $k11 (keys %File_M1){
	if(@{$File_M1{$k11}}>0){
         foreach my $f11 (@{$File_M1{$k11}}){
                open (IN,$f11)||die $!;
                open (OUT,">${f11}_new.xls")||die $!;
                while(<IN>){
                        chomp;
                        my @tmp = split /\t/,$_;
                        if(@tmp>2){
                                if(/^#|^RNA/){my ($id1,$id2,$info)=split /\t/,$_,3;print OUT "$id1\tSymbol\t$id2\tSymbol\t$info\n";next;}
                                my ($g1,$g2,$info)=split /\t/,$_,3;
                                if(exists $hash{$g1}){
                                        print OUT "$g1\t$hash{$g1}\t";
                                }else{
                                        print OUT "$g1\t$g1\t";
                                }
                                if(exists $hash{$g2}){
                                        print OUT "$g2\t$hash{$g2}\t";
                                }else{
                                        print OUT "$g2\t$g2\t";
                                }
                                print OUT "$info\n";
                        }else{
                                if(/^#|^RNA/){my ($id1,$id2)=split /\t/,$_,2;print OUT "$id1\tSymbol\t$id2\tSymbol\n";next;}
                                my ($g1,$g2)=split /\t/,$_,2;
                                if(exists $hash{$g1}){
                                        print OUT "$g1\t$hash{$g1}\t";
                                }else{
                                        print OUT "$g1\t$g1\t";
                                }
                                if(exists $hash{$g2}){
                                        print OUT "$g2\t$hash{$g2}\n";
                                }else{
                                        print OUT "$g2\t$g2\n";
                                }
                        }
                }
                close (IN);
                close (OUT);
                `mv ${f11}_new.xls ${f11}`;
        }
   }
}

foreach my $k3 (keys %File_D){
	if(@{$File_D{$k3}}){
          foreach my $f3 (@{$File_D{$k3}}){
		open (IN,$f3)||die $!;
		open (OUT,">${f3}_Symbol.xls") ||die $!;
		while (<IN>){
			chomp;
			if(/^#/){my($i1,$p,$i2,$info)=split /\t/,$_,4;print OUT "$i1\tSymbol\t$p\t$i2\tSymbol\t$info\n";next;}
			my ($g1,$pp,$g2,$info) = split /\t/,$_,4;
			my ($num1,$id1)=split /\./,$g1;
			my ($num2,$id2)=split /\./,$g2;
			if(exists $hash{$g1}){
				print OUT "$g1\t$hash{$g1}\t";
			}else{
				print OUT "$g1\t$g1\t";
			}
			print OUT "$pp\t";
			if(exists $hash{$g2}){
				print OUT "$g2\t$hash{$g2}\t";
			}else{
				print OUT "$g2\t$g2\t";
			}
			print OUT "$info\n";
		}
		close(IN);
		close(OUT);
		`mv ${f3}_Symbol.xls ${f3}`;
	}
    }
}

foreach my $k4 (keys %File_S){
	if(@{$File_S{$k4}}>0){
         foreach my $f4 (@{$File_S{$k4}}){
                open (IN,$f4)||die $!;
                open (OUT,">${f4}_Symbol.xls")||die $!;
                while (<IN>){
                        chomp;
                        if(/^#/){my($i1,$p,$i2)=split /\t/,$_,3;print OUT "$i1\tSymbol\t$p\t$i2\tSymbol\n";next;}
                        my ($g1,$pp,$g2) = split /\t/,$_,3;
                        my ($num1,$id1)=split /\./,$g1;
                        my ($num2,$id2)=split /\./,$g2;
                        if(exists $hash{$g1}){
                                print OUT "$hash{$g1}\t";
                        }else{
                                print OUT "$g1\t";
                        }
                        print OUT "$pp\t";
                        if(exists $hash{$g2}){
                                print OUT "$hash{$g2}\n";
                        }else{
                                print OUT "$g2\n";
                        }
                }
		close(IN);
		close(OUT);
		`mv ${f4}_Symbol.xls ${f4}`;
        }
    }
}

foreach my $k1 (keys %File_M){
	if(@{$File_M{$k1}}>0){
         foreach my $f1 (@{$File_M{$k1}}){
                open (IN,$f1)||die $!;
                open (OUT,">${f1}_new.xls")||die $!;
                while(<IN>){
                        chomp;
                        my @tmp = split /\t/,$_;
                        if(@tmp>3){
                                if(/^#|^RNA/){my ($id1,$id2,$info)=split /\t/,$_,3;print OUT "$id1\tSymbol\t$id2\tSymbol\t$info\n";next;}
                                my ($g1,$g2,$info)=split /\t/,$_,3;
                                my ($type1,$gene1)=split /:/,$g1,2;
                                my ($type2,$gene2)=split /:/,$g2,2;
                                if(exists $hash{$g1}){
                                        print OUT "$g1\t$hash{$g1}\t";
                                }else{
                                        print OUT "$g1\t$gene1\t";
                                }
                                if(exists $hash{$g2}){
                                        print OUT "$g2\t$hash{$g2}\t";
                                }else{
                                        print OUT "$g2\t$gene2\t";
                                }
                                print OUT "$info\n";
                        }else{
                                if(/^#|^RNA/){my ($id1,$id2)=split /\t/,$_,2;print OUT "$id1\tSymbol\t$id2\tSymbol\n";next;}
                                my ($g1,$g2)=split /\t/,$_,2;
                                my ($type1,$gene1)=split /:/,$g1,2;
                                my ($type2,$gene2)=split /:/,$g2,2;
                                if(exists $hash{$g1}){
                                        print OUT "$g1\t$hash{$g1}\t";
                                }else{
                                        print OUT "$g1\t$gene1\t";
                                }
                                if(exists $hash{$g2}){
                                        print OUT "$g2\t$hash{$g2}\n";
                                }else{
                                        print OUT "$g2\t$gene2\n";
                                }
                        }
                }
                close (IN);
                close (OUT);
                `mv ${f1}_new.xls ${f1}`;
        }
    }
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

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

sub max{#&max(lists or arry);
	#求列表中的最大值
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

################################################################################################################

sub min{#&min(lists or arry);
	#求列表中的最小值
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

################################################################################################################

sub revcom(){#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
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
     Contact:	Simon Young <yangxh\@biomarker.com.cn> 
Program Date:	2012.07.02
      Modify:	
 Description:	This program is used to ......
       Usage:
		Options:
		-idir <dir>	BMK_Result
		-name2list <file>	geneid2name.list
		-h		help

USAGE
	print $usage;
exit;
}
