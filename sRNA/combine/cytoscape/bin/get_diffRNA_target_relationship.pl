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
my ($idir,$odir,$data_cfg,$detail_cfg,$l2mc,$l2mt,$mi2l,$mi2c,$mi2g,$chostgene,$l2m);

GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$odir,
				"in:s"=>\$idir,
				"data_cfg:s"=>\$data_cfg,
				"detail_cfg:s"=>\$detail_cfg,
				"l2mc:s"=>\$l2mc,
				"l2mt:s"=>\$l2mt,
				"l2m:s"=>\$l2m,
				"mi2l:s"=>\$mi2l,
				"mi2c:s"=>\$mi2c,
				"mi2g:s"=>\$mi2g,
				"chostgene:s"=>\$chostgene,
				) or &USAGE;
&USAGE unless ($idir and $odir and $data_cfg and $detail_cfg);
# ------------------------------------------------------------------

$data_cfg = &ABSOLUTE_DIR($data_cfg);
$detail_cfg = &ABSOLUTE_DIR($detail_cfg);
$idir = &ABSOLUTE_DIR($idir);
`mkdir $odir` unless (-d "$odir");
$odir = &ABSOLUTE_DIR($odir);

$l2m ||= "$idir/../Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/novel_target_gene.xls" if (-e "$idir/../Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/novel_target_gene.xls");
$l2mc ||= "$idir/../Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/Cis_target_gene.xls" if(-e "$idir/../Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/Cis_target_gene.xls");
$l2mt ||= "$idir/../Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/Trans/Trans_target_gene.xls" if(-e "$idir/../Basic_Analysis/LncRNA_Analysis/Lnc_target_predict/Trans/Trans_target_gene.xls");
$mi2l ||= "$idir/../miRNA_Target/lncRNA/mir2target.list";
$mi2c ||= "$idir/../miRNA_Target/circRNA/mir2target.list";
$mi2g ||= "$idir/../miRNA_Target/gene/mir2target.list";
$chostgene = "$idir/../Basic_Analysis/circRNA_analysis/new_name/Circ_source_gene.xls";

my %sample = &relate($data_cfg);
my @diff;my %vs;
open (CFG,"$detail_cfg") or die $!;
print ("read cfgi and diff group is:");
while (<CFG>){
	chomp;
	next if(/^#|^$/);
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

my (%l2m,%m2l,%mi2l,%l2mi,%mi2c,%c2mi,%mi2g,%g2mi,%c2g,%g2c);

###read lncRNA2mRNA file and then get %l2m and %m2l
print ("1:L2M...... %l2m,%m2l\n");
if (-e "$l2m"){
	open (L2M,"$l2m");
	while(<L2M>){
		chomp;
		next if(/^#|^$/);
		my @A = split /\t/,$_;
		my @B = split /\;|,/,$A[1];
		my $nm = join (";",@B);
		$l2m{$A[0]}=$nm;
		for (my $i=0;$i<@B;$i++){
			if(exists $m2l{$B[$i]}){
				$m2l{$B[$i]} .= ";".$A[0];
			}else{
				$m2l{$B[$i]} = $A[0];
			}
		}
	}
	close (L2M);
}else {
	open (L2MC,"$l2mc");
	while(<L2MC>){
		chomp;
		next if(/^#|^\s$/);
		my @A = split /\t/,$_;
		my @B = split /\;|,/,$A[1];
		my $nm = join (";",@B);
		$l2m{$A[0]}{Cis}=$nm;
		for (my $i=0;$i<@B;$i++){
			if(exists $m2l{$B[$i]}{Cis}){
				$m2l{$B[$i]}{Cis} .= ";".$A[0];
			}else{
				$m2l{$B[$i]}{Cis} = $A[0];
			}
		}
	}
	if(-e $l2mt){
		open (L2MT,"$l2mt");
		while (<L2MT>){
			next if(/^#|^\s$/);
			my @A = split /\t/,$_;
	                my @B = split /\;|,/,$A[1];
                	my $nm = join (";",@B);
        	        $l2m{$A[0]}{Trans}=$nm;
	                for (my $i=0;$i<@B;$i++){
                	        if(exists $m2l{$B[$i]}{Trans}){
        	                        $m2l{$B[$i]}{Trans} .= ";".$A[0];
	                        }else{
	                               	$m2l{$B[$i]}{Trans} = $A[0];
                        	}
                	}
		}
	}
}


###read miRNA2lncRNA file and then get %mi2l and %l2mi
print "2:MI2L...... %mi2l,%l2mi\n";
open (MI2L,"$mi2l");
while(<MI2L>){
	chomp;
	next if(/^#|^$/);
	my @A = split /\t/,$_;
	my @B = split /\;|,/,$A[1];
	my $nm = join (";",@B);
	$mi2l{$A[0]}=$nm;
	for (my $i=0;$i<@B;$i++){
		if(exists $l2mi{$B[$i]}){
			$l2mi{$B[$i]} .= ";".$A[0];
		}else{
			$l2mi{$B[$i]} = $A[0];
		}
	}
}
close (MI2L);

###read miRNA2circRNA file
print ("3:MI2C...... %mi2c,%c2mi\n");
open(MI2C,"$mi2c");
while(<MI2C>){
	chomp;
	next if(/^#|^$/);
	my @A = split /\t/,$_;
	my @B = split /\;|,/,$A[1];
	my $nm = join (";",@B);
	$mi2c{$A[0]}=$nm;
	for (my $i=0;$i<@B;$i++){
		if(exists $c2mi{$B[$i]}){
			$c2mi{$B[$i]} .= ";".$A[0];
		}else{
			$c2mi{$B[$i]} = $A[0];
		}
	}
}
close (MI2C);

###read miRNA2gene file
print ("4:MI2G...... %mi2g,%g2mi\n");
open (MI2G,"$mi2g");
while(<MI2G>){
	chomp;
	next if(/^#|^$/);
	my @A = split /\t/,$_;
	my @B = split /\;|,/,$A[1];
	my $nm = join (";",@B);
	$mi2g{$A[0]}=$nm;
	for (my $i=0;$i<@B;$i++){
		if(exists $g2mi{$B[$i]}){
			$g2mi{$B[$i]} .= ";".$A[0];
		}else{
			$g2mi{$B[$i]} =$A[0];
		}
	}
}
close (MI2G);

###read circRNA hostgene file
print ("5:C2G...... %c2g,%g2c\n");
open (C2G,"$chostgene");
while (<C2G>){
	chomp;
	next if(/^#|^$/);
	my @A = split /\t/,$_;
	my @B = split /\;|,/,$A[1];
	my $nm = join (";",@B);
	$c2g{$A[0]}=$nm;
	for (my $i=0;$i<@B;$i++){
		if(exists $g2c{$B[$i]}){
			$g2c{$B[$i]} .= ";".$A[0];
		}else{
			$g2c{$B[$i]} = $A[0];
		}
	}
}
close(C2G);

print ("read target file progress finish!!!\n\n");

print "get diff result file and creat hash\n\n";
print "#group\tRNA_type\tdiff_files\n";
my (%h,%hash);
foreach my $vs(@diff){
	&switch($vs);
	foreach my $k ("sRNA","lncRNA","circRNA","gene"){
		$h{$vs}{$k} = "$idir/$k/$vs{$vs}{$k}/$vs{$vs}{$k}.DEG_final.xls";
		print "$vs\t$k\t$h{$vs}{$k}\n";
	}
}

print "\n%h===>\$hash{\$deg}{\$trans}{\$id}=\$regulate\n\n";
foreach my $deg (keys %h){
	foreach my $trans(keys %{$h{$deg}}){
		open (IN,"$h{$deg}{$trans}");
		while (<IN>){
			chomp;
			next if(/^#/);
			my @a = split /\t/,$_;
			my $info = $a[0]."(".$a[-1].")";
			$hash{$deg}{$trans}{$a[0]}=$info;
		}
		close(IN);
	}
}
###################lncRNA core#############################
print "core-lncRNA\n";
`mkdir $odir/lncRNA` unless (-d "$odir/lncRNA");
my $lncdir = "$odir/lncRNA";
if(-e $l2m){
	foreach my $test (keys %h){
		open (IN,"$h{$test}{lncRNA}");
		open (OUT,">$lncdir/$test.DEG.xls");
		my $lminter = "$lncdir/$test.lm.Interaction.list";
		my $lmatt = "$lncdir/$test.lm.Attribution.list";
		my $lsinter = "$lncdir/$test.ls.Interaction.list";
		my $lsatt = "$lncdir/$test.ls.Attribution.list";
		my (%lmCyto,%lsCyto);
		while(<IN>){
			chomp;
			if(/^#/){print OUT "$_\tdiff_target_gene\tdiff_target_miRNA\n";next;}
			my @A = split /\t/,$_;
			print OUT "$_\t";
			if(exists $l2m{$A[0]}){
				my @lm = split /\;/,$l2m{$A[0]};
				my $ginfo;my $flag =0;
				for (my $i=0;$i<@lm;$i++){
					if(exists $hash{$test}{gene}{$lm[$i]}){
						$flag++;
						if($flag==1){
							$ginfo = $hash{$test}{gene}{$lm[$i]};
						}else{
							$ginfo .= ";".$hash{$test}{gene}{$lm[$i]};
						}
					}	
				}
			
				if($flag >0){
					print OUT "$ginfo\t";
					$lmCyto{$A[0]}{$A[-1]}=$ginfo;
				}else{
					print OUT "--\t";
				}
			}else{
                	        print OUT "--\t";
                	}

			if(exists $l2mi{$A[0]}){
				my @lmi = split /\;/,$l2mi{$A[0]};
				my $sinfo;my $flag =0;
				for (my $i=0;$i<@lmi;$i++){
					if(exists $hash{$test}{sRNA}{$lmi[$i]}){
						$flag++;
						if($flag ==1){
							$sinfo = $hash{$test}{sRNA}{$lmi[$i]};
						}else{
							$sinfo .= ";".$hash{$test}{sRNA}{$lmi[$i]};
						}
					}
				}
				if($flag >0){
					print OUT "$sinfo\n";
					$lsCyto{$A[0]}{$A[-1]}=$sinfo;
				}else{
					print OUT "--\n";
				}
			}else{
				print OUT "--\n";
			}
		}
		&getCyto($lminter,$lmatt,"lm",\%lmCyto);
		&getCyto($lsinter,$lsatt,"ls",\%lsCyto);
		close (IN);
		close (OUT);	
		}
}else{
	foreach my $test (keys %h){
		open (IN,"$h{$test}{lncRNA}");
		open (OUT,">$lncdir/$test.DEG.xls");
		my $lmcinter = "$lncdir/$test.lm.Cis.Interaction.list";
		my $lmcatt = "$lncdir/$test.lm.Cis.Attribution.list";
        	my $lsinter = "$lncdir/$test.ls.Interaction.list";
        	my $lsatt = "$lncdir/$test.ls.Attribution.list";
		my $lmtinter = "$lncdir/$test.lm.Trans.Interaction.list";
        	my $lmtatt = "$lncdir/$test.lm.Trans.Attribution.list";
	        my (%lmcCyto,%lmtCyto,%lsCyto);
	        while(<IN>){
         	       chomp;
			if(-e $l2mt){
	                	if(/^#/){print OUT "$_\tdiff_Cis_target_gene\tdiff_Trans_target_gene\tdiff_target_miRNA\n";next;}
			}else{
				if(/^#/){print OUT "$_\tdiff_Cis_target_gene\tdiff_target_miRNA\n";next;}
			}
		    	my @A = split /\t/,$_;
		        print OUT "$_\t";
        	        if(exists $l2m{$A[0]}{Cis}){
                	        my @lm = split /\;/,$l2m{$A[0]}{Cis};
                        	my $ginfo;my $flag =0;
	                        for (my $i=0;$i<@lm;$i++){
        	                        if(exists $hash{$test}{gene}{$lm[$i]}){
                	                        $flag++;
                        	                if($flag==1){
                                	                $ginfo = $hash{$test}{gene}{$lm[$i]};
	                                        }else{
        	                                        $ginfo .= ";".$hash{$test}{gene}{$lm[$i]};
                	                        }
                        	        }
	                        }

        	                if($flag >0){
                	                print OUT "$ginfo\t";
                        	        $lmcCyto{$A[0]}{$A[-1]}=$ginfo;
	                        }else{
        	                        print OUT "--\t";
                	        }
	                }else{
        	                print OUT "--\t";
                	}
		
			if(-e $l2mt){
				if(exists $l2m{$A[0]}{Trans}){
                	        	my @lm = split /\;/,$l2m{$A[0]}{Trans};
	                        	my $ginfo;my $flag =0;
		                        for (my $i=0;$i<@lm;$i++){
        		                        if(exists $hash{$test}{gene}{$lm[$i]}){
                		                        $flag++;
                        		                if($flag==1){
                                		                $ginfo = $hash{$test}{gene}{$lm[$i]};
	                                        	}else{
        	                                        	$ginfo .= ";".$hash{$test}{gene}{$lm[$i]};
	                	                        }
        	                	        }
	        	                }

	        	                if($flag >0){
        	        	                print OUT "$ginfo\t";
                	        	        $lmtCyto{$A[0]}{$A[-1]}=$ginfo;
	                	        }else{
        	                	        print OUT "--\t";
                	        	}
		
	                	}else{
        	               		print OUT "--\t";
	                	}
			}

                if(exists $l2mi{$A[0]}){
                        my @lmi = split /\;/,$l2mi{$A[0]};
                        my $sinfo;my $flag =0;
                        for (my $i=0;$i<@lmi;$i++){
                                if(exists $hash{$test}{sRNA}{$lmi[$i]}){
                                        $flag++;
                                        if($flag ==1){
                                                $sinfo = $hash{$test}{sRNA}{$lmi[$i]};
                                        }else{
                                                $sinfo .= ";".$hash{$test}{sRNA}{$lmi[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$sinfo\n";
                                $lsCyto{$A[0]}{$A[-1]}=$sinfo;
                        }else{
                                print OUT "--\n";
                        }
                }else{
                        print OUT "--\n";
                }
        }
        &getCyto($lmcinter,$lmcatt,"lm",\%lmcCyto);
	&getCyto($lmtinter,$lmtatt,"lm",\%lmtCyto);
        &getCyto($lsinter,$lsatt,"ls",\%lsCyto);
        close (IN);
        close (OUT);
	}
}
print "core-circRNA\n";
`mkdir $odir/circRNA` unless (-d "$odir/circRNA");
my $circdir = "$odir/circRNA";
foreach my $circ (keys %h){
        open (IN,"$h{$circ}{circRNA}");
        open (OUT,">$circdir/$circ.DEG.xls");
	my $cginter = "$circdir/$circ.cm.Interaction.list";
        my $cgatt = "$circdir/$circ.cm.Attribution.list";
        my $csinter = "$circdir/$circ.cs.Interaction.list";
        my $csatt = "$circdir/$circ.cs.Attribution.list";
        my (%cgCyto,%csCyto);
        while(<IN>){
                chomp;
                if(/^#/){print OUT "$_\tdiff_host_gene\tdiff_target_miRNA\n";next;}
                my @A = split /\t/,$_;
                print OUT "$_\t";
                if(exists $c2g{$A[0]}){
                        my @cg = split /\;/,$c2g{$A[0]};
                        my $ginfo;my $flag =0;
                        for (my $i=0;$i<@cg;$i++){
                                if(exists $hash{$circ}{gene}{$cg[$i]}){
                                        $flag++;
                                        if($flag==1){
                                                $ginfo = $hash{$circ}{gene}{$cg[$i]};
                                        }else{
                                                $ginfo .= ";".$hash{$circ}{gene}{$cg[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$ginfo\t";
				$cgCyto{$A[0]}{$A[-1]}=$ginfo;
                        }else{
                                print OUT "--\t";
                        }
                }else{
                        print OUT "--\t";
                }

                if(exists $c2mi{$A[0]}){
                        my @cmi = split /;/,$c2mi{$A[0]};
                        my $sinfo;my $flag =0;
                        for (my $i=0;$i<@cmi;$i++){
                                if(exists $hash{$circ}{sRNA}{$cmi[$i]}){
                                        $flag++;
                                        if($flag ==1){
                                                $sinfo = $hash{$circ}{sRNA}{$cmi[$i]};
                                        }else{
                                                $sinfo .= ";".$hash{$circ}{sRNA}{$cmi[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$sinfo\n";
				$csCyto{$A[0]}{$A[-1]}=$sinfo;
                        }else{
                                print OUT "--\n";
                        }
                }else{
                        print OUT "--\n";
                }
        }
        close (IN);
	close (OUT);
	&getCyto($cginter,$cgatt,"cm",\%cgCyto);
        &getCyto($csinter,$csatt,"cs",\%csCyto);
}

print "core-miRNA\n";
`mkdir $odir/sRNA` unless (-d "$odir/sRNA");
my $sdir = "$odir/sRNA";
foreach my $mi (keys %h){
        open (IN,"$h{$mi}{sRNA}");
        open (OUT,">$sdir/$mi.DEG.xls");
	my $sginter = "$sdir/$mi.sm.Interaction.list";
        my $sgatt = "$sdir/$mi.sm.Attribution.list";
        my $slinter = "$sdir/$mi.sl.Interaction.list";
        my $slatt = "$sdir/$mi.sl.Attribution.list";
	my $scinter = "$sdir/$mi.sc.Interaction.list";
        my $scatt = "$sdir/$mi.sc.Attribution.list";
        my (%sgCyto,%slCyto,%scCyto);
        while(<IN>){
                chomp;
                if(/^#/){print OUT "$_\tdiff_target_gene\tdiff_target_lncRNA\tdiff_target_circRNA\n";next;}
                my @A = split /\t/,$_;
                print OUT "$_\t";
                if(exists $mi2g{$A[0]}){
                        my @sg = split /\;/,$mi2g{$A[0]};
                        my $ginfo;my $flag =0;
                        for (my $i=0;$i<@sg;$i++){
                                if(exists $hash{$mi}{gene}{$sg[$i]}){
                                        $flag++;
                                        if($flag==1){
                                                $ginfo = $hash{$mi}{gene}{$sg[$i]};
                                        }else{
                                                $ginfo .= ";".$hash{$mi}{gene}{$sg[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$ginfo\t";
				$sgCyto{$A[0]}{$A[-1]}=$ginfo;
                        }else{
                                print OUT "--\t";
                        }
                }else{
                        print OUT "--\t";
                }

                if(exists $mi2l{$A[0]}){
                        my @sl = split /;/,$mi2l{$A[0]};
                        my $linfo;my $flag =0;
                        for (my $i=0;$i<@sl;$i++){
                                if(exists $hash{$mi}{lncRNA}{$sl[$i]}){
                                        $flag++;
                                        if($flag ==1){
                                                $linfo = $hash{$mi}{lncRNA}{$sl[$i]};
                                        }else{
                                                $linfo .= ";".$hash{$mi}{lncRNA}{$sl[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$linfo\t";
				$slCyto{$A[0]}{$A[-1]}=$linfo;
                        }else{
                                print OUT "--\t";
                        }
                }else{
                        print OUT "--\t";
                }
		
		if(exists $mi2c{$A[0]}){
                        my @cl = split /;/,$mi2c{$A[0]};
                        my $cinfo;my $flag =0;
                        for (my $i=0;$i<@cl;$i++){
                                if(exists $hash{$mi}{circRNA}{$cl[$i]}){
                                        $flag++;
                                        if($flag ==1){
                                                $cinfo = $hash{$mi}{circRNA}{$cl[$i]};
                                        }else{
                                                $cinfo .= ";".$hash{$mi}{circRNA}{$cl[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$cinfo\n";
				$scCyto{$A[0]}{$A[-1]}=$cinfo;
                        }else{
                                print OUT "--\n";
                        }
                }else{
                        print OUT "--\n";
                }
        }
        close (IN);
	close (OUT);
	&getCyto($sginter,$sgatt,"sm",\%sgCyto);
	&getCyto($slinter,$slatt,"sl",\%slCyto);
	&getCyto($scinter,$scatt,"sc",\%scCyto);
}

print "core-gene\n";
`mkdir $odir/gene` unless (-d "$odir/gene");
my $gdir = "$odir/gene";
if(-e "$l2m"){
foreach my $g (keys %h){
        open (IN,"$h{$g}{gene}");
        open (OUT,">$gdir/$g.DEG.xls");
	my $gsinter = "$gdir/$g.ms.Interaction.list";
        my $gsatt = "$gdir/$g.ms.Attribution.list";
        my $glinter = "$gdir/$g.ml.Interaction.list";
        my $glatt = "$gdir/$g.ml.Attribution.list";
        my $gcinter = "$gdir/$g.mc.Interaction.list";
        my $gcatt = "$gdir/$g.mc.Attribution.list";
        my (%gsCyto,%glCyto,%gcCyto);
        while(<IN>){
                chomp;
                if(/^#/){print OUT "$_\tdiff_target_miRNA\tdiff_target_lncRNA\tdiff_target_circRNA\n";next;}
                my @A = split /\t/,$_;
                print OUT "$_\t";
                if(exists $g2mi{$A[0]}){
                        my @gs = split /\;/,$g2mi{$A[0]};
                        my $sinfo;my $flag =0;
                        for (my $i=0;$i<@gs;$i++){
                                if(exists $hash{$g}{sRNA}{$gs[$i]}){
                                        $flag++;
                                        if($flag==1){
                                                $sinfo = $hash{$g}{sRNA}{$gs[$i]};
                                        }else{
                                                $sinfo .= ";".$hash{$g}{sRNA}{$gs[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$sinfo\t";
				$gsCyto{$A[0]}{$A[-1]}=$sinfo;
                        }else{
                                print OUT "--\t";
                        }
                }else{
                        print OUT "--\t";
                }

                if(exists $m2l{$A[0]}){
                        my @gl = split /;/,$m2l{$A[0]};
                        my $linfo;my $flag =0;
                        for (my $i=0;$i<@gl;$i++){
                                if(exists $hash{$g}{lncRNA}{$gl[$i]}){
                                        $flag++;
                                        if($flag ==1){
                                                $linfo = $hash{$g}{lncRNA}{$gl[$i]};
                                        }else{
                                                $linfo .= ";".$hash{$g}{lncRNA}{$gl[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$linfo\t";
				$glCyto{$A[0]}{$A[-1]}=$linfo;
                        }else{
                                print OUT "--\t";
                        }
                }else{
                        print OUT "--\t";
                }

                if(exists $g2c{$A[0]}){
                        my @gc = split /;/,$g2c{$A[0]};
                        my $cinfo;my $flag =0;
                        for (my $i=0;$i<@gc;$i++){
                                if(exists $hash{$g}{circRNA}{$gc[$i]}){
                                        $flag++;
                                        if($flag ==1){
                                                $cinfo = $hash{$g}{circRNA}{$gc[$i]};
                                        }else{
                                                $cinfo .= ";".$hash{$g}{circRNA}{$gc[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$cinfo\n";
				$gcCyto{$A[0]}{$A[-1]}=$cinfo;
                        }else{
                                print OUT "--\n";
                        }
                }else{
                        print OUT "--\n";
                }
        }
        close (IN);
	close (OUT);
	&getCyto($gsinter,$gsatt,"ms",\%gsCyto);
	&getCyto($glinter,$glatt,"ml",\%glCyto);
	&getCyto($gcinter,$gcatt,"mc",\%gcCyto);
	}
}else{
	foreach my $g (keys %h){
 		open (IN,"$h{$g}{gene}");
        	open (OUT,">$gdir/$g.DEG.xls");
	        my $gsinter = "$gdir/$g.ms.Interaction.list";
        	my $gsatt = "$gdir/$g.ms.Attribution.list";
	        my $glcinter = "$gdir/$g.ml.Cis.Interaction.list";
	        my $glcatt = "$gdir/$g.ml.Cis.Attribution.list";
		my $gltinter = "$gdir/$g.ml.Trans.Interaction.list";
		my $gltatt = "$gdir/$g.ml.Trans.Attribution.list";
	        my $gcinter = "$gdir/$g.mc.Interaction.list";
	        my $gcatt = "$gdir/$g.mc.Attribution.list";
	        my (%gsCyto,%glcCyto,%gltCyto,%gcCyto);
        	while(<IN>){
        	        chomp;
			if(-e "$l2mt"){
	                	if(/^#/){print OUT "$_\tdiff_target_miRNA\tCis.diff_target_lncRNA\tTrans.diff_target_lncRNA\tdiff_target_circRNA\n";next;}
			}else{
				if(/^#/){print OUT "$_\tdiff_target_miRNA\tCis.diff_target_lncRNA\tdiff_target_circRNA\n";next;}
			}
                my @A = split /\t/,$_;
                print OUT "$_\t";
                if(exists $g2mi{$A[0]}){
                        my @gs = split /\;/,$g2mi{$A[0]};
                        my $sinfo;my $flag =0;
                        for (my $i=0;$i<@gs;$i++){
                                if(exists $hash{$g}{sRNA}{$gs[$i]}){
                                        $flag++;
                                        if($flag==1){
                                                $sinfo = $hash{$g}{sRNA}{$gs[$i]};
                                        }else{
                                                $sinfo .= ";".$hash{$g}{sRNA}{$gs[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$sinfo\t";
                                $gsCyto{$A[0]}{$A[-1]}=$sinfo;
                        }else{
                                print OUT "--\t";
                        }
                }else{
                        print OUT "--\t";
		}
		
		if(exists $m2l{$A[0]}{Cis}){
                        my @gl = split /;/,$m2l{$A[0]}{Cis};
                        my $linfo;my $flag =0;
                        for (my $i=0;$i<@gl;$i++){
                                if(exists $hash{$g}{lncRNA}{$gl[$i]}){
                                        $flag++;
                                        if($flag ==1){
                                                $linfo = $hash{$g}{lncRNA}{$gl[$i]};
                                        }else{
                                                $linfo .= ";".$hash{$g}{lncRNA}{$gl[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$linfo\t";
                                $glcCyto{$A[0]}{$A[-1]}=$linfo;
                        }else{
                                print OUT "--\t";
                        }
                }else{
                        print OUT "--\t";
                }
		if(-e "$l2mt"){
			if(exists $m2l{$A[0]}{Trans}){
                        my @gl = split /;/,$m2l{$A[0]}{Trans};
                        my $linfo;my $flag =0;
                        for (my $i=0;$i<@gl;$i++){
                                if(exists $hash{$g}{lncRNA}{$gl[$i]}){
                                        $flag++;
                                        if($flag ==1){
                                                $linfo = $hash{$g}{lncRNA}{$gl[$i]};
                                        }else{
                                                $linfo .= ";".$hash{$g}{lncRNA}{$gl[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$linfo\t";
                                $gltCyto{$A[0]}{$A[-1]}=$linfo;
                        }else{
                                print OUT "--\t";
                        }
                }else{
                        print OUT "--\t";
                }
	}
		
		if(exists $g2c{$A[0]}){
                        my @gc = split /;/,$g2c{$A[0]};
                        my $cinfo;my $flag =0;
                        for (my $i=0;$i<@gc;$i++){
                                if(exists $hash{$g}{circRNA}{$gc[$i]}){
                                        $flag++;
                                        if($flag ==1){
                                                $cinfo = $hash{$g}{circRNA}{$gc[$i]};
                                        }else{
                                                $cinfo .= ";".$hash{$g}{circRNA}{$gc[$i]};
                                        }
                                }
                        }
                        if($flag >0){
                                print OUT "$cinfo\n";
                                $gcCyto{$A[0]}{$A[-1]}=$cinfo;
                        }else{
                                print OUT "--\n";
                        }
                }else{
                        print OUT "--\n";
                }
        }
	
        close (IN);
        close (OUT);
        &getCyto($gsinter,$gsatt,"ms",\%gsCyto);
        &getCyto($glcinter,$glcatt,"ml",\%glcCyto);
	&getCyto($gltinter,$gltatt,"ml",\%gltCyto);
        &getCyto($gcinter,$gcatt,"mc",\%gcCyto);
	}
}


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
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
	open (OUT,">$name.veen.r");
	
	
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
	open (VENN, ">$name.set.xls");
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
		for my $deg_id ($list[$i]){
			$deg_id = "'".$deg_id."'";
			push @{$DEG{$label[$i-1]}}, $deg_id;
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
		$more_opts.= "    cat.cex = 0.4,\n";
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
	fill = c($color_content),
	height=800,
	width=800,
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
 Description:	This program is used to .... 
       Usage:
		Options:
		-in <dir>	input dir ,DEG_Analysis,force
		-od <dir>	output dir , Combine/Cytoscape,force
		-data_cfg <file>	data_cfg,force
		-detail_cfg <file>	detail_cfg,force
		-l2m <file>	lncRNA2mRNA target file, l2m can exist with (l2mc and l2mt) at the same time.
		-l2mc <file>	Cis target file, (l2mc and l2mt) can exist with l2m at the same time.
		-l2mt <file>	Trans target file, l2mt can not exist 
		-mi2l <file>	miRNA2lncrNA target file,
		-mi2c <file>	miRNA2circRNA target file,
		-mi2g <file>	miRNA2gene target file,
		-chostgene <file>	circRNA hostgene file
		-h		help

USAGE
	print $usage;
	exit;
}
