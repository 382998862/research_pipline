#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd qw(abs_path);
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
my ($idir,$odir,$cfg,$l2mc,$l2mt,$mi2l,$mi2c,$mi2g,$chostgene,$nl);

GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$odir,
				"in:s"=>\$idir,
				"cfg:s"=>\$cfg,
				"l2mc:s"=>\$l2mc,
				"l2mt:s"=>\$l2mt,
				"mi2l:s"=>\$mi2l,
				"mi2c:s"=>\$mi2c,
				"mi2g:s"=>\$mi2g,
				"chostgene:s"=>\$chostgene,
				"name2list:s"=>\$nl,
				) or &USAGE;
&USAGE unless ($idir and $odir and $cfg);
# ------------------------------------------------------------------
$cfg = abs_path($cfg);
$idir = abs_path($idir);
$odir = abs_path($odir);
`mkdir $odir` unless (-d "$odir");
$l2mc||="$idir/Cis_target_gene.xls";
$l2mt||="$idir/Trans_target_gene.xls";
$mi2l||="$idir/lncRNA.mir2target.list";
$mi2c||="$idir/circRNA.mir2target.list";
$mi2g||="$idir/gene.mir2target.list";
$chostgene||="$idir/circRNA.source.xls";
my %name;
if(defined $nl){
	open (IN,$nl);
	while(<IN>){
		chomp;
		next if(/^#/);
		my @t = split /\t/,$_;
		$name{$t[0]}=$t[1];
		my $up = $t[1]."(up)";
		my $down = $t[1]."(down)";
		$name{"$t[0](up)"}=$up;
		$name{"$t[0](down)"}=$down;
	}
	close(IN);
}
my %sample = &relate($cfg);
my @diff;my %vs;
open (CFG,"$cfg") or die $!;
print ("read cfg and diff group is:");
while (<CFG>){
	chomp;
	next if(/^#|^$|^\s$/);
	my @tmp = split /\s+/,$_;
	if($tmp[0] eq "Diff"){
		push @diff,$tmp[1];
		print "$tmp[1]\n";
	}
}
close(CFG);

my (%l2m,%m2l,%mi2l,%l2mi,%mi2c,%c2mi,%mi2g,%g2mi,%c2g,%g2c);

###read lncRNA2mRNA file and then get %l2m and %m2l
print ("1:L2M...... %l2m,%m2l\n");
if(-e $l2mc){
	open (L2MC,"$l2mc") or die $!;
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
}
if(-e $l2mt){
	open (L2MT,"$l2mt") or die $!;
	while (<L2MT>){
		chomp;
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


###read miRNA2lncRNA file and then get %mi2l and %l2mi
if(-e $mi2l){
	print "2:MI2L...... %mi2l,%l2mi\n";
	open (MI2L,"$mi2l") or die $!;
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
}
###read miRNA2circRNA file
if(-e $mi2c){
	print ("3:MI2C...... %mi2c,%c2mi\n");
	open(MI2C,"$mi2c") or die $!;
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
}
###read miRNA2gene file
if(-e $mi2g){
	print ("4:MI2G...... %mi2g,%g2mi\n");
	open (MI2G,"$mi2g") or die $!;
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
}
###read circRNA hostgene file
if(-e $chostgene){
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
}
print ("read target file progress finish!!!\n\n");

print "get diff result file and creat hash\n\n";
print "#group\tRNA_type\tdiff_files\n";
my (%h,%hash);
foreach my $vs(@diff){
	&switch($vs);
	foreach my $k ("sRNA","lncRNA","circRNA","gene"){
		$h{$vs}{$k} = "$idir/$k.$vs{$vs}{$k}.DEG_final.xls";
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
                	if(/^#/){
				my ($i1,$d1)=split /\t/,$_,2;
				if(defined $nl){
					print OUT "$i1\tSymbol\t$d1\tdiff_Cis_target_gene\tSymbol\tdiff_Trans_target_gene\tSymbol\tdiff_target_miRNA\n";
				}else{
					print OUT "$i1\t$d1\tdiff_Cis_target_gene\tdiff_Trans_target_gene\tdiff_target_miRNA\n";
				}
				next;
			}
		}else{
			if(/^#/){
				my ($i1,$d1)=split /\t/,$_,2;
				if(defined $nl){
					print OUT "$i1\tSymbol\t$d1\tdiff_Cis_target_gene\tSymbol\tdiff_target_miRNA\n";
				}else{
					print OUT "$i1\t$d1\tdiff_Cis_target_gene\tdiff_target_miRNA\n";
				}
				next;
			}
		}
	    	my @A = split /\t/,$_;
		my ($n,$m)=split /\t/,$_,2;
		if(defined $nl){
                        if(exists $name{$n}){
                                print OUT "$n\t$name{$n}\t$m\t";
                        }else{
                                print OUT "$n\t$n\t$m\t";
                    	}
		}else{
			print OUT "$n\t$m\t";
		}
        	        if(exists $l2m{$A[0]}{Cis}){
                	        my @lm = split /\;/,$l2m{$A[0]}{Cis};
                        	my $ginfo;my $flag =0;my $ninfo;
	                        for (my $i=0;$i<@lm;$i++){
        	                        if(exists $hash{$test}{gene}{$lm[$i]}){
                	                        $flag++;
						if (defined $nl){
						  if(exists $name{$hash{$test}{gene}{$lm[$i]}}){
                        		                if($flag==1){
                                		                $ginfo = $hash{$test}{gene}{$lm[$i]};
								$ninfo = $name{$hash{$test}{gene}{$lm[$i]}};
	                                	        }else{
        	                                	        $ginfo .= ";".$hash{$test}{gene}{$lm[$i]};
								$ninfo .= ";".$name{$hash{$test}{gene}{$lm[$i]}};
                	                        	}
						  }else{
							if($flag==1){
								$ginfo = $hash{$test}{gene}{$lm[$i]};
								$ninfo = $hash{$test}{gene}{$lm[$i]};
							}else{
								$ginfo .= ";".$hash{$test}{gene}{$lm[$i]};
								$ninfo .= ";".$hash{$test}{gene}{$lm[$i]};
							}
						  }
						}else{
							if($flag==1){
								$ginfo = $hash{$test}{gene}{$lm[$i]};
							}else{
								$ginfo .= ";".$hash{$test}{gene}{$lm[$i]};
							}
						}
                        	        }
	                        }

        	                if($flag >0){
					if(defined $nl){
                	                	print OUT "$ginfo\t$ninfo\t";
					}else{
						print OUT "$ginfo\t";
					}
                        	        $lmcCyto{$A[0]}{$A[-1]}=$ginfo;
	                        }else{
					if(defined $nl){
        	                        	print OUT "--\t--\t";
					}else{
						print OUT "--\t";
					}
                	        }
	                }else{
				if(defined $nl){
        	                	print OUT "--\t--\t";
				}else{
					print OUT "--\t";
				}
                	}
		
			if(-e $l2mt){
				if(exists $l2m{$A[0]}{Trans}){
                	        	my @lm = split /\;/,$l2m{$A[0]}{Trans};
	                        	my $ginfo;my $flag =0;my $ninfo;
		                        for (my $i=0;$i<@lm;$i++){
        		                        if(exists $hash{$test}{gene}{$lm[$i]}){
                		                        $flag++;
							if(defined $nl){
							  if(exists $name{$hash{$test}{gene}{$lm[$i]}}){
	                        		                if($flag==1){
        	                        		                $ginfo = $hash{$test}{gene}{$lm[$i]};
									$ninfo = $name{$hash{$test}{gene}{$lm[$i]}};
	        	                                	}else{
        	        	                                	$ginfo .= ";".$hash{$test}{gene}{$lm[$i]};
									$ninfo .= ";".$name{$hash{$test}{gene}{$lm[$i]}};
	                		                        }
							  }else{
								if($flag==1){
                                                                        $ginfo = $hash{$test}{gene}{$lm[$i]};
									$ninfo = $hash{$test}{gene}{$lm[$i]};
                                                                }else{
                                                                        $ginfo .= ";".$hash{$test}{gene}{$lm[$i]};
									$ninfo .= ";".$hash{$test}{gene}{$lm[$i]};
                                                                }
							  }
							}else{
								if($flag==1){
									$ginfo = $hash{$test}{gene}{$lm[$i]};
								}else{
									$ginfo .= ";".$hash{$test}{gene}{$lm[$i]};
								}
							}
        	                	        }
	        	                }
	        	                if($flag >0){
						if(defined $nl){
	        	        	                print OUT "$ginfo\t$ninfo\t";
						}else{
							print OUT "$ginfo\t";
						}
                	        	        $lmtCyto{$A[0]}{$A[-1]}=$ginfo;
	                	        }else{
						if(defined $nl){
        	                	        	print OUT "--\t--\t";
						}else{
							print OUT "--\t";
						}
                	        	}
	                	}else{
					if(defined $nl){
	        	               		print OUT "--\t--\t";
					}else{
						print OUT "--\t";
					}
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
        &getCyto($lmcinter,$lmcatt,"lm",\%lmcCyto,"1");
	&getCyto($lmtinter,$lmtatt,"lm",\%lmtCyto,"1");
        &getCyto($lsinter,$lsatt,"ls",\%lsCyto,"1");
        close (IN);
        close (OUT);
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
                if(/^#/){
			if(defined $nl){
				print OUT "$_\tdiff_host_gene\tSymbol\tdiff_target_miRNA\n";
			}else{
				print OUT "$_\tdiff_host_gene\tdiff_target_miRNA\n";
			}
			next;
		}
                my @A = split /\t/,$_;
                print OUT "$_\t";
                if(exists $c2g{$A[0]}){
                        my @cg = split /\;/,$c2g{$A[0]};
                        my $ginfo;my $flag =0;my $ninfo;
                        for (my $i=0;$i<@cg;$i++){
                                if(exists $hash{$circ}{gene}{$cg[$i]}){
                                        $flag++;
					if(defined $nl){
					  if(exists $name{$hash{$circ}{gene}{$cg[$i]}}){
	                                        if($flag==1){
        	                                        $ginfo = $hash{$circ}{gene}{$cg[$i]};
							$ninfo = $name{$hash{$circ}{gene}{$cg[$i]}};
                	                        }else{
                        	                        $ginfo .= ";".$hash{$circ}{gene}{$cg[$i]};
							$ninfo .= ";".$name{$hash{$circ}{gene}{$cg[$i]}};
                                	        }
					  }else{
						if($flag==1){
                                                        $ginfo = $hash{$circ}{gene}{$cg[$i]};
							$ninfo = $hash{$circ}{gene}{$cg[$i]};
                                                }else{
                                                        $ginfo .= ";".$hash{$circ}{gene}{$cg[$i]};
							$ninfo .= ";".$hash{$circ}{gene}{$cg[$i]};
                                                }
					  }
					}else{
						if($flag==1){
							$ginfo = $hash{$circ}{gene}{$cg[$i]};
						}else{
							$ginfo .= ";".$hash{$circ}{gene}{$cg[$i]};
						}
					}
                                }
                        }
                        if($flag >0){
				if(defined $nl){
	                                print OUT "$ginfo\t$ninfo\t";
				}else{
					print OUT "$ginfo\t";
				}
				$cgCyto{$A[0]}{$A[-1]}=$ginfo;
                        }else{
				if(defined $nl){
	                                print OUT "--\t--\t";
				}else{
					print OUT "--\t";
				}
                        }
                }else{
			if(defined $nl){
	                        print OUT "--\t--\t";
			}else{
				print OUT "--\t";
			}
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
	&getCyto($cginter,$cgatt,"cm",\%cgCyto,"1");
        &getCyto($csinter,$csatt,"cs",\%csCyto,"1");
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
                if(/^#/){
			if(defined $nl){
				print OUT "$_\tdiff_target_gene\tSymbol\tdiff_target_lncRNA\tSymbol\tdiff_target_circRNA\n";
			}else{
				print OUT "$_\tdiff_target_gene\tdiff_target_lncRNA\tdiff_target_circRNA\n";
			}
			next;
		}
                my @A = split /\t/,$_;
                print OUT "$_\t";
                if(exists $mi2g{$A[0]}){
                        my @sg = split /\;/,$mi2g{$A[0]};
                        my $ginfo;my $flag =0;my $ninfo;
                        for (my $i=0;$i<@sg;$i++){
                                if(exists $hash{$mi}{gene}{$sg[$i]}){
                                        $flag++;
					if(defined $nl){
					  if(exists $name{$hash{$mi}{gene}{$sg[$i]}}){
	                                        if($flag==1){
        	                                        $ginfo = $hash{$mi}{gene}{$sg[$i]};
							$ninfo = $name{$hash{$mi}{gene}{$sg[$i]}};
                	                        }else{
                        	                        $ginfo .= ";".$hash{$mi}{gene}{$sg[$i]};
							$ninfo .= ";".$name{$hash{$mi}{gene}{$sg[$i]}};
                                	        }
					  }else{
						if($flag==1){
                                                        $ginfo = $hash{$mi}{gene}{$sg[$i]};
							$ninfo = $hash{$mi}{gene}{$sg[$i]};
                                                }else{
                                                        $ginfo .= ";".$hash{$mi}{gene}{$sg[$i]};
							$ninfo .= ";".$hash{$mi}{gene}{$sg[$i]};
                                                }
					  }
					}else{
						if($flag==1){
							$ginfo = $hash{$mi}{gene}{$sg[$i]};
						}else{
							$ginfo .= ";".$hash{$mi}{gene}{$sg[$i]};
						}
					}
                                }
                        }
                        if($flag >0){
				if(defined $nl){
	                                print OUT "$ginfo\t$ninfo\t";
				}else{
					print OUT "$ginfo\t";
				}
				$sgCyto{$A[0]}{$A[-1]}=$ginfo;
                        }else{
				if(defined $nl){
	                                print OUT "--\t--\t";
				}else{
					print OUT "--\t";
				}
                        }
                }else{
			if(defined $nl){
	                        print OUT "--\t--\t";
			}else{
				print OUT "--\t";
			}
                }

                if(exists $mi2l{$A[0]}){
                        my @sl = split /;/,$mi2l{$A[0]};
                        my $linfo;my $flag =0;my $ninfo;
                        for (my $i=0;$i<@sl;$i++){
                                if(exists $hash{$mi}{lncRNA}{$sl[$i]}){
                                        $flag++;
					if(defined $nl){
					  if(exists $name{$hash{$mi}{lncRNA}{$sl[$i]}}){
	                                        if($flag ==1){
        	                                        $linfo = $hash{$mi}{lncRNA}{$sl[$i]};
							$ninfo = $name{$hash{$mi}{lncRNA}{$sl[$i]}};
                	                        }else{
                        	                        $linfo .= ";".$hash{$mi}{lncRNA}{$sl[$i]};
							$ninfo .= ";".$name{$hash{$mi}{lncRNA}{$sl[$i]}};
                                	        }
					  }else{
						if($flag ==1){
                                                        $linfo = $hash{$mi}{lncRNA}{$sl[$i]};
							$ninfo = $hash{$mi}{lncRNA}{$sl[$i]};
                                                }else{
                                                        $linfo .= ";".$hash{$mi}{lncRNA}{$sl[$i]};
							$ninfo .= ";".$hash{$mi}{lncRNA}{$sl[$i]};
                                                }
					  }
					}else{
						if($flag ==1){
							$linfo = $hash{$mi}{lncRNA}{$sl[$i]};
						}else{
							$linfo .= ";".$hash{$mi}{lncRNA}{$sl[$i]};
						}
					}
                                }
                        }
                        if($flag >0){
				if(defined $nl){
	                                print OUT "$linfo\t$ninfo\t";
				}else{
					print OUT "$linfo\t";
				}
				$slCyto{$A[0]}{$A[-1]}=$linfo;
                        }else{
				if(defined $nl){
	                                print OUT "--\t--\t";
				}else{
					print OUT "--\t";
				}
                        }
                }else{
			if(defined $nl){
                        	print OUT "--\t--\t";
			}else{
				print OUT "--\t";
			}
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
	&getCyto($sginter,$sgatt,"sm",\%sgCyto,"1");
	&getCyto($slinter,$slatt,"sl",\%slCyto,"1");
	&getCyto($scinter,$scatt,"sc",\%scCyto,"1");
}

print "core-gene\n";
`mkdir $odir/gene` unless (-d "$odir/gene");
my $gdir = "$odir/gene";
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
                	if(/^#/){
				my ($i,$d)=split /\t/,$_,2;
				if(defined $nl){
					print OUT "$i\tSymbol\t$d\tdiff_target_miRNA\tCis.diff_target_lncRNA\tSymbol\tTrans.diff_target_lncRNA\tSymbol\tdiff_target_circRNA\n";
				}else{
					print OUT "$i\t$d\tdiff_target_miRNA\tCis.diff_target_lncRNA\tTrans.diff_target_lncRNA\tdiff_target_circRNA\n";
				}
				next;
			}
		}else{
			if(/^#/){
				my ($i,$d)=split /\t/,$_,2;
				if(defined $nl){
					print OUT "$i\tSymbol\t$d\tdiff_target_miRNA\tCis.diff_target_lncRNA\tSymbol\tdiff_target_circRNA\n";
				}else{
					print OUT "$i\t$d\tdiff_target_miRNA\tCis.diff_target_lncRNA\tdiff_target_circRNA\n";
				}
				next;
			}
		}
                my @A = split /\t/,$_;
		my ($m,$p)=split /\t/,$_,2;
		if(defined $nl){
 		  if(exists $name{$m}){
			print OUT "$m\t$name{$m}\t$p\t";
		  }else{
			print OUT "$m\t$m\t$p\t";
		  }
		}else{
			print OUT "$m\t$p\t";
		}
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
                        my $linfo;my $flag =0;my $ninfo;
                        for (my $i=0;$i<@gl;$i++){
                                if(exists $hash{$g}{lncRNA}{$gl[$i]}){
                                        $flag++;
					if(exists $name{$hash{$g}{lncRNA}{$gl[$i]}}){
                                        	if($flag ==1){
                                                	$linfo = $hash{$g}{lncRNA}{$gl[$i]};
							$ninfo = $name{$hash{$g}{lncRNA}{$gl[$i]}};
	                                        }else{
        	                                        $linfo .= ";".$hash{$g}{lncRNA}{$gl[$i]};
							$ninfo .= ";".$name{$hash{$g}{lncRNA}{$gl[$i]}};
                	                        }
					}else{
						if($flag ==1){
                                                        $linfo = $hash{$g}{lncRNA}{$gl[$i]};
							$ninfo = $hash{$g}{lncRNA}{$gl[$i]};
                                                }else{  
                                                        $linfo .= ";".$hash{$g}{lncRNA}{$gl[$i]};
							$ninfo .= ";".$hash{$g}{lncRNA}{$gl[$i]};
                                                }
					}
                                }
                        }
                        if($flag >0){
				if(defined $nl){
	                                print OUT "$linfo\t$ninfo\t";
				}else{
					print OUT "$linfo\t";
				}
                                $glcCyto{$A[0]}{$A[-1]}=$linfo;
                        }else{
				if(defined $nl){
	                                print OUT "--\t--\t";
				}else{
					print OUT "--\t";
				}
                        }
                }else{
			if(defined $nl){
                        	print OUT "--\t--\t";
			}else{
				print OUT "--\t";
			}
                }
		if(-e "$l2mt"){
			if(exists $m2l{$A[0]}{Trans}){
                        my @gl = split /;/,$m2l{$A[0]}{Trans};
                        my $linfo;my $flag =0;my $ninfo;
                        for (my $i=0;$i<@gl;$i++){
                                if(exists $hash{$g}{lncRNA}{$gl[$i]}){
                                        $flag++;
					if(defined $nl){
					  if(exists $name{$hash{$g}{lncRNA}{$gl[$i]}}){
	                                        if($flag ==1){
        	                                        $linfo = $hash{$g}{lncRNA}{$gl[$i]};
							$ninfo = $name{$hash{$g}{lncRNA}{$gl[$i]}};
                	                        }else{
                        	                        $linfo .= ";".$hash{$g}{lncRNA}{$gl[$i]};
							$ninfo .= ";".$name{$hash{$g}{lncRNA}{$gl[$i]}};
                                	        }
					  }else{
						if($flag ==1){
                                                        $linfo = $hash{$g}{lncRNA}{$gl[$i]};
							$ninfo = $hash{$g}{lncRNA}{$gl[$i]};
                                                }else{
                                                        $linfo .= ";".$hash{$g}{lncRNA}{$gl[$i]};
							$ninfo .= ";".$hash{$g}{lncRNA}{$gl[$i]};
                                                }
					  }
					}else{
						if($flag ==1){
							$linfo = $hash{$g}{lncRNA}{$gl[$i]};
						}else{
							$linfo .= ";".$hash{$g}{lncRNA}{$gl[$i]};
						}
					}
                                }
                        }
                        if($flag >0){
				if(defined $nl){
	                                print OUT "$linfo\t$ninfo\t";
				}else{
					print OUT "$linfo\t";
				}
                                $gltCyto{$A[0]}{$A[-1]}=$linfo;
                        }else{
				if(defined $nl){
	                                print OUT "--\t--\t";
				}else{
					print OUT "--\t";
				}
                        }
                }else{
			if(defined $nl){
	                        print OUT "--\t--\t";
			}else{
				print OUT "--\t";
			}
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
        &getCyto($gsinter,$gsatt,"ms",\%gsCyto,"1");
        &getCyto($glcinter,$glcatt,"ml",\%glcCyto,"1");
	&getCyto($gltinter,$gltatt,"ml",\%gltCyto,"1");
        &getCyto($gcinter,$gcatt,"mc",\%gcCyto,"1");
	}


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub relate{
        my $cfg=shift;
        my %sample=();
        my @rnas=();
        open(CFG,$cfg)||die $!;
        while(<CFG>){
                chomp;next if($_!~/^Sample/);
                my @tmp=split(/\s+/,$_);shift @tmp;
                if($tmp[0] eq "ID" && scalar(@rnas)==0){
                        shift @tmp;@rnas=@tmp;
                }else{
                        my $id=shift @tmp;
                        for(my $i=0;$i<@tmp;$i++){
                                $sample{$id}{$rnas[$i]}=$tmp[$i];
                        }
                        $sample{$id}{gene}=  $sample{$id}{lncRNA} if(exists $sample{$id}{lncRNA} && !exists $sample{$id}{gene});
                }
        }
        close(CFG);
        return %sample; #relation{W01}{lncRNA}="L01";
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
	my $flag = shift;
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
				if($flag ==1){
					my ($id,$ud,$k) = split /\(|\)/,$p[$i];
					print O1 "$rna\t$id\t$line\n";
					if (!exists $ind{$id}){
						print O2 "$id\t$type\t$ud\n";
						$ind{$id}=1;
					}
				}else{
					print O1 "$rna\t$p[$i]\t$line\n";
					if(!exists $ind{$p[$i]}){
						print O2 "$p[$i]\t$type\ttarget\n";
						$ind{$p[$i]}=1;
					}
				}
			}
		}
	}
	close(O1);close(O2);
	return ($inter,$att);
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
		-cfg <file>	cfg,force
		-l2mc <file>	Cis target file, (l2mc and l2mt) can exist with l2m at the same time.
		-l2mt <file>	Trans target file, l2mt can not exist 
		-mi2l <file>	miRNA2lncrNA target file,
		-mi2c <file>	miRNA2circRNA target file,
		-mi2g <file>	miRNA2gene target file,
		-chostgene <file>	circRNA hostgene file
		-name2list <file>	name2list file
		-h		help

USAGE
	print $usage;
	exit;
}
