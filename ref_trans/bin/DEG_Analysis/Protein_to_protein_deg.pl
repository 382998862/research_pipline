#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $BEGIN_TIME=time();
my $version="1.0.1";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($flinks,$unigene,$ftest,$odir,$protein,$blast,$id,$clade_or_taxid,$queues,$cfg,$idir,$ppi);

GetOptions(
				"help|?" =>\&USAGE,
				"od:s"=>\$odir,
				"ppi:s"=>\$ppi,
				"idir:s"=>\$idir,
			) or &USAGE;
&USAGE unless ( $idir);
mkdir $odir unless -d "$odir";
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------


    my @list = glob "$idir/*_vs_*/*DEG*final.xls";
    for my $l (@list) {
        my $name = (split /\./,basename($l))[0];
        &abstract_interacts($l, "$ppi", "$odir/$name.DEG.detail.txt");
        &make_sif("$odir/$name.DEG.detail.txt", "$odir/$name.ppi.cytoscapeInput.sif");
    }


	` cat $idir/*_vs_*/*DEG*final.xls |grep -v "#" |cut -f 1|sort |uniq >$odir/used_gene.list `;
    &abstract_interacts("$odir/used_gene.list", "$ppi", "$odir/ppi_qurey.ppi.detail.txt");
    &make_sif("$odir/ppi_qurey.ppi.detail.txt", "$odir/ppi_qurey.ppi.cytoscapeInput.sif");
    


########################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################
sub abstract_interacts {
    my ($fId, $fIn, $fOut) = @_;

    my %ids;
    open (ID,$fId) or die $!;
    while (<ID>) {
        chomp;
        next if (/^\s*$/ or /^#/);
        $ids{(split /\t/)[0]}=1;
    }
    close ID;

    open (IN,$fIn) or die $!;
    open (OUT,">$fOut") or die $!;
    while (<IN>) {
        chomp;
        next if (/^\s*$/);
        if (/^#/) {
            print OUT "$_\n";
        } else {
            my ($k1, $k2) =(split /\t/)[0,2];
            print OUT "$_\n" if ( exists $ids{$k1} && exists $ids{$k2} );
        }
    }
    close IN;
    close OUT;
}

sub make_sif {
    my ($fIn, $fOut) = @_;
    #print "$fIn,$fOut\n";
    my %ppi;
    open (PPI,$fIn) or die $!;
    while (<PPI>) {
        chomp;
        next if (/^\s*$/ or /^#/);
        my ($a,$b,$c) = (split /\t/)[0,1,2];
        ($a,$c) = ($a lt $c) ? ($a,$c) : ($c,$a);
        $ppi{$a}{$c} = $b;
    }
    close PPI;

    open (OUT,">$fOut") or die $!;
    for my $i (sort keys %ppi) {
        for my $j (sort keys %{$ppi{$i}}) {
            print OUT "$i\tpp\t$j\n";
        }
    }
    close OUT;
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
     Contact:	Wang Yajing <wangyj\@biomarker.com.cn> 
Program Date:	2015.09.06
      Modify:	
 Description:	This program is used to ......
       Usage:
		Options:
		-h		help
		-idir		input deg dir
		-od		output deg ppi dir
		-ppi		input PPI.txt
USAGE
	print $usage;
	exit;
}

sub detail_cfg_read {
    #&log_current_time("detail config check:");
    my ($cfg_file, $detail_cfg) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+$/ or /^#/);
        my ($key, $value) = (split /\s+/)[0,1];
	next unless(defined $key);
        $detail_cfg->{$key} = $value;
        if ($key eq 'Project_name' or $key eq 'Customer_info' or $key eq 'Project_id' or $key eq 'Project_key') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Known_unigene' or $key eq 'Known_pep') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }

        if ($key eq 'Known_anno') {
            #die "$key: $value is not illegal!\n" unless (-e "$value/02.gene-annotation" and -e "$value/Result");
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Ref_seq' or $key eq 'Ref_ann') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }
        if ($key =~/^SNP_/) {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'nr' or $key eq 'Swissprot' or $key eq 'Kegg' or $key eq 'Pfam' or $key eq 'Cog' or $key eq 'Kog'  or $key eq 'eggNOG'  ) {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{anno}{$key} = $value;
			die "Kegg database is wrong ,should  be kobas data\n " if ($key eq 'Kegg' and $value!~/kobas/);
        }
        if ($key eq 'Queue_type1') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'Queue_type2') {
            $detail_cfg->{$key} = $value;
        }
        print "$key: $value\n" if (exists $detail_cfg->{$key});
    }
    close CFG;
    die "Must choose Queue_type1 and Queue_type2 !\n" unless (exists $detail_cfg->{Queue_type1} or exists $detail_cfg->{Queue_type2});
   # &log_current_time("detail config check done.");
}
