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
my $Rscript="/share/nas2/genome/biosoft/R/3.1.1/bin/Rscript";
my $fBam;
my $input_file_format = "depth"; ## can be sam, bam or depth
my $notename = `hostname`; chomp $notename;
my $fKey;
my $dOut;

my $x_title;
my $y_title;
my $title;
my $chrnum=0;
my $nchro;
my $window = 10000;
my $fai;
my $queue;
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fBam,
				"a:s"=>\$fai,
				"f:s"=>\$input_file_format,
				"o:s"=>\$dOut,
				"k:s"=>\$fKey,
				"n:s"=>\$nchro,
				"x:s"=>\$x_title,
				"y:s"=>\$y_title,
				"t:s"=>\$title,
				"q:s"=>\$queue,				
		
				) or &USAGE;
&USAGE unless ($fBam and $fKey );
&USAGE if ($input_file_format eq 'bam' and !$fai );

$dOut||="./";
mkdir $dOut unless (-d $dOut) ;
mkdir "$dOut/work_sh" unless (-d "$dOut/work_sh") ;
$x_title||= "Chromsome position";
$y_title||= "Median of read deinsity(log(2))";
$title||= "Genomewide distribution of read coverage";
$queue||="middle.q";

$fBam = &ABSOLUTE_DIR($fBam) unless ($input_file_format eq 'depth') ;
$fai = &ABSOLUTE_DIR($fai) if (defined $fai);
$dOut= &ABSOLUTE_DIR($dOut);
#-------------------------------------------------------------------
# Global value
#-------------------------------------------------------------------

#my %readDensity = ();  ## read density
#my %chroLength = ();  ## chromosome length

#-------------------------------------------------------------------
# Get depth info
#-------------------------------------------------------------------

my $prefix = "$dOut/$fKey";
my ($plus,$minus);

if ($input_file_format eq 'depth') {
        ($plus,$minus)=split(/:/,$fBam);
}elsif($input_file_format eq 'bam'){
        ## sort bam file
		$plus="$dOut/plus.txt";
		$minus="$dOut/minus.txt";
		open (OUT,">$dOut/work_sh/plot_Read_$fKey.sh") or die $!;
		print OUT "bedtools genomecov -ibam  $fBam  -g $fai -d -strand +  |awk '\$3!=0' >$plus && \n";
		print OUT "bedtools genomecov -ibam  $fBam  -g $fai -d -strand -  |awk '\$3!=0' >$minus && \n";
		close OUT;
		&Cut_shell_qsub("$dOut/work_sh/plot_Read_$fKey.sh",  18, "15G", $queue);

}elsif($input_file_format eq 'sam'){
        my $outbamfile = "$prefix.bam";
		$plus="$dOut/plus.txt";
		$minus="$dOut/minus.txt";
        system("samtools view $fBam > $outbamfile");
		open (OUT,">$dOut/work_sh/plot_Read_$fKey.sh") or die $!;        
		print OUT "bedtools genomecov -ibam  $outbamfile  -g $fai -d -strand +  |awk '\$3!=0' >$plus && \n";
		print OUT "bedtools genomecov -ibam  $outbamfile  -g $fai -d -strand -  |awk '\$3!=0' >$minus && \n";
		close OUT;
		&Cut_shell_qsub("$dOut/work_sh/plot_Read_$fKey.sh",  18, "15G", $queue);
}else{
        die "the input file format is wrong\n";
}

&chrsort("$plus","$plus.xls");
&chrsort("$minus","$minus.xls");
&readdensity("$plus.xls","$plus.stat");
&readdensity("$minus.xls","$minus.stat");


#------------------------------------------------------------------
# plot calling r scripts 
#------------------------------------------------------------------

#$nchro = 20 if ($nchro == 0 && $total_chro > 100);

if (defined $nchro) {
	my $cmd = "$Rscript $Bin/plot_coverage_distribution_v2.r  --input1 $plus.stat --input2 $minus.stat --output $dOut/$fKey.coverage.png  -w $window -n $nchro -x \"$x_title\" -y \"$y_title\" -t \"$title\"";
	print $cmd,"\n";
	process_cmd($cmd);
}
else {
	$nchro=($chrnum>20)?10:$chrnum;
	my $cmd = "$Rscript $Bin/plot_coverage_distribution_v2.r  --input1 $plus.stat --input2 $minus.stat --output $dOut/$fKey.coverage.png  -w $window -n $nchro -x \"$x_title\" -y \"$y_title\" -t \"$title\"";
	print $cmd,"\n";
	process_cmd($cmd);
}


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub chrsort{
	my ($input,$output)=@_;
	my %chrname;
	my %cover;
	open (IN,"$input") or die $!;

	$chrnum=0;
	while (<IN>) {
		next if $_=~/^\s*$|^#/;
		chomp;
		my($chr,$inform)=split(/\s+/,$_,2);
		my $name=$chr;
		#$name=~s/.*(\d+)$/$1/;
		if ($name=~/(\d*\.*\d+)$/) {
			$name=$1;
		}else{
			$name=100000000;
		}

		if (exists $chrname{$name}{$chr}) {
			$chrname{$name}{$chr}.="$chr\t$inform\n";
			$cover{$chr}++;
		}
		else {
			$chrnum++;
			$chrname{$name}{$chr}="$chr\t$inform\n";
			$cover{$chr}=1;
			last if $chrnum>20;
		}
	}
	close(IN);
	#my $mark=0;
	if ($chrnum<=20) {
		open (OUT,">$output") or die $!;
		foreach my $name (sort {$a<=>$b} keys  %chrname) {
			foreach my $chr (sort  keys  %{$chrname{$name}}) {
				print OUT $chrname{$name}{$chr};
			}
		}
		close(OUT);
		my @chrlen=sort {$b<=>$a} values %cover;
		$window=1000 if ($chrlen[0]<10000000);		
	}
	else {
		system"cp $input $output";
		my @chrlen=sort {$b<=>$a} values %cover;
	#	foreach my $chr (sort {$cover{$b}<=>$cover{$a}} keys  %cover) {
	#		$mark=$cover{$chr} if ($mark==0) ;
	#		foreach my $name (sort {$a<=>$b} keys  %chrname) { 
	#			print OUT $chrname{$name}{$chr} if (exists $chrname{$name}{$chr}) ;
	#		}
	#	}
		$window=1000 if ($chrlen[0]<10000000);
	}


}

sub readdensity {
	my ($input,$output)=@_;

	my @tmp = ();
	my $cur_pos_index = 0;
	my $last_chr = "";
	my $total_chro = 0;

	open (IN,"$input") or die $!;
	open (OUT,">$output") or die $!;
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^$/ || /^;/ || /^\#/) ;
		
		my ($chr, $pos, $depth) = split;
		
		if ($chr ne $last_chr) {
			## output last chromosome's info any more
			if ($last_chr ne "") {
				my $midian_cov = mid(\@tmp);

				print OUT join("\t",
					$last_chr,
					$cur_pos_index,
					log($midian_cov+1)/log(2),
				),"\n";
			}

			$cur_pos_index = 1;
			@tmp = ();

			push @tmp, $depth;
			$last_chr = $chr;
			$total_chro++;
			
			next;
		}

		if ($pos > $cur_pos_index * $window) {

			my $midian_cov = mid(\@tmp);
			
			print OUT join("\t",
				$chr,
				$cur_pos_index,
				log($midian_cov+1)/log(2),
			),"\n";

			while ((++$cur_pos_index) * $window < $pos) {
				
				print OUT join("\t",
					$chr,
					$cur_pos_index,
					0,
				),"\n";
			}

			@tmp = ();
		}
		
		push @tmp, $depth;

		$last_chr = $chr;
		
	}
	close (IN) ;
	close (OUT) ;
}
sub Cut_shell_qsub {#Cut shell for qsub 1000 line one file
	# &Cut_shell_qsub($shell,$cpu,$vf,$queue);
	my $shell = shift;
	my $cpu = shift;
	my $vf = shift;
	my $queue = shift;

	my $line = `less -S $shell |wc -l `;
	if ($line<=1000) {
		if ($notename=~/cluster/) {
			system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
		else
		{
			system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $shell --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
		}
	}
	if ($line>1000) {
		my @div=glob "$shell.div*";
		foreach (@div) {
			if (-e $_) {
				system "rm $_";
			}
		}
		@div=();
		my $div_index=1;
		my $line_num=1;
		open IN,"$shell" || die;
		while (<IN>) {
			chomp;
			open OUT,">>$shell.div.$div_index" || die;
			if ($line_num<1000) {
				print OUT "$_\n";
				$line_num++;
			}
			else {
				print OUT "$_\n";
				$div_index++;
				$line_num=1;
				close OUT;
			}
		}
		if ($line_num!=0) {
			close OUT;
		}
		@div=glob "$shell.div*";
		foreach my $div_file (@div) {
			if ($notename=~/cluster/) {
				system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
			else
			{
				system "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $div_file --queue $queue --maxproc $cpu --resource vf=$vf --reqsub --independent";
			}
		}
	}
}
sub mid {
	my $arr = shift;

	return 0 if (@$arr == 0) ;
	my @list = sort @$arr;
	return $list[(@list-@list%2)/2]; 
} 

sub process_cmd {#
	my ($cmd) = @_;

	print STDERR "CMD: $cmd\n";
	my $start_time = time;

	my $ret = system($cmd);
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}

	my $end_time = time;
	print STDERR "CMD Done, elapsed time ( ", $end_time - $start_time, "s )\n\n";

	return;
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
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Li Tiancheng <litc\@biomarker.com.cn> <ltc_gs\@qq.com>
Modified: Wang Yajing <wangyj\@biomarker.com.cn>

	画染色体reads分布图。
Usage:
  Options:
  -i <file>  input sam/bam file, result of samtools, forced.
  -f <str>   input file format, sam|bam, default [depth].
  -a <str>   input genome.fa.fai, forced .
  -k <str>   key of output file , forced.
  -o <dir>   direction where restult produced.
  -n <int>   if the gonome is in scaffold, choose the top n chromsomes to draw, 0 means plot all chromsomes, default [0].
  -x <str>   title of x axis, optional.
  -y <str>   titile of y axis, optional.
  -t <str>   titile of graph, optional.
  -q <str>   the queue is used for qsub jobs 
USAGE
	print $usage;
	exit;
}
