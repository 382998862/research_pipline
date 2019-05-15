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
my %CFG=%{readconf("$Bin/../../CFG")};
my $Rscript	 = $CFG{Rscript};
my $bedtools = $CFG{bedtools};
my $samtools = $CFG{samtools};

my $fBam;
my $input_file_format = "depth"; ## can be sam, bam or depth
my $notename = `hostname`; chomp $notename;
my $fKey;
my $dOut;
my $queue;
my $medical;

my $x_title;
my $y_title;
my $title;
my $chrnum=0;
my $nchro;
my $window = 10000;
my $fai;
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
				"queue:s"=>\$queue,
				"medical:s"=>\$medical,
				) or &USAGE;
&USAGE unless ($fBam and $fKey );
&USAGE if ($input_file_format eq 'bam' and !$fai );

$dOut||="./";
mkdir $dOut unless (-d $dOut) ;
mkdir "$dOut/work_sh" unless (-d "$dOut/work_sh") ;
$x_title||= "Chromsome position";
$y_title||= "Median of read deinsity(log(2))";
$title||= "Genomewide distribution of read coverage";
$queue||="medical.q";
$medical||="1";

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
		$plus="$dOut/$fKey.plus.txt";
		$minus="$dOut/$fKey.minus.txt";
		open (OUT,">$dOut/work_sh/plot_Read_$fKey.sh") or die $!;
		print OUT "$bedtools genomecov -ibam  $fBam  -g $fai -d -strand +  |awk '\$3!=0' >$plus && \n";
		print OUT "$bedtools genomecov -ibam  $fBam  -g $fai -d -strand -  |awk '\$3!=0' >$minus && \n";
		close OUT;
		qsubOrDie("$dOut/work_sh/plot_Read_$fKey.sh", "$queue", 18, "15G");

}elsif($input_file_format eq 'sam'){
        my $outbamfile = "$prefix.bam";
		$plus="$dOut/$fKey.plus.txt";
		$minus="$dOut/$fKey.minus.txt";
        runOrDie("$samtools view $fBam > $outbamfile");
		open (OUT,">$dOut/work_sh/plot_Read_$fKey.sh") or die $!;
		print OUT "$bedtools genomecov -ibam  $outbamfile  -g $fai -d -strand +  |awk '\$3!=0' >$plus && \n";
		print OUT "$bedtools genomecov -ibam  $outbamfile  -g $fai -d -strand -  |awk '\$3!=0' >$minus && \n";
		close OUT;
		qsubOrDie("$dOut/work_sh/plot_Read_$fKey.sh","$queue", 18, "15G");
}else{
        die "the input file format is wrong\n";
}

#&chrsort("$plus","$plus.xls");
#&chrsort("$minus","$minus.xls");
#&readdensity("$plus.xls","$plus.stat");
#&readdensity("$minus.xls","$minus.stat");
if(defined $medical){
        &readdensity("$plus","$plus.stat.tmp");
        &readdensity("$minus","$minus.stat.tmp");
        &sort_by_chr("$plus.stat.tmp","$plus.stat");
        &sort_by_chr("$minus.stat.tmp","$minus.stat");
        `rm $plus.stat.tmp $minus.stat.tmp`;
}else{
        &chrsort("$plus","$plus.xls");
        &chrsort("$minus","$minus.xls");
        &readdensity("$plus.xls","$plus.stat");
        &readdensity("$minus.xls","$minus.stat");
}

#------------------------------------------------------------------
# plot calling r scripts
#------------------------------------------------------------------

#$nchro = 20 if ($nchro == 0 && $total_chro > 100);

#if (defined $nchro) {
#	my $cmd = "$Rscript $Bin/plot_coverage_distribution.r  --input1 $plus.stat --input2 $minus.stat --output $dOut/$fKey.png  -w $window -n $nchro -x \"$x_title\" -y \"$y_title\" -t \"$title\" --col \"steelblue1,palegreen\"";
#	print $cmd,"\n";
#	process_cmd($cmd);
#}
#else {
#	$nchro=($chrnum>20)?10:$chrnum;
#	my $cmd = "$Rscript $Bin/plot_coverage_distribution.r  --input1 $plus.stat --input2 $minus.stat --output $dOut/$fKey.png  -w $window -n $nchro -x \"$x_title\" -y \"$y_title\" -t \"$title\" --col \"steelblue1,palegreen\"";
#	print $cmd,"\n";
#	process_cmd($cmd);
#}

if (defined $medical) {
        my $cmd = "$Rscript $Bin/plot_coverage_distribution_v2.r  --input1 $plus.stat --input2 $minus.stat --output $dOut/$fKey.png -C $Bin/medical.chr.list -w $window -x \"$x_title\" -y \"$y_title\" -t \"$title\"";
        print $cmd,"\n";
        process_cmd($cmd);
}else {
        $nchro=($chrnum>20)?10:$chrnum  unless(defined $nchro);
        my $cmd = "$Rscript $Bin/plot_coverage_distribution_v2.r  --input1 $plus.stat --input2 $minus.stat --output $dOut/$fKey.png  -w $window -n $nchro -x \"$x_title\" -y \"$y_title\" -t \"$title\"";
        print $cmd,"\n";
        process_cmd($cmd);
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub sort_by_chr{  ###sort for medical project
        my ($input,$output)=@_;
        my $dir=dirname $output;
        `mkdir -p $dir` if(!-d $dir);
        `grep '^[1-9]' $input|sort -n -k 1 -k 2 >$output`;
        `grep '^X' $input|sort -n -k 2 >>$output`;
        `grep '^Y' $input|sort -n -k 2 >>$output`;
}


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
			$name=10000000000;
		}

		if (exists $chrname{$name}{$chr}) {
			$chrname{$name}{$chr}.="$chr\t$inform\n";
			$cover{$chr}++;
		}
		else {
			$chrnum++;
			$chrname{$name}{$chr}="$chr\t$inform\n";
			$cover{$chr}=1;
			last if $chrnum>100;
		}
	}
	close(IN);
	#my $mark=0;
	if ($chrnum<=100) {
		open (OUT,">$output") or die $!;
		foreach my $name (sort {$a<=>$b} keys  %chrname) {
			foreach my $chr (sort  keys  %{$chrname{$name}}) {
				print OUT $chrname{$name}{$chr};
			}
		}
		close(OUT);
	}
	else {
		print "no sorting \n";
		system"cp $input $output";
		#my @chrlen=sort {$cover{$b}<=>$cover{$a}} keys  %cover;
		my @chrlen=sort {$b<=>$a} values  %cover;

	#	foreach my $chr (sort {$cover{$b}<=>$cover{$a}} keys  %cover) {
	#		$mark=$cover{$chr} if ($mark==0) ;
	#		foreach my $name (sort {$a<=>$b} keys  %chrname) {
	#			print OUT $chrname{$name}{$chr} if (exists $chrname{$name}{$chr}) ;
	#		}
	#	}
		print "$chrlen[0]\n";
		$window=500 if ($chrlen[0]<10000000);
		$window=100  if ($chrlen[0]<5000000);
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
  -queue 
  -medical 
USAGE
	print $usage;
	exit;
}
