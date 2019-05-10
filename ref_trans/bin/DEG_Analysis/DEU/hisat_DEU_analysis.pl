# Writer:         Mengf <mengf@biomarker.com.cn>
# Program Date:   2012.
# Modifier:       Mengf <mengf@biomarker.com.cn>
# Modifier:       lium <lium@biomarker.com.cn>
# Last Modified:  2012-7-28
my $ver="1.0.0";

use strict;
use Cwd;
use Getopt::Long;
my $BEGIN=time();
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="LncRNA";  
my $version="v2.0.3";  
my %config=%{readconf("$Bin/../../../config/db_file.cfg")};
$config{qsub_sh}="/share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh"	if(!exists $config{qsub_sh});$config{python}="/share/nas2/genome/biosoft/Python/2.7.8/bin/python"		if(!exists $config{python});
$config{Queue_type}="medical.q"			if(!exists $config{Queue_type});

my $programe_dir=basename($0);
my $path=dirname($0);
######################����д����֮ǰ��һ��д��ʱ�䡢������;������˵����ÿ���޸ĳ���ʱ��Ҳ������ע�͹���


our %opts;
GetOptions(\%opts,"cfg=s","od=s","h","hisat=s","gtfdir=s");


if(!defined($opts{cfg}) || !defined($opts{od}) || !defined($opts{hisat}) || defined($opts{h}) || !defined($opts{gtfdir})){
	&help();
	exit;
}


###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $programe_dir Time :[$Time_Start]\n\n";




my $Rscript=$config{Rscript};
my $cfg=&ABSOLUTE_DIR($opts{cfg});

&MKDIR($opts{od});
&MKDIR("$opts{od}/work_sh");
my $outdir=&ABSOLUTE_DIR($opts{od});
my $sh_dir=&ABSOLUTE_DIR("$opts{od}/work_sh");
my $hisat=&ABSOLUTE_DIR($opts{hisat});
my $notename=`hostname`;chomp $notename;
my %DEG_type;
&LOAD_PARA($cfg,\%DEG_type);

###### read config file into hash
my %CONTROL = ();
read_config_into_hash($cfg, \%CONTROL);
my $gtf_dir=&ABSOLUTE_DIR($opts{gtfdir});
my $finish=(glob"$gtf_dir/*gtf")[0];

open (SH,">$sh_dir/DEU_analysis.sh") or die $!;
control_DE_analysis($outdir, \%CONTROL);


###### do DE analysis according to control
sub control_DE_analysis{
	# get para
	my ($outdir, $control) = @_;
	# do compare with replicates
	if( exists $$control{Sep} ) {
		my $sep_groups = $$control{Sep};
		my $len = @$sep_groups;
		my $gtfdir=&ABSOLUTE_DIR($opts{gtfdir});
		my $gtf=glob("$gtfdir/*.gtf");
		my $deu;
		my $deu=&read_gtf_for($gtf);
		if ($deu==1){
        		print "deu analysis continue\n";
		}else {
		        print "the $gtf do not contain the feature \"exon\"\n";
		        exit;
		}
		my $hisat=$opts{hisat};
		my $strand;
		if ($$control{Lib_type} eq "fr-firststrand"){
			$strand="yes";
		}
		elsif ($$control{Lib_type} eq "fr-secondstrand"){
			$strand="reverse";
		}else{
			$strand="no";
		}
	
		for(my $i=0; $i<$len; $i++) {
			compare_with_replicates($$sep_groups[$i], $outdir, $gtf, $hisat,$strand);
		}
	}
	close SH;
	&Shell_qsub ("$sh_dir/DEU_analysis.sh",$DEG_type{'Queue_type'},10);
}

sub read_gtf_for {
	my ($gtf,$deu)=@_;
	my $dir=dirname($gtf);
	my @tem=readpipe( "cut -f 3 $gtf|sort|uniq ");
	foreach my $v (@tem){
		chomp $v;
		$v=~s/\s$//;
		if ($v eq 'exon'){
			$deu=1;
			last;
		}
	}
	return $deu;	
}

###### reading config file into hash
sub read_config_into_hash{
	# get config file name
	my ($in, $control) = @_;

	# open config file
	open(IN,"$in") || die "open file $in is failed: $!";

	# reading every line
	while( my $line = <IN> ) {
		chomp $line;
		next if($line!~/\w/); 	# ignore the empty line
		next if($line=~/^#/); 	# ignore the comment line

		# get parameters according the first word
		my @str = $line=~/(\S+)/g;
        next unless ($str[0] eq "Com" or $str[0] eq "Sep" or $str[0] eq "fold" or $str[0] eq "FDR" or $str[0] eq "Ref_ann" or $str[0] eq "Lib_type");

		# check the number of word in one line
		my $len = @str;
		if( $len != 2 ){ print "$in: $line: the number of world != 2\n"; exit; }

		# store the group info
		if( exists $$control{$str[0]} ){
			# check sample analysis
			if( ($str[0] ne "Com") && ($str[0] ne "Sep") ) {
				print "$in: repeat program word $str[0]\n"; exit;
			}
			# push
			my $tmp = $$control{$str[0]};
			push(@$tmp, $str[1]);
		}else{
			if( ($str[0] eq "Com") || ($str[0] eq "Sep") ) {
				my @tmp = ();
				$tmp[0] = $str[1];
				$$control{$str[0]} = \@tmp;
			} else {
				$$control{$str[0]} = $str[1];
			}
		}
	}

	# check key
	#if( ! exists $$control{fold} ) { die "ERROR: key fold not exists in file $in";}
	#if( ! exists $$control{FDR} ) { die "ERROR: key FDR not exists in file $in";}
	#if( (! exists $$control{Com}) && (! exists $$control{Sep}) ) {
	#	die "ERROR: key Com and Sep not exists in file $in";
	#}

	# close file handle
	close IN;
}


# here will do one compare which contain all levels
sub compare_with_replicates {
	# get para
	my ($sample_id, $outdir,$gtf,$hisat,$strand) = @_;
################create the path.txt used for DEU_analysis
	my @bams=glob("$hisat/*/*HISAT_aln.sorted.bam");
        mkdir "$outdir/htseq_count"unless -d "$outdir/htseq_count";
	`cp $gtf $outdir/htseq_count`;
        open CONFIG,">$outdir/htseq_count/bam_to_exon_count.config" ||die $!;
        foreach  (@bams) {
                print CONFIG "$_\n";
        }
        close CONFIG;
	$DEG_type{'Queue_type'}="medical.q"	if(!exists $DEG_type{'Queue_type'});
	print "$Rscript $Bin/bin/bam_count.R  qsub=$config{qsub_sh} python=$config{python} infile=$outdir/htseq_count/bam_to_exon_count.config gtf=$gtf type=bam od=$outdir/htseq_count queue=$DEG_type{'Queue_type'} strand=$strand";
        `$Rscript $Bin/bin/bam_count.R  qsub=$config{qsub_sh} python=$config{python} infile=$outdir/htseq_count/bam_to_exon_count.config gtf=$gtf type=bam od=$outdir/htseq_count queue=$DEG_type{'Queue_type'} strand=$strand`;
	$sample_id =~ s/\s//g;
	my @groups = split(';', $sample_id);
	if($#groups < 1) {die "ERROR: the number of group < 2\n"; }
	my @sample = ();
	for(my $i=0; $i<=$#groups; $i++) { 
		my @tmp = split(',', $groups[$i]);
		$sample[$i] = \@tmp; 
	}

	# combinate compare
	#
	for (my $i=0;$i<=$#groups-1;$i++) {
		for (my $j=$i+1;$j<=$#groups;$j++) {
			my $v1 = $sample[$i];
			my $v2 = $sample[$j];
			my $v1_name = join("_", @$v1);
			my $v2_name = join("_", @$v2);
			# get output prefix
			my $vs_out = "$v1_name" . "_vs_" . "$v2_name";
			# create dir for output
			&MKDIR("$outdir/$vs_out");

			my $cmd="$Rscript $Bin/bin/DEU_DEXSeq.R  infile=$outdir/htseq_count/path.txt group=$vs_out od=$outdir/$vs_out  ";
                        print SH "$cmd\n";
		}
	}
}







###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd $programe_dir Time :[$Time_End]\n\n";
&Runtime($BEGIN);






sub Shell_qsub
{ # &Shell_qsub($sh,$qeue,$cpu);
	my $sh = shift;
	my $qeue = shift;
	my $cpu = shift;

	if ($notename=~/login\-0\-4/)
	{
		`sh $config{qsub_sh} --queue $qeue --reqsub -maxproc $cpu --independent $sh `;
	}
	else
	{
		`sh $config{qsub_sh} --queue $qeue --reqsub --maxproc $cpu --independent $sh `;
	}
}




###########subs
sub LOAD_PARA {
	my $para_file= shift;
	my $para = shift;

	my $error_status = 0;
	open(IN, $para_file) || die "fail open: $para_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if (/^$/ or /^\#/) ;
		my ($para_key,$para_value) = split(/\s+/,$_);
		if ($para_key=~/_G/) {
			push @{$para->{$para_key}},$para_value;
			next;
		}
		$para->{$para_key} = $para_value;
		if (!defined $para_value) {
			warn "Non-exist: $para_key\n";
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error Undefined parameter\n" if ($error_status) ;
}

sub parse_config
{ # load config file
	my $config_file= shift;
	my $DataBase= shift;
	
	my $error_status = 0;
	open(IN,$config_file) || die "fail open: $config_file";
	while (<IN>) {
		chomp;
		s/^\s+//;s/\s+$//;s/\r$//;
		next if(/^$/ or /^\#/);
		my ($software_name,$software_address) = split(/\s+/,$_);
		$DataBase->{$software_name} = $software_address;
		if (! -e $software_address){
			warn "Non-exist:  $software_name  $software_address\n"; 
			$error_status = 1;
		}
	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
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
#############################################################################################################
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

# &show_log1("cmd")
sub show_log1()
{
    open LOG, ">$outdir/../../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}

# &show_log2("cmd")
sub show_log2()
{
    open LOG, ">>$outdir/../../bioinfor_pipeline.log";
	my ($txt) = @_ ;
    my $time = time();
    my $Time = &sub_format_datetime(localtime($time));
    print LOG "$Time:\t$txt\n" ;
    return ($time) ;
	close LOG;
}
#############################################################################################################
sub help{
print <<"Usage End.";
Description: Differencial Express Analysis pipeline;
Version: $ver

Usage:
-i                All_gene_counts.list                            must be given;
-cfg              soft parameter to DE miRNA Analysis             must be given;
-od               Out dir                                         must be given;
-hisat		  dir for DEU analysis                            must be given;
-gtfdir	          dir for ref_gtf  				  must be given; 
-h                help document
Usage End.
exit;
}
