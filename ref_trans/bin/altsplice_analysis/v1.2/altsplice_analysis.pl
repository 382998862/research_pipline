#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use Cwd 'abs_path';
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my %config=%{readconf("$Bin/../../../config/db_file.cfg")}; 
my $USAGE = "

#*******************************************************************************************
Options:
-cufflinks   <dir>         Cufflinks OUT dir                               [required]
-od          <dir>         AS output dir                                   [default: (./)]
-queue       <str>         the queue is used for qsub jobs                 [default: general.q]
-Ref_Genome  <dir>         genome dir,*.hdrs and *.gtf should include      [required]
-step        <int>         The step of program                             [default 1]
              1             Extract AS
              2             Stat AS Result


#******************************************************************************************

";

my $extract_as="$Bin/ASprofile.b-1.0.4/extract-as";
my $summary_as="$Bin/ASprofile.b-1.0.4/summarize_as.pl";
my $as_fpkm="$Bin/ASprofile.b-1.0.4/extract-as-fpkm";
my $collect_fpkm="$Bin/ASprofile.b-1.0.4/collect_fpkm.pl";

my ($cufflinks, $fout, $step,$Ref_Genome,$queues);

GetOptions (  "cufflinks:s"  =>\$cufflinks, "Ref_Genome:s"   =>\$Ref_Genome, "queue:s"   =>\$queues, "od:s"   =>\$fout,    "step:n" =>\$step, ) || die $USAGE;

die $USAGE unless ($cufflinks and $Ref_Genome);

$fout ||= "./";
mkdir $fout unless (-d $fout);
$fout       = abs_path($fout);
$Ref_Genome = abs_path($Ref_Genome);
$cufflinks  = abs_path($cufflinks);
my $hdrs=(glob("$Ref_Genome/*.hdrs"))[0];
my $gtf =(glob("$Ref_Genome/*.gtf"))[0];

$step ||= 1;
$queues||="general.q";

=c
my $note_name = `hostname`;  chomp $note_name;
if ($note_name =~ /cluster/)
{
	print "This jobs cannot be down on cluster host!\n";
	exit;
}
=cut
my @cufflinks=glob("$cufflinks/*");


my $work_sh = "$fout/work_sh";       mkdir $work_sh unless (-d $work_sh);
my $Step_1_sh = "$work_sh/Step_1_Extract_AS.sh";
my $Step_2_sh = "$work_sh/Step_2_Extract_fpkm.sh";
my $Step_3_sh = "$work_sh/Step_3_fpkm_compare.sh";

### Step 1 : Extract Alt gene list from cuffcmp.tracking file

if ($step==1)
{
	open (OUT1, ">", $Step_1_sh) || die "Open $Step_1_sh failed!\n";
	open (OUT2, ">", $Step_2_sh) || die "Open $Step_2_sh failed!\n";
	foreach my $cufflink (sort @cufflinks) 
	{
		my $sam=basename$cufflink;
		next if ($sam!~/\w+/);
		mkdir "$fout/$sam" unless (-d "$fout/$sam");
		print OUT1 "cd $fout/$sam && perl $Bin/bin/changename.pl $cufflink/cuffcmp.StringTie_asm.gtf.tmap $cufflink/StringTie_asm.gtf $fout/$sam && ";
		print OUT1 "$extract_as transcripts.gtf $hdrs  >$sam.tmap.as && ";
		print OUT1 "perl $summary_as  transcripts.gtf $sam.tmap.as  -p $sam && \n";
		print OUT2 "cd $fout/$sam && $as_fpkm  transcripts.gtf $hdrs $sam.as.nr  -W 9 >$sam.W9.fpkm \n";
	}
	close OUT1;
	close OUT2;
	&qsubOrDie($Step_1_sh, $queues,20, "2G");
	#&Check_qsub_error($Step_1_sh);
	print "\n\nStep 1 :";
	print "\n\tExtract_AS is Done!\n\n";
	print "sh $Step_2_sh ";
	`sh $Step_2_sh `;
	print "\n\tExtract_fpkm is Done!\n\n";
	$step++;
}

### Step 2 : AS stat and draw
my %stat;
my %sample;
if ($step==2) 
{

	foreach my $cufflink (sort @cufflinks) 
	{
		my $sam=basename$cufflink;
		next if ($sam!~/\w+/);
		$sample{$sam}=1;
		my $stats="$fout/$sam/$sam.W9.fpkm";
		open (IN,"$stats")||die $!;
		while (<IN>) {
			chomp;
			next if ($_=~/event|_OFF/);
			my $type=(split(/\t+/,$_))[1];
			$type=~s/_ON|_OFF//g;
			$stat{$type}{$sam}++;
		}
        close(IN);
		`mv $fout/$sam/$sam.W9.fpkm $fout `;

	}
	
#	open (OUT, ">", $Step_3_sh) || die "Open $Step_3_sh failed!\n";
#	print OUT "cd $fout && perl $collect_fpkm ";
#	@samples=sort keys %sample;
#	print OUT join(",",@samples),"  -s W9.fpkm > All_sample.W9.diff-exp";
#	close OUT;
#	`sh $Step_3_sh `;
#	print "\n\nStep 3 :";
#	print "\n\t Sample Fpkm Compare is Done!\n\n";
	
	open (OUT, ">$fout/As_Result.stat") || die "Open $fout/As_Result.stat  failed!\n";
	print OUT "AS\t",join("\t",sort {$a cmp $b} keys %sample),"\n";
	foreach my $type (sort {$b cmp $a} keys %stat) {
		print OUT "$type";
		foreach  my $sam(sort {$a cmp $b} keys %sample) {
			print OUT  "\t$stat{$type}{$sam}";
		}
		print OUT "\n";
	}
	close OUT;
	my $cmd="$config{Rscript} $Bin/bin/plot_stat_AS_events.r -i $fout/As_Result.stat -o $fout/As_event_stat -b ";
	print $cmd,"\n";
	system$cmd;

	print "\n\nStep 4 :";
	print "\n\t Stat AS Result is Done!\n\n";

}
