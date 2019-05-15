#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;
use Cwd 'getcwd';

################################# Prioritization ##########################################

################################## USAGE ##################################################
my $usage=
"Usage:

$0 -g go_annotation_file -v validated.csv -p predicted.csv -q

-h      print this usage
-g      file containing GO annotation
-v      validated miRNA-target file in degradome data
-p      predicted miRNA-target file using TAPIR or other tools
-q      Quiet mode .. no log/progress information to STDERR
";

###########################################################################################

################################### INPUT #################################################
my %options=();
getopts('g:v:p:hq',\%options);

###########################################################################################

############################## GLOBAL VARIABLES ###########################################

#hashes
my %hash_prediction;
my %hash_validate;
my %hash_map_go;

#other variables;
my $cur_path=getcwd(__FILE__);
###########################################################################

############################## MAIN #######################################################
#Check options and file format
check_options();
check_file_format();
unless($options{'q'}){
    print STDERR "###############Running the prioritization step###############\n\n";
}
#Read the validated file from degradome analysis
unless($options{'q'}){
    print STDERR "read the validated MTI file...\n";
}
read_validate($options{'v'});
#Read the prediction file
unless($options{'q'}){
    print STDERR "read the prediction MTI file...\n";
}
read_prediction($options{'p'});

#parse the gene_association file
unless($options{'q'}){
    print STDERR "parse the gene_association file...\n";
}
read_gene_asso($options{'g'});

#Prepare for the function similarity computing
unless($options{'q'}){
    print STDERR "Prepare for the function similarity computing...\n";
}
my $file_for_fs="pre_for_FS";
my $file_for_fs_r="pre_for_FS.R";
gene2go($file_for_fs);
pre_cal_fs_R($file_for_fs_r);

#invoke R for calculating function similarity
unless($options{'q'}){
    print STDERR "invoke R for calculating function similarity...\n";
}
system("R <$file_for_fs_r --vanilla 1>rlog.txt 2>&1");
clean();
unless($options{'q'}){
    print STDERR "Prioritization finished\n\n";
}
exit;

##############################################################################################

############################## SUBROUTINES ###################################################
sub check_options{
    if($options{'h'}) {
        die "$usage";
    }
    unless($options{'g'} and $options{'v'} and $options{'p'}){
        die "FATAL: please specify all files\n$usage";
    }
    unless(-f $options{'g'}){
        die "FATAL: can not find file $options{'g'}\n$usage";
    }
    unless(-f $options{'v'}){
        die "FATAL: can not find file $options{'v'}\n$usage";
    }
    unless(-f $options{'p'}){
        die "FATAL: can not find file $options{'p'}\n$usage";
    }
}
sub check_file_format{
    open IN,$options{'v'};
    my $vline=<IN>;
    close IN;
    open IN,$options{'p'};
    my $pline=<IN>;
    close IN;
    unless($vline =~ /^\S+\s+\S+/){
        die "Please make sure the first and second column of $options{'v'} should be miRNA and target, and separated by space or TAB\n";
    }
    unless($pline =~ /^\S+\s+\S+/){
        die "Please make sure the first and second column of $options{'p'} should be miRNA and target, and separated by space or TAB\n";
    }
    open IN,$options{'g'};
    my $gline=<IN>;
    close IN;
    unless($gline =~ /^\S+\s+GO:\d+/){
        die "The format of each line should be 'Gene/Protein <tab/space> GO:1234567 <tab/space> GO:1234567'\n";
    }
}
sub read_prediction{
    my $file=shift;
    open PRE,$file or die("Can not open $file");
    while(<PRE>){
           chomp;
           next if /^#/;
           my($mir,$gene)=(split /\s+/)[0,1];
           $gene=~s/\.\d+$//;
           $gene=~s/_T\d+$//;
           push(@{$hash_prediction{$mir}{$gene}},$_);
    }
    close PRE;
}
sub read_validate{
    my $file=shift;
    open VAL,$file or die("Can not open $file");
    while(<VAL>){
           chomp;
           next if /^#/;
           next if /^SiteID/;
           my($mir,$gene)=(split /\s+/)[1,2];
           $gene=~s/\.\d+$//;
           $gene=~s/_T\d+$//;
           push(@{$hash_validate{$mir}{$gene}},$_);
    }
    close VAL;
}

sub read_gene_asso{
    my $file=shift;
    open IN, $file or die("Can not open $file");
    while(<IN>){
           s/\r?\n$//;
           my @arr=split /\t/;
           my $gene=shift @arr;
           #next if $gene=~/\.\d+$/;
           $hash_map_go{$gene}=join(" ",@arr);
    }
    close IN;
}
sub gene2go{
    my $outfile=shift;
    open OUT,">$outfile";
    print OUT join("\t",qw(Mir Gene GO Validated))."\n";
    foreach my $mir(sort keys %hash_validate){
        my $exists_validate=0;
        next unless exists $hash_prediction{$mir};
        foreach my $gene(keys %{$hash_validate{$mir}}){
            delete $hash_prediction{$mir}{$gene} if exists $hash_prediction{$mir}{$gene};
            if(exists $hash_map_go{$gene}){
               my $gos=$hash_map_go{$gene};
               $exists_validate=1;
               print OUT "$mir\t$gene\t$gos\tV\n";
            }
            # else{
                 # print OUT "$mir\t$gene\t\tV\n";
            # }
        }
        next unless $exists_validate;
        foreach my $gene(keys %{$hash_prediction{$mir}}){
            if(exists $hash_map_go{$gene}){
               my $gos=$hash_map_go{$gene};
               print OUT "$mir\t$gene\t$gos\tP\n";
            }
            # else{
                 # print OUT "$mir\t$gene\t\tP\n";
            # }
        }
    }
    close OUT;
}
sub pre_cal_fs_R{
    my $outfile=shift;
    open OUTR,">$outfile" or die("Can not create $outfile.Please make sure you have the permission!");
    my $out_file=<<TMP;
library(csbl.go);
library(plyr);
library(reshape2);
setwd("$cur_path");
cal_prob_table <- function(ann.filename) {
    ent <- entities.from.text(ann.filename);
    # Similarity table
    prob.table <- entities.to.prob.table(ent);
    filename <- paste(ann.filename,"sim",sep=".");
    write.prob.table(prob.table, filename);
}
self.entity.sim<-function(x,metric,multiple="max"){
    validate<-x[x\$Validated=="V",];
    prediction<-x[x\$Validated=="P",];
    all_score_mf<-numeric();
    all_score_cc<-numeric();
    all_score_bp<-numeric();
    all_score<-numeric();
    len<-length(validate\$GO);
    for(i in prediction\$GO){
        score_mf<-numeric();
        score_bp<-numeric();
        score_cc<-numeric();
        for(j in validate\$GO){
            score_mf<-c(score_mf,entity.sim(i,j,"MF",metric,multiple));
            score_bp<-c(score_bp,entity.sim(i,j,"BP",metric,multiple));
            score_cc<-c(score_cc,entity.sim(i,j,"CC",metric,multiple));
        }
        score_mf[is.na(score_mf)]<-0;
        score_cc[is.na(score_cc)]<-0;
        score_bp[is.na(score_bp)]<-0;
        all_score_mf<-c(all_score_mf,mean(score_mf));
        all_score_cc<-c(all_score_cc,mean(score_cc));
        all_score_bp<-c(all_score_bp,mean(score_bp));
    }
    prediction<-transform(prediction,MF_score=all_score_mf,BP_score=all_score_bp,CC_score=all_score_cc,
                            all_score=all_score_mf+all_score_cc+all_score_bp
                          );
    prediction<-prediction[order(prediction\$all_score,decreasing=T),];
    return(prediction);
}

if(!file.exists("$options{'g'}.sim")){
    cal_prob_table("$options{'g'}");
}
set.prob.table(filename="$options{'g'}.sim");

raw_file<-read.table("$file_for_fs",header=T,sep="\\t",stringsAsFactors=F);
raw_file_list<-dlply(raw_file,.(raw_file[["Mir"]]));
result_list<-lapply(raw_file_list,function(x) self.entity.sim(x,"Resnik"));
out_result<-data.frame();
out_result<-ldply(result_list)[,-1];
write.table(out_result,"$options{'p'}.go",quote=F,sep="\\t",row.names=F);
TMP
   print OUTR $out_file;
   close OUTR;
}

sub clean{
    unlink $file_for_fs,$file_for_fs_r;
}