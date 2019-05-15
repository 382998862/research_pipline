#!/usr/bin/perl

use strict;
use LWP::Simple;




print STDERR "\n
###############################################################################################
#
# This is the miRDeep2 installer. 
# It will work under a bash,csh and ksh shell.
# It will try to download all necessary files and install them. 
# Please restart your shell to make changes take effect
#
###############################################################################################


";

my $time =time;
my $new = $ARGV[0];

if($new !~ /no/){
    print STDERR "making backup of .bashrc,.bash_profile and .cshrc and removing all entries of mirdeep in those files\n";
    rem_mirdeep(".bashrc");
    rem_mirdeep(".bash_profile");
    rem_mirdeep(".cshrc");
}

my $gcc=`gcc --version 2>&1`;
if($gcc !~ /(GCC)/i){
    die "\nError:\n\tno gcc compiler installed. Please install a gcc compiler\n";
}else{
    if($gcc =~ /^gcc\s*\S*\s*(\d+\S+)\s*/){
        print STDERR "gcc version: $1                                      already installed, nothing to do ...\n";
    }
} 


my $wget=`wgetx`;
my $curl=`curl 2>&1`;

my $dtool='';
my $dopt='';

if($wget =~ /URL/i){
	$dtool ="wget";
	
}elsif($curl){
	$dtool ="curl";
	$dopt=" -O";
}else{
	die "No commandline download tool found on your system. Please install wget or curl on your machine\n";
}

my $dir=`pwd 2>&1`;

chomp $dir;

my $err;

if(not -d 'mirdeep2' and $dir !~ /mirdeep2/){
    $err=system('tar --no-same-permissions -xvvzf mirdeep2.tar.gz');
    if($err){
        die "tar xvvzf mirdeep2.tar.gz could not be executed\n\n";
    }
}elsif($dir =~ /mirdeep2/i){
    chdir("..");
    $dir=`pwd 2>&1`;
    chomp $dir;
}

my $dfile='';

##only attach to config file if not yet existing
my $in=`grep  "$dir/mirdeep2:*" ~/.bashrc`;
if(not $in){
    `echo 'export PATH=\$PATH:$dir/mirdeep2' >> ~/.bashrc`;
}

$in=`grep  "$dir/mirdeep2:*" ~/.bash_profile`;
if(not $in){
    `echo 'export PATH=\$PATH:$dir/mirdeep2' >> ~/.bash_profile`;
}

my $in2=`grep  "$dir/mirdeep2:*" ~/.cshrc`;
if(not $in2){
    #`echo 'setenv PATH \$PATH:$dir/mirdeep2' >> ~/.cshrc`;
}


chdir("mirdeep2");
if(not -d "essentials"){
    `mkdir essentials`;
}

chdir("essentials");

my $a=`uname -a`;

my $bowtie;

my $bowtie_version="0.12.7";

my $binst=`bowtie --version 2>&1`;

if($binst =~ /version\s*(\S*)/i){
    print STDERR "Bowtie version: $1                                      already installed, nothing to do ...\n";
}else{
    if(not -d "bowtie-$bowtie_version"){

        print STDERR "Downloading bowtie $bowtie_version binaries\n\n";
        if($a =~ /Darwin/i){ ## download mac version
            $bowtie = "bowtie-$bowtie_version-macos-10.5-x86_64.zip";
        }elsif($a =~ /x86_64/i){
            $bowtie = "bowtie-$bowtie_version-linux-x86_64.zip";
        }else{
            $bowtie = "bowtie-$bowtie_version-linux-i386.zip";
        }
        
        if(not -f $bowtie){
			 
			if(check("http://ovh.dl.sourceforge.net/project/bowtie-bio/bowtie/$bowtie_version/$bowtie")){
				$err=system("$dtool http://ovh.dl.sourceforge.net/project/bowtie-bio/bowtie/$bowtie_version/$bowtie $dopt");
				if($err){
					die "\nError:\n\t$bowtie could not be downloaded\n\n\n";
				}
            }elsif(check("http://ovh.dl.sourceforge.net/project/bowtie-bio/bowtie/old/$bowtie_version/$bowtie")){
				$err=system("$dtool http://garr.dl.sourceforge.net/project/bowtie-bio/bowtie/old/$bowtie_version/$bowtie $dopt");
				if($err){
					die "\nError:\n\t$bowtie could not be downloaded\n\n\n";
				}
			}else{
				die "\nError:\n\t$bowtie not found on server http://garr.dl.sourceforge.net/project/bowtie-bio/bowtie/ \n\n\n";
			}
        }
        
        if(not -f "$bowtie"){
            die "$bowtie download failed \n";
        }
        
        
        print STDERR "Installing bowtie binaries\n\n";
        $err=system("unzip $bowtie");
        
        if($err){
            die "unzip $bowtie was not successful\n";
        }
    }

    $in = `grep  "$dir/mirdeep2/essentials/bowtie-$bowtie_version:*" ~/.bashrc`;
    if(not $in){
        `echo 'export PATH=\$PATH:$dir/mirdeep2/essentials/bowtie-$bowtie_version' >> ~/.bashrc`;
    }

    $in = `grep  "$dir/mirdeep2/essentials/bowtie-$bowtie_version:*" ~/.bash_profile`;
    if(not $in){
        `echo 'export PATH=\$PATH:$dir/mirdeep2/essentials/bowtie-$bowtie_version' >> ~/.bash_profile`;
    }


    $in2 = `grep  "$dir/mirdeep2/essentials/bowtie-$bowtie_version:*" ~/.cshrc`;

    if(not $in2){
        #`echo 'setenv PATH \$PATH:$dir/mirdeep2/essentials/bowtie-$bowtie_version' >> ~/.cshrc`;
    }

}

my $ret = checkBIN("RNAfold -h 2","usage");
#my $rna_inst=`RNAfold 2>&1`;

#if($rna_inst !~ /no\s*RNAfold/i){
if($ret == 0){
   print STDERR "RNAfold                                                   already installed, nothing to do ...\n";
}else{
    if(not -d "ViennaRNA-1.8.4"){
        $dfile="ViennaRNA-1.8.4.tar.gz";
        if(not -f $dfile){
            print STDERR "Downloading Vienna package now\n\n";
			if(check("http://www.tbi.univie.ac.at/~ivo/RNA/ViennaRNA-1.8.4.tar.gz")){
				$err=system("$dtool http://www.tbi.univie.ac.at/~ivo/RNA/ViennaRNA-1.8.4.tar.gz $dopt");
				if($err){
					die "Download of Vienna package not successful\n\n";
				}
			}else{
				die "Vienna package not found at  http://www.tbi.univie.ac.at/~ivo/RNA/ViennaRNA-1.8.4.tar.gz\n";
			}
    }
        
        
        if(not -f "ViennaRNA-1.8.4.tar.gz"){
            die "Vienna package download failed\n";
        }
        
        print STDERR "Installing Vienna package now \n\n";
        `tar xvvzf ViennaRNA-1.8.4.tar.gz`;
        chdir("ViennaRNA-1.8.4");
        `./configure --prefix=$dir/mirdeep2/essentials/ViennaRNA-1.8.4/install_dir`;
        `make`;
        `make install`;
        
        
        chdir("..");
    }
}
$in = `grep "$dir/mirdeep2/essentials/ViennaRNA-1.8.4/install_dir/bin:*" ~/.bashrc`;
if(not $in){
    print STDERR "Vienna package path has been added to \$PATH variable\n"; 
    `echo 'export PATH=\$PATH:$dir/mirdeep2/essentials/ViennaRNA-1.8.4/install_dir/bin' >> ~/.bashrc`;
}

$in = `grep "$dir/mirdeep2/essentials/ViennaRNA-1.8.4/install_dir/bin:*" ~/.bash_profile`;
if(not $in){
    print STDERR "Vienna package path has been added to \$PATH variable\n"; 
    `echo 'export PATH=\$PATH:$dir/mirdeep2/essentials/ViennaRNA-1.8.4/install_dir/bin' >> ~/.bash_profile`;
}



$in2 = `grep  "$dir/mirdeep2/essentials/ViennaRNA-1.8.4/install_dir/bin:*" ~/.cshrc`;
if(not $in2){
    #`echo 'setenv PATH \$PATH:$dir/mirdeep2/essentials/ViennaRNA-1.8.4/install_dir/bin' >> ~/.cshrc`;
}

$ret = checkBIN("randfold","let7");

#my $randf = `randfold -h`;

#if($randf =~ /no\s*randfold/i){ ## this should work
if($ret == 0){
    print STDERR "randfold\t\t\t\talready installed, nothing to do ...\n";
}else{

	$dfile="squid-1.9g.tar.gz";
    if(not -f $dfile){
        print STDERR "Downloading SQUID library now\n\n";
        `$dtool ftp://selab.janelia.org/pub/software/squid/squid-1.9g.tar.gz $dopt`;
    }

    if(not -f "squid-1.9g.tar.gz"){
        die "squid could not be downloaded\n";
    }

    if(not -d "squid-1.9g"){
        print STDERR "Extracting squid and configuring it now\n\n"; 
        `tar xxvzf squid-1.9g.tar.gz`;
        chdir("squid-1.9g");
        `./configure`;
        `make`;
        chdir("..");
    }

	$dfile="randfold-2.0.tar.gz";
    if(not -f $dfile ){
        print STDERR "Downloading randfold now\n\n";
        `$dtool http://bioinformatics.psb.ugent.be/supplementary_data/erbon/nov2003/downloads/randfold-2.0.tar.gz $dopt`;
    }


    if(not -f "randfold-2.0.tar.gz"){
        die "randfold could not be downloaded\n\n";
    }

    if(not -d "randfold-2.0"){
        print STDERR "Installing randfold now\n\n";
        `tar xvvzf randfold-2.0.tar.gz`;
        
        chdir("randfold-2.0");

        open IN,"<Makefile" or die "File Makefile not found\n\n";
        
        open OUT,">Makefile_new" or die "Makefile_new could not be created\n\n";

        while(<IN>){
            if(/INCLUDE=-I\.\s*/i){
                print OUT "INCLUDE=-I. -I$dir/mirdeep2/essentials/squid-1.9g -L$dir/mirdeep2/essentials/squid-1.9g/\n";
            }else{
                print OUT;
            }
        }
        close IN;
        close OUT;
        
        `mv Makefile Makefile.orig`;
        `mv Makefile_new Makefile`;

        `make`;

        chdir("..");
    }

    $in = `grep  "$dir/mirdeep2/essentials/randfold-2.0:*" ~/.bashrc`;
    if(not $in){ 
        print STDERR "Randfold path has been added to \$PATH variable\n"; 
        `echo 'export PATH=\$PATH:$dir/mirdeep2/essentials/randfold-2.0' >> ~/.bashrc`;
    }
    
    $in = `grep  "$dir/mirdeep2/essentials/randfold-2.0:*" ~/.bashrc_profile`;
    if(not $in){ 
        print STDERR "Randfold path has been added to \$PATH variable\n"; 
        `echo 'export PATH=\$PATH:$dir/mirdeep2/essentials/randfold-2.0' >> ~/.bash_profile`;
}


    $in2 = `grep  "$dir/mirdeep2/essentials/randfold-2.0:*" ~/.cshrc`;
    if($in2){
        #`echo 'setenv PATH \$PATH:$dir/mirdeep2/essentials/randfold-2.0' >> ~/.cshrc`;
    }

}


##check for zlib perl
my $zlib=`perl -e 'use Compress::Zlib;' 2>&1`;

if(not $zlib){
    print STDERR "Compress::Zlib\t\t\t\talready installed, nothing to do ...\n";
}else{
    die "please install Compress::Zlib by using CPAN before you proceed\n";

}

#my $pdfapi=`perl -e 'use PDF::API2;' 2>&1`;

$ret = checkBIN("perl -e \'use PDF::API2; print \"installed\";\'","installed");


#if(not $pdfapi){
if($ret == 0){
    print STDERR "PDF::API2                                           already installed, nothing to do ...\n";
}else{


	$dfile="PDF-API2-0.73.tar.gz";
    if(not -f $dfile){
        print STDERR "Downloading PDF-API2 now\n\n";
        `$dtool http://search.cpan.org/CPAN/authors/id/A/AR/AREIBENS/PDF-API2-0.73.tar.gz $dopt`;
    }

    if(not -f "PDF-API2-0.73.tar.gz"){
        die "Download of PDF-API2 failed\n\n";
    }

    #if(not -d "PDF-API2-0.73"){
        print STDERR "Installing PDF-API2 now\n\n";
        `tar xvvzf PDF-API2-0.73.tar.gz`; 
        chdir("PDF-API2-0.73");
        
        `perl Makefile.PL PREFIX=$dir/mirdeep2/lib LIB=$dir/mirdeep2/lib`;
        
        `make`;
        `make install`;
    #}
    chdir("..");
}

#my $perl = "$], $/";
#chomp $perl;
#$perl =~ s/0/\./g;
#$perl =~ s/\.+/\./g;
#$perl =~ s/\,+//g;

#$in = `grep "$dir/mirdeep2/lib/perl5/site_perl/$perl:*" ~/.bashrc`;

#if(not $in){ 
   print STDERR "PDF::API2 package path has been added to \$PERL5LIB variable\n"; 
    `echo 'export PERL5LIB=\$PERL5LIB:$dir/mirdeep2/lib/' >> ~/.bashrc`;
#}

#$in = `grep "$dir/mirdeep2/lib/perl5/site_perl/$perl:*" ~/.bash_profile`;

#if(not $in){
    print STDERR "PDF::API2 package path has been added to \$PERL5LIB variable\n"; 
    `echo 'export PERL5LIB=\$PERL5LIB:$dir/mirdeep2/lib/' >> ~/.bash_profile`;
#}



#$in2 = `grep  "$dir/mirdeep2/lib/perl5/site_perl/$perl:*" ~/.cshrc`;

#if(not $in2){
    #`echo 'setenv PERL5LIB \$PERL5LIB:$dir/mirdeep2/lib/' >> ~/.cshrc`;
#}

print STDERR "\n\nInstallation successful\n\n\nPlease start a new shell\n\n\n\n";



sub rem_mirdeep{
    my ($file)=@_;
    if(-f "$ENV{'HOME'}/$file"){
        `cp ~/$file ~/${file}_$time`;
        print STDERR "~/$file backup is ~/${file}_$time\n";
        `mv ~/$file ~/$file.bak`;
        open OUT,">$ENV{'HOME'}/$file" or die "Cannot create file $ENV{'HOME'}/$file\n";
        open IN,"$ENV{'HOME'}/$file.bak";
        my @line;
        my $tmp;
        while(<IN>){
            $tmp="";
            if(/mirdeep2/){
                @line = split(/:/);
                foreach(@line){
                    if(/mirdeep2/){}else{
                        $tmp.="$_";
                    }
                }
                if($file !~ /.cshrc/){
                    
                    if($tmp !~ /PATH=\$PATH$/ and $tmp !~ /PERL5LIB=\$PERL5LIB\s*$/){
                        print OUT "$tmp\n";
                    }
                }else{
                    if($tmp !~ /PATH \$PATH$/ and $tmp !~ /PERL5LIB \$PERL5LIB\s*$/){
                        print OUT "$tmp\n";
                    }
                }
            }else{
                print OUT;
            }
        }
    }else{
        print "file $ENV{'HOME'}/$file.bak not found\n";
    }
}


sub checkBIN{
    my ($a,$b) = @_;    
	my $e = system("$a> tmp 2>tmp2");
	
    open IN,"<tmp";
    my $found =1;
    while(<IN>){
		if(/$b/){
			$found =0;
			
		}
	}
	close IN;
	if($found){
		open IN,"<tmp2";
		while(<IN>){
			if(/$b/){
				$found =0;
				
			}
		}
	}
    close IN;
    return $found;
}


## this routine checks if files exists on remote servers
sub check{
	my ($url) = @_;
	if (head($url)) {
		return 1;
	}else{
		return 0;
	}
}
