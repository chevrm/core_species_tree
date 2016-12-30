#!/bin/env perl

use strict;
use warnings;

## Check to make sure genpro tsv present
my ($gp, $gpsite) = ('genprop0799.tsv', 'http://www.jcvi.org/cgi-bin/genome-properties/GenomePropDefinition.cgi?prop_acc=GenProp0799');
die "ERROR:  please create $gp from $gpsite\n" unless(-e $gp);

## Download and unpack the full tigr database
unless(-e 'TIGRFAMs_14.0_HMM.tar.gz'){
	system("wget ftp://ftp.tigr.org/pub/data/TIGRFAMs/14.0_Release/TIGRFAMs_14.0_HMM.tar.gz");
	system("tar zxvf TIGRFAMs_14.0_HMM.tar.gz");
	system("rm -r tigrfam") if(-d 'tigrfam');
	system("mkdir tigrfam");
	system("mv TIGR*.HMM tigrfam");
}

## Download and unpack the info
unless(-e 'TIGRFAMs_14.0_INFO.tar.gz'){
	system("wget ftp://ftp.tigr.org/pub/data/TIGRFAMs/14.0_Release/TIGRFAMs_14.0_INFO.tar.gz");
	system("tar zxvf TIGRFAMs_14.0_INFO.tar.gz");
	system("rm -r info") if(-d 'info');
	system("mkdir info");
	system("mv TIGR*.INFO info");
}

## Read in the genprop list, create cutoff tsv, press hmmdb
system("rm -r genprop0799.hmmdb*") if(-e 'genprop0799.hmmdb');
open my $gph, '<', $gp or die $!;
open my $ioh, '>', 'genprop0799.cutoffs.tsv' or die $!;
print $ioh join("\t", '#HMM', 'Trust_Cutoff', 'Noise_Cutoff') . "\n";
while(<$gph>){
	chomp;
	if($_ =~ m/.+(TIGR\d+).+/){
		my $t = $1;
		system("cat tigrfam/$t.HMM >> genprop0799.hmmdb");
		open my $ifh, '<', "info/$t.INFO" or die $!;
		my ($tc, $nc) = (0,0);
		while(<$ifh>){
			chomp;
			if($_ =~ m/^TC\s+(\d+\.\d+)/){
				$tc = $1;
			}elsif($_ =~ m/^NC\s+(\d+\.\d+)/){
				$nc = $1;
				last;
			}
		}
		close $ifh;
		print $ioh join("\t", $t, $tc, $nc) . "\n";
	}
}
close $gph;
close $ioh;
system("hmmpress genprop0799.hmmdb");
