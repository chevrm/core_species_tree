#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Data::Dumper;
use Cwd 'abs_path';
use Bio::TreeIO;

## Setup script path
my $script_dir = abs_path($0);
$script_dir =~ s/\/core_species_tree\.pl//;

## Params
our $cpu = 8; ## default=8
our $bootstraps = 100; ## default=100

## Paths
our $tigrdir = "$script_dir/tigrfam";
our $tigrdb = $tigrdir . '/genprop0799.hmmdb';
our $tigrcut = $tigrdir . '/genprop0799.cutoffs.tsv';
our $raxml = '/home/mchevrette/builds/raxml-8.2.4/raxmlHPC';
our $astral = '/home/mchevrette/builds/Astral/astral.4.10.12.jar';

my $skip = 0; ## SET TO 1 to skip to raxml and astral

## Call genes
foreach my $fna (@ARGV){
    prodigal($fna) unless($skip == 1);
    
}

## Make leaf map
my %leaf2num = ();
my %num2leaf = ();
my $n = 1;
foreach my $faa (glob("*.prod.faa")){
    my $pref = $faa;
    $pref =~ s/\.prod\.faa$//;
    $leaf2num{$pref} = $n;
    $num2leaf{$n} = $pref;
    $n += 1;
}

## Read in all cutoffs
my %cut = ();
open my $tch, '<', $tigrcut or die $!;
while(<$tch>){
	unless($_ =~ m/^#/){
		my ($model, $tc, $nc) = split(/\t/, $_);
		$cut{$model} = {
			'trust'	=> $tc,
			'noise'	=> $nc
		};
	}
}
close $tch;

## Read in all hmm lengths
my %hmmlen = ();
open my $tdbh, '<', $tigrdb or die $!;
my $cur = '';
while(<$tdbh>){
	chomp;
	if($_ =~ m/^NAME\s+(\S+)/){
		$cur = $1;
	}elsif($_ =~ m/^LENG\s+(\d+)/){
		$hmmlen{$cur} = $1;
	}
}
close $tdbh;

## Run HMMs and assemble ML
unless($skip == 1){
    my %famseq = ();
    my %multilocus = ();
    foreach my $faa (glob("*.prod.faa")){
	my ($fsr, $mlr) = runhmm($faa, \%famseq, \%multilocus);
	%famseq = %$fsr;
	%multilocus = %$mlr;
    }
    
    ## Align all genes (protein)
    alignall(\%famseq);
}
    
## Convert alignments to nucl and save as phylip
foreach my $afa ( glob("*.afa") ){
    nucphy($afa) unless($skip == 1);
}

## RaxML all
foreach my $phy ( glob("*.phy") ){
    my $t = $phy;
    $t =~ s/\.phy$//;
    raxmlall($phy) unless(-e "RAxML_bootstrap.$t.raxml");
}

## ASTRAL
open my $bfh, '>', 'bs-files' or die $!;
foreach my $bs (glob("RAxML_bootstrap.*.fs.raxml")){
    print $bfh "$bs\n";
}
close $bfh;
system("cat *.fs.tre > allfam.tre");
system("java -jar $astral -i allfam.tre -b bs-files -r $bootstraps -o astral.species.out 2> /dev/null");

## Translate leaves
my @aso = ();
open my $afh, '<', 'astral.species.out' or die $!;
while(<$afh>){
    chomp;
    push @aso, $_;
}
close $afh;
my $bsto = new Bio::TreeIO(-file=>'>astral.species.bs', -format=>'newick');
my $consto = new Bio::TreeIO(-file=>'>astral.species.consensus', -format=>'newick');
my $treto = new Bio::TreeIO(-file=>'>astral.species.tre', -format=>'newick');
for(my $i=0;$i<scalar(@aso);$i+=1){
    open my $tfh, '>', 'tmp.tree' or die $!;
    print $tfh "$aso[$i]\n";
    close $tfh;
    my $tf = new Bio::TreeIO(-file=>'tmp.tree', -format=>'newick');
    my $tree = $tf->next_tree;
    foreach my $leaf ($tree->get_leaf_nodes){
	$leaf->id($num2leaf{$leaf->id});
    }
    if($i < $bootstraps){
	$bsto->write_tree($tree);
    }elsif($i == $bootstraps){
	$consto->write_tree($tree);
    }elsif($i == $bootstraps + 1){
	$treto->write_tree($tree);
    }else{
	die "Indexing error: i=$i\n";
    }
}

## Cleanup
my @torm = ('tmp*', '*.reduced', 'RAxML_*', 'bs-files', '*.fs.*', '*.prod.*', 'allfam.tre');
foreach (@torm){
    #system("rm $_");
}




sub prodigal{
    my ($fna, $a) = (shift, shift);
    my $pref = '';
    if($fna =~ m/(.+)\.f(n)?a(sta)?$/){
	$pref = $1;
	if($pref =~ m/\//){
	    my @parr = split(/\//, $pref);
	    $pref = $parr[-1];
	}
    }else{
	die "\"$fna\" does not have proper extension (.fasta, .fna, or .fa)\n";
    }
    $pref =~ s/_/-/g;
    print STDERR "$fna\tCalling genes with prodigal...";
    system("prodigal -t $pref.prod.train -c -i $fna > /dev/null 2>&1") unless(-e "$pref.prod.train");
    system("prodigal -c -i $fna -a $pref.prod.faa -d $pref.prod.orf -t $pref.prod.train > /dev/null 2>&1") unless(-e "$pref.prod.orf");
    print STDERR "DONE!\n";
}

sub runhmm{
    my ($faa, $fsr, $mlr) = (shift, shift, shift);
    my %fs = %$fsr;
    my %ml = %$mlr; 
    print STDERR "$faa\tScanning with TIGRFAM HMMs...";
    ## HMMscan of all prodigal faas
    my $pref = $faa;
    $pref =~ s/\.prod\.faa//;
    system("hmmscan -o tmp.hmmscan.out --tblout tmp.hmmtbl.out -E 1e-5 --cpu $cpu --noali $tigrdb $faa"); ## COMMENT OUT TO SKIP
    ## Parse hmmtbl
    open my $sfh, '<', "tmp.hmmtbl.out" or die "Died in hmmscanner: $!";
    my %ann = ();
    my %bit = ();
    while(<$sfh>){
	unless($_ =~ m/^#/ || $_ =~ m/^\W/){
	    chomp;
	    my ($hit, $mod, $query, $dash, $evalue, $score, @rest) = split(/\s+/, $_);
	    $bit{$query}{$hit} = $score;
	}
    }
    close $sfh;
    ## Assign annotations
    my %grab = ();
    #die Dumper(%bit) . "\n"; 
    foreach my $q (sort keys %bit){
	foreach my $h (sort keys %{$bit{$q}}){
	    if( $bit{$q}{$h} >= $cut{$h}{'trust'}){
		if(exists $ann{$h}{'q'}){
		    if($ann{$h}{'b'} < $bit{$q}{$h}){
			$ann{$h}{'q'} = $q;
			$ann{$h}{'b'} = $bit{$q}{$h};
			$grab{$q} = 1;
		    }
		}else{
		    $ann{$h}{'q'} = $q;
		    $ann{$h}{'b'} = $bit{$q}{$h};
		    $grab{$q} = 1;
		}
	    }
	}
    }
    ## Grab protein seqs
    my %s = ();
    my $fa = new Bio::SeqIO(-file=>$faa, -format=>'fasta');
    while(my $seq = $fa->next_seq){
	if(exists $grab{$seq->id}){
	    my $ss = $seq->seq;
	    $ss =~ s/\*//;
	    $s{$seq->id} = $ss;
	}
    }
    ## Assemble ML seq 
    my $mls = '';
    foreach my $tf (sort keys %cut){
	if(exists $ann{$tf}{'q'}){
	    $mls .= $s{$ann{$tf}{'q'}};
	    $fs{$tf}{$pref}{'seq'} = $s{$ann{$tf}{'q'}};
	    $fs{$tf}{$pref}{'orf'} = $ann{$tf}{'q'};
	}else{
	    foreach(1..$hmmlen{$tf}){
		$mls .= 'X';
	    }
	}
    }
    $ml{$pref} = $mls;
    print "DONE!\n";
    return(\%fs, \%ml);
}

sub alignall{
    my $fsr = shift;
    my %fs = %$fsr;
    my $f = scalar(keys %fs);
    foreach my $t (keys %fs){
	print STDERR "$t\tAligning individual gene family...";
	open my $fsh, '>', "$t.fs.faa" or die $!;
	foreach my $g (sort keys %{$fs{$t}}){
	    print $fsh '>' . join('--', $fs{$t}{$g}{'orf'}, $g) . "\n" . $fs{$t}{$g}{'seq'} . "\n";
	}
	close $fsh;
	system("mafft --quiet --namelength 70 $t.fs.faa > $t.fs.afa");
	print "DONE!\n";
    }
}

sub nucphy{
    my $afa = shift;
    my $phy = $afa;
    $phy =~ s/\.afa/\.phy/;
    my $a = new Bio::SeqIO(-file=>$afa, -format=>'fasta');
    my $p = Bio::AlignIO->new(-format =>'phylip', -file => ">$phy");
    my $newaln = new Bio::SimpleAlign;
    while(my $seq = $a->next_seq){
	my ($orf, $genome) = split(/--/, $seq->id);
	open my $gfh, '<', "$genome.prod.orf" or die "$!\t$genome\n";
	my $flag = 0;
	my $nuc = '';
	while(<$gfh>){
	    chomp;
	    if($_ =~ m/^>(.+)/){
		my $id = $1;
		$id =~ s/^(\S+).+/$1/;
		if($id eq $orf){
		    $flag = 1;
		}else{
		    $flag = 0;
		}
	    }elsif($flag == 1){
		$nuc .= $_;
	    }
	}
	close $gfh;
	my @codon = $nuc =~ /\w{3}/g;
	my $nucaln = '';
	my $i = 0;
	foreach my $aa (split //, $seq->seq){
	    if($aa eq '-'){
		$nucaln .= '---';
	    }else{
		$nucaln .= $codon[$i];
		$i += 1;
	    }
	}
	my $newseq = new Bio::LocatableSeq(-seq=>$nucaln, -id=>$leaf2num{$genome}, -start=> 1, -end=>length($nuc)-3);
	$newaln->add_seq($newseq);
    }
    $p->write_aln($newaln);
}

sub raxmlall{
    my $phy = shift;
    my $pref = $phy;
    $pref =~ s/\.phy//;
    print STDERR "$phy\traxml...";
    system("$raxml -m GTRGAMMA -n $pref.raxml -s $phy -f a -x 897543 -N $bootstraps -p 345232 --silent > /dev/null");
    system("$raxml -f b -m GTRGAMMA -z RAxML_bootstrap.$pref.raxml -t RAxML_bestTree.$pref.raxml -n $pref.BS_TREE > /dev/null");
    system("mv RAxML_bipartitions.$pref.BS_TREE $pref.tre") if(-e "RAxML_bipartitions.$pref.BS_TREE");
    print STDERR "DONE!\n";
}
