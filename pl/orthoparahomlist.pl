#!/usr/bin/perl
# Compare gene sets of two species with BLAST and compile lists with
# orthologs, paralogs and other homologs
# Mario Stanke, 6.6.2011

use strict;
use warnings;

use Getopt::Long; # for parameter specification on the command line
use Bio::Tools::Run::StandAloneBlast;

my $usage = <<'ENDUSAGE';
orthoparahomlist.pl       compare gene sets of two species with BLAST and
                          compile lists with orthologs, paralogs and other homologs

SYNOPSIS

orthoparahomlist.pl --protA=proteinsA.fa --protB=proteinsB.fa 

    protA and protB specify a set of genes from two different species A and B, respectively.
    Every protein sequence is compared to every other protein sequence with BLASTP. This 
    script requires an installation of NCBI blastall and that databases have been created with
    fomatdb.
    
OPTIONS

    --help             output this help message
    --ortho=ortho.lst  Create an output file ortho.lst that contains for each gene a in A the
                       putative ortholog b in B as inferred by reciprocal best BLAST hit. List
                       only genes which have a reciprocal best hit. Format:
                       geneA.id <TAB> geneB.id
    --para=para.lst    Create an output file para.lst that contains for each gene a in A all other
                       genes in A that meet the E-value threshold. Include the E-value of the hit 
                       and the percentage of identical amino acids of the first HSP. Format:
                       geneA.id1 <TAB> geneB.id1 <TAB> E-value <TAB> percent identical
                       geneA.id1 <TAB> geneB.id2 <TAB> E-value <TAB> percent identical
                       geneA.id1 <TAB> geneB.id3 <TAB> E-value <TAB> percent identical
                       ...
                       geneA.id1 <TAB> geneB.idM <TAB> E-value <TAB> percent identical
                       geneA.id2 <TAB> geneB.id1 <TAB> E-value <TAB> percent identical
                       geneA.id2 <TAB> geneB.id2 <TAB> E-value <TAB> percent identical
                       ...
    --hom=hom.lst      Create an output file para.lst that contains for each gene a in A all
                       genes in B that meet the E-value threshold. Format like para.lst.

    --E                minimal E-value threshold for BLAST search
    --delBlastOutput   delete BLAST output after finishing

DESCRIPTION
      
  Example:

    orthoparahomlist.pl -E 1e-5 --protA=protA.fa --protB=protB.fa --ortho=ortho.lst --para=para.lst --hom=hom.lst

ENDUSAGE

my ($protAfname, $protBfname, $orthofname, $parafname, $homfname, $help, $delBlastOutput); # file names and options
my $E = 1e-10; # default E-value
sub run_blast;

GetOptions('protA=s'=>\$protAfname,
	   'protB=s'=>\$protBfname,
	   'ortho=s'=>\$orthofname,
	   'para=s'=>\$parafname,
	   'hom=s'=>\$homfname,
	   'E:f' =>\$E,
           'help!'=>\$help,
           'delBlastOutput!'=>\$delBlastOutput);

if ($help){
    print $usage;
    exit(1);
}

if (!defined($protAfname) || !defined($protBfname)){
    print "Missing input filename. Must specify both protA and protB.\n$usage";
    exit(1);
}

# check whether the two input files exist
if (! -f "$protAfname"){
    print "Input file $protAfname does not exist. Please check.\n";
    exit(1);
}

if (! -f "$protBfname"){
    print "Input file $protBfname does not exist. Please check.\n";
}

# check whether BLAST databases have been created for the two input files
if (! -f "$protAfname.pin" || ! -f "$protBfname.pin"){
    print "Must create BLAST database first, e.g. with\nformatdb -i $protAfname\nand/or\nformatdb -i $protBfname\n";
    exit(1);
}

if (!defined($orthofname) && !defined($parafname) && !defined($homfname)){
    print "You have not specified any of the three output file names. Will do nothing.\n";
    exit(0);
}

# run all 3 required BLAST combinations (A versus A, A versus B, B versus A)
my $resAvsA;
$resAvsA = run_blast($protAfname, $protAfname, "AvsA.blast.out") if ($parafname); # only required if paralogs are requested
my $resAvsB = run_blast($protBfname, $protAfname, "AvsB.blast.out");
my $resBvsA = run_blast($protAfname, $protBfname, "BvsA.blast.out");

# hash that stores for each gene in B the best hit in A (for reciprocity)
my %bestHitB = (); # keys: seq ids in B, values: id of best hit in A (if any)

# retrieve best hits of genes in B for testing reciprocity
while(my $result = $resBvsA->next_result()){
    my $firsthit = $result->next_hit();
    if ($firsthit){
	$bestHitB{$result->query_name()} = $firsthit->name();
    }
}

# open output files
if ($orthofname){
    open (ORTHO,">$orthofname") or die ("Could not open $orthofname for writing.\n");
}
if ($parafname){
    open (PARA,">$parafname") or die ("Could not open $parafname for writing.\n");
}
if ($homfname){
    open (HOM,">$homfname") or die ("Could not open $homfname for writing.\n");
}

# go through A and output orthologs and homologs in B
while(my $result = $resAvsB->next_result()){
    my $hit = $result->next_hit();
    if ($hit && defined($bestHitB{$hit->name()}) && $bestHitB{$hit->name()} eq $result->query_name()){ # reciprocal best hit
	print ORTHO $result->query_name(). "\t" . $hit->name() . "\n";
    }
    if ($homfname){
	while ($hit && (my $firsthsp = $hit->next_hsp())){
	    print HOM $result->query_name(). "\t" . $hit->name() . "\t". $hit->significance . "\t" . $firsthsp->percent_identity() . "\n";
	    $hit = $result->next_hit();
	}
    }
}

if ($parafname){
    # go through A and output paralogs
    while(my $result = $resAvsA->next_result()){
	while (my $hit = $result->next_hit()){
	    if ($result->query_name() ne $hit->name()){
		my $firsthsp = $hit->next_hsp();
		if ($firsthsp){
		    print PARA $result->query_name(). "\t" . $hit->name() . "\t". $hit->significance . "\t" . $firsthsp->percent_identity() . "\n";
		}
	    }
	}
    }
}

# delete BLAST output files if requested
if ($delBlastOutput){
    foreach my $fname ("AvsA.blast.out", "AvsB.blast.out", "BvsA.blast.out"){
	unlink $fname if (-f $fname);
    }
}

# run_blast
# 
# parameters: database name, query file name, output file name
# returns a reference to a BLAST object
sub run_blast {
    my ($database, $queryfname, $outfname) = @_;
    
    if (-f $outfname) { # BLAST done before, retrieve results from file
	return Bio::SearchIO->new( -format => "blast", -file   => $outfname);
    } else {  # run BLAST
	my $blastfactory = Bio::Tools::Run::StandAloneBlast->new( -program => "blastp",
								  -database => $database,
								  -expect => $E,
								  -outfile => $outfname);
	my $resultio = $blastfactory->blastall($queryfname);
	return $resultio; 
    }
}
