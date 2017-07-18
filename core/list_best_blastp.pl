#!/usr/local/bin/perl

use strict;
use warnings;

die ("Usage: $0 <deflines1> <deflines2> <best>") unless (scalar(@ARGV)==3);

my ($file_primary_fasta)=$ARGV[0];
my ($file_secondary_fasta)=$ARGV[1];
my ($file_best_hits)=$ARGV[2];

my (@PRIMARY_ID);
my (%PRIMARY_TO_NAME);
my (%SECONDARY_TO_NAME);
my (%PRIMARY_TO_SECONDARY);

print STDERR "Read primary...\n";
&read_primary();
print STDERR "Read secondary...\n";
&read_secondary();
print STDERR "Read best...\n";
&read_pairs();
print STDERR "Output...\n";
&print_pairs();
print STDERR "Done\n";

# Assumptions.
# Primary has lines like this. GenBank ID. No greater-than and no pipe characters.
# XP_019538164.1,"Description of protein"
# Secondary has lines like this. VectorBase ID. No suffix like -PA. 
# AGAP005712,"Description of protein"
# Protein descriptions may have leading spaces to indicate missing fields (unknown proteins).
# Pairs file is ordered: genbank vectorbase

sub read_primary () {
    my ($line);
    my ($accession,$descrip);
    my ($count)=0;
    open (PRIMARY, "<$file_primary_fasta") or die ("Cannot open primary");
    while (<PRIMARY>) {
	chomp;
	$line = $_;
	# Extract entire description including internal and leading blanks.
	if ($line =~ /(\S+),"(.+)"/) {
	    $accession=$1;
	    $descrip=$2;	    
	    $PRIMARY_ID[$count++]=$accession;
	    $PRIMARY_TO_NAME{$accession}=$descrip;
	}	
    }
    close (PRIMARY);
}
sub read_secondary () {
    my ($line);
    my ($accession,$descrip);
    open (SECONDARY, "<$file_secondary_fasta") or die ("Cannot open secondary");
    while (<SECONDARY>) {
	chomp;
	$line = $_;
	# Extract entire description including internal and leading blanks.
        if ($line =~ /^(\S+),"(.+)"/) {
	    $accession=$1;
	    $descrip=$2;
	    $accession = &strip_suffix($accession);
	    $SECONDARY_TO_NAME{$accession}=$descrip;
	    #print STDERR "DEBUG $line\nDEBUG ($accession) $descrip\n\n";
	}
    }
    close (SECONDARY);
}
sub read_pairs () {
    my ($line);
    my ($aid,$bid);
    my ($adesc,$bdesc);
    open (ALIGN, "<$file_best_hits") or die ("Cannot open best");
    while (<ALIGN>) {
	chomp;
	$line = $_;
	undef($bdesc);
	undef($adesc);
	# Assume lines are like this.
	# AGAP000002-PA XP_019525794.1
        if ($line =~ /^(\S+)\s(\S+)/) {
	    $aid=$1;
	    $bid=$2;
	    $aid=&strip_suffix($aid);
	    $bid=&strip_suffix($bid);
	    # Look up the description string for each ID
	    $adesc=$PRIMARY_TO_NAME{$aid};
	    $bdesc=$SECONDARY_TO_NAME{$bid};
	} else {
	    die ("Cannot parse pair: $line");
	}
	die ("Unrecognized ID b($bid) from line $line") unless (defined($bdesc));
	die ("Unrecognized ID a($aid) from line $line") unless (defined($adesc));
	$PRIMARY_TO_SECONDARY{$aid}=$bid;
    }
    close (ALIGN);
}
sub print_pairs () {
    my ($index,$acc,$adesc,$best,$bdesc);
    # Use the input index rather than the keys of the hash.
    # This ensures that output order is same as input order.
    for ($index=0; $index < scalar(@PRIMARY_ID); $index++) {
	$acc = $PRIMARY_ID[$index];
	$adesc = $PRIMARY_TO_NAME{$acc};
	$best = $PRIMARY_TO_SECONDARY{$acc};
	if (defined($best)) {
	    $bdesc=$SECONDARY_TO_NAME{$best};
	} else {
	    $best = "NONE";
	    $bdesc = "NONE";
	}
	die ("Undefined descrip for $acc") unless (defined($adesc));
	die ("Undefined best for $acc") unless (defined($best));
	die ("Undefined best descrip for $acc") unless (defined($bdesc));
	print STDOUT "${acc},\"${adesc}\",${best},\"${bdesc}\"\n";
    }
}
sub strip_suffix() {
    my($str)=shift;
    my ($pos);
    $pos=index($str,"-");	    
    $str=substr($str,0,$pos) if ($pos>0);
    return ($str);
}
