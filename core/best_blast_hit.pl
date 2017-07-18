#!/usr/local/bin/perl

use strict;
use warnings;

die ("Usage: $0 <blastin >bestout") unless (scalar(@ARGV)==0);

my (%ID_TO_NDX);
my (@BEST_LINE);
my (@BEST_EXPECT);
my (@BEST_SCORE);

print STDERR "Reading...\n";
&read_inputs();
print STDERR "Printing...\n";
&print_array();
print STDERR "Done.\n";

sub read_inputs () {
    my ($line);
    my ($accession,$expect,$score,$index);
    my ($prev_score,$prev_expect);
    my ($base)=0;
    while (<STDIN>) {
	chomp;
	$line = $_;
	if ($line =~ /^(\S+).*(\S+)\s+(\S+)$/) {
	    $accession=$1;
	    $expect=0+$2;
	    $score=0+$3;
	    $index=$ID_TO_NDX{$accession};
	    if (!defined($index)) {
		# This is a new accession
		$index=$base;
		$ID_TO_NDX{$accession}=$index;
		$prev_expect=100;
		$prev_score=0;
		$base++;
	    } else {
		# This is a repeat accession
		$prev_expect = $BEST_EXPECT[$index];
		$prev_score = $BEST_SCORE[$index];
	    }
	    if ($expect < $prev_expect ||
		($expect == $prev_expect && $score > $prev_score)) {
		$BEST_LINE[$index]=$line;
		$BEST_EXPECT[$index]=$expect;
		$BEST_SCORE[$index]=$score;
	    }
	} else {
	    print STDERR "WARNING: Could not parse!\n>>$line<<\n";
	}
    }
}
sub print_array () {
    my ($X);
    foreach $X (@BEST_LINE) {
	print "$X\n";
    }
}
