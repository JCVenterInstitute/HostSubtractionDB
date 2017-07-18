#! /usr/local/bin/perl -w
use strict;
use Getopt::Std;

my %options=();
getopts("v", \%options);
my ($NAMEFILE) = $ARGV[0];
die ("Usage: $0 [-v] <nameFile> <inFile >outFile") unless (defined($NAMEFILE));

print STDERR "Reading STDIN, Writing STDOUT\n";
print STDERR "Read FASTQ and write FASTQ\n";
print STDERR "Read names from = $NAMEFILE\n";
my ($INVERSE) = defined($options{v});
if ($INVERSE) {
    print STDERR "Inverse operation (exclude named reads).\n";
} else {
    print STDERR "Normal operation (include named reads).\n";
}

my ($PREFIX_SEQID) = '@';
my ($PREFIX_QVID) = '+';
my ($PREFIX_FASTA) = '>';
my (%NAMESINCORE);
my ($ins,$outs);
my ($defline);
my ($readname);
my ($sequence);
my ($quality);

sub slurp () {
    open (NAMES, "<$NAMEFILE") or die ("Cannot read nameFile");
    my ($name);
    while (<NAMES>) {
	chomp;
	$name=$_;
	$NAMESINCORE{$name}=1;
    }
    close NAMES;
}

sub process () {

    my ($line);
    my ($state);
    my ($prefix);
    
    $state = 0;
    $ins=0;
    $outs=0;
    while (<STDIN>) {
	chomp;
	$line = $_;
	if (++$state == 5) {$state=1;}
	if ($state == 1) {
	    $prefix = substr ($line, 0, 1);
	    die ("Expected $PREFIX_SEQID, got $line") unless ($prefix eq $PREFIX_SEQID);
	    $defline = substr ($line, 1);
	    ($readname) = split (/\s+/, $defline, 2);
	} elsif ($state == 2) {
	    $sequence=$line;
	} elsif ($state == 3) {
	    $prefix = substr ($line, 0, 1);
	    die ("Expected $PREFIX_QVID, got $line") unless ($prefix eq $PREFIX_QVID);
	} elsif ($state == 4) {
	    $quality = $line;
	    &output();
	} else {
	    die ("Unexpected: state=$state");
	}
    }
    die ("Nothing read from STDIN") if ($state==0);
    die ("Unexpected termainal state: $state") unless ($state==4);
}

sub output {
    ++$ins;
    if ($INVERSE && !defined($NAMESINCORE{$readname}) ||
	!$INVERSE && defined($NAMESINCORE{$readname})) {
	&printout();
	++$outs;
    }
}


sub printout () {
    # my ($hold_len) = length($sequence);
    print STDOUT "${PREFIX_SEQID}${defline}\n";
    print STDOUT "$sequence\n";
    print STDOUT "${PREFIX_QVID}\n";
    print STDOUT "$quality\n";
}

print STDERR "Loading names...\n";
slurp();
print STDERR "Loaded " . scalar(keys(%NAMESINCORE)) . " unique names.\n";
print STDERR "Filtering...\n";
process();
print STDERR "Read in $ins, wrote out $outs.\n";
print STDERR "Done\n";

