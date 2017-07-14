#! /usr/local/bin/perl -w
use strict;
use Getopt::Std;

print STDERR "Reading STDIN, Writing STDOUT\n";
print STDERR "Read FASTQ and write FASTQ\n";
die ("Usage: $0 <inFile >outFile") unless (scalar(@ARGV)==0);

#my %options=();
#getopts("v", \%options);
#my ($NAMEFILE) = $ARGV[0];
#die ("Usage: $0 <nameFile> <inFile >outFile") unless (defined($NAMEFILE));
#my ($INVERSE) = defined($options{v});
#if ($INVERSE) {
#    print STDERR "Inverse operation (exclude named reads).\n";
#} else {
#    print STDERR "Normal operation (include named reads).\n";
#}

my ($PREFIX_SEQID) = '@';
my ($PREFIX_QVID) = '+';
my ($PREFIX_FASTA) = '>';
my (%NAMESINCORE);
my ($ins,$outs);
my ($readname);
my ($sequence);
my ($quality);


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
	    $readname = $line;
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
    die ("Unexpected termainal state: $state") unless ($state==4);
}

sub output {
    ++$ins;
    if (!defined($NAMESINCORE{$readname})) {
	&printout();
	++$outs;
    }
    $NAMESINCORE{$readname}=1;	
}


sub printout () {
    print STDOUT "$readname\n";
    print STDOUT "$sequence\n";
    print STDOUT "${PREFIX_QVID}\n";
    print STDOUT "$quality\n";
}

print STDERR "Filtering...\n";
process();
my ($size) = scalar(keys(%NAMESINCORE));
print STDERR "Reads_input: $ins\n";
print STDERR "Unique_names:: $size\n";
print STDERR "Reads_output: $outs.\n";
print STDERR "Done\n";

