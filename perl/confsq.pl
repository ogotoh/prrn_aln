#!/usr/bin/perl

##############################################################################
#
#	confirm consistency between an MSA and its components
#
#       Osamu Gotoh, ph.D.      (-2001)
#	Saitama Cancer Center Research Institute
#	818 Komuro, Ina-machi, Saitama 362-0806, Japan
#
#	Osamu Gotoh, Ph.D.      (2001-)
#	National Institute of Advanced Industrial Science and Technology
#	Computational Biology Research Center (CBRC)
#	2-41-6 Aomi, Koutou-ku, Tokyo 135-0064, Japan
#
#	Osamu Gotoh, Ph.D.      (2003-2012)
#	Department of Intelligence Science and Technology
#	Graduate School of Informatics, Kyoto University
#	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
#
#	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
#
#___________________________________________________________________________
#
#	Usage: confsq.pl MSA dir/
#
#	MSA: multiple sequence alignment in prrn native format
#	dir: directory that contains aa sequences 
#
##############################################################################

use strict;

while ($_ = $ARGV[0], /^-/) {
	shift;
	if (/^-\?/ || /^-h/) {&usage;}
}

my $msa;
my $pas;
my $temp = "/tmp/confsq_tmp$$";

foreach my $arg (@ARGV) {
	$msa = $arg if (-s $arg && -f $arg);
	$pas = $arg if (-d $arg);
}

$pas = $pas . "/" if ($pas && substr($pas, -1, 1) ne "/");

open(CAT, "utp -l1 $msa |") || die "Can't run utp!\n";

while (<CAT>) {
	my ($snum, $snam) = split;
	$snam = $pas . $snam;
	unless (-e $snam) {
	    printf("%-39s\tnot found!\n", $snam);
	    next;
	}
	exit 1 if system("rdn -csd $msa $snum > $temp");
	my @msg = split(' ', `iden -n $temp $snam`);
	if (!@msg || $msg[2]) {
	    system("iden -w50 $temp $snam");
	} else {
	    printf("%-39s\tis OK\n", $snam);
	}
}
continue {
	unlink($temp);
}
close(CAT);
exit(0);

sub usage {
	print STDERR "Usage: confsq.pl MSA dir/ \n";
	exit(1);
}
