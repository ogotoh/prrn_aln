#!/usr/bin/perl

##############################################################################
#
#	Construct GSA-MPSA from aa sequences in Src using 'upg' as the guide tree
#
#
#	Osamu Gotoh, ph.D.      (-2001)
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
#	Usage:	ls *.upg.* | upg2prrn.pl -sSrc -L [-O]
#	or	upg2prrn.pl -sSrc [-O] *.upg.*
#
#	'upg' represents a UPGMSA tree file constructed by upg command
#	-O option forces overwrite the output MSA file if it exists
#############################################################################

my $opt = "-yl2 -KP";
my $lst;
my $vbs;
my $ovw;
my $pas;

sub System {
	my $cmd = shift;
	print STDERR "$cmd\n" if ($vbs & 1);
	system($cmd) == 0 or print STDERR "$cmd failed !\n";
}

while ($_ = $ARGV[0], /^-/) {
    shift;
    if (/^-L/) {$lst = 1;}
    elsif (/^-V/) {$vbs = 1;}
    elsif (/^-O/) {$ovw = 1;}
    elsif (/^(-s\S+)/) {$pas = $1;}
    else {$opt .= " $_";}
}

sub a2p {
    $_ = shift;
    my @a = split('\.', $_);
    my $n = pop(@a);
    my $tgt = "MSA.$n";
    next if (-s $tgt && !$ovw);
    my $cmd = "prrn5 -b$_ $opt -o$tgt";
    $cmd .= " $pas" if ($pas);
    &System($cmd);
}

if ($lst) {
    while (<>) {
	chop;
	&a2p($_);
    }
} else {
    while ($_ = $ARGV[0]) {
	shift;
	&a2p($_);
    }
}

