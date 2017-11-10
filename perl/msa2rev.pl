#!/usr/bin/perl

##############################################################################
#
#	rename last-refined MSA to REV
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
#____________________________________________________________________________
#
#	Usage: msa2rev.pl [-fFractionalPlace] [-rRename] MSA.*
#
#############################################################################

my $list;
my $nsb = 10;
my $nbdy = "REV";
my $Debug;

while ($_ = $ARGV[0], /^-/) {
	shift;
	/^-f(\d+)/ && ($nsb = $1);
	/^-r(\S+)/ && ($nbdy = $1);
	/^-D/ && ($Debug = 1);
}

foreach my $arg (@ARGV) {
	my ($head, $num, @sub) = split('\.', $arg);
	my $bdy = $head . '.' . $num;
	splice(@sub, $nsb);
	my $sbn = join('.', @sub);
	if (undef($list{$bdy}) || $sbn > $list{$bdy}) {
	    $list{$bdy} = $sbn;
	}
}

foreach my $sq (sort keys(%list)) {
	my $sqn = $sq;
	my @a = split('\.', $sq, 2);
	my $rev = $nbdy. '.' . $a[1];
	if ($list{$sq}) {
	    $sqn .= '.' . $list{$sq};
	    if ($Debug) {print "mv $sqn $rev\n";}
	    else {rename($sqn, $rev);}
	} else {
	    if ($Debuf) {print "link $sqn $rev\n";}
#	    else {link($sqn, $rev);}
	    system("cp $sqn $rev");
	}
}
