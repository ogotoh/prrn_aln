##############################################################################
#
#	Extract information about exon-intron structure from a extended
#	fasta protein/DNA/RNA sequence file
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
#############################################################################

package	Exin;
require	Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(exin);

use lib "$ENV{HOME}/perl";	# perl module
use SspTab;			# species specific table
use strict;

sub exin {
    my ($file, $Copt, $apos) = @_;
    return unless (open(FILE, $file));
    my $cmp = 0;
    my @exons = ();
    my @exlen = ();
    my %emssg = ();
    my $tcode = 0;
    my $mcode = 0;
    my @ret = ();

    my ($gene, $sgn, $rngl, $rngr, $o15);
    while (<FILE>) {
	if (/^(>|\;D)/) {
	    chop();
	    ($gene, $sgn, $rngl, $rngr, $o15) = &dbentry($_);
	} elsif (/^;C /) {
	    $cmp += s/complement//g;
	    my @a = split(/\D+/);
	    shift(@a);		# Remove leading blank
	    push(@exons, @a);
	} elsif (/^;M (\S+) (\d+)/) {
	    my $fs = ($1 eq "Insert")? "D": "I";
	    $fs .= $2;
	    ++$emssg{$fs};
	} elsif (/^\s*\d/) {
	    my @a = split;
	    $mcode++ if ($a[3] == 0);	# initiation
	    if ($a[5] < 0) {
		if ($a[3] == 0) {
		    ++$emssg{"M-"};	# missed initiation
		} elsif ($mcode) {
		    ++$emssg{"T-"};	# missed termination
		}
	    }
	} elsif (/^\s*\w/) {
	    $tcode += tr/*ouxOUX/*ouxOUX/;	# termination code
	}
    }
    close(FILE);
    unless ($gene) {
	my @a = split(/\//, $file);
	$gene = pop(@a);
    }
    my $gname = $gene;
    if ($Copt =~ /^d/) {
	return unless ($tcode || keys(%emssg));
	push(@ret, $tcode . 'O ') if ($tcode);
	foreach (sort keys(%emssg)) {
	    push(@ret, "$emssg{$_}$_");
	}
	return @ret;
    } elsif ($Copt =~ /^g/) {
	my $noexon = @exons / 2;
	@ret = ($file, $gene);
	if ($o15 eq 'Y') {
	    if ($sgn eq '-') {my $tmp = $rngl; $rngl = $rngr; $rngr = $tmp;}
	    push(@ret, $rngl, $rngr, $sgn, $noexon);
	} else {
	    push(@ret, shift(@exons), pop(@exons), $cmp? "-": "+", $noexon);
	}
	push(@ret, $o15);
	return (@ret) unless ($tcode || keys(%emssg));
	push(@ret, $tcode . 'O') if ($tcode);
	foreach (sort keys(%emssg)) {
	    push(@ret, "$emssg{$_}$_");
	}
	return @ret;
    } elsif ($Copt =~ /^b/) {
	push(@ret, $gname);
	push(@ret, $rngl, $rngr);
	push(@ret, $cmp? "-": "+", $o15);
	my $parity = 1;
	while (my $bound = shift(@exons)) {
	    push(@ret, $bound - $parity);
	    $parity = 1 - $parity;
	}
	return @ret;
    } elsif ($Copt =~ /^B/) {
	push(@ret, $gname, $cmp? "-": "+");
	my $base = $exons[0];
	my $parity = 0;
	while (my $bound = shift(@exons)) {
	    push(@ret, $bound - $base + $parity);
	    $parity = 1 - $parity;
	}
	return @ret;
    }
    my $repe = $Copt =~ /^e/;
    my $repi = $Copt =~ /^i/;
    return unless ($repe || $repi || defined($apos));
    my $glen = $exons[-1] - $exons[0] + 1;
    push(@ret, $file);
    if (!($repe || $repi || $apos)) {	# 2nd and 2nd last
	my $lpos = $exons[2] - $exons[0];
	my $rpos = $exons[-1] - $exons[-3];
	push(@ret, $cmp? $rpos: $lpos, $cmp? $lpos: $rpos)
	    if (@exons > 2);
	return (@ret);
    }
    $apos *= 3 if ($apos);
    my ($slen, $nbound) = (0, 0);
    if ($cmp == 1) {
	my $gbegin = $exons[-1];
	while (my $rbound = pop(@exons)) {
	    my $lbound = pop(@exons);
	    push(@exlen, $rbound - $lbound + 1) if ($repe);
	    push(@exlen, $nbound - $rbound - 1) if ($nbound && $repi);
	    $nbound = $lbound;
	    my $plen = $slen;
	    $slen += $rbound - $lbound + 1;
	    if ($apos && $plen <= $apos && $apos < $slen) {
		my $gpos = $gbegin - $rbound + $apos - $plen;
	        push(@ret, $gpos);
		push(@ret, $glen - $gpos);
	    }
	}
    } else {
	my $gbegin = $exons[0];
	while (my $lbound = shift(@exons)) {
	    my $rbound = shift(@exons);
	    push(@exlen, $rbound - $lbound + 1) if ($repe);
	    push(@exlen, $lbound - $nbound - 1) if ($nbound && $repi);
	    $nbound = $rbound;
	    my $plen = $slen;
	    $slen += $rbound - $lbound + 1;
	    if ($apos && $plen <= $apos && $apos < $slen) {
		my $gpos = $lbound - $gbegin + $apos - $plen;
		push(@ret, $gpos);
		push(@ret, $glen - $gpos);
	    }
	}
    }
    push(@ret, @exlen) if ($repe || $repi);
    return @ret;
}

1;

