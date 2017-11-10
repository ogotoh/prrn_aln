#!/usr/bin/perl

##############################################################################
#
#	Find transcript sequences that are derived from the same genomic region
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

sub usage {
	print STDERR "samegr.pl [-ON] [-Hthr] [-spas] seqs\n";
	exit 1;
}

# Examples:
# samegr.pl -U -O1 -r cyp/* famaln/prd/* > rentab.txt
# samegr.pl -U -O2 ../gnel2/* pas/* > rentab.txt
# samegr.pl -F * | clade

my $pas;
my $ovrthr = 1;
my $omode = 0;
my $summary = 0;
my $list;
my $uniq;
my $ngsid = 8;
my $lwrlmt = 400;
my $uprlmt = 800;
my $outform;
my $rev;

while ($_ = $ARGV[0], /^-/) {
	shift;
	last if /^-$/;
	/^-D/ && ($debug = 1);
	/^-F(\d*)/ && ($outform = $1? $1: 1);
	/^-H(\S+)/ && ($ovrthr = $1);
	/^-L/ && ($list = 1);
	/^-N(\d+)/ && ($ngsid = $1);
	/^-O(\S+)/ && ($omode = $1);
	/^-S/	&& ($summary = 1);
	/^-U/ && ($uniq = 1);
	/^-l(\d+)/ && ($lwrlmt = $1);
	/^-r/ && ($rev = 1);
	/^-s(\S+)/ && ($pas = $1);
	/^-u(\d+)/ && ($uprlmt = $1);
	/^\?/ && (&usage);
}

$pas .= '/' if ($pas && substr($pas, -1, 1) ne '/');

my @args = ();

if ($list) {
	$list = shift;
	open(LST, $list) or die "Can't open $list !\n";
	while (<LST>){
	    my @a = split;
	    push(@args, @a);
	}
	close(LST);
}
push(@args, @ARGV);

my %genspc = ();
foreach my $fn (@args) {
	my $sl = rindex($fn, '/');
	my $gs = substr($fn, $sl + 1, $ngsid);
	$genspc{$gs} .= "$fn ";
}
foreach my $gs (sort keys(%genspc)) {
	my @a = split(' ', $genspc{$gs});
	&samegr(\@a, $gs);
}

my @data;
my @loci;

sub byCpos {
	return ($loci[$a] cmp $loci[$b]) if ($loci[$a] ne $loci[$b]);
	my ($al, $ar) = split(' ', $data[$a], 2);
	my ($bl, $br) = split(' ', $data[$b], 2);
	return ($ar <=> $br) if ($ar != $br);
	return ($al <=> $bl);
}

sub injection {
	my $tru = shift;
	foreach my $rec (sort keys(%$tru)) {
	    my @a = split(':', $$tru{$rec});
	    next unless(@a);
	    my ($cor, $ovr) = ("", 0);
	    foreach my $snd (@a) {
		my ($cor2, $ovr2) = split(' ', $snd);
		if ($ovr2 > $ovr) {
		    $cor = $cor2;
		    $ovr = $ovr2;
		}
	    }
	    if ($rev) {
		print $cor,"\t",$rec,"\n";
	    } else {
		print $rec,"\t",$cor,"\n";
	    }
	}
}

sub samegr {
	my ($sqs, $gs) = @_;
	@data = ();
	@loci = ();
	my $npath = 0;
	my %paths = ();
	my %trunit_homl = ();
	my %trunit_homr = ();
	my %trunit_hetl = ();
	my %trunit_hetr = ();

	foreach my $fn (@$sqs) {
	    $fn = $pas . $fn;
	    my ($gname, $rngl, $rngr, $sns, $dmy, @b) = split(' ', `ExonLen -b $fn`);
	    $gname = substr($gname, 1) if (substr($gname, 0, 1) eq '$');
	    my $l = $b[0];			# left boundary
	    my $r = $b[$#b];			# right boundary
	    my $n = rindex($fn, '/');		# path
	    push(@loci, "$gname $sns");
	    my $dr = substr($fn, 0, $n);
	    my $fn = substr($fn, $n + 1);
	    push(@data, "$l $r $dr $fn");
	    my $np = $paths{$dr};
	    $paths{$dr} = $np = ++$npath unless ($np);
	    if ($np == 1)	{$trunit_hetl{$fn} = "";}
	    else		{$trunit_hetr{$fn} = "";}
	}

	my @odr = sort byCpos (0..$#loci);

	for (my $i = 0; $i < $#loci; ++$i) {
	    my $oi = $odr[$i];
	    my ($li, $ri, $dri, $fni) = split(' ', $data[$oi]);
	    for (my $j = $i + 1; $j <= $#loci; ++$j) {
		my $oj = $odr[$j];
		last if ($loci[$oi] ne $loci[$oj]);	# different chromosome
		my ($lj, $rj, $drj, $fnj) = split(' ', $data[$oj]);
		last if ($ri < $lj || $rj < $li);	# no overlap
		my $rr = ($ri < $rj)? $ri: $rj;
		my $ll = ($li > $lj)? $li: $lj;
		my $ovr = $rr - $ll;
		next if ($ovr < $ovrthr);		# under threshold
		my ($drl, $drr, $fnl, $fnr) = ($dri, $drj, $fni, $fnj);
		if ($paths{$dri} == 2) {
		    $drl = $drj; $fnl = $fnj;
		    $drr = $dri; $fnr = $fni;
		}
		my $hom = $drl eq $drr;
		if (!$uniq && !$summary) {
		    if (!$omode || ($hom && $omode & 1) || (!$hom && $omode & 2)) {
			if ($outform == 1) {
			    printf("%d %s %s\n", $ovr,  $fnl,  $fnr);
			} elsif ($outform == 2) {
			    printf("%-15s %-15s %s %d\n", $fnl, $fnr, $loci[$oi], $ovr);
			} else {
			    printf("%-15s %-15s %-15s %-15s %s %d\n", 
				$drl, $fnl, $drr, $fnr, $loci[$oi], $ovr);
			}
		    }
		} elsif ($summary) {
		    if ($hom) {
			if ($paths{$drl} == 1)	{$trunit_homl{$fnr} = $fnl;}
			else			{$trunit_homr{$fnr} = $fnl;}
		    } else {
			$trunit_hetl{$fnl} .= "$fnr ";
			$trunit_hetr{$fnr} .= "$fnl ";
		    }
		}  else {
		    if ($hom) {
			if ($paths{$drl} == 1)	{$trunit_homl{$fnr} = "$fnl $ovl:";}
			else			{$trunit_homr{$fnr} = "$fnl $ovl:";}
		    } else {
			$trunit_hetl{$fnl} .= "$fnr $ovr:";
			$trunit_hetr{$fnr} .= "$fnl $ovr:";
		    }
		}
	    }
	}

	if ($uniq) {
	    if ($outform == 0) {
		if ($summary) {
		    my $count = @$sqs - %trunit_homl;
		    print "$gs\t$count\n";
		} elsif ($omode == 1) {
		    &injection(\%trunit_hetl);
		} else {
		    &injection(\%trunit_hetr);
		}
	    } else {
		my ($total, $good) = (0, 0);
		foreach my $fn (@$sqs) {
		    my $glen = 0;
		    if ($trunit_homl{$fn}) {
			my @a = split(' ', `utp -n $trunit_homl{$fn}`);
			$glen = $a[5] if ($lwrlmt <= $a[5] && $a[5] <= $uprlmt);
		    }
		    unless ($glen) {
			my @a = split(' ', `utp -n $fn`);
			$glen = $a[5] if ($lwrlmt <= $a[5] && $a[5] <= $uprlmt);
		    }
		    ++$total;
		    ++$good if ($glen);
		    print "$fn\t$glen\n" if ($omode == 2);
		}
		print "$gs\t$total\t$good\n" if ($omode == 1);
	    }
	} elsif ($summary) {
	    if ($omode == 1) {
		my ($tl, $ul, $cl, $sl) = (0, 0, 0, 0);
		foreach my $bnl (keys(%trunit_hetl)) {
		    unless ($trunit_hetl{$bnl}) {
			++$sl;
		    } elsif (!$trunit_homl{$bnl}) {
			++$cl;
		    }
		    ++$ul unless ($trunit_homl{$bnl});
		    ++$tl;
		}
		my ($tr, $ur, $cr, $sr) = (0, 0, 0, 0);
		foreach my $bnr (keys(%trunit_hetr)) {
		    unless ($trunit_hetr{$bnr}) {
			++$sr;
		    } elsif (!$trunit_homr{$bnr}) {
			++$cr;
		    }
		    ++$ur unless ($trunit_homr{$bnr});
		    ++$tr;
		}
		printf("%s\t%7d\t%7d\t%7d\t%7d\t%7d\t%7d\t%7d\t%7d\n", $gs, 
		    $tl, $tr, $ul, $ur, $cl, $cr, $sl, $sr);
	    }
	    if ($omode & 2) {
		foreach my $bnr (sort keys(%trunit_homl)) {
		    print "$bnr\t$trunit_homl{$bnr}\n";
		}
		foreach my $bnr (sort keys(%trunit_homr)) {
		    print "$bnr\t$trunit_homr{$bnr}\n";
		}
		print "\n" if ($omode & 12);
	    }
	    if ($omode & 4) {
		foreach my $bnl (sort keys(%trunit_hetl)) {
		    print "$bnl\t$trunit_hetl{$bnl}\n";
		}
		print "\n" if ($omode & 8);
	    }
	    if ($omode & 8) {
		foreach my $bnr (sort keys(%trunit_hetr)) {
		    print "$bnr\t$trunit_hetr{$bnr}\n";
		}
	    }
	}
}
