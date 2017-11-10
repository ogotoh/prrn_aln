#!/usr/bin/perl

##############################################################################
#
#	Divide the given UPGMA tree into subtrees with a cut-off height
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
#___________________________________________________________________________
#
#	Usage: divtree.pl [-HMaxHeight] [-XMaxCluster] [-NMinCluster] upgma_tree
#	Options:
#	  MaxHeight: cut-off height
#	  MaxCluster: maximum number of OTUs in a subtree
#	  MinCluster: minimum number of OTUs in a subtree
#
#############################################################################

my $Dist = 110;		# Threshold between subtrees
my $MinMem = 3;		# Min members
my $MaxMem = 0x7fffffff;	# Max members
my $database;		# Blast DB
my $dir = "./";		# Directory of query
my ($No, $Rep, $Sub, $Akin, $Uniq) = (0, 1, 2, 3, 4);
my $OutNo;		# Output members in no rather than name
my $repre = $No;	# Type of output
my $eof = 0;
my $inc_outlier;	# include outliers to subtree

my %blasuite = (
	'p', 'blp',
	'n', 'bln',
	't', 'tbln',
	'x', 'blx'
);

sub usage {
    print STDERR "Usage: divtree.pl [-Hthr -Nmin_member -Xmax_number [-A|-R|-U|-B[pntx]]] tree\n";
    print STDERR "Output\t\ttree.N\t(N= 1..)\n";
    print STDERR "\t-R\tRepresentatives\n";
    print STDERR "\t-B\tblast for repls\n";
    exit;
}

while ($_ = $ARGV[0], /^-/) {
	shift;
	last if /^-$/;
	/^-D/ && ($Debug = 1);
	/^-d(\S+)/ && ($database = $1);
	/^-s(\S+)/ && ($dir = $1);
	/^-H(\S+)/ && ($Dist = $1);
	/^-I/	   && ($inc_outlier = 1);
	/^-N(\d+)/ && ($MinMem = $1);
	/^-X(\d+)/ && ($MaxMem = $1);
	/^-A(\S*)/ && ($repre = $Akin, $OutNo = $1);
	/^-R(\S*)/ && ($repre = $Rep, $OutNo = $1);
	/^-S(\S*)/ && ($repre = $Sub, $OutNo = $1);
	/^-U/ && ($repre = $Uniq);
	/^-B(\S*)/ && ($blast = $1? $blasuite{$1}: 'tbln');
	/^-\?/ && &usage;
}

my $sf = 0;	# subtree no.
my $nr = 0;	# non-blank line no.
my $delim = 1;
my @buf;
my $pline;
$MinMem *= 2;
$MaxMem = $MinMem + 1 if ($MaxMem <= $MinMem);
$MaxMem *= 2;

if ($dir && substr($dir, -1, 1) ne "/") {$dir .= "/";}

sub writetree {
	my ($n, $r) = @_;
	++$sf;
	my $name = ($#ARGV? $ARGV: "tree");
	$name = '_' if ($name eq '-');
	$name .= ".$sf";
	open(TREE, ">$name") || return;
	for ( ; $n < $r; ++$n) {
	   print TREE $buf[$n];
	}
	close(TREE);
}

sub wtree {
	my ($n, $r) = @_;
        for ( ; $n < $r; ++$n) {
           print STDERR $buf[$n];
        }
}

sub divtree {
    my ($lnode, $rnode) = @_;
#    printf STDERR ("%d %d %d\n", $lnode, $rnode, $rnode - $lnode);
    return writetree(@_) if ($rnode - $lnode < $MaxMem);
    my $x = my $mx = $lnode + $MinMem - 1;
    while (($x += 2) <= $rnode - $MinMem) {
	if ($buf[$x] > $buf[$mx]) {
	    $mx = $x;
	}
    }
    &divtree($lnode, $mx);
    &divtree($mx + 1, $rnode);
}

my $otu;
my $MemNo = 0;

while (<>) {
again:
    my ($clmn1) = split;	# 1st column
    next unless ($clmn1);	# Skip blank line
    $odd = $nr++ % 2;		# Branch height
    $otu = $clmn1 unless($odd);	# OTU
    $MemNo++ unless ($odd);
    print $OutNo? $MemNo: $otu, "\n" if ($repre == $Rep && $delim);
    if ($delim && $blast) {
	$query = $dir . $otu;
	$cbst = $otu . ".cbst";
	my $cmd =  "$blast -H0 $database $query | blasort.pl > $cbst";
	if ($Debug) {
	    print "$cmd\n";
	} else {
	    system("$cmd");
	    exit if ($? >> 8);
	}
    }
    $delim = $odd && $clmn1 > $Dist;
    if ($repre == $Akin || $repre == $Sub) {
	if ($delim) {
	    shift(@buf) if ($repre == $Sub);
	    if (@buf) {
		if ($dir) {
		    my @slen = ();
		    my @sbuf = ();
		    foreach my $mem (@buf) {
			my $fn = $dir . $mem;
			if (-s $fn && index($fn, ':') < 0) {
			    my @a = split(' ', `utp -n $fn`);
			    push(@slen, pop(@a));
			} else {
			    push(@slen, 0);
			}
		    }
		    my @odr = sort {$slen[$b] <=> $slen[$a];} (0 .. $#buf);
		    foreach my $i (0 .. $#buf) {
			push(@sbuf, $buf[$odr[$i]]);
		    }
		    @buf = @sbuf;
		}
		foreach my $mem (sort @buf) {
		    print $mem;
		    print (($repre == $Sub)? "\n": " ");
		}
		print "\n" if ($repre == $Akin);
		@buf = ();
	    }
	} elsif (!$odd) {
	    push(@buf, $OutNo? $MemNo: $otu);
	}
    } elsif ($repre == $Uniq) {
	if ($delim) {
	    $pline = $_;
	} elsif ($odd) {
	    $pline = "";
	} elsif ($pline) {
	    print $pline, $_;
	} elsif (!defined($pline)) {
	    print $_;
	}
    }
    next if ($repre || $blast);
    push(@buf, $_);
    if ($delim) {
	if ($#buf > $MinMem) {	# write to files
	    &divtree(0, $#buf);
	} else {
	    while ($_ = shift(@buf)) {
		if (/^\s*(\S+)/ && $1 !~ /^\d/) {
		    print STDERR $1, "\n";
		}
	    }
	    print STDERR "\n";
	}
	@buf = ();
    }
    last if ($eof);
}
continue {
    if (eof() && @buf > 0) {
	$eof++;
	$_ = 999999999;
	goto again;
    }
}
