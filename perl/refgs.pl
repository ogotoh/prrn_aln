#!/usr/bin/perl

##############################################################################
#
#	refgs.pl: refine gene structure prediction by iteration
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
#____________________________________________________________________________
#
#	Usage:	refgs.pl -sSrc -tTgt -I# [-R|-R1] MSA
#	 or	refgs.pl -sSrc -tTgt template seq1 seq2 ...
#	Options:
#	-sSrc:	directory containing extended fasta files before refinement
#	-tTgt:	directory containing extended fasta files after refinement
#	-I#:	maximum number of iterations
#	-R:	refinement with CR mode (closest reliable)
#	-R1:	refinement with PR mode (profile reliable)
#	none:	refinement with M1 mode (minus one)
#	MSA	GSA-MPSA constructed from aa sequences in Src directory
#
#############################################################################

use lib "$ENV{HOME}/perl";	# perl module
use Util;			# utilities
use SspTab;			# species specific table
use Dixon;			# Dixon's test of outliers
use Exin;			# Exon-intron boundaries
use File::Copy;			# equivalent to system cp
use strict;

my $INT_MAX = 2147483647;
my $tmpdir = "tmp/";
mkdir $tmpdir unless (-d $tmpdir);
my $TmpRefGs = $tmpdir . "RefGs";
my $cmngnmdb;
#my $yopt = "-yh1 -ys1";
my $yopt;
my $aln = "aln -yl2 -L -pq -H35 -yh1 ";
my $spaln = "spaln -pq -H35 ";
my $opt;
my $prrn = "prrn5 -KP -U -ye0.5 -H25 -I1 -Q1 -pq ";
my $reference;	# reference sequence or MSA
my $nas = "";	# directory wherein genomic DNA sequences (optional)
my $org = "";	# directory wherein initial prediced sequences
my $dst = "";	# directory to which the refined aa sequences are written
my $list;	# list of files to be tested
my $lmrgn;	# left margin
my $rmrgn;	# right margin
my $ave = 100;	# average exon length;
my $maxN = 10;	# max allowed N's
my $alpha = 0.1;	# error limit for outlier
my $allowedN = 50;	# allowable N's
my $allowedX = 1;	# allowable X's
my $allowedQ = 2;	# allowable defects
my $LongGap;		# Remove sequence with long gap
my $reCalcWt = 2;# 0: even weight, 1: seq weight, 2: pair weight 3: recalc
my $mem;	# number of members in current MSA
my $sitr = 0;	# status of iteration
my $vbs = 0;	# Verbose output
my $repeat = 0;	# Max number of iterations
my $out = "";	#
my $log;	# log file
my @fn;		# seq. name
my @pw;		# pair weight
my $GrpDist;	# Threshold between groups
my $reportrel;	# Retain references if ON
my $quiet;
my $relyonpp = 1;	# rely on the present prediction
my %tmpfiles;
my $cogonly;
my $ncommon = 2;
my $tlmt = 15;
my $ilmt = 20;
my $AmbN_thr = 10;
my %locrng;
my %leftbnd;
my %lb;
my %rb;
my $keep_unrel = 1;
my $max_plen = 1000;
my $min_plen = 0;
my $in_genetic = 0;
my $mklink;
my $defectP = 2;
my $omode;
my $TempPrf;
my $TempPrf2;
my $remove_processed_pg = 1;
my $ref_bka;
my $lendev = 2.0;
my $min_len_var = 2.0;

sub usage {
	print "Usage:\n";
	print "\trefgs [options] -In MSA\n";
	print " or\trefgs [options] MSA [n1 n2 ..| n-m]\n";
	print " or\trefgs [options] template seq1 seq2 ...\n";
	print "\nOptions:\n";
	print "\t-ln:	left mergin\n";
	print "\t-rn:	right mergin\n";
	print "\t-In:	Maximum number of iterations\n";
	print "\t-RS:	Reliable reference set\n";
	print "\t-R1:	refinement with CR mode (closest reliable)\n";
	print "\t-R:	refinement with PR mode (profile reliable)\n";
	print "\tnoR:	refinement with M1 mode (minus one)\n";
	print "\t-sS:	directory of present protein sequences\n";
	print "\t-tS:	directory of refined protein sequences\n";
	print "\t-SS:	directory of genomic DNA sequence\n";
	print "\t-On:	Output mode of aln\n";
	print "\t-V:	Verbose output\n";
	exit(1);
}

sub System {
	my $cmd = shift;
	print STDERR "$cmd\n" if ($vbs & 1);
	system($cmd) == 0 or die "$cmd failed !\n";
}

sub lcscmp {
	return ($$a[1] <=> $$b[1]) if ($$a[0] == $$b[0]);
	return ($$a[0] <=> $$b[0]);
}

sub setboundary {
    foreach my $sq (@_) {
	my $fnm = $org . $sq;
	next unless (-s $fnm);
	my ($seq, $rngl, $rngr, $sns, $o15, @bnd) = &exin($fnm, 'b');
	next unless (@bnd && ($o15 eq 'N' || $o15 eq 'Q'));
	my $ll = shift(@bnd);
	my $rr = pop(@bnd);
	my $sqs = $seq . $sns;
        $locrng{$sq} = [ $sqs, $ll, $rr ];
	push @{ $leftbnd{$sqs} }, [$ll, $rr];
    }
    foreach my $sqs (keys(%leftbnd)) {
	my @ra = sort lcscmp @{ $leftbnd{$sqs} };
	$leftbnd{$sqs} = [ @ra ];
    }
    foreach my $sq (@_) {
	next unless ($locrng{$sq});
	my $sqs = $locrng{$sq}[0];
	my $ll = $locrng{$sq}[1];
	my $rr = $locrng{$sq}[2];
	$lb{$sq} = 0;
	$rb{$sq} = $INT_MAX;
	my $na = @{ $leftbnd{$sqs} };
	for (my $i = 0; $i < $na; ++$i) {
	    my $l = $leftbnd{$sqs}[$i]->[0];
	    my $r = $leftbnd{$sqs}[$i]->[1];
	    if ($l < $ll) {$lb{$sq} = ($r < $ll)? $r: $ll;}
	    if ($r > $rr) {
		$rb{$sq} = ($l > $rr)? $l: $rr;
		last;
	    }
	}
    }
}

while ($_ = $ARGV[0], /^-/) {
    shift;
    last if /^--$/;
    if (/^-\?/ || /^-h/)	{&usage;} 
    elsif (/^-a(\S*)/)	{&Util::getoptarg(\$ref_bka, $1);}
    elsif (/^-c(\S*)/)	{&Util::getoptarg(\$ncommon, $1);}
    elsif (/^-f(\S*)/)	{&Util::getoptarg(\$relyonpp, $1);}
    elsif (/^-g(\S*)/)	{&Util::getoptarg(\$cmngnmdb, $1);}
    elsif (/^-i(\S*)/)	{&Util::getoptarg(\$ilmt, $1);}
    elsif (/^-k(\S*)/)	{&Util::getoptarg(\$keep_unrel, $1, 9999);}
    elsif (/^-j(\S*)/)	{&Util::getoptarg(\$tlmt, $1);}
    elsif (/^-l(\S*)/)	{&Util::getoptarg(\$lmrgn, $1);}
    elsif (/^-m/)	{$mklink = 1;}
    elsif (/^-n(\S*)/)	{&Util::getoptarg(\$in_genetic, $1);}
    elsif (/^-o(\S*)/)	{&Util::getoptarg(\$out, $1);}
    elsif (/^-p(\S*)/)	{&Util::getoptarg(\$allowedQ, $1);}
    elsif (/^-q(\S*)/)	{&Util::getoptarg(\$quiet, $1, 1);}
    elsif (/^-r(\S*)/)	{&Util::getoptarg(\$rmrgn, $1);}
    elsif (/^-s(\S*)/)	{&Util::getoptarg(\$org, $1);}
    elsif (/^-t(\S*)/)	{&Util::getoptarg(\$dst, $1);}
    elsif (/^-u(\S*)/)	{&Util::getoptarg(\$reportrel, $1);
	 $reference = 0; $keep_unrel = -1; $LongGap = 1;}	 # diagonize
    elsif (/^-w(\S*)/)	{&Util::getoptarg(\$reCalcWt, $1, 3);}
    elsif (/^-A(\S+)/)	{&Util::getoptarg(\$alpha, $1);}
    elsif (/^-C(\S*)/)	{$cogonly = 1;}
    elsif (/^-G(\S*)/)	{&Util::getoptarg(\$GrpDist, $1);}
    elsif (/^-I(\S*)/)	{&Util::getoptarg(\$repeat, $1);}
    elsif (/^-L(\S+)/)	{if ($1 eq 'S') {$yopt .= " -LS";} else {$list = $1;}}
    elsif (/^-N(\S+)/)	{&Util::getoptarg(\$maxN, $1);}
    elsif (/^-O(\S*)/)	{&Util::getoptarg(\$omode, $1, 1);}
    elsif (/^-P(\S*)/)	{$keep_unrel = 0 if (!$1 || $1 & 1); $LongGap = 1 if ($1 & 2);}
    elsif (/^-R(\S*)/)	{&Util::getoptarg(\$reference, $1);
	$reference = 0 unless(defined($reference));}
    elsif (/^-S(\S*)/)	{&Util::getoptarg(\$nas, $1);}
    elsif (/^-V(\S*)/)	{&Util::getoptarg(\$vbs, $1, 1);}
    elsif (/^(-y\S*)/)	{$yopt .= " $1";}
    else {$opt .= " " . $_;}
}

$nas .= '/' if ($nas && substr($nas, -1, 1) ne '/');
$org .= '/' if ($org && substr($org, -1, 1) ne '/');
$dst .= '/' if ($dst && substr($dst, -1, 1) ne '/');
my $pas = $org;
my $wrk = ($dst && $dst ne $pas)? $dst: $tmpdir; 
my $prof = $reference == 1;	# use profile reference
my $glmt = $tlmt + $ilmt if (defined($tlmt) && defined($ilmt));
$repeat = 1 if (($reference || $ref_bka) && !$prof);

if ($yopt) {$prrn .= $yopt; $opt .= $yopt;}

########## setup ##########

&opengnmdb($cmngnmdb);	# initialize gnm2tab

if ($list) {
	open(LST, $list) or die "Can't open $list !\n";
	while (<LST>) {
	    my ($fn, $rst) = split(' ', $_, 2);
	    push(@ARGV, $fn) if ($fn =~ /^[A-Z]/);
	}
	close(LST);
}

# Use the most similar reference in the list
if ($reference =~ /^[a-zA-Z]/) {
	my %defs = ();
	&mostsim(\%defs, \@ARGV);
	&rmtemp;
}
	
my ($msa, $attr) = split(' ', shift(@ARGV), 2);
$reportrel = "PQRS" if ($vbs & 4 || $reportrel =~ /A/);

my @seqlist;

if (@ARGV) {
	@seqlist = @ARGV;
} else {
	@seqlist = split(' ', `utp -l0 $msa`);
}

&setboundary(@seqlist);

########## iterative refinement ##########

my $msaorg = $msa;
my $iter = 0;

while (!$repeat || $iter++ < $repeat) {
	my $upd = $out? $out: "$msaorg.$iter";
	my ($rev, $prv) = (undef, undef);
	if (defined($reference)) {
	    ($rev, $prv) = &msa2ref($iter, $upd);
	} elsif ($iter == 1 && $ref_bka) {
	    ($rev, $prv) = &external_reference();
	    $rev = "" unless (defined($rev));
	} else {
	    my @seqs = ();
	    my @ndefs = ();
	    my @cdefs = ();
	    &toolongorshort($msa, \@seqs, \ @ndefs, \@cdefs);
	}
	my $now = defined($prv)? $prv: $msa;
	$rev = &onecycle($now) unless (defined($rev));
	last unless ($repeat);
	if ($rev eq "" && !$prv) {	# no more refinement
	    unlink($now) if (defined($prv));
	    last;
	}
	$now = "./$now" if (index($now, "/") < 0);
	&System("$prrn -s$wrk -o$upd $now $rev");
	unlink($now) if (defined($prv));
	my $temp = $wrk;
	$wrk = ($wrk eq $tmpdir)? $pas: $tmpdir;
	$pas = $temp;
	$msa = $upd;
	$sitr = 0;
}
&rmtemp;

########## end of main ##########

# find reliable references from MSA

sub msa2ref()
{
	my ($iter, $upd) = @_ ;
	my $rmvp = $iter >= $keep_unrel;			# remove P type
	my @samepat = ();
	my @lone = ();
	my @seqs = ();
	my @ndefs = ();
	my @cdefs = ();
	&toolongorshort($msa, \@seqs, \@ndefs, \@cdefs);	# total length
	return ("") unless (@seqs);
	&outliers($msa, \@ndefs, \@cdefs);
	my ($avrneij, $maxns) = &eijpattern($msa, \@seqs, \@samepat, \@lone);
	&pseudogene(\@seqs, \@ndefs, \@cdefs) if ($rmvp);
	&processed_pg($msa, \@ndefs, \@cdefs) if ($rmvp && $remove_processed_pg);
	my @rel = ();
	my @top = ();
	my @seed = ();
	my @unc = ();
	my @rmv = ();
	my $prv = undef;
	my $conspat = 2 * $maxns > @seqs;
	my $max_seq_name;
	for (my $ns = 0; $ns < @seqs; ++$ns) {
	    my $sq = $seqs[$ns];
	    my $ln = length($sq);
	    $max_seq_name = $ln if ($ln > $max_seq_name);
	    my $fn = $org . $sq;
	    if (!(-w $fn)) {
		push(@seed, $ns);
		next;
	    }
	    $fn = $pas . $sq;
	    if (!(-s $fn) && !$quiet && $mklink) {
		print STDERR "$fn not found or empty!\n";
	    }
	    my $ndef = $ndefs[$ns];
	    --$ndef if ($samepat[$ns] >= $ncommon && !$reportrel);
	    if ($ndef <= 0) {
		push(@rel, $ns);
		push(@top, $ns) if ($conspat && $samepat[$ns] == $maxns);
	    } elsif ($ndef < $allowedQ || !$rmvp) {
		push(@unc, $ns);
	    } else {
		push(@rmv, $ns);
	    }
	}
	if ($reportrel =~ /[A-Z]/) {
	    my $fmt = '%-' . sprintf("%d", $max_seq_name);
	    $fmt .= "s\t%s\t%d\t%s %s\n";
	    if ($reportrel =~ /T/) {
		foreach my $ns (@top) {
		    printf($fmt, $seqs[$ns], $msa, $ns + 1, 'T', $cdefs[$ns]);
		}
	    }
	    if ($reportrel =~ /S/) {
		foreach my $ns (@seed) {
		    printf($fmt, $seqs[$ns], $msa, $ns + 1, 'S', $cdefs[$ns]);
		}
	    }
	    if ($reportrel =~ /R/) {
		foreach my $ns (@rel) {
		    printf($fmt, $seqs[$ns], $msa, $ns + 1, 'R', $cdefs[$ns]);
		}
	    }
	    if ($reportrel =~ /Q/) {
		foreach my $ns (@unc) {
		    printf($fmt, $seqs[$ns], $msa, $ns + 1, 'Q', $cdefs[$ns]);
		}
	    }
	    if ($reportrel =~ /P/) {
		foreach my $ns (@rmv) {
		    my $c = ($cdefs[$ns] =~ /N/)? 'N': 'P';
		    printf($fmt, $seqs[$ns], $msa, $ns + 1, $c, $cdefs[$ns]);
		}
	    }
	    return ("") unless ($vbs & 4);
	} elsif ($vbs & 2) {
	    print "S:", join(' ', @seed), "\n";
	    print "R:", join(' ', @rel), "\n";
	    print "Q:", join(' ', @unc), "\n";
	    print "P:", join(' ', @rmv), "\n";
	}
	if ($reference || $relyonpp) {
		@rel = (@seed, @rel);
	} else {		# first time
	    if (@seed) {	# exam non-seed
		@unc = (@unc, @rel);
		@rel = @seed;
	    }
	}
	return (undef) unless (@rel);	# no reliable reference

	my $rel_n = "";
	if ($prof) {
	    foreach my $r (@rel) {
		my $r1 = $r + 1;
		$rel_n .= "$r1 ";
	    }
	}
	$reference = "Rel$$";
	if ($prof && $rel_n) {
	    my $cmd = "rdn -csd $msa $rel_n > $reference";
	    &System($cmd);
	} else {
	    open(CAT, "> $reference") or die "Can't write to $reference !\n";
	    foreach my $ns (@rel) {
		print CAT $seqs[$ns], "\n";
	    }
	    close(CAT);
	}

########## convert seq_no to seq_name ########## 

	for (my $i = 0; $i < @rel; ++$i) {$rel[$i] = $seqs[$rel[$i]];}
	for (my $i = 0; $i < @unc; ++$i) {$unc[$i] = $seqs[$unc[$i]];}
	my %defs = ();
	for (my $i = 0; $i < @seqs; ++$i) {$defs{$seqs[$i]} = $cdefs[$i];}

########## update questionable predictions  ########## 

	my $revised = &mostsim(\%defs, \@rel, \@unc, \@ARGV);
	unlink($reference);

########## remove defective genes ########## 

	if (@ARGV) {
	    @unc = sort(split(' ', $revised), @ARGV);
	    for (my $i = 1; $i <= $#unc; ) {	# unique
		if ($unc[$i] eq $unc[$i-1]) {
		    splice(@unc, $i, 1);
		} else {
		    $i++;
		}
	    }
	    @ARGV = ();
	    $revised = join(' ', @unc);
	}
	if (@rmv) {
	    my %rmv;
	    my $rnm;
	    foreach my $ns (@rmv) {
	        $rmv{$seqs[$ns]}++;		# to unique
		my $r = $ns + 1;
		$rnm .= "$r ";
		print STDERR "Remove $seqs[$ns]\n" unless($quiet);
	    }
	    for (my $i = 0; $i < @unc; ) {
		if ($rmv{$unc[$i]}) {
		    splice(@unc, $i, 1);
		} else {
		    ++$i;
		}
	    }
	    $revised = join(' ', @unc);
	    $prv = $upd.$$;
	    &System("rdn -ced $msa $rnm > $prv");
	    unless (-s $prv) {
		unlink($prv);
		$prv = undef;
	    }
	}
	$revised = "" unless (defined($revised));
	return ($revised, $prv);
}

sub onecycle {
    my $msa = shift;
    my @oklist = ("", "");	# list of remain unchanged
    unless ($ARGV[0]) {
	&initilize unless $sitr;
	if ($GrpDist) {
	    my @groups;
	    my $gn = 0;
	    open(SS, "phyln $msa | divtree.pl -AN -H$GrpDist |")
		or die "Can't run phyln!";
	    while (<SS>) {
		my @subfmem = split;
		foreach my $m (@subfmem) {
		    $oklist[&conf($fn[$m-1], &reref(@subfmem), $m)] .= "$fn[$m-1] ";
		}
	    }
	    close(SS);
	} else {
	    for (my $n = 0; $n < $mem; ++$n) {
		$oklist[&conf($fn[$n], &reref($n+1), $n+1)] .= "$fn[$n] ";
	    }
	}
    }

    my @retry;
    my $refgp;
    foreach my $arg (@ARGV) {
	if ($arg =~ /^(\d+)/) {
	    &initilize unless $sitr;
	    if ($arg =~ /-/) {
		my ($i, $j) = split(/-/, $arg);
		for (my $k = $i; $k <= $j; ++$k) {
		    push(@retry, $k);
		}
	    } else {
		push(@retry, $1);
	    }
	    $refgp = &reref(@retry);
	} else {
	    unless (-e $pas.$arg) {
		print STDERR "Can't found $arg\n";
		next;
	    }
	    $oklist[&conf($arg, $msa, 0)] .= "$arg ";
	}
    }
    if ($refgp) {
	for (my $n = 0; $n <= $#retry; ++$n) {
	    my $m = $retry[$n];
	    $oklist[&conf($fn[$m-1], $refgp, $m)] .= "$fn[$m-1] " if ($m);
	}
    }
    unlink($TempPrf)  if ($TempPrf);
    unlink($TempPrf2) if ($TempPrf2);
    $sitr = 0 if ($reCalcWt > 2); # recalculate pair weights
    return ($oklist[1]);
}

sub elem {
    my ($i, $j) = @_;
    ($i > $j)? $i * ($i-1)/2 + $j: $j * ($j-1)/2 + $i;
}

sub initilize {
    $sitr = 1;
    open(PHYL, "phyln -w -W \'$msa\' |") or
	 die "Can't exec phyln!\n";
    $mem = 0;
    while (<PHYL>) {
	my @a = split;
	$fn[$mem] = shift(@a);
	$pw[$mem++] = shift(@a);
    }
    close(PHYL);
    $TempPrf = $TmpRefGs . $$;
    $TempPrf2 = $TempPrf . "_";
    open(BODY, $msa) or die "Can't open $msa!\n";
}

sub reref {
    my @ref;
    my $delseq;

    foreach my $i (0 .. $#fn) {
	$ref[$i] = 0;
    }
    foreach my $i (@_) {
	$ref[$i - 1] = 1;
    }
    open(REF, ">$TempPrf") or die "Can't open $TempPrf!\n";
    seek(BODY, 0, 0);
    while (<BODY>) {
	print REF unless /^%/;
	if ($reCalcWt && /^>/) {
	    print REF "\n";
	    for (my $i = 0; $i <= $#fn; ) {
		my $wt;
		print REF "%" unless $i % 5;
		$wt = $ref[$i]? 0: $pw[$i];
		printf(REF " %14.7e", $wt);
		print REF "\n" unless (++$i % 5);
	    }
	    print REF "\n";
	}
    }
    close(REF);
    foreach my $i (@_) {
	$delseq .= "$i ";
    }
    return (undef) unless (-s $TempPrf);
    &System("rdn -ced $TempPrf $delseq > $TempPrf2");
    unless (-s $TempPrf2) {
	unlink($TempPrf2);
	return (undef);
    } else {
	return ($TempPrf2);
    }
}

sub skip_aln {
    my ($src, $dst) = @_;
    unless (-s $src) {
	if ($mklink) {
	    print STDERR "$src can't open!\n" unless ($quiet);
	    return (3);
	}
	return (1);
    }
    if (-e $dst && !(-w $dst)) {
	return (1);	# exist and protected
    }
    if (!(-w $src)) {	# protected
	if ($mklink && (!-e $dst || -w $dst)) {
	    print STDERR "$src is protected!\n" unless ($quiet);
	    unlink($dst) if (-e $dst);
	    link($src, $dst);
	}
	return (1);
    }
    return (0);
}

sub copy_or_link {
    my ($src, $dst) = @_;
    return (0) if (-e $dst && !-w $dst);
    if (-w $src) {copy($src, $dst);}
    elsif ($mklink) {
	unlink($dst) if (-e $dst);
	link($src, $dst);
    } else {return (0);}
    return (1);
}

sub conf {	# refine each gene structure
    my ($name, $refp, $nn, $qlm, $qrm) = @_;
    my $src = $pas . $name;
    my $temp = $wrk . $name;
    my $lm = $lmrgn;
    my $rm = $rmrgn;
    my $avi = &avrintlen($name) unless ($lm && $rm);
    unless ($lm) {
	$lm = $avi;
	$lm *= sqrt($lm) if ($qlm);
	$lm += $ave;
    }
    unless ($rm) {
	$rm = $avi;
	$rm *= sqrt($rm) if ($qrm);
	$rm += $ave;
    }
    my $skip = &skip_aln($src, $temp);
    return ($skip - 1) if ($skip);
    my ($seq, $rngl, $rngr, $sns, $o15, @bnd) = &exin($src, 'b');
    unless (($seq =~ /^\$/ && ($sns eq '+' || $sns eq '-')) 
	&& (@bnd && ($o15 eq 'N' || $o15 eq 'Q'))) {
	&copy_or_link($src, $temp) unless (-s $temp);	# not for subject of refinement
	return (1);
    }
    my $alnargs;
    unless (-s $refp) {
	print STDERR "$refp can't open!\n";
	return (2);
    }
    if (substr($seq, 0, 1) eq "\$") {
	my $ssp = &sspopt($seq);
	$alnargs = $ssp . &gnmdb($seq);
    } else {
	$alnargs = &sspopt($seq);
	$seq = $nas . $seq;
    }
    my $cmp = $sns eq "-";
    $sns = $cmp? " <": "";
    my $left = shift(@bnd) - ($cmp? $rm: $lm);
    if ($left < $lb{$name}) {$left = $lb{$name};}
    elsif ($left < 0) {$left = 0;}
    my $right = pop(@bnd) + ($cmp? $lm: $rm);
    $right = $rb{$name} if ($rb{$name} && $right > $rb{$name});
    my $gene = sprintf("\'%s %d %d %s G\'", $seq, $left, $right, $sns);
    if ($attr) {
	$refp = sprintf("\'%s %s\'", $refp, $attr);
    }
    my ($id, $lb, $qmany) = split(' ', `utp -n $refp`);
    my $alncmd = $qmany == 1? $spaln: $aln;
    $alnargs .= " $opt $gene '$refp P'";
    if (defined($omode)) {
	&System("$alncmd -O$omode $alnargs");
	return (2);
    }
    $tmpfiles{$temp}++ if ($wrk eq $tmpdir);
    &System("$alncmd -O7 -o $temp $alnargs");
    chmod(0666, $temp);
    my @u = split(' ', `utp -n $temp`);
    my $ulen = pop(@u);
    if ($ulen < $min_plen || $max_plen < $ulen) {
	printf STDERR "$temp has abnormal length: %7.1f < %d < %7.1f\n", 
		$min_plen, $ulen, $max_plen;
	return (2);
    }
    my @msg = split(' ', `iden -n -KP $temp $src`);
    if ((!@msg || $msg[2]) && $wrk eq $tmpdir && -e $src) {
	@msg = split(' ', `iden -n -KP $temp $src`);
    }
    if (!@msg || $msg[2]) {
	unless ($quiet) {
	    print STDERR $src;
	    print STDERR "\t[", $nn, "]" if $nn;
	    print STDERR "\thas a problem!\n";
	    &System("iden -w50 -KP $temp $src");
	}
	return (1);
    } else {
	print STDERR $src, "\tis OK\n" unless ($quiet);
	return (0);
    }
}

sub external_reference() {
	my $revised;
	foreach my $sq (@seqlist) {
	    my $trst = $wrk . $sq;
	    my $trss = $pas . $sq;	# backward copy
	    if (!-w $trss) {
		&copy_or_link($trss, $trst);
		next;
	    }
	    my $cmd = "spaln -Q4 -a$ref_bka -n $trss";
	    print STDERR "$cmd\n" if ($vbs & 1);
	    my ($dist, @rest) = split(' ', `$cmd`);
	    next if (!@rest);
	    my $rn = pop(@rest);
	    my $ref_f = $pas . $rn;
	    my $cf = &conf($sq, $ref_f, 0);
	    if ($cf == 1) {
		$revised .= "$sq ";
	    } elsif ($cf == 2) {
		next;
	    }
	    my $cf = &conf($sq, $ref_f, 0);
	    if ($cf == 1) {
		$revised .= "$sq ";
	    } elsif ($cf == 2) {
		next;
	    }
	}
	return ($revised)
}
	    
sub mostsim {
	my ($dels, $rel, $unc, $arg) = @_;
	my %rlist = ();
	my $revised;
	foreach my $sq (@$rel) {
	    $rlist{$sq} = 1;
	    my $trst = $wrk . $sq;
	    my $trss = $pas . $sq;
	    &copy_or_link($trss, $trst);
	    $tmpfiles{$trst}++ if ($wrk eq $tmpdir);
	}
	foreach my $sq (@$unc, @$arg) {
	    next if $rlist{$sq};
	    my $trst = $wrk . $sq;
	    my $trss = $pas . $sq;	# backward copy
	    &copy_or_link($trst, $trss) if (!-e $trss && -e $trst);
	    if ($cogonly) {	# skip if remedy mode
		my @a = &exin($trss, 'g');
		my $t = rindex($a[0], '/') + 1;
	        my $trst = $wrk . $sq;
	        my $trss = $pas . $sq;
	        next unless (&copy_or_link($trss, $trst));
		next if ($a[6] eq 'X');
		my $s = substr($a[1], 0, 1) eq '$'? 1: 0;
		next unless (substr($a[0], $t, 8) eq substr($a[1], $s, 8));
	    }
	    my $tlm = (index($$dels{$sq}, 'n') >= 0);
	    my $trm = (index($$dels{$sq}, 'c') >= 0);
	    my $ref_f = "./$reference";
	    my $rn = $reference;
	    unless ($prof) {
#		my $cmd = "qdiv -f$reference -s$pas $sq | sort -g | head -1";
		my $cmd = "aln -KP -ig:$reference -s$pas -n $sq | sort -g | head -1";
		print STDERR "$cmd\n" if ($vbs & 1);
		my ($dist, @rest) = split(' ', `$cmd`);
		my $tn = pop(@rest);
		$rn = pop(@rest);
		$ref_f = $pas . $rn;
	    }
	    my $cf = &conf($sq, $ref_f, 0, $tlm, $trm);
	    if ($cf == 1) {
		$revised .= "$sq ";
	    } elsif ($cf == 2) {
		next;
	    }
	    next if ($quiet > 1);	# report tested seqs
	    my $tn = $wrk . $sq;
	    my @a = split(' ', `utp -n $tn`);
	    my $len = pop(@a);
	    my $peudo = join(' ', &exin($tn, 'd')) . "\n";
	    my $nx = 0;
	    open(UTP, "utp -c $tn |") or die "Can't run utp \n";
	    while (<UTP>) {
		if (/X = (\d+)/) {
		    $nx = $1;	# No of Xs
		    last;
		}
	    }
	    close (UTP);
	    my $cmd = "aln -yJ10 -s$wrk $ref_f $sq";
	    print STDERR "$cmd\n" if ($vbs);
	    open(ALN, "$cmd |") or die "Can't run $cmd !\n";
	    my ($m, $npos, $nint, $ae, $be);
	    my ($p, $g, $u);
	    while (<ALN>) {
		if (/^;B/) {
		   ($npos, $nint) = split;
		} elsif (/^;b/) {
		    my @a = split;
		    for (shift(@a); @a; ) {
			shift(@a);
			my $num = shift(@a);
			$m++ if (chomp($num) == 2);
		    }
		} elsif (/^;m/) {
		    my @a = split;
		    for (shift(@a); @a; ) {
			if (shift(@a) == 1)	{++$ae;}
			else			{++$be;}
		    }
		} elsif (/^Score/) {
		    my @a = split;
		    $p = substr($a[12], 1);
		    $g = $a[8];
		    $u = $a[10];
		} elsif (/^ALIGNMENT/) {
		    last;
		}
	    }
	    close(ALN);
	    printf("%-16s %-16s %7d %3d %3d %3d %3d %3d %5.2f %s", 
		$rn, $sq, $len, $nx, $m, $be, $g, $u, $p, $peudo);
	}
	return ($revised);
}

sub rmtemp {
	foreach my $fn (keys(%tmpfiles)) {
	    unlink($fn);
	}
	exit (0);
}

sub toolongorshort {
    my ($msa, $seqs, $ndefs, $cdefs) = @_;
    my @len = ();	# seq. length
    my $cmd = "utp -l2 $msa";
    print STDERR "$cmd\n" if ($vbs & 1);
    open(UTP, "$cmd |") or die "Can't run $cmd !\n";
    while (<UTP>) {
	my ($sd, $ln) = split;
	push(@$seqs, $sd);
	push(@len, $ln);
    }
    close(UTP);
    foreach my $ol (&Dixon(0.01, \@len)) {
	$$cdefs[$ol] .= 'l';
	$$ndefs[$ol]++;
    }
    my ($n, $a, $v) = (0, 0, 0);
    my ($nn, $av, $vr) = (0, 0, 0);
    foreach (0 .. $#len) {
	++$n;
	$a += $len[$_];
	$v += $len[$_] * $len[$_];
	unless ($$ndefs[$_]) {
	    ++$nn;
	    $av += $len[$_];
	    $vr += $len[$_] * $len[$_];
	}
    }
    if ($nn <= 2) {
	$nn = $n; $av = $a; $vr = $v;
    }
    if ($nn == 1) {$vr = sqrt($av);}
    elsif ($nn > 1) {
	$av /= $nn;
	$vr = sqrt(($vr - $nn * $av * $av) / ($nn - 1));
    }
#    $min_plen = $av - 3 * $vr;
#    $max_plen = $av + 3 * $vr;
    $min_plen = $av / 2;
    $max_plen = $av * 2;
}

sub pseudogene {	# unfinished or putative pseudogene
    my ($seqs, $ndefs, $cdefs) = @_;
    for (my $no = 0; $no < @$seqs; ++$no) {
	my $sq = $$seqs[$no];
	my $fnm = "$pas$sq";
	next unless (-w $fnm);	# skip if write-protected
	my ($fn, $chr, $gl, $gr, $sns, $no_exon, $o15, @def) = &exin($fnm, 'g');
	next unless ($o15 eq 'N' || $o15 eq 'Q');
	foreach (@def) {
	    if (/C/) {$$cdefs[$no] .= 'C'; $$ndefs[$no] += $defectP;}
	    if (/D/) {$$cdefs[$no] .= 'D'; $$ndefs[$no] += $defectP;}
	    if (/I/) {$$cdefs[$no] .= 'I'; $$ndefs[$no] += $defectP;}
	    if (/M/) {$$cdefs[$no] .= 'M'; $$ndefs[$no] += $defectP;}
	    if (/O/) {$$cdefs[$no] .= 'O'; $$ndefs[$no] += $defectP;}
	    if (/T/) {$$cdefs[$no] .= 'T'; $$ndefs[$no] += $defectP;}
	}
	next unless ($in_genetic);
	my $cmd = &gnmdb($chr);
	unless ($cmd) {
	    print STDERR "Warning: $chr not founs !\n" if ($vbs & 4);
	    next;
	}
	$chr = '$' . $chr if (substr($chr, 0, 1) ne '$');
	$cmd = "utn -c -A1 $cmd'$chr $gl $gr'";
	open (UTN, "$cmd |") or die "Can't run utn !\n";
	while (<UTN>) {
	    if (/^N/) {
		my @a = split;
		if (pop(@a) >= $AmbN_thr) {$$cdefs[$no] .= 'N'; $$ndefs[$no] += $defectP;}
		last;
	    }
	}
	close(UTN);
    }
}

sub processed_pg {	# possible processed pseudogene
    my ($msa, $ndefs, $cdefs) = @_;
    my %genspc;
    my @introns;
    my @members;
    my $cmd = "utp -B4 $msa";
    open(UTP, "$cmd |") or die "Can't run $cmd !\n";
    while (<UTP>) {
	my ($mem, $no, $ni, $nu) = split;
	push(@introns, $ni);		# # introns
	push(@members, $mem);		# seq name
	$genspc{substr($mem, 0, 8)} += $ni;
    }
    my $with = 0;
    foreach (keys(%genspc)) {
	$with += ($genspc{$_}? +1: -1);
    }
    for (my $no = 0; $no < @introns; ++$no) {
	next if ($introns[$no]);	# has introns
	my $gs = substr($members[$no], 0, 8);
	if ($with > 0 || $genspc{$gs}) {
	    $$cdefs[$no] .= 'c';	# paralog with introns
	    $$ndefs[$no] += $defectP;
	}
    }
}

sub eijpattern {
	my ($msa, $seqs, $sameij, $lone) = @_;
	my %pat = ();
	my $nint = 0;

	my $cmd = "utp -B6 $msa |";
	print STDERR "$cmd\n" if ($vbs & 1);
	open(UTP, $cmd) or die "Can't run $cmd !\n";
	while (<UTP>) {
	    chop;
	    my ($name, $no, $ni, $li, $cln, $bp) = split(' ', $_, 6);
	    --$no;
	    $nint += $ni;
	    $pat{$bp} .= "$no ";
	    $$lone[$no] = $li;
	}
	close(UTP);
	my $maxns = 0;
	foreach my $bp (keys(%pat)) {
	    my @a = split(' ', $pat{$bp});
	    my $ns = @a;
	    $maxns = $ns if ($ns > $maxns);
	    foreach my $no (@a) {
		$$sameij[$no] = $ns;
	    }
	}
	return ($nint / @$seqs, $maxns);
}

sub outliers {
	my ($msa, $ndefs, $cdefs) = @_;
	my $cmd = "prrn5 -O2 -I0 -c$alpha $msa";
	print STDERR "$cmd\n" if ($vbs & 1);
	unless (open(PRRN, "$cmd |")) {
	    print STDERR "Can't run $cmd !\n";
	    return;
	}
	my $no = 0;
	while (<PRRN>) {
	    my @a = split;
	    for (my $i = 4; $i <= 10; ++$i) {
		next unless ($a[$i]);
		$$ndefs[$no] += $a[$i];			# total # of outliers
		if ($i == 6) {$$cdefs[$no] .= 'n';}	# N-terminal deletion
		elsif ($i == 10) {$$cdefs[$no] .= 'c';}	# C-terminal deletion
		else {$$cdefs[$no] .= 'd';}		# others
	    }
	    ++$no;
	}
	close(PRRN);
	return if ($no < 6);
	$cmd = "utp -l4 $msa";
	print STDERR "$cmd\n" if ($vbs & 1);
	unless (open(UTP, "$cmd |")) {
	    print STDERR "Can't run $cmd !\n";
	    return;
	}
	my $no = 0;
	my @idl = ();
	my @avr = ();
	my @var = ();
	while (<UTP>) {
	    my ($sid, @a) = split;
	    for (my $i = 0; $i < @a; ++$i) {
		$idl[$no][$i] = $a[$i];
		$avr[$i] += $a[$i];
		$var[$i] += $a[$i] * $a[$i];
	    }
	    ++$no;
	}
	my @llt = ();
	my @ult = ();
	for (my $i = 0; $i < 3; ++$i) {
	    $avr[$i] /= $no;
	    $var[$i] = sqrt(($var[$i] - $no * $avr[$i]) / ($no - 1));
	    $llt[$i] = $avr[$i] - $lendev * $var[$i];
	    $ult[$i] = $avr[$i] + $lendev * $var[$i];
	}
	my $np = 0; my @pavr = (); my @pvar = ();
	my $nn = 0; my @navr = (); my @nvar = ();
	for (my $n = 0; $n < $no; ++$n) {
	    for (my $i = 0; $i < 3; $i += 2) {
		if ($idl[$n][$i] >= $avr[$i]) {
		    $pavr[$i] += $idl[$n][$i];
		    $pvar[$i] += $idl[$n][$i] * $idl[$n][$i];
		    ++$np;
		} else {
		    $navr[$i] += $idl[$n][$i];
		    $nvar[$i] += $idl[$n][$i] * $idl[$n][$i];
		    ++$nn;
		}
	    }
	}
	for (my $i = 0; $i < 3; $i += 2) {
	    next if ($np < 3 || $nn < 3);
	    $pavr[$i] /= $np;
	    $pvar[$i] = sqrt(($pvar[$i] - $np * $pavr[$i]) / ($np - 1));
	    $navr[$i] /= $nn;
	    $nvar[$i] = sqrt(($nvar[$i] - $nn * $navr[$i]) / ($nn - 1));
	    $llt[$i] = $pavr[$i] - $lendev * $pvar[$i]
		if (($pavr[$i] - $navr[$i]) > $lendev * sqrt($pvar[$i] * $nvar[$i]));
	}
	for (my $n = 0; $n < $no; ++$n) {
	    if (!($$cdefs[$n] =~ /n/) && ($llt[0] > $idl[$n][0] || $idl[$n][0] > $ult[0]))
		{++$$ndefs[$n]; $$cdefs[$n] .= 'n';} 
	    if (!($$cdefs[$n] =~ /d/) && ($llt[1] > $idl[$n][1] || $idl[$n][1] > $ult[1]))
		{++$$ndefs[$n]; $$cdefs[$n] .= 'd';} 
	    if (!($$cdefs[$n] =~ /c/) && ($llt[2] > $idl[$n][2] || $idl[$n][2] > $ult[2]))
		{++$$ndefs[$n]; $$cdefs[$n] .= 'c';} 
	}
}
