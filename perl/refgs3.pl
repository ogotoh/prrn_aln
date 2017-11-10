#!/usr/bin/perl

##############################################################################
#
#	Refgs3.pl: Refinement of gene prediction with 3 template selection modes
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
#	Usage:	refgs3.pl [-m...] [-e] [-HN] [-sSrc] upg 
# 	or	refgs3.pl [-f1] [-m...] [-e] [-HN] [-sSrc] MSA
#
#	The aa sequences in extended fasta format before refinement must be 
#	in './prd' directory
#
#	The default set of template selection modes is CR-M1-PR
#
#############################################################################

use lib "$ENV{HOME}/perl";      # perl module
use Util;
use Exin;
use strict;

my @body = ("UNQ", "FST", "SND", "TRD");
my @pcap = ("Unq", "Fst", "Snd", "Trd");
my @dirs = ("prd", "prd1", "prd2", "prd3");
my @popt = ("", "-R", "-R1");
my $ordrstr = "102";
my $refgsopt = "-E -I5 -q2";

my $Debug;
my $samegr_thr = 100;
my $from_which = 0;
my $vrbs = 0;
my $redundant = "redundant.txt";
my $equalmem;
my $keep_im_file = 1;	# 0:no if; 1: keep 1st; 2: keep body; 3: keep all
my $pas = $dirs[0];
my $max_allowed_unp = 5;
my $referencelst;
my %reference = ();

while ($_ = $ARGV[0], /^-/) {
	shift;
	my $each;
	if (/^-e/)	{$equalmem = 1;}	# equalize the sizes of UNQ, FST, SND, and TRD
	elsif (/^-f(\S*)/) {			# startf from phase from_which
	    &Util::getoptarg(\$from_which, $1);
	} elsif (/^(-k\S*)/) {			# keep psuegogenes
	    &Util::getoptarg(\$each, $1);
	    $refgsopt .= $each if ($each);
	} elsif (/^-m(\S*)/) {			# a series of template selection modes
	     &Util::getoptarg(\$ordrstr, $1);
	} elsif (/^(-n\S*)/) {			# exam presence of 'N' in gene
	    &Util::getoptarg(\$each, $1);
	    $refgsopt .= " $each" if ($each);
	} elsif (/^-s(\S*)/) {			# directory wherein originals
	     &Util::getoptarg(\$pas, $1);
	} elsif (/^-D(\S*)/) {
	     &Util::getoptarg(\$Debug, $1, 1);
	} elsif (/^-H(\S*)/) {			# threshold of overlap range (nt)
	    &Util::getoptarg(\$samegr_thr, $1);
	} elsif (/^(-I\S*)/) {			# number of iterations
	    &Util::getoptarg(\$each, $1);
	    $refgsopt .= " $each" if ($each);
	} elsif (/-K(\S*)/) {
	    &Util::getoptarg(\$keep_im_file, $1, 2);	# keep intermediate files
	} elsif (/^-R(\S*)/) {			# external reference
	    &Util::getoptarg(\$referencelst, $1);
	} elsif (/^-V(\S*)/) {
	     &Util::getoptarg(\$vrbs, $1, 1);
	}
}

my $infile = shift;
my $dot = rindex($infile, '.');
my $ext = ($dot > 0)? substr($infile, $dot): "";
my $ref_db;

if ($referencelst) {
	open(REF, $referencelst) or die "Can't open $referencelst !\n";
	while (<REF>) {
	    chomp;
	    ++$reference{$_};
	}
	close(REF);
}

my @ordr = ();
for (my $i = 0; $i < length($ordrstr); ++$i) {
	my $c = substr($ordrstr, $i, 1);
	next unless ($c =~ /\d/);
	push(@ordr, $c + 0);
}

$dirs[3] = $dirs[1] if ($keep_im_file < 3);		# write back

unless (-d $pas) {die "Error! No $pas directory!\n";}
$pas .= '/' if (substr($pas, -1, 1) ne '/');

my $err = 0;
for (my $i = 1; $i <= @ordr; ++$i) {
	if (-e $dirs[$i] && -f $dirs[$i]) {
	    print STDERR "$dirs[$i] is an oridnay file !\n";
	    ++$err;
	} else {
	    mkdir $dirs[$i] unless (-d $dirs[$i]);
	}
}
exit(1) if ($err);

my @smemb = ();
my %nmemb = ();
my @sleng = ();

if ($from_which == 0) {

##### construct MSA with the given guide tree #####

	my $msa = "MSA" .  $ext;
	my $unq = $body[0] . $ext;
	unless (-s $msa) {
	    my $cmd = "upg2prrn.pl -s$pas -O $infile";
	    Util::System($cmd, $vrbs);
	}
	exit(0) if (lc($ordrstr) eq 'm');

	my $cmd = "utp -l3 $msa |";
	print STDERR $cmd, "\n" if ($vrbs);
	open(UTP, $cmd) or die "Can't run $cmd !\n";
	while (<UTP>) {
	    my ($no, $sn, $ln) = split;
	    $nmemb{$sn} = $no;
	    push(@smemb, $sn);
	    push(@sleng, $ln);
	}
	close(UTP);

##### make reference db if %reference #####

	my @firm = ();
	&extendedseed(\@smemb, \%reference, \@firm) if (keys(%reference));
	&make_bka(\@firm) if (@firm);		# custom reference

##### remove redundancy #####

	my $memb = join(' ', @smemb);
	my $lsttmpf = "samegr.lst$$";
	$cmd = "samegr.pl -F1 -H$samegr_thr -s$pas $memb > $lsttmpf";
	Util::System($cmd, $vrbs);
	if (-s $lsttmpf) {
	    open(REDUNDANT, ">> $redundant") or 
		print STDERR "Cant't write to $redundant !\n";
	    &rmsamegr($msa, $unq, \%nmemb, $lsttmpf);
	    close(REDUNDANT);
	}
	link $msa, $unq unless (-s $unq);
	unlink($lsttmpf);
	++$from_which;
}

my $msalast = $body[$#body] . $ext;
for (my $phs = $from_which; $phs <= @ordr; ++$phs) {
	my $now = $phs - 1;
	my $opt = $popt[$ordr[$now]];
	$opt .= " -m" if ($now == 0);	# copy/link protected seqs.
	my $msanow = $body[$now] . $ext;
	my @a = split(' ', `utp -n $msanow`);
	if (!@a || $a[2] < 3) {
	    rename $msanow, $msalast;
	    last;
	}
	my $cmd = "refgs.pl -s$dirs[$now] -t$dirs[$phs] $opt $refgsopt $msanow";
	Util::System($cmd, $vrbs);
	if (-s "$msanow.1") {
	    $cmd = "msa2rev.pl -r$body[$phs] $msanow.*";
	    Util::System($cmd, $vrbs);
	    unlink $msanow if ($keep_im_file == 0 || ($keep_im_file == 1 && $now));
	    Util::System("rm $msanow.*", $vrbs) if ($keep_im_file < 2);
	} else {
	    if ($keep_im_file < 2) {
		rename $msanow, "$body[$phs]$ext";
	    } else {
		link $msanow, "$body[$phs]$ext";
	    }
	}
	if ($keep_im_file < 2) {	# remove intermediates
	    if ($now == 0) {
		my $msa = "MSA" .  $ext;
		my $unq = $body[0] . $ext;
		if (-e $unq) {
		    unlink $unq unless (`diff -q $msa $unq`);
		}
	    } elsif ($dirs[$now] ne $dirs[1]) {
		foreach my $mem (@smemb) {
		    unlink("$dirs[$now]/$mem");
		}
	    }
	}
}
&delete_bka if ($ref_db);

##### unify the sizes of UNQ, FTS, SND, and TRD #####

&samememb("$body[@ordr]$ext") if ($equalmem);

##### end of main #####

##### make and delete reference.bka #####

sub make_bka {
	$ref_db = "reference" . $ext;
	$ref_db =~ tr/./_/; 
	$popt[1] .= " -a$ref_db";
	my $ref = shift;
	my $faa = $ref_db . ".faa";
	unlink($faa) if (-e $faa);
	foreach (@$ref) {
	    my $fn = $pas . $_;
	    my $cmd = "cat $fn >> $faa";
	    Util::System($cmd, $vrbs, 1);
	}
	my $cmd = "makeidx.pl -a $faa";
	Util::System($cmd, $vrbs);
}

sub delete_bka {
	my $cmd = "rm $ref_db.*";
	Util::System($cmd, $vrbs);
}

#### remove minor isoforms #####

sub rmsamegr {
	my $msa = shift;
	my $unq = shift;
	my $nmemb = shift;
	my $lsttmpf = shift;
	my @rmv = ();
	my $cmd = "clade.pl -R -H0 $lsttmpf |";
	open(CLD, $cmd) or return;
	while (<CLD>) {
	    my @srmlist = split;
	    my @nrmlist = ();
	    my @scr = ();
	    foreach my $mem (@srmlist) {
		push(@nrmlist, $$nmemb{$mem});
		my @a = split(' ', `head -1 $pas/$mem`);
		push(@scr, $a[-1]) if ($a[-2] eq 'N' || $a[-2] eq 'Q');
	    }
	    unless (@scr == @srmlist) {
	        my $rml = join(' ', @nrmlist);
	        my $tmp = "tmp$$";
	        my $cmd = "rdn -ced $msa $rml > $tmp";
	        Util::System($cmd, $vrbs, 1);
	        @scr = ();
	        foreach my $mem (@srmlist) {
		    my @aln = split(' ', `aln -n $tmp $pas/$mem`);
		    splice(@aln, 1, 2) if (@aln == 10);
		    push(@scr, $aln[1]);
	        }
	        unlink($tmp);
	    }
	    my $maxi = 0;		# highest score
	    for (my $i = 1; $i <= $#scr; ++$i) {
		$maxi = $i if ($scr[$i] > $scr[$maxi]);
	    }
	    my ($master) = splice(@srmlist, $maxi, 1);	# retain one
	    push(@rmv, @srmlist);
	    foreach (@rmv) {
		print REDUNDANT $_,"\t",$master,"\n";
	    }
	}
	close(CLD);
	my $srmv = "";
	foreach (@rmv) {
	    $srmv .= "$$nmemb{$_} ";
	}
	my $cmd = "rdn -ced $msa $srmv > $unq";
	Util::System($cmd, $vrbs);
}

sub samememb {
	my $trd = shift;
	my $dot = index($trd, '.');
	my $ext = substr($trd, $dot);
	my %rsv = ();
	my $ntrd = 0;
	open(UTP, "utp -l0 $trd |") or die "Can't run utp !\n";
	while (<UTP>) {
	    chop;
	    $rsv{$_} = 1;
	    ++$ntrd;
	}
	close(UTP);

	for (my $i = 0; $i < $#body; ++$i) {
	    my $org = $body[$i] . $ext;
	    last if ($org eq $trd);
	    my $sel = $pcap[$i] . $ext;
	    my @mems = ();
	    my $norg = 0;
	    open(UTP, "utp -l1 $org |") or die "Can't run utp !\n";
	    while (<UTP>) {
	        chop;
	        my ($no, $mem) = split;
	        push(@mems, $no) if ($rsv{$mem});
	        ++$norg;
	    }
	    close(UTP);
	    if ($norg > $ntrd) {
	        my $slist = join(' ', @mems);
	        my $cmd = "rdn -csd $org $slist > $sel";
	        Util::System($cmd, $vrbs);
	    } else {
		link($org, $sel);
	    }
	}
}

# change to write-protected if 
# same gene organization no gap
# compared with an already protected sequence

sub extendedseed {
	my ($seqs, $seed, $firm) = @_;
	my @seedpat = ();
	my %othrpat = ();
	foreach my $mem (@$seqs) {
	    my $fn = $pas . $mem;
	    my @a = &exin($fn, 'e');
	    shift(@a);
	    my $pat = join(' ', @a);
	    if ($$seed{$mem}) {
		push(@seedpat, $pat);
	    } else {
		$othrpat{$mem} = $pat;
	    }
	}
	foreach my $sdpat (@seedpat) {
	    foreach my $mem (keys(%othrpat)) {
		if ($sdpat eq $othrpat{$mem}) {			# same gene structure
		    my $fn = $pas . $mem;
		    chmod 0444, $fn if (-w $fn && !&exin($fn, 'd'));		# no defect
		    push(@$firm, $mem) if (!-w $fn);
		}
	    }
	}
}

sub extendedseed_1 {
	my ($msa, $slen, $seed) = @_;
	my $cmd = "utp -B11 $msa |";
	open (DVN, $cmd) or die "Can't run $cmd !\n";
	my ($i, $j) = (0, 1);
	while (<DVN>) {
	    my @a = split;
	    my $dlen = $$slen[$i] - $$slen[$j];
	    next if (shift(@a) > 0 || abs($dlen) > $max_allowed_unp);
	    if ($$seed{$a[3]}) {
		my $fn = "$pas/$a[4]";
		my @exlng = &exin($fn, 'd');
		chmod 0444, $fn if (-w $fn && !@exlng);
	    }
	    if ($$seed{$a[4]}) {
		my $fn = "$pas/$a[3]";
		my @exlng = &exin($fn, 'd');
		chmod 0444, $fn if (-w $fn && !@exlng);
	    }
	} continue {
	    if (++$i == $j) {
		++$j; $i = 0;
	    }
	}
	close(DVN);
}

