#!/usr/bin/perl

#	Input	.dvn/.dvp
#	output	< Thr clades
#	Option	-M (don't include Master) -S (Master only)

my $Thr = 100;
my $Mstr = 1;
my $Slav = 1;
my $Uniq = 1;
my $Sim = 0;
my $avoid;
my %member;
my %master;
my %subsid;

while ($_ = $ARGV[0], /^-/) {
	shift;
	last if /^-$/;
	/^-D/ && ($Debug = 1);
	/^-H(\d+)/ && ($Thr = $1);
	/^-A(\S+)/ && ($avoid = $1);
	/^-M/ && ($Mstr = 0);	# Slave only, Don't output Master
	/^-S/ && ($Slav = 0);	# Master only, Don't output Slave
	/^-U/ && ($Uniq = 0);	# Don't output unique members
	/^-R/ && ($Sim = 1);	# Resemblance
}

while (<>) {
	my (@a) = split;
	my ($scr) = shift(@a);
	my ($b) = pop(@a);
	my ($a) = pop(@a);
	if ($Sim && $scr < $Thr || !$Sim && $scr > $Thr) {
	    $member{$a} = 1;
	    $member{$b} = 1;
	    next;
	}
	if ($avoid && $a =~ /$avoid/ && !($b =~ /$avoid/)) {
	    my($tmp) = $a;
	    $a = $b;
	    $b = $tmp;
	}
	unless ($subsid{$b}) {
	    $subsid{$a} = $a unless $subsid{$a};
	    $mstr = $subsid{$b} = $subsid{$a};
	    $master{$mstr} .= $b . " ";
	}
}

sub bymemb {
	($master{$b} =~ tr/ / /) <=> ($master{$a} =~ tr/ / /);
}

foreach $mstr (sort bymemb keys(%master)) {
	print $mstr," " if ($Mstr);
	print $master{$mstr} if ($Slav);
	print "\n";
}

exit(0) unless ($Uniq);

foreach $mem (sort keys(%member)) {
	print $mem,"\n" unless($subsid{$mem});
}

