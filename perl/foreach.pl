#!/usr/bin/perl

##############################################################################
#
#       Execute command in parallel
#
#
#       Osamu Gotoh, ph.D.      (-2001)
#       Saitama Cancer Center Research Institute
#       818 Komuro, Ina-machi, Saitama 362-0806, Japan
#
#       Osamu Gotoh, Ph.D.      (2001-)
#       National Institute of Advanced Industrial Science and Technology
#       Computational Biology Research Center (CBRC)
#       2-41-6 Aomi, Koutou-ku, Tokyo 135-0064, Japan
#
#       Osamu Gotoh, Ph.D.      (2003-2012)
#       Department of Intelligence Science and Technology
#       Graduate School of Informatics, Kyoto University
#       Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
#
#       Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
#___________________________________________________________________________
#
#	Usage: foreach.pl [-tN] [-ppost] [-oS] command [-command options]
#	Options:
#	-t#	Number of threads
#	-pporst	Commnad for post process
#	-oS	output file name (_n is added for each input)
#
#############################################################################

use lib "$ENV{HOME}/perl";      # perl module
use strict;
use threads;
use Thread::Queue;

my $opt;
my $vbs = 1;
my $n_thread = 1;
my $out;
my $postarg;

sub System {
	my $cmd = shift;
	$cmd .= " $postarg" if ($postarg);
	print STDERR "$cmd\n" if ($vbs & 2);
	if (($vbs & 1) && system($cmd)) {
	    print STDERR "$cmd failed !\n";
	}
}

while ($_ = $ARGV[0], /-/) {
	shift;
	if (/^-V(\S*)/) {$vbs =  ($1 ne "")? $1: 3;}
	elsif (/^-D/) {$vbs = 2;}
	elsif (/^-t(\d+)/) {$n_thread = $1;}
	elsif (/^-o(\S*)/) {$out = $1? $1: "";}
	elsif (/^-p(.+)/) {$postarg = $1;}
}

my $cmd = shift;
$out = $cmd if (defined($out) && !$out);

while ($_ = $ARGV[0], /^-/) {
	shift;
	$cmd .= " $_";
}

print STDERR $cmd, "\n";

if ($n_thread < 2) {
	while ($_ = $ARGV[0]) {
	    shift;
	    &System("$cmd $_");
	}
	exit (0);
}
	
my $queue = new Thread::Queue;

while ($_ = $ARGV[0]) {
	shift;
	$queue->enqueue($_);
}

my @threads;
foreach (1 .. $n_thread) {
	my $thread = threads->new(\&my_thread, $_);
	push(@threads, $thread);
	$queue->enqueue(undef);
}

foreach (@threads) {
	$_->join;
}

sub my_thread {
        my $i = shift;
	while (my $q = $queue->dequeue) {
	    my $com = "$cmd $q";
	    $com .= "> out_$i" if ($out);
	    &System($com);
	    threads->yield();
	}
	return ($i);
}
