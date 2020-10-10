#!/usr/bin/perl 

##############################################################################
#
#       anno2gsiseq.pl	Last updated <2019-10-31>
#
#	Convert transcript sequences in FASTA format into 
#	gene-structure-informed (extended FASTA) format
#       based on gff/gtf and sequence files provided by Genome projects
#
#
#	Osamu Gotoh, Ph.D.      (2001-)
#	National Institute of Advanced Industrial Science and Technology
#	Artificial Intelligence Research Center (AIRC)
#	2-4-7 Aomi, Koutou-ku, Tokyo 135-0064, Japan
#
#	Osamu Gotoh, Ph.D.      (2003-)
#	Department of Intelligence Science and Technology
#	Graduate School of Informatics, Kyoto University
#	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
#
#	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
#
#___________________________________________________________________________
#
# 	Usage: 
#
#	anno2gsiseq.pl -f GffStyle [-KD] [-o out_file] gff_file(.gz) sequence_file(.gz)
#
##############################################################################

use strict;

my @style_list = ('NONE', 'ALN','BROAD', 'EMBL', 'EUPATH', 'HYMGD',
	'INSECTBASE','I5K', 'JGI', 'KAZUSA', 'NCBI', 'PHYTOZOME','WORMBASE');
my %style_list;

foreach my $i (0..$#style_list) {
	$style_list{$style_list[$i]} = $i;
}

sub usage {
	print "Usage: anno2gsiseq.pl -f GffStyle [-KD]  [-o out_file] gff_file(.gz) sequence_file(.gz)\n";
	print "\tSequence_file shold be in multi-fasta format, each entry of which corresponds to\n";
        print "\tthe pid (protein-id) or tid (transcript-id) speceified in the gff_file.\n";
	print "\nOptions: (default)\n";
	print "\t-f GffStyle \tSupported GffStyles:";
	foreach my $i (1..$#style_list) {
	    print ("\n\t\t") if (($i - 1) % 6 == 0);
	    print "$style_list[$i] ";
	}
	print "\n" if (@style_list % 6);
	print "\t-KD: nucleotide (aa) sequence\n";
	print "\t-l linewidth (60)\n";
	print "\t-o out_file (stdout)\n";
	exit(1);
}

my $outfile = "";
my $gff_style;
my $separator = '_';	
my $linewidth = 60;
my $gid_field = 0;
my $pid_term;
my $verb;
my $isaa = 1;

while ($_ = $ARGV[0], /^-/) {
	shift;
	/^-f(\S+)/	&& ($gff_style = uc($1));
	/^-f$/		&& ($gff_style = uc(shift));
	/^-h/		&& &usage;
	/^-KD/		&& ($isaa = 0);
	/^-l(\d+)/	&& ($linewidth = $1);
	/^-l$/		&& ($linewidth = shift);
	/^-o(\S+)/	&& ($outfile = $1);
	/^-o$/		&& ($outfile = shift);
	/^-p(\S+)/	&& ($pid_term = $1);
	/^-p$/		&& ($pid_term = shift);
	/^-s(\S+)/	&& ($separator = $1);
	/^-s$/		&& ($separator = shift);
	/^^V/		&& ($verb = 1);
}

&usage unless ($style_list{$gff_style});
my $gff_file = shift;
&usage unless ($gff_file);
my $seq_file = shift;
&usage unless ($seq_file);

$separator = '.'
    if ($gff_style eq 'WORMBASE' || $gff_style eq 'INSECTBASE' ||
	$gff_style eq 'EMBL' || $gff_style eq 'I5K');
if ($isaa) {
	$pid_term = 'CDS' unless (defined($pid_term));
	$pid_term = 'cds' if ($gff_style eq 'ALN');
	$pid_term = 'mRNA' if ($gff_style eq 'I5K');
	$pid_term = 'gene' if ($gff_style eq 'HYMGD');
} else {
	$pid_term = 'exon' unless (defined($pid_term));
}

my $genome_id = substr($gff_file, rindex($gff_file, '/') + 1);
my $sep = index($genome_id, $separator);
$genome_id = substr($genome_id, 0, $sep) if ($sep > 0);

if ($outfile) {
	open(OFD, "> $outfile") or die "Can't write to $outfile !\n";
} else {
	open(OFD, ">-") or die "Can't write to STDOUT !\n";
}

#####################################################################
#
#	Main

my %keyh = ();
&readgff($gff_file, \%keyh);
&read_write($seq_file, \%keyh);
close(OFD);
exit(0);

#	End of Main
#
#####################################################################

sub getpid {
	my ($rest) = @_;
	if ($gff_style eq 'ALN') {
	    if ($rest =~ /Name=(\S+)/) {
		my $pid = $1;
		return ($pid =~ tr/_/./);
	    } else {
		return (undef);
	    }
	} elsif ($gff_style eq 'BROAD') {
	    return ($rest =~ /transcript_id\s*(\S+)/)? $1: undef;
	} elsif ($gff_style eq 'EMBL') {
	    if ($rest =~ /ID=CDS\s+(\S+)/) {return $1;}
	    if ($rest =~ /(PAC \d+)/) {return $1;}
	    if ($rest =~ /ID=(\S+)/) {return $1;}
	    if ($rest =~ /protein_id=(\S+)/) {return $1;}
	    return ($rest =~ /protein_id\s+(\S+)/)? $1: undef;
	} elsif ($gff_style eq 'I5K') {
	    if ($rest =~ /ID=(\S+)/) {return $1;}
	    return undef;
	} elsif ($gff_style eq 'JCVI') {
	    if ($rest =~ /ID=(\S+)/) {
		$rest = $1;
		if ($rest =~ /(\S+)_cds/) {$rest = $1};
		return ($rest);
	    } else {
		return undef;
	    }
	} elsif ($gff_style eq 'JGI') {
	    if ($rest =~ /proteinId=(\d+)/) {return ($1);}
	    elsif ($rest =~ /protein\S+\s+(\S+)/) {return ($1);}
	    return (undef);
	} elsif ($gff_style eq 'KAZUSA') {
	    return ($rest =~ /ID= (\S+)/)? $1: undef;
	} elsif ($gff_style eq 'MISC') {
	    if ($rest =~ /Parent=(\S+)/) {
		$rest = $1;
		if ($rest =~ /(\S+).mRNA/) {$rest = $1};
		return ($rest);
	    } else {
		return $rest;
	    }
	} elsif ($gff_style eq 'NCBI') {
	    return ($rest =~ /protein_id=(\S+)/)? $1: undef;
	} elsif ($gff_style eq 'PHYTOZOME') {
	    return ($rest =~ /pacid=(\d+)/)? $1: undef;
	} elsif ($gff_style eq 'WORMBASE') {
	    if ($rest =~ /cds (\S+)/) {return $1;}
	    return ($rest =~ /CDS (\S+)/)? $1: undef;
	} elsif ($gff_style eq 'INSECTBASE') {
	    if ($rest =~ /TCOGS2 (\S+)/) {return $1;}	# Tribolium_castaneum
	    if ($rest =~ /ID=(\S+)/) {return $1;}
	    if ($rest =~ /Parent=(\S+)/)  {return $1;}
	    if ($rest =~ /Name=(\S+)/) {return $1;}
	} elsif ($gff_style eq 'EUPATH') {
	    if ($rest =~ /ID=cds_(\S+)/) {return $1;}
	} elsif ($gff_style eq 'HYMGD') {
	    if ($rest =~ /ID=(\S+)/) {return $1;}
	}
	return undef;
}

sub gettid {
	my ($rest) = @_;
	if ($gff_style eq 'ALN') {
	    if ($rest =~ /Name=(\S+)/) {
		my $pid = $1;
		return ($pid =~ tr/_/./);
	    } else {
		return (undef);
	    }
	} elsif ($gff_style eq 'BROAD') {
	    return ($rest =~ /transcript_id\s*(\S+)/)? $1: undef;
	} elsif ($gff_style eq 'EMBL') {
	    if ($rest =~ /Parent=transcript\s+(\S+)/) {return $1;}
	} elsif ($gff_style eq 'I5K') {
	    if ($rest =~ /ID=(\S+)/) {return $1;}
	    return undef;
	} elsif ($gff_style eq 'JCVI') {
	    if ($rest =~ /ID=(\S+)/) {
		$rest = $1;
		if ($rest =~ /(\S+)_cds/) {$rest = $1};
		return ($rest);
	    } else {
		return undef;
	    }
	} elsif ($gff_style eq 'JGI') {
	    if ($rest =~ /proteinId=(\d+)/) {return ($1);}
	    elsif ($rest =~ /protein\S+\s+(\S+)/) {return ($1);}
	    return (undef);
	} elsif ($gff_style eq 'KAZUSA') {
	    return ($rest =~ /ID= (\S+)/)? $1: undef;
	} elsif ($gff_style eq 'MISC') {
	    if ($rest =~ /Parent=(\S+)/) {
		$rest = $1;
		if ($rest =~ /(\S+).mRNA/) {$rest = $1};
		return ($rest);
	    } else {
		return $rest;
	    }
	} elsif ($gff_style eq 'NCBI') {
	    return ($rest =~ /transcript_id=(\S+)/)? $1: undef;
	} elsif ($gff_style eq 'PHYTOZOME') {
	    return ($rest =~ /pacid=(\d+)/)? $1: undef;
	} elsif ($gff_style eq 'WORMBASE') {
	    if ($rest =~ /cds (\S+)/) {return $1;}
	    return ($rest =~ /CDS (\S+)/)? $1: undef;
	} elsif ($gff_style eq 'INSECTBASE') {
	    if ($rest =~ /TCOGS2 (\S+)/) {return $1;}	# Tribolium_castaneum
	    if ($rest =~ /ID=(\S+)/) {return $1;}
	    if ($rest =~ /Parent=(\S+)/)  {return $1;}
	    if ($rest =~ /Name=(\S+)/) {return $1;}
	} elsif ($gff_style eq 'EUPATH') {
	    if ($rest =~ /ID=cds_(\S+)/) {return $1;}
	} elsif ($gff_style eq 'HYMGD') {
	    if ($rest =~ /ID=(\S+)/) {return $1;}
	}
	return undef;
}

sub readgff() {
	my ($gff, $keyh) = @_;
	my @pos = ();
	$pid_term = 'CDS' if ($gff_style eq 'JGI');
	if (substr($gff, -3, 3) eq '.gz') {
	    open (GFF, "gzip -cd $gff |") or die "CANNOT OPEN $gff\n";
	} else {
	    open(GFF, $gff) or die "CANNOT OPEN $gff\n";
	}
	my ($pchr, $psgn, $pid);
	while (<GFF>) {
	    if (/^#/) {
		$pid_term = 'gene' if ($gff_style eq 'JGI' && /version 3/);	# gff3
		next;
	    }
	    my ($chr, $dmy, $term, $st, $ed, $dot, $sgn, $phs, $rest) 
		= split(' ', $_, 9);
	    next unless ($term eq $pid_term);
	    $rest =~ tr/:;/  /;
	    if (my $cpid = $isaa? &getpid($rest): &gettid($rest)) {
		if ($pid eq $cpid) {
		    if ($st < $ed) {
			push(@pos, $st); push(@pos, $ed);
		    }
		} else {
		    $$keyh{$pid} = "$pchr $psgn " . join(' ', sort {$a <=> $b} @pos) if (@pos);
		    @pos = ();
		    $pchr =  $chr; $psgn = $sgn; 
		    push(@pos, $st); push(@pos, $ed);
		    $pid = $cpid;
		}
	    } elsif ($verb) {
		print STDERR $chr,' ',$rest;
	    }
	} 
	$$keyh{$pid} = "$pchr $psgn " . join(' ', @pos) if (@pos);
	close (GFF);
}

######## read the first seq #########

sub getpid_field {
	my ($seq, $keyh) = @_;
	if (substr($seq, -3, 3) eq '.gz') {
	    open(SEQ, "gzip -cd $seq |") or die "CANNOT OPEN $seq !\n";
	} else {
	    open(SEQ, $seq)  or die "CANNOT OPEN $seq !\n";
	}
	my $fld = -1;
LINE:	while (<SEQ>) {
	    next unless (/^>/);
	    tr/|:=/   /;
	    my @a = split;
	    $a[0] = substr($a[0], 1);
	    while (++$fld < @a) {
		last LINE if ($$keyh{$a[$fld]});
	    }
	    $fld = -1;
	}
	close(SEQ);
	return ($fld);
}

######## main subroutine #########

sub read_write
{
	my ($seq, $keyh) = @_;
	my $nn = 0;
	my $pid_field = &getpid_field($seq, $keyh);
	if (substr($seq, -3, 3) eq '.gz') {
	    open(SEQ, "gzip -cd $seq |") or die "CANNOT OPEN $seq !\n";
	} else {
	    open(SEQ, $seq)  or die "CANNOT OPEN $seq !\n";
	}
	while (<SEQ>) {
	    if (/^>/) {
		print OFD $_;
		++$nn;
		tr/|:=/   /;
		my @a = split;
		$a[0] = substr($a[0], 1) if ($pid_field <= 0);	# remove '>'
		my $pid = $pid_field >= 0? $a[$pid_field]: $nn;
		my ($coms, $sgn, @pos) = split(' ', $$keyh{$pid});
		&ncbiexon($sgn, \@pos) if (@pos >= 2);
	    } else {
		if (substr($_, -2, 1) eq '*') {
		    chop; chop;
		    $_ .= "\n";
		}
		print OFD unless (/^$/);
	    }
	}
	close(SEQ);
}

######## write extended FASTA #########

sub ncbiexon
{
	my ($sgn, $pos) = @_;
	my $line = ";C ";
	my $i = 0;
	my $n = $#$pos;
	if ($sgn eq '-') {
	    my $len = $$pos[1] - $$pos[0];		# 3' end exon
	    if ($len >= 3) {
	        $$pos[0] += 3;
	    } elsif (!$len) {
		$i += 2;
	    }
	} else {
	    my $len = $$pos[$n] - $$pos[$n - 1];	# 5' end exon
	    if ($len >= 3) {
		$$pos[$n] -= 3;
	    } elsif (!$len) {
		$n -= 2;
	    }
	}
	$line .=  "complement(" if ($sgn eq '-');
	$line .= "join(";
	for ( ; $i < $n - 1; $i += 2) {
	    my $exon = "$$pos[$i]..$$pos[$i+1],";
	    if (length($line) + length($exon) > $linewidth) {
		print OFD $line, "\n";
		$line = ";C ";
	    }
	    $line .= $exon;
	}
	my $exon = "$$pos[$i]..$$pos[$i+1])";
	if (length($line) + length($exon) > $linewidth) {
	    print OFD $line, "\n";
	    $line = ";C ";
	}
	print OFD $line . $exon;
	print OFD ")" if ($sgn eq '-');
	print OFD "\n";
}
