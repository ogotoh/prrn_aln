# PRRN information

## 

*   [Overview](#Ov)
*   [Install](#Inst)
*   [User Operation](#Usr)
*   [Examples](#Exam)
*   [References](#Ref)

* * *

## <a name="Ov">Overview</a>

The program **prrn** is an implementation of the randomized iterative strategy for multiple sequence alignment. **Prrn** accepts either nucleotide or protein sequences. **Prrn** repeatedly uses pairwise group-to-group alignment [[1]](#Ref1)[[2]](#Ref2) to improve the overall weighted sum-of-pairs score at each iterative step, where the pair weights are introduced to correct for uneven representations of the sequences to be aligned [[3]](#Ref3). The doubly-nested randomized iterative (DNR) method described in Ref. [[4]](#Ref4)[[5]](#Ref5) tries to make alignment, phylogenetic tree and pair weights mutually consistent. As this is a hill-climbing strategy, the true optimum may not be attained. Although earlier versions of **prrn** required that the input sequences must be pre-aligned by other methods, e.g. progressive alignment method, **prrn** now internally generate a provisional alignment from unaligned sequences. From version 4, gene-structure-aware multiple protein sequence alignment (GSA-MPSA) is obtained if each input sequence possesses information about the exon-intron organization of the parental gene [[6]](#Ref6). Version 5 implements a faster algorithm for constructing the guide tree from a set of unaligned sequences. The internal data structure that holds the information of parental gene structures has been modified from that of Version 4. In addition, Version 5 has enabled multi-thread computation. These improvements have made Version 5 considerably faster and more accurate than the previous versions. 

## <a name="Inst">Install</a>

```
% cd src
% ./configure [--help]

% make prrn
  or
% make prrnall

% make install
% make clearall
```

`% make prrn` generates the latest version (Ver5.1.0), whereas
`% make prrnall` generates both Ver4.1.0 and Ver5.1.0

If you have changed the location of the table directory after installation, set the env variable ALN_TAB:
  * % setenv ALN_TAB New_Aln_Tab (csh/tsh)
  * $ export ALN_TAB=New_Aln_Tab (sh/bsh)

## <a name="Usr">User Operation</a>

**prrn** is invoked from command line as follows.

```
	prrn [options1] [options2] seq1 [seq2 seq3 ...]
```

Options1 is specific to **prrn**, whereas options2 is common to other programs, and default values will be used if it is omitted. The sequence data in the command line arguments, seq_n_, may or may not be pre-aligned. When two or more sequences or an unaligned multiple sequence file are given, they are combined and pre-aligned by the progressive method according to an internally calculated guide tree. The `-b Tree` option described below generates the initial alignment by the progressive method using the given 'Tree' as the guide tree. The result is reported to SO (stdout) or the file specified by option -o (see below).

### Notice

At the beginning of each process, **prrn** reads the score files. Follow the instructions above to tell the program where these files are located.

### Options

Each option is specified by a dash followed by one or two character(s). Some options take an argument as shown below. (N is a number, S is a word, A space between the option character(s) and its argument is optional).

<pre>Options1  Default         Meaning

-DN      35     Threshold height that divides subtrees.
-ES      SE     Destination of supplementary messages.
-GN       2     Grouping.
-IN[:M]  10     N: Maximum number of iterations in the outer loop.
	        M: Maximum number of iterations in an inner loop.
-JN       1     Mode of division. N=1: 1:N-1; N=2 branch of tree; N=3: all.
-ON       1     Set output mode (omode) to N.
		  1: output alignment
		  2: output outlier indel information
		  4: output normalized alignment scores.
-Q        0     Use faster but less rigorous algorithm B (ref [[1]](#Ref1))
-RN	        Seed of the series of random numbers.
-SN       1     Number of iteration series (recalculation of tree and weights).
-U              Update mode. Members in seq1 are replaced by sequences of the
                same names in seq_n_ (n>1).  Seq1 must be a multiple alignment
                whose members should have unique names.
-XN     inf     Max member size of each subtree (subalignment).
-YHN     20     Threshold value used in finding conserved regions.
-bS    tree     Initial alignment is calculated by a progressive method using S
                as the guide tree.
-cN     0.1	Alpha value for outlier detection.
-eS             Prefix to subalignments.
-rN       1     Grain size for iterative refinement. N=1: serial.
</pre>

<pre>Options2  Default       Meaning

-FN       1     Output format. N=1-5: native; N=6: Phylip; N=7: GCG; 
                  N=8: CLUSTA; N=9: FASTA.
-FS             S=c: CLUSTAL; S=f: FASTA; S=h: Phylip; S=n: NEXUS
-lN      60     Set lpw (# of residues per line) = N > 8.
-KS             S=A or A=P: amino acid; S=D or S=N: nucleotide sequence
-LN       0     N=0: global; N=embpty or N=15: semiglobal
-h or -?        Show help.
-mS             Amino acid exchange matrix.
-oS      SO     Output resultant alignment to file.
-pi	  0	Display alignment with intron positions highlighted by color.
-ph	  	Same as above but output an html file.
-ps             The order of sequences in the output file is rearranged 
                  according to the calculated phylogenetic relationship.
-pp             Add taxonomic label to each sequence ID (Note 2).
-pq             Quiet mode.
-sS      ./     Set the default path to sequence files to path.
-tN	  1     Number of threads.
-uN       2     Set gap-extension penalty u = N.
-vN     4/9     Set gap-opening penalty v = N.
-wN     100     Set shldr = N.
-w-N            Set shldr = N % of shorter sequence
-yJN	 25	>=0	Bonus given to a matched intron positions.
-yeN	  0	Background gap extension penalty.
-yhN	  0	Include similarity in hydrophobicity profile into score.
-ylN	  2	N = 2: affine, N = 3: double affine gap penalty.
-ymN	  2	Set nucleotide match score.
-ynN	 -2	Set nucleotide mismatch score.
-ypN    250     Set pam = N.
-ysN	  0	Include similarity in predicted 2D propensities into score.
-ytN      1.    Factor of reducing terminal gap penalty. 1.: global; 
-yJN	 25	>=0	Bonus given to a matched intron positions.
</pre>

(Note1) SE and SO mean standard error (stderr) and standard output (stdout), respectively.  
(Note2) To make this option effective, 'gnm2tab' file must be present in the 'table' directory. In addition, the first eight characters of each sequence identifier must match one of the 'genspcID' of the first column in gnm2tab in case-insensitive manner.

## <a name="Exam">Examples</a>
```
% prrn -pi multi_fasta.faa
	: This produces MSA or GSA-MPSA of sequences in the input file.
	  Intron positions are indicated by color on the screen for GSA-MPSA.
% prrn -s dir -KA -m vtml200 a b c
	: Align amino acid sequences a, b, and c in 'dir' directory. Use 'vtml200'
          as the amino acid substitution matrix.
% prrn -U -ps -s dir -o Updated.msa ./Current.msa memb3 memb6
	: Replace sequences memb3 and memb6 in Current.msa with those in 'dir' 
	  directory, and then refine the alignment. The output is sorted on 
	  mutual similarity of members and stored in 'Updated.msa'. The file name
	  (memb3 or memb6) must be identical to the sequence identifier that
	  follows '>' in the FASTA file.
```

## <a name="Ref">References</a>

<a name="Ref1">[[1]](https://www.semanticscholar.org/paper/Optimal-alignment-between-groups-of-sequences-and-Gotoh/4a5eda88df350fae5fff64dcb790fe3acadd3825) Gotoh, O. (1993) "Optimal alignment between groups of sequences and its application to multiple sequence alignment." _CABIOS_ **9**, 361-370.

<a name="Ref2">[[2]](https://academic.oup.com/bioinformatics/article-abstract/10/4/379/230674?redirectedFrom=PDF) Gotoh, O. (1994) "Further improvement in group-to-group sequence alignment with generalized profile operations." _CABIOS_ **10**, 379-387\.

<a name="Ref3">[[3]](https://academic.oup.com/bioinformatics/article-abstract/11/5/543/236157) Gotoh, O. (1995) "A weighting system and algorithm for aligning many phylogenetically related sequences." _CABIOS_, **11**, 543-551.

<a name="Ref4">[[4]](http://www.sciencedirect.com/science/article/pii/S0022283696906798?via%3Dihub) Gotoh, O. (1996) "Significant improvement in accuracy of multiple protein sequence alignments by iterative refinement as assessed by reference to structural alignments." _J. Mol. Biol._ **264**, 823-838.

<a name="Ref5">[[5]](http://www.sciencedirect.com/science/article/pii/S0022283696906798?via%3Dihub) Gotoh, O. (1999) "Multiple sequence alignment: algorithms and applications." _Adv. Biophys._ **36**, 159-206.

<a name="Ref6">[[6]](http://www.biomedcentral.com/1471-2105/15/189) Gotoh, O., Morita, M. Nelson, D.R. (2014) "Assessment and refinement of eukaryotic gene structure prediction with gene-structure-aware multiple protein sequence alignment", _BMC Bioinformatics_ **15**:189.


* * *

Copyright (c) 1997-2017 Osamu Gotoh (o.gotoh@aist.go.jp) All Rights Reserved.

