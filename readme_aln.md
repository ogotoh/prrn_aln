# ALN Information

## 

*   [Overview](#Ov)
*   [Install](#Inst)
*   [User Operation](#Usr)
*   [Gene-structure prediction](#GeneStructure)
*   [Examples](#Exam)
*   [References](#Ref)

* * *

## <a name="Ov">Overview</a>

**Aln** is a program for aligning a pair of nucleotide or amino acid sequences or groups of sequences. **Aln** accepts various combinations of distinct types of inputs: a single nucleotide sequence (SN), a single amino acid sequence (SA), a pre-aligned group of nucleotide sequences (MN), and a group of amino acid sequences (MA). **Aln** can also perform mixed alignment between a nucleotide sequence and a protein sequence or a group of protein sequences [[5]](#Ref5). This feature of **aln** enables us to predict eukaryotic gene structures (protein-coding exons) based on sequence homology with known protein sequence(s). The details of this procedure will be described in a later section entitled '[Gene-structure prediction](#GeneStructure)'.

In sequence alignment, the hardest problem is to develop an efficient optimization algorithm under realistic constraints on deletions and insertions (indels or gaps). Currently the most popular methods use affine functions (AF) of the form, w(_k_) = _uk_ + _v_, to penalize an indel [[1]](#Ref1). A slightly more general penalty function is a piecewise linear function [[2]](#Ref2), a special and the simplest form of which is a double affine function (DA), w(_k_) = _uk + v_ for _k_ <= _K_ and w(_k_) = w(_K_) + _u'_(_k-K_) for (_k_ > _K_), where _u_, _u'_ (< _u_), _v_, and _K_ are non-negative constants. With the default settings, **aln** automatically selects seemingly the best algorithm, depending on the types of input sequences, as follows.

(1A) SA x SA : DP-DA  
(1N) SN x SN : DP-DA  
(2A) SA x MA or MA x MA : DP-AF  
(2N) SN x MN or MN x MN : DP-AF  
(3) SN x SA or SN x MA : DP-AF

where DP implies dynamic programming algorithm. More rigorous candidate-list (CL) algorithm may be used as an option. CL is a generalized version of DP, and was introduced to rigorously optimize a group-to-group alignment score with an affine gap-penalty function [[3]](#Ref3). **Aln** converts the input MA or MN into a generalized profile [[4]](#Ref4), if the performance is expected to be improved compared to the ordinary character-based method.

## <a name="Inst">Install</a>

```
% cd src
% ./configure [--help]
% make
% make install
% make clearall
```
If you have changed the location of the table directory after installation, set the env variable ALN_TAB:
  * % setenv ALN_TAB New_Aln_Tab (csh/tsh)
  * $ export ALN_TAB=New_Aln_Tab (sh/bsh)

## <a name="Usr">User Operation</a>

There are two methods to invoke a program. 

```
[1]	aln [option2] seq1 seq2
[2]	aln option1 [option2] [seq2]
```

Option2 is common to other programs (e.g. **Prrn**), and default values will be used if omitted. Method [1] is an 'instant mode'; the calculation starts immediately using the sequence data in files seq1 and seq2\. Method [2] is a 'batch mode', which will be discussed after explanation of the instant mode.

### Notice

At the beginning of each process, several files in the 'table' directory are read. Follow the instructions under Install to properly locate these files.

Usually, **aln** automatically judges the type of a sequence (amino acids or nucleotides) according to the composition of characters read in. You may explicitly specify the sequence type by adding an 'attribute' as described below.

### [1] Instant mode

This mode is the simplest to use. Two sequence files to be compared are given in the command line arguments, and the result is reported to stdout (standard output). Command line options, option2, are common to all modes.

#### Option2

Each option is specified by a dash followed by one or two 'option' character(s). Some options take a second argument as shown below, where _N_ is a number (real or integer), and _S_ is a word (file name in most cases). The options with an asterisk are meaningful only for DNA vs. protein sequence alignment.

<pre>
option	default	range	description
-FN      1      1-9     Output format. N=1-5: native; N=6: Phylip; N=7: GCG; 
                        N=8: CLUSTA; N=9: FASTA.
-FS             c,f,h,n S=c: CLUSTAL; S=f: FASTA; S=h: Phylip; S=n: NEXUS
-KS	none	A,D,P	A: amino acid, D: nucleotide, P: second is amino acid 
                          sequence (default: infer).
-L[N]	 0	0-15	Control penalty values given to terminal gaps.
-ON			Output mode.  See the 'Gene structure prediction'.
-RN	 0	>=0	The number of random sequence pairs used for a shuffle test.
-TS	ALN_TAB		The 'table' directory.
-lN	60	>8	Set line width = N.
-mS			Amino acid exchange matrix.
-n			Report only statistics rather than an alignment.
-oS	SO		Output file.
-pi	 	        Display alignment with intron positions highlighted by color.
-ph	  	        Same as above but output an html file.
-pp                     Add taxonomic label to each sequence ID (Note 2).
-pq                     Quiet mode.
-r			Don't remove intermediate files (with -a or -b).
-sS			Set the default path where sequence files reside.
-uN	2	>=0	Set the gap-extension penalty u = N.
-vN	9	>=0	Set the gap-opening penalty v = N.
-wN	100	>=0	Set shldr = N.
-ybN	0		Add N to each element of an amino acid exchange matrix.
-yiN*	15	>0	Intron penalty.
-yjN	0	>=0	The u' value in a DA gap penalty. Must be 0 <= u' < u.
-ykN	10	>0	The K value in a DA gap penalty.
-ylN*	1	1-3	N=1: assume the absence of introns; N=2: assume the presence of
			introns, and use the AF penalty: N=3: assume the presence of
			introns, and use the DA penalty.
-ymN	2	>=0	Similarity measure between identical nucleotides.
-ynN	2	>=0	Dissimilarity measure between different nucleotides.
-yoN	30	<=0	Penalty for a premature termination codon.
-ypN	250	0-300	The 'PAM' value of the mutation data matrix.
-yuN	2	>=0	Same as -uN.
-yvN	9	>=0	Same as -vN.
-ywN	100	>=0	Same as -wN.
-yxN*	30	>=0	Penalty for a frameshift.
-yyN*	8	>=0	Relative contribution rate of translational and splicing
			boundary signals to the total alignment scores.
-yzN*	2	>=0	Relative contribution rate of coding potential to the total
			alignment scores.
-yJN	25	>=0	Bonus given to a matched intron positions.
</pre>

### [2] Batch mode

This mode enables several calculations with a single command. Given a list of sequence names in a catalog file, you can compare various combinations of sequences according to option1\. Option1's are mutually exclusive; you should choose one of them. In a special case, you get a multiple-sequence alignment by recurrently applying pairwise alignment to groups of sequences. In this case, you must prepare a special format of tree file instead of a catalog file. If :catalog in the following list is omitted, individual sequences given in the arguments are compared with one another.

#### Option1

<pre>-acatalog    Construct multiple-sequence alignment simply by adding on in the
             order of the file names in the catalog.
-btree       Construct multiple-sequence alignment by the progressive method.
-ia:catalog  (2i-1)-th sequence vs. 2i-th sequence (i > 0).
-ie:catalog  All pairs of sequences in the catalog.
-if:catalog  First sequence vs. all other sequences in the catalog.
-ig:catalog  Every sequence in group1 vs. every sequence in  group2\.  The two
             groups are separated by a blank line.
-ii:catalog  Self comparisons.  Not useful with aln.
-ij:catalog  (2i-1)-th sequence vs. 2i-th sequence (i > 0) (same as -ia).
-il:catalog  Last sequence vs. all other sequences in the catalog.
-ip File1 File2	parallel read from File1 (e.g. genomic) and File2 (e.g. protein)
</pre>

### Tree file

A tree file appropriate for an input to **aln** can be made by upg/nj included in the distribution version. You can also perform the equivalent calculations with the **prrn** -b option without generating intermediate files.

The tree topology may be specified in the 'New Hampshire Standard' format or 'Newick' format, in which the OTU names correspond to the sequence file names to be aligned. Alternatively, you may directly indicate the order of calculation as shown in the following example.

<pre>
AA     ---------
1               |
BB              |
0               |    Content of a tree file.
CC              |    AA,...EE show sequence file names.
3               |    Numbers indicate the order of calculations.
DD              |
2               |
EE     --------
</pre>

In this example, BB and CC are aligned first, and the result is stored in a file named _s0 in the current directory. _s0 is then aligned with AA to give _s1\. Alignment between DD and EE is stored in _s2, which is further aligned with _s1, to give _s3, the final result. The number between sequence names indicates the distance to the internal node from the leaf furthest from the root. Since the present implementation uses only topological information, you need not care about precise distance values. These intermediate files _s.. are deleted automatically unless you give the '-r' command line option, and the final result is written in either the file specified by -o option or './ALN' if not specified. A second run of `aln -b tree` may replace the previous files, so that the final result (ALN or _s3 with -r option in this example) should be copied or renamed before another run of **aln**.

The strategy adopted here is the so-called 'progressive method' similar to those proposed by Hogeweg and Hesper (1984) _J. Mol. Evol._ **20**, 175-186, Feng and Doolittle (1987) _J. Mol. Evol._ **25**, 351-360, and others. This strategy rapidly generates a reasonably good alignment, which may be further refined by **prrn**. Please refer to the [accompanying document](./readme_prrn.md) for further details.

### Catalog file

Each line in a catalog file has the same format as that used for input of a sequence name under the Instant mode:

<pre>[[path][filename]] [start] [end] [attribute]
</pre>

With the -g option, a blank line separates two groups, otherwise a blank line (optional) indicates the end of the list.

## <a name="GeneStructure">Gene-structure prediction</a>

**Aln** also supports spliced alignment between a nucleotide sequence and a protein sequence or a protein profile. To run the program in this mode, you must confirm that the following five files with asterisk are present in the 'table' directory. Species specific (ss) files may be present in each subdirectory.

```
AlnParam (ss)
Branch (ss)
CodePotTab*
Intron53 (optional)
IntronPotTab (ss)
Splice3* (ss)
Splice5* (ss)
TransInit*
TransTerm*
```

These files are included in the distribution package. Subdirectories in the 'table' directory contain species-specific parameters. You can select the species closest to your genome by the -Txxxxxx option. Otherwise, the default parameter set is used.

The option -yl2 or -yl3, together with G attribute, tells the program that the nucleotide sequence may contain introns which are skipped at translation. As a consequence, the exon/intron structure of a genomic sequence can be predicted [[5]](#Ref5)[[6]](#Ref6). An affine (-yl2) or a double affine (-yl3) function may be chosen to penalize a gap (indel) at the translated amino-acid sequence level. The -L option should also usually be set to perform semi-global alignment. The standard command would look like:

<pre>% aln -yl2 (or -yl3) -L -ON 'DNA N1 N2 G' reference	(sense strand)

% aln -yl2 (or -yl3) -L -ON 'DNA N1 N2 < G' reference	(antisense strand)
</pre>

where N1 and N2 specify the region in the genomic DNA sequence within which the objective gene may reside, and 'reference' indicates a protein sequence or profile used as the reference. The character 'G' indicates that the DNA sequence can contain introns. The 'G' is optional when the reference is an amino acid sequence or profile. The -ON option controls the output:
```
-O0:	Gff3 gene format
-O1:	alignment between the DNA sequence and the reference.
-O2:	Gff match format
-O3:	Bed format
-O4:	exon-oriented native format
-O5:	intron-oriented native format
-O6:	concatenated exonic sequence i.e. 'cDNA' sequence.
-O7:	translated amino acid sequence.
```

With -O6 and -O7, boundary signal strengths are also reported after the end of the sequence. In addition, translation is conducted after correction for potential frameshift errors with -O7.

The accompanying file 'Exam' contains some examples for the application of **aln** to test data with several different option settings.

### **aln** and [**spaln**](https://github.com/ogotoh/spaln.git)

The following commands are nearly equivalent:

<pre>% spaln -Q0 'DNA N1 N2' reference
% aln -yl2 -L -O4 -yp150 'DNA N1 N2 G' reference
</pre>

However, there are several differences between **spaln** and spliced alignment mode of **aln**.

1.  **Aln** requires -yl2 (or -yl3) -L options to perform spliced alignment.
2.  **Spaln** does not need of these options or 'G' attribute to the genomic sequence.
3.  **Spaln** cannot use multiple alignment as the reference while **aln** can.
4.  **Spaln** has the genome mapping phase but **aln** does not.
5.  **Spaln** has quick modes whereas the quick modes of **aln** are currently disabled.
6.  Score values are integer in **spaln** whereas floating point in **aln**.
7.  The default PAM level is 150 in **spaln** whereas 250 in **aln**.

## <a name="Exam">Examples</a>
```
% aln -s dir -KA a b		: Align amino acid sequences a and b in dir directory
% aln -ie:catalog -s dir -n	: Calculate alignment scores for all combination
                                  of sequences in dir listed catalog
% aln -yl2 -L 'genomic 1001 3000 G' 'cDNA D'
                                : spliced alignment with cDNA reference
% aln -yl2 -L -O7 -T NematodC 'C_ele_Chr_II 9798260 9800250 <' aa.msa
                                : output predicted amino acid sequence
```

## <a name="Ref">References</a>

<a name="Ref1">[[1]](http://www.sciencedirect.com/science/article/pii/0022283682903989?via%3Dihub) Gotoh, O. (1982) "An improved algorithm for matching biological sequences." _J. Mol. Biol._ **162**, 705-708.

<a name="Ref2">[[2]](https://link.springer.com/article/10.1007/BF02458577) Gotoh, O. (1990) "Optimal sequence alignment allowing for long gaps." _Bull. Math. Biol._ **52**, 359-373.

<a name="Ref3">[[3]](https://www.semanticscholar.org/paper/Optimal-alignment-between-groups-of-sequences-and-Gotoh/4a5eda88df350fae5fff64dcb790fe3acadd3825) Gotoh, O. (1993) "Optimal alignment between groups of sequences and its application to multiple sequence alignment." _CABIOS_ **9**, 361-370.

<a name="Ref4">[[4]](https://academic.oup.com/bioinformatics/article-abstract/10/4/379/230674?redirectedFrom=PDF) Gotoh, O. (1994) "Further improvement in group-to-group sequence alignment with generalized profile operations." _CABIOS_ **10**, 379-387.

<a name="Ref5">[[5]](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/16.3.190) Gotoh, O. (2000) " Homology-based gene structure prediction: simplified matching algorithm by the use of translated codon (tron) and improved accuracy by allowing for long gaps." _Bioinformatics_ **16**, 190-202.

<a name="Ref6">[[6]](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btn460?ijkey=XajuzvyHlcQZoQd&keytype=ref) Gotoh, O. (2008) "Direct mapping and alignment of protein sequences onto genomic sequence." _Bioinformatics_ **24** (21) 2438-2444.

* * *

Copyright (c) 1997-2017 Osamu Gotoh (o.gotoh@aist.go.jp) All Rights Reserved.

