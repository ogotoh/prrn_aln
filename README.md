# [ALN](./readme_aln.md), [PRRN](./readme_prrn.md), and [Refgs.pl](./readme_refgs.md)

## Pairwise sequence alignment, multiple sequence alignment, and homology-based gene prediction
### Present Version 5.1.0
### Last Updated 2017-02-15

*   [Overview](#Ov)
*   [System](#Sys)
*   [Sequence File](#Sf)
*   [List of Related Publications](#Ref)

* * *

## <a name="Ov">Overview</a>

This repository contains sequence alignment and homology-based gene prediction programs developed by Osamu Gotoh and associates. [**Aln**](./readme_aln.md) performs pairwise alignment between sequences or pre-aligned groups of sequences, whereas [**Prrn**](./readme_prrn.md) performs multiple sequence alignment (MSA) of DNA or protein sequences. These programs use common internal codes and the same sequence file format described below. [**Refgs.pl**](./readme_refgs.md) is a perl script that performs concerted gene prediction and multiple sequence alignment, taking advantage of the unique feature of **prrn** that generates gene-structure-aware multiple protein sequence alignment (GSA-MPSA). A related program [**spaln**](https://github.com/ogotoh/spaln.git) that seamlessly performs genome mapping and spliced alignment is presented in a separate repository. For further details, crick on the individual tags.

## <a name="Sys">System</a>

Except for **refgs.pl**, the programs are written in C++ (ISO/C++), and distributed as source code. Users must compile the programs on their own system. Although the programs have been tested only on Linux and Cygwin, they are likely to be portable to most UNIX systems with little or no modifications. The distribution package also contains source codes of several related programs used by **refgs.pl**.

## <a name="Sf">Sequence File</a>

A sequence file may contain either a single sequence or a set of pre-aligned sequences. Avoid file names starting with a numeral.

### Single sequence

Almost all popular formats for nucleotide and amino acid sequences, FASTA, GenBank, EMBL, SwissProt, PIR, and ProDB, are acceptable. A simple text file without any additional information is also fine. In any format, a line beginning with a semicolon (;) is regarded as a comment line. The specific format is recognized automatically by the first word in the first non-comment line of each file. Only single-letter codes are accepted for an amino acid. The nomenclature recommended by NC-IUB (Eur. J. Biochem. (1985) 150, 1-5) is followed to represent ambiguous nucleotides. The programs are not case-sensitive.

To represent the exon-intron structure of the parental gene, the format of FASTA file should be extended. A line starting with ';C' shows the exon boundaries (inclusive). More than one line may be used if necessary. The format after ';C ' is essentially the same as that of Feature field of a GenBank file. Start and end positions of each exon are separated by two dots. Individual exons are delimited by a comma. The term 'complement' indicates that the corresponding gene lies in the complementary strand of the genomic sequence. Two examples are as follows:

```
>ce13a1 C. elegans chromosome II positive strand  
;C join(9803525..9803710,9803766..9804097,9804152..9804251,  
;C 9804299..9804855,9804926..9805069,9805115..9805349)  
MSFSILIAIAIFVGIISYYLWIWSFWIRKGVKGPRGLPFLGVIHKFTNYENPGALKFSEW  
TKKYGPVYGITEGVEKTLVISDPEFVHEVFVKQFDNFYGRKLTAIQGDPNKNKRVPLVAA  
QGHRWKRLRTLASPTFSNKSLRKIMGTVEESVTELVRSLEKASAEGKTLDMLEYYQEFTM  
DIIGKMAMGQEKSLMFRNPMLDKVKTIFKEGRNNVFMISGIFPFVGIALRNIFAKFPSLQ  
MATDIQSILEKALNKRLEQREADEKAGIEPSGEPQDFIDLFLDARSTVDFFEGEAEQDFA  
KSEVLKVDKHLTFDEIIGQLFVFLLAGYDTTALSLSYSSYLLATHPEIQKKLQEEVDREC  
PDPEVTFDQLSKLKYLECVVKEALRLYPLASLVHNRKCLKTTNVLGMEIEAGTNINVDTW  
SLHHDPKVWGDDVNEFKPERWESGDELFFAKGGYLPFGMGPRICIGMRLAMMEMKMLLTN  
ILKNYTFETTPETVIPLKLVGTATIAPSSVLLKLKSRF  
[EOF]

>ce13a2 C. elegans chromosome II negative strand  
;C complement(join(9798263..9798503,9798584..9798727,9798905..9799461,  
;C 9799519..9799618,9799680..9800011,9800058..9800243))  
MSLSILIAGASFIGLLTYYIWIWSFWIRKGVKGPRGFPFFGVIHEFQDYENPGLLKLGEW  
TKEYGPIYGITEGVEKTLIVSNPEFVHEVFVKQFDNFYGRKTNPIQGDPNKNKRAHLVSA  
QGHRWKRLRTLSSPTFSNKNLRKIMSTVEETVVELMRHLDDASAKGKAVDLLDYYQEFTL  
DIIGRIAMGQTESLMFRNPMLPKVKGIFKDGRKLPFLVSGIFPIAGTMFREFFMRFPSIQ  
PAFDIMSTVEKALNKRLEQRAADEKAGIEPSGEPQDFIDLFLDARANVDFFEEESALGFA  
KTEIAKVDKQLTFDEIIGQLFVFLLAGYDTTALSLSYSSYLLARHPEIQKKLQEEVDREC  
PNPEVTFDQISKLKYMECVVKEALRMYPLASIVHNRKCMKETNVLGVQIEKGTNVQVDTW  
TLHYDPKVWGEDANEFRPERWESGDELFYAKGGYLPFGMGPRICIGMRLAMMEKKMLLTH  
ILKKYTFETSTQTEIPLKLVGSATTAPRSVMLKLTPRHSN  
[EOF]
```

The line starting with ';C' may be simplified as follows:
```
>ce13a1 C. elegans chromosome II positive strand  
;C + 9803525 9803710 9803766 9804097 9804152 9804251  
;C 9804299 9804855 9804926 9805069 9805115 9805349  
...

>ce13a2 C. elegans chromosome II negative strand  
;C - 9798263 9798503 9798584 9798727 9798905 9799461  
;C 9799519 9799618 9799680 9800011 9800058 9800243  
...
```

When the parental gene contains one or more frame shifts, an output from **aln** or **spaln** contains the corresponding number of lines starting with ';M' such that:
  * ;M Deleted _n_ chars at _p_
  * ;M Insert _n_ chars at _p_

The first line indicates that _n_ (_n_ = 1 or 2) nucleotides have been deleted from the parental genomic sequence beginning at site _p_. Likewise, the second line indicates that _n_ blank characters are inserted after site _p_ to maintain the open reading frame. Such kinds of information are used to properly juxtapose intron positions along the alignment.

### Sequential format of multiple sequences

The concatenation of multiple sequences in one of the above single-sequence formats, except for the naked sequence-alone format, may be used as a program input. Extended FASTA format sequences as described above should be used to indicate the exon-intron organizations of the parental genes. You can indicate the number of sequences and the length of the longest sequence in the first line of the file by two numbers separated by one or more spaces. However, this is no longer mandatory and a file in the ordinary multi-fasta format is fine. Deletion characters are automatically padded at the end of shorter sequences. If the sequences are pre-aligned, internal deletion characters are preserved. Two examples of sequential formats are shown below. A multi-fasta file with gene structure information may be found in the sample directory of the distribution package.

<table>

<tbody>

<tr>

<td align="left">File 1:</td>

<td align="center" width="40">|</td>

<td align="left">File 2:</td>

</tr>

<tr>

<td></td>

<td align="center" width="40">|</td>

<td></td>

</tr>

<tr>

<td align="left">>Seq1</td>

<td align="center" width="40">|</td>

<td align="left">2 16</td>

</tr>

<tr>

<td align="left">aaatttcccggg</td>

<td align="center" width="40">|</td>

<td align="left">LOCUS Seq1</td>

</tr>

<tr>

<td align="left">>Seq2</td>

<td align="center" width="40">|</td>

<td align="left">ORIGIN</td>

</tr>

<tr>

<td align="left">atcgatcgatcgat</td>

<td align="center" width="40">|</td>

<td align="left">1 AAATTTCCCGGG</td>

</tr>

<tr>

<td align="left">>Seq3</td>

<td align="center" width="40">|</td>

<td align="left">//</td>

</tr>

<tr>

<td align="left">tcgatctcgagaat</td>

<td align="center" width="40">|</td>

<td align="left">LOCUS NCODE</td>

</tr>

<tr>

<td align="left">>Seq4</td>

<td align="center" width="40">|</td>

<td align="left">ORIGIN</td>

</tr>

<tr>

<td align="left">aaatgacacggg</td>

<td align="center" width="40">|</td>

<td align="left">1 -ACMGRSVTWYHKDBN</td>

</tr>

<tr>

<td align="left">ggatcaggctagg</td>

<td align="center" width="40">|</td>

<td align="left">//</td>

</tr>

</tbody>

</table>

### Interleaved (native) format of multiple sequences

This native format is designed for multiple-sequence alignment to be naturally recognized by human eyes. The alignment produced by **aln** can be used as an input to **aln** or **prrn**, and this is the most common way to have access to sets of pre-aligned sequences. Thus, the format of an aligned sequence file is the same as the default output format of **aln**. The first non-blank line in a file must indicate the number of sequences, _N_, involved in the alignment. This number is obtained as the sum of numbers in square brackets, e.g., when the first line is

Seq1[3] - Seq2[4]

_N_ is calculated to be 7\. Subsequent lines up to the first blank line are ignored. The rest of the file is composed of one or more blocks of a fixed column width of less than 254 characters. Each block is composed of _N_ 'sequence lines' and other (optional) 'non-sequence lines'. The general format of a sequence line is:

\<Position\> \<Sequence|\> \<Name\>

where \<Position\> is a numeral that indicates the sequence position of the first letter in \<Sequence|\> (Usually all \<Position\>s in the first block are 1, but it is not a prerequisite. Negative values are also appropriate). A line lacking the \<Position\> field is regarded as a 'non-sequence line' and ignored upon reading. The _i_-th sequence line in the second block is concatenated to the _i_-th sequence line in the first block, and so on. There is no particular limit on _N_ or the length, but the total number of characters to be stored is limited by MAXAREA defined in src/sqio.h. Several examples of native format are provided in the sample directory.

### Special Characters

Three characters, dash '-', tilde '\~' and caret '^' have special meanings in a multiple-sequence file of the native format. A dash indicates a 'deletion' introduced by some alignment procedure. Be careful not to use a space or dot instead of a dash. Spaces and dots are simply ignored, so that the file may be interpreted in a totally unexpected way. A tilde means the same residue as that in the first sequence line on that column in the block. On the other hand, a caret means the same residue as that in the previous sequence line on that column. These ditto characters are convenient in representing an alignment of closely related sequences. Neither '\~' nor '^' is allowed in the first sequence line in each block.

### Information on gene structures

For native format of multiple sequences, lines starting with ";B ", ";b ", and ";m " represent the information about organizations of corresponding genes. The first number that follows ";B ", _NP_, indicates the number of alignment positions where an intron intervenes at least one of the genes corresponding to the aligned protein or cDNA sequences. For a protein sequence, the position means that of coding nucleotide sequence so that phase as well as alignment column is also significant. The second number in this line, _NI_, indicates the total number of introns. The numbers that follow ";b " indicate the number of sequences that contain the intron at the 1st, 2nd, ..., _NP_-th position. The numbers that follow ";m " present the list of sequences that contain the intron at the 1st, 2nd, ..., _NP_-th position in this order. See the example ce13a.mfa in sample/pas directory.

## <a name="Ref">List of Related Publications</a>

[[1]](http://www.sciencedirect.com/science/article/pii/0022283682903989?via%3Dihub) Gotoh, O. (1982) "An improved algorithm for matching biological sequences." _J. Mol. Biol._ **162**, 705-708.

[[2]](https://ac.els-cdn.com/S0022519386801126/1-s2.0-S0022519386801126-main.pdf?_tid=7676197a-c463-11e7-b427-00000aacb361&acdnat=1510131893_4413edc236323e8ba199e83d25d90423) Gotoh, O. (1986) Alignment of three biological sequences with an efficient traceback procedure. _J. Theor. Biol._ **121**, (3) 327-337.

[[3]](https://link.springer.com/article/10.1007/BF02458577) Gotoh, O. (1990) "Optimal sequence alignment allowing for long gaps." _Bull. Math. Biol._ **52**, 359-373.

[[4]](http://www.sciencedirect.com/science/article/pii/S0092824005803593) Gotoh, O. (1990) Consistency of optimal sequence alignments. _Bull. Math. Biol_. **52**, (4) 509-525.

[[5]](https://www.semanticscholar.org/paper/Optimal-alignment-between-groups-of-sequences-and-Gotoh/4a5eda88df350fae5fff64dcb790fe3acadd3825) Gotoh, O. (1993) "Optimal alignment between groups of sequences and its application to multiple sequence alignment." _CABIOS_ **9**, 361-370.

[[6]](https://www.jstage.jst.go.jp/article/gi1990/4/0/4_0_109/_pdf) Gotoh, O. (1993) "Extraction of conserved or variable regions from a multiple sequence alignment." _Proceedings of Genome Informatics Workshop_ **IV**, pp. 109-113\.

[[7]](https://academic.oup.com/bioinformatics/article-abstract/10/4/379/230674?redirectedFrom=PDF) Gotoh, O. (1994) "Further improvement in group-to-group sequence alignment with generalized profile operations." _CABIOS_ **10**, 379-387\.

[[8]](https://academic.oup.com/bioinformatics/article-abstract/11/5/543/236157) Gotoh, O. (1995) "A weighting system and algorithm for aligning many phylogenetically related sequences." _CABIOS_, **11**, 543-551.

[[9]](http://www.sciencedirect.com/science/article/pii/S0022283696906798?via%3Dihub) Gotoh, O. (1996) "Significant improvement in accuracy of multiple protein sequence alignments by iterative refinement as assessed by reference to structural alignments." _J. Mol. Biol._ **264**, 823-838.

[[10]](https://academic.oup.com/mbe/article/15/11/1447/973602/Divergent-structures-of-Caenorhabditis-elegans) Gotoh, O. (1998) "Divergent structures of Caenorhabditis elegans cytochrome P450 genes suggest the frequent loss and gain of introns during the evolution of nematodes." _Mol. Biol. Evol._ **15**, 1447-1459.

[[11]](http://www.sciencedirect.com/science/article/pii/S0022283696906798?via%3Dihub) Gotoh, O. (1999) "Multiple sequence alignment: algorithms and applications." _Adv. Biophys._ **36**, 159-206.

[[12]](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/16.3.190) Gotoh, O. (2000) " Homology-based gene structure prediction: simplified matching algorithm by the use of translated codon (tron) and improved accuracy by allowing for long gaps." _Bioinformatics_ **16**, 190-202.

[[13]](http://www.crcnetbase.com/doi/abs/10.1201/9781420036275.ch3) Gotoh, O., Yamada, S., and Yada, T. (2006) Multiple Sequence Alignment, in _Handbook of Computational Molecular Biology,_ (Aluru, S. ed.) Chapman & Hall/CRC, Computer and Information Science Series, Vol. 9, pp. 3.1-3.36.

[[14]](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-524) Yamada, S., Gotoh, O., Yamana, H. (2006) Improvement in accuracy of multiple sequence alignment using novel group-to-group sequence alignment algorithm with piecewise linear gap cost, _BMC Bioinformatics_, **7**, 524.

[[15]](http://nar.oxfordjournals.org/cgi/content/abstract/gkn105?ijkey=N2yLVza41RuShAg&keytype=ref) Gotoh, O. (2008) "A space-efficient and accurate method for mapping and aligning cDNA sequences onto genomic sequence" _Nucleic Acids Research_ **36** (8) 2630-2638.

[[16]](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btn460?ijkey=XajuzvyHlcQZoQd&keytype=ref) Gotoh, O. (2008) "Direct mapping and alignment of protein sequences onto genomic sequence." _Bioinformatics_ **24** (21) 2438-2444.

[[17]](https://www.jstage.jst.go.jp/article/ipsjtbio/1/0/1_0_2/_article) 110.	Yamada, S., Gotoh, O., Yamana, H. (2008) Improvement in speed and accuracy of multiple sequence alignment program PRIME, _IPSJ Trans. Bioinformat._, **1**, 2-12.

[[18]](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-224) Nakato, R. and Gotoh, O. (2010) Cgaln: fast and space-efficient whole-genome alignment, _BMC Bioinformatics_, **11**, 224.

[[19]](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-45) Iwata, H. and Gotoh, O. (2011) Comparative analysis of information contents relevant to recognition of introns in many species, _BMC Genomics_, **12**, 45.

[[20]](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gks708) Iwata, H. and Gotoh, O. (2012) Benchmarking spliced alignment programs including  Spaln2, an extended version of Spaln that incorporates additional species-specific features, _Nucleic Acids Res._, **40** (20) e161.

[[21]](http://link.springer.com/protocol/10.1007%2F978-1-62703-646-7_2) Gotoh, O. (2013) "Heuristic Alignment Methods, in _Multiple Sequence Alignment Methods, (Russel, D. ed.), Methods in Molecular Biology_, Vol. **1079**, pp. 29-44, Humana Press.

[[22]](http://www.biomedcentral.com/1471-2105/15/189) Gotoh, O., Morita, M. Nelson, D.R. (2014) "Assessment and refinement of eukaryotic gene structure prediction with gene-structure-aware multiple protein sequence alignment", _BMC Bioinformatics_ **15**:189.

* * *

Copyright (c) 1997-2017 Osamu Gotoh (o.gotoh@aist.go.jp) All Rights Reserved.

