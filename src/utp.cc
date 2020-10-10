/*****************************************************************************
*
*	Utilities for nucleotide/amino-acid sequence management
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-41-6 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#include "aln.h"
#include "mseq.h"
#include "autocomp.h"
#include "prs.h"
#include "pattern.h"
#include "utilseq.h"

inline	double	protmwt(double w, int l) {return (w - 18.015 * (l - 1));}

class Mapsite {
	int	memno;
	const	char*	memid;
	int*	memlst;
	int	lst_no;
public:
	Mapsite(int& argc, const char**& argv);
	~Mapsite() {delete[] memlst;}
	void	mapsite(Seq* sd, bool s2n);
	void	output(Seq* sd, bool s2n);
};

static	char	pattern[MAXL];
static	char	spat[] = "Pattern : ";
static	char	pros[] = "PS#/all : ";

static	void	composit(FILE* fd, Seq* sd);
static	void	proutseq(FILE* fd, Seq* sd);
static	void	list(FILE* fd, Seq* sd);
static	int	icmpf(int* a, int* b);
static	void	mute(Seq* seqs[]);
static	void	muteseq(FILE* fd, Seq* seqs[]);
static	void	printskl(FILE* fd, Seq* sd);
static	char*	pgetstr(char* pat, char* prmpt);
static	void	find(FILE* fd, Seq* sd);
static	void	prosite(Seq* sd);
static	void	nameseq(FILE* fd, Seq* sd);
static	void	invtrnsl(FILE* fd, Seq* sd);
static	Seq*	forgeseq(FILE* fd, Seq* sd);
static	void	shufseq(Seq* sd);

static const	FTYPE	aamwt[24] = {
	  0.00,  89.09, 174.20, 132.12, 133.10, 121.15, 146.15, 147.13,
	 75.07, 155.16, 131.17, 131.17, 146.19, 149.21, 165.19, 115.13,
	105.09, 119.12, 204.23, 181.19, 117.15, 132.61, 146.64, 128.16
};
static	const	double	DefThr = 50.;
static	double	NGP = 0.5;
static	Mapsite*	mapsite = 0;

Mapsite::Mapsite(int& argc, const char**& argv)
	: memno(-2), memid(0), memlst(0)
{
	if (isdigit(argv[0][2])) {
	    memno = atoi(argv[0] + 2) - 1;
	} else if (argv[0][2]) {
	    memid = argv[0] + 2;
	}
	for (lst_no = 1; lst_no < argc; ++lst_no)
	    if (!isdigit(argv[lst_no][0])) break;
	if (--lst_no) {
	    memlst = new int[lst_no + 1];
	    for (int i = 0; i < lst_no; ++i) {
		memlst[i] = atoi(*++argv);
		--argc;
	    }
	    memlst[lst_no] = -1;
	}
}

void usage()
{
	fputs("Specific Options:\n", stderr);
	fputs("\t-c:\tComposition\n", stderr);
	fputs("\t-f[pattern]:\t Find pattern\n", stderr);
	fputs("\t-l:\tPrint sequence\n", stderr);
	fputs("\t-n:\tName and range\n", stderr);
	fputs("\t-BN [-s dir] ./msa: Exon junctions\n", stderr);
	fputs("\t\tN=0:  presense/absence bits\n", stderr);
	fputs("\t\tN=1:  position in nt\n", stderr);
	fputs("\t\tN=2:  phase\n", stderr);
	fputs("\t\tN=3:  position in aa + phase\n", stderr);
	fputs("\t\tN=4:  total and lonesome introns + bits\n", stderr);
	fputs("\t\tN=5:  number of mates\n", stderr);
	fputs("\t\tN=6:  total and lonesome introns + phase\n", stderr);
	fputs("\t\tN=7:  total and lonesome introns + position in aa + phase\n", stderr);
	fputs("\t\tN=8:  intron length with -s dir\n", stderr);
	fputs("\t\tN=9:  intron length:phase with -s dir\n", stderr);
	fputs("\t\tN=10: position on CDS with -s dir\n", stderr);
	fputs("\t\tN=11: exon legnth (difference with -s dir)\n", stderr);
	fputs("\t\tN=12: CDS position\n", stderr);
	fputs("\t\tN=13: aa position + phase\n", stderr);
	fputs("\t\tN=14: lower left distance matrix\n", stderr);
	fputs("\t\tN=15: pairwise distances\n", stderr);
	exit (0);
}

template <class seq_t>
int AlnServer<seq_t>::localoption(int& argc, const char**& argv)
{
static	const	char*   Localcmds = "cfglmnrtxBM";
	const	char*	opt = argv[0] + 1;
	const	char*	val = argv[0] + 2;
	int	c = *opt;
	int	n = 1;

	if (!c) return (0);
	const	char*	al = strchr(Localcmds, c);
	if (!al) {
	  switch (c) {
	    case 'G': if (*val)	NGP = atof(val); break;
	    case 'E': algmode.mlt = atoi(val); break;
	    case 'X': setdefmolc('X'); break;
	    default: n = 0; break;
	  }
	} else if (!this->jobcode) {	// prime option
	    this->jobcode = c;
	    c = argv[1]? argv[1][0]: 0;
	    switch (this->jobcode) {
	      case 'f': case 'r':
		if (*val)	strcpy(pattern, val);
		else if (c && c != OPTCHAR) {
		    ++argv; --argc;
		    strcpy(pattern, *argv);
		} else if (!pgetstr(pattern, this->jobcode == 'f'? spat: pros)) {
		    exit (1);
		}
		break;
	      case 'B': 
		if (isdigit(*val))	algmode.nsa = atoi(val);
		break;
	      case 'm': case 'M':
		mapsite = new Mapsite(argc, argv);
		break;
	      default: break;
	    }
	    if (const char* col = strchr(val, ':')) this->set_catalog(col + 1);
	    if (this->jobcode == 'l') n = 0;
	}
	return (n);
}

template <class seq_t>
void AlnServer<seq_t>::setparam(int lvl)
{
	if (lvl & 1) {
	    setlpw(QUERY);
	    setform(QUERY);
	    setprmode(QUERY, QUERY, SILENT);
	}
	if (lvl & 2) {
	    setCodonUsage(QUERY);
	}
}

static void proutseq(FILE* fd, Seq* sd)
{
	if (OutPrm.lpw > 8) {
	    sd->typeseq(fd);
	} else {
	    algmode.nsa = OutPrm.lpw;
	    sd->fpmem_len(fd);
	}
}

static void list(FILE* fd, Seq* sd)
{
	RANGE	temp;
	int 	from = sd->left + 1;

	if (!fd) fd = qout("");
	if (!fd) return;
	sd->saverange(&temp);
	promptin("from (%d) - to (%d): ", &from, &sd->right);
	sd->left = from - 1;

	proutseq(fd, sd);
	qclose();
	sd->restrange(&temp);
}

static int icmpf(int* a, int* b)
{
	return (*a - *b);
}

static void mute(Seq* seqs[])
{
	int	positions[MAXLOC];

	int	nn = fgetiarray(positions, MAXLOC, stdin);
	if (nn == 0 || positions[0] < 0) return;
	positions[nn] = -1;
	qsort((UPTR) positions, (INT) nn, sizeof(int), (CMPF) icmpf);
	seqs[0]->copyseq(seqs[1], CPY_ALL);
	CHAR*	s = seqs[1]->at(0);
	for (int n = 0; n < nn; n += 2) positions[n]--;
	for (int n = 0, i = 0, parity = 0; i < seqs[1]->len; ++i) {
	    while (i == positions[n]) {
		++n; parity = 1 - parity;
	    }
	    if (parity) {
		for (int j = 0; j < seqs[1]->many; ++j) 
		    *s++ = seqs[1]->code->amb_code;
	    } else
		s += seqs[1]->many;
	}
}

static void muteseq(FILE* fd, Seq* seqs[])
{
	mute(seqs);
	proutseq(stdout, seqs[1]);
}

void Mapsite::output(Seq* sd, bool s2n)
{
	if (s2n)
	    sd->pos2num(memno, memlst);
	else
	    sd->num2pos(memno, memlst);
	printf("%-15s\t", (*sd->sname)[memno]);
	for (int i = 0; i < lst_no; ++i) printf("%d ", memlst[i]);
	putchar('\n');
}

void Mapsite::mapsite(Seq* sd, bool s2n)
{
	if (memno < 0 && memid) memno = sd->sname2memno(memid);
	if (memno >= -1 && memlst) {
	    qsort((UPTR) memlst, (INT) lst_no, sizeof(int), (CMPF) icmpf);
	    if (memno >= 0) output(sd, s2n);
	    else {
		int*	resv = new int[lst_no + 1];
		vcopy(resv, memlst, lst_no + 1);
		for (memno = 0; memno < sd->many; ++memno) {
		    if (memno) vcopy(memlst, resv, lst_no);
		    output(sd, s2n);
		}
		delete[] resv;
		memno = -1;
	    }
	    return;
	}

	int	positions[MAXLOC];
	memlst = positions + 1;
	while (true) {
	    lst_no = fgetiarray(positions, MAXLOC, stdin);
	    if (lst_no == 0 || *memlst == 0) break;
	    memno = positions[0] - 1;
	    --lst_no;
	    qsort((UPTR) memlst, (INT) lst_no, sizeof(int), (CMPF) icmpf);
	    output(sd, s2n);
	}
	lst_no = 0; memlst = 0;
}

static void composit(FILE* fd, Seq* sd)
{
	if (!fd) fd = qout("");
	if (!fd) return;
	sd->composition();
	sd->fphseq(fd);
	FTYPE	w = 0.;
	int	l = sd->right - sd->left;
	for (int i = 1; i < AAS; i++) w += aamwt[i] * sd->cmps[i];
	fprintf(fd, "  %4d aa\'s, Mw = %9.2f\n", l, protmwt(w, l));
	int 	j = 0;
	for (int i = 0; i < sd->code->max_code; ++i) {
	   if (sd->cmps[i] > 0) {
		fprintf(fd, "   %c = %2.0lf", 
		    sd->code->decode[i], (double) sd->cmps[i]);
		if (!(++j % 8)) putc('\n', fd);
	    }
	}
	if (j % 8) putc('\n', fd);
	qclose();
}

static void printskl(FILE* fd, Seq* sd)
{
	if (!fd) fd = qout("");
	if (!fd) return;
	CHAR*	s1 = sd->at(sd->left);
	CHAR*	s2 = s1 + 1;
	CHAR*	ts = sd->at(sd->right);
	int	i = 0, j = 0;
	while (s1 < ts) {
	    while (IsGap(*s1)) {s1 += sd->many; ++i;}
	    while (IsGap(*s2)) {s2 += sd->many; ++j;}
	    fprintf(fd, "%d\t%d\n", i++, j++);
	    s1 += sd->many;
	    s2 += sd->many;
	}
	qclose();
}

static char* pgetstr(char* pat, char* prmpt)
{
	*pat = '\0';
	progets(pat, prmpt);
	return (*pat? pat: 0);
}

static void find(FILE* fd, Seq* sd)
{
	if (!fd) fd = qout("");
	if (!fd) return;
	while (pgetstr(pattern, spat)) findpat(fd, sd, pattern);
	qclose();
}

static void prosite(Seq* sd)
{
	FILE*	fd = qout("");
	ProSite	prs(PAT, 0);

	if (!fd) return;
	while (pgetstr(pattern, pros)) {
	    if (!strncmp(pattern, "all", 3))
		prs.allprs(fd, sd);
	    else if (isdigit(*pattern))
		prs.prosites(fd, sd, parspat(pattern));
	}
	qclose();
}

static void nameseq(FILE* fd, Seq* sd)
{
	RANGE	gate;

	fprintf(fd, "%-16s [ %3d ] %5d %5d", sd->sqname(), 
	    sd->many, sd->SiteNo(sd->left), sd->SiteNo(sd->right - 1));
	if (algmode.rng && sd->findGate(&gate))
	    fprintf(fd, " %5d %5d", gate.left, gate.right);
	fputc('\n', fd);
}

static Seq* forgeseq(FILE* fd, Seq* sd)
{
	promptin("Length [%d] : ", &sd->len);
	double*	pcmp = get_mdmcmp();
	sd->randseq(pcmp);
	sd->typeseq(fd);
	return (sd);
}

static void shufseq(Seq* sd)
{
	FILE*	fd = qout("");

	if (!fd) return;
	sd->rndseq();
	sd->typeseq(fd);
	qclose();
}

static void invtrnsl(FILE* fd, Seq* sd)
{
extern	int*	invtranslate[];
extern	float	CodonUsage[];
extern	CHAR	gencode[];
	if (!fd) fd = qout("");
	if (!fd) return;
	float	nfrq[3][4];
	CHAR*	ps = sd->at(sd->left);
	FTYPE	f0 = 1./ sd->many;
	CHAR	nc[4];

	mkinvtab();
	for (int l = sd->left; l < sd->right; ++l) {
#if USE_WEIGHT
	    FTYPE*	w = sd->weight;
#else
	    FTYPE*	w = 0;
#endif
	    for (int n = 0; n < 3; ++n)
		for (int i = 0; i < 4; ++i)
		    nfrq[n][i] = 0;
	    for (int i = 0; i < sd->many; ++i) {
		int*	invtab = invtranslate[*ps++];
		if (w) f0 = *w++;
		while (*invtab >= 0) {
		    FTYPE	f = CodonUsage[*invtab] * f0;
		    de_codon_4(nc, *invtab++);
		    for (int n = 0; n < 3; ++n)
			nfrq[n][nc[n]] += f;
		}
	    }
	    int	c = 0;
	    for (int n = 0; n < 3; ++n)
		c = 4 * c + vmax(nfrq[n], 4) - nfrq[n];
	    for (int n = 0; n < 3; ++n) {
		switch (n) {
		    case 0:	fprintf(fd, "%5d", l + 1); break;
		    case 1:	fprintf(fd, "%5s", amino3[gencode[c]]); break;
		    case 2:	fputs("     ", fd); break;
		}
	        for (int i = 0; i < 4; ++i) 
		    fprintf(fd, "\t%6.2f", 100. * nfrq[n][i]);
		putc('\n', fd);
	    }
	}
	qclose();
}

int utp_main(CalcServer<Seq>* svr, Seq* seqs[], ThQueue<Seq>* q)
{
	Seq*&	sd = seqs[0];
	switch (svr->jobcode) {
	    case 'B':	fouteij(stdout, sd); break;
	    case 'c':
	    case 'C':	composit(stdout, sd); break;
	    case 'f':	findpat(stdout, sd, pattern); break;
	    case 'g':	forgeseq(stdout, sd); break;
	    case 't':
	    case 'T':	invtrnsl(stdout, sd); break;
	    case 'l':
	    case 'L':	proutseq(stdout, sd); break;
	    case 'n':
	    case 'N':	nameseq(stdout, sd); break;
	    case 'M':	mapsite->mapsite(sd, 0); break;
	    case 'm':	mapsite->mapsite(sd, 1); break;
	    case 'x':	muteseq(stdout, seqs); break;
	    default:	usage();
	}
	return (OK);
}

template <class seq_t>
void AlnServer<seq_t>::menu_job()
{
static	const	char prmssg[] =
		"\n#(+)/Comp/Find/List/Map/pRosite/Skl/invT/Quit : ";
	char	cmd[CONLINE];
	algmode.mlt = 2;
	int	c = *restsq(cmd, 0, 1);
	setprompt(1, 0);
	Seq*&	sd = this->in_face[0];
	do {
	  switch (c) {
	    case '1': inputseq(this->in_face, cmd); break;
	    case 'B': fouteij(0, sd); break;
	    case 'c': 
	    case 'C': composit(0, sd); break;
	    case 'j': printskl(0, sd); break;
	    case 'f': find(0, sd); break;
	    case 'g': forgeseq(stdout, sd); break;
	    case 'l':
	    case 'L': list(0, sd); break;
	    case 'p': setparam(1); break;
	    case 'P': setparam(2); break;
	    case 'M': mapsite->mapsite(sd, 0); break;
	    case 'm': mapsite->mapsite(sd, 1); break;
	    case 'r':
	    case 'R': prosite(sd); break;
	    case 'S': shufseq(sd); break;
	    case 't':
	    case 'T': invtrnsl(0, sd); break;
	    case 'x': mute(this->in_face);
		list(0, this->in_face[1]); break;
	    default:  break;
	  }
	  c = *progets(cmd, prmssg);
	} while (c != 'q');
}

int main(int argc, const char** argv)
{
	setdefmolc(PROTEIN);
	alprm2.spb = 1;
	alprm.thr = DefThr;
	setprmode(Row_Last, 'L', SILENT);
	AlnServer<Seq>	svr(argc, argv, IM_SNGL, IM_SNGL,(void*) 0, &utp_main);
	if (svr.autocomp() == 1) svr.menu_job();
	EraDbsDt();
	eraStrPhrases();
	delete mapsite;
	return (0);
}
