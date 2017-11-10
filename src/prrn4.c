/*****************************************************************************
*
*	Refinement of alignment by iteration
*
*	See prrn.doc for details
*
**	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-41-6 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-2012)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#define DEBUG	0

#include "aln.h"
#include "mseq.h"
#include "maln.h"
#include "autocomp.h"
#include "gfreq.h"
#include "mgaps.h"
#include "css.h"
#include "consreg.h"
#include "randiv.h"
#include "fspscore.h"
#include "prrn4.h"

#define MONIT	1

inline	int	icmpf(const int* a, const int* b) {return (*a - *b);}

static	void	sortseq(mSeq* dst, mSeq* src, Subset* ss, int* odr);
static	mSeq*	replseq(mSeq** sqs);
static	void	childfact(Knode* node, double fact);
static	void	calcfact(FTYPE* w, Knode* node);
static	int	n_common(int* a, int* b);
static	void	printlst(FILE* fd, int* lstbuf[]);
static	int*	udmember(int* dl, mSeq* dest, mSeq* ud);
static	int	lcmp(const int* a, const int* b);
static	void	setdefparam();
#if DEBUG
static	void	printset(Subset* ss);
static	void	printwt(mSeq* sd);
#endif

static	int	maxitr = 10;
static	int	grp_code = 2;
static	const	char*	anchor = 0;
static	const	char*	block = 0;
static	char*	groups = 0;
static	int	nseries = 1;
static	INT	countaln = 0;
static	int	rejectup = 1;
static	int	shuffle_odr = 0;
static	int	locpam = 100;
static	float	local_u = 3.;
static	float	local_v = 10.;
static	int	thrrng = 0;
static	int	randseed = 1;
static	int	recycle = 10;
static	FILE*	fmessg = stderr;
static	int	reAlign = ThisAln;
static	int	update = 0;
static	int*	g_locked = 0;
static	int*	u_locked = 0;
static	int	n_locked = 0;
static	float	olr_alpha = 0.1;
static	float	olr_thr = 20.;
static	int	MaxNoGaps = 10;
static	const	char*	guidetree = 0;
//static	const	char*	grouping = 0;

void usage()
{
	fputs("*** prrn version 4.2.0 <160815> ***\n\n", stderr);
	fprintf(stderr, "Usage: ---%s\n", progname);
	fprintf(stderr, "\t%s [-option ...] [input_sequences]\n", progname);
	fputs("Options:\n", stderr);
	fputs("\t-A#\tAlgorithm [0-4]\n", stderr);
	fputs("\t-a\tSequential Prealignment\n", stderr);
	fputs("\t-b\tProgressive Prealignment\n", stderr);
	fputs("\t-c#\tSet significant level of outlier test\n", stderr);
	fputs("\t-E[file name]\tOutput message to\n", stderr);
	fputs("\t-F#\tOutput format [6:PHYLIP;7:GCG;8:CLUSTAL;9:GDE]\n", stderr);
	fputs("\t-H#\tThreshold of conserved region\n", stderr);
	fputs("\t-I#\tMax # of outer loop >= 0\n", stderr);
	fputs("\t-L\tLocal alignment mode\n", stderr);
	fputs("\t-m[file name]\tSubstitution matrix\n", stderr);
	fputs("\t-O#\tSet printmode [0-15]\n", stderr);
	fputs("\t-po\tOnly sumary of outlier tests\n", stderr);
	fputs("\t-ps\tRearrange members according to the tree\n", stderr);
	fputs("\t-o[file name]\tOutput result to\n", stderr);
	fputs("\t-R#\tSeed of random number > 0\n", stderr);
	fputs("\t-s[directory]\tDirectory of sequence files\n", stderr);
	fputs("\t-S#\t# of series of iterations\n", stderr);
	fputs("\t-u#\tGap extension penalty\n", stderr);
	fputs("\t-v#\tGap opening penalty\n", stderr);
	fputs("\t-w#\tBand width off diagonal\n", stderr);
	fputs("\t-U#:\tUpdate mode\n", stderr);
	fputs("\t-Y#\tPAM value for conserved region\n", stderr);
	fputs("\t-?\tThis\n", stderr);
	exit (1);
}

template <class seq_t>
void AlnServer<seq_t>::setparam(int level)
{
	char	str[MAXL];

	if (level & 1) {
	    setalprm();
	    setstrip(QUERY);
	    promptin("# of series ? (%d) : ", &nseries);
	    promptin("Reject up-going cycle ? [1/0] (%d) : ", &rejectup);
	    setalgmode(QUERY, QUERY);
	    promptin("# of iter. cycles (%d) : ", &maxitr);
	    if (setrandiv(QUERY)) setdivseq(QUERY, SILENT, QUERY);
	    promptin("Groups old[0]/new[1]/all[2] (%d) : ", &grp_code);
	    switch (grp_code) {
	      case 1:	progets(str, ": ");
		if (!*str) grp_code = 0;
		else if (*str == '2') {
		    delete[] groups;
		    groups = 0;
		}
		else groups = strrealloc(groups, str);
		break;
	      case 2: delete[] groups;
		groups = 0; break;
	      default:	break;
	    }
	}
	if (level & 2) {
	    setthr(DQUERY);
	    if (alprm.thr) {
		promptin("Local pam (%d) u (%.1f), v(%.1f) : ", 
			&locpam, &local_u, &local_v);
		promptin("Max to go silent (%d) : ", &thrrng);
	    }
	    if (setlpw(QUERY)) {
		setform(QUERY);
		setprmode(QUERY, QUERY, QUERY);
	    }
	}
}

template <class seq_t>
int AlnServer<seq_t>::localoption(int& argc, const char**& argv)
{
	const	char*	opt = argv[0] + 1;
	const	char*	val = argv[0] + 2;
	const	char*	s2;
	int	c = *opt;
	int	n = 1;
extern	int	ndesc_thr;

	switch (c) {
	  case 'E':
	    val = getarg(argc, argv);
	    if (!val || !*val || *val == '-') fmessg = 0;
	    else if (!(fmessg = fopen(val, "w")))
		fatal(no_file, val);
	    break;
	  case 'B': block = val; break;
	  case 'D': if (*val) rejectup = atoi(val); break;
	  case 'G': groups = strrealloc(groups, val); break;
	  case 'I': recycle = *val? atoi(val): 0; break;
	  case 'J':
	    switch (tolower(*val)) {
		case '\0': case 'n': setphyl(NJ_METHOD); break;
		case 'u':setphyl(UPG_METHOD); break;
		case 't':	reAlign = ThisAln; break;
		case 'c':	reAlign = Composition; break;
		case 'a':
		case 'd':	reAlign = DynAln; break;
	    }
	    break;
	  case 'L': algmode.lcl = 15; break;
	  case 'N': if (*val) setpwfact(atof(val)); break;
	  case 'R': randseed = (*val)? atoi(val): time(0); break;
	  case 'S':
	    if (!*val) break;
	    nseries = atoi(val);
	    s2 = strchr(val, '.');
	    if (s2)	maxitr = atoi(s2+1);
	    break;
	  case 'X': 
	    switch (*val) {
		case 'n': ndesc_thr = atoi(val); break;
		case 'o': olr_thr = atof(val); break;
	    }
	    break;
	  case 'Y': 
	    if (isdigit(*val)) {locpam = atoi(val); break;}
	    switch (*val) {
		case 'u':
		case 'U': local_u = atof(val+1); break;
		case 'v':
		case 'V': local_v = atof(val+1); break;
		case 'p':
		case 'P': locpam = atoi(val+1); break;
		default: break;
	    }
	    break;
	  case 'Z': anchor = val; break;
	  case 'U': 
	    update = 1; 
	    this->input_mode = IM_ADON;
	    break;
	  case 'a':
	    this->input_mode = IM_ADON;
	    this->catalog = getarg(argc, argv);
	    break;
	  case 'b':
	    this->calc_mode = IM_TREE;
	    guidetree = getarg(argc, argv);
	    break;
	  case 'c':
	    olr_alpha = atof(getarg(argc, argv)); break;
	  case 'u': case 'v': case 'w': readalprm(opt); break;
	  case 'y': readalprm(val); break;
	  case '?': case 'h': usage();
	  default: n = 0; break;
	}
	return (n);
}

static void sortseq(mSeq* dst, mSeq* src, Subset* ss, int* odr)
{
	int 	n = 0;
	int*	which = new int[src->many + 1];

	for (int i = 0; i < ss->num; ++i) {
	    int*	gp = ss->group[odr? odr[i]: i];
	    while (*gp >= 0) which[*gp++] = n++;
	}
	which[n] = -2;
	dst = src->extseq(dst, which, CPY_SEQ + CPY_SBI + (odr? CPY_LBL: 0));
	delete[] which;
}

static mSeq* replseq(mSeq** sqs)
{
	mSeq*	old = sqs[0];
	mSeq*	nnw = sqs[1];
	mSeq*	tmp = sqs[2];
	RANGE	tail;

	if (old->many != nnw->many) return (0);
	old->saverange(&tail);
	tmp->refresh(old->many);
	old->right = old->left;
	old->left = 0;
	tmp = old->cutseq(tmp, CPY_SEQ | CPY_NBR | CPY_SBI);
	tmp = nnw->catseq(tmp);
	old->left = tail.right;
	old->right = old->len;
	tmp = old->catseq(tmp);
	old->restrange(&tail);
	tmp->left  = tail.left;
	tmp->right = tail.left + nnw->len;
	if (tmp->weight) delete[] tmp->weight;
	if (tmp->pairwt) delete[] tmp->pairwt;
	tmp->weight = old->weight; old->weight = 0;
	tmp->pairwt = old->pairwt; old->pairwt = 0;
	old->copylbl(tmp);
	gswap(sqs[0], sqs[2]);
	return (sqs[0]);
}

void Msa::prntgap(FILE* fd)
{
	GapsList	glist(ss->num);
	takeapart(&glist, msd, false);
	glist.print(fd);
}

void Msa::readgap(FILE* fd)
{
	GapsList	glist(fd);
	if (!glist.num) return;
	if (glist.num != ss->num) {
	    prompt("Improper gaps file!\n");
	    return;
	}

	mSeq**	smem = takeapart(0, msd, true);
	for (int i = 0; i < ss->num; ++i) smem[i]->elim_column(DEL_GAP);
	gatherseq(seqs[1], smem, &glist);
	replseq(seqs);
	clearseq(smem, ss->num);
	delete[] smem;
}

int Msa::checkss()
{
	if (ss == 0 || grp_code != 0) {
	    if (ss) delete ss;
	    ss = new Subset(msd->many, groups);
	    if (!ss || msd->many != ss->elms) 
		return (ERROR);
	    if (block) {
		delete[] g_locked;
		g_locked = new int[msd->many];
		n_locked = getiarray(g_locked, msd->many, block);
		if (n_locked > 1) {
		    qsort((void*) g_locked, (INT) n_locked, sizeof(int), (CMPF) icmpf);
		    int*	pi = g_locked;
		    for (int n = 0; n++ < n_locked; ++pi) (*pi)--;
		    *pi = -1;
		} else {
		    delete[] g_locked;
		    n_locked = 0; g_locked = 0;
		}
	    }
	} else if (ss) {
	    int**	mdf = ss->group;
	    int		n = 0;

	    for (int** grp = mdf; *grp; ++grp) {
		if (**grp >= 0) {
		    *mdf++ = *grp;
		    ++n;
		}
	    }
	    *mdf = 0;
	    ss->num = n;
	}
	return (OK);
}

void Msa::updatesrl(Ssrel* srl)
{
	if (!srl) return;
	gswap(ss, srl->ss); gswap(ktree, srl->ktree);
	gswap(spscr, srl->spscr); gswap(pairwt, srl->pairwt);
	delete srl;
}

VTYPE Msa::phyl_pwt(int dyn_aln)
{
	VTYPE	scr = 0;

        if (!msd->weight) msd->weight = new FTYPE[msd->many];
	delete[] pairwt; delete ktree;
	if (ss->num < ss->elms) {
	    mSeq*	tmp = 0;
	    int**	gr = ss->group;
	    Knode*	lead = new Knode[2 * ss->num - 1];
	    FTYPE*	swt = new FTYPE[ss->elms];
	    for (int k = 0; *gr; ++gr, ++k) {
		tmp = msd->extseq(tmp, *gr, CPY_SEQ + CPY_LBL);
		if (tmp->many > 1) {
		    Ktree subtree(tmp);
		    subtree.calcwt(swt);
		    lead[k].ndesc = subtree.root->ndesc;
		    if (update) lead[k].height = lead[k].res = 0;
		    else {
			lead[k].height = subtree.root->height;
			lead[k].res = subtree.root->res;
		    }
		    int*	gi = *gr;
		    for (int i = 0; *gi >= 0; )
			msd->weight[*gi++] = swt[i++];
		} else {
		    msd->weight[**gr] = 1.;
		    lead[k].ndesc = 1;
		    lead[k].height = lead[k].res = 0;
		}
		lead[k].tid = k;
		lead[k].left = lead[k].right = 0;
	    }
	    delete[] swt;
	    delete tmp;
	    if (dyn_aln) setdivseq((int) reAlign, SILENT, SILENT);
	    ktree = new Ktree(msd, ss, UPG_METHOD, lead);
	    pairwt = ktree->calcpw();
	    if (dyn_aln) setdivseq(POPUP, SILENT, SILENT);
	    gr = ss->group;
	    for (int k = 0; *gr; ++gr, ++k) {
		for (int* gi = *gr; *gi >= 0; ++gi)
		   msd->weight[*gi] *= lead[k].vol;
	    } 
	} else {
	    if (dyn_aln) setdivseq((int) reAlign, SILENT, SILENT);
	    ktree = new Ktree(msd);
	    pairwt = ktree->calcpw();
	    if (dyn_aln) setdivseq(POPUP, SILENT, SILENT);
	    for (int i = 0; i < msd->many; ++i)
		msd->weight[i] = ktree->lead[i].vol;
	}
	return (scr);
}

static FTYPE* wfact;
static void childfact(Knode* node, double fact)
{
	if (node->isleaf()) 
	    wfact[node->tid] = node->vol * fact;
	else {
	    childfact(node->left, fact);
	    childfact(node->right, fact);
	}
}

static void calcfact(FTYPE* w, Knode* node)
{
	Knode*	father = node->parent;
	double	fact = 1. / father->vol;

	wfact = w;
	childfact(node, fact);
	do {
	    if (father->left != node)
		childfact(father->left, fact);
	    else
		childfact(father->right, fact);
	    node = father;
	    fact *= node->cur * node->cur;
	} while ((father = node->parent));
}

static int n_common(int* a, int* b)
{
	int	n = 0;

	while (*a >= 0 && *b >= 0) {
	    if (*a < *b)	++a;
	    else if (*a > *b)	++b;
	    else     {++n; ++a; ++b;}
	}
	return (n);
}

int Prrn::totalNoSeq(int* lst)
{
	int	many = 0;
	while (*lst >= 0)
	    many += slst[*lst++]->many;
	return (many);
}

GAPS* Prrn::gather(int k, int* gr, GapsList* gbuf)
{
	GAPS*	hh = 0;
	mSeq*&	sd = seqs[k];
	if (gr) {
	    glst->copy_to(gbuf, gr);
	    hh = gbuf->delcommongap();
	    if (gbuf->num == 1) {
		mSeq*&	src = slst[*gr];
		if (src->many == 1) src->aliaseq(sd);
		else {
		    src->copyseq(sd);
		    for (int j = 0; j < sd->many; ++j)
			sd->weight[j] *= wlst[*gr];
		}
		sd->sumwt = wlst[*gr];
	    } else {
		mSeq**	sq = sbuf;
		VTYPE*	ww = wbuf;
		for ( ; *gr >= 0; ++gr) {
		    *sq++ = slst[*gr];
		    if (wlst) *ww++ = wlst[*gr];
		}
		*sq = 0;
		aggregate_sb(sd, sbuf, gbuf, wbuf);
	    }
	} else {
	    aggregate_sb(sd, slst, glst, wlst);
	}
	sd->exg_seq(sd->inex.exgl, sd->inex.exgr);
	return (hh);
}

SKL* Prrn::divideseq(Randiv* rdiv, int* lst[], GapsList* gsubs[])
{
	LRAND	mark = rdiv->nextrandiv();
	MCTYPE	mcran = rdiv->mcr->mcrand_now();

	rdiv->bin2lst2(lst, mark);
	calcfact(wlst, srl->ktree->lead + mcran);
	if (totalNoSeq(lst[0]) < totalNoSeq(lst[1]))
	    gswap(lst[0], lst[1]);
	GAPS*	gp[2] =
	    {gather(1, lst[0], gsubs[0]), gather(2, lst[1], gsubs[1])};
	SKL*	skl = gap2skl(gp[0], gp[1]);
	return (skl);
}

VTYPE Prrn::profscore(FTYPE* pw)
{
	FTYPE	scr = 0;
	mSeq*	sqs[2];
	float	u0 = alprm.u0;
	alprm.u0 = 0.;

#if DEBUG
	printf("prof\n");
#endif
	for (int j = 1; slst[j]; ++j) {
	    for (int i = 0; i < j; ++i) {
		sqs[0] = slst[i];
		sqs[1] = slst[j];
		SKL*	skl = gap2skl((*glst)[i], (*glst)[j]);
		PreSpScore	pss(sqs, skl);
		scr += *pw++ * pss.calcSpScore();
		delete[] skl;
	    }
	}
	if (u0 == 0.) return ((VTYPE) scr);
	VTYPE	blank = 0;
	for (int j = 0; slst[j]; ++j) {
	    int	n = 0;
	    mSeq*&	msd = slst[j];
	    mSeqItr	ssi(msd);
	    for (GAPS* gp = (*glst)[j]; gaps_intr(++gp); ) {
		ssi.reset(gp->gps - n);
		blank += ssi.dns->efq * gp->gln;
		n += gp->gln;
	    }
	}
	scr -= u0 * blank;
	gswap(alprm.u0, u0);
	return ((VTYPE) scr);
}

static void printlst(FILE* fd, int* lstbuf[])
{
	for (int i = 0; lstbuf[0][i] >= 0; ++i)
	    fprintf(fd, " %d", lstbuf[0][i]);
	putc(':', fd);
	for (int i = 0; lstbuf[1][i] >= 0; ++i)
	    fprintf(fd, " %d", lstbuf[1][i]);
	putc('\n', fd);
}

VTYPE Prrn::rir(VTYPE prv)
{
	SKL*	skl[2];
	INT	nrep = 0;
	int	maxk = 0;
	int&	num = srl->ss->num;
	VTYPE	scr = 0;
	VTYPE	maxsp = NEG_INT;
	int*	lstbuf[2];
	Randiv	rdiv(seqs[0], srl, TREEDIV, randseed);
static	const	char frm2[] = "!! %3d %4d %8.1lf < %8.1lf (%8.1lf) %2d %4d !!\n";

	lstbuf[0] = new int[num + 1];
	lstbuf[1] = new int[num + 1];

	VTYPE	sps = prv;
	int	k = 0;

	for (int i = 0; i < maxitr; ++i) {
	  for (INT j = 0; j < rdiv.cycle || nrep > 1; ++j) {
	    ++k;
	    glst->copy_to(grsv);  //	reserve	previous glist
// to be parallelized
	    do {
		skl[0] = divideseq(&rdiv, lstbuf, glists + 1);
	    } while (skl[0] == 0);
	    PwdM	pwd(seqs + 1);
	    if (pwd.swp) swapskl(skl[0]);
	    PreSpScore	pss(seqs + 1, &pwd);
	    sps = prv - pss.calcSpScore(skl[0]);
	    Gsinfo	alninf;
	    skl[1] = align2(seqs + 1, &pwd, &scr, &alninf);
	    sps += pss.calcSpScore(skl[1]);
	    if (pwd.swp) {
		if (alninf.sigII && alninf.sigII->lst)
		    alninf.sigII->swaplst(seqs[1]->many, seqs[2]->many);
		swapskl(skl[1]);
		gswap(seqs[1], seqs[2]);
	    }
	    synthgap(glists, skl[1], lstbuf);
	    nrep = lt(prv, sps)? 1: nrep + 1;
	    if (rejectup && lt(sps, prv)) {
		if (rejectup > 1) {
		    if (fmessg) fprintf(fmessg, 
			frm2, k, nrep, (double) prv, (double) sps, 
			(double) scr, seqs[1]->many, seqs[2]->many);
		    if (rejectup > 2) {
			printlst(stdout, lstbuf);
			if (rejectup == 3) {
			    seqs[1]->typeseq(stdout);
			    seqs[2]->typeseq(stdout);
			} else {
#if DEBUG
			    printseqgfq(stdout, seqs[1], 1);
			    printseqgfq(stdout, seqs[2], 1);
#endif
			}
			exit(1);
		    }
		}
		sps = prv;
		grsv->copy_to(glst);
	    } else {
		prv = sps;
		if (!rejectup && i > 1 && lt(maxsp, sps)) {
		    maxsp = sps;
		    maxk = k;
		    glst->copy_to(gopt);
		}
	    }
	    delete[] skl[0]; delete skl[1];
	    if (nrep == rdiv.cycle) goto exitloop;
	  }
/*	  rdiv.mcr->smcrand();*/
	}
exitloop:
	if (!rejectup && lt(sps, maxsp)) {
	    sps = maxsp;
	    nrep = maxk;
	    gopt->copy_to(glst);
	}
	countaln += k;
	delete[] lstbuf[0]; delete[] lstbuf[1];
	return (sps);
}

Prrn::Prrn(mSeq** sqs, Ssrel* trl, VTYPE& ref) 
	: seqs(sqs), srl(trl), slst(0), sbuf(0), wlst(0), wbuf(0), 
	glst(glists[0]), gorg(glists[3]), 
	gopt(glists[4]), gmax(glists[5]), grsv(glists[6])
{
#if MONIT
	long	start = time(0);
#endif
	int&	num = srl->ss->num;
	if (num < 2) return;
	VTYPE	maxsp = VABORT + 1;

	int	kk = num + 1;
	countaln = 0;
	int	k = 0;
	glst = new GapsList(num);
	wlst = new FTYPE[kk];
	wbuf = new FTYPE[kk];
	slst = srl->takeapart(glst, seqs[0], true);
	sbuf = new mSeq*[kk];
	int	maxgaps = seqs[0]->len;
	for (int i = 0; i <  num; ++i) {
	    mSeq*&	sd = slst[i];
#if SSHP
	    sd->inex.sshp = sshpprm? 1: 0;
#endif
	    sd->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
	    sd->convseq(RAWSEQ);
	    if (sd->len < maxgaps) maxgaps = sd->len;
	}
	maxgaps = max(seqs[0]->len - maxgaps + 3, MaxNoGaps);
	glst->resize(maxgaps);
	for (int i = 1; i < NGLIST; ++i)
	    glists[i] = new GapsList(num, maxgaps, true);
	if (g_locked) {
	    u_locked = new int[kk];
	    for (int i = k = 0; i < num; ++i)
   		if (n_common(srl->ss->group[i], g_locked))
		    u_locked[k++] = i;
	    u_locked[k] = -1;
	    n_locked = k;
	}
	glst->copy_to(gorg);
	glst->copy_to(gmax);
	VTYPE	initsp = srl->pairsum_ss(seqs[0], true);
	for (int i = 0; i <  num; ++i) slst[i]->clean();
	for (int i = 0; i < nseries; ++i) {
	    gorg->copy_to(glst);
	    VTYPE	sps = rir(initsp);
	    if (lt(maxsp, sps)) {
		maxsp = sps;
		glst->copy_to(gmax);
	    }
	}
	for (int i = 0; i < num; ++i) slst[i]->elim_column(DEL_GAP);
	delete[] gatherseq(seqs[2], slst, gmax);
	clearseq(slst, num);
	for (int i = 0; i < NGLIST; ++i) delete glists[i];
	if (u_locked) {
	    delete[] u_locked;
	    u_locked = 0;
	}
	sortseq(seqs[1], seqs[2], srl->ss, 0);
	replseq(seqs);
	ref += maxsp - initsp;
#if MONIT
	if (fmessg) {
	  long	stop = time(0);
	  seqs[0]->fphseq(2, fmessg);
	  fprintf(fmessg, "  %8.1lf <-- %8.1lf, %2d grp, %4d rep, %2ld sec\n", 
	    (double) maxsp, (double) initsp, num, countaln, stop - start);
	}
#endif
}

#if DEBUG
static void printset(Subset* ss)
{
	printf("%3d %3d\n", ss->num, ss->elms);
	for (int i = 0; i < ss->num; ++i) {
	    for (int* g = ss->group[i]; *g >= 0; ++g)
		printf("%2d ", *g);
	    putchar('\n');
	}
}

static void printwt(mSeq* sd)
{
	if (sd->weight)
	    for (int i = 0; i < sd->many; ++i)
		printf("%6.3f ", sd->weight[i]);
	putchar('\n');
}

#endif

void Msa::swapprm()
{
	gswap(local_u, alprm.u);
	gswap(local_v, alprm.v);
	gswap(u0, alprm.u0);
	gswap(mtx_no, alprm.mtx_no);
}

int Msa::do_job()
{
	int	oldlen = msd->len;
	int	exgl = msd->inex.exgl;
	int	exgr = msd->inex.exgr;
	VTYPE	refined = 0;
static	const char lclwarn[] = "Error: -L option is not feasible to this msa !\n";
        phyl_pwt(0);
	bool	lclsp = algmode.lcl == 15;
	double	tpwt = fmessg? sumofpwt(): 0;
	if (alprm.thr) {
	    swapprm();
	    Ssrel*	rrl = consregss(seqs);
	    Ssrel*	srl = rrl? rrl: this;
	    if (srl->ss->num < 2) {
		delete rrl;
		return (false);
	    }
	    FTYPE*	tmpwt = msd->saveseqwt();
	    RANGE	svrng;
	    msd->saverange(&svrng);
	    RANGE*	atkrng = srl->consreg(msd, DISSIM);
	    if (fmessg) fprintrng(fmessg, atkrng);
	    INEX	inex = msd->inex;
	    RANGE* rng = lastrng(atkrng);
	    if (lclsp && rng->left == 0 && rng->right == msd->len)
		 fatal(lclwarn);
	    for ( ; rng > atkrng; --rng) {
		if (thrrng && rng->right - rng->left > thrrng) {
		    fprintf(stderr, "Range = %d: ", rng->right - rng->left);
		    int	c = progetc("Are you sure to go ahead? [y/n] : ");
		    if (tolower(c) != 'y') continue;
		}
		msd->restrange(rng);
		Ssrel*	trl = srl->consregss(seqs);
		Ssrel*	url = trl? trl: srl;
		if (!url->pairwt) url->recalcpw();
		swapprm();
		msd->inex.exgl = (svrng.left == rng->left)? exgl: 0;
		msd->inex.exgr = (svrng.right == rng->right)? exgr: 0;
		Prrn	prrn(seqs, url, refined);
		msd->restseqwt(tmpwt);
		swapprm();
		delete trl;
	    }
	    swapprm();
	    delete rrl;
	    msd->restrange(&svrng);
	    msd->right += msd->len - oldlen;
	    msd->inex = inex;
	    delete[] tmpwt;
	    delete[] atkrng;
	} else {
	    if (lclsp && msd->left == 0 && msd->right == msd->len)
		fatal(lclwarn);
	    Prrn	prrn(seqs, this, refined);
	}
	if (fmessg) {
	    VTYPE	newsp = pairsum_ss(msd, true);
	    msd->fphseq(2, fmessg);
	    fprintf(fmessg, "  %8.1lf <-- %8.1lf, %2d grp, %8.2lf %8.4lf\n", 
	    (double) newsp, (double) refined,
	    ss->num, (double) newsp / tpwt, (double) tpwt);
	}
	return (gt(refined, 0));
}

void Msa::individuallen(int* leng)
{
	int*	last = leng + msd->many;

	if (ss) {
	    mSeq*	tmp = 0;
	    for (int** grp = ss->group; *grp; ++grp) {
		if (!(tmp = msd->extseq(tmp, *grp, CPY_SEQ))) continue;
		tmp->elim_column(DEL_GAP);
		*leng++ = tmp->len;
	    }
	    delete tmp;
	    return;
	} 
	CHAR*	sp = msd->at(msd->left);
	CHAR*	t = msd->at(msd->right);
	for ( ; leng < last; ++leng) {
	    *leng = 0;
	    for (CHAR* s = sp++; s < t; s += msd->many)
		if (IsntGap(*s)) (*leng)++;
	}
}

static	int*	seq_lengths = 0;
static	bool	ascend = false;

static	int lcmp(const int *a, const int *b) {
	int	d = seq_lengths[*a] - seq_lengths[*b];
	return (ascend? d: -d);
}

void Msa::stuckupseq(mSeq** argsq)
{
	int	m = 0;
	int	pick[2];
	int	nn  = ss? ss->num: msd->many;
	int*	grp = pick;
	int*	odr = new int[nn];

	for (int i = 0; i < nn; ++i) odr[i] = i;
	ascend = shuffle_odr == SHORTER;
	switch (shuffle_odr) {
	    case AS_INPUT:	break;
	    case SHORTER: case LONGER:
		seq_lengths = new int[nn];
		individuallen(seq_lengths);
		qsort(odr, nn, sizeof(int), (CMPF) lcmp);
		delete[] seq_lengths;
		break;
	    default:
		srand(shuffle_odr);
		for (int i = 0; i < nn; ++i) 
		    gswap(odr[i], odr[rand() % nn]);
		break;
	}
	int i = 0;
	if (argsq) {
	    gswap(seqs[1], argsq[0]);
	    m = seqs[1]->many;
	} else {
	    for ( ; i < nn; ++i) {
		if (ss) grp = ss->group[odr[0]];
		else	{pick[0] = odr[0]; pick[1] = -1;}
		if (*grp >= 0) {	// not empty
		    msd->extseq(seqs[1], grp, CPY_ALL);
		    seqs[1]->elim_column(DEL_GAP);
		    break;
		}
	    }
	}

	for (i = 0; ++i < nn; ) {
	    if (argsq) {
		gswap(seqs[2], argsq[i]);
		m += seqs[2]->many;
	    } else {
		if (ss) grp = ss->group[odr[i]];
		else	pick[0] = odr[i];
		if (*grp < 0) continue;
		msd->extseq(seqs[2], grp, CPY_ALL);
		seqs[2]->elim_column(DEL_GAP);
	    }
	    msd = align2(seqs + 1, msd);
	    if (argsq) {
	        if (i == 1) gswap(seqs[1], argsq[0]);
		gswap(seqs[2], argsq[i]);
	    }
	    gswap(seqs[0], seqs[1]);
	}
	if (shuffle_odr && OutPrm.sortodr == AS_INPUT) {
	    int*	rdo = new int[msd->many + 1];
	    for (int j = i = 0; i < nn; ++i) {
		if (ss) {
		    for (grp = ss->group[odr[i]]; *grp >= 0; ++grp)
			rdo[*grp] = j++;
		} else	rdo[odr[i]] = i;
	    }
	    rdo[msd->many] = -1;
	    seqs[1]->extseq(msd, rdo, CPY_ALL);
	    delete[] rdo;
	} else
	    gswap(seqs[0], seqs[1]);
	msd->inex.algn = 1;
	seqs[1]->refresh();
	delete[] odr;
}

mSeq* Msa::prog_up(Knode* node, mSeq** argsq)
{
	mSeq*	sd = 0;

	if (node->isleaf()) {
	    if (argsq) sd = argsq[node->tid];
	    else {
		sd = msd->extseq(0, ss->group[node->tid], CPY_ALL + RLT_SBI);
		sd->elim_column(DEL_GAP);
	    }
	    lst_odr[lst_idx++] = node->tid;
	} else {
	    mSeq*	sqs[3];
	    sqs[0] = prog_up(node->left, argsq);
	    sqs[1] = prog_up(node->right, argsq);
	    sqs[2] = 0;
	    sd = align2(sqs);
	    clearseq(sqs, 2);
	}
	return (sd);
}

void Msa::progressive()
{
	int	w = alprm.sh;

	lst_idx = 0;
	lst_odr = new int[ss->num + 1];
	alprm.sh += alprm.sh;
	phyl_pwt(reAlign);
	delete seqs[2];
	seqs[2] = prog_up(ktree->root, 0);
	lst_odr[lst_idx] = -1;
	if (seqs[2]->sigII) seqs[2]->sigII->renumlst(lst_odr);
	sortseq(seqs[0], seqs[2], ss, lst_odr);
	delete[] lst_odr;
	lst_odr = 0;
	alprm.sh = w;
	msd->inex.algn = 1;
}

template<>
mSeq* ProgMsa<Knode>::getseq(Knode* node)
{
	lst_odr[lst_idx++] = node->tid;
	return (seqs[node->tid]->aliaseq(0));
}

template<>
mSeq* ProgMsa<Tnode>::getseq(Tnode* node)
{
	mSeq*	sd = new mSeq;
	return (sd->getseq(node->tname));
}

void Msa::lstodr(Knode* node)
{
	if (node->isleaf())
	    lst_odr[lst_idx++] = node->tid;
	else {
	    lstodr(node->left);
	    lstodr(node->right);
	}
}

static float* flength = 0;
static int fpcmp(const int* a, const int* b)
{
	if (flength[*a] == flength[*b]) return (0);
	return (flength[*a] > flength[*b])? 1: -1;
}

Outlier* Msa::findoutliers(float* stddev)
{
	Outlier*	outlier = new Outlier[msd->many];
	vclear(outlier, msd->many);

//	lone eijunction
	if (msd->sigII) {
	    msd->sigII->mkeijtab(msd->many);
	    for (int i = 0; i < msd->many; ++i)
		outlier[i].eij = msd->sigII->lone[i];
	}

//	unusual amount of indels
	if (!ktree) phyl_pwt(0);
	Ssrel*	rrl = consregss(seqs);
	Ssrel*	srl = rrl? rrl: this;
	if (srl->ss->num > 1) {
	    RANGE	tmp;
	    msd->saverange(&tmp);
	    FTYPE*	lclsod = 0;
	    FTYPE*	glbsod = calcdistsum((Seq*) msd);
	    flength = new float[msd->many];
	    int*	order = new int[msd->many];
	    int*	olist = new int[msd->many];
	    swapprm();
	    gswap(olr_thr, alprm.thr);
	    RANGE*	atkrng = srl->consreg(msd, DISSIM);
	    RANGE*	rng = atkrng;
	    RANGE*	lst = atkrng + atkrng->left - 1;
	    int	last = atkrng->left - 3;
	    Dixon	dxn(olr_alpha);
	    for (int i = 0; ++rng < lst; ++i) {
		CHAR*	sp = msd->at(rng->left);
		CHAR*	tt = msd->at(rng->right);
		for (int j = 0; j < msd->many; ++j, ++sp) {
		    order[j] = j;
		    flength[j] = 0;
		    for (CHAR* s = sp; s < tt; s += msd->many)
			if (IsntGap(*s)) ++flength[j];
		}
		qsort(order, msd->many, sizeof(int), (CMPF) fpcmp);
		int	nol = dxn.dixon(olist, flength, order, msd->many, 2);
		for (int k = 0; k < nol; ++k) {
		    if (i == 0) {
			if (olist[k] >= 0) ++outlier[olist[k]].ins_f;
			else	++outlier[-olist[k] - 1].del_f;
		    } else if (i == last) {
			if (olist[k] >= 0) ++outlier[olist[k]].ins_l;
			else	++outlier[-olist[k] - 1].del_l;
		    } else {
			if (olist[k] >= 0) ++outlier[olist[k]].ins_m;
			else	++outlier[-olist[k] - 1].del_m;
		    }
		}
		float	x;
		vavsd(x, flength, msd->many);
		stddev[0] += 1;
		stddev[1] += x * x;
		if (last) {
		    if (i == 0) stddev[2] = x;
		    if (i == last) stddev[3] = x;
		} else
		    stddev[2] = stddev[3] = -1.;

//	unusually divergent seq
		msd->restrange(rng);
		lclsod = calcdistsum((Seq*) msd, 0, lclsod);
		for (int j = 0; j < msd->many; ++j) {
		    order[j] = j;
		    flength[j] = lclsod[j] / glbsod[j];
		}
		qsort(order, msd->many, sizeof(int), (CMPF) fpcmp);
		nol = dxn.dixon(olist, flength, order, msd->many);
		for (int k = 0; k < nol; ++k)
		    if (olist[k] >= 0) ++outlier[olist[k]].match;
	    }
	    gswap(olr_thr, alprm.thr);
	    swapprm();
	    delete[] atkrng;
	    delete[] olist;
	    delete[] order;
	    delete[] flength; flength = 0;
	    delete[] lclsod;
	    delete[] glbsod;
	    msd->restrange(&tmp);
	}
	delete rrl;
	return outlier;
}

VTYPE Msa::nomal_pairwt()
{
	int	nn = ncomb(msd->many);
	VTYPE	f = 0.;
	if (pairwt)
	    for (int i = 0; i < nn; ++i)
		f += pairwt[i];
	else
	    f = (VTYPE) nn;
	return (f);
}

FTYPE Msa::sp_score(FSTAT* fst)
{
	vclear(fst, 2);
	FTYPE*	pwt = pairwt;
	FTYPE	sumpwt = 0;
	Simmtx*	sm = getSimmtx(0);
	for (int j = 1; j < msd->many; ++j) {
	    for (int i = 0; i < j; ++i) {
		CHAR*	ss = msd->at(msd->left);
		CHAR*	ts = msd->at(msd->right);
		int	agl = 0;
		int	bgl = 0;
		for ( ; ss < ts; ss += msd->many) {
		    bool	ag = IsGap(ss[i]);
		    bool	bg = IsGap(ss[j]);
		    if (ag && bg) continue;
		    if (ag && !bg) {
			++fst->unp;
			if (pwt) fst[1].unp += *pwt;
			if (agl <= bgl) {
			    ++fst->gap;
			    if (pwt) fst[1].gap += *pwt;
			}
			++agl; bgl = 0;
		    } else if (!ag && bg) {
			++fst->unp;
			if (pwt) fst[1].unp += *pwt;
			if (agl >= bgl) {
			    ++fst->gap;
			    if (pwt) fst[1].gap += *pwt;
			}
			++bgl; agl = 0;
		    } else {
			if (ss[i] == ss[j]) {
			    ++fst->mch;
			    if (pwt) fst[1].mch += *pwt;
			} else {
			    ++fst->mmc;
			    if (pwt) fst[1].mmc += *pwt;
			}
			fst->val += sm->mtx[ss[i]][ss[j]];
			if (pwt) fst[1].val += *pwt * sm->mtx[ss[i]][ss[j]];
			agl = bgl = 0;
		    }
		}
		if (pwt) sumpwt += *pwt++;
	    }
	}
	fst->val -= alprm.v * fst->gap + alprm.u * fst->unp;
	if (pwt) fst[1].val -= alprm.v * fst[1].gap + alprm.u * fst[1].unp;
	return (sumpwt);
}
		    
void Msa::prrn_main(AlnServer<mSeq>* svr)
{
	RANGE*	givrng = 0;
	RANGE	svrng;
	int	exgl = msd->inex.exgl;
	int	exgr = msd->inex.exgr;

	if (alprm.thr == INT_MAX) {
	    if (msd->isdrna()) {
		locpam = 2;
		setthr(5.);
	    }  else
		setthr(20.);
	}
	if (checkss() == ERROR)
	    fatal("Error -- Improper grouping!");
	if (ss->num < 2) 
	    fatal("Error -- Give multiple sequences!");
	msd->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
	if (!msd->inex.algn) {
	  switch (svr->input_mode) {
	    case IM_ADON:
		if (reAlign >= 0) stuckupseq(0); break;
	    case IM_TREE:
		if (reAlign >= 0) progressive(); break;
	    default:
		if (!msd->test_aligned()) {
		    if (reAlign != DynAln) reAlign = Composition;
		    progressive();
		}
		break;
	  }
	}
	if (anchor) givrng = getrng(anchor);
	Strlist*	spath_rsv = 0;
	gswap(msd->spath, spath_rsv);
	while (update-- >= 0) {
	  for (int rc = 0; rc < recycle; ++rc) {
	    clear();
	    bool	rv = false;		// revised?
	    if (givrng) {
		INEX	inex = msd->inex;
		int	oldlen = msd->len;
		msd->saverange(&svrng);
		for (RANGE* rng = lastrng(givrng); rng > givrng; --rng) {
		    int	tmplen = msd->len;
		    msd->restrange(rng);
		    msd->inex.exgl = (svrng.left == rng->left)? exgl: 0;
		    msd->inex.exgr = (svrng.right == rng->right)? exgr: 0;
		    rv = do_job() || rv;
		    msd->saverange(rng);
		    tmplen -= msd->len;
		    for (RANGE* tmp = rng + 1; tmp <= lastrng(givrng); ++tmp) {
			tmp->left -= tmplen;
			tmp->right -= tmplen;
		    }
		}
		msd->restrange(&svrng);
		msd->inex = inex;
		msd->right += msd->len - oldlen;
	    } else 
		rv = do_job();
	    if (!rv || ss->num == 2) break;
	  }
	  if (update == 0) {
	    grp_code = 2;
	    if (checkss() == ERROR)
		fatal("Error -- Improper grouping!");
	  }
	}
	gswap(msd->spath, spath_rsv);
	if (givrng) {
	    if (fmessg) fprintrng(fmessg, givrng);
	    delete[] givrng;
	}
	if (out_fd) {
	    if (OutPrm.sortodr == BY_TREE) {
		clear();
		ktree = new Ktree(msd);
		lst_odr = new int[msd->many + 1];
		lst_idx = 0;
		lstodr(ktree->root);
		lst_odr[lst_idx] = -1;
		msd->extseq(seqs[2], lst_odr, CPY_ALL + RLT_SBI);
		gswap(seqs[0], seqs[2]);
		delete[] lst_odr;
	    }
	    if (algmode.nsa & 1) {
		if (getlpw()) msd->typeseq();
		else if (ss || checkss() != ERROR) prntgap(out_fd);
	    }
	    float stdiv[4] = {0, 0, 0, 0};
	    if (algmode.nsa & 2) {
		Outlier*	outlier = findoutliers(stdiv);
		Outlier*	olr = outlier;
		FILE*	fo = stdout;
		int*	acc = new int[10];
		vclear(acc, 10);
		if (OutPrm.out_file) {
		    char	str[MAXL];
		    strcpy(str, OutPrm.out_file);
		    strcat(str, OLR_EXT);
		    fo = fopen(str, "w");
		    if (!fo) fo = stdout;
		}
		for (int i = 0; i < msd->many; ++i, ++olr) {
		    bool	def = olr->match || olr->ins_f || olr->del_f ||
			olr->ins_m || olr->del_m || olr->ins_l || olr->del_l;
		    acc[0] += def; acc[1]+= olr->eij; acc[2] += olr->match;
		    acc[3] += olr->ins_f; acc[4] += olr->del_f; acc[5] += olr->ins_m;
		    acc[6] += olr->del_m; acc[7] += olr->ins_l; acc[8] += olr->del_l;
		    if (!OutPrm.olrsum) fprintf(fo, "%5d %-23s\t%3d %2d %d %d %d %d %d %d %d\n",
			i + 1, (*msd->sname)[i], def, 
			olr->eij, olr->match, olr->ins_f, olr->del_f, 
			olr->ins_m, olr->del_m, olr->ins_l, olr->del_l);
		}
		if (OutPrm.olrsum) {
		    fprintf(fo, "%-15s %3d", msd->sqname(), msd->many);
		    for (int i = 0; i < 9; ++i) fprintf(fo, " %3d", acc[i]);
		    if (stdiv [0] > 1) stdiv[1] /= stdiv [0];
		    if (stdiv [1] > 0) stdiv[1] = sqrt(stdiv[1]);
		    fprintf(fo, " %2d", int(stdiv[0]));
		    for (int i = 1; i < 4; ++i) fprintf(fo, " %5.1f", stdiv[i]);
		    fputc('\n', fo);
		}
		delete[] outlier; delete[] acc;
		if (fo != stdout) fclose(fo);
	    }
	    int	span = msd->right - msd->left;
	    if (algmode.nsa & 4 && span > 0) {
		phyl_pwt(0);
		VTYPE	spscr = pairsum_ss(msd, false);
		VTYPE	wspscr = pairsum_ss(msd, true);
		fprintf(out_fd, "%s [ %d ] %d\t%7.1f %7.3f %7.1f %7.3f\n", 
		    msd->sqname(), msd->many, span, 
		    float(spscr), float(100. * spscr / ncomb(msd->many) / span),
		    float(wspscr), float(100. * wspscr / nomal_pairwt() / span));
	    }
	}
}

static int* udmember(int* dl, mSeq* dest, mSeq* ud)
{
	for (int j = 0; j < ud->many; ++j) {
	    for (int i = 0; i < dest->many; ++i) {
		if (!strcmp((*ud->sname)[j], (*dest->sname)[i]))
		    *dl++ = i;
	    }
	}
	*dl = -1;
	return (dl);
}

void Msa::makemsa(const char* gtree)
{
	Btree<Tnode> ltree(gtree);
	Seq	sd(1);
	sd.getseq(ltree.mname[0]);
	setdefmolc(sd.inex.molc);
	prePwd(sd.inex.molc);
	ProgMsa<Tnode>	promsa(ltree, gtree);
	msd = promsa.prog_up(ltree.root);
}

Msa::Msa(AlnServer<mSeq>& svr) : Ssrel(), msd(seqs[0])
{
	u0 = 0.; mtx_no = 1;
	initseq(seqs, 3); seqs[3] = 0;
	if (svr.calc_mode == IM_TREE) {
	    makemsa(guidetree);
	} else if (svr.sql1->nextseq(seqs) != IS_OK) usage();
	prePwd(msd);
	INT	molc = seqs[0]->inex.molc;
	if (svr.memb->membersize <= 1) return;	// single input

	mSeq**	tmp = new mSeq*[svr.memb->membersize + 1];
	initseq(tmp, svr.memb->membersize + 1);
	ss = new Subset(svr.memb->membersize, (const char*) 0);
	int**	grp = ss->group;
	gswap(seqs[0], tmp[0]);
	int	maxmem = ss->elms = tmp[0]->many;
	int	n = 0;
	int	nerror = 0;
	InSt	ist = IS_OK;
	while ((ist = svr.sql1->nextseq(tmp + ++n)) != IS_END) {
	    if (ist == IS_ERR || tmp[n]->inex.molc != molc || tmp[n]->len == 0) {
		++nerror;
		--ss->num;
		delete tmp[n];
		tmp[n] = 0;
	    } else {
		ss->elms += tmp[n]->many;
		if (tmp[n]->many > maxmem) maxmem = tmp[n]->many;
	    }
	}
	if (update) {				// remove older sequences
	    setdefmolc(seqs[0]->inex.molc);
	    int*	updlist = new int[maxmem + 1];
	    int*	rsvlist = new int[maxmem + 1];
	    int	i = n;
	    while (--i >= 0 && !tmp[i]) ;
	    if (tmp[i]) ss->elms = tmp[i]->many;
	    while (--i >= 0) {
		if (!tmp[i]) continue;
		int*	upd = updlist;
		for (int j = n - 1; j > i; --j)	// find member to be updated
		    if (tmp[j]) upd = udmember(upd, tmp[i], tmp[j]);
		if (int k = upd - updlist) {	// found
		    *upd = -1;
		    if (k > 1) qsort((UPTR) updlist, k, sizeof(int), (CMPF) icmpf);
		    int* rsv = upd = updlist;
		    while (*++upd >= 0)
			if (*rsv != *upd) *++rsv = *upd; // remove redundcy
		    *++rsv = -1;
		    rsv = rsvlist;
		    upd = updlist;
		    for (int j = 0; j < tmp[i]->many; ++j) {
			if (j == *upd) ++upd;
			else	*rsv++ = j;
		    }
		    *rsv = -1;
		    if (int k = rsv - rsvlist) {
			ss->elms += k;
			tmp[n] = tmp[i]->extseq(tmp[n], rsvlist, CPY_ALL + RLT_SBI);
			tmp[n]->elim_column(DEL_GAP);
			gswap(tmp[i], tmp[n]);
			delete tmp[n]; tmp[n] = 0;
		    } else {
			delete tmp[i]; tmp[i] = 0;
			ss->num--;
			nerror++;
		    }
		} else {
		    ss->elms += tmp[i]->many;
		}
	    }
	    delete[] rsvlist; delete[] updlist;
	}
	if (nerror) {
	    int	j = 0;
	    for (int i = 0; i < n; ++i)
		if (tmp[i]) tmp[j++] = tmp[i];
	    tmp[n = j] = 0;
	}
	delete[] ss->pool;
	int*	pi = ss->pool = new int[ss->elms + ss->num];
	for (int i = 0, k = 0; i < n; ++i) {
	    *grp++ = pi;
	    for (int j = 0; j < tmp[i]->many; ++j) *pi++ = k++;
	    *pi++ = EOTAB;
	}
	*grp = 0;
	stuckupseq(tmp);
	clearseq(tmp, svr.memb->membersize + 1);
	delete[] tmp;
	grp_code = 0;
}

static void setdefparam()
{
	treemode.n_edge = 1;
	optimize(GLOBAL, MAXIMUM);
	setlsegs(1);
	setalgmode(0, 0);
	setdefPprm(250, 2., 9.);
	setdefNprm(-2., 2., 4.);	// n, u, v
	setprmode(Row_Last, 'L', SILENT);
	algmode.mns = 1;
	algmode.nsa = 1;
	OutPrm.SkipLongGap = 0;	// long gaps are reported
}

int main(int argc, const char** argv)
{
	progname = *argv;
	setdefparam();
	AlnServer<mSeq>	svr(argc, argv, IM_MULT, IM_MULT);
	setdefPprm(SILENT, alprm.u, alprm.v, 0);
	setdefPprm(locpam, local_u, local_v, 1);
#if SSHP
	initSsHpPrm();
#endif
	spb_fact();
	Msa	msa(svr);
	msa.prrn_main(&svr);
	delete[] g_locked;
#if SSHP
	eraseSsHpPrm();
#endif
	EraDbsDt();
	resetSimmtxes();
	return (0);
}

