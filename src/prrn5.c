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

#include "mseq.h"
#include "sltree.h"
#include "maln.h"
#include "autocomp.h"
#include "gfreq.h"
#include "mgaps.h"
#include "css.h"
#include "consreg.h"
#include "randiv.h"
#include "fspscore.h"
#include "prrn5.h"

#define MONIT	1

inline	int	icmpf(const int* a, const int* b) {return (*a - *b);}

static	void	sortseq(mSeq* dst, mSeq* src, Subset* ss, int* odr);
static	mSeq*	replseq(mSeq** sqs);
static	void	setdefparam();
static	int	slfprrn_main(CalcServer<Sltree>* svr, Sltree** nodes, ThQueue<Sltree>* q);
static	mSeq*	makemsa(mSeq** seqs, int nn, ALPRM* alp, bool sglthrd);
static	mSeq*	makemsa(int argc, const char** argv, int num);
static	mSeq*	cut_in(mSeq* mom, mSeq* dau);
static	void	print_msa(mSeq* sd);

static	const	int	eff_paral_num = 4;
static	const	int	ndescthr = 1000;

static	int	maxitr = 10;
static	char*	groups = 0;
static	int	DefPrrnPam = 250;
static	int	DefScndPam = 100;
static	float	DefMmcPen = -2.;
static	float	local_u = 3.;
static	float	local_v = 10.;
static	float	local_thr = 35.;
static	int	dist_sh = -60;
static	int	randseed = 1;
static	int	recycle = 1;
static	int	update = 0;
static	float	olr_alpha = 0.1;
static	float	olr_thr = 20.;
static	int	min_host_size = 2;
static	const	char*	header = 0;
static	const	char*	catalog = 0;
static	const	char*	guidetree = 0;
static	RunStat	runstat; 
static	int	paral_num = -1;
static	int	MaxNoGaps = 10;
static	int	max_memb = INT_MAX;
static	int	min_memb = 2;
static	int	msa_id = 0;
static	int	min_ndesc = 5;
static	int	min_seqs = 16;
static	DivMode	divmode = TREEDIV;
static	const	char*	fail_to_align = "Fail to align!\n";
static	bool	incl_outsider = true;

static	int	fget(FILE* fd, const char* fn = 0) {return 1;}

void usage()
{
	fputs("*** prrn version 5.1.0 <170215> ***\n\n", stderr);
	fprintf(stderr, "Usage: ---%s\n", progname);
	fprintf(stderr, "\t%s [-option ...] [input_sequences]\n", progname);
	fputs("Options:\n", stderr);
	fputs("\t-D#\tMaxmum distance to connect nodes\n", stderr);
	fputs("\t-F#\tOutput format:\n", stderr);
	fputs("\t\t[f:FASTA;m:PIR;n:NEXUS;6:PHYLIP;7:GCG;8:CLUSTAL;9:GDE]\n", stderr);
	fputs("\t\t[C:column number;M:column marker;N:sequence number at right\n", stderr);
	fputs("\t-H#\tThreshold of conserved region\n", stderr);
	fputs("\t-I#\tMax # of outer loop >= 0\n", stderr);
	fputs("\t-K[A|D] specify amino acid or DNA sequence\n", stderr);
	fputs("\t-L\tSemi global alignment mode\n", stderr);
	fputs("\t-O#\tSet printmode [0-15]\n", stderr);
	fputs("\t-R#\tSeed of random number > 0\n", stderr);
	fputs("\t-S#\t# of series of iterations\n", stderr);
	fputs("\t-U\tUpdate mode\n", stderr);
	fputs("\t-X#\tClump size\n", stderr);
	fputs("\t-Y[u|v|p]#\t alignment parameters for conserved regions\n", stderr);
	fputs("\t-b\tGiven guide tree\n", stderr);
	fputs("\t-c#\tSet significant level of outlier test\n", stderr);
	fputs("\t-d#\tDivide mode [1:tree leaves;2:tree edges]\n", stderr);
	fputs("\t-h[Header]\tPrefix of sub-alignments\n", stderr);
	fputs("\t-m[file name]\tSubstitution matrix\n", stderr);
	fputs("\t-o[file name]\tOutput result to\n", stderr);
	fputs("\t-ph\tOutput html format with intron position markers\n", stderr);
	fputs("\t-pi\tDisplay intron position on terminal\n", stderr);
	fputs("\t-po\tOnly sumary of outlier tests\n", stderr);
	fputs("\t-pq\tSuppress run time messages\n", stderr);
	fputs("\t-ps\tRearrange members according to the tree\n", stderr);
	fputs("\t-r#\t# of divisions to refine in parallel\n", stderr);
	fputs("\t-s[directory]\tDirectory of sequence files\n", stderr);
	fputs("\t-t#\t# of threads\n", stderr);
	fputs("\t-u#\tGap extension penalty\n", stderr);
	fputs("\t-v#\tGap opening penalty\n", stderr);
	fputs("\t-w#\tBand width off diagonal\n", stderr);
	fputs("\t-yp#\t PAM level\n", stderr);
	fputs("\t-?\tThis\n", stderr);
	exit (1);
}

template <class seq_t>
int AlnServer<seq_t>::localoption(int& argc, const char**& argv)
{
extern	int	ndesc_thr;
	const	char*	opt = argv[0] + 1;
	const	char*	val = argv[0] + 2;
	int	n = 1;
	int	temp;

	switch (*opt) {
	  case 'A':
	    if ((val = getarg(argc, argv, true)))
		algmode.any = atoi(val);
		break;
	  case 'D':
	    if ((val = getarg(argc, argv, true))) 
		set_max_height((FTYPE) atof(val));
	    break;
	  case 'E':
	    runstat.setfmessg(argc, argv);
	    break;
	  case 'G':
	    if ((val = getarg(argc, argv)))
		groups = strrealloc(groups, val);
	    break;
	  case 'I':
	    recycle = (val = getarg(argc, argv, true))? atoi(val): 0;
	    break;
	  case 'J':
	    temp = (val = getarg(argc, argv, true))? atoi(val): 1;
	    if (0 <= temp && temp < 4) divmode = (DivMode) temp;
	    break;
	  case 'L': 
	    algmode.lcl = (val = getarg(argc, argv, true))? atoi(val): 15;
	    break;
	  case 'R':
	    randseed = (val = getarg(argc, argv, true))? atoi(val): time(0);
	    break;
	  case 'S':
	    if ((val = getarg(argc, argv, true))) maxitr = atoi(val);
	    break;
	  case 'U': update = 1; break;
	  case 'X': 
	    switch (*val) {
		case 'n': ndesc_thr = atoi(val); break;
		case 'o': olr_thr = atof(val); break;
		default: 
		    if ((val = getarg(argc, argv, true))) max_memb = atoi(val);
		    break;
	    }
	    break;
	  case 'Y': 
	    switch (*val) {
		case 'H': local_thr = atof(val+1); break;
		case 'u': local_u = atof(val+1); break;
		case 'v': local_v = atof(val+1); break;
		case 's': dist_sh = atoi(val+1); break;
		default: break;
	    }
	    break;
	  case 'b':
	    guidetree = getarg(argc, argv);
	    break;
	  case 'c':
	    olr_alpha = atof(getarg(argc, argv, true)); break;
	  case 'e':
	    header = (val = getarg(argc, argv))? val: def_header;
	    break;
	  case 'l':
	    if ((val = getarg(argc, argv, true))) setlpw(atoi(val));
	    break;
	  case 'o':
	    OutPrm.out_file = ((val = getarg(argc, argv)))? val: "";
	    break;
	  case 'r':
	    paral_num = (val = getarg(argc, argv, true))? atoi(val): INT_MAX;
	    break; 
	  case 'h': case '?': usage();
	  default: n = 0; break;
	}
	return (n);
}

void RunStat::stamp(int val)
{
	time_t	t = time(0);
	long	intvl = timepoint? t - previoust: 0;
	prompt("%d\t%.2f\t%d\n", timepoint, float(intvl) / 60, val);
	previoust = t;
	if (timepoint < N_Stamp) {
	    values[timepoint] = val;
	    timestamp[timepoint] = t;
	}
	++timepoint;
}

void RunStat::conclude()
{
	if (fmessg && timepoint) {
	    for (int i = 1; i < timepoint; ++i)
		fprintf(fmessg, "%ld\t", timestamp[i] - timestamp[i - 1]);
	    long	secs = timestamp[--timepoint] - timestamp[0];
	    float	mins = float(secs) / 60;
	    fprintf(fmessg, "%ld secs %.2f mins\n", secs, mins);
	}
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
#if USE_WEIGHT
	if (tmp->weight) delete[] tmp->weight;
	if (tmp->pairwt) delete[] tmp->pairwt;
	tmp->weight = old->weight; old->weight = 0;
	tmp->pairwt = old->pairwt; old->pairwt = 0;
#endif	// USE_WEIGHT
	old->copylbl(tmp);
	std::swap(sqs[0], sqs[2]);
	return (sqs[0]);
}

void IterMsa::prntgap(FILE* fd)
{
	GapsList	glist(ss->num);
	takeapart(&glist, msd, false);
	glist.print(fd);
}

void IterMsa::readgap(FILE* fd)
{
	GapsList	glist(fd);
	if (!glist.num) return;
	if (glist.num != ss->num) {
	    prompt("Improper gaps file!\n");
	    return;
	}

	mSeq**	smem = takeapart(0, msd, true);
	for (int i = 0; i < ss->num; ++i) smem[i]->elim_column(DEL_GAP);
	gatherseq(infc[1], smem, &glist);
	replseq(infc);
	clearseq(smem, ss->num);
	delete[] smem;
}

int IterMsa::checkss()
{
	if (ss == 0) {
	    ss = new Subset(msd->many, groups);
	    if (!ss || msd->many != ss->elms) 
		return (ERROR);
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

void IterMsa::phyl_pwt()
{
#if USE_WEIGHT
	if (!msd->weight) msd->weight = new VTYPE[msd->many];
	if (msd->many == 1) {
	    msd->weight[0] = 1;
	    return;
	}
	delete[] pairwt; pairwt = 0;
	delete ktree; ktree = 0;
	if (ss->num < ss->elms) {
	    int**	gr = ss->group;
	    Knode*	lead = new Knode[2 * ss->num - 1];
	    mSeq*	tmp;
	    VTYPE*	swt = new VTYPE[ss->elms];
	    for (int k = 0; *gr; ++gr, ++k) {
		int	l = lst_odr? lst_odr[k]: k;
		tmp = seqs? seqs[l]: msd->extseq(tmp, *gr, CPY_SEQ + CPY_LBL);
		if (tmp->many > 1) {
		    if (!tmp->weight) {
			Ktree subtree(tmp);
			subtree.calcwt(swt);
			tmp->height = subtree.root->height;
			tmp->resist = subtree.root->res;
		    }
		    VTYPE*	wt = tmp->weight? tmp->weight: swt;
		    for (int* gi = *gr; *gi >= 0; )
			msd->weight[*gi++] = *wt++;
		} else {
		    msd->weight[**gr] = 1.;
		}
		if (update) lead[k].height = lead[k].res = 0;
		else {
		    lead[k].height = tmp->height;
		    lead[k].res = tmp->resist;
		}
		lead[k].tid = k;
		lead[k].ndesc = tmp->many;
		lead[k].left = lead[k].right = 0;
	    }
	    delete[] swt;
	    if (!seqs) delete tmp;
	    ktree = new Ktree(msd, ss, UPG_METHOD, lead);
	    pairwt = ktree->calcpw();
	    gr = ss->group;
	    for (int k = 0; *gr; ++gr, ++k) {
		for (int* gi = *gr; *gi >= 0; ++gi)
		   msd->weight[*gi] *= lead[k].vol;
	    } 
	    msd->height = ktree->root->height;
	    msd->resist = ktree->root->res;
	} else {
	    ktree = new Ktree(msd);
	    pairwt = ktree->calcpw(msd->weight);
	    msd->height = ktree->root->height;
	    msd->resist = ktree->root->res;
	}
	msd->sumwt = 0;
	for (int i = 0; i < msd->many; ++i) msd->sumwt += msd->weight[i];
#else
	if (msd->many == 1) return;
	delete ktree; ktree = 0;
	if (ss->num < ss->elms) {
	    int**	gr = ss->group;
	    Knode*	lead = new Knode[2 * ss->num - 1];
	    for (int k = 0; *gr; ++gr, ++k) {
		lead[k].height = lead[k].res = 0;
		lead[k].tid = k;
		int	many = 0;
		for (int* gi = *gr; *gi >= 0; ++gi) ++many;
		lead[k].ndesc = many;
		lead[k].left = lead[k].right = 0;
	    }
	    ktree = new Ktree(msd, ss, UPG_METHOD, lead);
	} else {
	    ktree = new Ktree(msd);
	    msd->height = msd->resist = 0;
	}
#endif // USE_WEIGHT
}

#if USE_WEIGHT
void Prrn::childfact(Knode* node, double fact)
{
	if (node->isleaf()) 
	    wfact[node->tid] = node->vol * fact;
	else {
	    childfact(node->left, fact);
	    childfact(node->right, fact);
	}
}

VTYPE Prrn::calcfact(VTYPE* w, Knode* node)
{
	wfact = w;
	VTYPE	rv = node->cur;
	childfact(node, 1. / node->vol);
	FTYPE	fact = 1.;

	while (Knode* father = node->parent) {
	    if (father->left != node)
		childfact(father->left, fact / father->vol);
	    else
		childfact(father->right, fact / father->vol);
	    node = father;
	    fact *= node->cur;
	}
	return (rv);
}
#endif

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
#if USE_WEIGHT
	    if (gbuf->num == 1) {
		slst[*gr]->aliaseq(sd);
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
#else
	    if (gbuf->num == 1) {
		mSeq*&	src = slst[*gr];
		if (src->many == 1) src->aliaseq(sd);
		else src->copyseq(sd);
	    } else {
		mSeq**	sq = sbuf;
		for ( ; *gr >= 0; ++gr) *sq++ = slst[*gr];
		*sq = 0;
		aggregate_sb(sd, sbuf, gbuf, wbuf);
	    }
#endif // USE_WEIGHT
	} else {
	    aggregate_sb(sd, slst, glst, wlst);
	}
	sd->exg_seq(sd->inex.exgl, sd->inex.exgr);
	return (hh);
}

SKL* Prrn::divideseq(Randiv* rdiv, int* lst[], GapsList* gsubs[], 
	VTYPE* pwt, MCTYPE* bch)
{
	LRAND	mark = rdiv->nextrandiv();
	MCTYPE	mcran = rdiv->mcr->mcrand_now();
	if (bch) *bch = mcran;

	rdiv->bin2lst2(lst, mark);
#if USE_WEIGHT
	*pwt = calcfact(wlst, srl->ktree->lead + mcran);
#endif
	if (totalNoSeq(lst[0]) < totalNoSeq(lst[1]))
	    std::swap(lst[0], lst[1]);
	GAPS*	gp[2] =
	    {gather(1, lst[0], gsubs[0]), gather(2, lst[1], gsubs[1])};
	SKL*	skl = gap2skl(gp[0], gp[1]);
	return (skl);
}

SKL* Prrn::divideseq(MCTYPE mcran, int* lst[], GapsList* gsubs[], VTYPE* pwt)
{
#if USE_WEIGHT
	*pwt = calcfact(wlst, srl->ktree->lead + mcran);
#endif
	if (totalNoSeq(lst[0]) < totalNoSeq(lst[1]))
	    std::swap(lst[0], lst[1]);
	GAPS*	gp[2] =
	    {gather(1, lst[0], gsubs[0]), gather(2, lst[1], gsubs[1])};
	SKL*	skl = gap2skl(gp[0], gp[1]);
	return (skl);
}

VTYPE Prrn::onecycle(SKL** skl, VTYPE pwt)
{
	if (!skl[0]) {
	    skl[1] = 0;
	    return (NEVSEL);
	}
	PwdM	pwd(seqs + 1, alnprm);
	if (pwd.swp) swapskl(skl[0]);
	PreSpScore	pss(seqs + 1, &pwd);
	VTYPE	sps = pss.calcSpScore(skl[0]);
	Gsinfo	alninf;
	VTYPE	scr;
//	seqs[1]->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
//	seqs[2]->exg_seq(algmode.lcl & 4, algmode.lcl & 8);
	skl[1] = align2(seqs + 1, &pwd, &scr, &alninf);
	if (!skl[1]) {
	    if (pwd.swp) std::swap(seqs[1], seqs[2]);
	    return (NEVSEL);
	}
	sps = sameskl(skl[0], skl[1])? 0: pwt * (alninf.fstat.val- sps);
	if (pwd.swp) {
	    if (alninf.sigII && alninf.sigII->lst)
		alninf.sigII->swaplst(seqs[1]->many, seqs[2]->many);
	    swapskl(skl[1]);
	    std::swap(seqs[1], seqs[2]);
	}
	return (sps);
}

VTYPE Prrn::onecycle(SKL** skl, int** lb, Randiv* rdiv)
{
	VTYPE	pwt;
	skl[0] = divideseq(rdiv, lb, glists + 1, &pwt);
	if (!skl[0]) {skl[1] = 0; return (NEVSEL);}
	return (onecycle(skl, pwt));
}

VTYPE Prrn::onecycle(SKL** skl, int** lb, int k)
{
	if (k == best_k || !lt(0, thargs[k].scr)) {
	    skl[0] = skl[1] = 0;
	    return (0);
	}
	lb[0] = thargs[k].lst[0];
	lb[1] = thargs[k].lst[1];
	VTYPE	pwt;
	skl[0] = divideseq(thargs[k].bch, lb, glists + 1, &pwt);
	if (!skl[0]) {skl[1] = 0; return (NEVSEL);}
	return (onecycle(skl, pwt));
}
 
static void* thread_onecycle(void* arg)
{
	ThreadArg*	targ = (ThreadArg*) arg;
	if (!targ->skl[0]) {
	    targ->skl[1] = 0;
	    targ->scr = NEVSEL;
	    return (void*) 0;
	}
	PwdM	pwd(targ->sqs, targ->alp);
	if (pwd.swp) swapskl(targ->skl[0]);
	PreSpScore	pss(targ->sqs, &pwd);
	VTYPE sps = pss.calcSpScore(targ->skl[0]);
	Gsinfo	alninf;
	VTYPE	scr;
	targ->skl[1] = align2(targ->sqs, &pwd, &scr, &alninf);
	if (!targ->skl[1]) {
	    targ->scr = NEVSEL;
	    return (void*) 0;
	}
	targ->scr = targ->pwt * (scr - sps);
	if (pwd.swp) {
	    if (alninf.sigII && alninf.sigII->lst)
		alninf.sigII->swaplst(targ->sqs[0]->many, targ->sqs[1]->many);
	    swapskl(targ->skl[1]);
	    std::swap(targ->sqs[0], targ->sqs[1]);
	}
	return (void*) 0;
}

VTYPE Prrn::best_of_n(SKL** skl, int** lstbuf, Randiv* rdiv)
{
	for (int i = 0, j = 0; j < num; ++j) {
	    thargs[i].skl[0] = divideseq(rdiv, thargs[i].lst, thargs[i].gsub, 
		&thargs[i].pwt, &thargs[i].bch);
	    if (thargs[i].skl[0]) {
		std::swap(seqs[1], thargs[i].sqs[0]);
		std::swap(seqs[2], thargs[i].sqs[1]);
		if (++i == no_thread) break;
	    }
	}

	for (int i = 0; i < no_thread; ++i) {
	    pthread_create(&handle[i], NULL, thread_onecycle,
		(void*) &thargs[i]);
	}
	for (int i = 0; i < no_thread; ++i) {
	    pthread_join(handle[i], NULL);
	}
	VTYPE	sps = NEVSEL;
	best_k = 0;
	for (int i = 0; i < no_thread; ++i) {
	    if (thargs[i].scr > sps)
		sps = thargs[best_k = i].scr;
	}
	for (int i = 0; i < 2; ++i) {
	    skl[i] = thargs[best_k].skl[i];
	    lstbuf[i] = thargs[best_k].lst[i];
	    thargs[best_k].gsub[i]->copy_to(glists[i + 1]);
	}
	for (int i = 0; i < no_thread; ++i) {
	    if (i != best_k) {
		delete[] thargs[i].skl[0];
		delete[] thargs[i].skl[1];
	    }
	}
	return (sps);
}

VTYPE Prrn::rir(VTYPE sps)
{
	SKL*	skl[2];
	int*	lstbuf[2] = {thlsts[0], thlsts[1]};
	INT	nrep = 0;
	Randiv	rdiv(seqs[0], srl, dm, randseed);

	INT	i = 0;
	int	k = 0;
	INT	maxi = maxitr * rdiv.cycle;
	for ( ; i < maxi; ++i) {
	    VTYPE	delta_scr = 0;
	    if (!thargs) {
		delta_scr = onecycle(skl, lstbuf, &rdiv);
	    } else if (!k) {
		delta_scr = best_of_n(skl, lstbuf, &rdiv);
		k = no_thread;
	    } else if (--k != best_k) {
		delta_scr = onecycle(skl, lstbuf, k);
	    } else	continue;
	    if (lt(0, delta_scr)) {	// improved
		if (synthgap(glists, skl[1], lstbuf)) {
		    nrep = 1;
		    sps += delta_scr;
		} else {
		    fatal("Alignment error !\n");
		}
	    } else	++nrep;
	    delete[] skl[0]; delete[] skl[1];
	    if (nrep >= rdiv.cycle) break;
	}
	countaln += i;
	return (sps);
}

Prrn::~Prrn()
{
	if (thargs) {
	    clearseq(thseqs, 3 * no_thread);
	    delete[] thseqs;
	    delete[] thskls;
	    ThreadArg*	wrk = thargs;
	    for (int i = 0; i < no_thread; ++i, ++wrk) {
		delete wrk->gsub[0];
		delete wrk->gsub[1];
	    }
	    delete[] thargs;
	    delete[] handle;
	}
	if (thlsts) {
	    delete[] thlbuf;
	    delete[] thlsts;
	}
}

Prrn::Prrn(mSeq** sqs, ALPRM* alp, Ssrel* trl, VTYPE& ref, bool sglthrd, DivMode dm_) 
	: seqs(sqs), msd(seqs[0]), alnprm(alp), srl(trl), num(trl->ss->num), dm(dm_),
	handle(0), thargs(0), thseqs(0), thskls(0), thpwds(0), thlsts(0), thlbuf(0),
	slst(0), sbuf(0), wlst(0), wbuf(0)
{
#if MONIT
	long	start = time(0);
#endif
	if (num <= 1) return;
	int	kk = num + 1;

	no_thread = (sglthrd && num >= paral_num)? min(paral_num, thread_num): 0;
	glst = new GapsList(num); 
#if USE_WEIGHT
	wlst = new VTYPE[kk];
	wbuf = new VTYPE[kk];
#endif
	slst = srl->takeapart(glst, msd, true);
	int	maxgaps = msd->len;
	sbuf = new mSeq*[kk];
	for (int i = 0; i <  srl->ss->num; ++i) {
	    mSeq*&	sd = slst[i];
#if SSHP
	    sd->inex.sshp = sshpprm? 1: 0;
#endif
	    sd->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
	    sd->convseq(RAWSEQ);
	    if (sd->len < maxgaps) maxgaps = sd->len;
	}
	maxgaps = std::max(msd->len - maxgaps + 3, MaxNoGaps);
	glists[0] = glst->resize(maxgaps);
	for (int i = 1; i < NGLIST; ++i) 
	    glists[i] = new GapsList(num, maxgaps, true);
	if (no_thread) {
	    handle = new pthread_t[no_thread];
	    thargs = new ThreadArg[no_thread];
	    ThreadArg*	wrk = thargs;
	    ThreadArg*	prv = thargs;
	    prv->sqs = thseqs = new mSeq*[3 * no_thread];
	    prv->skl = thskls = new SKL*[2 * no_thread];
	    prv->lst = thlsts = new int*[2 * no_thread];
	    prv->lst[0] = thlbuf = new int[2 * no_thread * kk];
	    prv->alp = alnprm;

	    for (int i = 0; i < no_thread; ++i, ++wrk) {
		wrk->sqs = prv->sqs + (i? 3: 0);
		initseq(wrk->sqs, 2); wrk->sqs[2] = 0;
		wrk->skl = prv->skl + (i? 2: 0);
		wrk->lst = prv->lst + (i? 2: 0);
		wrk->lst[0] = prv->lst[0] + (i? 2 * kk: 0);
		wrk->lst[1] = wrk->lst[0] + kk;
		wrk->gsub[0] = new GapsList(num, maxgaps, true);
		wrk->gsub[1] = new GapsList(num, maxgaps, true);
		wrk->scr = 0;
		wrk->alp = prv->alp;
		if (i) ++prv;
	    }
	} else {
	    thlsts = new int*[2];
	    thlsts[0] = thlbuf = new int[2 * kk];
	    thlsts[1] = thlsts[0] + kk;
	    thargs = 0; thseqs = 0; thskls = 0;
	}

#if USE_WEIGHT
	VTYPE	initsp = srl->pairsum_ss(msd, true);
#else
	VTYPE	initsp = srl->pairsum_ss(msd, false);
#endif
	countaln = 0;
	VTYPE	sps = rir(initsp);
	for (int i = 0; i <  srl->ss->num; ++i)
	    slst[i]->elim_column(DEL_GAP);
	gather(2);	// result is stored in seqs[2]
	clearseq(slst, srl->ss->num);
	for (int i = 0; i < NGLIST; ++i) delete glists[i];
	delete[] slst; delete[] sbuf;
#if USE_WEIGHT
	delete[] wlst; delete[] wbuf;
#endif

	sortseq(seqs[1], seqs[2], srl->ss, 0);
	replseq(seqs);
	ref += sps - initsp;
#if MONIT
	if (getprompt()) {
	  long	stop = time(0);
	  msd->fphseq(2, stderr);
	  prompt("  %8.1lf <-- %8.1lf, %2d grp, %4d rep, %2ld sec\n", 
	    (double) sps, (double) initsp, srl->ss->num, 
	    countaln, stop - start);
	}
#endif
}

static	const char lclwarn[] = "Error: -L option is not feasible to this msa !\n";

bool	IterMsa::preprrn(bool sglthrd)
{
	if (ss->num <= 2) return (false);
	int	oldlen = msd->len;
	int	exgl = msd->inex.exgl;
	int	exgr = msd->inex.exgr;
	bool	progmessg = getprompt();
	phyl_pwt();
	VTYPE	refined = 0;
	if (local_thr) {
	    Ssrel*	rrl = consregss(infc);
	    Ssrel*	srl = rrl? rrl: this;
	    if (srl->ss->num < 2) {
		delete rrl;
		return (false);
	    }
	    RANGE	svrng;
	    msd->saverange(&svrng);
	    RANGE*	atkrng = srl->consreg(msd, DISSIM);
	    if (progmessg) fprintrng(stderr, atkrng);
	    INEX	inex = msd->inex;
	    RANGE* rng = lastrng(atkrng);
	    for ( ; rng > atkrng; --rng) {
		msd->restrange(rng);
		Ssrel*	trl = srl->consregss(infc);
		Ssrel*	url = trl? trl: srl;
#if USE_WEIGHT
		if (!url->pairwt) url->recalcpw();
#endif
		msd->inex.exgl = (svrng.left == rng->left)? exgl: 0;
		msd->inex.exgr = (svrng.right == rng->right)? exgr: 0;
		Prrn	prrn(infc, alnprm, url, refined, sglthrd, dm);
		delete trl;
	    }
	    delete rrl;
	    msd->restrange(&svrng);
	    msd->right += msd->len - oldlen;
	    msd->inex = inex;
	    delete[] atkrng;
	} else {
	    Prrn	prrn(infc, alnprm, this, refined, sglthrd, dm);
	}
	if (progmessg) {
	    FTYPE	tpwt = 1.;
#if USE_WEIGHT
	    tpwt = sumofpwt();
#endif
	    VTYPE       newsp = pairsum_ss(msd, true);
	    msd->fphseq(2, stderr);
	    prompt("  %8.1lf <-- %8.1lf, %2d grp, %8.2lf %8.4lf\n",
	    (double) newsp, (double) (newsp - refined),
	    ss->num, newsp / tpwt,  tpwt);
	}
	return (gt(refined, 0));
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

template<>
mSeq* ProgMsa<Slnode>::getseq(Slnode* node)
{
	mSeq*	msd = 0;
	if (node->tid & artp_bit) {
	    msd = seqs[node->tid & bit_mask];
	    seqs[node->tid & bit_mask] = 0;
	} else {
	    msd = new mSeq;
	    msd->read_dbseq(wdbf, verdid[node->tid]);
	}
	if (ss) {
	    leaves[ss->num] = msd;
	    ss->group[ss->num++] = pl;
	    for (int i = 0; i < msd->many; ++i)
		*pl++ = ss->elms++;
	    *pl++ = EOTAB;
	    ss->group[ss->num] = 0;
	}
	return (msd);
}

IterMsa::IterMsa(mSeq* sd, mSeq** sqs, int nn, ALPRM* alp, int* lstodr, DivMode dm_)
	: Ssrel(), msd(infc[0]), seqs(sqs), no_seqs(nn), 
	  rsv(sd), alnprm(alp), lst_odr(lstodr), dm(dm_)
{
	initseq(infc, 3); infc[3] = 0;
	std::swap(msd, rsv);
	if (nn) {
	    ss = new Subset(nn, (const char*) 0);
	    int**	grp = ss->group;
	    ss->elms = 0;
	    for (int i = 0; i < nn; ++i)
		ss->elms += seqs[i]->many;
	    delete[] ss->pool;
	    int*	pi = ss->pool = new int[ss->elms + ss->num];
	    for (int i = 0, k = 0; i < nn; ++i) {
		*grp++ = pi;
		int	l = lst_odr? lst_odr[i]: i;
		for (int j = 0; j < seqs[l]->many; ++j) *pi++ = k++;
		*pi++ = EOTAB;
	    }
	    *grp = 0;
	}
}

IterMsa::IterMsa(mSeq* sd, mSeq** sqs, Subset* ss, ALPRM* alp, DivMode dm_)
	: Ssrel(ss), msd(infc[0]), seqs(sqs), no_seqs(ss? ss->num: 0), 
	  rsv(sd), alnprm(alp), lst_odr(0), dm(dm_)
{
	initseq(infc, 3); infc[3] = 0;
	std::swap(msd, rsv);
}

mSeq*	IterMsa::msa(bool sglthrd)
{
	if (checkss() == ERROR)
	    fatal("Error -- Improper grouping!");
	for (int rc = 0; rc < recycle; ++rc)
	    if (!preprrn(sglthrd)) break;
	std::swap(rsv, msd);
	return (rsv);
}

mSeq* SlfPrrn::make_msa(Slnode* node, bool sglthrd)
{
	Subset*	ss = node->ndesc > node->nunit? 
		new Subset(node->nunit, node->ndesc) : 0;
	ProgMsa<Slnode>	promsa(wdbf, verdid, seqs, ss, &alprm);
	mSeq*	msd = promsa.prog_up(node);
	if (!recycle) return (msd);
	IterMsa	itrmsa(msd, promsa.get_leaves(), ss, &alprm, divmode);
	msd = itrmsa.msa(sglthrd);
	msd->left = 0; msd->right = msd->len;
	return (msd);
}

static mSeq* makemsa(mSeq** seqs, int nn, ALPRM* alp, bool sglthrd)
{
	if (nn == 2) {
	    mSeq*	msd =  align2(seqs[0], seqs[1]);
#if USE_WEIGHT
	    if (msd) calcweight((Seq*) msd);
#endif // USE_WEIGHT
	    return (msd);
	}
	for (int i = 0; i < nn; ++i) seqs[i]->sid = i;
	DistCal	realin = distPrm.contbt_spj > 0? DynAln: (DistCal) algmode.any;
	DistMat	dm(seqs, nn, realin);
	Ktree	dt(&dm);
	ProgMsa<Knode>	promsa(seqs, nn, alp);
	mSeq*	msd = promsa.prog_up(dt.root);
	if (!recycle) return (msd);
	IterMsa	itrmsa(msd, seqs, nn, alp, promsa.lstodr(), divmode);
	msd = itrmsa.msa(sglthrd);
	msd->left = 0; msd->right = msd->len;
	return (msd);
}

static int getseq_main(CalcServer<mSeq>* svr, mSeq** sqs, ThQueue<mSeq>* q)
{
	mSeq**	seqs = (mSeq**) svr->prm;
	int	sid = sqs[0]->sid;
	std::swap(sqs[0], seqs[sid]);
	return (OK);
}

static mSeq* makemsa(int argc, const char** argv, int num)
{
	mSeq**	seqs = new mSeq*[num];
	initseq(seqs, num);
	AlnServer<mSeq> svr(argc, argv, IM_SNGL, catalog, (void*) seqs, &getseq_main);
	svr.autocomp(false);
	mSeq*	msd = makemsa(seqs, num, &alprm, true);
	clearseq(seqs, num);
	delete[] seqs;
	return (msd);
}

static mSeq* cut_in(mSeq* mom, mSeq* dau)
{
#if USE_WEIGHT
	if (!mom->weight) calcweight((Seq*) mom);
#endif // USE_WEIGHT
	mSeq*	msd = align2(mom, dau, 0);
	if (!msd) return (0);
#if USE_WEIGHT
	msd->weight = new VTYPE[mom->many + dau->many];
	vcopy(msd->weight, mom->weight, mom->many);
	for (int i = mom->many; i < msd->many; ++i) msd->weight[i] = 1;
#endif // USE_WEIGHT
	delete mom;
	return msd;
}

static	int sscmpf(const mSeq** a, const mSeq** b) {
	return ((*a)->right - (*a)->left - (*b)->right + (*b)->left);
}

void SlfPrrn::findap(Slnode* node, std::vector<Slnode*>& subtrees)
{
	if (node->nunit >= max_memb) {
	    node->left->root = node->right->root = node;
	    if (node->left->nunit < max_memb && node->right->nunit < max_memb) {
		if (!node->root) return;
		node->tid = (no_seqs + subtrees.size()) | artp_bit;
		subtrees.push_back(node);
	    } else {
		findap(node->left, subtrees);
		findap(node->right, subtrees);
	    }
	}
}

void SlfPrrn::move_to_outsiders(Slnode*& node)
{
	int	nn = node->ndesc;
	if (nn >= min_ndesc) return;
	if (node->isleaf()) {
	    outsiders.push_back(verdid[node->tid]);
	    ++no_outsider;
	    node = 0;
	    return;
	}
	FTYPE   hthr = alprm.thr * (1. - 1. / float(nn));
	if (node->dist > hthr) {
	    move_to_outsiders(node->left);
	    move_to_outsiders(node->right);
	    if (node->left && node->right) {
		trees->push_back(node->right);
		node = node->left;
	    } else if (node->left) node = node->left;
	    else	node = node->right;
	}
}

void SlfPrrn::prune_subtrees()
{
	if (!no_trees) return;
	for (TreePtrItr st = trees->begin(); st < trees->end(); ++st) {
	    Slnode*	node = *st;
	    move_to_outsiders(node);
	    if (!node) {
		trees->erase(st);
		--no_trees;
	    } else	*st = node;
	}
}

int SlfPrrn::restruct(Slnode* node)
{
	int	redc = 0;
	if (node->nunit >= max_memb) {
	    if (node->left->nunit < max_memb && node->right->nunit < max_memb) {
		if (!node->root) return (redc);
		redc = node->nunit - 1;
		node->left = node->right = 0;
	    } else
		redc = restruct(node->left) + restruct(node->right);
	    node->nunit -= redc;
	}
	return (redc);
}

static void print_msa(mSeq* sd)
{
	char	str[MAXL];
	sprintf(str, "%s.%d", header, ++msa_id);
	FILE*	fd = fopen(str, "w");
	if (!fd) fatal("Can't write to %s !\n", str);
	sd->typeseq(fd);
	fclose(fd);
}

//	articulation points

int SlfPrrn::make_artps(std::vector<Slnode*>& subtrees)
{
	int	no_atps = subtrees.size();
	if (!header && no_atps) {
	    mSeq**	atps = new mSeq*[no_seqs + no_atps];
	    if (seqs) {
		vcopy(atps, seqs, no_seqs);
		delete[] seqs;
	    }
	    seqs = atps;
	}
	Sltree* 	sltrees = new Sltree[no_atps];
	Sltree**	psltrees = new Sltree*[no_atps];
	Sltree* 	slt = sltrees;
	Sltree**	pslt = psltrees;
	int	tid = no_seqs;
	int	no_large = 0;
	int	nunitthr = no_atps < eff_paral_num? paral_num: INT_MAX;
	for (TreePtrItr st = subtrees.begin(); st < subtrees.end(); ++st, ++tid) {
	    if ((*st)->nunit >= nunitthr || (*st)->ndesc >= ndescthr) {
		mSeq*	sd = make_msa(*st, true);
//fprintf(stderr, "SlfPrrn: %d %d %d\n", (*st)->ndesc, (*st)->nunit, sd->len);
		if (!sd) fatal(fail_to_align);
		if (header) {
		    print_msa(sd);
		    delete sd;
		} else
		    seqs[tid] = sd;
		++no_large;
	    } else {
		slt->root = *st;
		slt->tid = tid;
		slt->vrtl = 0;
		slt->fget = fget;
		*pslt++ = slt++;
	    }
	}
	int	no_small = no_atps - no_large;
	runstat.stamp(no_large);
	if (no_small) {
	    CalcServer<Sltree> svr(IM_SNGL, (void*) this, 
		slfprrn_main, 0, 0, psltrees, no_small);
	    svr.autocomp();
	    runstat.stamp(no_small);
	}
	no_seqs = tid;
	delete[] sltrees;
	delete[] psltrees;
	return (no_atps);
}

mSeq* SlfPrrn::print_submsas()
{
	for (TreePtrItr st = trees->begin(); st < trees->end(); ++st) {
	    if ((*st)->nunit > 1) {
		mSeq*	sd = make_msa(*st, true);
		if (sd) print_msa(sd);
		delete sd;
	    }
	}
	return (0);
}

mSeq* SlfPrrn::make_msa(int molc)
{
	if (no_trees) {
	    std::vector<Slnode*>	subtrees;
	    while (true) {
		for (TreePtrItr st = trees->begin(); st < trees->end(); ++st)
		    findap(*st, subtrees);
		if (subtrees.size() == 0) break;
		no_seqs += make_artps(subtrees);
		for (TreePtrItr st = trees->begin(); st < trees->end(); ++st)
		    restruct(*st);
		subtrees.clear();
	    }
	    make_artps(*trees);
	}
	int	no_leaf = no_trees + no_outsider;
	mSeq**	tmp = (no_seqs < no_leaf)? new mSeq*[no_leaf]: 0;
	if (seqs) {
	    mSeq**	dst = tmp? tmp: seqs;
	    mSeq**	trm = seqs + no_seqs;
	    for (mSeq** src = seqs; src < trm; ++src) {
		if (*src) {
		    if (dst == src) ++dst;
		    else {*dst++ = *src; *src = 0;}
		}
	    }
	}
	if (tmp) {
	    delete[] seqs;
	    seqs = tmp;
	    initseq(seqs + no_trees, no_outsider);
	}
	no_seqs = no_leaf;

// multiple profiles alignmnet

	if (no_outsider) {	// read from db
	    for (int i = 0; i < no_outsider; ++i) {
		if (!seqs[i + no_trees]) seqs[i + no_trees] = new mSeq;
		seqs[i + no_trees]->read_dbseq(wdbf, outsiders[i]);
	    }
	}
	mSeq*	msd = 0;
	if (no_seqs == 1)
	    std::swap(msd, seqs[0]);
	else if (no_trees == 0) 
	    msd = makemsa(seqs, no_seqs, &alnprm, true);
	else {
	    if (no_trees == 1)
		std::swap(msd, seqs[0]);
	    else
		msd = makemsa(seqs, no_trees, &alnprm, true);
	    runstat.stamp(no_trees);
	    if (incl_outsider && msd && no_outsider) {
		if (no_outsider > 1)
		    qsort((UPTR) (seqs + no_trees), (INT) no_outsider, 
			sizeof(mSeq*), (CMPF) sscmpf);
		if (recycle) {
		    int	k = no_trees - 1;
		    seqs[k] = msd->aliaseq(seqs[k], true);
		    Subset	ss(no_outsider + 1, msd->many + no_outsider);
		    int*	pl = ss.pool;
		    int**	gr = ss.group;
		    *gr++ = pl;
		    int	i = 0;
		    while (i < msd->many) *pl++ = i++;
		    *pl++ = EOTAB;
		    for (int g = no_trees; g < no_seqs; ++g) {
			*gr++ = pl;
			*pl++ = i++;
			*pl++ = EOTAB;
			msd = cut_in(msd, seqs[g]);
		    }
		    *gr = 0;
		    IterMsa	itrmsa(msd, seqs + k, &ss, &alnprm, ONE_DIV);
		    msd = itrmsa.msa(true);
		} else {
		    for (int g = no_trees; g < no_seqs; ++g)
			msd = cut_in(msd, seqs[g]);
		}
	    }
	}
	clearseq(seqs, no_seqs);
	delete[] seqs;
	runstat.stamp(no_outsider);
	return (msd);
}

static void setdefparam()
{
	treemode.n_edge = 1;
	optimize(GLOBAL, MAXIMUM);
	setlsegs(1);
	setalgmode(0, 0);
	setdefPprm(DefPrrnPam, 2., 9., 0);
	setdefNprm(-2., 2., 4.);	// n, u, v
	setprmode(Row_Last, 'L', SILENT);
//	setprompt(false, false);
	alprm.sh = -60;		// band shoulder = 1/2 of shorter
	alprm.thr = 70;		// max distance to connect edges
	algmode.any = DynScr;	// alignment-score-based distance
	algmode.mns = 1;
	algmode.nsa = 1;
	OutPrm.SkipLongGap = 0;	// long gaps are reported
}

static	int slfprrn_main(CalcServer<Sltree>* svr, Sltree** nodes, ThQueue<Sltree>* q)
{
	SlfPrrn*	slp = (SlfPrrn*) svr->prm;
	mSeq*	sd = slp->make_msa((*nodes)->root, !q);
//fprintf(stderr, " Slfprrn_main: %d %d %d\n", (*nodes)->root->ndesc, (*nodes)->root->nunit, sd->len);
	if (header) {
	    m_thread_Lock(q);
	    print_msa(sd);
	    m_thread_Unlock(q);
	    delete sd;
	} else {
	    m_thread_Lock(q);
	    slp->seqs[(*nodes)->tid] = sd;
	    m_thread_Unlock(q);
	}
	return (OK);
}

//	construct MSA with slf-guide trees and then run prrn

static mSeq* de_novo_prrn(mSeq* before, int argc = 0, const char** argv = 0)
{
	int	molc = setdefmolc();
	ALGMODE	rsv_algmode = algmode;	//	reserve params;
	ALPRM	rsv_alprm = alprm;
	if (molc == PROTEIN) {
	    DefPrrnPam = getpam(0);
	    DefScndPam = getpam(1);
	} else	DefMmcPen = getsmn(4);

	sltree_defparam(molc);
	alprm.sh = dist_sh;

	prePwd(molc, true);
	Slforest*	slf = before?
	    new Slforest(before):
	    new Slforest(argc, argv, molc, catalog, min_seqs);
	runstat.stamp(slf->no_edges());
	resetSimmtxes();
	algmode = rsv_algmode;		//	restore pramas;
	alprm = rsv_alprm;
	if (molc == PROTEIN) {
	    setpam(DefPrrnPam, 0);
	    setpam(DefScndPam, 1);
	} else	setNpam(4, DefMmcPen);
	prePwd(molc);
	SlfPrrn	slfprrn(*slf);
	slfprrn.prune_subtrees();
	mSeq*	msd = header? slfprrn.print_submsas(): slfprrn.make_msa(molc);
	if (header) slf->print_singletons();
	delete slf;
	return (msd);
}

static void replacemem(mSeq* msd, int m, mSeq* upd)
{
	int*	which = new int[msd->many + 1];
	int	nn = msd->many - 1;
	which[nn] = which[msd->many] = -1;
#if USE_WEIGHT
	if (!msd->pairwt) {
	    Ktree	ktree(msd);
	    msd->pairwt = ktree.calcpw();
	}
	if (!msd->weight) msd->weight = new VTYPE[msd->many];
#endif // USE_WEIGHT
	for (int i = 0, k = 0; i < msd->many; ++i) {
	    if (i != m) which[k++] = i;
#if USE_WEIGHT
	    msd->weight[i] = (i == m)? 0: msd->pairwt[elem(i, m)];
#endif // USE_WEIGHT
	}
	mSeq*	nsd = msd->extseq(0, which, CPY_ALL + RLT_SBI);
	nsd->elim_column(DEL_GAP);
	mSeq*	asd = align2(nsd, upd);
	delete nsd;
	if (!asd) fatal(fail_to_align);
	which[m] = nn;
#if USE_WEIGHT
	FTYPE*	wt = msd->weight; msd->weight = 0;
	FTYPE*	pw = msd->pairwt; msd->pairwt = 0;
#endif // USE_WEIGHT
	for (int i = m; i < nn; ++i) which[i + 1] = i;
	msd = asd->extseq(msd, which, CPY_ALL + RLT_SBI);
#if USE_WEIGHT
	msd->weight = wt;
	msd->pairwt = pw;
#endif // USE_WEIGHT
	delete asd;
	delete[] which;
}

static mSeq* update_prrn(mSeq*& host, int argc = 0, const char** argv = 0)
{
	if (update) min_host_size = 2;
	int	no_hosts = 1;
	int	no_guests = 0;
	for (int c = 0; c < argc; ++c) {
	    Seq*	sd = new Seq(argv[c]);
	    if (!sd || !sd->many) continue;
	    if (sd->many < min_host_size) ++no_guests;
	    else	++no_hosts;
	    delete sd;
	}
	int	no_seqs = no_hosts + no_guests;
	mSeq**	seqs = new mSeq*[no_seqs];
	seqs[0] = host;
	int	h = 1;
	int	g = 0;
	int	hmany = host->many;
	int	maxmem = hmany;
	for (int c = 0; c < argc; ++c) {
	    mSeq* sd = new mSeq(argv[c]);
	    if (!sd || !sd->many) continue;
	    if (sd->many < min_host_size) seqs[no_hosts + g++] = sd;
	    else {
		seqs[h++] = sd;
		hmany += sd->many;
		if (sd->many > maxmem) maxmem = sd->many;
	    }
	}
	for (int h = 0; h < no_hosts; ++h) {
	    mSeq*	msd = seqs[h];
	    if (!msd->inex.algn) {
		if (!msd->test_aligned()) {
		    seqs[h] = de_novo_prrn(msd);
		    if (!seqs[h]) fatal(fail_to_align);
		    delete msd;
		    if (!h)	host = seqs[0];
		}
	    }
	}
	if (update && no_hosts && no_guests) {
	    StrHash<int>	hostmem(hmany);
	    for (int h = 0, k = 0; h < no_hosts; ++h, k += maxmem) {
		mSeq*&	sd = seqs[h];
		for (int i = 0; i < sd->many; ++i)
		    hostmem.assign((*sd->sname)[i], k + i + 1);
	    }
	    for (int g = 0; g < no_guests; ++g) {
		mSeq*&	sd = seqs[no_hosts + g];
		KVpair<INT, int>*	kv = hostmem.find((*sd->sname)[0]);
		if (kv) {
		    int	code = kv->val - 1;
		    replacemem(seqs[code / maxmem], code % maxmem, sd);
		    delete sd; sd = 0;
		}
	    }
	}
	runstat.stamp(no_seqs);
	mSeq*	msd = 0;
	if (no_hosts == 1) {
	    msd = host->aliaseq(0, true);
	} else {
	    msd = makemsa(seqs, no_hosts, &alprm, true);
	    clearseq(seqs + 1, no_hosts - 1);
	}
	runstat.stamp(no_hosts);
	for (int i = no_hosts; i < no_seqs; ++i)
	    if (seqs[i]) msd = cut_in(msd, seqs[i]);
	clearseq(seqs + no_hosts, no_guests);
	delete[] seqs;
	runstat.stamp(no_guests);
	if (update || no_hosts == 1) {
	    IterMsa	itrmsa(msd);
	    msd = itrmsa.msa(true);
	};
	return (msd);
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

void Msa::phylsort()
{
	if (!pairwt) calcpw(msd);
	lst_odr = new int[msd->many + 1];
	lst_idx = 0;
	lstodr(ktree->root);
	lst_odr[lst_idx] = -1;
	mSeq*	nsd = msd->extseq(0, lst_odr, CPY_ALL + RLT_SBI);
	delete[] lst_odr; lst_odr = 0;
	std::swap(msd, nsd);
	delete nsd;
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

//      lone eijunction
	if (msd->sigII) {
	    msd->sigII->mkeijtab(msd->many);
	    for (int i = 0; i < msd->many; ++i)
		outlier[i].eij = msd->sigII->lone[i];
	}

//      unusual amount of indels
	RANGE	tmp;
	msd->saverange(&tmp);
	FTYPE*	lclsod = 0;
	FTYPE*	glbsod = calcdistsum((Seq*) msd);
	flength = new float[msd->many];
	int*	order = new int[msd->many];
	int*	olist = new int[msd->many];
	std::swap(olr_thr, local_thr);
	RANGE*	atkrng = consreg(msd, DISSIM);
	RANGE*	rng = atkrng;
	RANGE*	lst = atkrng + atkrng->left - 1;
	int last = atkrng->left - 3;
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
		    else    ++outlier[-olist[k] - 1].del_f;
		} else if (i == last) {
		    if (olist[k] >= 0) ++outlier[olist[k]].ins_l;
		    else    ++outlier[-olist[k] - 1].del_l;
		} else {
		    if (olist[k] >= 0) ++outlier[olist[k]].ins_m;
		    else    ++outlier[-olist[k] - 1].del_m;
		}
	    }
	    float   x;
	    vavsd(x, flength, msd->many);
	    stddev[0] += 1;
	    stddev[1] += x * x;
	    if (last) {
	  	if (i == 0) stddev[2] = x;
		if (i == last) stddev[3] = x;
	    } else
		stddev[2] = stddev[3] = -1.;

//      unusually divergent seq
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
	std::swap(olr_thr, local_thr);
	delete[] atkrng;
	delete[] olist;
	delete[] order;
	delete[] flength; flength = 0;
	delete[] lclsod;
	delete[] glbsod;
	msd->restrange(&tmp);
	return outlier;
}

FTYPE Msa::nomal_pairwt()
{
	int	nn = ncomb(ss->num);
#if USE_WEIGHT
	if (pairwt) {
	    if (ss->num > thr_gfq_23)
		return (ktree->sum_of_pairwt());
	    FTYPE   f = 0;
	    for (int i = 0; i < nn; ++i)
		f += pairwt[i];
	    return (f);
	} else
#endif // USE_WEIGHT
	return ((FTYPE) nn);
}

mSeq* Msa::output()
{
	setup_output(algmode.nsa);
	if (OutPrm.sortodr == BY_TREE) phylsort();
	if (algmode.nsa & 1) {
	    if (OutPrm.out_file && !(out_fd = fopen(OutPrm.out_file, "w")))
		fatal("Can't write to %s !\n", OutPrm.out_file);
	    if (getlpw()) msd->typeseq();
	}
	if (algmode.nsa & 2) {
	    float stdiv[4] = {0, 0, 0, 0};
	    Outlier*	outlier = findoutliers(stdiv);
	    Outlier*	olr = outlier;
	    FILE*	fo = stdout;
	    int	acc[10];
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
		char	fmt[MAXL];
		sprintf(fmt, "%%5d %%-%ds\t", msd->sname->longest());
		strcat(fmt + strlen(fmt), "%3d %2d %d %d %d %d %d %d %d\n");
		if (!OutPrm.olrsum) fprintf(fo, fmt, 
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
	    delete[] outlier;
	    if (fo != stdout && fo != out_fd) fclose(fo);
	}
	int span = msd->right - msd->left;
	if (algmode.nsa & 4 && span > 0) {
#if USE_WEIGHT
	    FTYPE*   nopwt = 0;
	    std::swap(nopwt, pairwt);
	    VTYPE   spscr = pairsum_ss(msd, false);
	    std::swap(nopwt, pairwt);
	    VTYPE   wspscr = pairsum_ss(msd);
#else
	    VTYPE   spscr = pairsum_ss(msd, false);
	    VTYPE   wspscr = spscr;
#endif	// USE_WEIGHT
	    fprintf(out_fd, "%s [ %d ] %d\t%7.1f %7.3f %7.1f %7.3f\n",
		msd->sqname(), msd->many, span,
		float(spscr), float(100. * spscr / ncomb(msd->many) / span),
		float(wspscr), float(100. * wspscr / nomal_pairwt() / span));
	}
	return (msd);
}

int main(int argc, const char** argv)
{
	progname = *argv;
	setdefparam();
	InputSeqTest	tv = zero_InputSeqTest;
	AlnServer<Seq>* ssvr = new AlnServer<Seq>(argc, argv, 
	    IM_SNGL, IM_MULT, (void*) &tv, &testInputSeq);
	catalog = ssvr->catalog;
	ssvr->autocomp(false);
	if (tv.bad) fatal("%d seqs were not found !\n", tv.bad);
	delete ssvr;
	if (!argc && !catalog && !guidetree) usage();
	if (paral_num < 0) paral_num = thread_num;
	set_localprm(local_u, local_v, local_thr);	// for consreg
	spb_fact();					// set spb factor

	runstat.stamp();
	mSeq*	host = 0;
	set_max_memb(max_memb);
	set_min_memb(min_memb);
	if (tv.molc) {
	    setdefmolc(tv.molc);
#if SSHP
	    if (tv.molc == PROTEIN) initSsHpPrm();
#endif
	    prePwd(tv.molc);
	}
	if (guidetree) {				// given guide tree
	    Btree<Tnode> ltree(guidetree);
	    if (ltree.members < 2) fatal("No or Too small guide tree !\n");
	    if (!tv.molc) {
		Seq     sd(1);
		sd.getseq(ltree.mname[0]);
		tv.molc = sd.inex.molc;
		setdefmolc(tv.molc);
#if SSHP
		if (tv.molc == PROTEIN) initSsHpPrm();
#endif
		prePwd(tv.molc);
	    }
	    ProgMsa<Tnode>	promsa(ltree, guidetree);
	    mSeq* prgsd = promsa.prog_up(ltree.root);
	    host = update_prrn(prgsd);
	    delete prgsd;
	} else if (tv.many < 2) {			// sequence files
	    fatal("Too small number of input sequences: %d !\n", tv.many);
	} else if (tv.maxmany >= min_host_size) {
	    host = new mSeq(*argv);
	}
	mSeq*	msd = 0;
	if (host) {				// aligned
	    if (!guidetree) {--argc; ++argv;}
	    msd = update_prrn(host, argc, argv);
	} else if (tv.num >= min_seqs)		// sl-forest
	    msd = de_novo_prrn(host, argc, argv);
	else 					// all-by-all
	    msd = makemsa(argc, argv, tv.many);
	delete host;
	if (msd) {
	    Msa	msa(msd);
	    msd = msa.output();
	    delete msd;
	}
#if SSHP
	if (tv.molc == PROTEIN) eraseSsHpPrm();
#endif
	runstat.conclude();
	eraStrPhrases();
	EraDbsDt();
	resetSimmtxes();
	return (0);
}

