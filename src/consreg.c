/*****************************************************************************
*
*	Locate significantly conserved or significantly divergent
*	regions in a multiple sequence alignment
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
*	Osamu Gotoh, Ph.D.	(2003-)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@i.kyoto-u.ac.jp>>
*****************************************************************************/

#include "aln.h"
#include "mseq.h"
#include "maln.h"
#include "css.h"
#include "mgaps.h"
#include "phyl.h"
#include "consreg.h"
#include "randiv.h"
#include "gfreq.h"

#ifndef MAIN
#define MAIN	0
#endif

static	const	int	RECS = 4;
static	const	int	lwrlmt = 3;
static	ALPRM	localprm = 
	{3., 10, 0., 0., 0, 1.0, 35., 1., 8., 0.5, 7, 2, 100, 1};


static	RANGE*	consself(mSeq* sd);

static	int	monit = 0;

void set_localprm(float u, float v, float t)
{
	if (u >= 0.) localprm.u = u;
	if (v >= 0.) localprm.v = v;
	if (t >= 0.) localprm.thr = t;
}

Conserved1::Conserved1(mSeq* sd_) : sd(sd_)
{
	msp = new Msap(sd, 0, &localprm);
	vthr = (VTYPE) (localprm.thr * sd->sumwt);
	fthr = (VTYPE) -((2 * localprm.u + localprm.v) * sd->sumwt);
	length = sd->right - sd->left;
	results = (monit & 4)? new float[length * RECS]: 0;
}

Conserved1::~Conserved1()
{
	delete msp;
	delete[] results;
}

Conserved2::Conserved2(mSeq** sqs) : seqs(sqs)
{
	pwd = new PwdM(seqs, &localprm);
	a = seqs[0]; b = seqs[1];
	fthr = (VTYPE) -((2 * localprm.u + localprm.v)) * a->sumwt * b->sumwt;
	length = a->right - a->left;
	results = (monit & 4)? new float[length * RECS]: 0;
}

Conserved2::~Conserved2()
{
	if (pwd->swp) swapseq(seqs, seqs + 1);
	delete	pwd;
	delete[] results;
}

RANGE* Conserved1::cons1_nv()
{
	VTYPE	scr = 0;
	VTYPE	mxv = 0;
	VTYPE	sum = 0;
	RANGE*	cng = new RANGE[sd->right - sd->left];
	RANGE*	wng = cng;
	int*	gla = new int[sd->many];

	sd->pregap(gla);
	mSeqItr	asi(sd, sd->left);
	float*	fwk = results;
	RANGE	rng = {0, 0};
	*wng++ = rng;

	for (int i = sd->left; i < sd->right; ++i, ++asi) {
	    VTYPE	s = (msp->*msp->sim1)(asi) + (msp->*msp->crg1)(gla, asi);
	    incrgap(gla, asi.res, sd->many);
	    sum += s;
	    if (scr == 0 && s > 0) rng.left = i;
	    scr += s;
	    if (scr < 0) scr = 0;
	    else if (scr >= vthr && scr > mxv) {
		mxv = scr;
		rng.right = i + 1;
	    }
	    if (mxv > 0 && (scr <= 0 || scr < mxv - vthr)) {
		*wng++ = rng;
		mxv = scr = sum = 0;
	    }
	    if (monit & 2) 
		printf("%3d %6.1f %6.1f %6.1f %6.1f\n", i, 
		(double) s, (double) sum, (double) scr, (double) mxv);
	    if (monit & 4) {
		*fwk++ = s; *fwk++ = sum; *fwk++ = scr; *fwk++ = mxv;
	    }
	}
	if (scr >= vthr) *wng++ = rng;
	delete[] gla;
	*wng++ = endrng;
	cng->left = cng->right = wng - cng;
	return (cng);
}

RANGE* Conserved1::cons1_pf()
{
	VTYPE	scr = 0;
	VTYPE	mxv = 0;
	VTYPE	sum = 0;
	RANGE*	cng = new RANGE[sd->right - sd->left];
	RANGE*	wng = cng;
	float*	fwk = results;

	mSeqItr	asi(sd, sd->left);
	RANGE	rng = {0, 0};
	*wng++ = rng;

	for (int i = sd->left; i < sd->right; ++i, ++asi) {
	    VTYPE	s = (msp->*msp->sim1)(asi) + 
		msp->newgap3(*asi.sfq, *asi.tfq);
	    sum += s;
	    if (scr == 0 && s > 0) rng.left = i;
	    scr += s;
	    if (scr < 0) scr = 0;
	    else if (scr >= vthr && scr > mxv) {
		mxv = scr;
		rng.right = i + 1;
	    }
	    if (mxv > 0 && (scr <= 0 || scr < mxv - vthr)) {
		*wng++ = rng;
		mxv = scr = sum = 0;
	    }
	    if (monit & 2) 
		printf("%3d %6.1f %6.1f %6.1f %6.1f\n", i, 
		(double) s, (double) sum, (double) scr, (double) mxv);
	    if (monit & 4) {
		*fwk++ = s; *fwk++ = sum; *fwk++ = scr; *fwk++ = mxv;
	    }
	}
	if (scr >= vthr) *wng++ = rng;
	*wng++ = endrng;
	cng->left = cng->right = wng - cng;
	return (cng);
}

RANGE* Conserved2::cons2_ng()
{
	VTYPE	scr = 0;
	VTYPE	mxv = 0;
	VTYPE	sum = 0;
	RANGE*	cng = new RANGE[a->right - a->left];
	RANGE*	wng = cng;

	mSeqItr	asi(a, a->left);
	mSeqItr	bsi(b, b->left);
	float*	fwk = results;
	RANGE	rng = {0, 0};
	*wng++ = rng;

	for (int i = a->left; i < a->right; ++i, ++asi, ++bsi) {
	    VTYPE	s = (pwd->*pwd->sim2)(asi, bsi);
	    sum += s;
	    if (scr == 0 && s > 0) rng.left = i;
	    scr += s;
	    if (scr < 0) scr = 0;
	    else if (scr >= pwd->Vthr && scr > mxv) {
		mxv = scr;
		rng.right = i + 1;
	    }
	    if (mxv > 0 && (scr <= 0 || scr < mxv - pwd->Vthr)) {
		*wng++ = rng;
		mxv = scr = sum = 0;
	    }
	    if (monit & 2) 
		printf("%3d %6.1f %6.1f %6.1f %6.1f\n", i, 
		(double) s, (double) sum, (double) scr, (double) mxv);
	    if (monit & 4) {
		*fwk++ = s; *fwk++ = sum; *fwk++ = scr; *fwk++ = mxv;
	    }
	}
	if (scr >= pwd->Vthr) *wng++ = rng;
	*wng++ = endrng;
	cng->left = cng->right = wng - cng;
	return (cng);
}

bool Conserved2::may_fuse_ng()
{
	VTYPE	scr = 0;

	mSeqItr	asi(a, a->left);
	mSeqItr	bsi(b, b->left);
	mSeqItr	tsi(a, a->right);

	for ( ; asi < tsi; ++asi, ++bsi) {
	    scr += (pwd->*pwd->sim2)(asi, bsi);
	    if (scr > 0) scr = 0;
	    if (scr < fthr) break;
	}
	return (asi == tsi);
}

RANGE* Conserved2::cons2_nv()
{
	VTYPE	scr = 0;
	VTYPE	mxv = 0;
	VTYPE	sum = 0;
	RANGE*	cng = new RANGE[a->right - a->left];
	RANGE*	wng = cng;
	int*	gla = new int[a->many];
	int*	glb = new int[b->many];

	a->pregap(gla);
	b->pregap(glb);
	mSeqItr	asi(a, a->left);
	mSeqItr	bsi(b, b->left);
	float*	fwk = results;
	RANGE	rng = {0, 0};
	*wng++ = rng;

	for (int i = a->left; i < a->right; ++i, ++asi, ++bsi) {
	    VTYPE	s = (pwd->*pwd->sim2)(asi, bsi) 
			+ (pwd->*pwd->crg2)(gla, glb, asi, bsi, 0);
	    incrgap(gla, asi.res, a->many);
	    incrgap(glb, bsi.res, b->many);
	    sum += s;
	    if (scr == 0 && s > 0) rng.left = i;
	    scr += s;
	    if (scr < 0) scr = 0;
	    else if (scr >= pwd->Vthr && scr > mxv) {
		mxv = scr;
		rng.right = i + 1;
	    }
	    if (mxv > 0 && (scr <= 0 || scr < mxv - pwd->Vthr)) {
		*wng++ = rng;
		mxv = scr = sum = 0;
	    }
	    if (monit & 2) 
		printf("%3d %6.1f %6.1f %6.1f %6.1f\n", i, 
		(double) s, (double) sum, (double) scr, (double) mxv);
	    if (monit & 4) {
		*fwk++ = s; *fwk++ = sum; *fwk++ = scr; *fwk++ = mxv;
	    }
	}
	if (scr >= pwd->Vthr) *wng++ = rng;
	delete[] gla; delete[] glb;
	*wng++ = endrng;
	cng->left = cng->right = wng - cng;
	return (cng);
}

bool Conserved2::may_fuse_nv()
{
	VTYPE	scr = 0;
	int*	gla = new int[a->many];
	int*	glb = new int[b->many];

	mSeqItr	asi(a, a->left);
	mSeqItr	bsi(b, b->left);
	mSeqItr	tsi(a, a->right);
	a->pregap(gla);
	b->pregap(glb);

	for ( ; asi < tsi; ++asi, ++bsi) {
	    scr += (pwd->*pwd->Sim2)(asi, bsi) 
		+ (pwd->*pwd->crg2)(gla, glb, asi, bsi, 0);
	    incrgap(gla, asi.res, a->many);
	    incrgap(glb, bsi.res, b->many);
	    if (scr > 0) scr = 0;
	    if (scr < fthr) break;
	}
	delete[] gla; delete[] glb;
	return (asi == tsi);
}

RANGE* Conserved2::cons2_hf()
{
	VTYPE	scr = 0;
	VTYPE	mxv = 0;
	VTYPE	sum = 0;
	RANGE*	cng = new RANGE[a->right - a->left];
	RANGE*	wng = cng;
	float*	fwk = results;
	mSeqItr	asi(a, a->left);
	mSeqItr	bsi(b, b->left);
	int*	glb = new int[b->many];
	RANGE	rng = {0, 0};
	*wng++ = rng;

	b->pregap(glb);
	for (int i = a->left; i < a->right; ++i, ++asi, ++bsi) {
	    VTYPE s = (pwd->*pwd->sim2)(asi, bsi) 
		+ pwd->newgap2(asi, bsi.res, glb);
	    incrgap(glb, bsi.res, b->many);
	    sum += s;
	    if (scr == 0 && s > 0) rng.left = i;
	    scr += s;
	    if (scr < 0) scr = 0;
	    else if (scr >= pwd->Vthr && scr > mxv) {
		mxv = scr;
		rng.right = i + 1;
	    }
	    if (mxv > 0 && (scr <= 0 || scr < mxv - pwd->Vthr)) {
		*wng++ = rng;
		mxv = scr = sum = 0;
	    }
	    if (monit & 2) 
		printf("%3d %6.1f %6.1f %6.1f %6.1f\n", i, 
		(double) s, (double) sum, (double) scr, (double) mxv);
	    if (monit & 4) {
		*fwk++ = s; *fwk++ = sum; *fwk++ = scr; *fwk++ = mxv;
	    }
	}
	if (scr >= pwd->Vthr) *wng++ = rng;
	delete[] glb;
	*wng++ = endrng;
	cng->left = cng->right = wng - cng;
	return (cng);
}

bool Conserved2::may_fuse_hf()
{
	VTYPE	scr = 0;

	mSeqItr	asi(a, a->left);
	mSeqItr	tsi(a, a->right);
	mSeqItr	bsi(b, b->left);
	int*	glb = new int[b->many];

	b->pregap(glb);
	for ( ; asi < tsi; ++asi, ++bsi) {
	    scr += (pwd->*pwd->sim2)(asi, bsi) 
		+ pwd->newgap2(asi, bsi.res, glb);
	    incrgap(glb, bsi.res, b->many);
	    if (scr > 0) scr = 0;
	    if (scr < fthr) break;
	}
	delete[] glb;
	return (asi == tsi);
}

RANGE* Conserved2::cons2_pf()
{
	VTYPE	scr = 0;
	VTYPE	mxv = 0;
	VTYPE	sum = 0;
	RANGE*	cng = new RANGE[a->right - a->left];
	RANGE*	wng = cng;
	float*	fwk = results;
	mSeqItr	asi(a, a->left);
	mSeqItr	bsi(b, b->left);
	RANGE	rng = {0, 0};
	*wng++ = rng;

	for (int i = a->left; i < a->right; ++i, ++asi, ++bsi) {
	    VTYPE	s = (pwd->*pwd->sim2)(asi, bsi)
		+ pwd->newgap3(*asi.sfq, *bsi.tfq)
		+ pwd->newgap3(*bsi.sfq, *asi.tfq);
	    sum += s;
	    if (scr == 0 && s > 0) rng.left = i;
	    scr += s;
	    if (scr < 0) scr = 0;
	    else if (scr >= pwd->Vthr && scr > mxv) {
		mxv = scr;
		rng.right = i + 1;
	    }
	    if (mxv > 0 && (scr <= 0 || scr < mxv - pwd->Vthr)) {
		*wng++ = rng;
		mxv = scr = sum = 0;
	    }
	    if (monit & 2) 
		printf("%3d %6.1f %6.1f %6.1f %6.1f\n", i, 
		(double) s, (double) sum, (double) scr, (double) mxv);
	    if (monit & 4) {
		*fwk++ = s; *fwk++ = sum; *fwk++ = scr; *fwk++ = mxv;
	    }
	}
	if (scr >= pwd->Vthr) *wng++ = rng;
	*wng++ = endrng;
	cng->left = cng->right = wng - cng;
	return (cng);
}

bool Conserved2::may_fuse_pf()
{
	VTYPE	scr = 0;

	mSeqItr	asi(a, a->left);
	mSeqItr	bsi(b, b->left);
	mSeqItr	tsi(a, a->right);

	for ( ; asi < tsi; ++asi, ++bsi) {
	    scr += (pwd->*pwd->sim2)(asi, bsi)
		+ pwd->newgap3(*asi.sfq, *bsi.tfq)
		+ pwd->newgap3(*bsi.sfq, *asi.tfq);
	    if (scr > 0) scr = 0;
	    if (scr < fthr) break;
	}
	return (asi == tsi);
}

static RANGE* consself(mSeq* sd)
{
	Conserved1 csv1(sd);

	switch (csv1.msp->calc_mode) {
	    case BAD_ALN: return 0;
	    case NTV_ALB: return csv1.cons1_nv();
	    case HLF_ALB: 
	    default:	return csv1.cons1_pf();
	}
}

RANGE* Ssrel::constwo(mSeq** seqs)
{
	Conserved2 csv2(seqs);

	switch (csv2.pwd->alnmode) {
	    case NGP_ALB: return csv2.cons2_ng();
	    case NTV_ALB: return csv2.cons2_nv();
	    case RHF_ALB:
	    case HLF_ALB: return csv2.cons2_hf();
	    case GPF_ALB:
	    default:	return csv2.cons2_pf();
	}
}

bool ConservedT::mayfuse(mSeq** seqs)
{
	Conserved2 csv2(seqs);

	switch (csv2.pwd->alnmode) {
	    case NGP_ALB: return csv2.may_fuse_ng();
	    case NTV_ALB: return csv2.may_fuse_nv();
	    case RHF_ALB:
	    case HLF_ALB: return csv2.may_fuse_hf();
	    case GPF_ALB:
	    default:	return csv2.may_fuse_pf();
	}
}

void Ssrel::dividseq(mSeq* sons[], mSeq* fath, Subset* ss,
	 Randiv& rdv, int* lstbuf[])
{
	LRAND	rnbr = rdv.nextrandiv();
	rdv.bin2lst2g(lstbuf, rnbr, ss);
	fath->extseq(sons[0], lstbuf[0], CPY_SEQ + CPY_NBR);
	fath->extseq(sons[1], lstbuf[1], CPY_SEQ + CPY_NBR);
}

void Ssrel::shiftrng(RANGE* rng, int by)
{
	if (by == 0) return;
	while (neorng(++rng)) {
	    rng->left += by;
	    rng->right += by;
	}
}

RANGE*	Ssrel::consreg(mSeq* sd, int dissim)
{
	if (!ktree) ktree = new Ktree(sd);
	RANGE*	united = 0;
	Randiv	rdv(sd, this, TREEDIV, 0);

	if (!ss) ss = new Subset(sd->many);
	if (ss->num != 4) monit &= ~4;
	mSeq*	seqs[2];
	initseq(seqs, 2);
	RANGE*	prv = singlerng(0, sd->left, sd->right);
	RANGE*	sdrng = dissim? copyrng(prv): 0;
	int*	lstbuf[2];
	lstbuf[0] = new int[ss->elms + 1];
	lstbuf[1] = new int[ss->elms + 1];
	for (INT i = 0; i < rdv.cycle; ++i, prv = united) {
	    dividseq(seqs, sd, ss, rdv, lstbuf);
	    RANGE*	rng = constwo(seqs);
	    shiftrng(rng, sd->left);
	    united = cmnrng(prv, rng, 0);
	    delete[] rng;
	    delete[] prv;
	    if (emptyrng(united)) break;
	}
	delete[] lstbuf[0];
	delete[] lstbuf[1];
	clearseq(seqs, 2);
	if (dissim) {
	    prv = complerng(sdrng, united);
	    delete[] united;
	    united = prv;
	}
	delete[] sdrng;
	return (united);
}

int ConservedT::grcopy(int i)
{
	*gg = pp;
	if (gg_par) {
	    int*	aa = gg_par[i];
	    while ((*pp++ = *aa++) >= 0) ;
	} else {
	    *pp++ = i;
	    *pp++ = -1;
	}
	return (gg++ - gg_son);
}

int* ConservedT::testcons(Knode* node)
{
	if (node->isleaf())
	    return (gg_son[grcopy(node->tid)]);

	int*	lst1 = testcons(node->left);
	int*	lst2 = testcons(node->right);
	if (lst1 && lst2) {
	    seqs[0]->extseq(seqs[1], lst1, CPY_SEQ);
	    seqs[0]->extseq(seqs[2], lst2, CPY_SEQ);
	    if (mayfuse(seqs + 1)) {
		while ((lst2[-1] = *lst2) >= 0) lst2++;
		gg--; pp--;
		return (lst1);
	    }
	}
	return (0);
}

ConservedT::ConservedT(mSeq** sqs, Subset* ss, Subset* tt, bool lf)
	: seqs(sqs)
{
	gg_par = ss? ss->group: 0;
	gg = gg_son = tt->group;
	pp = tt->pool;
	ntid = 0; lead = 0;
	if (lf) {
	    int	n = ss? ss->num: sqs[0]->many;
	    int	k = 2 * n - 1;
	    leaf = new int[k];
	    for (int i = 0; i < k; ++i) leaf[i] = -1;
	} else {
	    leaf = 0;
	}
}

int ConservedT::testssrel(Knode* node)
{
	if (node->isleaf())
	    return (leaf[node->tid] = grcopy(node->tid));

	int	lt = testssrel(node->left);
	int	rt = testssrel(node->right);
	if (lt >= 0 && rt >= 0) {
	    int*	lst1 = gg_son[lt];
	    int*	lst2 = gg_son[rt];
	    seqs[0]->extseq(seqs[1], lst1, CPY_SEQ);
	    seqs[0]->extseq(seqs[2], lst2, CPY_SEQ);
	    if (mayfuse(seqs + 1)) {
		while ((lst2[-1] = *lst2) >= 0) lst2++;
		--pp; --gg;
		leaf[node->tid] = lt;
		return (lt);
	    }
	}
	return (-1);
}

void Ssrel::clear()
{
	delete ktree; ktree = 0;
	delete[] spscr; spscr = 0;
	delete[] pairwt; pairwt = 0;
}

mSeq** Ssrel::takeapart(GapsList* glst, mSeq* sd, bool save)
{
	int	snl = save && glst? CPY_SEQ + CPY_SBI: CPY_SEQ;
	int	nn = ss? ss->num: sd->many;
	mSeq**	slist = 0;
	mSeq*	sq = 0;

	if (save) {
	    slist = new mSeq*[nn + 1];
	    initseq(slist, nn);
	    slist[nn] = 0;
	}
	for (int i = 0; i < nn; ++i) {
	    if (slist) sq = slist[i];
	    sq = sd->extseq(sq, ss->group[i], snl, ktree->lead[i].vol);
	    if (glst || (sq->sigII && sq->len)) {
		GAPS*   gps = sq->elim_column(DEL_GAP | RET_GAP);
		if (glst) {
		    delete[] (*glst)[i];
		    (*glst)[i] = gps;
		} else
		    delete[] gps;
	    }
	}
	if (glst)	glst->num = nn;
	if (!save)	delete sq;
	return slist;
}

#if USE_WEIGHT

void Ssrel::calcpw(Seq* sd)
{
	int	nn = ss? ss->num: sd->many;
	if (nn < 2) return;
	delete[] pairwt;
	ktree = new Ktree(sd, ss);
	pairwt = ktree->calcpw();
}

FTYPE Ssrel::sumofpwt()
{
	FTYPE	twt = 0;
	FTYPE*	pw = pairwt;
	int	n = ncomb(ss->num);

	while (n--) 	twt += *pw++;
	return twt;
}

#endif	// USE_WEIGHT

void Ssrel::fprint(FILE* fd)
{
	fprintf(fd, "%d\t%d\n", ss->num, ss->elms);
	for (int** gg = ss->group; *gg; ++gg) {
	    int*	g = *gg;
	    while (*g >= 0) fprintf(fd, "%d ", *g++);
	    fputc('\n', fd);
	}
	for (int i = 0; i < ss->num * 2 - 1; ++i)
	    fprintf(fd, "%2d %10.4f %10.4f\n", 
		ktree->lead[i].tid, ktree->lead[i].cur, ktree->lead[i].vol);
#if USE_WEIGHT
	for (int i = 0; i < ncomb(ss->num); ++i)
	    fprintf(fd, "%10.4f\n", pairwt[i]);
#endif	// USE_WEIGHT
}

Knode* ConservedT::rearrange(Knode* node)
{
	Knode*	lft = 0;
	Knode*	rht = 0;
	bool	ilf = leaf[node->tid] >= 0;

	if (!ilf) {
	    lft = rearrange(node->left);
	    rht = rearrange(node->right);
	    leaf[node->tid] = ntid++;
	}
	Knode*	targ = lead + leaf[node->tid];
	*targ = *node;
	targ->tid = leaf[node->tid];
	if (ilf) {
	    targ->left = targ->right = 0;
	} else {
	    targ->left = lft;
	    targ->right = rht;
	    targ->left->parent = targ->right->parent = targ;
	}
	return (targ);
}

/*  Dim. of sqs must >= 3   */
Ssrel*	Ssrel::consregss(mSeq** sqs)
{
	Ssrel*	trl = new Ssrel(ss, sqs[0]->many);
	ConservedT	ct(sqs, ss, trl->ss, true);
	ct.testssrel(ktree->root);
	ct.endstamp();
	int	k = ss? ss->num: sqs[0]->many;
	int	n = trl->ss->num = ct.grnum();

	ct.setntid(n);
	if (n < k) {
	    trl->ktree = new Ktree(n);
	    ct.setlead(trl->ktree->lead);
	    trl->ktree->root = ct.rearrange(ktree->root);
	} else {
	    delete trl;
	    trl = 0;
	}
	return (trl);
}

RANGE*	bluntreg(Seq* sd, RANGE* rng)
{
	RANGE*	wk = rng;
	RANGE*	br = rng;
	int	lmt = sd->left;

	br->right = -(lwrlmt + 1);
	while (neorng(++wk)) {
	    int	i = wk->left;
	    CHAR*	ss = sd->at(i);
	    if (!sd->nogap(ss))
		for ( ; i > lmt; --i)
		    if (sd->nogap(ss -= sd->many)) break;
	    if (i - br->right > lwrlmt)
		(++br)->left = i;
	    lmt = neorng(wk + 1)? (wk + 1)->left: sd->right;
	    i = wk->right;
	    ss = sd->at(i - 1);
	    if (!sd->nogap(ss))
		for ( ; i < lmt; i++)
		    if (sd->nogap(ss += sd->many)) break;
	    lmt = br->right = i;
	}
	++br;
	*br = endrng;
	rng->left = rng->right = ++br - rng;
	return (rng);
}

#if MAIN
#include "consreg.inc"
#endif
