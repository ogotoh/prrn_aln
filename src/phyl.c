/****************************************************************************
*
*	Subroutines for nucleotide/amino-acid sequence divergence
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
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japang
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@i.kyoto-u.ac.jp>>
*****************************************************************************/

#if NJ || UPG || WARD || PHYL
#define MAIN	1
#endif

#include "aln.h"
#include "wln.h"
#include "mseq.h"
#include "maln.h"
#include "mgaps.h"
#include "autocomp.h"
#include "phyl.h"
#include "qdiv.h"
#include "consreg.h"
#include "fspscore.h"

struct	Svr_Prm	{
	PwdM*	pwd;
	FTYPE*	dist;
	ScrRng*	scrrng;
	Svr_Prm(FTYPE* d, ScrRng* s = 0) 
	    : pwd(0), dist(d), scrrng(s) {}
	~Svr_Prm() {delete pwd;}
};

static	int	dvx(CalcServer<Seq>* svr, Seq* sqs[], ThQueue<Seq>* q);
static 	int	dpaln_main(CalcServer<mSeq>* svr, mSeq* sqs[], ThQueue<mSeq>* q);
static	int	dpscore(CalcServer<mSeq>* svr, mSeq* sqs[], ThQueue<mSeq>* q);
static	int	selfscr(CalcServer<mSeq>* svr, mSeq* sqs[], ThQueue<mSeq>* q);
static	void	scr2dist(Svr_Prm* svp, int num);
static	PwdM*	dpaln_reset(PwdM* pwd0, mSeq** sqs, ThQueue<mSeq>* q);
static	void	dpaln_setup(CalcServer<mSeq>* svr);
static	void	divseq2(FSTAT* stt, Seq* sd, int i, int j);
static	int	grpmem(Subset* ss, int i);
static	void	mins(int& ii, int& jj, FTYPE* d, FTYPE* r, int m);
static	void	minh(int& ii, int& jj, FTYPE* d, FTYPE* r, Knode** nodes);

TOUTMODE	treemode = {UPG_METHOD};

static	double	wfact = 0.;
static	bool	cfact = true;

static void divseq2(FSTAT* stt, Seq* sd, int i, int j)
{
	CHAR*	cs = sd->at(sd->left);
	CHAR*	ts = sd->at(sd->right);
	int	ga = 0;
	int	gb = 0;
	int	mch = 0;
	int	mmc = 0;
	int	unp = 0;
	int	gap = 0;

	for ( ; cs < ts; cs += sd->many) {
	    int	a = cs[i];
	    int	b = cs[j];
	    if (IsntGap(a)) {
		if (IsntGap(b)) {
		    ga = gb = 0;
		    if (a == b)	++mch;
		    else	++mmc;
		} else {
		    if (ga >= gb) ++gap;
		    ga = 0; ++gb; ++unp;
		}
	    } else {
		if (IsntGap(b)) {
		    if (ga <= gb) ++gap;
		    gb = 0; ++ga; ++unp;
		} else {
		    ++ga; ++gb;
		}
	    }
	}
	stt->mch = mch; stt->mmc = mmc;
	stt->gap = gap; stt->unp = unp;
}

void divseq(FSTAT* stat, Seq* sd, int* group1, int* group2)
{
	mSeq*	sqs[2];

	sqs[0] = new mSeq(sd, group1, CPY_SEQ);
	sqs[1] = new mSeq(sd, group2, CPY_SEQ);
	calcscore_grp(sqs, stat);
	clearseq(sqs, 2);
}

static int grpmem(Subset* ss, int i)
{
	if (ss && ss->elms > ss->num) {
	    int*	gr = ss->group[i];
	    for ( ; *gr >= 0; gr++) ;
	    int	n = gr - ss->group[i];
	    if (n == 1)
		return ss->group[i][0];
	    else if (distPrm.pickup_one)
		return ss->group[i][rand() % n];
	    else
		return -1;
	}
	return i;
}

static PwdM* dpaln_reset(PwdM* pwd0, mSeq** sqs, ThQueue<mSeq>* q)
{
	mSeq*&	a = sqs[0];
	mSeq*&	b = sqs[1];
	PwdM*	pwd = pwd0;
	if (!pwd) {
	    prePwd((Seq**) sqs);
	    pwd = new PwdM(sqs);
	} else if (a->many > 1 || b->many > 1) {
	    m_thread_Lock(q);
	    pwd = new PwdM(sqs);
	    m_thread_Unlock(q);
	} else {
	    a->mkthick();
	    b->mkthick();
	}
	return pwd;
}

static	void	dpaln_setup(CalcServer<mSeq>* svr)
{
	prePwd((Seq**) svr->in_face);
	Svr_Prm*	svp = (Svr_Prm*) svr->prm;
	svp->pwd = new PwdM(svr->in_face);
}

static int dpaln_main(CalcServer<mSeq>* svr, mSeq* sqs[], ThQueue<mSeq>* q)
{
	Svr_Prm*	svp = (Svr_Prm*) svr->prm;

	PwdM*	pwd = dpaln_reset(svp->pwd, sqs, q);

	Gsinfo	gsi;
	VTYPE	scr;
	SKL* skl = align2(sqs, pwd, &scr, &gsi);
	FTYPE	dst = pamcorrect(degdiv(&gsi), sqs[0]->inex.molc);
	mSeq*&   a = sqs[0];
	mSeq*&   b = sqs[1];
	int	al = a->left;
	int	bl = b->left;
	int	ar = a->right;
	int	br = b->right;
	if (algmode.lcl) {
	    skl = trimskl((Seq**) sqs, skl);
	    al = skl[1].m;
	    bl = skl[1].n;
	    ar = skl[skl->n].m;
	    br = skl[skl->n].n;
	}
	delete[] skl;

	ScrRng*&	srg = svp->scrrng;
	m_thread_Lock(q);
	svp->dist[elem(sqs[0]->sid, sqs[1]->sid)] = dst;
	if (srg) {
	    if (al < srg[a->sid].left) srg[a->sid].left = al;
	    if (bl < srg[b->sid].left) srg[b->sid].left = bl;
	    if (ar > srg[a->sid].right) srg[a->sid].right = ar;
	    if (br > srg[b->sid].right) srg[b->sid].right = br;
	}
	m_thread_Unlock(q);

	if (pwd != svp->pwd) delete pwd;
	return (OK);
}

static int dpscore(CalcServer<mSeq>* svr, mSeq* sqs[], ThQueue<mSeq>* q)
{
	Svr_Prm*	svp = (Svr_Prm*) svr->prm;

	int	lends[2];
	int*	ends = algmode.lcl? lends: 0;
	FTYPE	dst;
	if (sqs[0]->many == 1 && sqs[1]->many == 1)
	    dst = alnscore2dist((Seq**) sqs, svp->pwd, ends);
	else {
	    PwdM*	pwd = dpaln_reset(svp->pwd, sqs, q);
	    dst = alnscore2dist(sqs, pwd, ends);
	    if (pwd != svp->pwd) delete pwd;
	}
	mSeq*&   a = sqs[0];
	mSeq*&   b = sqs[1];
	int	al = a->left;
	int	bl = b->left;
	int	ar = a->right;
	int	br = b->right;
	if (ends) {
	    if (ends[0] > 0) bl += ends[0];
	    else if (ends[0] < 0) al -= ends[0];
	    if (ends[1] > 0) ar -= ends[1];
	    else if (ends[1] < 0) br += ends[1];
	}
	int	cmb = elem(a->sid, b->sid);
	ScrRng*&	srg = svp->scrrng;
	m_thread_Lock(q);
	svp->dist[cmb] = 100. * dst;
	if (srg) {
	    if (al < srg[a->sid].left) srg[a->sid].left = al;
	    if (bl < srg[b->sid].left) srg[b->sid].left = bl;
	    if (ar > srg[a->sid].right) srg[a->sid].right = ar;
	    if (br > srg[b->sid].right) srg[b->sid].right = br;
	}
	m_thread_Unlock(q);
	return (OK);
}

static int selfscr(CalcServer<mSeq>* svr, mSeq* sqs[], ThQueue<mSeq>* q)
{
	Svr_Prm*	svp = (Svr_Prm*) svr->prm;

	mSeq*&	a = sqs[0];
	ScrRng*	srg = svp->scrrng + a->sid;
	a->left = srg->left;
	a->right = srg->right;
	FTYPE	denom = (FTYPE) selfAlnScr(a);
	m_thread_Lock(q);
	srg->scr = denom;
	m_thread_Unlock(q);
	return (OK);
}

static int dvx(CalcServer<Seq>* svr, Seq* seqs[], ThQueue<Seq>* q)
{
	FSTAT	fst;

	vclear(&fst);
	CHAR*	ss = svr->msd->at(svr->msd->left);
	CHAR*	tt = svr->msd->at(svr->msd->right);
	int	agl = 0;
	int	bgl = 0;
	int	i = seqs[0]->sid;
	int	j = seqs[1]->sid;

	for ( ; ss < tt; ss += svr->msd->many) {
	    if (ss[i] == nil_code || ss[j] == nil_code) continue;
	    bool	ag = IsGap(ss[i]);
	    bool	bg = IsGap(ss[j]);
	    if (ag && bg) continue;
	    if (ag && !bg) {
		++fst.unp;
		if (agl <= bgl) ++fst.gap;
		++agl; bgl = 0;
	    } else if (!ag && bg) {
		++fst.unp;
		if (agl >= bgl) ++fst.gap;
		++bgl; agl = 0;
	    } else {
		if (ss[i] == ss[j]) ++fst.mch;
		else	++fst.mmc;
		agl = bgl = 0;
	    }
	}
	FTYPE	dst = pamcorrect(degdiv(&fst), seqs[0]->inex.molc);
	FTYPE*	dist = (FTYPE*) svr->prm;
	m_thread_Lock(q);
	dist[elem(i, j)] = dst;
	m_thread_Unlock(q);
	return (OK);
}

FTYPE* calcdist(Seq* sd, Subset* ss)
{
	int	nn = ss? ss->num: sd->many;
	if (nn < 2) return (0);
	FTYPE*	dist = 0;
	bool	nmult = distPrm.realign == DynAln;
	mSeq**	sbuf = nmult? new mSeq*[nn]: 0;
	if (nmult) {
	    initseq(sbuf, nn);
	    int	grp[2] = {0, -1};
	    for (int i = 0; i < nn; ++i) {
		grp[0] = i;
		sbuf[i]->extract(sd, ss? ss->group[i]: grp, CPY_SEQ);
	    }
	}
	if (distPrm.realign == ThisAln) {
	    dist = new FTYPE[ncomb(nn)];
	    AlnServer<Seq> svr(sd, IM_EVRY, ss, 1, (void*) dist, &dvx);
	    svr.auto_comp();
	} else if (distPrm.realign == Composition) {
	    dist = calcdist_kmer(sd, IM_EVRY, ss);
	} else if (distPrm.realign == GeneOrg){
	    dist = eijdmx(sd);
	} else if (distPrm.realign == DynAln) {
	    dist = new FTYPE[ncomb(nn)];
	    Svr_Prm	svp(dist);
	    AlnServer<mSeq> svr(sbuf, nn, IM_EVRY, (void*) &svp, 
		&dpaln_main, &dpaln_setup);
	    svr.auto_comp();
	} 
	if (nmult) {clearseq(sbuf, nn); delete[] sbuf;}
	return (dist);
}

// calculate the sum of distances from each member to others

FTYPE* calcdistsum(Seq* sd, Subset* ss, FTYPE* dist)
{
	int	nn = ss? ss->num: sd->many;
	if (nn < 2) return (0);
	bool	nmult = ss && ss->num < ss->elms && !distPrm.pickup_one;
	mSeq*	sbuf[2];
	if (nmult) initseq(sbuf, 2);
	if (!dist) dist = new FTYPE[nn];
	FTYPE*	dd = dist;
	for (int i = 0; i < nn; ++i, ++dd) {
	    int	grpi = 0;
	    if (nmult)	sbuf[0]->extract(sd, ss->group[i], CPY_SEQ);
	    else	grpi = grpmem(ss, i);
	    *dd = 0;
	    for (int j = 0; j < nn; ++j) {
		if (i == j) continue;
		FSTAT	stt;
		if (nmult) {
		    sbuf[1]->extract(sd, ss->group[j], CPY_SEQ);
		    calcscore_grp(sbuf, &stt);
		} else {
		    int	grpj = grpmem(ss, j);
		    divseq2(&stt, sd, grpi, grpj);
		}
		*dd += pamcorrect(degdiv(&stt), sd->inex.molc);
	    }
	}
	if (nmult) clearseq(sbuf, 2);
	return (dist);
}

FTYPE*	calcdist_i(Seq* sd, int k, Subset* ss, FTYPE* dist)
{
	int	nn = ss? ss->num: sd->many;
	if (nn < 2) return (0);
	bool	nmult = ss && ss->num < ss->elms && !distPrm.pickup_one;
	int	grpk = 0;
	mSeq*	sbuf[2];
	if (nmult) {
	    initseq(sbuf, 2);
	    sbuf[1]->extract(sd, ss->group[k], CPY_SEQ);
	} else
	    grpk = grpmem(ss, k);
	if (!dist) dist = new FTYPE[nn];
	FTYPE*	dd = dist;
	for (int i = 0; i < nn; ++i) {
	    if (i == k) {
		*dd++ = 0;
		continue;
	    }
	    FSTAT	stt;
	    if (nmult) {
		sbuf[0]->extract(sd, ss->group[i], CPY_SEQ);
		calcscore_grp(sbuf, &stt);
	    } else {
		int	grpi = grpmem(ss, i);
		divseq2(&stt, sd, grpi, grpk);
	    }
	    *dd++ = pamcorrect(degdiv(&stt), sd->inex.molc);
	}
	if (nmult) clearseq(sbuf, 2);
	return (dist);
}

/*
static void scr2dist(Svr_Prm* svp, int num)
{
	FTYPE*	d = svp->dist;
	ScrRng*	srg = svp->scrrng;
	for (int j = 1; j < num; ++j) {
	    for (int i = 0; i < j; ++i, ++d) {
		FTYPE denom = max(srg[i].scr, srg[j].scr);
		*d = 100. * (1. - *d / denom);
	    }
	}
}
*/

static void scr2dist(Svr_Prm* svp, int num)
{
	int	nn = ncomb(num);
	FTYPE*  d = svp->dist;
	FTYPE	dmax = d[vmax(d, nn) - d];
	FTYPE*	till = d + nn;
	for ( ; d < till; ++d) *d = dmax - *d;
}

DistMat::DistMat(int argc, const char** argv, const char* catalog, DistCal realin)
	: realign(realin), bias(0), scrrng(0)
{
	sname = new Strlist();
	InputSeqTest	tv = {sname, 0, 0, 0, 0, 0, 0};
	AlnServer<Seq> ssvr(argc, argv, IM_SNGL, catalog, (void*) &tv, &testInputSeq);
	ssvr.auto_comp(false);
	if (tv.bad) fatal("%d seqs were not found !\n", tv.bad);
	if (tv.num < 2) fatal("not enough data !\n");
	numb = tv.num;
	if (realign == DynAln || realign == DynScr) {
	    dist = new FTYPE[ncomb(numb)];
	    if (algmode.lcl) {
		scrrng = new ScrRng[numb];
		vset(scrrng, nullScrRng, numb);
	    }
	    Svr_Prm	svp(dist, scrrng);
	    AlnServer<mSeq> svr(argc, argv, IM_EVRY, catalog, (void*) &svp,
		realign == DynAln? &dpaln_main: &dpscore, &dpaln_setup);
	    svr.auto_comp();
	    if (algmode.lcl) {
		AlnServer<mSeq> svr2(argc, argv, IM_SNGL, catalog, 
		    (void*) &svp, &selfscr);
		svr2.auto_comp();
		scr2dist(&svp, numb);
	    }
	} else {
	    resetqdiv(ssvr.in_face[0]);
	    dist = calcdist_kmer(ssvr.in_face[0], argc, argv, catalog, numb, IM_EVRY);
	}
}

DistMat::DistMat(Seq** seqs, int nn) 
	: numb(nn), sname(0), bias(0)
{
	dist = calcdist_kmer(seqs, nn, IM_EVRY, sname);
}

DistMat::DistMat(mSeq** seqs, int nn, DistCal realin)
	: numb(nn), realign(realin), sname(0), bias(0), scrrng(0)
{
	if (realign == DynAln || realign == DynScr) {
	    convertseqs(seqs, nn, getSimmtx(0));
	    dist = new FTYPE[ncomb(numb)];
	    if (algmode.lcl) {
		scrrng = new ScrRng[numb];
		vset(scrrng, nullScrRng, numb);
	    }
	    Svr_Prm	svp(dist, scrrng);
	    CalcServer<mSeq> svr(IM_EVRY, (void*) &svp, 
		realign == DynAln? &dpaln_main: &dpscore, 
		&dpaln_setup, 0, seqs, nn);
	    svr.auto_comp();
	    if (algmode.lcl) {
		CalcServer<mSeq> svr2(IM_SNGL, (void*) &svp, &selfscr,
		    0, 0, seqs, nn);
		svr2.auto_comp();
		scr2dist(&svp, numb);
	    }
	} else {
	    dist = calcdist_kmer((Seq**) seqs, nn, IM_EVRY, sname);
	}
}

DistMat::~DistMat()
{
	delete[] dist;
	delete[] scrrng;
	delete sname;
}

/*****************************************************
*	class Knode
*****************************************************/

double Knode::calres()		// combined resistance
{
	if (isleaf()) return (res = 0.);
	double	rr = left->calres() + left->length;
	double	rl = right->calres() + right->length;
	return (res = (rr > 0. && rl > 0.)? rr * rl / (rr + rl): 0.);
}

Knode* Knode::findroot(double brl)
{
	Knode*	chng;
	Knode*	keep;

	height = (left->height + right->height + brl) / 2.;
	left->length = height - left->height;
	right->length = height - right->height;
	if (left->length < 0.) {
	    chng = left;
	    keep = right;
	} else if (right->length < 0.) {
	    chng = right;
	    keep = left;
	} else return (this);

	keep->length = brl;
	if (chng->left->height + chng->left->length >
	    chng->right->height + chng->right->length) {
	    left = chng->left;
	    chng->left = chng->right;
	} else {
	    left = chng->right;
	}
	brl = left->length;
	chng->right = keep;
	right = chng;
	double	hl = chng->left->height + chng->left->length;
	double	hr = chng->right->height + chng->right->length;
	chng->height = max(hl, hr);
	return (findroot(brl));
}

int Knode::teachparent()
{
	if (isleaf()) return (tid);
	left->parent = right->parent = this;
	int	l = left->teachparent();
	int	r = right->teachparent();
	ndesc = left->ndesc + right->ndesc;
	if (l > r) {
	    l = r;
	    gswap(left, right);
	}
	return (l);
}

Knode* Knode::findcenter()
{
	int	ndl = left? left->ndesc: 0;
	int	ndr = right? right->ndesc: 0;
	Knode*	node = ndl > ndr? left: right;
	Knode*	sist = ndl > ndr? right: left;
	int	nn = sist? sist->ndesc: 0;

	int	unb0 = abs(nn - ndl - ndr);
	int	unb1 = abs(nn + ndl - ndr);
	int	unb2 = abs(nn - ndl + ndr);
	if (unb0 <= unb1 && unb0 <= unb2) return (this);
	sist->length = 2 * height - sist->height - height;
	height = max(height, sist->height + sist->length);
	length = 0;
	sist->cur *= cur;
	cur = 1.;
	if (unb1 < unb2) {
	    left = node;
	    right = left;
	    left = sist;
	    ndesc = nn + ndl;
	} else {
	    right = node;
	    left = right;
	    right = sist;
	    ndesc = nn + ndr;
	}
	return (findcenter());
}

Knode* Knode::balanced()
{
	if (ndesc <= 3) return (this);
	Knode*	root = findcenter();
	root->left = root->left->findcenter();
	root->right = root->right->findcenter();
	return (root);
}

void Knode::kirchhof()
{
	if (isroot()) {
	    vol = res;
	    cur = 1.;
	} else {
	    double	pres = res + length;
	    cur = pres > 0.? parent->vol / pres: parent->cur / 2.;
	    vol = parent->vol - length * cur;
	}
	if (left) left->kirchhof();
	if (right) right->kirchhof();
}

/*****************************************************
*	class Ktree
*****************************************************/

Ktree::Ktree(int mem, Knode* led) : members(mem)
{
	lead = led;
	if (!led) {
	    lead = mem? new Knode[2 * members]: 0;
	    if (lead) vclear(lead, 2 * members);
	}
	leaf_no = 0; leaves = 0;
#if USE_WEIGHT
	 pwt = 0; bwt = 0.;
#endif
}

void Ktree::plant(DistTree& dt)
{
	lead = dt.lead; dt.lead = 0;
	root = dt.root;
	members = dt.members;
	leaf_no = 0; leaves = 0;
#if USE_WEIGHT
	pwt = 0; bwt = 0.;
#endif
}

Ktree::Ktree(DistMat* dmat, TreeMet tmethod, Knode* led)
{
	DistTree	dt(dmat, tmethod, led);
	plant(dt);
}

Ktree::Ktree(Seq* sd, Subset* ss, TreeMet tmethod, Knode* led) {
	DistTree	dt(sd, ss, tmethod, led);
	plant(dt);
}

#if USE_WEIGHT

FTYPE* Ktree::calcwt(FTYPE* wt)
{
	bwt = F_UNIF_WEIGHT;
	root->kirchhof();
	if (wt) {
	    FTYPE	nf = 1. / (1. + bwt);
	    for (int i = 0; i < members; ++i)
		wt[i] = (root->ndesc * lead[i].cur + lead[i].ndesc * bwt) * nf;
	}
	return (wt);
}

LEAF* Ktree::pairwt(Knode* node, double ros_)
{
	double	wbc, wac, wab;
	double	wfa = wfact;
	double	wfb = wfact;
	LEAF*	leaf = 0;
	LEAF*	reaf;

	node->ros = ros_;
	if (node->isleaf()) {
	    node->vol = node->parent->vol * node->cur;
	    if (leaves) {
		leaf = leaves + node->tid;
		leaf->wheight = node->vol + node->ndesc * bwt;
		leaf->tid = node->tid;
		leaf->lnk = 0;
	    }
	} else {
	    node->left->parent = node->right->parent = node;
	    double	a = node->left->res + node->left->length;
	    double	b = node->right->res + node->right->length;
	    if (node->isroot()) {
		node->cur = node->left->cur = node->right->cur = 1;
	    } else if (node->ros <= fepsilon || a + b <= fepsilon) {
		wab = 0.25; a = b = 0;
		wbc = wac = node->left->cur = node->right->cur = 0.5;
		node->vol = node->cur * node->parent->vol;
	    } else {
		if (a <= 0.) {b += a; a = fepsilon;}
		if (b <= 0.) {a += b; b = fepsilon;}
		double	c = node->length + node->ros;
 		wab = a * b / (a + b);
		wbc = a * (b + c);
		if (cfact) {
		    wfa = 1. + a * node->ros / ((wab + c) * (a + c));	
		    wfb = 1. + b * node->ros / ((wab + c) * (b + c));	
		}
		wab = wbc + b * c;
		wbc /= wab * wfb;
		wac = b * (a + c) / (wab * wfa);
		wab = c * (a + b) / wab;
		if (!cfact) wab /= wfa;
		a *= node->ros / (a + node->ros);
		b *= node->ros / (b + node->ros);
		node->cur *= sqrt(wac * wbc / wab);
		node->vol = node->cur * node->parent->vol;
		node->left->cur = sqrt(wab * wac / wbc);
		node->right->cur = sqrt(wab * wbc / wac);
	    }
	    leaf = pairwt(node->left, b);
	    reaf = pairwt(node->right, a);
	    if (!pwt || !leaves) return (leaf);
	    wab = 1. / (node->vol * node->vol);
	    LEAF*	lwk = leaf;
	    for (;; lwk = lwk->lnk) {
		for (LEAF* rwk = reaf; rwk; rwk = rwk->lnk)
		    pwt[elem(lwk->tid, rwk->tid)] = 
			wab * lwk->wheight * rwk->wheight;
		if (!lwk->lnk) break;
	    }
	    lwk->lnk = reaf;
	}
	return (leaf);
}

LEAF* Ktree::repairwt(Knode* node)
{
	LEAF*	leaf;
	LEAF*	reaf;

	if (node->isleaf()) {
	    leaf = leaves + leaf_no++;
	    leaf->wheight = node->vol + node->ndesc * bwt;
	} else {
	    leaf = repairwt(node->left);
	    reaf = repairwt(node->right);
	    double	wab = 1. / (node->vol * node->vol);
	    LEAF*	lwk = leaf;
	    for (;; lwk = lwk->lnk) {
		for (LEAF* rwk = reaf; rwk; rwk = rwk->lnk)
		    pwt[elem(lwk->tid, rwk->tid)] = 
			wab * lwk->wheight * rwk->wheight;
		if (!lwk->lnk) break;
	    }
	    lwk->lnk = reaf;
	}
	return (leaf);
}

FTYPE* Ktree::recalcpw(Knode* node, int num)
{
	if (!node) node = root;
	if (num == 0) num = node->ndesc;
	leaf_no = 0;

	leaves = new LEAF[num];
	for (int i = 0; i < num; ++i) {
	    leaves[i].tid = i;
	    leaves[i].lnk = 0;
	}
	FTYPE*	rsv = pwt;	// save old pwt
	pwt = new FTYPE[ncomb(num)];
	repairwt(node);
	delete[] leaves; leaves = 0;
	FTYPE*	rtn = pwt;
	pwt = rsv;		// restore old pwt
	return (rtn);		// return new pwt
}

FTYPE* Ktree::calcpw(FTYPE* wt)
{
	leaves = new LEAF[members];
	if (!pwt) pwt = new FTYPE[ncomb(members)];
	bwt = F_UNIF_WEIGHT / root->ndesc;
	root->vol = 1.;
	pairwt(root, fInfinit);
	if (wt) {
	    for (int i = 0; i < members; ++i)
		wt[i] = lead[i].vol;
	}
	delete[] leaves; leaves = 0;
	return (pwt);
}

FTYPE Ktree::sum_of_pairwt(Knode* node)
{
	if (!node) node = root;
	if (node->isleaf()) {
	    node->ros = node->vol;
	    return (0);
	} else {
	    VTYPE spw = sum_of_pairwt(node->left) + sum_of_pairwt(node->right);
	    node->ros = node->left->ros + node->right->ros;	// sum of weight
	    return (spw + node->left->ros *  node->right->ros / node->vol / node->vol);
	}
}

void calcweight(Seq* sd)
{
	if (!sd->weight) sd->weight = new FTYPE[sd->many];
	if (sd->many == 1) sd->weight[0] = 1.;
	else if (sd->many == 2) sd->weight[0] = sd->weight[1] = 0.5;
	else {
	    Ktree	tree(sd);
	    tree.calcwt(sd->weight);
	}
}

#endif	// USE_WEIGHT

/****************************************************************************
*
*	Dollo parsynomy
*
****************************************************************************/

void Ktree::dollo_1st(Knode* node, int j)
{
	if (node->isleaf()) {
	    node->stat = sgi->eijtab[node->tid][j]? 1: 0;
	} else {
	    dollo_1st(node->left, j);
	    dollo_1st(node->right, j);
	    node->stat = node->left->stat? 1: 0;
	    if (node->right->stat) ++node->stat;
	}
}

void Ktree::dollo_2nd(Knode* node, bool peist)
{
	bool	leaf = node->isleaf();
	if (!peist && ((leaf && node->stat == 1) || node->stat == 2)) {
	    ++node->gain;
	    peist = true;
	} else if (peist && node->stat == 0) {
	    ++node->loss;
	    peist = false;
	}
	if (leaf) return;
	dollo_2nd(node->left, peist);
	dollo_2nd(node->right, peist);
}

void Ktree::dollo(Seq* sd)
{
	sgi = sd->sigII;
	if (!sgi) return;
	sgi->mkeijtab(sd->many);
	for (int j = 0; j < sgi->pfqnum; ++j) {
	    dollo_1st(root, j);
	    dollo_2nd(root, false);
	}
}

/****************************************************************************
*
*	UPG Method to construct a Dendrogram
*
*	See DMX.DOC for detail.
*
*	Osamu Gotoh
*	Department of Biochemistry
*	Saitama Cancer Center Research Institute
*	Ina-machi, Saitama 362, Japan
*
*****************************************************************************/

FTYPE DistTree::dminidx(int& ii, int members)
{
	ii = row[0];
	FTYPE	dmin = distance(ii, nnbr[ii]);

	for (int j = 0; ++j < members; ) {
	    int	jj = row[j];
	    FTYPE	dij = distance(jj, nnbr[jj]);
	    if (dij < dmin) {
		ii = jj;
		dmin = dij;
	    }
	}
	return (dmin);
}

int DistTree::dminrow(int ii, int members)
{
	FTYPE	dmin = FLT_MAX;
	int	jj = ii;
	for (int k = 0; k < members; ++k) {
	    int	kk = row[k];
	    if (kk == ii) continue;
	    FTYPE	dij = distance(ii, kk);
	    if (dij < dmin) {
		dmin = dij;
		jj = kk;
	    }
	}
	return (jj);
}

Knode* DistTree::upg_method()
{
	if (members < 2) fatal("No of OTUs must be >= 2!\n");
	Knode**	nodes = new Knode*[members];
	row = new int[members];
	nnbr = new int[members];

	nnbr[0] = 1;
	vclear(nnbr + 1, members - 1);
	int	m = 0;
	for (int k = 0; m < members; ++m) {
	    nodes[m] = lead + m;
	    nodes[m]->tid = row[m] = m;
	    nodes[m]->ndesc = max(nodes[m]->ndesc, 1);
	    for (int n = 0; n < m;  ++n, ++k) {
		if (dist[k] < distance(m, nnbr[m])) nnbr[m] = n;
		if (dist[k] < distance(n, nnbr[n])) nnbr[n] = m;
	    }
	}

	root = *nodes;
	for (int n = members; n-- > 1; ) {
	    int	ii = 0;
	    FTYPE dij = dminidx(ii, n);
	    int	jj = nnbr[ii];
	    root = lead + m;
	    root->tid = m++;
	    Knode*&	lnode = root->left = nodes[ii];
	    Knode*&	rnode = root->right = nodes[jj];
	    root->height = dij / 2;
	    lnode->length = root->height - lnode->height;
	    if (lnode->length < 0) lnode->length = 0;
	    rnode->length = root->height - rnode->height;
	    if (rnode->length < 0) rnode->length = 0;
	    FTYPE	rl = lnode->res + root->height - lnode->height;
	    FTYPE	rr = rnode->res + root->height - rnode->height;
	    root->res = (rl > fepsilon && rr > fepsilon)? (rl * rr) / (rl + rr): fepsilon;
	    root->ndesc = lnode->ndesc + rnode->ndesc;
	    int	j = 0;
	    for (int k = 0; k <= n; ++k) {
		nnbr[ii] = -1;
		int	kk = row[k];
		if (kk == ii) continue;
		if (kk == jj) {j = k; continue;}
		FTYPE&	x = distance(kk, ii);
		FTYPE	y = distance(kk, jj);
		Knode*&	knode = nodes[kk];
		switch (treemode.method) {
		  case UPG_METHOD:	// upgma
		    x *= lnode->ndesc;
		    x += rnode->ndesc * y;
		    x /= root->ndesc;
		    break;
		  case SLINK_METHOD:	// single linkage
		    if (y < x) x = y; break;
		  case CLINK_METHOD:	// complete linkage
		    if (y > x) x = y; break;
		  case WARD_METHOD:	// Ward's minimum variance
		    x = x * (lnode->ndesc + knode->ndesc) +
			y * (rnode->ndesc + knode->ndesc) -
			distance(ii, jj) * knode->ndesc;
		    x /= (lnode->ndesc + rnode->ndesc + knode->ndesc);
		  default:		// simple mean
		    x = (x + y) / 2; break;
		}			// indicator for recalculation
		if (treemode.method != SLINK_METHOD && 
			(nnbr[kk] == ii || nnbr[kk] == jj))
		    nnbr[kk] = -1;
		else if (nnbr[kk] == jj) nnbr[kk] = ii;
	    }
	    nodes[ii] = root;
	    row[j] = row[n];
	    for (int k = 0; k < n; ++k) {
		int	kk = row[k];
		if (nnbr[kk] < 0) nnbr[kk] = dminrow(kk, n);
	    }
	}
	delete[] nodes;
	root->teachparent();
	root->parent = 0;
	return (root);
}

/*********************************************************************
*
*	Neighbor-Join Method to construct a Phylogenic Tree
*	Ref. (1):	Saitou, N. and Nei, M. (1987) 
*					Mol. Biol. Evol. 4, 406-425
*	Ref. (2):	Studier, J.A. and Keppler, K.L. (1988)
*					Mol. Biol. Evol. 5, 729-731
*	Comment:
*		R[k] in Ref. 2 is now gotten by recursion
*		R[k] = R[k] - (D[i][k] + D[j][k] + D[i][j]) / 2  (k != i, j)
*		R[k] = (R[i] + R[j] - m * D[i][j]) / 2		 (k = new node)
*
*	See DMX.DOC for detail.
*
*	Osamu Gotoh
*	Department of Biochemistry
*	Saitama Cancer Center Research Institute
*	Ina-machi, Saitama 362, Japan
*
***********************************************************************/

void DistTree::lowesthi(Knode* node, double hi)
{
	node->height = hi -= node->length;
	if (treemode.n_edge || (!ntimes && treemode.unrooted)) {
	    if (node->left && node->left->length < 0.) {
		node->right->length -= node->left->length;
		node->left->length = 0.;
	    }
	    if (node->right && node->right->length < 0.) {
		node->left->length -= node->right->length;
		node->right->length = 0.;
	    }
	}
	ntimes++;
	if (node->height < lwhi) lwhi = node->height;
	if (node->left) lowesthi(node->left, hi);
	if (node->right) lowesthi(node->right, hi);
}

double DistTree::recalhi(Knode* node, double hi)
{
	lwhi = ntimes = 0;
	lowesthi(node, hi);
	return (lwhi);
}

static void mins(int& ii, int& jj, FTYPE* d, FTYPE* r, int m)
{
	double	fm2 = (double) (m - 2);
	double	fmin = *d++ * fm2 - r[0] - r[1];

	ii = 0;
	jj = 1;
	for (int j = 2; j < m; j++) {
	    for (int i = 0; i < j; i++, d++) {
		double	tmp = *d * fm2 - r[i] - r[j];
		if (tmp < fmin) {
		    ii = i;
		    jj = j;
		    fmin = tmp;
		}
	    }
	}
}

static void minh(int& ii, int& jj, FTYPE* d, FTYPE* r, Knode** nodes)
{
	r += 2;
	nodes += 2;
	double	hmax = 2. * (*nodes)->height + *r - *d;
	ii = 0;
	jj = 1;
	for (int i = 0; i < 2; i++) {
	    ++d; --r; --nodes;
	    double	tmp = 2. * (*nodes)->height + *r - *d;
	    if (tmp > hmax) {
		ii = i;
		jj = 2;
		hmax = tmp;
	    }
	}
}

Knode* DistTree::nj_method()
{
	if (members < 3) fatal("No of OTUs must be >= 3!\n");
	FTYPE*	sum = new FTYPE[members];
	Knode**	nodes = new Knode*[members];

	int m = 0;
	for ( ; m < members; m++) {
	    nodes[m] = lead + m;
	    nodes[m]->tid = m;
	    nodes[m]->ndesc = max(nodes[m]->ndesc, 1);
	    for (int k = (int) (sum[m] = 0); k < members; k++)
		if (k != m) sum[m] += dist[elem(k, m)];
	}

	for (int n = members; n >= 3; ) {
	    int	i, j;
	    if (n > 3)	mins(i, j, dist, sum, n);	//	i < j
	    else 	minh(i, j, dist, sum, nodes);
	    double	dd  = (sum[i] - sum[j]) / (double) (n - 2);
	    double	dij = dist[elem(i, j)];
	    double	hl = (dij + dd) / 2.;
	    double	hr = (dij - dd) / 2.;
	    sum[i] += sum[j] - n * dij;
	    sum[i] /= 2.;
	    root = lead + m;
	    root->tid = m++;
	    root->left = nodes[i];
	    root->right = nodes[j];
	    root->left->length = hl;
	    root->right->length = hr;
	    root->ndesc = nodes[i]->ndesc + nodes[j]->ndesc;
	    hl += root->left->height;
	    hr += root->right->height;
	    root->height = max(hl, hr);
	    nodes[i] = root;
	    for (int k = 0; k < n; k++) {
		if (k == i || k == j) continue;
		dd = dist[elem(k, i)] + dist[elem(k, j)];
		dist[elem(k, i)] = (dd - dij) / 2.;
		sum[k] -= (dd + dij) / 2.;
	    }
	    if (j != --n) {		//	Move last node to j-th
		sum[j] = sum[n];
		nodes[j] = nodes[n];
		for (int k = 0; k < n; k++)
		    if (k != j) dist[elem(k, j)] = dist[elem(k, n)];
	    }
	}
	root = lead + m;
	root->tid = m;
	root->left = nodes[0];
	root->right = nodes[1];
	root->length = 0.;
	root->ndesc = members;
	root = root->findroot(dist[0]);
	root->parent = 0;
	root->teachparent();
	double	hl = recalhi(root, root->height);
	if (hl < 0.) recalhi(root, root->height - hl);
	root->calres();
	delete[] sum;
	delete[] nodes;
	return (root);
}

/****************************************************************************
*
*	The parts commonly used by distance matrix methods for
*	constructing phylogenetic trees
*
*	See DMX.DOC for detail.
*
*	Osamu Gotoh
*	Department of Biochemistry
*	Saitama Cancer Center Research Institute
*	Ina-machi, Saitama 362, Japan
*
*****************************************************************************/


TreeMet setphyl(int method)
{
retry:
	switch (method) {
	    case UPG_METHOD: case SLINK_METHOD: case CLINK_METHOD: 
	    case NJ_METHOD: case WARD_METHOD:
		treemode.method = TreeMet(method);
		break;
	    case QUERY:	
		method = (int) treemode.method;
		promptin("Tree method UPG[0]/NJ[1]/SLINK[2]/CLINK[3]/WARD[4] (%d) : ", &method);
		goto retry;
	    case SILENT: break;
	    default:	goto retry;
	}
	return (treemode.method);
}

void setpwfact(double wf)
{
	wfact = wf;
	cfact = wf <= 0.;
}

void DistTree::maketree(TreeMet met)
{
	lwhi = ntimes = 0;
	mtrx = treemode.ddmx? vcopy((FTYPE*) 0, dist, ncomb(members)): 0;
	if (met == Def_METHOD) met = treemode.method;
	root = (met == NJ_METHOD)? nj_method(): upg_method();
}

DistTree::DistTree(DistMat* dmat, TreeMet tmethod, Knode* led)
	: Ktree(dmat->numb, led), dist(dmat->dist),
	  row(0), nnbr(0), newdmat(false), mname(0), mnbuf(0), mtrx(0) 
{
	if (dmat->sname) {
	    mname = new char*[members + 1];
	    for (int i = 0; i < dmat->numb; ++i)
		mname[i] = (*dmat->sname)[i];
	}
	maketree(tmethod);
}

DistTree::DistTree(Seq* sd, Subset* ss, TreeMet tmethod, Knode* led, bool nm)
	: Ktree(ss? ss->num: sd->many, led), 
	  row(0), nnbr(0), newdmat(true), mname(0), 
	  mnbuf(0), ntimes(0), lwhi(0), mtrx(0)
{
	if (nm) {
	    mname = new char*[members + 1];
	    char**	mnm = mname;
	    if (ss) {
		for (int** gr = ss->group; *gr; ++gr)
		   *mnm++ = (*sd->sname)[**gr];
	    } else {
		for (int i = 0; i < sd->many; ++i)
		    *mnm++ = (*sd->sname)[i];
	    }
	}
	dist = calcdist(sd, ss);
	maketree(tmethod);
}

DistTree::DistTree(const char* user_tree)
	: Ktree(), dist(0), row(0), nnbr(0), newdmat(false), mname(0), 
	 mnbuf(0), ntimes(0), lwhi(0), mtrx(0)
{
	Btree<Knode>	st(user_tree);
	members = st.members;
	mname = st.mname; st.mname = 0;
	mnbuf = st.mnbuf; st.mnbuf = 0;
	lead = st.lead; st.lead = 0;
	root = st.root;
}
	
LEAF* DistTree::estimatedist(Knode* node)
{
	LEAF*	leaf;
	LEAF*	reaf;

	if (node->isleaf()) {
	    leaf = leaves + node->tid;
	    leaf->tid = node->tid;
	    leaf->lnk  = 0;
	    leaf->wheight = node->height;
	} else {
	    leaf = estimatedist(node->left);
	    reaf = estimatedist(node->right);
	    LEAF*	lwk = leaf;
	    while (lwk) {
		for (LEAF* rwk = reaf; rwk; rwk = rwk->lnk)
		    mtrx[elem(lwk->tid, rwk->tid)] = 
		    2 * node->height - lwk->wheight - rwk->wheight;
		lwk = lwk->lnk;
	    }
	    lwk->lnk = reaf;
	}
	return (leaf);
}

#if MAIN

static	int	clmwd =  CLMWD;
static	int	treewd = CLMWD - OTULEN - LMARGIN;
static	int	horline = HORLINE;
static	int	verline = VERLINE;

DistMat::DistMat(FILE* fd)
	: numb(0), bias(0)
{
	Mfile	mfd(sizeof(char *));

	prompt("Input members\n");
	sname = new Strlist();
	while (true) {
	    char    str[MAXL];
	    prompt("%3d: ", &numb);
	    if (!fgets(str, MAXL, fd) || isBlankLine(str)) break;
	    int	m = strlen(str);
	    str[m - 1] = '\0';
	    sname->push(str);
	    ++numb;
	}
	dist = new FTYPE[ncomb(numb)];
	FTYPE*  d = dist;
	prompt("Input distance\n");
	for (int i = 1; i < numb; i++) {
	    for (int j = 0; j < i; j++, d++) {
		prompt("%s - %s : ", (*sname)[i], (*sname)[j]);
		float	v;
		if (fscanf(fd, "%f", &v) <= 0)
		    fatal("Insuficient Dist. Values");
		*d = treemode.negate? -v: v;
		if (*d < bias) bias = *d / 2.;
	    }
	}
}

void DistTree::putdmx(FILE* fo)
{
	double	rms = 0.;
static	const	char* space8 = "		";

	putc('\n', fo);
	fputs(space8, fo);
	for (int i = 0; i < members; i++)
	    fprintf(fo, "%8s", mname[i]);
	putc('\n', fo);
	for (int i = 0; i < members; i++) {
	    fprintf(fo, "%8s", mname[i]);
	    for (int j = 0; j < members; j++) {
		int	e = elem(i, j);
		if (i < j) {
		    fprintf(fo, "%8.2f", dist[e]);
		    double	d = dist[e] - mtrx[e];
		    rms += d * d;
		} else if (i > j)
		    fprintf(fo, "%8.2f", mtrx[e]);
		else	fputs(space8, fo);
	    }
	    putc('\n', fo);
	}
	fprintf(fo, "\n%sRms = %8.3f\n", space8, sqrt(rms / ncomb(members)));
}

void PpPrm::prepare()
{
	minelem = 1;
	lastelem = 0;
	factor = 1.;
	maxhi = 0.;
	for (int i = 0; i < MAXELEM; ++i)
	    list[i].idx = EOL;
	list[0].data = 0.;
	linebuf = new char[clmwd + 1];
	baseline = linebuf;
}

PpPrm::PpPrm(FILE* fd, Seq* sd, Subset* ss, TreeMet tmethod, double bs) 
	: DistTree(sd, ss, tmethod, 0, true), fo(fd), bias(bs)
{
	prepare();
}

PpPrm::PpPrm(FILE* fd, DistMat* dmat, TreeMet tmethod, double bs) 
	: DistTree(dmat, tmethod), fo(fd), bias(bs)
{
	prepare();
}

PpPrm::PpPrm(FILE* fd, const char* user_tree, double bs)
	: DistTree(user_tree), fo(fd), bias(0)
{
	prepare();
}

void PpPrm::add(double data)
{
	if (minelem < MAXELEM) {
	    list[minelem].data = data;
	    lastelem = list[lastelem].idx = minelem;
	    while (++minelem < MAXELEM && list[minelem].idx != EOL) ;
	}
	else	fatal("Array Overflow !\n");
}

void PpPrm::del(double data)
{
	int	j;

	for (int i = 0; (j = list[i].idx) != EOL; i = j) {
	    if (list[j].data == data) {
		list[i].idx = list[j].idx;
		list[j].idx = EOL;
		if (j == lastelem) lastelem = i;
		minelem = min(minelem, j);
		break;
	    }
	}
}

double PpPrm::nextdat()
{
	if ((datp = list[datp].idx) == EOL)
	    return EOS;
	else	return list[datp].data;
}

double PpPrm::firstdat()
{
	datp = 0;
	return nextdat();
}

void PpPrm::fillstr(int cha, int from, int to)
{
	int i = 0;
	while (i < from)	baseline[i++] = ' ';
	while (i < to)		baseline[i++] = cha;
	while (i < treewd)	baseline[i++] = ' ';
}

void PpPrm::wdadjust()
{
	otulen = OTULEN;
	int	n = members;
	char**	nams = mname;
	while (n--) {
	    int	m = strlen(*nams++);
	    if (m >= otulen) otulen = m;
	}
	int 	m = (int) (otulen - log10(maxhi));
	if (m < 0) m = 0;
	if (m > 5) m = 5;
	if (treemode.lroot) {
	    sprintf(sfrmt, " %%%ds", otulen);
	    sprintf(ffrmt, " %%%d.%df", otulen, m);
	} else {
	    sprintf(sfrmt, "%%%ds ", otulen);
	    sprintf(ffrmt, "%%%d.%df ", otulen, m);
	}
	treewd = clmwd - otulen - LMARGIN;
	if (treemode.weight) treewd -= 10;
	factor = (FTYPE) treewd / (maxhi - bias);
}

void PpPrm::superimp(int a, int z, double len)
{
	char	buf[32];
	char*	p = buf;
	char*	str = baseline + ++a;

	sprintf(buf, "%.2f", len);
	z = (z - a - (int) strlen(buf)) / 2;
	if (z > 0) str += z;
	while (*p) *str++ = *p++;
	if (str >= baseline + treewd) *str = '\0';
}

void PpPrm::superimp(int a, int z, int gain, int loss)
{
	char	buf[32];
	char*	p = buf;
	char*	str = baseline + ++a;

	if (gain && loss)
	    sprintf(buf, "+%d -%d", gain, loss);
	else if (gain)
	    sprintf(buf, "+%d", gain);
	else if (loss)
	    sprintf(buf, "-%d", loss);
	z = (z - a - (int) strlen(buf)) / 2;
	if (z > 0) str += z;
	while (*p) *str++ = *p++;
	if (str >= baseline + treewd) *str = '\0';
}

void PpPrm::superimpnn(int a, int z, int nn)
{
	char	buf[32];
	char*	p = buf;
	char*	str; 

	sprintf(buf, "%d", nn);
	if (a > z && (a -= strlen(buf) + 2) < z) a = z;
	str = baseline + ++a;
	while (*p) *str++ = *p++;
	if (str >= baseline + treewd) *str = '\0';
}

void PpPrm::printcur(Knode* node)
{
	if (!node) return;
	fprintf(fo, "%3d\t%15.7e\n", node->tid, node->cur);
	printcur(node->left);
	printcur(node->right);
}

void PpPrm::printpar(Knode* node, bool newpar)
{
	int	l_gt_r = 0;

	if (node->isleaf()) {
	    if (treemode.branch)
		sprintf(str, "%s:%.3f", mname[node->tid], node->length);
	    else
		sprintf(str, "%s", mname[node->tid]);
	    fputs(str, fo);
	    clmpos += strlen(str);
	    return;
	}
	if (newpar) {
	    putc('(', fo);
	    clmpos++;
	}
	if (treemode.unrooted && node->isroot()) 
	    l_gt_r = node->left->length >= node->right->length? 1: -1;
	if (node->left) printpar(node->left, l_gt_r >= 0);
	putc(',', fo);
	clmpos++;
	if (clmpos >= treewd) {
	    putc('\n', fo);
	    clmpos = 0;
	}
	if (node->right) printpar(node->right, l_gt_r <= 0);
	if (newpar) {
	    putc(')', fo);
	    clmpos++;
	}
	if (treemode.branch && !node->isroot()) {
	    if (treemode.branch) sprintf(str, ":%.3f", node->length);
	    fputs(str, fo);
	    clmpos += strlen(str);
	}
}

void PpPrm::printnode(Knode* node, double father, int lr)
{
	int 	from = normal(0., treemode.lroot);
	int 	temp = normal(father, treemode.lroot);
static	const	char*	frmcur = "%7.5f ";

	if (node->left) printnode(node->left, node->height, LEFT);
	if (lr == RIGHT) del(father);
	if (!node->isleaf() || !treemode.fixedbase)
	    from = normal(node->height, treemode.lroot);
	if (treemode.prnode)	superimpnn(from, temp, node->tid);
	if (from > temp) gswap(from, temp);
	if (treemode.branch && node->length != 0.)
	    superimp(from, temp, node->length);
	else if (treemode.height && node->parent)
	    superimp(from, temp, root->height - node->height);
	else if (treemode.dollo && (node->gain || node->loss))
	    superimp(from, temp, node->gain, node->loss);
	fputs(linebuf, fo);
	putc('\n', fo);

	baseline[treewd] = baseline[treewd+1] = '\0';
	if (node->isleaf()) {
	    if (treemode.lroot)
		sprintf(baseline + treewd, sfrmt, mname[node->tid]);
	    else {
		sprintf(linebuf, sfrmt, mname[node->tid]);
		if (treemode.weight)
		    sprintf(linebuf + otulen + 1, frmcur, node->cur);
		baseline[-1] = ':';
	    }
	} else if (!treemode.lroot) {
	    sprintf(linebuf, ffrmt, node->height);
	    char*	ps = linebuf;
	    while (*ps)  ps++;
	    while (ps < baseline) *ps++ = ' ';
	}
	if (treemode.lroot) {from++; temp++;}
	fillstr(horline, from, temp);
	for (double data = firstdat(); data != EOS; data = nextdat()) {
	    int	temp = normal(data, treemode.lroot);
	    baseline[temp] = verline;
	}
	if (lr == LEFT) add(father);
	if (node->right) printnode(node->right, node->height, RIGHT);
}

void  PpPrm::puttree()
{
	if (treemode.paren) {
	    clmpos = 0;
	    printpar(root, true);
	    putc(';', fo); putc('\n', fo);
	} else {
	    maxhi = root->height;
	    wdadjust();
	    if (!treemode.lroot) {
		baseline += otulen + 2;
		if (treemode.weight) baseline += 8;
	    }
	    fputc('\n', fo);
	    memset(linebuf, ' ', baseline - linebuf + treewd);
	    baseline[treewd] = '\0';
	    printnode(root, root->height, 0);
	    fputs(linebuf, fo);
	    putc('\n', fo);
	}
}

void PpPrm::putheight(Knode* node)
{
	if (node->isleaf()) {
	    printf("%-15s\t%7.2f\n", mname[node->tid], root->height - node->height);
	} else {
	    putheight(node->left);
	    putheight(node->right);
	}
}

void PpPrm::putwt()
{
	FTYPE*	wt = new FTYPE[members];

	calcwt(wt);
	for (int i = 0; i < members; ++i) {
	    if (treemode.method == UPG_METHOD)
		fprintf(fo, "%15s %15.7e\n",  mname[i], wt[i]);
	    else
		fprintf(fo, "%-15s %15.7e %15.7e\n",
		     mname[i], wt[i], root->height - lead[i].height);
	}
	delete[] wt;
}

void PpPrm::putpw()
{
	int	k = ncomb(members);
	double	sum = 0;

	delete[] pwt; pwt = 0;
	calcpw();
	if (treemode.balance) root = root->balanced();
	if (treemode.branch) printcur(root);
	FTYPE*	pw = pwt;
	for (int i = 0; i < k; ++i) sum += *pw++;
	sum = k / sum;
	pw = pwt;
	for (int i = 0; i < k; ++i)
	    fprintf(fo, "%15.7e\n", *pw++ * sum);
	delete[] pwt; pwt = 0;
}

LEAF* PpPrm::tracepath(Knode* node)
{
	LEAF*	leaf = 0;
	LEAF*	reaf = 0;

	if (node->isleaf()) {
	    leaf = leaves + node->tid;
	    leaf->tid = node->tid;
	    leaf->lnk  = 0;
	} else {
	    node->left->parent = node->right->parent = node;
	    leaf = tracepath(node->left);
	    reaf = tracepath(node->right);
	    LEAF*	lwk = leaf;
	    while (lwk) {
		for (LEAF* rwk = reaf; rwk; rwk = rwk->lnk) {
		    Knode*	wnd = lead + lwk->tid;
		    fprintf(fo, "%d\t", elem(lwk->tid, rwk->tid));
		    for ( ; wnd->tid != node->tid; wnd = wnd->parent)
			fprintf(fo, "%d ", wnd->tid);
		    wnd = lead + rwk->tid;
		    for ( ; wnd->tid != node->tid; wnd = wnd->parent)
			fprintf(fo, "%d ", wnd->tid);
		    putc('\n', fo);
		}
		lwk = lwk->lnk;
	    }
	    lwk->lnk = reaf;
	}
	return (leaf);
}

static	const	char*	pronam;

#if PHYL

static	const	char*	groups = 0;
static	void	phyln(Seq* msd);
static	const	char*	tree_file = 0;

void usage()
{
	char	str[MAXL];

	fprintf(stderr, "Usage:\t%s [-options] [aligned_seq]\n", 
	    partfnam(str, pronam, "b"));
	fputs("Options:\n", stderr);
	fputs("  -aN:\tRealign pairwise [0-2].\t", stderr);
	fputs("  -b: \tPrint branch length.\n", stderr);
	fputs("  -cN:\tSet column width to N.\t", stderr);
	fputs("  -d: \tPrint distance matrix.\n", stderr);
	fputs("  -e: \tSuppress negative edge.\t", stderr);
	fputs("  -gN:\tGap equivalence [0-3].\n", stderr);
	fputs("  -hC:\tHorizontal line char=C.\t", stderr);
	fputs("  -kN:\tCorrect multi-hits [0-2].\n", stderr);
	fputs("  -n: \tPrint node number.\t", stderr);
	fputs("  -r: \tRoot is put left.\n", stderr);
	fputs("  -t: \tNew Hapmshire standard.\t", stderr);
	fputs("  -u: \tUnrooted tree.\n", stderr);
	fputs("  -vC:\tVertical line char=C.\t", stderr);
	fputs("  -w: \tWeight without tree.\n", stderr);
	fputs("  -wp: \tPair weight No tree.\t", stderr);
	fputs("  -y: \tWeight with tree.\n", stderr);
	fputs("  -B: \tGene organizaion.\t", stderr);
	fputs("  -D: \tDollo parsymony on Gene Organization.\n", stderr);
	fputs("  -C: \tComplete linkage method.", stderr);
	fputs("  -N: \tNJ method.\n", stderr);
	fputs("  -W: \tWard method.\n", stderr);
	fputs("  -S: \tSingle linkage method.\t", stderr);
	fputs("  -T: \tNewick format.\n", stderr);
	fputs("  -U: \tUPGMA method.\t\t", stderr);
	fputs("  -?: \tThis.\n", stderr);
	exit(0);
}

template <class seq_t>
int AlnServer<seq_t>::localoption(int& argc, const char**& argv)
{
	const	char*	opt= argv[0] + 1;
	const	char*	val = argv[0] + 2;
	int	c = *opt;
	DistMes	dm = Qpamd;
	int	n = 1;

	switch (c) {
	  case 'b':   
	    if (*val == '-') treemode.branch = false; else
	    if (*val == '+') treemode.branch = true; else
	    treemode.branch = !treemode.branch;
	    break;
	  case 'c':
	    clmwd = atoi(getarg(argc, argv, true));
	    break;
	  case 'd':
	    treemode.tree = (*val == 't');
	    treemode.ddmx = true; break;
	  case 'e': treemode.n_edge = true; break;
	  case 'f':
	    wfact = atof(getarg(argc, argv, true));
	    if (wfact > 0.) cfact = false;
	    break;
	  case 'h': treemode.tree = (*val == 't'); treemode.height = 1; break;
	  case 'j': treemode.balance = true; break;
	  case 'l': horline = *val? *val: ' '; break;
	  case 'm':
	    dm = (DistMes) atoi(getarg(argc, argv, true));
	    setqdiv(0, 0, 0, dm);
	    break;
	  case 'n': treemode.prnode = true; break;
	  case 'q': algmode.qck = 1; break; // calc dist by k-mer count
	  case 'r': treemode.lroot = true; break;
	  case 'u': treemode.unrooted = true; break;
	  case 'v': verline = *val? *val: ' '; break;
	  case 'w':
	    treemode.tree = (*val == 't');
	    if (*val == '2' || *val == 'p') treemode.pairwt = true;
	    else	treemode.weight = true;
	    break;
	  case 'z': setexprm_z(val); break;
	  case 'B': distPrm.realign = GeneOrg; break;
	  case 'C': setphyl(CLINK_METHOD); break;
	  case 'D': treemode.dollo = 1; break;
	  case 'G': groups = val; break;
	  case 'M':
	    switch(*val) {
	      case '\0':
	      case 'c': case 'C': setphyl(CLINK_METHOD); break;
	      case 'n': case 'N': setphyl(NJ_METHOD); break;
	      case 's': case 'S': setphyl(SLINK_METHOD); break;
	      case 'w': case 'W': setphyl(WARD_METHOD); break;
	      default:	break;
	    }
	  case 'N': setphyl(NJ_METHOD); break;
	  case 'P':
	   treemode.tree = (*val == 't');
	    treemode.path = true; break;
	  case 'S': setphyl(SLINK_METHOD); break;
	  case 'T':treemode.paren = treemode.branch = true; break;
	  case 'U': setphyl(UPG_METHOD); break;
	  case 'W': setphyl(WARD_METHOD); break;
	  case 'Y':		// user tree
	    if (*val) tree_file = val;
	    else if (argv[1][0] != OPTCHAR) {
		tree_file = *++argv; --argc;
	    }
	    break;
	  case '?': usage();
	  default:  n = 0; break;
	}
	return (n);
}

static void phyln(Seq* msd)
{
	Subset*	ss = groups? new Subset(msd->many, groups): 0;

	PpPrm*	dtree = 0;
	Seq*	tmp = 0;
	if (tree_file) {
	    dtree = new PpPrm(out_fd, tree_file);
	    int*	which = new int[msd->many + 1];
	    int*	badmn = new int[msd->many];
	    int		bnm = 0;
	    for (int i = 0; i < msd->many; ++i) {
		which[i] = msd->sname2memno(dtree->lname(i));
		if (which[i] < 0) badmn[bnm++] = i;
	    }
	    if (bnm) {
		fputs("Msa and tree are incompatible !!\n", stderr);
		for (int i = 0; i < bnm; ++i)
		    fprintf(stderr, "%s\n", dtree->lname(i));
		delete dtree; delete[] which; delete[] badmn;
		exit(1);
	    }
	    which[msd->many] = -1;
	    tmp = new Seq(msd, which, CPY_ALL | RLT_SBI);
	    delete[] which;
	    delete[] badmn;
	    gswap(msd, tmp);
	} else {
	    dtree = new PpPrm(out_fd, msd, ss, Def_METHOD);
	}
	delete ss;
	if (treemode.dollo) dtree->dollo(msd);
	if (treemode.weight) {
	    if (treemode.tree) {
		dtree->calcwt(0);
		dtree->puttree();
	    } else
		dtree->putwt();
	} else if (treemode.height && !treemode.tree) {
	    dtree->putheight(dtree->root);
	}
	if (treemode.pairwt) {
	    dtree->putpw();
	} 
	if (!treemode.weight && treemode.tree) {
	    if (!treemode.path) dtree->puttree();
	}
	if (treemode.ddmx) {
	    dtree->mkleaves();
	    dtree->estimatedist(dtree->root);
	    dtree->putdmx();
	}
	if (treemode.path) {
	    dtree->mkleaves();
	    dtree->tracepath(dtree->root);
	}
	if (tmp) {
	    gswap(msd, tmp);
	    delete tmp;
	}
	delete dtree;
}

template <class seq_t>
void AlnServer<seq_t>::setparam(int level)
{
	setphyl(QUERY);
	setdivseq(QUERY, QUERY, QUERY);
}

int main(int argc, const char** argv)
{
	pronam = argv[0];
	treemode.tree = true;
	AlnServer<Seq>	svr(argc, argv, IM_MULT, IM_MULT);
	if (svr.sql1->nextseq(svr.in_face) != IS_OK) usage();
	phyln(svr.in_face[0]);
	eraWlprms();
	EraDbsDt();
	return (0);
}

#else	// PHYL

#if NJ || UPG || WARD

void usage()
{
	char    str[NAMSIZ];

	fprintf(stderr, "Usage:\t%s [-options] [dist. matrix] [output]\n",
		partfnam(str, pronam, "b"));
	fputs("Options:\n", stderr);
	fputs("  -A[S]\t:Amino acid classification pattern\n", stderr);
	fputs("  -B[S]\t:Bit pattern for matches\n", stderr);
	fputs("  -b: \tPrint branch length.\t", stderr);
	fputs("  -cN:\tSet column width to \'N\'.\n", stderr);
	fputs("  -dS:\tSequence database name\n", stderr);
	fputs("  -D: \tPrint distance matrix.\t", stderr);
	fputs("  -e: \tSuppress negative edge.\t", stderr);
	fputs("  -hC:\tHorizontal line char = C.\n", stderr);
	fputs("  -i: \tInteractive mode.\t", stderr);
	fputs("  -kN:\tk-tuple size\n", stderr);
	fputs("  -m: \tSign is reversed.\n", stderr);
	fputs("  -n: \tPrint node #.\t", stderr);
	fputs("  -q[S]:\tRapid calculation sequence distance\n", stderr);
	fputs("  -r: \tRoot is put left.\n", stderr);
	fputs("  -T: \tTree topology in Newick format.\t", stderr);
	fputs("  -u: \tUnrooted tree.\n", stderr);
	fputs("  -vC:\tVertical line char = C.\t", stderr);
	fputs("  -w: \tPrint weight with tree.\n", stderr);
	fputs("  -y: \tWeight without tree.\t", stderr);
	fputs("  -?: \tThis.\n", stderr);
	fputs("\nFormat of dist. matrix: (Example)\n\n", stderr);
	fputs("Human\t\t: Name of Species 1\n", stderr);
	fputs("Gorilla\t\t: Name of Species 2\n", stderr);
	fputs("Chimpanzee\t: Name of Species 3\n", stderr);
	fputs("Orangutan\t: Name of Species 4\n", stderr);
	fputs("\t\t: A blank line\n", stderr);
	fputs("4\t\t: Dist[1-2]\n", stderr);
	fputs("3 5\t\t: Dist[1-3] Dist[2-3]\n", stderr);
	fputs("15 15 16\t: Dist[1-4] Dist[2-4] Dist[3-4]\n", stderr);
	exit(0);
}

int main(int argc, const char** argv)
{
const	char*   ap = 0;
const	char*   bp = 0;
const	char*   catalog = 0;
const	char*	seqpath = 0;
	int     k = 0;
	DistMes	dm = QJuCa;

	pronam = argv[0];
#if NJ
	treemode.method = NJ_METHOD;
#elif WARD
	treemode.method = WARD_METHOD;
#else
	treemode.method = UPG_METHOD;
#endif
	setprompt(false, false);
	treemode.tree = true;
	while(--argc > 0 && (++argv)[0][0] == '-') {
	  const char* val = argv[0] + 2;
	  switch(argv[0][1]) {
	    case 'b':
		if (*val == '-') treemode.branch = false; else
		if (*val == '+') treemode.branch = true; else
		treemode.branch = !treemode.branch;
		break;
	    case 'c':
		if (*val) clmwd = atoi(val);
		else if(--argc > 0)
		    clmwd = atoi((++argv)[0]);
		break;
	    case 'e': treemode.n_edge = true; break;
	    case 'f': 
		wfact = atof(getarg(argc, argv, true));
		if (wfact > 0.) cfact = false;
		break;
	    case 'h':
		treemode.tree = (*val == 't'); treemode.height = 1; break;
	    case 'i': setprompt(true, false); break;
	    case 'j': treemode.balance = true; break;
	    case 'k':       /* tuple size */
		k = atoi(getarg(argc, argv, true)); break;
	    case 'l': horline = *val? *val: ' '; break;
	    case 'm':       /* distance measure */
		dm = (DistMes) atoi(getarg(argc, argv, true));
		break;
	    case 'p':
		if (*val == 'p') OutPrm.taxoncode = 1;
		break;
	    case 'q':       /* calc dist by k-mer count */
		if (*val)
		    catalog = *val == ':'? val + 1: val;  /* seqs. in catalog */
		else if (argc > 2 && argv[1][0] != OPTCHAR && argv[2][0] == OPTCHAR)
		    {catalog = *++argv; argc--;}
		else
		    catalog = "";   /* seqs. in args */
		break;
	    case 'r': treemode.lroot = true; break;
	    case 's': 
		seqpath = getarg(argc, argv, false);
		if (seqpath) setdfn(seqpath);
		break;
	    case 't': thread_num = atoi(getarg(argc, argv, true)); break;
	    case 'u': treemode.unrooted = true; break;
	    case 'v': verline = *val? *val: ' '; break;
	    case 'w':
		treemode.tree = (*val == 't');
		if (*val == '2' || *val == 'p') treemode.pairwt = true;
		else	treemode.weight = true;
		break;
	    case 'A':       /* aa classification pattern */
		ap = getarg(argc, argv, false);
		break;
	    case 'B':       /* matching bit pattern */
		bp = getarg(argc, argv, true);
		if (!bp) bp = "";
		break;
	    case 'C': treemode.method = CLINK_METHOD; break;
	    case 'D': treemode.ddmx = true; break;
	    case 'E': treemode.method = AB_METHOD; break;
	    case 'K': treemode.prnode = true; break;
	    case 'N': treemode.method = NJ_METHOD; break;
	    case 'P': treemode.path = true; break;
	    case 'S': treemode.method = SLINK_METHOD; break;
	    case 'T': treemode.paren = treemode.branch = true; break;
	    case 'U': treemode.method = UPG_METHOD; break;
	    case 'W': treemode.method = WARD_METHOD; break;
	    case '-': treemode.negate = true; break;
	    case '?':
	    default:    usage();
	  }
	}

	FILE*   fo = stdout;
	FILE*   fd = stdin;
	DistMat*	dmat = 0;
	if (catalog) {
	    setqdiv(k, ap, bp, dm);
	    dmat = new DistMat(argc, argv, catalog);
	    if (!dmat->dist)  usage();
	} else {
	    if (argc > 0 && !(fd = fopen(*argv, "r")))
		fatal(not_found, *argv);
	    if (--argc > 0 && !(fo = fopen(*++argv, "w")))
		fatal("Can't create %s", *argv);
	    dmat = new DistMat(fd);
	    if (!dmat->dist) usage();
	    if (fd != stdin) fclose(fd);
	}
	PpPrm	dtree(fo, dmat, Def_METHOD, dmat->bias);
	if (treemode.pairwt) {
	    dtree.putpw();
	} else if (treemode.weight) {
	    if (treemode.tree) {
		dtree.calcwt(0);
		dtree.puttree();
	    } else
		dtree.putwt();
	} else if (treemode.height && !treemode.tree) {
	    dtree.putheight(dtree.root);
	} else {
	    if (!treemode.path)
		dtree.puttree();
	}
	if (treemode.ddmx) {
	    dtree.mkleaves();
	    dtree.estimatedist(dtree.root);
	    dtree.putdmx();
	}
	if (treemode.path) {
	    dtree.mkleaves();
	    dtree.tracepath(dtree.root);
	}
	delete dmat;
	eraStrPhrases();
	return (0);
}

#endif	// NJ || UPG || WARD
#endif	// PHYL
#endif	// MAIN
