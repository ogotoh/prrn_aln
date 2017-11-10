/*****************************************************************************
*
*	Alignment of two protein/nucleotide sequences 
*	5' and 3' splice site signals and coding potential
*	calculated from ditron frequnecies are considered.
*	Assumes there are internal gap(s) in the reference protein sequence.
*	The maximum number of candidates is limited to three.
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
*	Copyright(c) Osamu Gotoh <<o.gotoh@i.kyoto-u.ac.jp>>
*****************************************************************************/

#define MONITOR	1
#define DEBUG	1
#define TERMGOP	0

#include "aln.h"
#include "mseq.h"
#include "maln.h"
#include "gfreq.h"
#include "vmf.h"

#if MONITOR
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <sys/timeb.h>
#endif

static	const	int	NOLL = 3;
static	const	int	INTR = 2;

inline	bool	isEIJ(int phs) {return (algmode.lsg && phs > -2);}

struct DPUH {
	VTYPE	val;
	long	ptr;
	int	dir;
	int	jnc;
};

struct DPUH_hf : public DPUH {
	IDELTA*	dla;
	int	glb;
};

static	const	DPUH	blackDpUnit = {NEVSEL};
static	const	TraceBackDir	hori[4] = {HORI, HOR1, HOR2, HORI};
static	const	TraceBackDir	vert[4] = {VERT, SLA1, SLA2, VERT};

void reset(DPUH* rcd, bool all = true) {if (all) *rcd = blackDpUnit;}
void copy(DPUH* dst, DPUH* src) {*dst = *src;}

void reset(DPUH_hf* rcd, bool all = true) {
	if (all) reset((DPUH*) rcd);
	rcd->glb = 0;
	cleardelta(rcd->dla);
}

void copy(DPUH_hf* dst, DPUH_hf* src)
{
	copy((DPUH*) dst, (DPUH*) src);
	copydelta(dst->dla, src->dla);
	dst->glb = src->glb;
}

class Fwd2h3 {
	mSeq**	seqs;
	mSeq*	a;
	Seq*	b;
	WINDOW*	wdw;
	PwdM*	pwd;
	FTYPE	u2divu1;
	FTYPE	v2divv1;
	Vmf*	vmf;
	int	bbt;
	int	bbt3;
	int	recsize;
	Mfile*	mfd;
	int*	sigii;
	VTYPE	SpbFact;
	SpJunc*	spjcs;
	DPUH_hf*	buf;
	DPUH_hf*	ff[NOL - 1];	// F
	DPUH_hf*	hh[NOL + 1];	// H, G matrix +1 SPJ
	DPUH_hf*	hq;		// previous state
	DPUH_hf*	hl[NOLL][3];	// [DIA, HORI, HORL][phase][candidates]
	int*	glab;
	IDELTA*	gfab;
	void	initializeH();
	void	initH_hf();
	DPUH_hf*	lastH_hf();
	VTYPE	gapopen(DPUH_hf* rcd, SeqItr& asi, int d3);
	void	update(DPUH_hf* dst, DPUH_hf* src, SeqItr& asi, int d3);
public:
	Fwd2h3(mSeq** _seqs, PwdM* _pwd, WINDOW* _wdw);
	~Fwd2h3();
	VTYPE	forwardH_hf(long pp[]);
	SKL*	globalH_hf(VTYPE* scr);
};

static	VTYPE	SumCodePot(EXIN* bb, int i, CHAR* cs, PwdM* pwd);

void Fwd2h3::initializeH()
{
	int	apb = a->gfq->hetero + 1;
	gfab = new IDELTA[apb * recsize];
	IDELTA*	w = gfab - apb;
	DPUH_hf* r = buf;
	for (int n = 0; n < recsize; ++n, ++r) {
	    r->val = NEVSEL; r->ptr = r->dir = 0;
	    cleardelta(r->dla = w += apb);
	    r->glb = 0;
	}
	glab = 0;
}

VTYPE Fwd2h3::gapopen(DPUH_hf* rcd, SeqItr& asi, int d3)
{
	if (d3 == 0)
	    return pwd->newgap2(*asi.tfq, rcd->glb, rcd->dla);
	else if (d3 > 0)
	    return pwd->newgap1(*asi.sfq, rcd->dla, rcd->glb);
	else
	    return pwd->newgap2(*asi.rfq, rcd->glb, rcd->dla);
}

void Fwd2h3::update(DPUH_hf* dst, DPUH_hf* src, SeqItr& asi, int d3)
{
	if (d3 == 0) {
	    dst->dir = isdiag(src)? DIAG: NEWD;
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    dst->glb = 0;
	} else if (d3 > 0) {
	    dst->dir = VERT;
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    dst->glb = src->glb + 1;
	} else {
	    dst->dir = HORI;
	    incdelta(dst->dla, src->dla);
	    dst->glb = 0;
	}
	dst->ptr = src->ptr;
	dst->jnc = src->jnc;
}

Fwd2h3::Fwd2h3(mSeq* _seqs[], PwdM* _pwd, WINDOW* _wdw) :
	seqs(_seqs), wdw(_wdw), pwd(_pwd)
{
	a = seqs[0];
	b = (Seq*) seqs[1];
	bbt = b->many;
	bbt3 = 3 * bbt;
	sigii = pfq2loc(a);
	SpbFact = spb_fact();
	spjcs = new SpJunc(b);
	u2divu1 = pwd->BasicGEP < 0? pwd->LongGEP / pwd->BasicGEP: 0;
	v2divv1 = pwd->BasicGOP < 0? pwd->LongGOP / pwd->BasicGOP: 0;
	recsize = wdw->width * pwd->Nrow + 3 * pwd->Noll + NOLL * 3 * (INTR+1);
	buf = new DPUH_hf[recsize];
	DPUH_hf*	h = buf;
	for (int i = 0; i < pwd->Noll - 1; ++i, h += 3) ff[i] = h;
	for (int i = 0; i < pwd->Nrow; ++i, h += wdw->width)
	    hh[i] = h - wdw->lw + 3;
	hq = h; h += pwd->Noll;
	for (int i = 0; i < NOLL; ++i)
	    for (int j = 0; j < 3; ++j, h += INTR+1)
		hl[i][j] = h;
	initializeH();
	mfd = 0;
	vmf = 0;
}

Fwd2h3::~Fwd2h3()
{
	delete vmf;
	delete[] buf;
	delete[] glab;
	delete[] gfab;
	delete mfd;
	delete spjcs;
}

void Fwd2h3::initH_hf()
{
	int	m = a->left;
	int	n = b->left;
	int	r = n - 3 * m;
	int	rr = b->right - 3 * m;
	int	dir = a->inex.exgl? DEAD: DIAG;
	DPUH_hf*	h = hh[0];
	VTYPE	cand[4] = {0};

	if (wdw->up < rr) rr = wdw->up;
	SeqItr	asi(a, m - 1);
	GsqItr	bsi(b, n);
	for (int i = 0; r <= rr + 3; ++r, ++n, ++i, ++h, ++bsi) {
	    if (i == 0 || (a->inex.exgl && i < 3)) {
		h->dir = dir;
		h->val = (bsi.bb[1].sigS > 0)? bsi.bb[1].sigS: 0;
		h->ptr = vmf? vmf->add(m, n, 0): 0;
		h->jnc = n;
            } else if (a->inex.exgl) {	// semi-global
		cand[1] = h[-1].val + pwd->GapW1;
		cand[2] = h[-2].val + pwd->GapW2;
		int	k = n - h[-3].jnc;
		cand[3] = h[-3].val + pwd->GapExtPen3(k) + bsi.bb[-2].sigE;
		if ((algmode.lcl & 1)) cand[0] = bsi.bb[1].sigS; else
		if ((algmode.lcl & 4)) cand[0] = bsi.bb->sig3;
		k = vmax(cand, 4) - cand;
		h->val = cand[k];
		if (k) {		// extension
		    update(h, h - k, asi, -1);
		    h->dir = hori[k];
		} else {		// new start
		    reset(h, false);
		    h->ptr = vmf? vmf->add(m, n, 0L): 0L;
		    h->dir = dir;
		    h->jnc = n;
		}
	    }
	}

	r = b->left - 3 * a->left;
	rr = b->left - 3 * a->right;
	if (wdw->lw > rr) rr = wdw->lw;
	h = hh[0] + r;
	for (int i = 1; --r >= wdw->lw - 3; ++i) {
	    --h;
	    if (i < 3) {
		update(h, h + i, asi, 1);
		h->val += (i == 1? pwd->GapW1: pwd->GapW2);
		h->dir = VERT + i;
	    } else if (b->inex.exgl) {	// semi global
		h->val = 0;
		h->dir = DEAD;
		h->jnc = b->left + (i + 1) % 3;
		h->ptr = vmf? vmf->add(m, h->jnc, 0L): 0L;
	    } else {
		update(h, h + 3, asi, 1);
		h->val += (i == 3)? pwd->GapPenalty(3): pwd->GapExtPen3(i);
		h->dir = VERT;
	    }
	    if (i % 3 == 0) {
		 ++m; ++asi;
	    }
	}
}

DPUH_hf* Fwd2h3::lastH_hf()
{
	int	glen[3] = {0, 0, 0};
	int	rw = wdw->lw;
	int	m3 = 3 * a->right;
	int	rf = b->left - m3;
	if (rf > rw) rw = rf;
	else	rf = rw;
	int	r9 = b->right - m3;
	DPUH_hf*	h = hh[0] + rw;
	DPUH_hf*	h9 = hh[0] + r9;
	DPUH_hf*	f = h9;
	VTYPE	mxh = h9->val;
	GsqItr	bsi(b, rw + m3);

	if (a->inex.exgr) {
	    for (int p = 0; h <= h9; ++h, ++bsi, ++rf, ++p) {
		if (p == 3) p = 0;
		glen[p] += 3;
		VTYPE	cand[3] = {h->val, NEVSEL, NEVSEL};
		if (rf - rw >= 3 && h[-3].dir != DEAD) {
		    cand[1] = h[-3].val + bsi.bb[-2].sigE + pwd->GapExtPen3(glen[p]);
		    if ((algmode.lcl & 2) && !(h->dir & SPIN))
			cand[2] = h[-3].val + bsi.bb[-2].sigT;
		}
		int	k = vmax(cand, 3) - cand;
		if (k) {
		    copy(h, h - 3);
		    h->val = cand[k];
		} else if (!ishori(h)) glen[p] = 0;
		if (k == 2) {	// termination codon
		    h->dir = DEAD;
		    if (h->val > mxh) {
			f = h; mxh = h->val;
			h->ptr = vmf->add(a->right, rf + m3 - 3, h->ptr);
		    }
		} else {
		    if (k) h->dir = HORI;
		    if ((algmode.lcl & 16) && bsi.bb->sig5 > 0)
			cand[k] += bsi.bb->sig5;
		    if (cand[k] > mxh) {
			f = h; mxh = cand[k];
		    }
		}
	    }
	}
	f->val = mxh;
	if (b->inex.exgr) {
	    rw = wdw->up;
	    rf = b->right - 3 * a->left;
	    if (rf < rw) rw = rf;
	    h = hh[0] + rw;
	    for ( ; h > h9; --h, --rw) {
		VTYPE	x = h->val + (rw % 3? pwd->ExtraGOP: 0);
		if (f->val < x) {f = h; f->val = x;}
	    }
	}
	int p = f - h9;
	rf = a->right;	// m9
	rw = b->right;	// n9
	if (p > 0) {
	    rf -= (p + 2) / 3;
	    if (p %= 3) rw -= (3 - p);
	} else if (p < 0)	rw += p;
	if (vmf || a->inex.exgr) f->ptr = vmf->add(rf, rw, f->ptr);
	return (f);
}

VTYPE Fwd2h3::forwardH_hf(long* pp)
{
	if (!a->isprotein() || b->isprotein())
	    fatal("Inproper combination of seq. types !\n");
	bool	alcl = a->inex.exgl || a->inex.exgr;
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0L};
	bool	Local = algmode.lcl & 16;
	bool	LocalL = Local && a->inex.exgl && b->inex.exgl;
	bool	LocalR = Local && a->inex.exgr && b->inex.exgr;
	VTYPE	bthk[2] = {b->sumwt, 0};

	if (vmf) vmf->add(0, 0, 0L);	// Skip 0-th record
	initH_hf();

	int	m = a->left;
	if (!a->inex.exgl) --m;		// global
	int	n1 = 3 * m + wdw->lw - 1;
	int	n2 = 3 * m + wdw->up;
	VTYPE	pub = pwd->BasicGEP;
	SeqItr asip(a, m);
	for (SeqItr asi = asip; ++m <= a->right; asi = asip) {
	    ++asip;
	    bool	internal = !a->inex.exgr || m < a->right;
	    n1 += 3; n2 += 3;
	    int	n0 = max(n1, b->left);
	    int	n9 = min(n2, b->right);
	    int	n = n0;
	    int	r = n - 3 * m;
	    VTYPE	pua = alcl? 0: (pwd->*pwd->unpa)(asi.res, bthk);
	    DPUH_hf*	hf[NOLL];		// [DIAG, HORI, HORL]
	    DPUH_hf*	hg[NOLL + 1];		// [DIAG, VERT, VERL, SPJ]
	    int		nl[NOLL][3][INTR+1];	// [DIA, HORI, HORL][phase][candidates]
	    reset(hq);
	    for (int k = 0; k < pwd->Nrow; ++k) hg[k] = hh[k] + r;
	    for (int k = 0; k < pwd->Noll - 1; ++k) 
		for (int p = 0; p < 3; ++p) reset(ff[k] + p);
	    for (int k = 0; k < pwd->Noll; ++k) {
	      for (int p = 0; p < 3; ++p) {
		DPUH_hf*	phl = hl[k][p];
		for (int l = 0; l <= INTR; ++l, ++phl)
		  {nl[k][p][l] = l; reset (phl);}
	      }
	    }
#if DEBUG
	    if (algmode.nsa & 8) {
		printf("%2d %2d %2d", m, n, hg[0]->dir);
		putvar(hg[0]->val); putchar('\n');
	    }
#endif
	    GsqItr	bsi(b, n - 1);		// center of triplet
	    for (int q = 0; ++n <= n9; ++bsi) {
		VTYPE	cand[4];
		VTYPE&	x = cand[0];
		VTYPE&	y = cand[1];
		EXIN*	bb = bsi.bb;
		for (int k = 0; k < pwd->Nrow; ++k) ++hg[k];
		for (int k = 0; k < pwd->Noll - 1; ++k) hf[k+1] = ff[k] + q;

//	Diagonal
		VTYPE	sigE = bb->sigE;
		DPUH_hf*	src = hf[0] = hg[0];
		DPUH_hf*	mx = src;
		DPUH_hf*	from = src;
		DPUH_hf*	sj = hg[pwd->Noll];
		DPUH_hf*	tgt = hg[1];
		int		k = 0;
		if (m == a->left) goto Horizon;	// initialize column
		if (n > b->left + 2) {
		    copy(hq, src);	// reserve
		    if (sj->dir) {	// phase -1 3'
			src = sj;
			src->dir = DIAG;
			sj->dir = 0;
		    } else {
			src->val += (pwd->*pwd->Sim2)(asi.res, bsi.res) + sigE
				 + gapopen(src, asi, 0);
//			if (sndpro) src->val += sspsim(asi, bsi);
			src->dir = isdiag(src)? DIAG: NEWD;
		    }
		    update(mx, src, asi, 0);
		} else reset(mx);

//	vertical gap extention
		cand[0] = tgt[3].val + gapopen(tgt + 3, asi, 1);
//	1 nt deletion
		++from;
		cand[1] = from->val + (isvert(from)? pwd->GapE1: pwd->GapW1);
//	2 nt deletion
		++from;
		cand[2] = from->val + (isvert(from)? pwd->GapE2: pwd->GapW2);
//	normal deletion
		++from;
		cand[3] = from->val + gapopen(from, asi, 1);
		k = vmax(cand, 4) - cand;
		tgt->val = cand[k];
		tgt->dir = vert[k];
		src = k? src + k: tgt + 3;
		if (alcl) pua = (pwd->*pwd->unpa)(asi.res, bthk);
		tgt->val += pua;
		update(tgt, src, asi, 1);
		if (vmf && tgt->dir & SPIN)
		    tgt->ptr = vmf->add(m - 1, n, tgt->ptr);
		if (tgt->val > mx->val) mx = tgt;

//	long deletion
		if (pwd->Noll == 3) {
		    tgt = hg[2];
		    src = tgt + 3;
		    cand[3] = from->val + v2divv1 * gapopen(from, asi, 1);
		    if (cand[3] > src->val) {
			src = from;
			tgt->val = cand[3];
			tgt->dir = VERL;
		    } else {
			tgt->val = src->val;
		    }
		    tgt->val += u2divu1 * pua;
		    update(tgt, src, asi, 1);
		    if (vmf && tgt->dir & SPIN)
			tgt->ptr = vmf->add(m - 1, n, tgt->ptr);
		    if (tgt->val > mx->val) mx = tgt;
		}
Horizon:
		from = hg[0] - 3;
//	long insertion
		if (n < b->left + 3) continue;
		if (pwd->Noll == 3) {
		    tgt = src = hf[2];
		    copy(hq + 2, tgt);
		    x = from->val + v2divv1 * gapopen(from, asi, -1);
		    if (x > src->val) {
			src = from;
			tgt->val = x;
		    }
		    tgt->val += pwd->LongGEP;
		    if (!(src->dir & SPF2)) tgt->val += sigE;
		    tgt->dir = (src->dir & SPIN) + HORL;
		    update(tgt, src, asi, -1);
		    if (tgt->val > mx->val) mx = tgt;
		}

		src = from;
		tgt = hf[1];
		copy(hq + 1, tgt);
//	nomal extension
		cand[3] = tgt->val;
//	nomal insertion
		cand[0] = from->val + gapopen(from, asi, -1);
//	2 nt insertion
		cand[1] = (++from)->val + pwd->GapW2;
//	1 nt insertion
		cand[2] = (++from)->val + pwd->GapW1;
		k = vmax(cand, 4) - cand;
		src = k? src + k: tgt;
		tgt->val = cand[k];
		tgt->dir = hori[k]+ (src->dir & SPIN);
		tgt->val += pwd->BasicGEP;
		if (!(src->dir & SPF2)) tgt->val += sigE;
		update(tgt, src, asi, -1);
		if (tgt->val > mx->val) mx = tgt;
		if (++q == 3) q = 0;

//	3' boundary, assume no overlapping signals
		bb += 2;
		if (isEIJ(bb->phs3) && internal) {
		    int phs = (bb->phs3 == 2)? 1: bb->phs3;
Acceptor:
		    if (phs == -1)
			sj->val = mx->val + gapopen(mx, asi, 0)
			    + (pwd->*pwd->Sim2)(asip.res, bsi.res + 2 * bbt);
		    int	nb = n - phs;
		    VTYPE sigJ = sigii? SpbFact * sigii[3 * m - phs]: 0;
		    DPUH_hf* maxdhl = 0;
		    for (k = 0; k < pwd->Noll; ++k) {
			DPUH_hf*	maxphl = 0;
			int*	pnl = nl[k][phs + 1];
			tgt = hf[k];
			for (int l = 0; l < INTR; ++l) {
			    DPUH_hf*	phl = hl[k][phs + 1] + pnl[l];
			    if (!phl->dir) break;
			    y = phl->val + sigJ +
                                pwd->IntPen->Penalty(nb - phl->jnc) +
				b->exin->sig53(phl->jnc, nb, IE53);
			    if (phs == 1) {
				CHAR*	cs = spjcs->spjseq(phl->jnc, nb);
				y += pwd->pmt->prematT(cs);
				if (isdiag(phl)) y += (pwd->*pwd->Sim2)(asi.res, cs);
			    } else if (phs == -1) {
				CHAR*	cs = spjcs->spjseq(phl->jnc, nb) + bbt;
				y += pwd->pmt->prematT(cs);
				x = y + (pwd->*pwd->Sim2)(asip.res, cs)
				      + gapopen(phl, asip, 0);
				if (x > sj->val) {
				    maxdhl = phl;
 				    sj->val = x;
				}
			    }
			    if (y > tgt->val) {
 				tgt->val = y;
 				maxphl = phl;
			    }
			}
 			if (maxphl) {
			    copy(tgt, maxphl);
		 	    if (vmf) tgt->ptr = vmf->add(m, n, 
				vmf->add(m, maxphl->jnc + phs, maxphl->ptr));
			    tgt->jnc = nb;
			    tgt->dir = maxphl->dir | SPJCI;
			    if (phs == -1) tgt->dir |= SPF2;
			    if (isvert(maxphl)) copy(hg[1], tgt);
			    if (maxphl->dir == VERL) copy(hg[2], tgt);
			    if (tgt->val > mx->val) mx = tgt;
			}
		    }
		    if (maxdhl) {
			copy(sj, maxdhl);
			sj->jnc = nb;
			sj->dir = maxdhl->dir | SPALL;
		 	if (vmf) sj->ptr = vmf->add(m, n, 
			    vmf->add(m, maxdhl->jnc + phs, maxdhl->ptr));
		    }
		    if (bb->phs3 - phs == 1) {	// AGAG
			phs = -1;
			goto Acceptor;
		    }
		}

		if (vmf && mx->dir == NEWD)
		    mx->ptr = vmf->add(m - 1, n - 3, mx->ptr);

//	Find optimal path
		tgt = hf[0];
		y = tgt->val;
		if (mx != tgt) copy(tgt, mx);
		else if (Local && y > hq->val) {
		    if (LocalL && hq->dir == 0 && !(tgt->dir & SPJC))
			tgt->ptr = vmf? vmf->add(m - 1, n - 3, 0): 0;
		    else if (LocalR && y > maxh.val) {
			maxh.val = y;
			maxh.p = tgt->ptr;
			maxh.m = m;
			maxh.n = n;
		    }
		}
		if (LocalL && tgt->val <= 0) reset(tgt);
		else if (vmf && tgt->dir == NEWD)
		    tgt->ptr = vmf->add(m - 1, n - 3, tgt->ptr);

//	5' boundary
		if (isEIJ(bb->phs5) && internal) { 
		    int	phs = (bb->phs5 == 2)? 1: bb->phs5;
Donor:
		    int	nb = n - phs;
		    VTYPE	sigJ = b->exin->sig53(nb, 0, IE5);
		    for (int k = 0; k < pwd->Noll; ++k) {
 			src = (phs == 1)? hq + k: hf[k];
                               // An orphan exon is disallowed
			if (!src->dir || isvert(src) || (src->dir & SPIN))
			    continue;
			if (phs == 1 && nb - from->jnc == 2 && bb[-2].sigS > 0)
			    from->val -= bb[-2].sigS;
			x = src->val + sigJ;
			DPUH_hf*	phl = hl[k][phs + 1];
			int*	pnl = nl[k][phs + 1];
			phl[pnl[INTR]].dir = 0; // put x to the last
			int l = INTR;
			while (--l >= 0) {
			    if (x > phl[pnl[l]].val)
				gswap(pnl[l], pnl[l + 1]);
			    else
				break;
			}
			if (++l < INTR) {
			    phl[pnl[INTR]].dir = 0;
			    phl += pnl[l];
			    copy(phl, src);
			    phl->val = x;
			    phl->jnc = nb;
			}
		    }
		    if (bb->phs5 - phs == 1) {	//	GTGT..
			phs = -1;
			goto Donor;
		    }
		}
#if DEBUG
	if (algmode.nsa & 8) {
	    printf("%2d %2d %2d ", m, n, mx->dir);
	    putvar(mx->val); putvar(y);
	    for (int i = 1; i < pwd->Noll; ++i) {
		putvar(hg[i]->val); putvar(hf[i]->val);
	    }
	    if (algmode.lsg) {
		putvar(hl[0][2][nl[0][2][0]].val);
		putvar(hl[0][1][nl[0][1][0]].val);
		putvar(hl[0][0][nl[0][0][0]].val);
  		printf(" %2d %2d %5.2f", bb->phs5, bb->phs3, (double) sigE);
		printf(" %5.1f %5.1f", (double) bb->sig5, (double) bb->sig3);
	    }
	    putchar('\n');
	}
#endif
	    }	// end of n-loop
	}	// end of m-loop
	DPUH_hf*	tgt = lastH_hf();
	pub = tgt->val;
	*pp = tgt->ptr;
	return (pub);
}

static VTYPE SumCodePot(EXIN* bb, int i, CHAR* cs, PwdM* pwd)
{
	VTYPE	y = 0;

	for (++bb; i > 0; i -= 3, bb += 3) {
	    VTYPE	hvl;
	    if (cs) {
		hvl = pwd->pmt->prematT(cs);
		cs = 0;
	    } else
		hvl = bb->sigE;
	    y += hvl;
	}
	return (y);
}

VTYPE skl_rngH_hf(mSeq* seqs[], Gsinfo* gsi, PwdM* pwd)
{
	mSeq*	a = seqs[0];
	Seq*	b = (Seq*) seqs[1];
	int	abt = a->byte;
	int	bbt = b->byte;
	SKL*	wsk = gsi->skl;
	int	num = (wsk++)->n;
	int*	sigii = pfq2loc(a);
	VTYPE	SpbFact = spb_fact();
	SpJunc*	spjcs = new SpJunc(b);
	VTYPE	h = 0;
	VTYPE	ha = 0;
	VTYPE	hb = 0;
	VTYPE	hvl = 0;
	VTYPE	xi = 0;
	VTYPE	xc = 0;
	VTYPE	xc2 = NEVSEL;
	VTYPE	gop = 0;
	VTYPE	sig5 = 0;
	VTYPE	sig3 = 0;
	VTYPE	wgop = pwd->wgop(1);
	int	insert = 0;
	int	deletn = 0;
	int	intlen = 0;
	int	phs = 0;
	int	nb = 0;
	EXIN*	bbi = 0;
	EISCR	rbuf;
	FSTAT*  fst = &gsi->fstat;
	FSTAT   pst;
	Eijnc*	eijnc = gsi->eijnc = new Eijnc();
static	const char*	fswarn = "%s %s FrameShift: %d %d\n";

	vclear(fst);	// clear fst
	vclear(&pst);	// clear pst
	vclear(&rbuf);	// clear rbuf
	CHAR*	cs = b->at(wsk[num-1].n - 2);
	CHAR*	cs2 = 0;
	bool	termcodon = (*cs == TRM || *cs == TRM2);
	cs = 0;
	if (termcodon) wsk[num-1].n -= 3;
	if (wsk[1].n == wsk->n && b->inex.exgl) {++wsk; --num;}
	IDELTA*	dla = new IDELTA[a->gfq->hetero + 3];
	int	glb = 0;
	cleardelta(dla);
	int	m = wsk->m;
	int	n = wsk->n;
	int	i2 = 0;
	SeqItr	asi(a, m);
	GsqItr	bsi(b, n);
	GsqItr	brim(b, b->right - 1);

	if ((algmode.lcl & 1) && bsi.bb[1].sigS > 0) h = bsi.bb[1].sigS;
	if ((algmode.lcl & 4) && bsi.bb->sig3 > h) h = bsi.bb->sig3;
	rbuf.left = n;
	rbuf.rleft = m;
	rbuf.phs = phs;
	rbuf.iscr = NEVSEL;
	rbuf.sig3 = h;
	while (--num) {
	    ++wsk;
	    int	mi = (wsk->m - m) * 3;
	    if (mi && insert) {
		if (intlen) {
		    xi = rbuf.iscr + xc;
		    xi += pwd->GapPenalty3(insert - intlen, gop);
		} else	xi = NEVSEL;
		VTYPE	x = xc + pwd->GapPenalty3(insert, gop);
#if !TERMGOP
		if (a->inex.exgl && m == a->left) x -= pwd->BasicGOP;
#endif
		if (xi >= x) {			// intron
		    hb = ha;
		    insert -= intlen;
		    if (rbuf.right - rbuf.left > 1) eijnc->push(&rbuf);
		    rbuf.left = rbuf.right + intlen;
		    rbuf.rleft = m;
		    rbuf.sig3 = sig3;
		    h += xi;
		} else {			// gap
		    h += x;
		}
		fst->unp += insert;
		if (insert) fst->gap += gop / wgop;
		if ((phs = insert % 3)) {	// insertion frame shift
		    rbuf.right = n - phs;
		    rbuf.rright = m;
		    rbuf.iscr = NEVSEL;
		    prompt(fswarn, a->sqname(), b->sqname(), b->SiteNo(rbuf.right), phs);
		    eijnc->push(&rbuf);
		    rbuf.left = n;
		    rbuf.rleft = m;
                    h += (phs == 1? pwd->GapE1: pwd->GapE2);
		}
		insert = intlen = i2 = 0;
		gop = xc = 0;
		xc2 = rbuf.iscr = NEVSEL;
	    }
	    int	ni = wsk->n - n;
	    if (ni && deletn && mi) {
		if ((phs = deletn % 3)) {	// deletion frame shift
		    rbuf.right = n + phs;
		    rbuf.rright = m;
		    rbuf.iscr = NEVSEL;
		    prompt(fswarn, a->sqname(), b->sqname(), b->SiteNo(rbuf.right), -phs);
		    eijnc->push(&rbuf);
		    rbuf.left = n;
		    rbuf.rleft = m;
		    h += pwd->ExtraGOP;
		    ++asi;
		    phs = 3 - phs;
		    bsi += phs;
		}
		deletn = 0;
	    }
	    int	i = mi - ni;
	    int	d = (i >= 0)? ni: mi;
	    n += d;
	    m += d / 3;
	    VTYPE	x = 0;
	    for ( ; d > 2; d -= 3, ++asi, bsi += 3) {
		if (cs) {
		    hvl = (pwd->*pwd->Sim2)(asi.res, cs + bbt) + 
				pwd->pmt->prematT(cs + bbt);
		    (pwd->*pwd->stt2)(fst, asi.res, cs + bbt, 0);
		    cs = 0;
		} else {
		    hvl = (pwd->*pwd->Sim2)(asi.res, bsi.res + bbt) + bsi.bb[1].sigE;
		    (pwd->*pwd->stt2)(fst, asi.res, bsi.res, 0);
		}
		x += hvl += (bsi.dullend()? 0: pwd->newgap2(*asi.tfq, glb, dla));
		newdelta(dla, *asi.tfq, dla);
		glb = 0;
	    }
	    h += x;
	    if (i > 0) {
		deletn += i;
		for ( ; i > 2; i -= 3, ++asi) {
		    h += (pwd->*pwd->unpa)(asi.res, &b->sumwt) + 
			pwd->newgap1(*asi.sfq, dla, glb);
		    newdelta(dla, *asi.tfq, dla);
		    ++glb;
		}
	    } else if (i < 0 && wsk->n < b->right) {
		if (!insert) bbi = bsi.bb;
		EXIN*	b3 = bsi.bb + (i = -i);
		int	n3 = n + i;
		int	phs5 = (bsi.bb->phs5 == 2)? b3->phs3: bsi.bb->phs5;
		int	phs3  = (b3->phs3 == 2)? bsi.bb->phs5: b3->phs3;
		 	// potential intron
		x = xi = NEVSEL;
		if (phs3 == 2 && phs5 == 2) {
			// phs3 = phs5 = -1; 
		    nb = n + 1; n3 = nb + i;
		    if (isJunct(phs3, phs5)) {
			x = b->exin->sig53(nb, n3, IE5P3);
		        if (sigii) x += SpbFact * sigii[3 * m + 1];
		    }
		    phs3 = phs5 = 1;
		}
		nb = n - phs3; n3 = nb + i;
		if (isJunct(phs3, phs5)) {
		    xi = b->exin->sig53(nb, n3, IE5P3);
		    if (sigii) x += SpbFact * sigii[3 * m - phs3];
		    if (phs3 != 0) cs = spjcs->spjseq(nb, n3);
		    if (insert == 0 && phs3 == 1)
			xi += (pwd->*pwd->Sim2)(asi.res - abt, cs) - hvl 
			   + pwd->pmt->prematT(cs);
		    if (x > xi) {
			xi = x;
			phs3 = -1;
			nb = n - phs3; n3 = nb + i;
			cs = spjcs->spjseq(nb, n3);
		    }
		    xi += pwd->IntPen->Penalty(i);
		    if (phs3 != -1) cs = 0;
		}
		x = SumCodePot(bsi.bb, i, cs? cs + bbt: 0, pwd);
		if ((xi > x + pwd->GapPenalty3(i) && xi > rbuf.iscr)) {
		    if (i2) {
			gop += pwd->newgap2(asi.rfq[-1], glb, dla);
			incdelta(dla, dla, i2 / 3);
			glb = 0;
			xc = xc2;
		    }
		    i2 = i;
		    xc2 = xc + x;
		    cs2 = cs;
		    intlen = i;		// intron
		    rbuf.right = nb;
		    rbuf.rright = m;
		    rbuf.phs = phs3;
		    rbuf.iscr = xi;
		    rbuf.escr = h + xc + pwd->GapPenalty3(insert, gop);
		    rbuf.sig5 = sig5;
		    ha = rbuf.escr + xi - b3->sig3;
		    rbuf.escr += sig5 - hb;
		    rbuf.mch = (int) (fst->mch - pst.mch);
		    rbuf.mmc = (int) (fst->mmc - pst.mmc);
		    rbuf.gap = (int) (fst->gap - pst.gap);
		    rbuf.unp = (int) (fst->unp - pst.unp) / 3;
                    pst = *fst;
		} else {		// gap
		    gop += pwd->newgap2(asi.rfq[-1], glb, dla);
		    incdelta(dla, dla, i / 3);
		    glb = 0;
		    if (i2 && cs != cs2)
			x = SumCodePot(bsi.bb, i, cs2? cs2 + bbt: 0, pwd);
		    xc += x;
		    cs = 0;
		}
		bsi += i;
		insert += i;
	    }
	    m = wsk->m;
	    n = wsk->n;
	}
	xc2 = 0;
	if (n >= 2) {
	    if (insert) {
		h += xc + pwd->GapPenalty3(insert, gop);
		fst->gap += gop / wgop;
#if !TERMGOP
		if (a->inex.exgr && m == a->right && insert) h -= pwd->BasicGOP;
#endif
		if ((phs = insert % 3)) {	// insertion frame shift
		    rbuf.right = n - phs;
		    rbuf.rright = m;
		    prompt(fswarn, a->sqname(), b->sqname(), b->SiteNo(rbuf.right), phs);
		    eijnc->push(&rbuf);
		    rbuf.left = n;
		    rbuf.rleft = m;
		}
	    }
	    if (deletn && !(b->inex.exgr && n == b->right)) {
	        if ((phs = deletn % 3)) {	// deletion frame shift
		    rbuf.right = n + phs;
		    rbuf.rright = m;
		    prompt(fswarn, a->sqname(), b->sqname(), b->SiteNo(rbuf.right), -phs);
		    eijnc->push(&rbuf);
		    rbuf.left = n;
		    rbuf.rleft = m;
		    h += pwd->ExtraGOP;
		}
	    }
	    if (bsi < brim) {
		VTYPE	x = 0;
		if (algmode.lcl & 2 && bsi.bb[1].sigT > 0)
		    x = xc2 = bsi.bb[1].sigT;
		if (algmode.lcl & 8 && bsi.bb->sig5 > 0 
			&& bsi.bb[1].sig5 > bsi.bb[1].sigT)
		    x = xc2 = bsi.bb->sig5;
		h += x;
	    }
	}
	delete[] dla;
	rbuf.escr = h - hb;
	rbuf.iscr = 0;
	rbuf.sig5 = xc2;
	rbuf.right = n;
	rbuf.rright = m;
	rbuf.mch = (int) (fst->mch - pst.mch);
	rbuf.mmc = (int) (fst->mmc - pst.mmc);
	rbuf.gap = (int) (fst->gap - pst.gap);
	rbuf.unp = (int) (fst->unp - pst.unp) / 3;
	eijnc->push(&rbuf);
	rbuf.left = endrng.left;
	rbuf.right = endrng.right;
	eijnc->push(&rbuf);
	eijnc->flush();
	gsi->noeij = eijnc->size() - 1;
	delete spjcs;
	fst->mch /= a->sumwt;
	fst->mmc /= a->sumwt;
	fst->unp /= 3;
	fst->val = (FTYPE) h;
	if (termcodon) wsk->n += 3;
	return (h);
}

SKL* Fwd2h3::globalH_hf(VTYPE* scr)
{
	mfd = new Mfile(sizeof(SKL));
	SKL	wsk;
	mfd->write((UPTR) &wsk);	// dummy call
	vmf = new Vmf();
	long	ptr;
	*scr = forwardH_hf(&ptr);
	if (ptr) {
	    SKL*	lskl = vmf->traceback(ptr);
	    SKL*	lwsk = lskl;
	    while (lskl->n--) mfd->write((UPTR) ++lwsk);
	    delete[] lskl;
	}
	wsk.n = (int) mfd->size();
	SKL*	skl = (SKL*) mfd->flush();
	skl->n = wsk.n - 1;
	skl->m = 1;
	if (skl->n == 0) {
	    delete[] skl;
	    return (0);
	}
	return stdskl3(&skl);
}

VTYPE HomScoreH_hf(mSeq* seqs[], PwdM* pwd)
{
	WINDOW	wdw;

	stripe31((Seq**) seqs, &wdw, alprm.sh);
	Fwd2h3 alnv(seqs, pwd, &wdw);
	return alnv.forwardH_hf(0);
}

SKL* alignH_hf(mSeq* seqs[], PwdM* pwd, VTYPE* scr)
{
	WINDOW  wdw;

	stripe31((Seq**) seqs, &wdw, alprm.sh);	// Recacl. window boundaries
	Fwd2h3 alnv(seqs, pwd, &wdw);
	return (alnv.globalH_hf(scr));
}

