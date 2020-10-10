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
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#ifndef  _FWD2H_
#define  _FWD2H_

#define MONITOR	1
#define F2DEBUG	1
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

static	const	int	NCAND_H = 4;
static	const	int	NQUE = 3;
static	const	TraceBackDir	hori[4] = {HORI, HOR1, HOR2, HORI};
static	const	TraceBackDir	vert[4] = {VERT, SLA1, SLA2, VERT};

template <class recd_t>
class Fwd2h {
	mSeq**	seqs;
	mSeq*&	a;
	mSeq*&	b;
	mSeqItr	asi;
	mSeqItr	bsi;
	mSeqItr	bzsi;
const	WINDOW*	wdw;
const	PwdM*	pwd;
	Vmf*	vmf;
	Mfile*	mfd;
	IDELTA*	gfab;
	int	abt;
	int	bbt;
	int	bbt3;
	FTYPE	u2divu1;
	FTYPE	v2divv1;
	int	Nod;
	SpJunc*	spjcs;
	recd_t*	recbuf;
	recd_t*	hh[NOL + 1];	// H, G matrix +1 SPJ
	recd_t*	hl[3];		// [phase][candidates]
	recd_t*	e1, *e2;
	recd_t*	hq;
	void	initializeH(int recsize);
	void	initH();
	recd_t*	lastH();
	VTYPE	gapopen(const recd_t* rcd, const mSeqItr& asi, int d3);
	void	update(recd_t* rcd, const mSeqItr& asi, int d3);
	void	update(recd_t* dst, const recd_t* src, const mSeqItr& asi, 
		VTYPE gop, int d3);
	mSeqItr	blank_Seqitr;
public:
	Fwd2h(mSeq** _seqs, const PwdM* _pwd, const WINDOW* _wdw = 0);
	~Fwd2h();
	VTYPE	forwardH(long pp[]);
	SKL*	globalH(VTYPE* scr);
	VTYPE	verify(Gsinfo* gsi);
};

template <class recd_t>
Fwd2h<recd_t>::Fwd2h(mSeq* _seqs[], const PwdM* _pwd, const WINDOW* _wdw) :
	seqs(_seqs), a(seqs[0]), b(seqs[1]),
	asi(a, a->left), bsi(b, b->left),
	wdw(_wdw), pwd(_pwd), vmf(0), mfd(0), gfab(0),
	abt(a->byte), bbt(b->many), bbt3(3 * bbt),
	u2divu1(pwd->BasicGEP < 0? (FTYPE) pwd->LongGEP / pwd->BasicGEP: 0),
	v2divv1(pwd->BasicGOP < 0? (FTYPE) pwd->LongGOP / pwd->BasicGOP: 0),
	Nod(2 * pwd->Noll - 1)
{
	SeqThk	bthk = {b->sumwt, 0, b->sumwt};
	a->setsimmtx(pwd->simmtx);
	bzsi.dns = &bthk;
	spjcs = new SpJunc(b);
	size_t	recsize = wdw? wdw->width * pwd->Nrow + 
		3 * (NCAND_H + 1) + NQUE * (pwd->Noll - 1): 3;
	recbuf = new recd_t[recsize + 1];
	initializeH(recsize + 1);
	if (!wdw) return;
	recd_t*	h = recbuf - wdw->lw + 3;
	for (int i = 0; i < pwd->Nrow; ++i, h += wdw->width) hh[i] = h;
	h += (wdw->lw - 3);
	for (int i = 0; i < 3; ++i, h += (NCAND_H + 1)) hl[i] = h;
	e1 = h; h += NQUE;
	if (pwd->Noll == 3) {e2 = h; h += NQUE;}
	hq = h;
}

template <class recd_t>
Fwd2h<recd_t>::~Fwd2h()
{
	delete vmf;
	delete[] recbuf;
	delete[] gfab;
	delete mfd;
	delete spjcs;
}

template <class recd_t>
void Fwd2h<recd_t>::initH()
{
	int	m = a->left;
	int	n = b->left;
	int	r = n - 3 * m;
	int	rr = b->right - 3 * m;
	int	dir = a->inex.exgl? DEAD: DIAG;
	recd_t*	h = hh[0] + r;
	VTYPE	cand[4];
	VTYPE&	x = cand[0];

	if (wdw->up < rr) rr = wdw->up;
	asi.reset(m - 1);
	bsi.reset(n);
	for (int i = 0; r <= rr; ++r, ++n, ++i, ++h, ++bsi) {
	    if (i == 0 || (a->inex.exgl && i < 3)) {
		h->dir = dir;
		h->val = (bsi.bb[1].sigS > 0)? bsi.bb[1].sigS: 0;
		h->ptr = vmf? vmf->add(m, n, 0): 0;
		h->jnc = n;
            } else if (a->inex.exgl) {	// semi-global
		cand[1] = h[-1].val + pwd->GapW1;
		cand[2] = h[-2].val + pwd->GapW2;
		int	k = n - h[-3].jnc;
		cand[3] = h[-3].val + pwd->TermGapExtPen3(k) + bsi.bb[-2].sigE;
		x = 0;
		if ((algmode.lcl & 1) && bsi.bb[1].sigS > x) x = bsi.bb[1].sigS;
		if ((algmode.lcl & 4) && bsi.bb->sig3 > x) x = bsi.bb->sig3;
		k = vmax(cand, 4) - cand;
		recd_t*	src = h - k;
		if (k) {		// extension
		    update(h, src, asi, cand[k] - src->val, -k);
		    h->dir = hori[k];
		} else {		// new start
		    h->reseth();
		    h->val = x;
		    h->ptr = vmf? vmf->add(m, n, 0L): 0L;
		    h->dir = dir;
		    h->jnc = n;
		}
	    }
	}

	r = b->left - 3 * a->left;
	rr = b->left - 3 * a->right;
	if (wdw->lw> rr) rr = wdw->lw;
	h = hh[0] + r;
	for (int i = 1; --r >= wdw->lw; ++i) {
	    --h;
	    if (b->inex.exgl) {	// semi global
		h->val = 0;
		h->dir = DEAD;
		h->jnc = b->left + i % 3;
		h->ptr = vmf? vmf->add(m, h->jnc, 0L): 0L;
	    } else if (i < 3) {
		update(h, h + i, asi, i == 1? pwd->GapW1: pwd->GapW2, i);
		h->dir = VERT + i;
	    } else {
	  	VTYPE	pua = (pwd->*pwd->unpa)(asi, bzsi);
	  	VTYPE	gnp = gapopen(h + 3, asi, 3);
	 	gnp = (m - a->left < pwd->codonk1)? gnp + pua:
		    (VTYPE) (v2divv1 * gnp + u2divu1 * pua);
		update(h, h + 3, asi, gnp, 3);
		h->dir = VERT;
	    }
	    if (i % 3 == 0) {
		 ++m; ++asi;
	    }
	}
}

template <class recd_t>
recd_t* Fwd2h<recd_t>::lastH()
{
	int	glen[3] = {0, 0, 0};
	int	rw = wdw->lw;
	int	m3 = 3 * a->right;
	int	rf = b->left - m3;
	if (rf > rw) rw = rf;
	else	rf = rw;
	int	r9 = b->right - m3;
	recd_t*	h = hh[0] + rw;
	recd_t*	h9 = hh[0] + r9;
	recd_t*	f = h9;
	VTYPE	mxh = h9->val;
	bsi.reset(rw + m3);

	if (a->inex.exgr) {
	    for (int p = 0; h <= h9; ++h, ++bsi, ++rf, ++p) {
		if (p == 3) p = 0;
		glen[p] += 3;
		VTYPE	cand[3] = {h->val, NEVSEL, NEVSEL};
		if (rf - rw >= 3 && h[-3].dir != DEAD) {
		    cand[1] = h[-3].val + bsi.bb[-2].sigE + pwd->TermGapExtPen3(glen[p]);
		    if ((algmode.lcl & 2) && !(h->dir & SPIN))
			cand[2] = h[-3].val + bsi.bb[-2].sigT;
		}
		int	k = vmax(cand, 3) - cand;
		if (k) {
		    *h = h[-3];
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
	rw = b->right; 	// n9
	if (p > 0) {
	    rf -= (p + 2) / 3;
	    if (p %= 3) rw -= (3 - p);
	} else if (p < 0)	rw += p;
	if (vmf || a->inex.exgr) f->ptr = vmf->add(rf, rw, f->ptr);
	return (f);
}

template <class recd_t>
VTYPE Fwd2h<recd_t>::forwardH(long* pp)
{
	if (!a->isprotein() || b->isprotein())
	    fatal("Inproper combination of seq. types !\n");
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
	recd_t*	hf[NOD];		// [DIAG, HORI, VERT, HORL, VERTL]
	int	nx[3][NCAND_H + 1];
	recd_t*	maxphl[NOD + 1];	// +1 for sj
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0L};
	bool	Local = algmode.lcl & 16;
	bool	LocalL = Local && a->inex.exgl && b->inex.exgl;
	bool	LocalR = Local && a->inex.exgr && b->inex.exgr;
	mSeqItr	csi(b);

	if (vmf) vmf->add(0, 0, 0L);	// Skip 0-th record
	initH();

	int	m = a->left;
	PfqItr	api(a, m);
	int	api_size = api.size();
	if (!a->inex.exgl) --m;		// global
	int	n1 = 3 * m + wdw->lw - 1;
	int	n2 = 3 * m + wdw->up;
	mSeqItr apsi(a, m);
	for (mSeqItr asi = apsi; ++m <= a->right; asi = apsi) {
	    ++apsi;
	    bool	internal = !a->inex.exgr || m < a->right;
	    n1 += 3; n2 += 3;
	    int	n0 = max(n1, b->left);
	    int	n9 = min(n2, b->right);
	    int	n = n0;
	    int	r = n - 3 * m;
	    int		k = 0;
	    recd_t*&	h = hf[0] = hh[k++] + r;
	    recd_t*&	g = hf[2] = hh[k++] + r;
	    recd_t*&	g2 = hf[4] = dagp? hh[k++] + r: 0;
	    recd_t*	sj = hh[k] + r;
const	    VTYPE*	qprof = asi.prof();	// sim2(asi, .)
const	    VTYPE*	qprof1 = apsi.prof();	// sim2(asi + 1, .)
	    for (int p = 0; p < NQUE; ++p) {
		e1[p].reseth();
		if (dagp) e2[p].reseth();
	    }
	    if (!b->inex.exgl && m == a->left) {
		e1[2] = hh[0][r];
		e1[2].val = pwd->GapW3;
		if (dagp) {
		    e2[2] = hh[0][r];
		    e2[2].val = pwd->GapW3L;
		}
	    }
	    for (int p = 0; p < 3; ++p) {
		for (int l = 0; l <= NCAND_H; ++l) {
		    hl[p][l].reseth();
		    nx[p][l] = l;
	        }
	    }
	    VTYPE	pua = internal? (pwd->*pwd->unpa)(asi, bzsi): 0;
	    int	ncand[3] = {0, 0, 0};
#if F2DEBUG
	    if (OutPrm.debug) {
		printf("%2d %2d %2d", m, n, h->dir);
		putvar(h->val); putchar('\n');
	    }
#endif
	    bool	a_in_zone = api_size && api.eq(m);
	    bsi.reset(n - 1);			// center of triplet
	    mSeqItr	bsip2(b, n + 2);	// center of triplet
	    for (int q = 0; ++n <= n9; ++bsi, ++bsip2) {
		++h; ++g; if (g2) ++g2; ++sj;
		VTYPE	cand[4];
		VTYPE&	x = cand[0];
		VTYPE&	y = cand[1];
		EXIN*	bb = bsi.bb;
		recd_t*&	eq1 = hf[1] = e1 + q;
		recd_t*&	eq2 = hf[3] = dagp? e2 + q: 0;
	        *hq = *h;	// reserve previous state
		int	k = 0;

//	Diagonal
		VTYPE	gop = 0;
		STYPE&	sigE = bb->sigE;
		recd_t*	mx = h;
		recd_t*	from = h;
		recd_t*	src = h;
		if (m == a->left) goto Horizon;	// initialize column
		if (n > b->left + 2) {
		    if (sj->dir) {	// phase -1 3'
			*h = *sj;
			sj->dir = 0;
		    } else {
			update(h, h, asi, qprof[*bsi.res] + sigE + gapopen(h, asi, 0), 0);
//			if (sndpro) src->val += sspsim(asi, bsi);
		    }
		    h->dir = isdiag(h)? DIAG: NEWD;
		} else mx->reseth();

//	vertical gap extention
		cand[0] = g[3].val + gapopen(g + 3, asi, 3);
//	1 nt deletion
		++from;
		cand[1] = from->val + (isvert(from)? pwd->GapE1: pwd->GapW1);
//	2 nt deletion
		++from;
		cand[2] = from->val + (isvert(from)? pwd->GapE2: pwd->GapW2);
//	normal deletion
		++from;
		gop = gapopen(from, asi, 3);
		cand[3] = from->val + gop;
		k = vmax(cand, 4) - cand;
		src = k? h + k: g + 3;
		update(g, src, asi, cand[k] - src->val + pua, k? k: 3);
		g->dir = vert[k] + (src->dir & SPIN);
		if (g->val > mx->val) mx = g;

//	long deletion
		if (dagp) {
		    src = g2 + 3;
		    cand[k = 0] = src->val + (VTYPE) (v2divv1 * gapopen(g2 + 3, asi, 3));
		    cand[3] = from->val + (VTYPE) (v2divv1 * gop);
		    if (cand[3] > cand[0]) {src = from; k = 3;}
		    update(g2, src, asi, cand[k] - src->val + (VTYPE) (u2divu1 * pua), k);
		    if (vmf && g2->dir & SPIN)
			g2->ptr = vmf->add(m - 1, n, g2->ptr);
		    if (g2->val > mx->val) mx = g2;
		}
Horizon:
		from = h - 3;
//	long insertion
		if (n > b->left + 2) {
		  gop = gapopen(from, asi, -3);
		  if (dagp) {
		    x = (VTYPE) (v2divv1 * gop);
		    if (from->val + x > eq2->val) {
			src = from;
		    } else {
			src = eq2;
			x = 0;
		    }
		    x += pwd->LongGEP;
		    if (!(src->dir & SPF2)) x += sigE;
		    update(eq2, src, asi, x, -3);
		    eq2->dir = (src->dir & SPIN) + HORL;
		    if (eq2->val > mx->val) mx = eq2;
		  }

//	nomal extension
		  cand[0] = eq1->val;
//	nomal insertion
		  cand[3] = from->val + gop;
		} else {cand[0] = cand[3] = NEVSEL;}
		++from;
//	2 nt insertion
		cand[2] = (n > b->left + 1)? 
		    from->val + (ishori(from)? pwd->GapE2: pwd->GapW2): NEVSEL;
		++from;
//	1 nt insertion
		cand[1] = from->val + (ishori(from)? pwd->GapE1: pwd->GapW1);
		k = vmax(cand, 4) - cand;
		src = k? ++from - k: eq1;
		x = cand[k] - src->val + pwd->BasicGEP;
		if (!(src->dir & SPF2)) x += sigE;
		update(eq1, src, asi, x, k? -k: -3);
		eq1->dir = hori[k] + (src->dir & SPIN);
		if (eq1->val > mx->val) mx = eq1;
		if (++q == 3) q = 0;

//	3' boundary, assume no overlapping signals
		bb += 2;
		if (isEIJ(bb->phs3) && internal) {
		    int phs = (bb->phs3 == 2)? -1: bb->phs3;
Acceptor:
		    int	nb = n - phs;
		    VTYPE sigJ = a_in_zone? api.match_score(3 * m - phs): 0;
		    int*	pnl = nx[phs + 1];
		    vclear(maxphl, Nod + 1);
		    for (int l = 0; l < ncand[phs + 1]; ++l) {
			recd_t*	phl = hl[phs + 1] + pnl[l];
			x = phl->val + sigJ +
                            pwd->IntPen->Penalty(nb - phl->jnc) +
			    b->exin->sig53(phl->jnc, nb, IE53);
			if (phl->dir == 0 && phs) {
			    csi.sqset(spjcs->spjseq(phl->jnc, nb), phs == 1? 0: 1);
			    if (phs == 1) {
				x += pwd->pmt->prematT1(csi.res) + qprof[*csi.res];
			    } else {
				y = x + pwd->pmt->prematT1(csi.res) + qprof1[*csi.res]
				      + gapopen(phl, apsi, 0);
				if (y > mx->val + qprof1[*bsip2.res]) {
				    sj->val = y;
				    maxphl[Nod] = phl;
				}
			    }
			}
			from = hf[phl->dir];
			if (x > from->val) {
			    from->val = x;
			    maxphl[phl->dir] = phl;
			}
		    }
		    if (phs == -1) {
			if (maxphl[0]) sj->dir = 0;
			else if (recd_t* phl = maxphl[Nod]) {
			    sj->dir = NEWD;
			    if (vmf) 
				sj->ptr = vmf->add(m, phl->jnc + phs, phl->ptr);
			}
		    }
		    for (int d = 0; d < Nod; ++d) {
			recd_t*	phl = maxphl[d];
			if (!phl) continue;
			from = hf[d];
		 	if (vmf) {
			    from->ptr = vmf->add(m, n, 
			    vmf->add(m, phl->jnc + phs, phl->ptr));
			}
			from->dir |= SPJCI;
			from->jnc = nb;
			if (from->val > mx->val) mx = from;
		    }
		    if (bb->phs3 - phs == 3) {	// AGAG
			phs = 1;
			goto Acceptor;
		    }
		}

//	Find optimal path
		y = h->val;
		if (mx != h) *h = *mx;
		else if (Local && y > hq->val) {
		    if (LocalL && hq->dir == 0 && !(h->dir & SPJC))
			h->ptr = vmf? vmf->add(m - 1, n - 3, 0): 0;
		    else if (LocalR && y > maxh.val) {
			maxh.val = y;
			maxh.p = h->ptr;
			maxh.m = m;
			maxh.n = n;
		    }
		}
		if (LocalL && h->val <= 0) h->reseth();
		else if (vmf && h->dir == NEWD)
		    h->ptr = vmf->add(m - 1, n - 3, h->ptr);

//	5' boundary
		if (isEIJ(bb->phs5) && internal) { 
		    int	phs = (bb->phs5 == 2)? -1: bb->phs5;
Donor:
		    int	nb = n - phs;
		    VTYPE	sigJ = b->exin->sig53(nb, 0, IE5);
		    int		hd = dir2nod[mx->dir & 15];
		    for (int k = (hd == 0 || phs == 1)? 0: 1; k < Nod; ++k) {
			bool	crossspj = phs == 1 && k == 0;
			from = crossspj? hq: hf[k];
			if (!from->dir || (from->dir & SPIN)) continue;
			if (!crossspj && k != hd && hd >= 0) {
			    y = mx->val;
			    if (hd == 0 || (k - hd) % 2) y += pwd->GOP[k / 2];
			    if (from->val <= y) continue;	// prune
			}
			x = from->val + sigJ;
			recd_t*	phl = hl[phs + 1];
			int*	pnl = nx[phs + 1];
			int&	nc = ncand[phs + 1];
			int	l = nc < NCAND_H? ++nc: NCAND_H;
			while (--l >= 0) {
			    if (x > phl[pnl[l]].val)
				swap(pnl[l], pnl[l + 1]);
			    else
				break;
			}
			if (++l < INTR) {
			    phl += pnl[l];
			    phl->val = x;
			    phl->jnc = nb;
			    phl->dir = k;
			    phl->ptr = (vmf && crossspj && isntdiag(from))?
				vmf->add(m - 1, n - 3, from->ptr): from->ptr;
			} else --nc;
		    }
		    if (bb->phs5 - phs == 3) {	//	GTGT..
			phs = 1;
			goto Donor;
		    }
		}
#if F2DEBUG
	if (OutPrm.debug) {
	    printf("%2d %2d %2d ", m, n, mx->dir);
	    putvar(mx->val); putvar(y);
	    for (int i = 1; i < Nod; ++i) putvar(hf[i]->val); 
	    if (algmode.lsg) {
		putvar(hl[2][nx[2][0]].val);
		putvar(hl[1][nx[1][0]].val);
		putvar(hl[0][nx[0][0]].val);
  		printf(" %2d %2d %5.2f", bb->phs5, bb->phs3, (double) sigE);
		printf(" %5.1f %5.1f", (double) bb->sig5, (double) bb->sig3);
	    }
	    putchar('\n');
	}
#endif
	    }	// end of n-loop
	    if (a_in_zone) ++api;	// has exon-exon junctions
	}	// end of m-loop
	recd_t*	mx = lastH();
	if (maxh.val > mx->val) {
	    if (pp) *pp = vmf->add(maxh.m, maxh.n, maxh.p);
	} else {
	    maxh.val = mx->val;
	    if (pp) *pp = mx->ptr;
	}
	return (maxh.val);
}

template <class recd_t>
VTYPE Fwd2h<recd_t>::verify(Gsinfo* gsi)
{
static	const char*	fswarn = "%s %s FrameShift: %d %d\n";
	SKL*	wsk = gsi->skl;
	int	num = (wsk++)->n;
	recd_t*	h = recbuf;
	recd_t*	hi = 0;
	VTYPE	ha = 0;
	VTYPE	hb = 0;
	VTYPE	sig5 = 0;
	VTYPE	sig3 = 0;
	int	insert = 0;
	int	deletn = 0;
	int	intlen = 0;
	int	preint = 0;
	int	phs = 0;
	int	nb = 0;
	EISCR	rbuf;
	FSTAT*  fst = &gsi->fstat;
	FSTAT	lst;
	FSTAT   pst;
	Eijnc*	eijnc = gsi->eijnc = new Eijnc();
	Cigar*	cigar = gsi->cigar = 0;
	Vulgar*	vlgar = gsi->vlgar = 0;
	switch (algmode.nsa) {
	    case CIG_FORM: cigar = gsi->cigar = new Cigar(); break;
	    case VLG_FORM: vlgar = gsi->vlgar = new Vulgar(); break;
	    default: break;
	}
	vclear(fst);
	vclear(&lst);
	vclear(&pst);
	vclear(&rbuf);
	gsi->noeij = 0;
	mSeqItr	csi(b, wsk[num-1].n - 2);
	bool	termcodon = IsTerm(*csi.res);
	csi.reset();
	CHAR*	cs = 0;
	if (termcodon) wsk[num-1].n -= 3;
	if (wsk[1].n == wsk->n && b->inex.exgl) {++wsk; --num;}
	int	m = wsk->m;
	int	n = wsk->n;
	asi.reset(m);
	bsi.reset(n);
	PfqItr	api(a, m);
	int	api_size = api.size();
	bool	usespb = api_size && use_spb();
	mSeqItr	brim(b, b->right - 1);
	SeqThk	bthk = {b->sumwt, 0, b->sumwt};
	mSeqItr	bzsi;
	bzsi.dns = &bthk;

	h->val = 0;
	if ((algmode.lcl & (16 + 1)) && bsi.bb[bbt].sigS > 0)
	     h->val = bsi.bb[bbt].sigS;
	if ((algmode.lcl & (16 + 4)) && bsi.bb->sig3 > h->val)
	     h->val = bsi.bb->sig3;
	rbuf.left = n;
	rbuf.rleft = m;
	rbuf.phs = phs;
	rbuf.iscr = NEVSEL;
	rbuf.sig3 = h->val;
	VTYPE	hvl = 0;
	if (cigar && m) if (m) cigar->push('H', m);	// local alignment
	while (--num) {
	    ++wsk;
	    bool	tail = num == 1;		// tail gap?
	    int	mi = (wsk->m - m) * 3;
	    if (insert && (mi || tail)) {		// end of insertion
		h->val += pwd->UnpPenalty3(insert);	// not intron
		if (hi && insert > intlen)	// pre-intron gap
		    hi->val += pwd->UnpPenalty3(insert - intlen);
		if (hi && hi->val >= h->val) {		// intron
		    if (cigar) {
			if (preint) cigar->push('D', preint);
			cigar->push('N', intlen);
		    }
		    if (vlgar) {
			if (preint) vlgar->push('G', 0, preint);
			if (phs == -1) vlgar->push('S', 0, 1); else
			if (phs ==  1) vlgar->push('S', 0, 2); 
			vlgar->push('5', 0, 2);
			vlgar->push('I', 0, intlen - 4);
			vlgar->push('3', 0, 2);
			if (phs == -1) vlgar->push('S', 1, 2); else
			if (phs ==  1) vlgar->push('S', 1, 1); 
		    }
		    hb = ha;
		    if (eijnc && rbuf.right - rbuf.left > 1) eijnc->push(&rbuf);
		    rbuf.left = rbuf.right + intlen;
		    rbuf.rleft = m;
		    rbuf.sig3 = sig3;
		    *h = *hi;
		    hi = 0;
		    insert -= (preint + intlen);
		}
		if (insert) {
		    if (cigar) cigar->push('D', insert);
		    if (vlgar && insert) vlgar->push('G', 0, insert);
		    insert = intlen = preint = 0;
		}
	    }
	    int	ni = wsk->n - n;
	    if (ni && deletn) {		// end of delettion
		if (vlgar && deletn > 2) vlgar->push('G', deletn / 3, 0);
		deletn = 0;
	    }
	    int	i = mi - ni;
	    int	d = (i >= 0)? ni: mi;
	    if (d) {			// match
		if (cigar) cigar->push('M', d);
		if (vlgar) vlgar->push('M', d / 3, d);
		n += d;
		m += d / 3;
		for ( ; d > 2; d -= 3, ++asi, bsi += 3) {
		    lst = *fst;
		    if (cs) {		// exon boundary
			csi.sqset(++cs);
			hvl = (pwd->*pwd->sim2)(asi, csi) + 
				pwd->pmt->prematT(csi.res);
			(pwd->*pwd->stt2)(fst, asi, csi);
			cs = 0;
		    } else { 
			mSeqItr	bsip1(b, bsi.pos + 1);
			hvl = (pwd->*pwd->sim2)(asi, bsip1) + bsip1.bb->sigE;
			(pwd->*pwd->stt2)(fst, asi, bsip1);
		    }
		    VTYPE	gop = bsi.dullend()? 0: gapopen(h, asi, 0);
		    fst->gap += gop;
		    h->val += hvl + gop;
		    update(h, asi, 0);
		}
	    }
	    if (i > 0) {		// deletion
		if ((phs = i % 3)) {	// deletion frame shift
		    rbuf.right = n + phs;
		    rbuf.rright = m;
		    rbuf.iscr = NEVSEL;
		    prompt(fswarn, a->sqname(), b->sqname(), b->SiteNo(rbuf.right), -phs);
		    if (eijnc) eijnc->push(&rbuf);
		    if (vlgar) vlgar->push('F', phs, 0);
		    rbuf.left = n;
		    rbuf.rleft = m;
		    h->val += (phs == 1? pwd->GapE1: pwd->GapE2);
		    if (phs == 2) ++asi;
		    phs = 3 - phs;
		    bsi += phs;
		}
		for (int j = 0; j < i / 3; ++j, ++asi) { 
		    VTYPE	gop = bsi.dullend()? 0: gapopen(h, asi, 3);
		    fst->gap += gop;
		    (pwd->*pwd->stt2)(fst, asi, bzsi);
		    h->val += gop + pwd->GapExtPen3(d);
		    update(h, asi, 3);
		}
		deletn += i;
		if (cigar) cigar->push('I', i);
	    } else if (i < 0) {		// insertion
		int	n3 = n + (i = -i);
		EXIN*	b3 = bsi.bb + i;
		VTYPE	xi = NEVSEL;
		int	phs5 = (bsi.bb->phs5 == 2)? b3->phs3: bsi.bb->phs5;
		int	phs3  = (b3->phs3 == 2)? bsi.bb->phs5: b3->phs3;
		if (!hi && i >= IntronPrm.llmt && wsk->n < b->right) {	// intron?
		    VTYPE	xm = NEVSEL;
		 	// potential intron
		    if (phs3 == 2 && phs5 == 2) {	// GTGT....AGAG
			// phs3 = phs5 = -1; 
			nb = n + 1; n3 = nb + i;
			if (isJunct(phs3, phs5)) {
			    xm = b->exin->sig53(nb, n3, IE5P3);
			    if (usespb) xm += api.match_score(3 * m + 1);
			}
			phs3 = phs5 = 1;
		    }
		    nb = n - phs3; n3 = nb + i;
		    if (isJunct(phs3, phs5)) {
			sig5 = b->exin->sig53(nb, 0, IE5);
			sig3 = b->exin->sig53(nb, n3, IE53);
			xi = b->exin->sig53(nb, n3, IE5P3);
			if (usespb) xi += api.match_score(3 * m - phs3);
			preint = insert;
			if (phs3 != 0) cs = spjcs->spjseq(nb, n3);
			if (insert == 0 && phs3 == 1) {
			    mSeqItr asim1(a, asi.pos - 1);
			    csi.sqset(cs);
			    xi += (pwd->*pwd->sim2)(asim1, csi) - hvl
				+ pwd->pmt->prematT(cs);
			    *fst = lst;
			}
		    }
		    if (xm > xi) {
			xi = xm;
			phs3 = -1;
			nb = n - phs3; n3 = nb + i;
			sig5 = b->exin->sig53(nb, 0, IE5);
			sig3 = b->exin->sig53(nb, n3, IE53);
			cs = spjcs->spjseq(nb, n3);
		    }
		    if (xi > NEVSEL) xi += pwd->IntPen->Penalty(i);
		}
		phs = i % 3;
		if (xi > pwd->GapPenalty3(i - phs)) {
		    if (phs3 != -1) cs = 0;	// intron
		    hi = recbuf + 1;
		    *hi = *h;
		    hi->val += xi;
		    intlen = i;
		    rbuf.right = nb;		// 5' information
		    rbuf.rright = m;
		    rbuf.phs = phs3;
		    rbuf.iscr = xi;
		    rbuf.sig5 = sig5;
		    rbuf.escr = (VTYPE) h->val;
		    ha = rbuf.escr + xi - sig3;
		    rbuf.escr += sig5 - hb;
		    rbuf.mch = (int) (fst->mch - pst.mch);
		    rbuf.mmc = (int) (fst->mmc - pst.mmc);
		    rbuf.gap = (int) ((fst->gap - pst.gap) / alprm.v);
		    rbuf.unp = (int) (fst->unp - pst.unp);
                    pst = *fst;
		} else {			// ordinary insetion
		    if (phs) {			// insertion frame shift
			rbuf.right = n - phs;
			rbuf.rright = m;
			rbuf.iscr = NEVSEL;
			prompt(fswarn, a->sqname(), b->sqname(), 
			b->SiteNo(rbuf.right), phs);
			if (eijnc) eijnc->push(&rbuf);
			if (vlgar) vlgar->push('F', 0, phs);
			rbuf.left = n;
			rbuf.rleft = m;
			h->val += (phs == 1? pwd->GapE1: pwd->GapE2);
		    }
		    VTYPE gop = (bsi.dullend()? 0: gapopen(h, asi, -1));
		    fst->gap += gop;
		    fst->unp += a->sumwt * (i - phs) / 3;
		    mSeqItr	bsip1(b, bsi.pos + 1);
		    h->val += gop;
		    update(h, asi, -i);
		    for (int j = 0; j < i - phs; j += 3, bsip1 += 3)
			h->val += bsip1.bb->sigE;
		}
		bsi += i;
		insert += i;
	    }
	    m = wsk->m;
	    n = wsk->n;
	    if (usespb) while (!api.end() && api.lt(m)) ++api;
	}
	if (bsi < brim) {
	    sig5 = 0;
	    if (algmode.lcl & (16 + 2) && bsi.bb[1].sigT > 0)
		sig5 = bsi.bb[1].sigT;
	    if (algmode.lcl & (16 + 8) && bsi.bb->sig5 > 0 
		    && bsi.bb->sig5 > bsi.bb[1].sigT)
		sig5 = bsi.bb->sig5;
	    h->val += sig5;
	}
	if (eijnc) {
	    rbuf.escr = (VTYPE) (h->val + fst->gap - hb);
	    rbuf.iscr = 0;
	    rbuf.sig5 = sig5;
	    rbuf.right = n;
	    rbuf.rright = m;
	    rbuf.mch = (int) (fst->mch - pst.mch);
	    rbuf.mmc = (int) (fst->mmc - pst.mmc);
	    rbuf.gap = (int) ((fst->gap - pst.gap) / alprm.v);
	    rbuf.unp = (int) ((fst->unp - pst.unp));
	    eijnc->push(&rbuf);
	    rbuf.left = endrng.left;
	    rbuf.right = endrng.right;
	    eijnc->push(&rbuf);
	    eijnc->flush();
	    gsi->noeij = eijnc->size() - 1;
	}
	if (cigar) cigar->flush();
	if (vlgar) {
	    vlgar->push('E', 0, 0);	// dummy
	    vlgar->flush();
	    vlgar->postproc();		// correct match length
	}
	fst->mch /= a->sumwt;
	fst->mmc /= a->sumwt;
	fst->val = (FTYPE) h->val;
	fst->unp /= a->sumwt;
	fst->gap /= pwd->BasicGOP;
	if (termcodon) wsk->n += 3;
	return (gsi->scr = h->val);
}

template <class recd_t>
SKL* Fwd2h<recd_t>::globalH(VTYPE* scr)
{
	mfd = new Mfile(sizeof(SKL));
	SKL	wsk;
	mfd->write((UPTR) &wsk);	// dummy call
	vmf = new Vmf();
	long	ptr;
	*scr = forwardH(&ptr);
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

// interface 

template <class recd_t>
VTYPE HomScoreH(mSeq* seqs[], const PwdM* pwd)
{
	WINDOW	wdw;

	stripe31((const Seq**) seqs, &wdw, alprm.sh);
	Fwd2h<recd_t> alnv(seqs, pwd, &wdw);
	return alnv.forwardH(0);
}

template <class recd_t>
SKL* alignH(mSeq* seqs[], const PwdM* pwd, VTYPE* scr)
{
	WINDOW  wdw;

	stripe31((const Seq**) seqs, &wdw, alprm.sh);	// Recacl. window boundaries
	Fwd2h<recd_t> alnv(seqs, pwd, &wdw);
	return (alnv.globalH(scr));
}

template <class recd_t>
VTYPE skl_rngH(mSeq* seqs[], Gsinfo* gsi, const PwdM* pwd)
{
	Fwd2h<recd_t> fwd2h(seqs, pwd);
	return fwd2h.verify(gsi);
}

#endif // _FWD2H_
