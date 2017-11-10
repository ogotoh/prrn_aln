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

extern	VTYPE	SumCodePot(EXIN* bb, int i, CHAR* cs, PwdM* pwd);
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

static	const	DPUH	blackDPUH = {NEVSEL, 0, 0, 0};
static	const	TraceBackDir	hori[4] = {HORI, HOR1, HOR2, HORI};
static	const	TraceBackDir	vert[4] = {VERT, SLA1, SLA2, VERT};

template <class recd_t>
void	clear(recd_t* rcd);

template <class recd_t>
void	reseth(recd_t* rcd, bool all = true);

template <class recd_t>
void	copyh(recd_t* dst, recd_t* src, bool all = true);

template <class recd_t>
class Fwd2h {
	mSeq**	seqs;
	mSeq*	a;
	mSeq*	b;
	mSeqItr	asi;
	mSeqItr	bsi;
	mSeqItr	bzsi;
	WINDOW*	wdw;
	PwdM*	pwd;
	FTYPE	u2divu1;
	FTYPE	v2divv1;
	Vmf*	vmf;
	int	abt;
	int	bbt;
	int	bbt3;
	Mfile*	mfd;
	SpJunc*	spjcs;
	recd_t*	recbuf;
	recd_t*	ff[NOL - 1];	// F
	recd_t*	hh[NOL + 1];	// H, G matrix +1 SPJ
	recd_t*	hq;		// previous state
	recd_t*	hl[NOLL][3];	// [DIA, HORI, HORL][phase][candidates]
	IDELTA*	gfab;
	void	initializeH(int recsize);
	void	initH();
	recd_t*	lastH();
	VTYPE	gapopen(recd_t* rcd, mSeqItr& asi, int d3);
	VTYPE	ins_penalty(recd_t* rcd, mSeqItr& asi, int d3);
	VTYPE	del_penalty(recd_t* rcd, mSeqItr& asi, int d3, FSTAT* fst);
	void	update(recd_t* rcd, mSeqItr& asi, int d3);
	void	update(recd_t* dst, recd_t* src, mSeqItr& asi, int d3);
	mSeqItr	blank_Seqitr;
public:
	Fwd2h(mSeq** _seqs, PwdM* _pwd, WINDOW* _wdw = 0);
	~Fwd2h();
	VTYPE	forwardH(long pp[]);
	SKL*	globalH(VTYPE* scr);
	VTYPE	verify(Gsinfo* gsi);
};

template <class recd_t>
Fwd2h<recd_t>::Fwd2h(mSeq* _seqs[], PwdM* _pwd, WINDOW* _wdw) :
	seqs(_seqs), wdw(_wdw), pwd(_pwd)
{
	a = seqs[0];
	b = seqs[1];
	asi.reset(a->left, a);
	bsi.reset(b->left, b);
	SeqThk	bthk = {b->sumwt, 0, b->sumwt};
	bzsi.dns = &bthk;
	abt = a->byte;
	bbt = b->many;
	bbt3 = 3 * bbt;
	spjcs = new SpJunc(b);
	u2divu1 = pwd->BasicGEP < 0? (FTYPE) pwd->LongGEP / pwd->BasicGEP: 0;
	v2divv1 = pwd->BasicGOP < 0? (FTYPE) pwd->LongGOP / pwd->BasicGOP: 0;
	mfd = 0; vmf = 0; gfab = 0;
	int	recsize = wdw?
	    wdw->width * pwd->Nrow + 3 * pwd->Noll + NOLL * 3 * (INTR+1):
	    3;
	recbuf = new recd_t[recsize];
	initializeH(recsize);
	if (!wdw) return;
	recd_t*	h = recbuf;
	for (int i = 0; i < pwd->Noll - 1; ++i, h += 3) ff[i] = h;
	for (int i = 0; i < pwd->Nrow; ++i, h += wdw->width)
	    hh[i] = h - wdw->lw + 3;
	hq = h; h += pwd->Noll;
	for (int i = 0; i < NOLL; ++i)
	    for (int j = 0; j < 3; ++j, h += INTR+1)
		hl[i][j] = h;
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
		h->val = cand[k];
		if (k) {		// extension
		    update(h, h - k, asi, -k);
		    h->dir = hori[k];
		} else {		// new start
		    reseth(h, false);
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
		update(h, h + i, asi, i);
		h->val += (i == 1? pwd->GapW1: pwd->GapW2);
		h->dir = VERT + i;
	    } else {
		update(h, h + 3, asi, 3);
		h->val += (i == 3)? pwd->GapPenalty(3): pwd->GapExtPen3(i);
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
		    copyh(h, h - 3);
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
	bool	alcl = a->inex.exgl || a->inex.exgr;
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
	VTYPE	pub = pwd->BasicGEP;
	mSeqItr apsi(a, m);
	for (mSeqItr asi = apsi; ++m <= a->right; asi = apsi) {
	    ++apsi;
	    bool	internal = !a->inex.exgr || m < a->right;
	    n1 += 3; n2 += 3;
	    int	n0 = max(n1, b->left);
	    int	n9 = min(n2, b->right);
	    int	n = n0;
	    int	r = n - 3 * m;
	    VTYPE	pua = alcl? 0: (pwd->*pwd->unpa)(asi, bzsi);
	    recd_t*	hf[NOLL];		// [DIAG, HORI, HORL]
	    recd_t*	hg[NOLL + 1];		// [DIAG, VERT, VERL, SPJ]
	    int		nl[NOLL][3][INTR+1];	// [DIA, HORI, HORL][phase][candidates]
	    for (int k = 0; k < pwd->Nrow; ++k) hg[k] = hh[k] + r;
	    for (int k = 0; k < pwd->Noll - 1; ++k) 
		for (int p = 0; p < 3; ++p) reseth(ff[k] + p);
	    for (int k = 0; k < pwd->Noll; ++k) {
	      reseth(hq + k);
	      for (int p = 0; p < 3; ++p) {
		recd_t*	phl = hl[k][p];
		for (int l = 0; l <= INTR; ++l, ++phl)
		  {nl[k][p][l] = l; reseth(phl);}
	      }
	    }
#if F2DEBUG
	    if (algmode.nsa & 8) {
		printf("%2d %2d %2d", m, n, hg[0]->dir);
		putvar(hg[0]->val); putchar('\n');
	    }
#endif
	    bool	a_in_zone = api_size && api.eq(m);
	    bsi.reset(n - 1);			// center of triplet
	    mSeqItr	bsip2(b, n + 1);	// center of triplet
	    for (int q = 0; ++n <= n9; ++bsi, ++bsip2) {
		VTYPE	cand[4];
		VTYPE&	x = cand[0];
		VTYPE&	y = cand[1];
		EXIN*	bb = bsi.bb;
		for (int k = 0; k < pwd->Nrow; ++k) ++hg[k];
		for (int k = 0; k < pwd->Noll - 1; ++k) hf[k+1] = ff[k] + q;

//	Diagonal
		VTYPE	sigE = bb->sigE;
		recd_t*	src = hf[0] = hg[0];
		recd_t*	mx = src;
		recd_t*	from = src;
		recd_t*	sj = hg[pwd->Noll];
		recd_t*	tgt = hg[1];
		int		k = 0;
		if (m == a->left) goto Horizon;	// initialize column
		if (n > b->left + 2) {
		    copyh(hq, src);	// reserve
		    if (sj->dir) {	// phase -1 3'
			copyh(src, sj);
			sj->dir = 0;
		    } else {
			src->val += (pwd->*pwd->sim2)(asi, bsi) + sigE
				 + gapopen(src, asi, 0);
//			if (sndpro) src->val += sspsim(asi, bsi);
		    }
		    update(src, asi, 0);
		    src->dir = isdiag(src)? DIAG: NEWD;
		} else reseth(mx);

//	vertical gap extention
		cand[0] = tgt[3].val;
//	1 nt deletion
		++from;
		cand[1] = from->val + (isvert(from)? pwd->GapE1: pwd->GapW1);
//	2 nt deletion
		++from;
		cand[2] = from->val + (isvert(from)? pwd->GapE2: pwd->GapW2);
//	normal deletion
		++from;
		cand[3] = from->val + gapopen(from, asi, 3);
		k = vmax(cand, 4) - cand;
		tgt->val = cand[k];
		tgt->dir = vert[k];
		src = k? src + k: tgt + 3;
		if (alcl) pua = (pwd->*pwd->unpa)(asi, bzsi);
		tgt->val += pua;
		update(tgt, src, asi, k? k: 3);
		if (vmf && tgt->dir & SPIN)
		    tgt->ptr = vmf->add(m - 1, n, tgt->ptr);
		if (tgt->val > mx->val) mx = tgt;

//	long deletion
		if (pwd->Noll == 3) {
		    tgt = hg[2];
		    src = tgt + 3;
		    cand[3] = from->val + (VTYPE) (v2divv1 * gapopen(from, asi, 3));
		    if (cand[3] > src->val) {
			src = from;
			tgt->val = cand[3];
			tgt->dir = VERL;
		    } else {
			tgt->val = src->val;
		    }
		    tgt->val += (VTYPE) (u2divu1 * pua);
		    update(tgt, src, asi, 3);
		    if (vmf && tgt->dir & SPIN)
			tgt->ptr = vmf->add(m - 1, n, tgt->ptr);
		    if (tgt->val > mx->val) mx = tgt;
		}
Horizon:
		from = hg[0] - 3;
//	long insertion
		if (n > b->left + 2) {
		  if (pwd->Noll == 3) {
		    tgt = src = hf[2];
		    copyh(hq + 2, tgt);
		    x = from->val + (VTYPE) (v2divv1 * gapopen(from, asi, -3));
		    if (x > src->val) {
			src = from;
			tgt->val = x;
		    }
		    tgt->val += pwd->LongGEP;
		    if (!(src->dir & SPF2)) tgt->val += sigE;
		    tgt->dir = (src->dir & SPIN) + HORL;
		    update(tgt, src, asi, -3);
		    if (tgt->val > mx->val) mx = tgt;
		  }

		  src = from;
		  tgt = hf[1];
		  copyh(hq + 1, tgt);
//	nomal extension
		  cand[0] = tgt->val;
//	nomal insertion
		  cand[3] = from->val + gapopen(from, asi, -3);
		} else {cand[0] = cand[3] = NEVSEL;}
		++from;
//	2 nt insertion
		cand[2] = (n > b->left + 1)? 
		    from->val + (ishori(from)? pwd->GapE2: pwd->GapW2): NEVSEL;
		++from;
//	1 nt insertion
		cand[1] = from->val + (ishori(from)? pwd->GapE1: pwd->GapW1);

		k = vmax(cand, 4) - cand;
		src = k? ++from - k: tgt;
		tgt->val = cand[k];
		tgt->dir = hori[k] + (src->dir & SPIN);
		tgt->val += pwd->BasicGEP;
		if (!(src->dir & SPF2)) tgt->val += sigE;
		update(tgt, src, asi, k? -k: -3);
		if (tgt->val > mx->val) mx = tgt;
		if (++q == 3) q = 0;

//	3' boundary, assume no overlapping signals
		bb += 2;
		if (isEIJ(bb->phs3) && internal) {
		    int phs = (bb->phs3 == 2)? 1: bb->phs3;
Acceptor:
		    if (phs == -1)
			sj->val = mx->val + gapopen(mx, asi, 0)
			    + (pwd->*pwd->sim2)(apsi, bsip2);
		    int	nb = n - phs;
		    VTYPE sigJ = a_in_zone? api.match_score(3 * m - phs): 0;
		    recd_t* maxdhl = 0;
		    for (k = 0; k < pwd->Noll; ++k) {
			recd_t*	maxphl = 0;
			int*	pnl = nl[k][phs + 1];
			tgt = hf[k];
			for (int l = 0; l < INTR; ++l) {
			    recd_t*	phl = hl[k][phs + 1] + pnl[l];
			    if (!phl->dir) break;
			    y = phl->val + sigJ +
                                pwd->IntPen->Penalty(nb - phl->jnc) +
				b->exin->sig53(phl->jnc, nb, IE53);
			    if (phs == 1) {
				csi.sqset(spjcs->spjseq(phl->jnc, nb));
				y += pwd->pmt->prematT(csi.res);
				if (isdiag(phl)) y += (pwd->*pwd->sim2)(asi, csi);
			    } else if (phs == -1) {
				csi.sqset(spjcs->spjseq(phl->jnc, nb), 1);
				y += pwd->pmt->prematT(csi.res);
				x = y + (pwd->*pwd->sim2)(apsi, csi)
				      + gapopen(phl, apsi, 0);
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
			    copyh(tgt, maxphl, false);
		 	    if (vmf) tgt->ptr = vmf->add(m, n, 
				vmf->add(m, maxphl->jnc + phs, maxphl->ptr));
			    tgt->jnc = nb;
			    tgt->dir = maxphl->dir | SPJCI;
			    if (phs == -1) tgt->dir |= SPF2;
			    if (isvert(maxphl)) copyh(hg[1], tgt);
			    if (maxphl->dir == VERL) copyh(hg[2], tgt);
			    if (tgt->val > mx->val) mx = tgt;
			}
		    }
		    if (maxdhl) {
			copyh(sj, maxdhl, false);
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
		if (mx != tgt) copyh(tgt, mx);
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
		if (LocalL && tgt->val <= 0) reseth(tgt);
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
			x = src->val + sigJ;
			if (phs == 1 && nb - src->jnc == 2 && bb[-2].sigS > 0)
			    x -= bb[-2].sigS;
			recd_t*	phl = hl[k][phs + 1];
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
			    copyh(phl, src);
			    phl->val = x;
			    phl->jnc = nb;
			}
		    }
		    if (bb->phs5 - phs == 1) {	//	GTGT..
			phs = -1;
			goto Donor;
		    }
		}
#if F2DEBUG
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
	    if (a_in_zone) ++api;	// has exon-exon junctions
	}	// end of m-loop
	recd_t*	tgt = lastH();
	pub = tgt->val;
	*pp = tgt->ptr;
	return (pub);
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
	VTYPE	gop = 0;
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
	    case EXN_FORM: 
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
		h->val += ins_penalty(h, asi, insert);
		if (hi && insert > intlen)
		    hi->val += ins_penalty(hi, asi, insert - intlen);
		if (hi && hi->val > h->val) { 		// intron
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
		    copyh(h, hi);
		    insert -= (preint + intlen);
		}
		if (insert) {
		    if (cigar) cigar->push('D', insert);
		    phs = insert % 3;
		    insert -= phs;
		    if (!((a->inex.exgl && m == a->left) ||
			  (a->inex.exgr && m == a->right)))
			fst->gap += gop;
		    fst->unp += pwd->Vab * insert;
		    if (phs) {		// insertion frame shift
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
		    if (vlgar && insert) vlgar->push('G', 0, insert);
		}
		gop = insert = intlen = preint = 0; hi = 0;
	    }
	    if (deletn) {		// end of delettion
		h->val += del_penalty(h, asi, deletn, fst);
		if ((phs = deletn % 3)) {	// deletion frame shift
		    rbuf.right = n + phs;
		    rbuf.rright = m;
		    rbuf.iscr = NEVSEL;
		    prompt(fswarn, a->sqname(), b->sqname(), b->SiteNo(rbuf.right), -phs);
		    if (eijnc) eijnc->push(&rbuf);
		    if (vlgar) vlgar->push('F', phs, 0);
		    rbuf.left = n;
		    rbuf.rleft = m;
		    h->val += (phs == 1? pwd->GapE1: pwd->GapE2);
		    ++asi;
		    phs = 3 - phs;
		    bsi += phs;
		}
		if (vlgar && deletn > 2) vlgar->push('G', deletn / 3, 0);
		deletn = 0;
	    }
	    int	ni = wsk->n - n;
	    int	i = mi - ni;
	    int	d = (i >= 0)? ni: mi;
	    if (d) {
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
		    fst->gap += (bsi.dullend()? 0: gapopen(h, asi, 0));
		    h->val += hvl;
		    update(h, asi, 0);
		}
	    }
	    if (i > 0) {		// deletion
		deletn += i;
		if (cigar) cigar->push('I', i);
	    } else if (i < 0) {		// insertion
		if (!insert) gop = gapopen(h, asi, -3);
		i = -i;
		if (hi) {			// after intron
		    hi->val += SumCodePot(bsi.bb, i, cs? cs + bbt: 0, pwd);
		    update(hi, asi, i);
		} else if (i >= IntronPrm.llmt) {	// intron?
		    EXIN*	b3 = bsi.bb + i;
		    int	phs5 = (bsi.bb->phs5 == 2)? b3->phs3: bsi.bb->phs5;
		    int	n3 = n + i;
		    int	phs3  = (b3->phs3 == 2)? bsi.bb->phs5: b3->phs3;
		    VTYPE	xi = NEVSEL;
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
			    (pwd->*pwd->stt2)(fst, asim1, csi);
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
		    if (xi > NEVSEL) {
			xi += pwd->IntPen->Penalty(i);
			if (phs3 != -1) cs = 0;
			copyh(hi = recbuf + 1, h);
			hi->val += xi;
			intlen = i;
			rbuf.right = nb;			// 5' information
			rbuf.rright = m;
			rbuf.phs = phs3;
			rbuf.iscr = xi;
			rbuf.sig5 = sig5;
			rbuf.escr = (VTYPE) (h->val + fst->gap);
			if (insert) rbuf.escr += ins_penalty(hi, asi, insert);
			ha = rbuf.escr + xi - sig3;
			rbuf.escr += sig5 - hb;
			rbuf.mch = (int) (fst->mch - pst.mch);
			rbuf.mmc = (int) (fst->mmc - pst.mmc);
			rbuf.gap = (int) ((fst->gap - pst.gap) / pwd->BasicGOP);
			rbuf.unp = (int) ((fst->unp - pst.unp) / (3 * pwd->Vab));
                	pst = *fst;
		    }
		}
		h->val += SumCodePot(bsi.bb, i, 0, pwd);
		update(h, asi, -i);
		bsi += i;
		insert += i;
	    }
	    m = wsk->m;
	    n = wsk->n;
	    if (usespb) while (api.lt(m)) ++api;
	}
	if (bsi < brim) {
	    if (algmode.lcl & (16 + 2) && bsi.bb[1].sigT > 0)
		sig5 = bsi.bb[1].sigT;
	    if (algmode.lcl & (16 + 8) && bsi.bb->sig5 > 0 
		    && bsi.bb[1].sig5 > bsi.bb[1].sigT)
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
	    rbuf.gap = (int) ((fst->gap - pst.gap) / pwd->BasicGOP);
	    rbuf.unp = (int) ((fst->unp - pst.unp) / (3 * pwd->Vab));
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
	fst->unp /= 3;
	fst->val = (FTYPE) h->val + fst->gap;
	fst->gap /= pwd->BasicGOP;
	if (termcodon) wsk->n += 3;
	return (h->val);
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
VTYPE HomScoreH(mSeq* seqs[], PwdM* pwd)
{
	WINDOW	wdw;

	stripe31((Seq**) seqs, &wdw, alprm.sh);
	Fwd2h<recd_t> alnv(seqs, pwd, &wdw);
	return alnv.forwardH(0);
}

template <class recd_t>
SKL* alignH(mSeq* seqs[], PwdM* pwd, VTYPE* scr)
{
	WINDOW  wdw;

	stripe31((Seq**) seqs, &wdw, alprm.sh);	// Recacl. window boundaries
	Fwd2h<recd_t> alnv(seqs, pwd, &wdw);
	return (alnv.globalH(scr));
}

template <class recd_t>
VTYPE skl_rngH(mSeq* seqs[], Gsinfo* gsi, PwdM* pwd)
{
	Fwd2h<recd_t> fwd2h(seqs, pwd);
	return fwd2h.verify(gsi);
}

#endif // _FWD2H_
