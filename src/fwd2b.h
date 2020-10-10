/*****************************************************************************
*
*	Alignment of two protein or nucleotide sequence
*	 or groups of sequences.
*	Splicing is NOT considered.
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-4-7 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-2012)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#ifndef  _FWD2B_H_
#define  _FWD2B_H_

#define FDEBUG	0

#include "aln.h"
#include "vmf.h"
#include "mseq.h"
#include "maln.h"
#include "mgaps.h"
#include "gfreq.h"
#include "dpunit.h"

template <class recd_t>
class Fwd2b {
	mSeq**	seqs;
	mSeq*	a;
	mSeq*	b;
	WINDOW	wdw;
	PwdM*	pwd;
	FTYPE	u2divu1;
	FTYPE	v2divv1;
	Vmf*	vmf;
	recd_t*	buf;
	recd_t*	hdiag;
	recd_t*	f1;
	recd_t*	f2;
	recd_t*	hh[3];
	void	initializeB();
	VTYPE	gapopen(const recd_t* rcd, const mSeqItr& asi, const mSeqItr& bsi, const int d3);
	void	update(recd_t* dst, const recd_t* src, const mSeqItr& asi, const mSeqItr& bsi, 
		const VTYPE gpn, const int d3);
	void	initA();
	void	initB();
	void	store_ild(Mfile* mfd, recd_t* h, int r, ISLAND* mxd);
	recd_t*	lastB();
	size_t	bufsize;
	int	an;
	int	bn;
	int*	glab;
	IDELTA*	gfab;
public:
	Fwd2b(mSeq** _seqs, PwdM* _pwd, const bool trb, bool rectangle = false, WINDOW* pwdw = 0);
	~Fwd2b();
	VTYPE	forwardB(long pp[] = 0);
	VTYPE	forwardA(long pp[] = 0);
	SKL*	traceback() {return vmf->traceback(-1L);}
};

template <class recd_t>
Fwd2b<recd_t>::Fwd2b(mSeq* _seqs[], PwdM* _pwd, const bool trb, bool rectangle, WINDOW* pwdw) :
	seqs(_seqs), pwd(_pwd)
{
	a = seqs[0];
	b = seqs[1];
	u2divu1 = pwd->BasicGEP < 0? (FTYPE) pwd->LongGEP / pwd->BasicGEP: 0;
	v2divv1 = pwd->BasicGOP < 0? (FTYPE) pwd->LongGOP / pwd->BasicGOP: 0;
	if (pwdw) wdw = *pwdw;
	else if (!rectangle) stripe((Seq**) seqs, &wdw, pwd->alnprm.sh);
	vmf = trb? new Vmf(): 0;
	int	l = rectangle? b->right - b->left + 2: wdw.width;
	int	n = rectangle? b->left - 1: wdw.lw;
	bufsize = l * pwd->Noll + 3;
	buf = new recd_t[bufsize];
	hdiag = buf;
	f1 = buf + 1;
	f2 = buf + 2;
	hh[0] = buf + 3 - n;
	for (int k = 1; k < pwd->Noll; ++k) hh[k] = hh[k-1] + l;
	initializeB();
}

template <class recd_t>
Fwd2b<recd_t>::~Fwd2b() 
{
	delete vmf;
	delete[] buf;
	delete[] glab;
	delete[] gfab;
}

template <class recd_t>
void Fwd2b<recd_t>::initA()
{
	int     r = b->left - a->left;
	int     n = b->left - 1;
	recd_t* h = hh[0] + n;
	clear(h);

	if (vmf) {
	    vmf->add(0, 0, 0L); // Skip 0-th record
	    h->ptr = vmf->add(a->left, b->left, 0L);
	} else
	    h->ptr = r;
	*hdiag = *h;
	if (a->inex.exgl) {     // semi-global
	    while (++n < b->right) {
		clear(++h);
		if (!vmf) h->ptr = ++r;
	    }
	} else {
	    mSeqItr asi(a, a->left);
	    mSeqItr bsi(b, b->left);
	    for ( ; ++n < b->right; ++bsi) {
		(++h)->ptr = vmf? hdiag->ptr: ++r;
		VTYPE   gnp = gapopen(h - 1, asi, bsi, -1)
			 + (pwd->*pwd->unpb)(bsi, asi);
		update(h, h - 1, asi, bsi, (VTYPE) (alprm.tgapf * gnp), -1);
	    }
	}
}

template <class recd_t>
void Fwd2b<recd_t>::initB()
{
	int	m = a->left;
	int	n = b->left;
	int	r = n - m;
	int	rr = b->right - a->left;
	long	origin = vmf? vmf->add(m, n, 0): r;
	recd_t*	h = hh[0] + r;

	h->val = 0;
	h->dir = DIAG;
	h->ptr = origin;
	mSeqItr	asi(a, m - 1);
	mSeqItr	apsi(a, m);
	mSeqItr	bsi(b, n);
	if (wdw.up < rr) rr = wdw.up;
	for ( ; ++r <= rr; ++bsi) {
	    ++h; ++n;
	    if (a->inex.exgl) {		// semi-global
		h->val = 0;
		h->dir = HORI;
		h->ptr = vmf? origin: r;
	    } else {			// global
		VTYPE	pub = (pwd->*pwd->unpb)(bsi, apsi);
		VTYPE	gpn = gapopen(h - 1, asi, bsi, -1);
		gpn = (n < pwd->codonk1)? gpn + pub: 
		    (VTYPE) (v2divv1 * gpn + u2divu1 * pub);
		update(h, h - 1, asi, bsi, gpn, -1);
	    }
	}

	n = b->left;
	r = n - m;
	rr = b->left - a->right;
	h = hh[0 + r];
	bsi.reset(n - 1, b);
	mSeqItr bpsi(b, n);
	if (wdw.lw > rr) rr = wdw.lw;
	while (--r >= rr) {
	    --h; ++m; ++asi;
	    if (b->inex.exgl) {		// semi-global
		h->val = 0;
		h->dir = VERT;
		h->ptr = vmf? origin: r;
	    } else {			// global
		VTYPE	pua = (pwd->*pwd->unpa)(asi, bpsi);
		VTYPE	gpn = gapopen(h + 1, asi, bsi, 1);
		gpn = (m < pwd->codonk1)? gpn + pua:
		   (VTYPE) (v2divv1 * gpn + u2divu1 * pua);
		update(h, h + 1, asi, bsi, gpn, 1);
	    }
	}
}

template <class recd_t>
void Fwd2b<recd_t>::store_ild(Mfile* mfd, recd_t* h, int r, ISLAND* mxd)
{
	if (h->val > mxd->val) {
	    mxd->val = h->val;
	    mxd->lwr = h->ptr;
	    mxd->upr = r;
	}
	if (!mfd) return;
	++mxd;
	if (h->val > mxd->val) {
	    if (h->val > pwd->Vthr) {
		mxd->val = h->val;
		mxd->lwr = h->ptr;
		mxd->upr = r;
	    }
	} else if (mxd->val > pwd->Vthr + max(h->val, (VTYPE) 0)) {
	    mfd->write(mxd);
	    mxd->val = 0;
	}
}

template <class recd_t>
VTYPE Fwd2b<recd_t>::forwardA(long pp[])
{
	initA();		/* Initialize corners */
	int     m = a->left;
	ISLAND  mxd[2] = {{0, 0, 0}, {0, 0, 0}};
	Mfile*  mfd = pp? new Mfile(sizeof(ISLAND)): 0;
	PfqItr  api(a, m);
	PfqItr  bpi(b);
	mSeqItr bsi(b, 0);
	mSeqItr bpsi(b, 0);
	mSeqItr	apsi(a, m + 1);

	for (mSeqItr asi(a, m); m < a->right; ++asi, ++apsi, ++m) {
	    int	 n = b->left - 1;
	    recd_t*     h = hh[0] + n;
	    *hdiag = *h;
	    mSeqItr csi(b, b->left);
	    if (!b->inex.exgl) {
		VTYPE   gnp = gapopen(h, asi, csi, 1) +
			(pwd->*pwd->unpa)(asi, csi);
		update(h, h, asi, csi, (VTYPE) (alprm.tgapf * gnp), 1);
	    }
	    recd_t*     g = hh[1] + n;
	    recd_t*     g2 = (pwd->Noll == 3)? hh[2] + n: 0;
	    reset(f1);
	    reset(f2);
	    bsi.reset(++n);
	    bpsi.reset(n + 1);
	    bool	a_in_zone = api && m;
	    if (a_in_zone) bpi.reset(n);
	    for ( ; n < b->right; ++n, ++bsi, ++bpsi) {
		++h; ++g; if (g2) ++g2;
//      Diagonal
		VTYPE   dab = (pwd->*pwd->sim2)(asi, bsi);
		VTYPE   gop = gapopen(hdiag, asi, bsi, 0);
		update(hdiag, hdiag, asi, bsi, dab + gop, 0);
//      Vertical
		VTYPE   gnp = gapopen(g, asi, bsi, 1);
		if (h->dir != VERT &&
		    (h->val + (gop = gapopen(h, asi, bsi, 1)) > g->val + gnp))
			update(g, h, asi, bsi, gop, 1);
		else    update(g, g, asi, bsi, gnp, 1);
		VTYPE   pua = (pwd->*pwd->unpa)(asi, bpsi);
		g->val += pua;
		recd_t*     mx = g;
//      Vertical2
		if (g2) {
		    gnp = (VTYPE) (v2divv1 * gapopen(g2, asi, bsi, 1));
		    gop = (VTYPE) (v2divv1 + gop);
		    if (h->dir != VERT &&
			(h->val + gop > g2->val + gnp))
				update(g2, h,  asi, bsi, gop, 1);
		    else	update(g2, g2, asi, bsi, gnp, 1);
		    g2->val += (VTYPE) (u2divu1 * pua);
		    if (g2->val > mx->val) mx = g2;
		}
//      Horizontal
		VTYPE   pub = (pwd->*pwd->unpb)(bsi, apsi);
		gnp = gapopen(f1, asi, bsi, -1);
		if (h[-1].dir != HORI &&
		    (h[-1].val + (gop = gapopen(h - 1, asi, bsi, -1)) > f1->val + gnp))
			update(f1, h - 1, asi, bsi, gop, -1);
		else    update(f1, f1, asi, bsi, gnp, -1);
		f1->val += pub;
		if (f1->val >= mx->val) mx = f1;
//      Horizontal2
		if (g2) {
		    gnp = (VTYPE) (v2divv1 * gapopen(f2, asi, bsi, -1));
		    gop = (VTYPE) (v2divv1 * gop);
		    if (h[-1].dir != HORI &&
			(h[-1].val + gop > f2->val + gnp))
				update(f2, h - 1, asi, bsi, gop, -1);
		    else	update(f2, f2, asi, bsi, gnp, -1);
		    f2->val += (VTYPE) (u2divu1 * pub);
		    if (f2->val >= mx->val) mx = f2;
		}

//      Find optimal path
		if (a_in_zone) {
		    if (bpi && n) {
			hdiag->val += api.match_score(bpi, true);
			mx->val += api.match_score(bpi, false);
			++bpi;
		    }
		}
		if (mx->val > hdiag->val) *hdiag = *mx;
		else if (hdiag->dir == NEWD) {
		    hdiag->ptr = vmf? vmf->add(m, n, hdiag->ptr): n - m;
		}
#if FDEBUG
		if (OutPrm.debug) {
		    printf("%2d %2d %6.1lf %6.1lf %6.1lf %2d\n",
			m + 1, n + 1, (double) hdiag->val, (double) g->val,
			(double) f1.val, hdiag->dir);
		}
#endif
		swap(*hdiag, *h);
	    }   /* end of n-loop */
	    if (a_in_zone) ++api;
	    if (a->inex.exgr) store_ild(mfd, h - 1, --n - m, mxd);
	}       /* end of m-loop */
	if (b->inex.exgr) {
	    --m;
	    for (int n = b->left; n <= b->right; ++n)
		store_ild(mfd, hh[0] + n, n - m, mxd);
	} else if (!a->inex.exgr) {
	    mxd->val = hdiag->val;
	    mxd->lwr = hdiag->ptr;
	}
	if (vmf) mxd->lwr = vmf->add(a->right, b->right, mxd->lwr);
	if (pp) {
	    m = mfd->size();
	    ISLAND*     ild = (ISLAND*) mfd->flush();
	    delete mfd;
	    report_ild(seqs, ild, (INT) m, a->len > b->len);
	    delete[] ild;
	    pp[0] = mxd->lwr;
	    pp[1] = mxd->upr;
	}

#if MONITOR
	fprintf(stderr, "A: (%d-%d), Dist = %6.1f, Records = %ld\n",
	    a->inex.vect, b->inex.vect, x, mxd->lwr);
#endif
	return (mxd->val);
}

template <class recd_t>
recd_t* Fwd2b<recd_t>::lastB()
{
	recd_t*	h9 = hh[0] + b->right - a->right;
	recd_t*	mx = h9;

	if (b->inex.exgr) {
	    int	rw = wdw.up;
	    int	rf = b->right - a->left;
	    if (rf < rw) rw = rf;
	    for (recd_t* h = hh[0] + rw; h > h9; --h)
		if (h->val > mx->val) mx = h;
	}
	if (a->inex.exgr) {
	    int	rw = wdw.lw;
	    int	rf = b->left - a->right;
	    if (rf > rw) rw = rf;
	    else	 rf = rw;
	    for (recd_t* h = hh[0] + rw; h <= h9; ++h)
		if (h->val > mx->val) mx = h;
	}
	int	i = mx - h9;
	int	ar = a->right;	// m9
	int	br = b->right;	// n9
	if (i > 0) ar -= i;
	if (i < 0) br += i;
	if (vmf) {
	    if (i) mx->ptr = vmf->add(ar, br, mx->ptr);
	    mx->ptr = vmf->add(a->right, b->right, mx->ptr);
	}
	return (mx);
}

template <class recd_t>
VTYPE Fwd2b<recd_t>::forwardB(long pp[])
{
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0};
	int	Local = algmode.lcl & 16;
	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
	int	LocalR = Local && a->inex.exgr && b->inex.exgr;
	int	alcl = a->inex.exgl || a->inex.exgr;

	if (vmf) vmf->add(0, 0, 0);	// Skip 0-th record
	initB();
	int	m = a->left;
	int	n1 = m + wdw.lw + 1;
	int	n2 = m + wdw.up;
	PfqItr  api(a, m);
	PfqItr  bpi(b);
	mSeqItr bsi(b, 0);
	mSeqItr bpsi(b, 0);
	mSeqItr	apsi(a, m + 1);

	for (mSeqItr asi(a, m); m < a->right; ++asi, ++apsi, ++m, ++n1, ++n2) {
	    int		n = max(n1, b->left);
	    int		n9 = min(n2, b->right);
	    int		r = n - m;
	    int		k = 0;
	    bsi.reset(n);
	    bpsi.reset(n + 1);
	    bool	a_in_zone = api && m;
	    if (a_in_zone) bpi.reset(n);
	    VTYPE	pua = (pwd->*pwd->unpa)(asi, bsi);
	    recd_t*	h = hh[k] + r;
	    recd_t*	g = hh[++k] + r;
	    recd_t*	g2 = (pwd->Noll == 3)? hh[++k] + r: 0;
	    reset(f1);
	    reset(f2);
#if FDEBUG
	if (OutPrm.debug) {
	    printf("%2d %2d %2d", m + 1, n + 1, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    for ( ; n < n9; ++bsi, ++bpsi, ++n, ++h, ++g) {
//	Diagonal
		VTYPE	diag = h->val;
		VTYPE	dab = (pwd->*pwd->sim2)(asi, bsi);
		VTYPE	gop = gapopen(h, asi, bsi, 0);
		update(h, h, asi, bsi, dab + gop, 0);
//	Vertical
		recd_t*	mx = g;
		if (alcl) pua = (pwd->*pwd->unpa)(asi, bpsi);
		VTYPE	gnp = gapopen(g + 1, asi, bsi, 1);
		if (h[1].dir != VERT &&
		    (h[1].val + (gop = gapopen(h + 1, asi, bsi, 1)) > g[1].val + gnp))
			update(g, h + 1, asi, bsi, gop, 1);
		else	update(g, g + 1, asi, bsi, gnp, 1);
		g->val += pua;
//	Vertical2
		if (g2) {
		    gnp = (VTYPE) (v2divv1 * gapopen(g2 + 1, asi, bsi, 1));
		    gop = (VTYPE) (v2divv1 * gop);
		    if (h[1].dir != VERT && 
			(h[1].val + gop > g2[1].val + gnp))
				update(g2, h + 1, asi, bsi, gop, 1);
		    else	update(g2, g2 + 1, asi, bsi, gnp, 1);
		    g2->val += (VTYPE) (u2divu1 * pua);
		    if (g2->val > mx->val) mx = g2;
		}
//	Horizontal
		VTYPE	pub = (pwd->*pwd->unpb)(bsi, apsi);
		gnp = gapopen(f1, asi, bsi, -1);
		if (h[-1].dir != HORI &&
		    (h[-1].val + (gop = gapopen(h - 1, asi, bsi, -1)) > f1->val + gnp))
			update(f1, h - 1, asi, bsi, gop, -1);
		else	update(f1, f1, asi, bsi, gnp, -1);
		f1->val += pub;
		if (f1->val >= mx->val) mx = f1;
//	Horizontal2
		if (g2) {
		    gnp = (VTYPE) (v2divv1 * gapopen(f2, asi, bsi, -1));
		    gop = (VTYPE) (v2divv1 * gop);
		    if (h[-1].dir != HORI &&
			(h[-1].val + gop > f2->val + gnp))
				update(f2, h - 1, asi, bsi, gop, -1);
		    else	update(f2, f2, asi, bsi, gnp, -1);
		    f2->val += (VTYPE) (u2divu1 * pub);
		    if (f2->val >= mx->val) mx = f2;
		}

//	Find optimal path
		if (a_in_zone) {
		    if (bpi && n) {
			hdiag->val += api.match_score(bpi, true);
			mx->val += api.match_score(bpi, false);
			++bpi;
		    }
		}
		if (mx->vla < h->val) copy(h, mx);	// non-diagonal
		else if (Local && h->val > diag) {
		    if (LocalL && diag == 0)
			h->ptr = vmf? vmf->add(m - 1, n - 1, 0): 0;
		    else if (LocalR && h->val > maxh.val) {
			maxh.val = h->val;
			maxh.p = h->ptr;
			maxh.m = m;
			maxh.n = n;
		    }
		}
#if FDEBUG
		if (OutPrm.debug) {
		    printf("%2d %2d %2d ", m + 1, n + 1, mx->dir);
		    putvar(mx->val); putvar(diag); 
		    putvar(g->val); putvar(f1->val);
		    if (g2) {
			putvar(g2->val); putvar(f2->val);
		    }
		    putchar('\n');
		}
#endif
		if (LocalL && h->val <= 0) clear(h);
		else if (vmf && h->dir == NEWD)
		    h->ptr = vmf->add(m, n, h->ptr);
		else if (!vmf && m == a->left)
		    h->ptr = n - m;
		if (g2) ++g2;
	    }	/* end of n-loop */
	    if (a_in_zone) ++api;
	}	/* end of m-loop */

	if (LocalR) {
	    if (vmf) *pp = vmf->add(maxh.m, maxh.n, maxh.p);
	} else {
	    recd_t*	mx = lastB();
	    maxh.val = mx->val;
	    if (pp) {
		pp[0] = mx->ptr;
		pp[1] = b->left - a->left + mx - hh[0];
	    }
	}
	return (maxh.val);
}

// interface 

template <class recd_t>
VTYPE HomScoreB(mSeq* seqs[], PwdM* pwd, long rr[], bool rectangle = false, WINDOW* pwdw = 0)
{
	Fwd2b<recd_t>	pwa(seqs, pwd, false, rectangle, pwdw);
	return rectangle? pwa.forwardA(rr): pwa.forwardB(rr);
}

template <class recd_t>
SKL* alignB(mSeq* seqs[], PwdM* pwd, VTYPE* scr, bool rectangle = false, WINDOW* pwdw = 0)
{
	Fwd2b<recd_t>	pwa(seqs, pwd, true, rectangle, pwdw);

	*scr = rectangle? pwa.forwardA(0): pwa.forwardB(0);
	return pwa.traceback();
}

#endif	// _FWD2B_H_
