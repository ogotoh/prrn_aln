/*****************************************************************************
*
*	fwd2c.h
*
*	Alignment of two protein or nucleotide sequence
*	 or groups of sequences with Algorithm <<C>>
*
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

#ifndef  _FWD2C_H_
#define  _FWD2C_H_

#define FDEBUG	1

#include "aln.h"
#include "vmf.h"
#include "mseq.h"
#include "maln.h"
#include "mgaps.h"
#include "gfreq.h"
#include "dpunit.h"

template <class recd_t>
class Fwd2c {
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
	void	initializeC();
	VTYPE	gapopen(const recd_t* rcd, const mSeqItr& asi,
		    const mSeqItr& bsi, const int d3);
	void	update(recd_t* dst, const recd_t* src, const mSeqItr& asi, 
		    const mSeqItr& bsi, const VTYPE gnp, const int d3);
	void	initA();	// rectangle
	void	initB();	// banded
	void	initC();	// banded SWG
	recd_t*	lastC();
	void	store_ild(Mfile* mfd, recd_t* h, int r, ISLAND* mxd);
	size_t	bufsize;
	int	an;
	int	bn;
	int*	glab;
	IDELTA*	gfab;
public:
	Fwd2c(mSeq** _seqs, PwdM* _pwd, const bool trb, 
	    bool rectangle = false, WINDOW* pwdw = 0);
	~Fwd2c();
	VTYPE	forwardB(long pp[] = 0);	// banded
	VTYPE	forwardA(long pp[] = 0);	// square
	Colonies* forwardC();
	SKL*	traceback(long rr = -1L) {return vmf->traceback(rr);}
};

template <class recd_t>
Fwd2c<recd_t>::Fwd2c(mSeq* _seqs[], PwdM* _pwd, const bool trb, 
		bool rectangle, WINDOW* pwdw)
	: seqs(_seqs), a(seqs[0]), b(seqs[1]), pwd(_pwd)
{
	u2divu1 = pwd->BasicGEP < 0? (FTYPE) pwd->LongGEP / pwd->BasicGEP: 0;
	v2divv1 = pwd->BasicGOP < 0? (FTYPE) pwd->LongGOP / pwd->BasicGOP: 0;
	if (pwdw) wdw = *pwdw;
	else if (!rectangle) stripe((Seq**) seqs, &wdw, pwd->alnprm.sh);
	vmf = trb? new Vmf(): 0;
	int	l = rectangle? b->right - b->left + 2: wdw.width;
	int	n = (rectangle? b->left: wdw.lw) - 1;
	bufsize = l * pwd->Noll + 3;
	buf = new recd_t[bufsize];
	hdiag = buf;
	f1 = buf + 1;
	f2 = buf + 2;
	hh[0] = buf + 3 - n;
	for (int k = 1; k < pwd->Noll; ++k) hh[k] = hh[k-1] + l;
	initializeC();
}

template <class recd_t>
Fwd2c<recd_t>::~Fwd2c() 
{
	delete vmf;
	delete[] buf;
	delete[] glab;
	delete[] gfab;
}

template <class recd_t>
void Fwd2c<recd_t>::initA()
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
	mSeqItr asi(a, a->left - 1);
	mSeqItr bsi(b, b->left);
	for ( ; ++n < b->right; ++bsi) {
	    (++h)->ptr = vmf? hdiag->ptr: ++r;
	    VTYPE	pub = (pwd->*pwd->unpb)(bsi, asi);
	    VTYPE	gnp = gapopen(h - 1, asi, bsi, -1);
	    gnp = (n - b->left < pwd->codonk1)? gnp + pub: 
		(VTYPE) (v2divv1 * gnp + u2divu1 * pub);
	    update(h, h - 1, asi, bsi, gnp, -1);
	}
}

template <class recd_t>
void Fwd2c<recd_t>::initB()
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
	mSeqItr	bsi(b, n);
	if (wdw.up < rr) rr = wdw.up;
	for ( ; ++r <= rr; ++bsi) {
	    ++h; ++n;
	    VTYPE	pub = (pwd->*pwd->unpb)(bsi, asi);
	    VTYPE	gnp = gapopen(h - 1, asi, bsi, -1);
	    gnp = (n - b->left < pwd->codonk1)? gnp + pub: 
		(VTYPE) (v2divv1 * gnp + u2divu1 * pub);
	    update(h, h - 1, asi, bsi, gnp, -1);
	}

	n = b->left;
	r = n - m;
	rr = b->left - a->right;
	h = hh[0] + r;
	bsi.reset(n - 1);
	if (wdw.lw > rr) rr = wdw.lw;
	while (--r >= rr) {
	    --h; ++m; ++asi;
	    VTYPE	pua = (pwd->*pwd->unpa)(asi, bsi);
	    VTYPE	gnp = gapopen(h + 1, asi, bsi, 1);
	    gnp = (m - a->left < pwd->codonk1)? gnp + pua:
		(VTYPE) (v2divv1 * gnp + u2divu1 * pua);
	    update(h, h + 1, asi, bsi, gnp, 1);
	}
}

template <class recd_t>
void Fwd2c<recd_t>::initC()
{
	int	m = a->left;
	int	n = b->left;
	int	r = n - m;
	int	rr = b->right - a->left;

	if (wdw.up < rr) rr = wdw.up;
	for ( ; r <= rr; ++r) {
	    recd_t*	h = hh[0] + r;
	    h->val = 0;
	    h->dir = DEAD;
	    h->upr = h->lwr = r;
	    h->mlb = a->left;
	    h->nlb = n++;
	}

	r = b->left - a->left;
	rr = b->left - a->right;
	if (wdw.lw > rr) rr = wdw.lw;
	while (--r >= rr) {
	    recd_t*	h = hh[0] + r;
	    h->val = 0;
	    h->dir = DEAD;
	    h->upr = h->lwr = r;
	    h->mlb = ++m;
	    h->nlb = b->left;
	}
}

template <class recd_t>
void Fwd2c<recd_t>::store_ild(Mfile* mfd, recd_t* h, int r, ISLAND* mxd)
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
VTYPE Fwd2c<recd_t>::forwardA(long pp[])
{
	initA();		/* Initialize corners */
	int     m = a->left;
	ISLAND  mxd[2] = {{0, 0, 0}, {0, 0, 0}};
	Mfile*  mfd = pp? new Mfile(sizeof(ISLAND)): 0;
	PfqItr	api(a, m);
	int	api_size = api.size();
	PfqItr	bpi(b);
	int	bpi_size = bpi.size();
	mSeqItr	bsi(b, 0);

	for (mSeqItr asi(a, m); m < a->right; ++asi, ++m) {
	    int	 n = b->left - 1;
	    recd_t*     h = hh[0] + n;
	    *hdiag = *h;
	    VTYPE   gnp = gapopen(h, asi, bsi, 1) +
			(pwd->*pwd->unpa)(asi, bsi);
	    update(h, h, asi, bsi, gnp, 1);
	    recd_t*     g = hh[1] + n;
	    recd_t*     g2 = (pwd->Noll == 3)? hh[2] + n: 0;
	    reset(f1); reset(f2); reset(g); if (g2) reset(g2);
	    bsi.reset(++n);
	    bool	a_in_zone = api_size && bpi_size && (api && m);
	    if (a_in_zone) bpi.reset(n);
	    for ( ; n < b->right; ++n, ++bsi) {
		++h; ++g; if (g2) ++g2;
//      Diagonal
		VTYPE   dab = (pwd->*pwd->sim2)(asi, bsi);
		VTYPE   gop = gapopen(hdiag, asi, bsi, 0);
		update(hdiag, hdiag, asi, bsi, dab + gop, 0);
//      Vertical
		VTYPE   gnp = gapopen(g, asi, bsi, 1);
		gop = gapopen(h, asi, bsi, 1);
		if (isntvert(h) &&
		    (h->val + gop > g->val + gnp))
			update(g, h, asi, bsi, gop, 1);
		else    update(g, g, asi, bsi, gnp, 1);
		VTYPE   pua = (pwd->*pwd->unpa)(asi, bsi);
		g->val += pua;
		recd_t*	mx = g;
//      Vertical2
		if (g2) {
		    gnp = (VTYPE) (v2divv1 * gapopen(g2, asi, bsi, 1));
		    gop = (VTYPE) (v2divv1 + gop);
		    if (isntvert(h) &&
			(h->val + gop > g2->val + gnp))
				update(g2, h,  asi, bsi, gop, 1);
		    else	update(g2, g2, asi, bsi, gnp, 1);
		    g2->val += (VTYPE) (u2divu1 * pua);
		    if (g2->val > mx->val) mx = g2;
		}
//      Horizontal
		VTYPE   pub = (pwd->*pwd->unpb)(bsi, asi);
		gnp = gapopen(f1, asi, bsi, -1);
		if (isnthori(h - 1) &&
		    (h[-1].val + (gop = gapopen(h - 1, asi, bsi, -1)) > f1->val + gnp))
			update(f1, h - 1, asi, bsi, gop, -1);
		else    update(f1, f1, asi, bsi, gnp, -1);
		f1->val += pub;
		if (f1->val >= mx->val) mx = f1;
//      Horizontal2
		if (g2) {
		    gnp = (VTYPE) (v2divv1 * gapopen(f2, asi, bsi, -1));
		    gop = (VTYPE) (v2divv1 * gop);
		    if (isnthori(h - 1) &&
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
		if (mx->val > hdiag->val) {
		    copy(hdiag, mx);
		    if (vmf && (hdiag->dir == NEWV || hdiag->dir == NEWH))
			hdiag->ptr = vmf->add(m, n, hdiag->ptr);
		} else if (hdiag->dir == NEWD)
		    hdiag->ptr = vmf? vmf->add(m, n, hdiag->ptr): n - m;
#if FDEBUG
		if (algmode.nsa & 8) {
		    printf("%2d %2d %6.1lf %6.1lf %6.1lf %2d\n",
			m + 1, n + 1, (double) hdiag->val, (double) g->val,
			(double) f1->val, hdiag->dir);
		}
#endif
		gswap(*hdiag, *h);
	    }   /* end of n-loop */
	    if (a_in_zone) ++api;
	    if (a->inex.exgr) store_ild(mfd, h - 1, --n - m, mxd);
	}       /* end of m-loop */
	if (b->inex.exgr) {
	    --m;
	    for (int n = b->left; n <= b->right; ++n)
		store_ild(mfd, hh[0] + n, n - m, mxd);
	} else if (!a->inex.exgr) {
	    hdiag = hh[0] + b->right - 1;
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
VTYPE Fwd2c<recd_t>::forwardB(long pp[])
{
	if (vmf) vmf->add(0, 0, 0);	// Skip 0-th record
	initB();
	int	m = a->left;
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
	mSeqItr	bsi(b);
	PfqItr	api(a, m);
	int	api_size = api.size();
	PfqItr	bpi(b);
	int	bpi_size = bpi.size();

	for (mSeqItr asi(a, m); m < a->right; ++asi, ++m, ++n1, ++n2) {
	    int		n = max(n1, b->left);
	    int		n9 = min(n2, b->right);
	    int		r = n - m;
	    int		k = 0;
	    bsi.reset(n);
	    bool	a_in_zone = api_size && bpi_size && (api && m);
	    if (a_in_zone) bpi.reset(n);
	    VTYPE	pua = (pwd->*pwd->unpa)(asi, bsi);
	    VTYPE	gnp, pub;
	    recd_t*	h = hh[k] + r;
	    recd_t*	g = hh[++k] + r;
	    recd_t*	g2 = (pwd->Noll == 3)? hh[++k] + r: 0;
	    reset(f1);
	    reset(f2);
#if FDEBUG
	if (algmode.nsa & 8) {
	    printf("%2d %2d %2d", m + 1, n + 1, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    for (; n < n9; ++bsi, ++n, ++h, ++g) {
//	Diagonal
		VTYPE	diag = h->val;
		VTYPE	dab = (pwd->*pwd->sim2)(asi, bsi);
		VTYPE	gop = gapopen(h, asi, bsi, 0);
		update(h, h, asi, bsi, dab + gop, 0);
//	Vertical
		recd_t*	mx = g;
		if (m <= a->left) goto hori;
		if (a->inex.nils) pua = (pwd->*pwd->unpa)(asi, bsi);
		gnp = gapopen(g + 1, asi, bsi, 1);
		gop = gapopen(h + 1, asi, bsi, 1);
		if (isntvert(h + 1) &&
		    (h[1].val + gop > g[1].val + gnp))
			update(g, h + 1, asi, bsi, gop, 1);
		else	update(g, g + 1, asi, bsi, gnp, 1);
		g->val += pua;
//	Vertical2
		if (g2) {
		    gnp = (VTYPE) (v2divv1 * gapopen(g2 + 1, asi, bsi, 1));
		    gop = (VTYPE) (v2divv1 * gop);
		    if (m <= a->left || (isntvert(h + 1) && 
			(h[1].val + gop > g2[1].val + gnp)))
				update(g2, h + 1, asi, bsi, gop, 1);
		    else	update(g2, g2 + 1, asi, bsi, gnp, 1);
		    g2->val += (VTYPE) (u2divu1 * pua);
		    if (g2->val > mx->val) mx = g2;
		}
//	Horizontal
hori:		if (n <= b->left) goto opt;
		pub = (pwd->*pwd->unpb)(bsi, asi);
		gnp = gapopen(f1, asi, bsi, -1);
		gop = gapopen(h - 1, asi, bsi, -1);
		if (isnthori(h - 1) &&
		    (h[-1].val + gop > f1->val + gnp))
			update(f1, h - 1, asi, bsi, gop, -1);
		else	update(f1, f1, asi, bsi, gnp, -1);
		f1->val += pub;
		if (f1->val >= mx->val) mx = f1;
//	Horizontal2
		if (g2) {
		    gnp = (VTYPE) (v2divv1 * gapopen(f2, asi, bsi, -1));
		    gop = (VTYPE) (v2divv1 * gop);
		    if (n <= b->left || (isnthori(h - 1) &&
			(h[-1].val + gop > f2->val + gnp)))
				update(f2, h - 1, asi, bsi, gop, -1);
		    else	update(f2, f2, asi, bsi, gnp, -1);
		    f2->val += (VTYPE) (u2divu1 * pub);
		    if (f2->val >= mx->val) mx = f2;
		}

//	Find optimal path
opt:
		if (a_in_zone) {
		    if (bpi && n) {
			h->val += api.match_score(bpi, true);
			mx->val += api.match_score(bpi, false);
			++bpi;
		    }
		}
		if (mx->val > h->val) copy(h, mx);	// non-diagonal
#if FDEBUG
		if (algmode.nsa & 8) {
		    printf("%2d %2d %2d ", m + 1, n + 1, h->dir);
		    putvar(h->val); putvar(diag); 
		    putvar(g->val); putvar(f1->val);
		    if (g2) {
			putvar(g2->val); putvar(f2->val);
		    }
		    putchar('\n');
		}
#endif
		if (vmf) {
		    if (h->dir == NEWD || h->dir == NEWV || h->dir == NEWH)
			h->ptr = vmf->add(m, n, h->ptr);
		} else if (m == a->left)
		    h->ptr = n - m;
		if (g2) ++g2;
	    }	// end of n-loop
	    if (a_in_zone) ++api;
	}	// end of m-loop

	recd_t*	mx = hh[0] + b->right - a->right;
	if (vmf) mx->ptr = vmf->add(a->right, b->right, mx->ptr);
	if (pp) {
	    pp[0] = mx->ptr;
	    pp[1] = b->left - a->left + mx - hh[0];
	}
	return (mx->val);
}

template <class recd_t>
Colonies* Fwd2c<recd_t>::forwardC()
{
	Colonies* cls = new Colonies();
	COLONY*	cc0 = cls->at();

	initC();
	int	m = a->left;
	PfqItr	api(a, m);
	int	api_size = api.size();
	PfqItr	bpi(b);
	int	bpi_size = bpi.size();
	mSeqItr	bsi(b);
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;

	for (mSeqItr asi(a, m); m < a->right; ++asi, ++m, ++n1, ++n2) {
	    int 	n = max(n1, b->left);
	    int 	n9 = min(n2, b->right);
	    int 	r = n - m;
	    int		k = 0;
	    bsi.reset(n);
	    bool	a_in_zone = api_size && bpi_size && (api && m);
	    if (a_in_zone) bpi.reset(n);
	    VTYPE	gnp, pub, pua;
	    recd_t*	h = hh[k] + r;
	    recd_t*	g = hh[++k] + r;
	    recd_t*	g2 = (pwd->Noll == 3)? hh[++k] + r: 0;
	    reset(f1);
	    reset(f2);
	    recd_t*	hlb = hh[0] + r;
	    recd_t*	hrb = hh[0] + n9 - m;

#if DEBUG
	if (algmode.nsa & 8) {
	    printf("%2d %2d %2d", m+1, n+1, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    for ( ; n < n9; ++bsi, ++n, ++r, ++h, ++g) {
//	Diagonal
		VTYPE	diag = h->val;
		VTYPE	dab = (pwd->*pwd->sim2)(asi, bsi);
		VTYPE	gop = gapopen(h, asi, bsi, 0);
		update(h, h, asi, bsi, dab + gop, 0);
		recd_t*	mx = g;
//	Vertical
		if (m <= a->left) goto horiC;
		pua = (pwd->*pwd->unpa)(asi, bsi);
		gnp = gapopen(g + 1, asi, bsi, 1);
		gop = gapopen(h + 1, asi, bsi, 1);
		if (isntvert(h + 1) &&
		    (h[1].val + gop > g[1].val + gnp))
			update(g, h + 1, asi, bsi, gop, 1);
		else	update(g, g + 1, asi, bsi, gnp, 1);
		g->val += pua;
//	Vertical2
		if (g2) {
		    gnp = (VTYPE) (v2divv1 * gapopen(g2 + 1, asi, bsi, 1));
		    gop = (VTYPE) (v2divv1 * gop);
		    if (isntvert(h + 1) && 
			(h[1].val + gop > g2[1].val + gnp))
				update(g2, h + 1, asi, bsi, gop, 1);
		    else	update(g2, g2 + 1, asi, bsi, gnp, 1);
		    g2->val += (VTYPE) (u2divu1 * pua);
		    if (g2->val > mx->val) mx = g2;
		}
//	Horizontal
horiC:		if (n <= b->left) goto optC;
		pub = (pwd->*pwd->unpb)(bsi, asi);
		gnp = gapopen(f1, asi, bsi, -1);
		gop = gapopen(h - 1, asi, bsi, -1);
		if (isnthori(h - 1) &&
		    (h[-1].val + gop > f1->val + gnp))
			update(f1, h - 1, asi, bsi, gop, -1);
		else	update(f1, f1, asi, bsi, gnp, -1);
		f1->val += pub;
		if (f1->val >= mx->val) mx = f1;
//	Horizontal2	
		if (g2) {
		    gnp = (VTYPE) (v2divv1 * gapopen(f2, asi, bsi, -1));
		    gop = (VTYPE) (v2divv1 * gop);
		    if (isnthori(h - 1) &&
			(h[-1].val + gop > f2->val + gnp))
				update(f2, h - 1, asi, bsi, gop, -1);
		    else	update(f2, f2, asi, bsi, gnp, -1);
		    f2->val += (VTYPE) (u2divu1 * pub);
		    if (f2->val >= mx->val) mx = f2;
		}

//	Find optimal path
optC:
		if (a_in_zone) {
		    if (bpi && n) {
			h->val += api.match_score(bpi, true);
			mx->val += api.match_score(bpi, false);
			++bpi;
		    }
		}
		if (mx->val > h->val) {			// non-diagonal
		    copy(h, mx);
		    if (h->lwr > r) h->lwr = r;
		    if (h->upr < r) h->upr = r;
		} else if (h->val > diag) {
		    if (diag == 0) {		// new colony
			h->upr = h->lwr = r;
			h->mlb = m;
			h->nlb = n;
		    }
		    if (h->val > cc0->val) {	// max local score
			cc0->val = h->val;
			cc0->mrb = m + 1;
			cc0->nrb = n + 1;
			cc0->lwr = h->lwr;
			cc0->upr = h->upr;
			cc0->mlb = h->mlb;
			cc0->nlb = h->nlb;
		    }
		}
		if (h->val < 0) {		/* reset to blank */
		    clear(h); clear(f1); clear(f1);
		    if (g2) {clear(f2); clear(g2);}
		    h->clny = 0;
		}
		if (algmode.mlt > 1 && h->val >= pwd->Vthr && !h->clny) {
		    if (cls->size() >= (int) OutPrm.MaxOut) {	/* Gabage collecton */
			for (recd_t* hw = hlb; hw < hrb; ++hw)	/* mark active colony */
			    if (hw->clny) hw->clny->mark = 1;
			if (algmode.mlt == 2) cls->removeoverlap();
			if (cls->size() >= (int) OutPrm.MaxOut)
			    cls->removelowscore();
			for (recd_t* hw = hlb; hw < hrb; ++hw)
			    if (hw->clny) hw->clny = cls->at(hw->clny->clno);
			for (COLONY* cc = cls->at(1); cc <= cls->end(); ++cc) {
			    cc->mark = 0;
			    cc->clno = cls->index(cc);
			}
		    }
		    h->clny = cls->next();
		}
		COLONY*	cc = h->clny;
		if (cc) {
		    if (h->val > cc->val) {
			cc->val = h->val;
			cc->mrb = m + 1;
			cc->nrb = n + 1;
			cc->lwr = h->lwr;
			cc->upr = h->upr;
			cc->mlb = h->mlb;
			cc->nlb = h->nlb;
		    } else if (h->val <= cc->val - pwd->Vthr) {
			clear(h); clear(f1); clear(g);	/* X-drop */
			if (g2) {clear(f2); clear(g2);}
			h->clny = 0;
		    }
		}
#if DEBUG
	if (algmode.nsa & 8) {
		printf("%2d %2d %2d ", m + 1, n + 1, mx->dir);
		putvar(mx->val); putvar(diag); 
		putvar(g->val); putvar(f1->val);
		if (g2) {
		    putvar(g2->val); putvar(f2->val);
		}
		putchar('\n');
	}
#endif
		if (g2) ++g2;
	    }	/* end of n-loop */
	    if (a_in_zone) ++api;
	}	/* end of m-loop */

	if (algmode.mlt == 2) cls->removeoverlap();
	cls->sortcolonies();
	return (cls);
}

// interface 

template <class recd_t>
VTYPE HomScoreC(mSeq* seqs[], PwdM* pwd, long rr[], bool rectangle = false, WINDOW* pwdw = 0)
{
	Fwd2c<recd_t>	pwa(seqs, pwd, false, rectangle, pwdw);
	return rectangle? pwa.forwardA(rr): pwa.forwardB(rr);
}

template <class recd_t>
SKL* alignC(mSeq* seqs[], PwdM* pwd, VTYPE* scr, bool rectangle = false, WINDOW* pwdw = 0)
{
	Fwd2c<recd_t>	pwa(seqs, pwd, true, rectangle, pwdw);

	*scr = rectangle? pwa.forwardA(0): pwa.forwardB(0);
	return pwa.traceback();
}

template <class recd_t>
SKL* swg2ndC(mSeq* seqs[], PwdM* pwd, Gsinfo* gsi, COLONY* clny)
{
	mSeq*&	a = seqs[0];
	mSeq*&	b = seqs[1];
	gswap(a->left, clny->mlb);
	gswap(a->right, clny->mrb);
	gswap(b->left, clny->nlb);
	gswap(b->right, clny->nrb);
	gsi->skl = align2(seqs, pwd, &gsi->scr, gsi);
	gswap(a->left, clny->mlb);
	gswap(a->right, clny->mrb);
	gswap(b->left, clny->nlb);
	gswap(b->right, clny->nrb);
	return (gsi->skl);
}

template <class recd_t>
Colonies* swg1stC(mSeq* seqs[], PwdM* pwd, WINDOW* pwdw = 0)
{
	Fwd2c<recd_t>	pwa(seqs, pwd, false, false, pwdw);
	return pwa.forwardC();
}

#endif	// _FWD2C_H_
