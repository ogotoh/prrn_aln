/**********************************************************************
*
*	Subsequence matching of two protein/nucleotide sequences 
*	by means of Smisth-Waterman-Gotoh algorithm
*
*	See SWG.DOC for detail.
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

#include "aln.h"
#include "autocomp.h"
#include "fspscore.h"
#include "css.h"
#include "vmf.h"

#define MONITOR	    0
#define FDEBUG	    1
#define Epsilon	    (0.01)
#define HEMIPOS	    4
#define HEMINEG	    5

#if FVAL
#define EqualVal(a, b) (ABS((a) - (b)) < Epsilon)
#else
#define EqualVal(a, b) ((a) == (b))
#endif

typedef struct colony {VTYPE val; int ml, nl, mx, nx, clw, cup, clno, mark;} ILAND;
typedef struct {VTYPE val; int lw, up, ml, nl; ILAND* cl;} CRECD;
typedef struct {VTYPE val; long ptr;} RECD;

static	void	setparam(int level);
static	VTYPE	lhmseq_1(Seq* seqs[]);
static	VTYPE	lhmseq_2(Seq* seqs[]);
static	VTYPE	lhmseq_3(Seq* seqs[]);
static	VTYPE	lhmseq(Seq* seqs[]);
static	int	comp(ILAND* a, ILAND* b);
static	int	bymx(ILAND* a, ILAND* b);
static	void	detectoverlap(ILAND* cc);
static	ILAND*	removeoverlap(ILAND* cix);
static	int	forward_1(Seq* seqs[]);
static	int	forward_2(Seq* seqs[]);
static	int	forward_3(Seq* seqs[]);
static	void	swapcolony();
static	int	forward_lhm(Seq* seqs[]);
static	VTYPE	backward(Seq* seqs[], ILAND* cl);
static	void	rearrange(SKL* skl);
static	void	printout();
static	int	getlhm(ILAND* cl);
static	void	pickup(int n);
static	void	repmask(FILE* fd, Seq* src, RANGE* prg);
static	void	eachrng(FILE* fd, Seq* sd, RANGE* prg);
static	void	aminusb(RANGE* prg, ILAND* cl);
static	void	addrng(RANGE* rng, int* pn, RANGE* prg);
static	RANGE*	homrng(Seq* seqs[], RANGE* rng, int* pn);
static	int	hemiseq(Seq* seqs[], int calc);
static	void	others();
static	void	stt_nsm(Seq* seqs[], VTYPE val);
static	int	match_2(Seq* seqs[]);
static	int	match2(Seq* seqs[]);
extern	int	main(int argc, char** argv);

static	const	int	defpam = 100;
static	VMF*	vmf;
static	GAPS*	gaps[2];
static	Seq*	seqs[4];	/* seqs[2] and seqs[3] are dummy */
static	VTYPE	vthr;
static	INT 	top = 0;
static	INT	clno = 0;
static	ILAND*	clny = NULL;
static	CRECD	Zcrcd = {0, POS_INT, NEG_INT, 0, 0, NULL};
static	RECD	Zrcd = {0, 0};
static	Gsinfo	Alni;
static	VTYPE 	hml;
static	VTYPE	vab;

static void setparam(int level)
{
	if (level & 1) {
		setalprm();
		setthr(DQUERY);
	}
	if (level & 2) {
	    setalgmode(QUERY, SILENT);
	    if (algmode.mlt > 0 && setlpw(QUERY))
		setprmode(QUERY, QUERY, QUERY);
	    setminus(QUERY);
	    setshffl(QUERY, QUERY);
	} 
}

/*	a->byte == 1 && b->byte == 1	*/

static VTYPE lhmseq_1(Seq* seqs[])
{
	Seq*	a = seqs[0];
	Seq*	b = seqs[1];
	CHAR*	as;
	CHAR*	bs;
	CHAR*	ts;
	int	n9;
	int 	half = samerange(a, b);
	int 	m = b->right - b->left + 1;
	VTYPE	hh, f, x;
	VTYPE	*hm, *gm;
	VTYPE*	h = (VTYPE *) balloc(m, sizeof(VTYPE)) + 1;
	VTYPE*	g = (VTYPE *) balloc(m, sizeof(VTYPE)) + 1;
	VTYPE 	maxlhm = 0;
	VTYPE 	mvv;
	VTYPE	muu;
	VTYPE*	smtx;

	vab = setsim2(seqs);
	mvv = (VTYPE) (vab * -alprm.v);
	muu = (VTYPE) (vab * -alprm.u);
	as = a->seq + a->left;
	for (m = a->left; m < a->right; ++m, ++as) {
	    hh = f = 0;
	    hm = h - 1;
	    gm = g;
	    if (half == 0)	n9 = b->right;
	    else if (half > 0)  n9 = m;
	    else		n9 = b->right - m;
	    bs = b->seq + b->left;
	    ts = b->seq + n9;
	    smtx = simmtx[*as];
	    for ( ; bs < ts; ++bs, ++gm) {
		x = *hm + mvv;
		f = MAX(x, f) + muu;
		x = *++hm + mvv;
		*gm = MAX(x, *gm) + muu;
		x = hh + smtx[*bs];
		hh = MAX(f, *gm);
		if (x > hh) hh = x;
		if (hh < 0) hh = 0;
		if (hh > maxlhm) maxlhm = hh;
		gswap(*hm, hh);
	    }
	}
	free((void*) (h - 1));
	free((void*) (g - 1));
	return (maxlhm);
}

/*	a->byte > 1 && b->byte == 1	*/

static VTYPE lhmseq_2(Seq* seqs[])
{
	Seq*	a = seqs[0];
	Seq*	b = seqs[1];
	CHAR*	as;
	CHAR*	bs;
	CHAR*	ts;
	int 	half = samerange(a, b);
	int 	m = b->right - b->left + 1;
	int	n9;
	VTYPE*	va;
	VTYPE	hh, f, x;
	VTYPE	*hm, *gm;
	VTYPE*	h = (VTYPE *) balloc(m, sizeof(VTYPE)) + 1;
	VTYPE*	g = (VTYPE *) balloc(m, sizeof(VTYPE)) + 1;
	VTYPE 	maxlhm = 0;
	VTYPE 	mvv;
	VTYPE	muu;

	a->inex.prof = 1;
	vab = setsim2(seqs);
	mvv = (VTYPE) (vab * -alprm.v);
	muu = (VTYPE) (vab * -alprm.u);
	as = a->seq + a->left * a->byte;
	for (m = a->left; m < a->right; ++m, as += a->byte) {
	    hh = f = 0;
	    hm = h - 1;
	    gm = g;
	    if (half == 0)	n9 = b->right;
	    else if (half > 0)  n9 = m;
	    else		n9 = b->right - m;
	    bs = b->seq + b->left;
	    ts = b->seq + n9;
	    va = (VTYPE*) as;
	    for ( ; bs < ts; ++bs, ++gm) {
		x = *hm + mvv;
		f = MAX(x, f) + muu;
		x = *++hm + mvv;
		*gm = MAX(x, *gm) + va[gap_code];
		x = hh + va[*bs];
		hh = MAX(f, *gm);
		if (x > hh) hh = x;
		if (hh < 0) hh = 0;
		if (hh > maxlhm) maxlhm = hh;
		gswap(*hm, hh);
	    }
	}
	free((void*) (h - 1));
	free((void*) (g - 1));
	return (maxlhm);
}

/*	a->byte > 1 && b->byte > 1	*/

static VTYPE lhmseq_3(Seq* seqs[])
{
	Seq*	a = seqs[0];
	Seq*	b = seqs[1];
	CHAR*	as;
	CHAR*	bs;
	CHAR*	ts;
	int 	half = samerange(a, b);
	int 	m = b->right - b->left + 1;
	int	n9;
	VTYPE	hh, f, x;
	VTYPE	*hm, *gm;
	VTYPE*	h = (VTYPE *) balloc(m, sizeof(VTYPE)) + 1;
	VTYPE*	g = (VTYPE *) balloc(m, sizeof(VTYPE)) + 1;
	VTYPE 	maxlhm = 0;
	VTYPE 	mvv;
	VTYPE	muu;

	vab = setsim2(seqs);
	mvv = (VTYPE) (vab * -alprm.v);
	as = a->seq + a->left * a->byte;
	for (m = a->left; m < a->right; ++m, as += a->byte) {
		hh = f = 0;
		hm = h - 1;
		gm = g;
		if (half == 0)	    n9 = b->right;
		else if (half > 0)  n9 = m;
		else		    n9 = b->right - m;
		bs = b->seq + b->left * b->byte;
		ts = b->seq + n9 * b->byte;
		muu = (*unpa)(as, b->thk->bdy);
		for ( ; bs < ts; bs += b->byte, ++gm) {
			x = *hm + mvv;
			f = MAX(x, f) + (*unpb)(bs, a->thk->bdy);
			x = *++hm + mvv;
			*gm = MAX(x, *gm) + muu;
			x = hh + (*sim2)(as, bs);
			hh = MAX(f, *gm);
			if (x > hh) hh = x;
			if (hh < 0) hh = 0;
			if (hh > maxlhm) maxlhm = hh;
			gswap(*hm, hh);
		}
	}

	free((void*) (h - 1));
	free((void*) (g - 1));
	return (maxlhm);
}

#define seqswap(s) {Seq* t = s[0]; s[0] = s[1]; s[1] = t;}

static VTYPE lhmseq(Seq* seqs[])
{
	VTYPE	lhm;

	if (seqs[0]->byte == 1 && seqs[1]->byte == 1)
		return (VTYPE) lhmseq_1(seqs);
	if (seqs[1]->byte == 1)
		return (VTYPE) lhmseq_2(seqs);
	if (seqs[0]->byte == 1) {
		seqswap(seqs);
		lhm = lhmseq_2(seqs);
		seqswap(seqs);
		return (lhm);
	}
	return (VTYPE) lhmseq_3(seqs);
}

static int comp(ILAND* a, ILAND* b)
{
	if (b->val == a->val) return(0);
	return (b->val > a->val? 1: -1);
}

static int bymx(ILAND* a, ILAND* b)
{
	return (a->mx - b->mx);
}

static int byclno(ILAND* a, ILAND* b)
{
	return (a->clno - b->clno);
}

static void detectoverlap(ILAND* cc)
{
	ILAND*	cw = cc;
	int     mlb = cc->ml + OutPrm.AllowdOverlap;

	while ((--cw)->mx > mlb) {
	    if (cw->mark == 1) continue;
	    if (cc->mx - cw->ml > OutPrm.AllowdOverlap &&
		cc->nx - cw->nl > OutPrm.AllowdOverlap &&
		cw->nx - cc->nl> OutPrm.AllowdOverlap) {
		if (cc->val < cw->val)	cc->mark = -1;
		else	cw->mark = -1;
	    }
	}
}

static ILAND* removeoverlap(ILAND* cix)
{
	ILAND*	cw;
	ILAND*	cc;
	int	nc = cix - clny;

	clny->ml = clny->mx = 0;
	qsort((UPTR) (clny + 1), (INT) nc, sizeof(ILAND), (CMPF) bymx);
	for (cc = cix; cc > clny + 1; --cc)
	    if (cc->mark == 0) detectoverlap(cc);	/* non-active */
	qsort((UPTR) (clny + 1), (INT) nc, sizeof(ILAND), (CMPF) byclno);
	for (cc = cw = clny + 1; cw <= cix; ++cw) {
	    if (cw->mark >= 0) {		/* retained */
		nc = cc - clny;
		*cc++ = *cw;
	    } else	nc = 0;
	    cw->clno = nc;
	}
	for (cw = cc; cw <= cix; ) (cw++)->val = 0;
	return (--cc);
}

/*	a->byte == 1 && b->byte == 1	*/

static int forward_1(Seq* seqs[])
{
	register CRECD  *hm, *gm, *nd;
	register CHAR*	bs;
	register int n, r;
	Seq*	a = seqs[0];
	Seq*	b = seqs[1];
	CHAR*	as;
	CRECD	hh, f, f2;
	int 	n9;
	int 	half = samerange(a, b);
	ILAND*	cc;
	ILAND*	cix;
	VTYPE	x;
	VTYPE	mvv, lvv;
	VTYPE	muu, luu;
	VTYPE*	smtx;
	int 	m = b->right - b->left + 1;
	CRECD*	h = (CRECD *) balloc(m, sizeof(CRECD)) + 1;
	CRECD*	g = (CRECD *) balloc(m, sizeof(CRECD)) + 1;
	CRECD*	g2 = NULL;
	CRECD*  gm2 = NULL;
	CRECD*	hw;
	CRECD*	ht;

	vab = setsim2(seqs);
	mvv = (VTYPE) (vab * -alprm.v);
	muu = (VTYPE) (vab * -alprm.u);
	luu = (VTYPE) (vab * -alprm.u1);
	lvv = (VTYPE) (vab * -(alprm.v + alprm.k1 * (alprm.u - alprm.u1)));
	for (n = 0, hm = h - 1, gm = g - 1; n < m; n++)
	    (hm++)->val = (gm++)->val = 0;
	if (alprm.ls == 3) {
	    g2 = (CRECD *) balloc(m, sizeof(CRECD)) + 1;
	    for (n = 0, gm2 = g2 - 1; n < m; n++)
		(gm2++)->val = 0;
	} 
	top = OutPrm.MaxOut;
	if (clny) free(clny);
	clny = (ILAND *) balloc(top + 1, sizeof(ILAND));
	cix = clny;
	as = a->seq + a->left;

	for (m = a->left; m < a->right; ++m, ++as) {
	    hm = h - 1;
	    gm = g;
	    if (g2) gm2 = g2 - 1;
	    if (half == 0)	n9 = b->right;
	    else if (half > 0)  n9 = m;
	    else		n9 = b->right - m;
	    n = b->left;
	    bs = b->seq + n;
	    smtx = simmtx[*as];
	    hh = f = f2 = Zcrcd;
	    ht = h + n9;
	    for (r = n - m; n < n9; ++n, ++r, ++bs, ++gm) {
		if (g2) ++gm2;
		if ((x = hm->val + mvv) > f.val && x > 0)
		    { f = *hm; f.val = x;}
		f.val += muu;
		if ((x = (++hm)->val + mvv) > gm->val && x > 0)
		    { *gm = *hm; gm->val = x;}
		gm->val += muu;
		if (g2) {
		    if ((x = hm->val + lvv) > f2.val && x > 0)
			{ f2 = *hm; f2.val = x;}
		    f2.val += luu;
		    if ((x = (hm)->val + lvv) > gm2->val && x > 0)
			{ *gm2 = *hm; gm2->val = x;}
		    gm2->val += luu;
		}
		x = hh.val;
		hh.val += smtx[*bs];
		nd = (gm->val > f.val)? gm: &f;
		if (g2) {
		    if (f2.val > nd->val) nd = &f2;
		    if (gm2->val > nd->val) nd = gm2;
		}
		if (nd->val > hh.val) {
		    hh = *nd;
		    if (hh.lw > r) hh.lw = r;
		    if (hh.up < r) hh.up = r;
		} else if (hh.val > x) {
		    if (x == 0) {
			hh.lw = hh.up = r;
			hh.ml = m;
			hh.nl = n;
		    }
		    if (hh.val > clny->val) {
			clny->val = hh.val;
			clny->mx = m;
			clny->nx = n;
			clny->clw = hh.lw;
			clny->cup = hh.up;
			clny->ml = hh.ml;
			clny->nl = hh.nl;
		    }
		}
		if (hh.val < 0) {
			hh.val = f.val = gm->val = 0;
			if (g2) f2.val = gm2->val = 0;
			hh.cl = NULL;
		}
		if (hh.val >= vthr && !hh.cl) {
		    if (cix - clny >= OutPrm.MaxOut) {
			algmode.mlt = 2;
			for (hw = h; hw < ht; ++hw)
			    if (hw->cl) hw->cl->mark = 1;
			cix = removeoverlap(cix);
			for (hw = h; hw < ht; ++hw)
			    if (hw->cl) hw->cl = clny + hw->cl->clno;
			for (cc = clny + 1; cc <= cix; ++cc) {
			    cc->mark = 0;
			    cc->clno = cc - clny;
			}
		    }
		    if (cix - clny < OutPrm.MaxOut) {
			hh.cl = ++cix;
			cix->val = cix->mark = 0;
			cix->clno = cix - clny;
		    }
		}
		if ((cc = hh.cl)) {
		    if (hh.val > cc->val) {
			cc->val = hh.val;
			cc->mx = m;
			cc->nx = n;
			cc->clw = hh.lw;
			cc->cup = hh.up;
			cc->ml = hh.ml;
			cc->nl = hh.nl;
		    } else if (algmode.mlt > 1 && hh.val <= cc->val - vthr) {
			hh.val = f.val = gm->val = 0;
			hh.cl = NULL;
		    }
		}
#if FDEBUG
if (algmode.nsa & 8) {
printf("%3d %3d %4d: %5.1lf %5.1lf %5.1lf %5.1lf: %5ld %5d %5d\n", m, n, r, 
(double) x+simmtx[*as][*bs], (double) f.val, (double) gm->val, (double) hh.val,
hh.cl? hh.cl - clny: 0, hh.lw, hh.up);
}
#endif
		gswap(*hm, hh)
	    }
	}

	if (algmode.mlt == 2) cix = removeoverlap(cix);
	if ((top = cix - clny)) {
	    *clny = *cix;
	    qsort((UPTR) clny, (INT) top, sizeof(ILAND), (CMPF) comp);
	}
	if (top > OutPrm.NoOut) top = OutPrm.NoOut;
	clny = (ILAND *) realloc(clny, (top + 1) * sizeof(ILAND));
#if MONITOR
	fprintf(stderr, "LHM = %6.1f %d  lhm >= %d %-10s %-10s\n", 
	    (float) clny->val, top, alprm.thr, sqname(a), sqname(b));
#endif
	free((void*) (h - 1));
	free((void*) (g - 1));
	if (g2) free((void*) (g2 - 1));
	return (top);
}

/*	a->byte > 1 && b->byte == 1	*/

static int forward_2(Seq* seqs[])
{
	register CRECD  *hm, *gm, *nd;
	register CHAR*	bs;
	register int n, r;
	Seq*	a = seqs[0];
	Seq*	b = seqs[1];
	CHAR*	as;
	CRECD	hh, f;
	int 	n9;
	int 	half = samerange(a, b);
	ILAND*	cc;
	ILAND*	cix;
	VTYPE	x;
	VTYPE	mvv;
	VTYPE	muu;
	VTYPE*	va;
	int 	m = b->right - b->left + 1;
	CRECD*	h = (CRECD *) balloc(m, sizeof(CRECD)) + 1;
	CRECD*	g = (CRECD *) balloc(m, sizeof(CRECD)) + 1;

	a->inex.prof = 1;
	vab = setsim2(seqs);
	mvv = (VTYPE) (vab * -alprm.v);
	muu = (VTYPE) (vab * -alprm.u);
	for (n = 0, hm = h - 1, gm = g - 1; n < m; n++)
	    (hm++)->val = (gm++)->val = 0;
	top = OutPrm.MaxOut;
	if (clny) free(clny);
	clny = (ILAND *) balloc(top + 1, sizeof(ILAND));
	cix = clny;
	as = a->seq + a->left * a->byte;

	for (m = a->left; m < a->right; ++m, as += a->byte) {
	    hm = h - 1;
	    gm = g;
	    if (half == 0)	n9 = b->right;
	    else if (half > 0)  n9 = m;
	    else		n9 = b->right - m;
	    n = b->left;
	    bs = b->seq + n;
	    hh = f = Zcrcd;
	    va = (VTYPE*) as;
	    for (r = n - m; n < n9; ++n, ++r, ++bs, ++gm) {
		if ((x = hm->val + mvv) > f.val)
		    { f = *hm; f.val = x;}
		f.val += muu;
		if ((x = (++hm)->val + mvv) > gm->val)
		    { *gm = *hm; gm->val = x;}
		gm->val += va[gap_code];
		x = hh.val;
		hh.val += va[*bs];
		nd = (gm->val > f.val)? gm: &f;
		if (nd->val > hh.val) {
		    hh = *nd;
		    if (hh.lw > r) hh.lw = r;
		    if (hh.up < r) hh.up = r;
		} else if (hh.val > x) {
		    if (x == 0)
			hh.lw = hh.up = r;
		    if (hh.val > clny->val) {
			clny->val = hh.val;
			clny->mx = m;
			clny->nx = n;
			clny->clw = hh.lw;
			clny->cup = hh.up;
		    }
		}
		if (hh.val < 0) {
			hh.val = f.val = gm->val = 0;
			hh.cl = NULL;
		}
		if (hh.val >= vthr && !hh.cl && cix - clny < OutPrm.MaxOut)
			hh.cl = ++cix;
		if ((cc = hh.cl)) {
		    if (hh.val > cc->val) {
			cc->val = hh.val;
			cc->mx = m;
			cc->nx = n;
			cc->clw = hh.lw;
			cc->cup = hh.up;
		    } else if (algmode.mlt > 1 && hh.val <= cc->val - vthr) {
			hh.val = f.val = gm->val = 0;
			hh.cl = NULL;
		    }
		}
		gswap(*hm, hh)
	    }
	}

	if ((top = cix - clny)) {
	    *clny = *cix;
	    qsort((UPTR) clny, (INT) top, sizeof(ILAND), (CMPF) comp);
	}
	if (top > OutPrm.NoOut) top = OutPrm.NoOut;
	clny = (ILAND *) realloc(clny, (top + 1) * sizeof(ILAND));
#if MONITOR
	fprintf(stderr, "LHM = %6.1f %d  lhm >= %d %-10s %-10s\n", 
	    (float) clny->val, top, alprm.thr, sqname(a), sqname(b));
#endif
	free((void*) (h - 1));
	free((void*) (g - 1));
	return (top);
}

/*	a->byte > 1 && b->byte > 1	*/

static int forward_3(/* More Multiple Regions */ Seq* seqs[])
{
	register CRECD  *hm, *gm, *nd;
	register CHAR*	bs;
	register int n, r;
	Seq*	a = seqs[0];
	Seq*	b = seqs[1];
	CHAR*	as;
	CRECD	hh, f;
	int 	n9;
	int 	half = samerange(a, b);
	ILAND*	cc;
	ILAND*	cix;
	VTYPE	x;
	VTYPE	mvv;
	VTYPE	muu;
	int 	m = b->right - b->left + 1;
	CRECD*	h = (CRECD *) balloc(m, sizeof(CRECD)) + 1;
	CRECD*	g = (CRECD *) balloc(m, sizeof(CRECD)) + 1;

	vab = setsim2(seqs);
	mvv = (VTYPE) (vab * -alprm.v);
	for (n = 0, hm = h - 1, gm = g - 1; n < m; n++)
	    (hm++)->val = (gm++)->val = 0;
	top = OutPrm.MaxOut;
	if (clny) free(clny);
	clny = (ILAND *) balloc(top + 1, sizeof(ILAND));
	cix = clny;
	as = a->seq + a->left * a->byte;

	for (m = a->left; m < a->right; m++, as += a->byte) {
	    hm = h - 1;
	    gm = g;
	    if (half == 0)	n9 = b->right;
	    else if (half > 0)  n9 = m;
	    else		n9 = b->right - m;
	    n = b->left;
	    bs = b->seq + n * b->byte;
	    hh = f = Zcrcd;
	    muu = (*unpa)(as, b->thk->bdy);
	    for (r = n - m; n < n9; n++, r++, bs += b->byte, gm++) {
		if ((x = hm->val + mvv) > f.val) {
		    f = *hm;
		    f.val = x;
		}
		f.val += (*unpb)(bs, a->thk->bdy);
		if ((x = (++hm)->val + mvv) > gm->val) {
		    *gm = *hm;
		    gm->val = x;
		}
		gm->val += muu;
		x = hh.val;
		hh.val += (*sim2)(as, bs);
		nd = (gm->val > f.val)? gm: &f;
		if (nd->val > hh.val) {
		    hh = *nd;
		    if (hh.lw > r) hh.lw = r;
		    if (hh.up < r) hh.up = r;
		} else if (hh.val > x) {
		    if (x == 0)
			hh.lw = hh.up = r;
		    if (hh.val > clny->val) {
			clny->val = hh.val;
			clny->mx = m;
			clny->nx = n;
			clny->clw = hh.lw;
			clny->cup = hh.up;
		    }
		}
		if (hh.val < 0) {
			hh.val = f.val = gm->val = 0;
			hh.cl = NULL;
		}
		if (hh.val >= vthr && !hh.cl && cix - clny < OutPrm.MaxOut)
			hh.cl = ++cix;
		if ((cc = hh.cl)) {
		    if (hh.val > cc->val) {
			cc->val = hh.val;
			cc->mx = m;
			cc->nx = n;
			cc->clw = hh.lw;
			cc->cup = hh.up;
		    } else if (algmode.mlt > 1 && hh.val <= cc->val - vthr) {
			hh.val = f.val = gm->val = 0;
			hh.cl = NULL;
		    }
		}
		gswap(*hm, hh)
	    }
	}

	if ((top = cix - clny)) {
	    *clny = *cix;
	    qsort((UPTR) clny, (INT) top, sizeof(ILAND), (CMPF) comp);
	}
	if (top > OutPrm.NoOut) top = OutPrm.NoOut;
	clny = (ILAND *) realloc(clny, (top + 1) * sizeof(ILAND));
#if MONITOR
	fprintf(stderr, "LHM = %6.1f %d  lhm >= %d %-10s %-10s\n", 
	    (float) clny->val, top, alprm.thr, sqname(a), sqname(b));
#endif
	free((void*) (h - 1));
	free((void*) (g - 1));
	return (top);
}

static void swapcolony()
{
	INT	i = 0;
	ILAND*	cix = clny;

	while (i++ <= top) {
		gswap(cix->mx, cix->nx);
		gswap(cix->clw, cix->cup);
		cix->clw = -cix->clw;
		cix->cup = -cix->cup;
		++cix;
	}
}

static int forward_lhm(Seq* seqs[])
{
	int	lhm;

	if (seqs[0]->byte == 1 && seqs[1]->byte == 1)
		return forward_1(seqs);
	if (seqs[1]->byte == 1)
		return forward_2(seqs);
	if (seqs[0]->byte == 1) {
		seqswap(seqs);
		lhm = forward_2(seqs);
		swapcolony();
		seqswap(seqs);
		return (lhm);
	}
	return forward_3(seqs);
}

static VTYPE backward(Seq* seqs[], ILAND* cl)
{
	register CHAR* bs;
	Seq*	a = seqs[0];
	Seq*	b = seqs[1];
	RECD	f, f2;
	RECD*	nd;
	VTYPE	mvv;
	VTYPE	muu;
	VTYPE	nuu;
	VTYPE	lvv;
	VTYPE	luu;
	int 	r, n1;
	int 	r_lw = cl->clw - 1;
	int 	n = cl->cup - cl->clw + 3;
	RECD*	h = (RECD *) balloc(n, sizeof(RECD)) - r_lw;
	RECD*	g = (RECD *) balloc(n, sizeof(RECD)) - r_lw;
	RECD*   g2 = NULL;
	CHAR*	dir = (CHAR *) alloc(n) - r_lw;
	VTYPE 	x;
	double	ur = alprm.u1 / alprm.u;
	int 	m  = cl->mx;
	int 	n0 = m + cl->clw;
	int 	n9 = m + cl->cup;
	CHAR*	as;

	vab = setsim2(seqs);
	mvv = (VTYPE) (vab * -alprm.v);
	if (alprm.ls == 3) {
	    lvv = (VTYPE) (vab * -(alprm.v + alprm.k1 * (alprm.u - alprm.u1)));
	    g2  = (RECD *) balloc(n, sizeof(RECD)) - r_lw;
	}
	as = a->seq + m * a->byte;
	adr(0, 0, 0L);
	for ( ; m >= 0; m--, n0--, n9--, as -= a->byte) {
	    n1 = MAX(n0, 0);
	    n  = MIN(n9, cl->nx);
	    bs = b->seq + n * b->byte;
	    r = n - m;
	    h[r+1] = f = f2 = Zrcd;
	    muu = (*unpa)(as, b->thk->bdy);
	    luu = (VTYPE) (ur * muu);
	    for ( ; n >= n1; n--, r--, bs -= b->byte) {
		if ((x = h[r+1].val + mvv) >= f.val) {
		    f.val = x;
		    f.ptr = h[r+1].ptr;
		}
		nuu = (*unpb)(bs, a->thk->bdy);
		f.val += nuu;
		if ((x = h[r-1].val + mvv) >= g[r-1].val) {
		    g[r].val = x;
		    g[r].ptr = h[r-1].ptr;
		} else
		    g[r] = g[r-1];
		g[r].val += muu;
		nd = (g[r].val > f.val)? g+r: &f;
		if (g2) {
		    if ((x = h[r+1].val + lvv) >= f2.val) {
			f2.val = x;
			f2.ptr = h[r+1].ptr;
		    }
		    f2.val += (VTYPE) (ur * nuu);
		    if ((x = h[r-1].val + lvv) >= g2[r-1].val) {
			g2[r].val = x;
			g2[r].ptr = h[r-1].ptr;
		    } else
			g2[r] = g2[r-1];
		    g2[r].val += luu;
		    if (f2.val > nd->val) nd = &f2;
		    if (g2[r].val > nd->val) nd = g2 + r;
		}
		if (nd->val > (h[r].val += (*sim2)(as, bs))) {
		    h[r] = *nd;
		    dir[r] = 0;
		} else {
		    if (!dir[r] && h[r].val > 0) {
			h[r].ptr = adr(m, n, h[r].ptr);
			dir[r] = 1;
		    }
		    if (EqualVal(h[r].val, cl->val)) goto find;
		}
	    }
	}

find:
	adr(m, n, h[n-m].ptr);
	free((void*) (h + r_lw));
	free((void*) (g + r_lw));
	if (g2) free((void*) (g2 + r_lw));
	free((void*) (dir + r_lw));
	return (cl->val);
}

static void rearrange(SKL* skl)
{
	int	num = (skl++)->n;
	SKL*	nxt = skl + 1;
	int	i;

	while (--num > 0) {
	    i = nxt->m - nxt->n - skl->m + skl->n;
	    if (i > 0) skl->m += i;
	    if (i < 0) skl->n -= i;
	    nxt->m++; nxt->n++;
	    skl++; nxt++;
	}
}

static void printout()
{
	FILE*	prn = autoout("");
	if (algmode.mlt == HEMIPOS || algmode.mlt == HEMINEG)
		hemiseq(seqs, NO);
	print2(prn, seqs, gaps, (float) hml, &Alni, clno, top, 0);
	autoclose();
}

static int getlhm(ILAND* cl)
{
	FILE*	prn = autoout("");

	if (cl->val <= 0) return (ERROR);
	vmf = openvmf(sizeof(SAVE));
	alignvmf(vmf);
	hml = backward(seqs, cl);
	Alni.skl = trcback2(vmf, -1L);
	closevmf(vmf);
	rearrange(Alni.skl);
	calcscore_skl(seqs, &Alni);
	if (getlpw()) {
	    togaps2(gaps, Alni.skl);
	    toimage(gaps, 2);
	    printout();
	} else
	    repalninf(prn, seqs, &Alni);
	FreeGsInfo(&Alni);
	autoclose();
	return (OK);
}

static void pickup(int n)
{
	clno = n;
	if (n) getlhm(clny + n - 1);
	else
	    for (clno = 1; clno <= top; clno++) 
		if (getlhm(clny + clno - 1) == ERROR) break;
}

static void repmask(FILE* fd, Seq* src, RANGE* prg)
{
	int	i;
	int	j;
	CHAR*	s;
	Seq*	dst = copyseq(NULL, src, CPY_ALL);

	if (!dst) return;
	s = dst->seq + dst->left * dst->byte;
	prg->right = 0;
	for (i = dst->left; i < dst->right; ++i) {
	    if (i >= prg->right && !neorng(++prg)) break;
	    if (i >= prg->left)
		for (j = 0; j < dst->byte; j++)
		    *s++ = src->code->amb_code;
	    else
		s += dst->byte;
	}
	typeseq(fd, dst);
	eraseq(dst);
}

static void eachrng(FILE* fd, Seq* sd, RANGE* prg)
{
	RANGE	svr;

	saverange(&sd, &svr, 1);
	while (neorng(++prg)) {
	    restrange(&sd, prg, 1);
	    typeseq(fd, sd);
	}
	restrange(&sd, &svr, 1);
}

static void aminusb(RANGE* prg, ILAND* cl)
{
	SKL*	skl;

	vmf = openvmf(sizeof(SAVE));
	alignvmf(vmf);
	hml = backward(seqs, cl);
	skl = trcback2(vmf, -1L);
	rearrange(skl);
	closevmf(vmf);
	prg->left = skl[1].m;
	prg->right = skl[skl->n].m;
	free(skl);
}

static void addrng(RANGE* rng, int* pn, RANGE* prg)
{
	RANGE*	wkr = rng + *pn - 1;
	RANGE*	inp = rng;

	while (wkr >= rng && wkr->left > prg->right)
	    --wkr;
	if (wkr >= rng) {
	    if (wkr->right >= prg->left) {	/* overlap */
		wkr->left = MIN(wkr->left, prg->left);
		wkr->right = MAX(wkr->right, prg->right);
		return;
	    } else
		inp = ++wkr;
	}
	for (wkr = rng + (*pn)++; wkr > inp; --wkr)
	    *wkr = *(wkr - 1);
	*wkr = *prg;	/* insert */
}

static RANGE* homrng(Seq* seqs[], RANGE* rng, int* pn)
{
	RANGE	prg;
	static int ncl = 0;

	if (forward_lhm(seqs) == EOF) return(NULL);
	if (!rng && (ncl = 2 * top))
		rng = (RANGE*) balloc(ncl + 2, sizeof(RANGE));
	else if (top + *pn > ncl)
		rng = (RANGE*) realloc(rng, (top + *pn + 2) * sizeof(RANGE));
	for (clno = 1; clno <= top; clno++) {
		aminusb(&prg, clny + clno - 1);
		addrng(rng + 1, pn, &prg);
	}
	return (rng);
}

static int hemiseq(Seq* seqs[], int calc)
{
	static RANGE*	rng = NULL;
	int	n = 0;
	FILE*	fd = autoout("");

	vab = setsim2(seqs);
	if (calc || !rng) {
	    vthr = (VTYPE) (vab * alprm.thr);
	    if (rng ) free(rng);
	    rng = homrng(seqs, NULL, &n);
	    if (isdrna(seqs[0])) {
		comrev(seqs[1]);
		rng = homrng(seqs, rng, &n);
		comrev(seqs[1]);
	    }
	    if (rng) {
		rng->left = rng->right = ++n + 1;
		rng[n] = endrng;
	    }
	}
	if (rng) {
	    if (algmode.mlt == HEMINEG)
		repmask(fd, seqs[0], rng);
	    else
		eachrng(fd, seqs[0], rng);
	} else if (algmode.mlt == HEMINEG)
	    typeseq(fd, seqs[0]);
	autoclose();
	return (OK);
}
	
static void others()
{
	int 	n;

	do {
		n = top;
		promptin("Colony ( # <= %d [0: all/ -1: menu] ) : ", &n);
		if (n < 0) return;
	} while (n > top);
	pickup(n);
}

static void stt_nsm(Seq* seqs[], VTYPE val)
{
	double	aval = val / (thickness(seqs[0]) * thickness(seqs[1]));

	fpavsd(stdout, val);
	printf("%6.2lf %6.2lf ", val, aval);
	fprint_seq_mem(stdout, seqs, 2);
}

static int match_2(Seq* seqs[])
{
	VTYPE	val;

	vab = setsim2(seqs);
	vthr = (VTYPE) (vab * alprm.thr);
/*
	if (seqs[0]->many > 1)
	    vthr *= consenseq(NULL, seqs[0]) / selfAlnScr(seqs[0]);
	if (seqs[1]->many > 1)
	    vthr *= consenseq(NULL, seqs[1]) / selfAlnScr(seqs[1]);
*/
	if (algmode.mlt == 0) 
		val = lhmseq(seqs);
	else if (forward_lhm(seqs) == EOF)
		 return(ERROR);
	switch (algmode.mlt) {
		case 0: if (val > vthr) stt_nsm(seqs, val); break;
		case 1: pickup(1); break;
		case 2: pickup(0); break;
		case 3:
		default: pickup(top? 0: 1); break;
	}
	return(OK);
}

static int match2(Seq* seqs[])
{
	int	rv;
	int	wc;

	if (algmode.mlt == HEMIPOS || algmode.mlt == HEMINEG)
		return hemiseq(seqs, YES);
	shuffle(avsd, lhmseq, seqs);
	rv = match_2(seqs);
	if (rv == ERROR || !algmode.mns) return (rv);
	if (!isprotein(seqs[1]))	wc = 1;
	else if (!isprotein(seqs[1]))	wc = 0;
	else	return (rv);
	comrev(seqs[wc]);
	rv = match_2(seqs);
	comrev(seqs[wc]);
	return(rv);
}

int main(int argc, char** argv)
{
	int	 n;

	setpam(defpam, 0., 0);
	initseq(seqs, 2);
	optimize(LOCAL, MAXIMUM);
	setalgmode(2, 2);
/*	Mark between Every recode 	*/
/*	setprmode(Row_Every, 'N', SILENT); */
	n = getoption(argc, argv);
	getargseq(argc - n, argv + n, seqs, getseq, 2);
	if (autocomp(seqs, openseq, getseq, setparam, match2)) {
		setalgmode(1, SILENT);
		getoption(argc, argv);
		menucomp(seqs,inputseq, setparam, match2, others, printout);
	}
	if (gaps[0]) free(gaps[0]);
	if (gaps[1]) free(gaps[1]);
	if (clny) free(clny);
	clearseq(seqs, 2);
	freemtx();
	return (0);
}
