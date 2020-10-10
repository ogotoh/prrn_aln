/*****************************************************************************
*
*	Management of 'static' and 'dynamic' gap states
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

#define DEBUG	0

#include "aln.h"
#include "gfreq.h"
#include "mseq.h"
#if FVAL
#include <math.h>
#endif


int	warnrecovf = 0;

static	int	ipack(GFREQ* gf, GFREQ a, int endg);
static	void	accume(GFREQ* gf, GFREQ* kf);
static	int	cmpf(const GFREQ* a, const GFREQ* b);
static	void	mergegfq(GFREQ* syn, const GFREQ* add, const 
			IDELTA* dlt, int acc, GFREQ* wbuf, FTYPE wf);

#if DEBUG

static	void	printgfq(FILE* fd, GFREQ* w);
static	const	int	afq_dim = 20;
static	const	double	afq_ulmt = 1e-3;
static	const	double	afq_rate = sqrt(0.1l);
static	double	min_afq = FLT_MAX;
static	int	afq_distr[afq_dim];

void print_min_afq()
{
	fprintf(stderr, "Min_afq = %15.7e\n", min_afq);
	double	f = afq_ulmt;
	for (int i = 0; i < afq_dim; ++i, f *= afq_rate) {
	    fprintf(stderr, "%12.4le\t%d\n", f, afq_distr[i]);
	}
}

static void printgfq(FILE* fd, GFREQ* w)
{
	for ( ; neogfq(w); ++w)
#if FVAL
	    fprintf(fd, " %d-%.2f", w->glen, w->freq);
#else
	    fprintf(fd, " %d-%ld", w->glen, long(w->freq));
#endif
}

void printseqgfq(FILE* fd, mSeq* sd, int fb)
{
	for (int i = -1; i <= sd->len; i++) {
	    fprintf(fd, "%4d", i);
	    if (fb > 0) printgfq(fd, sd->gfq->sfrq[i]);
	    putc(':', fd);
	    if (fb > 0) printgfq(fd, sd->gfq->tfrq[i]);
	    putc('\n', fd);
	}
}

#endif // DEBUG

/*
	Insert 'a' into the bottom of 'gf' and eliminate
	those members with 'freq' = 0
*/

static int ipack(GFREQ* gf, GFREQ a, int endg)
{
	GFREQ	tmp;
	GFREQ*  kf = gf;
	GFREQ*  nf = gf;

	do {
	    tmp = *kf++;
	    if (endg || a.nres) *nf++ = a;
	    a = tmp;
	} while (neogfq(&tmp));
	*nf = delmgfq;
	return (nf - gf);
}

/*
	convert to the accumelated form.
	kf points to the end of gf array, which
	is searched if kf = NULL
*/

static void accume(GFREQ* gf, GFREQ* kf)
{
	VTYPE   s = 0;

	if (!kf) for (kf = gf; neogfq(kf++); ) ;
	while (--kf >= gf)
	    kf->freq = s += kf->freq;
}

/*
	copy a set of gfq.
*/

void copygfq(GFREQ* nf, const GFREQ* df)
{
	do {*nf++ = *df;} while (neogfq(df++));
}

static int cmpf(const GFREQ* a, const GFREQ* b)
{
	return (a->glen - b->glen);
}

/*
	conversion of a single column position
*/

int Gfq::seq2gfq(int kk, CHAR* ps)
{
	GFREQ   a = zerogfq;	/*  # of cons res   */
	GFREQ   b = delmgfq;	/*  # of new gaps   */
	GFREQ   key;
	GFREQ*  cf = sbuf;
	GFREQ*  df = tbuf;
	GFREQ*	rf = rbuf;
	VTYPE   w = 1;
	int	endg = 0;
	int*	ln = lbuf;
#if USE_WEIGHT
	FTYPE*	wt = msd->weight;
#else
	FTYPE*	wt = 0;
#endif

	while (neogfq(df)) {
	    df->glen += grain;
	    cf->glen = (df++)->glen;	// glen
	    cf->nres = 0;
	    (cf++)->freq = 0;
	}
	*cf = *df;
	*rf = zerogfq;
	for (int i = 0; i < msd->many; ++i, ++ln, ++ps) {
	    if (wt) w = *wt++;
	    VTYPE	bu = msd->gapdensity(ps - msd->many, i);
	    VTYPE	cu = msd->gapdensity(ps, i);
	    VTYPE	ru = msd->postgapdensity(ps, i);
	    if (cu > 0) {
		if (*ln == 0) {			// *-
		    b.glen = 0;
		    b.freq += w * cu;
		    ++b.nres;
		}
		*ln += grain;			// ?-
	    } else if (*ps == nil_code) {
		if (*ln == 0) ++endg;		// *.
		*ln += grain;			// ..
	    } else if (ps[-msd->many] == nil_code && bu == 0) {		// .*
		if (!neogfq(cf)) {
		    cf[1] = *cf;
		    cf->freq = cf->nres = 0;
		    cf->glen = *ln;		// glen
		}
		cf->freq += w; ++cf->nres;
		if (ru > 0) {rf->freq += w * ru; ++rf->nres;}
		*ln = 0;
	    } else {
		if (ru > 0) {rf->freq += w * ru; ++rf->nres;}		// ?*
		if (*ln) {			// -*
		    key.glen = *ln;		// glen
		    df = (GFREQ*) bsearch(&key, tbuf, kk, sizeof(GFREQ), (CMPF) cmpf);
		    if (df) {
			df->freq -= bu * w; --df->nres;
			cf = df - tbuf + sbuf;
			cf->freq += w; ++cf->nres;
#if DEBUG
#if FVAL
			VTYPE	afq = fabs(df->freq);
			if (afq && afq < afq_ulmt) {
			    int	i = int(2. * log(afq_ulmt / afq) / log(10));
			    if (i >= afq_dim) i = afq_dim - 1;
			    afq_distr[i]++;
			}
			if (!df->nres) df->freq = 0.; else 
			if (afq < min_afq) min_afq = afq;
			afq = fabs(cf->freq);
			if (afq && afq < afq_ulmt) {
			    int	i = int(2. * log(afq_ulmt / afq) / log(10));
			    if (i >= afq_dim) i = afq_dim - 1;
			    afq_distr[i]++;
			}
			if (!cf->nres) cf->freq = 0.; else 
			if (afq < min_afq) min_afq = afq;
#endif
#endif
		    }	/* else fatal("GFR error! %d %f", key.glen, key.freq); */
		} else  {a.freq += w; ++a.nres;}	// **
		*ln = 0;
	    }
	}
	kk = ipack(sbuf, a, 0);
	accume(sbuf, sbuf + kk);
	kk = ipack(tbuf, b, endg);

	if (rf->nres == 0) --rf;
	for (df = tbuf; neogfq(df); ++df) {
	    *++rf = *df;
	    rf->glen += grain;	// glen
	}
	*++rf = *df;
	return (kk);
}

void Gfq::store(const int& kk, const int* wk)
{
	GFREQ**	sf = new GFREQ*[3 * kk];
	GFREQ**	tf = sf + kk;
	GFREQ**	rf = tf + kk;
	sfrq = sf + 1;
	tfrq = tf + 1;
	rfrq = rf + 1;
	for (int k = 0; k++ < kk; ) {
	    *sf++ = gbuf + *wk++;
	    *tf++ = gbuf + *wk++;
	    *rf++ = gbuf + *wk++;
	}
}

//  Make the gap-profile for a character-based multiple-sequence

Gfq::Gfq(mSeq* sd, int gr) : msd(sd), grain(gr)
{
	FTYPE	ltgapf = sd->inex.exgl? 0: alprm.tgapf;
	int	kk = sd->inex.dels? sd->many + 2: 2;
	int	htr = 0;
	CHAR*	ts = sd->at(sd->len);
	CHAR*	ps = sd->at(0);
	CHAR*	ls = sd->inex.dels? 0: sd->at(1);
	CHAR*	rs = sd->inex.dels? 0: sd->at(sd->len - 1);
	sbuf = new GFREQ[3 * kk];
	tbuf = sbuf + kk;
	rbuf = tbuf + kk;
	int	thk_len = sd->inex.dels? sd->len: 2;
	int*	pbuf = new int[3 * thk_len + 7];
	int*	wk = pbuf;
	Mfile	mfd(sizeof(GFREQ));

	sbuf[0] = tbuf[0] = rbuf[0] = sbuf[1] = tbuf[1] = rbuf[1] = delmgfq;
	*wk = 0;
	mfd.write(&delmgfq);	/*  ---*  s  */
	if (ltgapf > 0) {	/*  ----  t  */
	    tbuf->glen = 0;
	    tbuf->freq = sd->sumwt * ltgapf;
	    tbuf->nres = sd->many;
	    *sbuf = *rbuf = *tbuf;
	    *++wk = mfd.size();
	    mfd.write(tbuf);
	    mfd.write(&delmgfq);
	    *++wk = mfd.size();
	    mfd.write(tbuf);	/*  ---?  r  */
	    mfd.write(&delmgfq);
	    tbuf[0] = delmgfq;
	} else {
	    *++wk = 0;		// -*
	    *++wk = 0;		// ??
	}
	lbuf = new int[sd->many];
	vclear(lbuf, sd->many);
	for (kk = 1; ps < ts; ps += sd->many) {
	    if (ls && ps == ls) {ps = rs; ls = 0;}
	    kk = seq2gfq(kk, ps);
	    if (kk > htr) htr = kk;
	    *++wk = mfd.size();
	    GFREQ*	wcd = sbuf;
	    if (neogfq(wcd)) {
		do {mfd.write(wcd);} while (neogfq(wcd++));
	    } else --(*wk);
	    *++wk = mfd.size();
	    wcd = tbuf;
	    if (neogfq(wcd)) {
		do {mfd.write(wcd);} while (neogfq(wcd++));
	    } else --(*wk);
	    *++wk = mfd.size();
	    wcd = rbuf;
	    if (neogfq(wcd)) {
		do {mfd.write(wcd);} while (neogfq(wcd++));
	    } else --(*wk);
	}
	for (int i = 0; i < 3; ++i) *++wk = mfd.size();
	delete[] lbuf;
	delete[] sbuf; 
	gbuf = (GFREQ*) mfd.flush();
	hetero = htr + 1;
	store(thk_len + 2, pbuf);
	delete[] pbuf;
}

/********************************************************
    Makeup gfq from a list of seqs/gaps
********************************************************/

Gfq::Gfq(mSeq** slist, FTYPE* wtlst, GAPS** glist)
	: msd(*slist), hetero(0), grain(0)
{
	int	kk = 0;
	int	k = 2;
	for (mSeq** sq = slist; *sq; ++sq, ++kk) {
	    if ((*sq)->gfq) {
		k += (*sq)->gfq->hetero;
		if (grain != (*sq)->gfq->grain) {
		    if (!grain) grain = (*sq)->gfq->grain;
		    else fatal("Gfq grain %d != %d !\n",
			grain, (*sq)->gfq->grain);
		}
	    } else ++k;
	}
	if (!grain) grain = 1;
	sbuf = new GFREQ[4 * k];
	tbuf = sbuf + k;
	rbuf = tbuf + k;
	GFREQ*	wbuf = rbuf + k;
	GFREQ***	cptr = new GFREQ**[3 * kk];
	GFREQ***	dptr = cptr + kk;
	GFREQ***	rptr = dptr + kk;
	IDELTA**	dlt = new IDELTA*[kk + 1];
	dlt[0] = new IDELTA[k + kk];
	GAPS**	gp = glist? new GAPS*[kk]: 0;
	int	len = glist? gaps_span(*glist): slist[0]->len;
	int	upp = 0;
	int*	pbuf = new int[3 * len + 6];
	int*    wk = pbuf - 1;
	Mfile	mfd(sizeof(GFREQ));

	mSeq**	sq = slist;
	for (kk = 0; *sq; ++sq, ++kk) {
	    cptr[kk] = (*sq)->gfq? (*sq)->gfq->sfrq - 1: 0;
	    dptr[kk] = (*sq)->gfq? (*sq)->gfq->tfrq - 1: 0;
	    rptr[kk] = (*sq)->gfq? (*sq)->gfq->rfrq - 1: 0;
	    if (glist) gp[kk] = glist[kk] + 1;
	    cleardelta(dlt[kk]);
	    dlt[kk + 1] = dlt[kk] + 1 + ((*sq)->gfq? (*sq)->gfq->hetero: 1);
	}
	for (int i = -1; i++ <= len; ) {
	    sbuf[0] = tbuf[0] = rbuf[0] = delmgfq;
	    FTYPE*	wt = wtlst;
	    mSeq**	sq = slist;
	    for (k = 0; k < kk; ++sq, ++k) {
		FTYPE	wf = wt? *wt++: 1;
		if (glist && gaps_intr(gp[k]) && (upp = gp[k]->gps + gp[k]->gln) < i)
		    ++gp[k];
		if (glist && gp[k]->gps < i) {
		    if (((*sq)->inex.exgr && upp == len)) continue;
		    mergegfq(tbuf, rptr[k]? rptr[k][-1]: unitgfq, dlt[k], 0, wbuf, wf);
		    mergegfq(rbuf, rptr[k]? rptr[k][-1]: unitgfq, dlt[k], 2, wbuf, wf);
		    incdelta(dlt[k], dlt[k]);
		} else {
		    mergegfq(sbuf, cptr[k]? *cptr[k]++: unitgfq, dlt[k], 1, wbuf, wf);
		    mergegfq(tbuf, dptr[k]? *dptr[k]: &zerogfq, dlt[k], 0, wbuf, wf);
		    newdelta(dlt[k], dptr[k]? *dptr[k]++:  &delmgfq, dlt[k]);
		    mergegfq(rbuf, rptr[k]? *rptr[k]++: unitgfq, dlt[k], 0, wbuf, wf);
		}
	    }
	    accume(sbuf, 0);		/* to acc form */
	    *++wk = mfd.size();
	    GFREQ*	wcd = sbuf;
	    if (neogfq(wcd)) {
		do {mfd.write(wcd);} while (neogfq(wcd++));
	    } else --(*wk);
	    *++wk = mfd.size();
	    wcd = tbuf;
	    if (neogfq(wcd)) {
		do {mfd.write(wcd);} while (neogfq(wcd++));
	    } else --(*wk);
	    *++wk = mfd.size();
	    wcd = rbuf;
	    if (neogfq(wcd)) {
		do {mfd.write(wcd);} while (neogfq(wcd++));
	    } else --(*wk);
	    k = wcd - rbuf - 1;
	    if (k > hetero) hetero = k;
	}
	delete[] sbuf; delete[] cptr;
	delete[] gp; delete[] dlt[0]; delete[] dlt;
	gbuf = (GFREQ*) mfd.flush();
	store(len + 2, pbuf);
	delete[] pbuf;
}

void mergegfq(GFREQ* syn, const GFREQ* add, const IDELTA* dlt, 
	int acc, GFREQ* wbuf, FTYPE wf)
{
	GFREQ*	wk = wbuf;
	GFREQ*	sn = syn;
const	GFREQ*	ad = add;
	int	alen;

	for ( ; neogfq(ad); ++ad) {
	    alen = ad->glen;
	    if (acc == 2) ++alen;		// glen
	    while (ad->glen >= dlt[1].glen) ++dlt;
	    alen += dlt->nins;
	    while (neogfq(sn) && sn->glen < alen) *wk++ = *sn++;
	    if (neogfq(sn) && sn->glen == alen) {
		wk->freq = sn->freq;
		wk->nres = (sn++)->nres;
	    } else {wk->freq = wk->nres = 0;}
	    if (acc == 1 && neogfq(ad + 1)) {
		wk->freq += wf * (ad->freq - ad[1].freq);
		wk->nres += ad->nres - ad[1].nres;
	    } else {
		wk->freq = wf * ad->freq;
		wk->nres = ad->nres;
	    }
	    wk->glen = alen;
	    if (wk->nres) ++wk;
	    else	wk->freq = wk->nres = 0;
	}
	do {*wk++ = *sn;} while (neogfq(sn++));
	copygfq(syn, wbuf);
}

/********************************************************
    release memory used by a gfq
********************************************************/

FTYPE	meanhetero(mSeq* sd)
{
	FTYPE	av = 0;

	if (!sd->gfq) sd->gfq = new Gfq(sd);
	GFREQ**	df = sd->gfq->tfrq + sd->left;
	for (int i = sd->left; i < sd->right; ++i) {
	    GFREQ*	ff = *df++;
	    int	n = 0;
	    while (neogfq(ff)) {
		if ((ff++)->nres) ++n;
	    }
	    av += n;
	}
	av /= (sd->right - sd->left);
	return (av + 1.); 
}

/***********************************************************
	Management of 'dynamic' gap states
***********************************************************/

void printdlt(FILE* fd, const IDELTA* dlt)
{
    do {
	fprintf(fd, " %d %d", dlt->glen, dlt->nins);
    } while (neodelta(++dlt));
}

void warnovfmsg()
{
	if (warnrecovf) {
	    fputs("Warning: May not be optimal!\n", stderr);
	    if (warnrecovf & DGSRECOVF) 
		fputs("\tDGS records overflow!\n", stderr);
	    if (warnrecovf & CANRECOVF)
		fputs("\tCandidate records overflow!\n", stderr);
	}
}

/********************************************************
	g = SUM_j {C(j) * d(j)}
	C(j) = SUM_i>=j c(i)
*********************************************************/

/********************************************************
	evaluate the gap-cost by the gap-profile algorithm with
	gap-states synthesized from static and dynamic gap-state
	variables
*********************************************************/

VTYPE newgap(const GFREQ* cf, const GFREQ* df)
{
	VTYPE   g = 0;

	for ( ; neogfq(df); ++df) {
	    for ( ; neogfq(cf); ++cf) {
		if (cf->glen >= df->glen) break;
	    }
	    if (!neogfq(cf)) break;
	    g += cf->freq * df->freq;
	}
	return (g);
}

VTYPE newgap(const GFREQ* cf, const IDELTA* dlc, const GFREQ* df, const IDELTA* dld)
{
	VTYPE   g = 0;

	for ( ; neogfq(df); ++df) {
	    int	j = GapLenSD(df, dld);
	    for ( ; neogfq(cf); ++cf) {
		int	i = GapLenSD(cf, dlc);
		if (i >= j) break;
	    }
	    if (!neogfq(cf)) break;
	    g += cf->freq * df->freq;
	}
	return (g);
}

VTYPE newgap(const GFREQ* cf, const IDELTA* dlc, int j)
{
//	int	j = dld->nins;

	for ( ; neogfq(cf); ++cf) {
	    int	i = GapLenSD(cf, dlc);
	    if (i >= j) return (cf->freq);
	}
	return (0);
}

VTYPE newgap(const GFREQ* df, int i, const IDELTA* dld)
{
//	i = dlc->nins;
	VTYPE   g = 0;

	for ( ; neogfq(df); ++df) {
	    while (df->glen >= dld[1].glen) ++dld;
	    if (i < df->glen + dld->nins) break;
	    g += df->freq;
	}
	return (g);
}

/********************************************************
	copy and clear IDELTA
*********************************************************/

void copydelta(IDELTA* dlt, const IDELTA* dln)
{
    do {
	*dlt++ = *dln;
    } while (neodelta(++dln));
    *dlt = *dln;
}

void cleardelta(IDELTA* dlt)
{
    *dlt++ = ZeroDelta;
    *dlt   = LastDelta;
}

/********************************************************
	filter dynamic gap-state variable through a sequence position
	    represented by static gap-state 'df'
********************************************************/

void newdelta(IDELTA* dlt, const GFREQ* df, const IDELTA* dln, int n)
{
	IDELTA*	dst = dlt;
	IDELTA	tmp = ZeroDelta;	// buffer in case of dlt == dln
	for ( ; neogfq(df); ++df) {
	    if (df->glen >= dln->glen) {
		while (df->glen >= dln[1].glen) ++dln;
		if (dln->nins > tmp.nins) {
		    int	dlnins = dln->nins;	// buffer in case of dlt == dln
		    *dst++ = tmp;
		    tmp.nins = dlnins;
		    tmp.glen = df->glen + n;
		}
	    }
	}
	*dst++ = tmp;
	*dst = LastDelta;
}

/********************************************************
	dynamic gap-state after insertion of a new gap
********************************************************/

void incdelta(IDELTA* dlt, int n)
{
	 while (neodelta(dlt)) (dlt++)->nins += n;
}

void incdelta(IDELTA* dlt, const IDELTA* dln, int n)
{
	do {
	    *dlt = *dln;
	    (dlt++)->nins += n;
	} while (neodelta(++dln));
	*dlt = *dln;
}

