/*****************************************************************************
*
*	Ditection of difference in two closely related sequences.
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
#include "mseq.h"
#include "autocomp.h"
#include "vmf.h"

#define NEVSELP (INT_MAX / 8 * 7)

struct RECD {VTYPE val; long ptr;};

class AlnIdn {
	Seq**	seqs;
	Seq*	a;
	Seq*	b;
	WINDOW	wdw;
	PwdB*	pwd;
	Vmf*	vmf;
	VTYPE	uu;
	VTYPE	vv;
	void	initdist(VTYPE* dd);
	VTYPE	enddist(VTYPE* dd);
	void	InitInfMtx(RECD* dd);
	RECD	FinitInfMtx(RECD* dd);
	VTYPE	forwardA();
public:
	AlnIdn(Seq** _seqs, PwdB* _pwd);
	~AlnIdn() {delete vmf;}
	SKL*	align_2(VTYPE* scr);
	VTYPE	distanceA();
	void	output(SKL* slctd, VTYPE scr);
};

class AlnOut {
	Seq**	seqs;
	char*	decode;
	CHAR*	image[2];
	int	left[2];
	int	nbr[2];
	void	indicate();
	void	oneblock(int k);
	bool	diffblock();
	void	prnt_idn();
public:
	AlnOut(Seq** seqs, GAPS** gps);
	~AlnOut() {delete *image;};
};

static	const	int	BLANK = 0xa0;
static	const	char	spc10[] = "          ";

static	int	calcnbr(Seq* a, int gp);

void usage()
{
}

template <class seq_t>
void AlnServer<seq_t>::setparam(int level)
{
	if (level & 1) {
	    setstrip(QUERY);
	    promptin("Thr = (%4.2f%%) : ", &alprm.thr);
	}
	if (level & 2) {
	    setlpw(QUERY);
	    setoutmode(QUERY);
	}
}

void AlnOut::indicate()
{
	fputs(spc10, out_fd);
	for (int i = 0; i < OutPrm.lpw; ++i) {
	    bool	mismatch = image[0][i] != image[1][i];
	    putc(mismatch? '*': ' ', out_fd);
	}
	putc('\n', out_fd);
}

void AlnOut::oneblock(int k)
{
	char	pline[MAXL];
	CHAR*	p1 = image[k];
	char*	ps = pline;

	for (int i = 0; i < OutPrm.lpw; ++i, ++p1)
	    *ps++ = (*p1 == BLANK)? ' ': decode[*p1];
	*ps = '\0';
	if (left[k] >= nbr[k]) 
	    fprintf(out_fd, "%s%s\n", spc10, pline);
	else
	    fprintf(out_fd, "%8d  %s%6d\n", left[k] + 1, pline, nbr[k]);
}

bool AlnOut::diffblock()
{
	CHAR*	p1 = image[0];
	CHAR*	p2 = image[1];
	for (int i = 0; i < OutPrm.lpw; ++i, ++p1, ++p2) {
	    if (*p1 == BLANK || *p2 == BLANK) continue;
	    if (*p1 != *p2) return (true);
	}
	return (false);
}

void AlnOut::prnt_idn()
{
	if (diffblock()) {
	    putc('\n', out_fd);
	    oneblock(0);
	    indicate();
	    oneblock(1);
	}
}

static int calcnbr(Seq* a, int gp)
{
	register int	n = *a->nbr;
	register CHAR*	s = a->at(0);

	for ( ; gp--; ++s) if (IsntGap(*s)) ++n;
	return (n);
}

AlnOut::AlnOut(Seq** sqs, GAPS** gaps) : seqs(sqs)
{
	int	active = 2;
	GAPS*	gp[2];
	CHAR*	seq[2];
	CHAR*	wrk[2];
	*image = new CHAR[2 * OutPrm.lpw];
	decode = seqs[0]->code->decode;
	for (int j = 0; j < 2; ++j) {
	    nbr[j] = calcnbr(seqs[j], gaps[j]->gps);
	    if (j) image[j] = image[j-1] + OutPrm.lpw;
	    wrk[j] = seq[j] = seqs[j]->at(gaps[j]->gps);
	    gp[j] = gaps[j];
	    if (!gp[j]->gln) ++gp[j];
	}

	int	z = 0;
	do {
	    for (int j = 0; j < 2; j++) left[j] = nbr[j];
	    for (int clm = 0; clm < OutPrm.lpw; ++clm, ++z) {
		if (active) active = 2;
		for (int j = 0; j < 2; ++j) {
	loop:	    int	pos = gp[j]->gps - gaps[j]->gps;
		    int	gap = gp[j]->gln;
		    if (!active) image[j][clm] = BLANK;
		    else if (pos > z) {
			image[j][clm] = *wrk[j];
			if (IsntGap(*wrk[j])) ++nbr[j];
			++wrk[j];
		    } else if (pos + gap > z) {
			image[j][clm] = gap_code;
		    } else if (gap < 0) {
			image[j][clm] = gap_code;
			if (!--active && !clm--) return;
		    } else {
			gp[j]++;
			goto loop;
		    }
		}
	    }
	    if (out_fd) prnt_idn();
	} while (active);
}

void AlnIdn::output(SKL* slctd, VTYPE scr)
{
	Gsinfo	gsi;
	gsi.skl = slctd;
	gsi.scr = skl_rngS_ng(seqs, &gsi, pwd);
	GAPS*	gaps[2] = {0, 0};
	skl2gaps(gaps, slctd);
	toimage(gaps, 2);
	FSTAT&	fst = gsi.fstat;
	int 	span = (int) (fst.mch + fst.mmc + fst.unp);
	if (!gaps[0] || !span) return;
	fphseqs(seqs, 2);
	FTYPE	percent = 100. * fst.mch / span;
	fprintf(out_fd, "Dist = %4d, Cons = %3d, Repl = %3d,",
	    (int) scr, (int) fst.mch, (int)  fst.mmc);
	fprintf(out_fd, "  Gaps = %2d, Unpairs = %3d, (%6.2f %%)\n",
	    (int) fst.gap, (int) fst.unp, percent);
	if (gsi.scr > 0 && OutPrm.lpw)
	    AlnOut	ao(seqs, gaps);
	fputs("\n", out_fd);
	delete[] gaps[0];
	delete[] gaps[1];
}

void AlnIdn::initdist(VTYPE* dd)
{
	int	r = b->left - a->left;
	VTYPE	g = a->inex.exgl? 0: vv;
	VTYPE	u = a->inex.exgl? 0: uu;
	VTYPE*	d = dd + r;
	CHAR*	as = a->at(a->left);
	CHAR*	bs = b->at(b->left);

	for (*d = 0; ++r < wdw.up; ++bs)
	    *++d = (g += u);
	*++d = NEVSELP;

	g = b->inex.exgl? 0: vv;
	u = b->inex.exgl? 0: uu;
	r = b->left - a->left;
	d = dd + r;
	for ( ; --r > wdw.lw; ++as)
	    *--d = (g += u);
	*--d = NEVSELP;
}

VTYPE AlnIdn::enddist(VTYPE* dd)
{
	int	rr = b->right - a->right;
	VTYPE	f = NEVSELP;
	VTYPE	g = NEVSELP;
	VTYPE*	d = dd + rr;

	if (a->inex.exgr) {
	    int		r = wdw.lw + 1;
	    d = dd + r;
	    CHAR*	ss = b->at(a->right + r);
	    for ( ; r < rr; ++r, ++d) {
		if (*d < g) g = *d;
		++ss;
	    }
	}
	if (b->inex.exgr) {
	    int	r = wdw.up - 1;
	    d = dd + r;
	    CHAR*	ss = a->at(b->right - r);
	    for ( ; r > rr; --r, --d) {
		if (*d < f) f = *d;
		++ss;
	    }
	}
	if (f < *d) *d = f;
	if (g < *d) *d = g;
	return (*d);
}

VTYPE AlnIdn::distanceA()
{
	VTYPE*	dd = new VTYPE[wdw.width] - wdw.lw;
	VTYPE*	gg = new VTYPE[wdw.width] - wdw.lw;
	VTYPE*	g = gg + wdw.lw;
	for (int r = wdw.lw; r <= wdw.up; r++) *g++ = NEVSELP;
	initdist(dd);
	int	m = a->left;
	int	n1 = m + wdw.lw + 1;
	int	n2 = m + wdw.up;
	CHAR*	as = a->at(m);
	for ( ; m < a->right; ++m, ++as, ++n1, ++n2) {
	    int		n = max(n1, b->left);
	    int		n9 = min(n2, b->right);
	    CHAR*	bs = b->at(n);
	    int		r = n - m;
	    VTYPE*	d = dd + r;
	    VTYPE		f = NEVSELP;
	    g = gg + r;
	    for ( ; n < n9; ++n, ++d, ++g, ++bs) {
		VTYPE	x  = d[-1] + vv;
		f  = min(x, f) + uu;
		x  = d[1] + vv;
		*g = min(x, g[1]) + uu;
		VTYPE	nd = (f < *g)? f: *g;
		*d += (*as != *bs);
		if (nd < *d) *d = nd;
	    }
	}
	VTYPE	dist = enddist(dd);

	delete[] (dd + wdw.lw);
	delete[] (gg + wdw.lw);
#if MONITOR
	fprintf(stderr, "A: (%d-%d), Dist = %6.1f\n", 
	    a->inex.vect, b->inex.vect, (FTYPE) dist);
#endif
	return (dist / uu);
}

void AlnIdn::InitInfMtx(RECD* dd)
{
	int	r = b->left - a->left;
	VTYPE	g = a->inex.exgl? 0: vv;
	VTYPE	u = a->inex.exgl? 0: uu;
	CHAR*	as = a->at(a->left);
	CHAR*	bs = b->at(b->left);
	RECD*	d = dd + r;
	long	origin = vmf->add(a->left, b->left, 0);

	d->val = 0; d->ptr = origin;
	for ( ; ++r < wdw.up; ++bs) {
	    (++d)->val = (g += u);
	    d->ptr = origin;
	}
	(++d)->val = NEVSELP;

	r = b->left - a->left;
	d = dd + r;
	g = b->inex.exgl? 0: vv;
	u = b->inex.exgl? 0: uu;
	for ( ; --r > wdw.lw; ++as) {
	    (--d)->val = (g += u);
	    d->ptr = origin;
	}
	(--d)->val = NEVSELP;
}

RECD AlnIdn::FinitInfMtx(RECD* dd)
{
	int	rr = b->right - a->right;
	RECD	f = {NEVSELP, 0};
	RECD	g = f;
	RECD*	d = dd + rr;

	if (a->inex.exgr) {
	    int	r = wdw.lw + 1;
	    d = dd + r;
	    CHAR*	ss = b->at(a->right + r);
	    for ( ; r < rr; ++r, ++d) {
		if (d->val < g.val) g = *d;
		++ss;
	    }
	}
	if (b->inex.exgr) {
	    int	r = wdw.up - 1;
	    d = dd + r;
	    CHAR*	ss = a->at(b->right - r);
	    for ( ; r > rr; --r, --d) {
		if (d->val < f.val) f = *d;
		++ss;
	    }
	}
	if (f.val < d->val) *d = f;
	if (g.val < d->val) *d = g;
	return (*d);
}

VTYPE AlnIdn::forwardA()
{
	CHAR*	dir = new CHAR[wdw.width] - wdw.lw;
	RECD*	dd = new RECD[wdw.width] - wdw.lw;
	RECD*	gg = new RECD[wdw.width] - wdw.lw;
	for (int r = wdw.lw; r <= wdw.up; r++) gg[r].val = NEVSELP;
	vmf->add(0, 0, 0);	// Skip 0-th record
	InitInfMtx(dd);			// Initialize corners

	int	m = a->left;
	CHAR*	as = a->at(m);
	int	n1 = m + wdw.lw + 1;
	int	n2 = m + wdw.up;
	for ( ; m < a->right; ++m, ++as, ++n1, ++n2) {
	    int	n = max(n1, b->left);
	    int	n9 = min(n2, b->right);
	    int	r = n - m;
	    CHAR*	bs = b->at(n);
	    RECD*	d = dd + r;
	    RECD*	g = gg + r;
	    CHAR*	e = dir + r;
	    RECD	f = {NEVSELP, 0};
	    for (int r = n - m; n < n9; ++r, ++n, ++d, ++g, ++e, ++bs) {
		VTYPE	x = d[-1].val + vv;
		if (x < f.val) {
		    f.val = x;
		    f.ptr = (d-1)->ptr;
		}
		f.val += uu;
		x  = d[1].val + vv;
		if (x < (g+1)->val) {
		    g->val = x;
		    g->ptr = (d+1)->ptr;
		} else {
		    g->val = (g+1)->val;
		    g->ptr = (g+1)->ptr;
		}
		g->val += uu;
		RECD	nd = (f.val < g->val)? f: *g;
		d->val += (*as != *bs);
		if (nd.val < d->val) {
		    *d = nd;
		    *e = 0;
		} else if (!*e) {
		    d->ptr = vmf->add(m, n, d->ptr);
		    *e = 1;
		}
#if FDEBUG
	    printf("%2d %2d %4.1lf %4.1lf %4.1lf\n", m, n, 
	(double) f.val, (double) g->val, (double) d->val);
#endif
	    }
	}
	RECD	nd = FinitInfMtx(dd);
	nd.ptr = vmf->add(a->right, b->right, nd.ptr);

#if MONITOR
	fprintf(stderr, "A: (%d-%d), Dist = %6.1f, Records = %ld\n",
	    a->inex.vect, b->inex.vect, (FTYPE) nd.val, *pp);
#endif
	
	delete[] (dir + wdw.lw);
	delete[] (dd + wdw.lw);
	delete[] (gg + wdw.lw);
	return (nd.val / uu);
}

SKL* AlnIdn::align_2(VTYPE* scr)
{
	*scr = forwardA();
	SKL*	skl = (*scr == VABORT)? 0: vmf->traceback(-1L);
	vreverse(skl + 1, skl->n);
	return (skl);
}

AlnIdn::AlnIdn(Seq* _seqs[], PwdB* prm) : seqs(_seqs), pwd(prm)
{
	stripe(seqs, &wdw, alprm.sh);
	a = seqs[0];
	b = seqs[1];
	vv = (VTYPE) alprm.v;
	uu = (VTYPE) alprm.u;
	vmf = new Vmf();
}

int idn_main(CalcServer<Seq>* svr, Seq* seqs[], ThQueue<Seq>* q = 0)
{
	Seq*&	a = seqs[0];
	Seq*&	b = seqs[1];

	if (a->many != 1 || b->many != 1)
	    fatal("Unable to compare MSAs !\n");
	int 	cut = int((a->right - a->left + b->right - b->left) * alprm.thr / 100);
	if (algmode.nsa) {	// alignment
	    AlnIdn	aln(seqs, (PwdB*) svr->prm);
	    VTYPE	score;
	    SKL* skl = aln.align_2(&score);
	    if (!skl) return (ERROR);
	    if (abs(int(score)) > cut || algmode.nsa & 1) {
		if (algmode.lcl) skl = trimskl(seqs, skl);
		m_thread_Lock(q);
		aln.output(skl, score);
		m_thread_Unlock(q);
	    }
	} else {		// score only
	    int lendiff = a->right - a->left - b->right + b->left;
	    if (abs(lendiff) > cut) return (OK); // Big diff in lengths
	    AlnIdn	aln(seqs, (PwdB*) svr->prm);
	    VTYPE	score = aln.distanceA();
	    if (score < cut) {
		m_thread_Lock(q);
		fprintf(out_fd, "%-12s %-12s %3d\n", 
		    a->sqname(), b->sqname(), int(score));
		m_thread_Unlock(q);
	    }
	}
	return (OK);
}

template <class seq_t>
int AlnServer<seq_t>::localoption(int& argc, const char**& argv)
{
	return 0;
}

void idn_setup(CalcServer<Seq>* svr)
{
	setup_output(algmode.nsa);
	prePwd(svr->in_face);
	PwdB*	pwd = new PwdB(svr->in_face);
	svr->prm = (void*) pwd;
}

void idn_cleanup(CalcServer<Seq>* svr)
{
	delete (PwdB*) svr->prm;
}

static void setdefparam()
{
	alprm.sh = 2;
	alprm.thr = 1;
	alprm.u = alprm.v = 1;
}

int main(int argc, const char** argv)
{
	setdefparam();
	AlnServer<Seq>	svr(argc, argv, IM_NONE, IM_ALTR, (void*) 0, 
	&idn_main, &idn_setup, &idn_cleanup);
	if (svr.autocomp() == 1) svr.menucomp();
	EraDbsDt();
	eraStrPhrases();
	return (0);
}
