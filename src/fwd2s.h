/*****************************************************************************
*
*	Alignment of two nucleotide sequences 
*	5' and 3' splice site signals are considered.
*	Assumes there are internal gap(s) in the reference nucleotide sequence.
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

#ifndef  _FWD2S_
#define  _FWD2S_

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

static	const	int	NCAND_S = 4;

template <class recd_t>
class Fwd2s {
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
	FTYPE	u2divu1;
	FTYPE	v2divv1;
	int	Nod;
	SpJunc*	spjcs;
	recd_t*	recbuf;
	recd_t*	hh[NOL];	// H, G matrix
	recd_t*	hl;	// [DIA, HORI, HORL][candidates]
	recd_t*	f1, *f2;
	void	initializeS(int recsize);
	void	initS();
	recd_t*	lastS();
	VTYPE	gapopen(const recd_t* rcd, const mSeqItr& asi, int d3);
	VTYPE	ins_penalty(const recd_t* rcd, const mSeqItr& asi);
	VTYPE	del_penalty(const recd_t* rcd, const mSeqItr& asi);
	void	update(recd_t* rcd, const mSeqItr& asi, int d3);
	void	update(recd_t* dst, const recd_t* src, const mSeqItr& asi, 
		VTYPE gop, int d3);
	mSeqItr	blank_Seqitr;
public:
	Fwd2s(mSeq** _seqs, const PwdM* _pwd, const WINDOW* _wdw = 0);
	~Fwd2s();
	VTYPE	forwardS(long pp[]);
	SKL*	globalS(VTYPE* scr);
	VTYPE	verify(Gsinfo* gsi);
};

template <class recd_t>
Fwd2s<recd_t>::Fwd2s(mSeq* _seqs[], const PwdM* _pwd, const WINDOW* _wdw) :
	seqs(_seqs), a(seqs[0]), b(seqs[1]), 
	asi(a, a->left), bsi(b, b->left),
	wdw(_wdw), pwd(_pwd), vmf(0), mfd(0), gfab(0),
	abt(a->byte), bbt(b->many),
	u2divu1(pwd->BasicGEP < 0? (FTYPE) pwd->LongGEP / pwd->BasicGEP: 0),
	v2divv1(pwd->BasicGOP < 0? (FTYPE) pwd->LongGOP / pwd->BasicGOP: 0),
	Nod(2 * pwd->Noll - 1)
{
	SeqThk	bthk = {b->sumwt, 0, b->sumwt};
	a->setsimmtx(pwd->simmtx);
	bzsi.dns = &bthk;
	spjcs = new SpJunc(b);
	size_t	recsize = wdw? wdw->width * pwd->Nrow + NCAND_S + pwd->Noll: 3;
	recbuf = new recd_t[recsize];
	initializeS(recsize);
	if (!wdw) return;
	recd_t*	h = recbuf - wdw->lw + 1;
	for (int i = 0; i < pwd->Nrow; ++i, h += wdw->width) hh[i] = h;
	hl = h += (wdw->lw - 1); h += (NCAND_S + 1);
	f1 = h++;
	f2 = (pwd->Noll == 3)? h: 0;
}

template <class recd_t>
Fwd2s<recd_t>::~Fwd2s()
{
	delete vmf;
	delete[] recbuf;
	delete[] gfab;
	delete mfd;
	delete spjcs;
}

template <class recd_t>
void Fwd2s<recd_t>::initS()
{
	int	m = a->left;
	int	n = b->left;
	int	r = n - m;
	int	rr = b->right - m;
	int	dir = a->inex.exgl? DEAD: DIAG;
	recd_t*	h = hh[0] + r;

	h->val = 0;
	h->dir = dir;
	h->ptr = vmf? vmf->add(a->left, n, 0): 0;
	h->jnc = n;
	if (a->inex.exgl) {		// semi-global
	    if (wdw->up < rr) rr = wdw->up;
	    while (++r <= rr) {
		++h;
		h->val = 0;
		h->dir = DIAG;
		h->glb = h->jnc = ++n;
		h->ptr = vmf? vmf->add(a->left, n, 0): 0;
	    }
	}

	mSeqItr	asi(a, m - 1);
	r = b->left - a->left;
	rr = b->left - a->right;
	if (wdw->lw> rr) rr = wdw->lw;
	h = hh[0] + r;
	while (--r >= rr) {
	    --h; ++m; ++asi;
	    if (b->inex.exgl) {	// semi global
		h->val = 0;
		h->dir = DEAD;
		h->jnc = b->left;
		h->ptr = vmf? vmf->add(m, b->left, 0L): 0L;
	    } else {
	  	VTYPE	pua = (pwd->*pwd->unpa)(asi, bzsi);
	  	VTYPE	gnp = gapopen(h + 1, asi, 1);
	 	gnp = (m - a->left < pwd->codonk1)? gnp + pua:
		    (VTYPE) (v2divv1 * gnp + u2divu1 * pua);
		update(h, h + 1, asi, gnp, 1);
		h->dir = VERT;
	    }
	}
}

template <class recd_t>
recd_t* Fwd2s<recd_t>::lastS()
{
	recd_t*	h9 = hh[0] + b->right - a->right;
	recd_t*	mx = h9;

	if (b->inex.exgr) {
	    int	rw = wdw->up;
	    int	rf = b->right - a->left;
	    if (rf < rw) rw = rf;
	    for (recd_t* h = hh[0] + rw; h > h9; --h)
		if (h->val > mx->val) mx = h;
	}

	int	rw = wdw->lw;
	int	rf = b->left - a->right;
	if (rf > rw) rw = rf;
	else	 rf = rw;
	if (a->inex.exgr)
	    for (recd_t* h = hh[0] + rw; h <= h9; ++h)
		if (h->val > mx->val) mx = h;
	int	i = mx - h9;
	rf = a->right;	// m9
	rw = b->right;	// n9
	if (i > 0) rf -= i;
	if (i < 0) rw += i;
	if (vmf) mx->ptr = vmf->add(rf, rw, mx->ptr);
	return (mx);
}

template <class recd_t>
VTYPE Fwd2s<recd_t>::forwardS(long* pp)
{
	if (a->isprotein() || b->isprotein())
	    fatal("Inproper combination of seq. types !\n");
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
	recd_t*	hf[NOD] = {0, f1, 0, f2, 0};
	int	nx[NCAND_S + 1];	// [DIA, HORI, HORL][candidates]
	recd_t*	maxphl[NOD];	// 
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0L};
	bool	Local = algmode.lcl & 16;
	bool	LocalL = Local && a->inex.exgl && b->inex.exgl;
	bool	LocalR = Local && a->inex.exgr && b->inex.exgr;

	if (vmf) vmf->add(0, 0, 0L);	// Skip 0-th record
	initS();

	int	m = a->left;
	PfqItr	api(a, m);
	int	api_size = api.size();
	if (!a->inex.exgl) --m;		// global
	int	n1 = m + wdw->lw - 1;
	int	n2 = m + wdw->up;
	mSeqItr apsi(a, m);
	for (mSeqItr asi = apsi; ++m <= a->right; asi = apsi) {
	    ++apsi;
	    bool	internal = !a->inex.exgr || m < a->right;
	    ++n1; ++n2;
	    int	n = max(n1, b->left);
	    int	n9 = min(n2, b->right);
	    int	r = n - m;
	    recd_t*&	h = hf[0] = hh[0] + r;
	    recd_t*&	g = hf[2] = hh[1] + r;
	    recd_t*&	g2 = hf[4] = dagp? hh[2] + r: 0;
const	    VTYPE*	qprof = asi.prof();	// sim2(asi, .)
	    f1->reseth(); if (dagp) f2->reseth();
	    for (int l = 0; l <= NCAND_S; ++l) {
		hl[l].reseth();
		nx[l] = l;
	    }
	    VTYPE	pua = internal? (pwd->*pwd->unpa)(asi, bzsi): 0;
	    int	ncand = 0;
#if F2DEBUG
	    if (OutPrm.debug) {
		printf("%2d %2d %2d", m, n, h->dir);
		putvar(h->val); putchar('\n');
	    }
#endif
	    bool	a_in_zone = api_size && api.eq(m);
	    bsi.reset(n);
	    for ( ; ++n <= n9; ++bsi) {
		++h; ++g; if (dagp) ++g2;
		VTYPE	x, y;

//	Diagonal
		VTYPE	gop = 0;
		VTYPE	gnp = 0;
		recd_t*	from = h;
		recd_t*	mx = h;
		VTYPE	diag = h->val;
		if (m == a->left) goto Horizon;	// initialize column
		gop = gapopen(h, asi, 0);
		update(h, h, asi, qprof[*bsi.res] + gop, 0);
		h->dir = isdiag(h)? DIAG: NEWD;

//	normal deletion
		++from;
		gop = gapopen(from, asi, 1);
		gnp = gapopen(g + 1, asi, 1);
		if (isntvert(from) &&
		    (from->val + gop > g[1].val + gnp))
			update(g, from, asi, gop, 1);
		else	update(g, g + 1, asi, gnp, 1);
		g->val += pua;
		g->dir = VERT;
		if (vmf && g->dir & SPIN)
		    g->ptr = vmf->add(m - 1, n, g->ptr);
		if (g->val > mx->val) mx = g;

//	long deletion
		if (dagp) {
		    gnp = (VTYPE) (v2divv1 * gapopen(g2 + 1, asi, 1));
		    gop = (VTYPE) (v2divv1 * gop);
		    if ((isntvert(from) && 
			(from->val + gop > g2[1].val + gnp)))
			update(g2, from, asi, gop, 1);
		    else
			update(g2, g2 + 1, asi, gnp, 1);
		    g2->val += (VTYPE) (u2divu1 * pua);
		    g2->dir = VERL;
		    if (vmf && g2->dir & SPIN)
			g2->ptr = vmf->add(m - 1, n, g2->ptr);
		    if (g2->val > mx->val) mx = g2;
		}
Horizon:
//	nomal insertion
		from = h - 1;
		gop = gapopen(from, asi, -1);
		if (isnthori(from) &&
		    (from->val + gop > f1->val))
			update(f1, from, asi, gop, -1);
		else	update(f1, f1, asi, 0, -1);
		f1->val += pwd->BasicGEP;
		f1->dir = (f1->dir & SPIN) + HORI;
		if (f1->val >= mx->val) mx = f1;
//	long insertion
		if (dagp) {
		    gop = (VTYPE) (v2divv1 * gop);
		    if ((isnthori(from) &&
			(from->val + gop > f2->val)))
				update(f2, from, asi, gop, -1);
		    else	update(f2, f2, asi, 0, -1);
		    f2->val += pwd->LongGEP;
		    f2->dir = (f2->dir & SPIN) + HORL;
		    if (f2->val > mx->val) mx = f2;
		}

//	3' boundary, assume no overlapping signals
		if (internal && b->exin->isAccpt(n)) {
		    VTYPE	sigJ = a_in_zone? api.match_score(m): 0;
		    vclear(maxphl, Nod);
		    for (int l = 0; l < ncand; ++l) {
			recd_t*	phl = hl + nx[l];
			x = phl->val + sigJ + 
				pwd->IntPen->Penalty(n - phl->jnc) +
				b->exin->sig53(phl->jnc, n, IE53);
			from = hf[phl->dir];
			if (x > from->val) {
			    from->val = x;
			    maxphl[phl->dir] = phl;
			}
		    }
		    for (int d = 0; d < Nod; ++d) {
			recd_t*	phl = maxphl[d];
			if (!phl) continue;
			from = hf[d];
			if (vmf) from->ptr = vmf->add(m, n, 
			    vmf->add(m, phl->jnc, phl->ptr));
			from->jnc = n;
			from->dir |= SPJCI;
			if (from->val > mx->val) mx = from;
		    }
		}

//	Find optimal path
		if (mx != h) *h = *mx;
		else if (Local && h->val > diag) {
		    if (LocalL && diag == 0 && !(h->dir & SPJC))
			h->ptr = vmf? vmf->add(m - 1, n - 1, 0): 0;
		    else if (LocalR && h->val > maxh.val) {
			maxh.val = h->val;
			maxh.p = h->ptr;
			maxh.m = m;
			maxh.n = n;
		    }
		}
		if (LocalL && h->val <= 0) h->reseth();
		else if (vmf && h->dir == NEWD)
		    h->ptr = vmf->add(m - 1, n - 1, h->ptr);

//	5' boundary
		if (internal && b->exin->isDonor(n)) { 
		    VTYPE	sigJ = b->exin->sig53(n, 0, IE5);
		    int		hd = dir2nod[mx->dir & 15];
		    for (int k = hd == 0? 0: 1; k < Nod; ++k) {
                               // An orphan exon is disallowed
			from = hf[k];
			if (!from->dir || (from->dir & SPIN)) continue;
			if (k != hd && hd >= 0) {
			    y = mx->val;
			    if (hd == 0 || (k - hd) % 2) y += pwd->GOP[k / 2];
			    if (from->val <= y) continue;	// prune
			}
			x = from->val + sigJ;
			int	l = ncand < NCAND_S? ++ncand: NCAND_S;
			while (--l >= 0) {
			    if (x > hl[nx[l]].val)
				swap(nx[l], nx[l + 1]);
			    else
				break;
			}
			if (++l < INTR) {
			    recd_t*	phl = hl + nx[l];
			    *phl = *from;
			    phl->val = x;
			    phl->jnc = n;
			    phl->dir = k;
			} else --ncand;
		    }
		}
#if F2DEBUG
		if (OutPrm.debug) {
		    printf("%2d %2d %2d ", m, n, mx->dir);
		    putvar(mx->val); putvar(diag); 
		    putvar(g->val); putvar(f1->val);
		    if (dagp) {
			putvar(g2->val); putvar(f2->val);
		    }
		    if (algmode.lsg) {
			putvar(hl[nx[0]].val);
			putvar(hl[nx[1]].val);
		    }
		    putchar('\n');
		}
#endif
	    }	// end of n-loop
	    if (a_in_zone) ++api;	// has exon-exon junctions
	}	// end of m-loop

	if (LocalR) {
	    if (pp) *pp = vmf->add(maxh.m, maxh.n, maxh.p);
	} else {
	    recd_t*	mx = lastS();
	    maxh.val = mx->val;
	    if (pp) *pp = mx->ptr;
	}
	return (maxh.val);
}

template <class recd_t>
VTYPE Fwd2s<recd_t>::verify(Gsinfo* gsi)
{
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
	EISCR	rbuf;
	FSTAT*  fst = &gsi->fstat;
	FSTAT   pst;
	Eijnc*	eijnc = gsi->eijnc = new Eijnc(true);
	Cigar*	cigar = gsi->cigar = 0;
	Vulgar*	vlgar = gsi->vlgar = 0;
	Samfmt*	samfm = gsi->samfm = 0;
	switch (algmode.nsa) {
	    case CIG_FORM: cigar = gsi->cigar = new Cigar(); break;
	    case VLG_FORM: vlgar = gsi->vlgar = new Vulgar(); break;
	    case SAM_FORM: samfm = gsi->samfm = new Samfmt(); break;
	    default: break;
	}
	vclear(fst);
	vclear(&pst);
	vclear(&rbuf);
	gsi->noeij = 0;
	if (wsk[1].n == wsk->n && b->inex.exgl) {++wsk; --num;}
	int	m = wsk->m;
	int	n = wsk->n;
	asi.reset(m);
	bsi.reset(n);
	PfqItr	api(a, m);
	int	api_size = api.size();
	bool	usespb = api_size && use_spb();
	SeqThk	bthk = {b->sumwt, 0, b->sumwt};
	mSeqItr	bzsi;
	bzsi.dns = &bthk;

	h->val = 0;
	rbuf.left = n;
	rbuf.rleft = m;
	rbuf.iscr = NEVSEL;
	rbuf.sig3 = sig3;
	if (cigar && m) cigar->push('H', m);	// local alignment
	if (samfm) {
	    samfm->left = m;
	    if (m) samfm->push('H', m);		// local alignment
	    if (b->inex.sens == 0) samfm->pos = n;
	}
	while (--num) {
	    ++wsk;
	    int	mi = wsk->m - m;
	    if (insert && mi) {			// end of insertion
		h->val += pwd->UnpPenalty(insert);
		if (hi && insert > intlen)	// pre-intron gap
		    hi->val += pwd->UnpPenalty(insert - intlen);
		if (hi && hi->val >= h->val) { 	// intron
		    if (cigar) {
			if (preint) cigar->push('D', preint);
			cigar->push('N', intlen);
		    }
		    if (samfm) {
			if (preint) samfm->push('D', preint);
			samfm->push('N', intlen);
		    }
		    if (vlgar) {
			if (preint) vlgar->push('G', 0, preint);
			vlgar->push('5', 0, 2);
			vlgar->push('I', 0, intlen - 4);
			vlgar->push('3', 0, 2);
		    }
		    hb = ha;
		    if (eijnc && rbuf.right - rbuf.left > 1) eijnc->push(&rbuf);
		    rbuf.left = rbuf.right + intlen;
		    rbuf.rleft = m;
		    rbuf.sig3 = sig3;
		    rbuf.iscr = NEVSEL;
		    *h = *hi;
		    hi = 0;
		    insert -= (preint + intlen);
		}
		if (insert) {
		    if (cigar) cigar->push('D', insert);
		    if (samfm) samfm->push('D', insert);
		    if (vlgar) vlgar->push('G', 0, insert);
		    insert = intlen = preint = 0;
		}
	    }
	    int	ni = wsk->n - n;
	    if (ni && deletn) {		// end of delettion
		if (vlgar) vlgar->push('G', deletn, 0);
		deletn = 0;
	    }
	    int	i = mi - ni;
	    int	d = (i >= 0)? ni: mi;
	    if (d) {
		if (cigar) cigar->push('M', d);
		if (vlgar) vlgar->push('M', d, d);
		if (samfm) samfm->push('M', d);
		m += d;
		n += d;
		if (a->many > 1) {
		  for ( ; d; --d, ++asi, ++bsi) {
		    fst->gap += (bsi.dullend()? 0: gapopen(h, asi, 0));
		    h->val += (pwd->*pwd->sim2)(asi, bsi);
		    (pwd->*pwd->stt2)(fst, asi, bsi);
		    if (eijnc) eijnc->shift(rbuf, *fst, pst,
			n - rbuf.left == alprm2.jneibr);
		    update(h, asi, 0);
		  }
		} else {
		  int	run = 0;
		  for ( ; d; --d, ++asi, ++bsi) {
		    h->val += (pwd->*pwd->sim2)(asi, bsi);
		    (pwd->*pwd->stt2)(fst, asi, bsi);
		    if (*asi.res == *bsi.res) {
		      if (run < 0) {
			if (samfm) samfm->push('X', -run);
			run = 0;
		      }
		      ++run;
		    } else {
		      if (run > 0) {
			if (samfm) samfm->push('=', run);
			run = 0;
		      }
		      --run;
		    }
		    if (eijnc) eijnc->shift(rbuf, *fst, pst,
			n - rbuf.left == alprm2.jneibr);
		    if (samfm) {
		      if (run > 0) samfm->push('=', run); else
		      if (run < 0) samfm->push('X', -run);
		    }
		  }
		}
	    }
	    if (i > 0) {		// deletion
		for (int j = 0; j < i; ++j, ++asi) {
		    VTYPE gop = bsi.dullend()? 0: gapopen(h, asi, 1);
		    fst->gap += gop;
		    (pwd->*pwd->stt2)(fst, asi, bzsi);
		    h->val += gop + (pwd->*pwd->unpa)(asi, bzsi);
		    update(h, asi, 1);
		}
		deletn += i;
		if (cigar) cigar->push('I', i);
		if (samfm) samfm->push('I', i);
		if (vlgar) vlgar->push('G', i, 0);
	    } else if (i < 0) {		// insertion
		int	n3 = n + (i = -i);
		VTYPE	xi = NEVSEL;
		if (!hi && i >= IntronPrm.llmt) {	// intron?
		    sig5 = b->exin->sig53(n, n3, IE5);
		    sig3 = b->exin->sig53(n, n3, IE3);
		    xi = b->exin->sig53(n, n3, IE5P3)
			+ pwd->IntPen->Penalty(i);
		    if (usespb) xi += api.match_score(m);
		}
		if (xi > pwd->GapPenalty(i) && xi > rbuf.iscr) {
		    preint = insert;
		    intlen = i;
		    rbuf.right = n;
		    rbuf.rright = m;
		    rbuf.iscr = xi;
		    rbuf.escr = h->val + sig5 - hb;
		    rbuf.sig5 = sig5;
		    hi = recbuf + 1;
		    *hi = *h;
		    hi->val += xi;
		    ha = h->val + xi - sig3;
		    if (eijnc) {
			eijnc->store(rbuf, *fst, pst,
			n - rbuf.left < alprm2.jneibr);
			pst = *fst;
		    }
		} else if (!a->inex.exgl || m != a->left) {
		    VTYPE gop = (bsi.dullend()? 0: gapopen(h, asi, -1));
		    fst->gap += gop;
		    fst->unp += i * a->sumwt;
		    h->val += gop;
		    update(h, asi, -i);
		}
		bsi += i;
		insert += i;
	    }
	    m = wsk->m;
	    n = wsk->n;
	    if (usespb) while (!api.end() && api.lt(m)) ++api;
	}
	if (insert && !(a->inex.exgr && m == a->right)) {
	    if (cigar) cigar->push('D', insert);
	    if (samfm) samfm->push('D', insert);
	    if (vlgar) vlgar->push('G', 0, insert);
	}
	if (deletn && !(b->inex.exgr && n == b->right)) {
	    if (cigar) cigar->push('I', deletn);
	    if (samfm) samfm->push('I', deletn);
	    if (vlgar) vlgar->push('G', deletn, 0);
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
	    rbuf.unp = (int) (fst->unp - pst.unp);
	    eijnc->push(&rbuf);
	    rbuf.left = endrng.left;
	    rbuf.right = endrng.right;
	    eijnc->push(&rbuf);
	    eijnc->flush();
	    gsi->noeij = eijnc->size() - 1;
	}
	if (cigar) cigar->flush();
	if (samfm) {
	    samfm->right = m;
	    if (m < a->len) samfm->push('H', a->len - m);
	    samfm->flush();
	    samfm->mapq = 30 + int(100 * (fst->mmc + fst->unp) / a->sumwt / a->len);
	    if (b->inex.sens) {
		samfm->flag |= 0x10;
		samfm->pos = n - 1;
		vreverse(samfm->rec, samfm->size());
	    }
	}
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
	return (gsi->scr = h->val);
}

template <class recd_t>
SKL* Fwd2s<recd_t>::globalS(VTYPE* scr)
{
	mfd = new Mfile(sizeof(SKL));
	SKL	wsk;
	mfd->write((UPTR) &wsk);	// dummy call
	vmf = new Vmf();
	long	ptr;
	*scr = forwardS(&ptr);
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
	return stdskl(&skl);
}

// interface 

template <class recd_t>
VTYPE HomScoreS(mSeq* seqs[], const PwdM* pwd)
{
	WINDOW	wdw;

	stripe((const Seq**) seqs, &wdw, alprm.sh);
	Fwd2s<recd_t> alnv(seqs, pwd, &wdw);
	return alnv.forwardS(0);
}

template <class recd_t>
SKL* alignS(mSeq* seqs[], const PwdM* pwd, VTYPE* scr)
{
	WINDOW  wdw;

	stripe((const Seq**) seqs, &wdw, alprm.sh);	// Recacl. window boundaries
	Fwd2s<recd_t> alnv(seqs, pwd, &wdw);
	return (alnv.globalS(scr));
}

template <class recd_t>
VTYPE skl_rngS(mSeq* seqs[], Gsinfo* gsi, const PwdM* pwd)
{
	Fwd2s<recd_t> fwd2s(seqs, pwd);
	return fwd2s.verify(gsi);
}

#endif // _FWD2S_
