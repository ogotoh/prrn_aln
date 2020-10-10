/*****************************************************************************
*
*	mseq.c:	multiple sequence
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
*	Osamu Gotoh, Ph.D.      (2003-2012)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#include "mseq.h"

mSeq* mSeq::clean()
{
	if (!vrtl || !wasmst.thk) {
	    if (thk)	delete[] (thk - 1);
	    if (dithk)	delete[] (dithk -  1);
	    delete[] internalres;
	}
	if (!vrtl || !wasmst.gfq) delete gfq;
	if (!vrtl || !wasmst.psq) delete[] pseq;
	if (!vrtl || !wasmst.cmp) {delete[] cmps; cmps = 0;}

#if SSHP
	if (!vrtl || !wasmst.snd) delete[] sshp_prof;
#endif
	clear();
	return this;
}

mSeq::~mSeq()
{
	refresh();
}

mSeq* mSeq::refresh(const int& m, const int& l)
{
	clean();
	Seq::refresh(m, l);
	return (this);
}

// alias of this

mSeq* mSeq::aliaseq(mSeq* dest, bool this_is_alias)
{
	if (!dest)	dest = new mSeq;
	else	dest->refresh();
	int	destid = dest->sid;
//	*dest = *this;
	memcpy(dest, this, sizeof(mSeq));
	dest->jxt = 0;
	if (this_is_alias) {
	    vrtl = sid;
	    wasmst.thk = (bool) thk;
	    wasmst.gfq = (bool) gfq;
	    wasmst.psq = (bool) pseq;
	    wasmst.cmp = (bool) cmps;
#if SSHP
	    wasmst.snd = (bool) sshp_prof;
#endif
#if USE_WEIGHT
	    inex.vtwt = !weight;
#endif
	} else {
	    dest->vrtl = destid;
	    dest->wasmst.thk = (bool) dest->thk;
	    dest->wasmst.gfq = (bool) dest->gfq;
	    dest->wasmst.psq = (bool) dest->pseq;
	    dest->wasmst.cmp = (bool) dest->cmps;
#if SSHP
	    dest->wasmst.snd = (bool) dest->sshp_prof;
#endif
#if USE_WEIGHT
	    dest->inex.vtwt = !dest->weight;
#endif
	}
	return (dest);
}

mSeq* mSeq::getseq(const char* str, DbsDt* dbf)
{
	if (!Seq::getseq(str, dbf)) return (0);
	clean();
	return (this);
}

mSeq* mSeq::copyseq(mSeq* dst, int srl)
{
	if (dst) dst->clean();
	else {
	    dst = new mSeq(many, len);
	    dst->clear();
	}
	Seq::copyseq(dst, srl);
	dst->mkthick();
	return (dst);
}

mSeq* mSeq::catseq(mSeq* head)
{
	clean();
	if (!head) head = new mSeq(many, 0);
	Seq::catseq(head);
	return (head);
}

mSeq* mSeq::cutseq(mSeq* dest, int snl)
{
	if (dest) dest->clean();
	else	dest = new mSeq(many, len);
	Seq::cutseq(dest, snl);
	return (dest);
}

mSeq* mSeq::extseq(mSeq* dst, int* grp, int srl, FTYPE nfact)
{
	if (dst) dst->clean();
	else	dst = new mSeq;
	Seq::extseq(dst, grp, srl, nfact);
	return (dst);
}

mSeq* mSeq::extract(Seq* src, int* grp, int srl)
{
	clean();
	src->extseq(this, grp, srl);
	return (this);
}

mSeq* mSeq::read_dbseq(DbsDt* dbf, long pos)
{
	clean();
	if (!Seq::read_dbseq(dbf, pos)) return (0);
	return (this);
}

void mSeq::mkthick()
{
	if (thk) return;
	bool	dd = algmode.qck & 1;
#if USE_WEIGHT
	bool	vwt1 = weight && unit_mode();
	if (vwt1 && sumwt == 0) {
	    for (int i = 0; i < many; ++i)
		sumwt += weight[i];
	} else	
#endif
	if (sumwt == 0)	sumwt = many;
#if SSHP
	if (many) evenwt = 1. / many;
#endif
	VTYPE	wt1 = 1;
	FTYPE	ltgapf = inex.exgl? 0: alprm.tgapf;
	FTYPE	rtgapf = inex.exgr? 0: alprm.tgapf;
	internalres = (inex.dels || inex.nils)? new CHAR*[many]: 0;
	int	thk_len = inex.dels? len: (internalres? 2: 0);
	SeqThk	bdy = {sumwt, 0, sumwt};
	SeqThk*	rthk = new SeqThk[thk_len + 2];
	if (dd) {
	    dithk = new DiThk[thk_len + 2];
	    vclear(dithk, thk_len + 2);
	}
	if (thk_len == 0) {				// no gap, full term gap 
	    *rthk = rthk[1] = bdy;
	    thk = ++rthk;
	    if (dd) {
		dithk->cc = dithk->cx = dithk[1].cc = dithk[1].cx = sumwt;	// **
		++dithk;
	    }
	    return;
	}
	CHAR*	ls = at(0);
	CHAR*	rs = at(len);
	CHAR*	fs = ls;
	CHAR*	bs = fs - many;
	bool*	egap = new bool[many];

	rthk->cfq = 0;				// thk[-1]
	rthk->dfq = rthk->efq = sumwt * ltgapf;
	thk = ++rthk;
	thk[thk_len].cfq = thk[thk_len].efq = 0;
	thk[thk_len].dfq = sumwt * rtgapf;	// thk[len]
	DiThk*	dwk = 0;
	if (dd) {
	    dithk->cc = dithk->cx = sumwt * ltgapf;
	    dwk = ++dithk;
	    dithk[thk_len].cc = dithk[thk_len].cx = sumwt * rtgapf;
	}

// find right-most left end gap

	if (ltgapf < 1.) {	// fill left end gap
	    vset(egap, true, many);
	    while (fs < rs) {
		VTYPE	w0 = 0;	// ..
		VTYPE	w1 = 0;	// @-	@ = [-*]
		VTYPE	w2 = 0;	// ?*	? = [.-*]
		VTYPE	w3 = 0;	// .*
		int	c = 0;
		for (int i = 0; i < many; ++i, ++fs, ++bs) {
#if USE_WEIGHT
		    if (vwt1) wt1 = weight[i];
#endif
		    if (*fs == nil_code) {	// ..
			w0 += wt1;
			if (dd) dwk->dd += wt1 * ltgapf;
			++c;
		    } else {
			if (egap[i]) {		// .*
			    egap[i] = false;
			    internalres[i] = fs;
			    w3 += wt1;
			    if (dd) dwk->dc += wt1 * ltgapf;
			}
			if (*fs == gap_code) {	// -
			    w1 += wt1;		// -- : *-
			    if (dd) {
				if (*bs == gap_code) dwk->dd += wt1;
				else		dwk->cd += wt1;
			    }
			} else {		// *
			    w2 += wt1;
			    if (dd && *bs == gap_code) dwk->dc += wt1;	// -*
			    if (dd && *bs  > gap_code) dwk->cc += wt1;	// **
			}
		    }
		}
		rthk->cfq = w2;
		rthk->dfq = w1 + w0 * ltgapf;
		rthk->efq = w1 + w2 + w0 * ltgapf;
		++rthk;
		if (dd) {
		    dwk[-1].cx = dwk->cc + dwk->cd;
		    dwk[-1].dx = dwk->cd + dwk->dd;
		    ++dwk;
		}
		if (c == 0) break;
	    }
	} else if (internalres) {
	    for (int i = 0; i < many; ++i)
		internalres[i] = ls + i;
	    *rthk = bdy;
	}
	SeqThk*	lthk = rthk;
	DiThk*	ldithk = dwk;

// find left-most right end gap

	if (rtgapf < 1.) {
	    vset(egap, true, many);
	    rthk = thk + thk_len;
	    if (dd) dwk = dithk + thk_len - 1;
	    bs = rs - many;
	    while (rs > ls) {
		VTYPE	w0 = 0;	// ..
		VTYPE	w1 = 0;	// -@	@ = [-*]
		VTYPE	w2 = 0;	// *?	? = [.-*]
		VTYPE	w3 = 0;	// *.
		VTYPE	w4 = 0;	// @@
		int	c = 0;
		for (int i = many; --i >= 0; ) {
		    --rs; --bs;
#if USE_WEIGHT
		    if (vwt1) wt1 = weight[i];
#endif
		    if (*rs == nil_code) {	// ..
			w0 += wt1;
			if (dd) {
			    if (*bs == nil_code) dwk->dd += wt1 * rtgapf;
			    else		 dwk->cd += wt1 * rtgapf;
			}
			++c;
		    } else {
			if (egap[i]) {		// *.
			    w3 += wt1;
			    egap[i] = false;
			} else	w4 += wt1;	// @@
			if (*rs == gap_code) {	// -
			    w1 += wt1;		// -- : *-
			    if (dd) {
				if (*bs == gap_code) dwk->dd += wt1;
				else		dwk->cd += wt1;
			    }
			} else {
			    w2 += wt1;		// -* : **
			    if (dd) {
				if (*bs == gap_code) dwk->dc += wt1;
				else		dwk->cc += wt1;
			    }
			}
		    }
		}
		--rthk;
		rthk->cfq = w2;
		rthk->dfq = w1 + w0 * rtgapf;
		rthk->efq = w4 + (w0 + w3) * rtgapf;
		if (dd) {
		    --dwk;
		    dwk->cx = dwk[1].cc + dwk[1].cd;
		    dwk->dx = dwk[1].dc + dwk[1].dd;
		}
		if (c == 0) break;
	    }
	}
	if (inex.dels) {
	    rthk = lthk;
	    dwk = ldithk;
	    bs = fs - many;
	    while (fs < rs) {
		VTYPE	w1 = 0;
		VTYPE	w2 = 0;
		for (int i = 0; i < many; ++i, ++fs, ++bs) {
#if USE_WEIGHT
		    if (vwt1) wt1 = weight[i];
#endif
		    if (*fs == gap_code) {
			w1 += wt1;		// -- : *-
			if (dd) {
			    if (*bs == gap_code) dwk->dd += wt1;
			    else		 dwk->cd += wt1;
			}
		    } else {
			w2 += wt1;		// -* : **
			if (dd) {
			    if (*bs == gap_code) dwk->dc += wt1;
			    else		 dwk->cc += wt1;
			}
		    }
		}
		rthk->cfq = w2;
		rthk->dfq = w1;
		rthk->efq = sumwt;
		++rthk;
		if (dd) {
		    dwk->cx = w2;
		    dwk->dx = w1;
		    ++dwk;
		}
	    }
	}
	delete[] egap;
}

void initseq(mSeq** seqs, int n)
{
	for ( ; n--; ++seqs) *seqs = new mSeq();
}

void clearseq(mSeq** seqs, int n)
{
	for ( ; n--; ++seqs) delete *seqs;
}

void mSeq::ntor(VTYPE* v, int k, VTYPE w)
{
	int     m;

	switch (k) {
	    case NIL: v[nil_code] += w; break;
	    case _: v[gap_code] += w; break;
	    case A: v[2] += w; break;
	    case C: v[3] += w; break;
	    case G: v[4] += w; break;
	    case T: v[5] += w; break;
	    default:
		m = nbits[k - gap_code];
		w /= m;
		const int*	j = amblist + ambaddr[k];
		for (int n = 0; n < m; ++n) v[*j++] += w;
		break;
	}
}

void mSeq::res2prof(VTYPE v[], CHAR r, VTYPE w)
{
	for (int i = nil_code; ++i < code->max_code; )
	    v[i] += simmtx->mtx[i][r] * w;
}

void mSeq::profile_n(VTYPE v[], VTYPE w[])
{
	v[nil_code] = 0;
	for (int i = nil_code; ++i < felm; ) {
	    int	k = decompact[i];
	    VTYPE*	m = simmtx->mtx[k];
	    v[k] = 0;
	    for (int j = 1; j < felm; ++j)
		v[k] += m[decompact[j]] * w[j];
	}
	for (int i = C; ++i < code->max_code; ) {
	    int	m = nbits[i - gap_code];
	    if (m == 1) continue;
	    v[i] = 0;
	    const int* j = amblist + ambaddr[i];
	    for (int n = 0; n < m; ++n)
		v[i] += v[*j++];
	    v[i] /= m;
	}
}

void mSeq::profile_p(VTYPE v[], VTYPE w[])
{
	v[nil_code] = 0;
	for (int i = nil_code; ++i < felm; ) {
	    VTYPE*	m = simmtx->mtx[i];
	    v[i] = 0;
	    for (int j = nil_code; ++j < felm; )
		v[i] += *++m * w[j];
	}
	v[ASX] = (v[ASN] + v[ASP]) / 2;
	v[GLX] = (v[GLN] + v[GLU]) / 2;
}

void mSeq::profile(VTYPE v[], VTYPE w[])
{
	v[nil_code] = 0;
	for (int i = nil_code; ++i < simmtx->dim; ) {
	    VTYPE*	m = simmtx->mtx[i];
	    v[i] = 0;
	    for (int j = nil_code; ++j < felm; )
		v[i] += *++m * w[j];
	}
}

VTYPE mSeq::unp_prof(VTYPE *w)
{
	VTYPE*	m = simmtx->mtx[gap_code];
	VTYPE*	t = m + felm;
	VTYPE	g = 0;

	while (++m < t) g += *m * *++w;
	return (g);
}

void mSeq::nuc2cvec(VTYPE* v, CHAR* s, FTYPE* w)
{
	VTYPE*	prf = v;
	vclear(prf, felm);
	for (int i = 0; i < many; ++i)
	    ntor(prf, *s++, w? *w++: 1);
}

void mSeq::aas2cvec(VTYPE* v, CHAR* s, FTYPE* w)
{
	VTYPE*	prf = v;
	vclear(prf, felm);
	for (int i = 0; i < many; ++i) {
	    int	k = *s++;
	    VTYPE	wt = w? *w++: 1;
	    switch (k) {
		case ASX:
		    prf[ASN] += wt / 2;
		    prf[ASP] += wt / 2;
		    break;
		case GLX:
		    prf[GLN] += wt / 2;
		    prf[GLU] += wt / 2;
		    break;
		default:
		    prf[k] += wt;
	    }
	}
}

void mSeq::seq2vec(VTYPE* v, CHAR* s, FTYPE* w)
{
	VTYPE*	frq = v;
	vclear(frq, felm);
	for (int i = 0; i++ < many; ) 
	    frq[*s++] += w? *w++: 1;
}

void mSeq::prepseq(bool mksp)
{
#if SSHP
	 if (mksp && inex.sshp && !sshp_prof) makesshpprof(false);
#endif
	switch (inex.molc) {
	    case DNA: case RNA: case GENOME: felm = 6; break;
	    case PROTEIN: felm = ASX; break;
	    case TRON:
	    default: felm = code->max_code; break;
	}
	code->gap_prof = felm + 1;
	nelm = felm + (simmtx? simmtx->dim: 0);
	eth_code = nelm++;
	delete[] pseq;
	pseq = new VTYPE[(len + 2) * nelm];
	vclear(pseq, (len + 2) * nelm);
}

mSeq* mSeq::convseq(INT vect, int step)
{
	if (!thk) mkthick();
	if (inex.dels && !gfq) gfq = new Gfq(this, step);
	if (inex.vect == RAWSEQ) {
#if SSHP
	    if (inex.sshp && !sshp_prof) makesshpprof();
#endif
	    if (vect == RAWSEQ) return (this);

	    void (mSeq::*stov)(VTYPE*, CHAR*, FTYPE*) = &mSeq::seq2vec;
	    switch (inex.molc) {
	      case DNA: case RNA: case GENOME: stov = &mSeq::nuc2cvec; break;
	      case PROTEIN: stov = &mSeq::aas2cvec; break;
	      case TRON: default: felm = code->max_code; break;
	    }
	    prepseq(false);

	    VTYPE*	wst = pseq;
	    CHAR*	soc = at(-1);
	    CHAR*	trm = at(len + 1);
	    FTYPE*	wtb = 0;
#if USE_WEIGHT
	    wtb = new FTYPE[many];
	    for (int i = 0; i < many; ++i) wtb[i] = weight? weight[i]: 1;
#endif	// USE_WEIGHT
	    for ( ; soc < trm; soc += many, wst += nelm) {
		(this->*stov)(wst, soc, wtb);
		if (simmtx) wst[code->gap_prof] = unp_prof(wst);
		VTYPE	eth = 0;
		for (int i = 0; i < many; ++i) {
		    if (soc[i] == gap_code) {
#if USE_WEIGHT
			if (weight) eth += weight[i]; else
#endif	// USE_WEIGHT
			eth += 1;
		    }
		}
		wst[eth_code] = (VTYPE) eth;
	    }
	    delete[] wtb;
	    inex.vect = VECTOR;
	}

	if (inex.vect == vect) return (this);
	if (vect == VECPRO && simmtx) {
	    void (mSeq::*prof)(VTYPE*, VTYPE*) = &mSeq::profile;
	    switch (inex.molc) {
	      case DNA: case RNA: case GENOME:
		prof = &mSeq::profile_n; break;
	      case PROTEIN:
		if (simmtx->dim == simmtx->rows)
		    prof = &mSeq::profile_p;
		else 
		    prof = &mSeq::profile;
		break;
	      default: break;
	    }
	    composition();
	    VTYPE*	wst = pseq;
	    VTYPE*	tst = fat(len + 1);
	    int		nnelm = felm + simmtx->dim;
	    eth_code = nnelm++;
	    if (nnelm == nelm) {
		for ( ; wst < tst; wst += nelm) {
		   (this->*prof)(wst + felm, wst);
		}
	    } else {
		VTYPE*	npseq = new VTYPE[(len + 2) * nnelm];
		VTYPE*	nst = npseq;
		for ( ; wst < tst; wst += nelm, nst += nnelm) {
		    vcopy(nst, wst, felm);
		    (this->*prof)(nst + felm, wst);
		    nst[eth_code] = wst[felm];
		}
		nelm = nnelm;
		delete[] pseq;
		pseq = npseq;
	    }
	    inex.vect = vect;
	}
	byte = inex.prof? nelm * sizeof(VTYPE): many;
	return (this);
}

Seq* mSeq::consenseq(Seq* cs)
{
	if (cs) cs->refresh(1, right - left + 2);
	else cs = new Seq(1, right - left + 2);
	copyattr(cs);
	cs->sname->push(sqname());
	CHAR*	cc = cs->at(0);
	if (!simmtx) simmtx = getSimmtx(0);
	convseq(VECPRO);
	VTYPE*	ss = pat(left);
	VTYPE*	tt = pat(right);
	for ( ; ss < tt; ss += nelm) {
	    VTYPE*	mxv = vmax(ss, code->max_code);
	    int	s = mxv - ss;
	    if (s > gap_code) *cc++ = s;
	    else {
		*cc++ = code->amb_code;
		inex.ambs = 1;
	    }
	}
	cs->postseq(cc);
	return (cs);
}

void antiseq(mSeq** seqs)
{
	mSeq**	as = (*seqs)->anti_;
	if (as) {
	    swapseq(seqs, as);
	    (*seqs)->anti_ = as;
	    (*as)->anti_ = seqs;
	} else (*seqs)->comrev();
}

mSeq* inputseq(mSeq** seqs, char* str)
{
static	const	char aa[] = "aa";
static	const	char nt[] = "nt";
static	const	char tc[] = "tc";
	int 	i = 0;
	const	char*	res;
	mSeq*	sd = 0;

	sscanf(str, "%d", &i);
	if (0 < i && --i < noseq) {
	    sd = seqs[i];
	    while (isdigit(*str)) str++;
	    if (*str == _APPN)	{
		sd->apndseq(0);
		str++;
	    } else
		sd->getseq(0);
	    if (sd->empty()) return (sd);
	    switch (sd->inex.molc) {
		case PROTEIN:	res = aa; break;
		case TRON:	res = tc; break;
		case DNA: case RNA: case GENOME:
		default:	res = nt; break;
	    }
	    fprintf(stderr, " %c%s", sd->senschar(), sd->sqname());
	    if (sd->inex.vect) putc('%', stderr);
	    fprintf(stderr, "[%d] ( %d %s\'s )  [ %d - %d ]\n", 
		sd->many, sd->len, res, sd->SiteNo(sd->left),
		sd->SiteNo(sd->right - 1));
	    sd->clean();
	}
	return (sd);
}

Gep1st::Gep1st(mSeq* sd, int k1, int iv) : many(sd->many)
{
#if USE_WEIGHT
	weight = sd->weight;
	pairwt = sd->pairwt;
#endif	// USE_WEIGHT
	gep1 = new Queue<int>*[many];
	for (int i = 0; i < many; ++i)
	    gep1[i] = new Queue<int>(k1, iv);
}

Gep1st::~Gep1st()
{
	for (int i = 0; i < many; ++i) delete gep1[i];
	delete[] gep1;
}

void Gep1st::shift(const CHAR* s, const int n)
{
	for (int i = 0; i < many; ++i) 
	     if (IsntGap(*s++)) gep1[i]->shift(n);
}

VTYPE Gep1st::longup(const CHAR* s, const int n, const int tgl, bool sft)
{
	VTYPE	lu = 0;
	for (int i = 0; i < many; ++i) {
	    if (IsntGap(*s++)) {
		int	cp = n - (sft? gep1[i]->shift(n): gep1[i]->oldest());
		if (tgl > cp) {
#if USE_WEIGHT
		    if (weight) lu += weight[i]; else
#endif
		    ++lu;
		}
	    }
	}
	return (lu);
}

VTYPE Gep1st::longup(const CHAR* s, const int n, const int* gla)
{
	VTYPE	tlu = 0;
#if USE_WEIGHT
	FTYPE*	pw = pairwt;
	for (int i = 1; i < many; ++i) {
	    if (IsGap(*s++)) {
		if (pw) pw += i;
		continue;
	    }
	    for (int j = 0; j < i; ++j) {
		int	cp = n - gep1[j]->oldest();
		if (gla[i] > cp) {
		    if (pw) tlu += *pw; else
		    ++tlu;
		}
		if (pw) ++pw;
	    }
	}
#else	// USE_WEIGHT
	for (int i = 1; i < many; ++i) {
	    if (IsGap(*s++)) continue;
	    for (int j = 0; j < i; ++j) {
		int	cp = n - gep1[j]->oldest();
		if (gla[i] > cp) ++tlu;
	    }
	}
#endif	// USE_WEIGHT
	return (tlu);
}

VTYPE Gep1st::longup(GFREQ* df, IDELTA* dld, int pos)
{				// half profile
	VTYPE   lunp = 0;
	GFREQ* 	tf = df;

	while (neogfq(tf)) ++tf;
	while (--tf > df) {
	    int	gi = GapLenSD(tf, dld);
	    if (longup(0, pos, gi)) lunp += tf->freq;
	    else break;
	}
	shift(0, pos);
	return (lunp);
}

VTYPE Gep1st::longup(GFREQ* df, IDELTA* dld, SeqItr& si)
{				// both profile
	VTYPE   lunp = 0;
	GFREQ* 	tf = df;

	while (neogfq(tf)) ++tf;
	while (--tf > df) {
	    int	gi = GapLenSD(tf, dld);
	    VTYPE	lu = longup(si.res, si.pos, gi, false) * tf->freq;
	    if (lu == 0) break;
	    lunp += lu;
	}
	shift(si.res, si.pos);
	return (lunp);
}

void mSeqItr::reset(int n, mSeq* sd)
{
	if (sd) {
	    SeqItr::reset(n, (Seq*) sd);
	    msd = sd;
	    nelm = sd->nelm;
	    felm = sd->felm;
	    vss = sd->inex.vect? sd->fat(n): 0;
	    thk_mode = sd->inex.dels? 2: (sd->inex.nils? 1: 0);
	    switch (thk_mode) {
		case 1:
		    sfq = tfq = rfq = 0;
		    repos();
		    break;
		case 2: 
		    dns = sd->thk? sd->thk + n: 0;
		    didns = sd->dithk? sd->dithk + n: 0;
		    if (sd->gfq) {
			sfq = sd->gfq->sfrq + n;
			tfq = sd->gfq->tfrq + n;
			rfq = sd->gfq->rfrq + n;
		    } else	sfq = tfq = rfq = 0;
		    break;
		default: 
		    sfq = tfq = rfq = 0;
		    dns = sd->thk? (sd->inex.sngl? &unit_dns: sd->thk - 1): 0;
		    didns = (sd->dithk && !sd->inex.sngl)? sd->dithk - 1: 0;
		    break;
	    }
	    eth_code = nelm? nelm - 1: 0;
#if SSHP
	    sshp = sd->sat(n);
#endif
	} else if (res) {
	    int	shft = n - pos;
	    SeqItr::reset(n);
	    if (vss) vss += shft * nelm;
	    if (sfq) {sfq += shft; tfq += shft; rfq += shft;}
	    if (thk_mode == 1) repos(); else
	    if (thk_mode == 2) {
		dns += shft;
		if (didns) didns += shft;
	    }
#if SSHP
	    if (sshp) sshp += shft * sshpprm->sshpelems;
#endif
	} else {
	    many = 0; res = 0; bb = 0;
	    nelm = felm = 0; msd = 0;
	    vss = 0; sfq = tfq = rfq = 0; dns = 0; didns = 0;
	    eth_code = 0; pos = n;
#if SSHP
	    sshp = 0;
#endif
	}
}

bool mSeqItr::dullend() {
	if (vss) {
	    if (vss[nil_code] > 0) return true;
	    return (vss[nil_code - nelm] > 0);
	} else {
	    return (res[0] == nil_code || res[-many] == nil_code);
	}
}

void mSeqItr::insertres(mSeqItr& rsi, int bias, VTYPE wt)
{
	memcpy(res + bias, rsi.res, rsi.many);
	if (msd && msd->thk) {
	    SeqThk*	tmp = (SeqThk*) dns;
	    tmp->cfq += wt * rsi.dns->cfq;
	    tmp->dfq += wt * rsi.dns->dfq;
	    tmp->efq += wt * rsi.dns->efq;
	}
	if (msd && msd->dithk) {
	    didns->cc += wt * rsi.didns->cc;
	    didns->cd += wt * rsi.didns->cd;
	    didns->dc += wt * rsi.didns->dc;
	    didns->dd += wt * rsi.didns->dd;
	    didns->cx += wt * rsi.didns->cx;
	    didns->dx += wt * rsi.didns->dx;
	}
#if SSHP
	if (sshp) {
	    for (int i = 0; i < sshpprm->sshpelems; ++i)
		sshp[i] += wt * rsi.sshp[i];
	}
#endif
}

void mSeqItr::insertgap(mSeqItr& rsi, int bias)
{
	CHAR*	ss = res + bias;
	CHAR*	rr = rsi.res;
	CHAR*	tt = rr + rsi.many;
	while (rr < tt)
	    *ss++ = (*rr++ == nil_code)? nil_code: gap_code;
}

#if SSHP

//	secondary structure 

void mSeq::ssprof(mSeqItr& asi)
{
	CHAR*	sq = asi.res;
	VTYPE*&	ssp = asi.sshp;
	VTYPE	ss[NOSS];
#if USE_WEIGHT
	FTYPE*	w = weight;
#else
	FTYPE*	w = 0;
#endif

	for (int i = 0; i < many; ++i, ++sq) {
	    if (IsGap(*sq)) continue;
	    CHAR*	sp = sq;
	    CHAR	rr = is_tron? tron2aa(*sp): *sp;
	    if (IsAA(rr)) {
		rr -= ALA;
		for (int s = 0; s < NOSS; ++s)
		    ss[s] = sshpprm->psstbl[s][SSWING][rr];
	    } else	vclear(ss, NOSS);
	    for (int j = SSWING; (sp += sshp_step) < seq_end && j < SSWIDTH - 1; ) {
		if (IsGap(*sp)) continue;
		if (is_tron && IsTerm(*sp)) break;
		++j;
		CHAR	rr = is_tron? tron2aa(*sp): *sp;
		if (IsAA(rr)) {
		    rr -= ALA;
		    for (int s = 0; s < NOSS; ++s)
			ss[s] += sshpprm->psstbl[s][j][rr];
		}
	    }
	    sp = sq;
	    for (int j = SSWING; (sp -= sshp_step) >= seq_ && j > 0; ) {
		if (IsGap(*sp)) continue;
		if (is_tron && IsTerm(*sp)) break;
		--j;
		CHAR	rr = is_tron? tron2aa(*sp): *sp;
		if (IsAA(rr)) {
		    rr -= ALA;
		    for (int s = 0; s < NOSS; ++s)
			ss[s] += sshpprm->psstbl[s][j][rr];
		}
	    }
	    FTYPE	realwt = w? *w++: evenwt;
	    for (int s = 0; s < NOSS; ++s)
		ssp[s] += realwt * ss[s];
	}
	for (int s = 0; s < NOSS; ++s)
	    ssp[s] = (ssp[s] - sshpprm->sshpav[s]) / sshpprm->sshpsd[s];
}

//	hydrophobicity

void mSeq::hyprof(mSeqItr& asi)
{
	CHAR*	sq = asi.res;
	VTYPE&	hp = asi.sshp[sshpprm->sndstates];
#if USE_WEIGHT
	FTYPE*	w = weight;
#else
	FTYPE*	w = 0;
#endif

	for (int i = 0; i < many; ++i, ++sq) {
	    if (IsGap(*sq)) continue;
	    CHAR*	sp = sq;
	    CHAR	rr = is_tron? tron2aa(*sp): *sp;
	    VTYPE	hh = IsAA(rr)? sshpprm->hphtbl[rr - ALA]: 0;
	    for (int j = 0; (sp += sshp_step) < seq_end && j < alprm3.hpwing; ) {
		if (IsGap(*sp)) continue;
		if (is_tron && IsTerm(*sp)) break;
		++j;
		CHAR	rr = is_tron? tron2aa(*sp): *sp;
		if (IsAA(rr)) hh += sshpprm->hphtbl[rr - ALA];
	    }
	    sp = sq;
	    for (int j = 0; (sp -= sshp_step) >= seq_ && j < alprm3.hpwing; ) {
		if (IsGap(*sp)) continue;
		if (is_tron && IsTerm(*sp)) break;
		++j;
		CHAR	rr = is_tron? tron2aa(*sp): *sp;
		if (IsAA(rr)) hh += sshpprm->hphtbl[rr - ALA];
	    }
	    hp += (w? *w++ * hh: hh);
	}
	if (!w) hp /= many;
	hp /= sshpprm->hpwidth;
	hp = (hp - sshpprm->sshpav[NOSS]) / sshpprm->sshpsd[NOSS];
}

//	hydrophobicity moment

void mSeq::hmprof(mSeqItr& asi, int aid)
{
	CHAR*	sq = asi.res;
	VTYPE&	hm = asi.sshp[sshpprm->sndstates + sshpprm->hphstates + aid];
	VTYPE	hhp[2] = {0, 0};
#if USE_WEIGHT
	FTYPE*	w = weight;
#else
	FTYPE*	w = 0;
#endif

	for (int i = 0; i < many; ++i, ++sq) {
	    if (IsGap(*sq)) continue;
	    CHAR*	sp = sq;
	    CHAR	rr = is_tron? tron2aa(*sp): *sp;
	    VTYPE	hh[2] = {0, IsAA(rr)? sshpprm->hmttbl[rr - ALA]: 0};
	    for (int j = 0; (sp += sshp_step) < seq_end && j < HMWING; ) {
		if (IsGap(*sp)) continue;
		if (is_tron && IsTerm(*sp)) break;
		++j;
		CHAR	rr = is_tron? tron2aa(*sp): *sp;
		if (IsAA(rr)) {
		    rr -= ALA;
		    hh[0] += sshpprm->sincrv[aid][j] * sshpprm->hmttbl[rr];
		    hh[1] += sshpprm->coscrv[aid][j] * sshpprm->hmttbl[rr];
		}
	    }
	    sp = sq;
	    for (int j = 0; (sp -= sshp_step) >= seq_ && j < HMWING; ) {
		if (IsGap(*sp)) continue;
		if (is_tron && IsTerm(*sp)) break;
		++j;
		CHAR	rr = is_tron? tron2aa(*sp): *sp;
		if (IsAA(rr)) {
		    rr -= ALA;
		    hh[0] -= sshpprm->sincrv[aid][j] * sshpprm->hmttbl[rr];
		    hh[1] += sshpprm->coscrv[aid][j] * sshpprm->hmttbl[rr];
		}
	    }
	    FTYPE	wt = w? *w++: evenwt;
	    hhp[0] += wt * hh[0];
	    hhp[1] += wt * hh[1];
	}
	aid += NOSS + 1;
	hm = sqrt(hhp[0] * hhp[0] + hhp[1] * hhp[1]);
	hm = (hm - sshpprm->sshpav[aid]) / sshpprm->sshpsd[aid];
}

void mSeq::sshp_normalize(int i)
{
	VTYPE	sm = 0.;
	VTYPE	sd = 0.;
	VTYPE*	sh = sshp_prof + sshpprm->sshpelems + i;
	for (int n = 0; n++ < sshp_size; sh += sshpprm->sshpelems) {
	    sm += *sh;
	    sd += *sh * *sh;
	}
	VTYPE	av = sm / sshp_size;
	sd = sqrt((sd - av * sm) / (sshp_size - 1));
	sh = sshp_prof + sshpprm->sshpelems + i;
	for (int n = 1; n <= sshp_size; ++n, sh += sshpprm->sshpelems)
	    *sh = (*sh - av) / sd;
}
	
void mSeq::makesshpprof(bool fill)
{
	if (!sshpprm) return;
	inex.sshp = 1;
	int	size = right - left;
	is_tron = istron();
	if (size < SSWIDTH * (is_tron? 3: 1)) return;
	int	bias = left - 1;
	if (sshp_prof && size == sshp_size && bias == sshp_bias) return;
	delete[] sshp_prof;
	sshp_step = many * (is_tron? 3: 1);
	sshp_bias = bias;
	sshp_size = size;
	seq_end = at(len);
	mSeqItr tsi(this, right);
	sshp_prof = new VTYPE[sshpprm->sshpelems * (sshp_size + 2)];
	vclear(sshp_prof, sshpprm->sshpelems * (sshp_size + 2));
	if (!fill) return;
	for (mSeqItr asi(this, left); asi < tsi; ++asi) {
	    if (sshpprm->sndstates) ssprof(asi);
	    if (sshpprm->hphstates) hyprof(asi);
	    for (int a = 0; a < sshpprm->hmtstates; ++a)
		hmprof(asi, a);
	}
}

#endif	// SSHP
