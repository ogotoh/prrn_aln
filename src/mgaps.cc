/*****************************************************************************
*
*	abstract structures of alingmnet
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

#include "aln.h"
#include "mseq.h"
#include "gfreq.h"
#include "mgaps.h"

#define DEBUG	0
#define CDEBUG	1

static const	char memerror[] = "Memory error at %s, %d > %d";
static const	char gapsidf[] = "Gaps structure: %d\n";

static	void	copyglist(GAPS** dest, GAPS** sorc, int n, int span = 0);

/*	Each GAPS* variable has a header containing the length of the
	alignment and the number of elements in that array.  These values
	are returned by the functions defined in "gaps.h".
*/

inline GAPS* copygaps(GAPS* d, GAPS* ss)
{
	vcopy(d, ss, ss->gln);
	return (d);
}

static void copyglist(GAPS** dest, GAPS** sorc, int n, int span)
{
	for (int i = 0; i < n; ++i) {
	    GAPS*&	g = *sorc++;
	    GAPS*&	d = *dest++;
	    if (!d) d = new GAPS[span? span: gaps_size(g)];
	    vcopy(d, g, g->gln);
	}
}

GapsList::GapsList(int n, int spn, bool wk)
	: size(n), span(spn), work(wk), num(0)
{
	glst = new GAPS*[n];
	vclear(glst, n);
	if (spn) {
	    GAPS**	gg = glst;
	    for (int i = 0; i < n; ++i) {
		GAPS*&	g = *gg++;
		g = new GAPS[spn];
		g->gps = g->gln = 0;
	    }
	}
	if (work && span) {
	    tmp = new GAPS[span];
	    res = new GAPS[span + span];
	} else {
	    tmp = res = 0;
	}
}
		
GapsList::GapsList(GapsList& src, bool wk) 
	: size(src.size), span(src.span), work(wk), num(src.num)
{
	glst = new GAPS*[size];
	if (work && span) {
	    tmp = new GAPS[span];
	    res = new GAPS[span + span];
	} else {
	    tmp = res = 0;
	}
	vclear(glst, size);
	copyglist(glst, src.glst, size, span);
}

GapsList::GapsList(FILE* fd)
	: span(0), work(false), tmp(0), res(0), glst(0)
{
	if (fscanf(fd, gapsidf, &num) < 1) goto input_error;
	size = num;
	glst = new GAPS*[size];
	for (int i = 0; i < num; ++i) {
	    GAPS	gap;
	    if (fscanf(fd, "%d %d", &gap.gps, &gap.gln) < 2)
		goto input_error;
	    GAPS*	gp = glst[i] = new GAPS[gap.gln];
	    *gp++ = gap;
	    for (int j = 0; ++j < gap.gln; ++gp) {
		if (fscanf(fd, "%d %d", &gp->gps, &gp->gln) < 2)
		    goto input_error;
	    }
	}
	return;
input_error:
	prompt("Not or incomplete gaps struture file!\n");
	delete[] glst; glst = 0; num = size = 0;
}

GapsList::~GapsList()
{
	for (int i = 0; i < size; ++i)
	    delete[] glst[i];
	delete[] glst; delete[] tmp; delete[] res;
}

void GapsList::copy_to(GapsList* dst)
{
	dst->num = num;
	copyglist(dst->glst, glst, num, span);
}

void GapsList::copy_to(GapsList* dst, int *gr)
{
	int	ngr = 0;
	while (gr[ngr] >= 0) ++ngr;
	dst->num = ngr;
	GAPS**	dest = dst->glst;
	for (int i = 0; i < ngr; ++i) {
	    GAPS*	g = glst[*gr++];
	    vcopy(*dest++, g, g->gln);
	}
}

void GapsList::print(FILE* fd)
{
	fprintf(fd, gapsidf, num);
	for (int i = 0; i < num; ++i) {
	    GAPS*	gp = glst[i];
	    do {
		fprintf(fd, " %d %d", gp->gps, gp->gln);
	    } while (gaps_intr(gp++));
	    putc('\n', fd);
	}
}

void GapsList::merge_from(GapsList* src, int *gr)
{
	GAPS**	sorc = src->glst;
	while (*gr >= 0) {
	    GAPS*	g = *sorc++;
	    vcopy(glst[*gr++], g, g->gln);
	}
}

GapsList* GapsList::resize(int spn)
{
	if (spn < span) return (this);
	span = spn;
	GAPS**	gg = glst;
	for (int i = 0; i < num; ++i) {
	    GAPS*&	g = *gg++;
	    if (g->gln > spn)
		fatal("Panic: Bad gaps list %d > %d !\n", g->gln, spn);
	    GAPS*	d = new GAPS[span];
	    vcopy(d, g, g->gln);
	    delete[] g;
	    g = d;
	}
	return (this);
}
	    
/****************************************
	removes blank columns
****************************************/

GAPS* GapsList::delcommongap()
{
#if DEBUG
	print(stdout);
#endif
	if (!tmp || !res)
	    fatal("set up work space for gaps list !\n");

	copygaps(res, *glst);
	for (int i = 1; i < num; ++i) {
	    copygaps(tmp, res);
	    GAPS*	g = glst[i] + 2;
	    GAPS*	h = res + 2;
	    GAPS*	f = tmp + 2;
	    while (gaps_intr(f) && gaps_intr(g)) {
		GAPS**	uu;
		GAPS**	ll;
		if (f->gps + f->gln <= g->gps + g->gln)
		    {ll = &f; uu = &g;}
		else
		    {ll = &g; uu = &f;}
		int	upp = (*ll)->gps + (*ll)->gln;
		if ((*uu)->gps < upp) {
		    h->gps = max((*ll)->gps, (*uu)->gps);
		    h->gln = upp - h->gps;
		    h++;
		}
		(*ll)++;
	    }
	    while (gaps_intr(f) && f->gps >= g->gps) *h++ = *f++;
	    while (gaps_intr(f)) ++f;
	    while (gaps_intr(g) && g->gps >= f->gps) *h++ = *g++;
	    while (gaps_intr(g)) ++g;
	    *h = *g;
	    res->gps = h->gps;
	    res->gln = ++h - res;
	}
	if (gaps_free(res)) goto funcend;
	if (gaps_size(res) > span) 
	    fatal(memerror, "delcommongap()", gaps_size(res), span);

	for (int i = 0; i < num; ++i) {
	    GAPS*	gg = glst[i];
	    GAPS*	g = gg + 1;
	    GAPS*	h = res + 2;
	    GAPS*	f = g + 1;
	    int		ln = 0;
	    while (gaps_intr(++g)) {
		int	upp = g->gps + g->gln;
		g->gps -= ln;
		while (gaps_intr( h) && h->gps + h->gln <= upp) {
		    g->gln -= h->gln;
		    ln += h->gln;
		    h++;
		}
		if (g->gln) *f++ = *g;
	    }
	    g->gps -= ln;
	    *f++ = *g;
	    gg->gps = g->gps - gg[1].gps;
	    gg->gln = f - gg;
	}
funcend:
	return (res);
}

void GapsList::insertgap(GAPS* hh)
{
	if (gaps_free(hh)) return;	//  no gap to be inserted
	for (int i = 0; i < num; ++i) {
	    int		len = 0;
	    GAPS*	ff = glst[i];
	    copygaps(tmp, ff);
	    GAPS*	g = tmp;
	    GAPS*	h = hh + 1;
	    GAPS*	f = ff + 1;
	    *f++ = *++g;
	    do {
		while (gaps_intr(h) && h->gps < g->gps) {
		    f->gps = h->gps + len;
		    (f++)->gln = h->gln;
		    len += (h++)->gln;
		}
		f->gps = g->gps + len;
		for (f->gln = 0; gaps_intr(g); ++g) {
		    f->gln += g->gln;
		    if (g[1].gps != g->gps) break;
		}
		int	upp = g->gps + f->gln;
		while (gaps_intr( h) && h->gps <= upp) {
		    f->gln += h->gln;
		    len += (h++)->gln;
		}
		if (f->gln) f++;
	    } while (gaps_intr(g++)); 
	    *f = *--g;
	    ff->gps = (f->gps += len) - (ff + 1)->gps;
	    ff->gln = ++f - ff;
	}
}

mSeq* aggregate(mSeq* sd, mSeq** slist, GAPS** glst, FTYPE* wtlst)
{
	int	many = 0;
	int	len = glst? gaps_span(*glst): (*slist)->len;

	for (mSeq** sq = slist; *sq; ++sq) many += (*sq)->many;
	if (!many) return (0);
	if (sd) sd->refresh(many, len);
	else	sd = new mSeq(many, len);
	sd->code = (*slist)->code;
	(*slist)->copyattr(sd);
	sd->inex.dels = sd->inex.ambs = 
	sd->inex.vect = sd->inex.gfrq = sd->inex.prof = 0;
#if USE_WEIGHT
	FTYPE*	wkwt = 0;
	int	k = 0;
	if (wtlst) {
	    if (!sd->weight) sd->weight = new FTYPE[many];
	    wkwt = wtlst;
	}
#endif
	mSeqItr	ssi(sd);
	mSeq** sq = slist;
	for (int bias = 0; *sq; bias += (*sq++)->many) {
	    if ((*sq)->inex.ambs) sd->inex.ambs = 1;
	    if ((*sq)->inex.dels) sd->inex.dels = 1;
	    FTYPE	wt = 1;
#if USE_WEIGHT
	    if (wtlst) {
		wt = *wkwt++;
		for (int j = 0; j < (*sq)->many; ++j)
		    if ((*sq)->weight)
			sd->weight[k++] = (*sq)->weight[j] * wt;
		    else
			sd->weight[k++] = wt;
	    }
#endif
	    ssi.reset(-1);
	    mSeqItr	rsi(*sq, -1);
	    if (glst) {
		GAPS*	gp = *glst;
		if (!gaps_free(gp++)) sd->inex.dels = 1;
		for (int i = -1; i++ < len; ++ssi) {
		    if (gaps_intr(gp) && (gp->gps + gp->gln) < i) ++gp;
		    if (gp->gps < i) ssi.insertgap(rsi, bias);
		    else {
			ssi.insertres(rsi, bias, wt);
			++rsi;
		    }
		}
		++glst;
	    } else 
		for (int i = -1; i++ < len; ++ssi, ++rsi)
		    ssi.insertres(rsi, bias, wt);
	}
	sd->postseq(ssi.res);
	sd->inex.nils = slist[0]->inex.nils;
	sd->mkthick();
	return (sd);
}

/*******************************************
       |             ||   |  intron (unfold -> fold)
_______V_____________VV___V__________________
           ||       |  |     gap (unfold)

*********************************************/

bool synthgap(GapsList* glists[], SKL* skl, int* lst[])
{
	GAPS*	gp[2];

	skl2gaps(gp, skl, true);
	int	len[2];
	for (int i = 0; i < 2; ++i) {
	    glists[i+1]->insertgap(gp[i]);
	    len[i] = glists[i+1]->glst[0][0].gps;
	}
	if (len[0] == len[1]) {
	    for (int i = 0; i < 2; ++i) {
		glists[0]->merge_from(glists[i+1], lst[i]);
	    }
	} else {
	    fprintf(stderr, "Synthgap Failed: %d != %d !\n", len[0], len[1]);
	}
	delete[] gp[0]; delete[] gp[1];
	return (len[0] == len[1]);
}

mSeq* aggregate_sb(mSeq* sd, mSeq** slist, GapsList* glist, FTYPE* wtlst)
{
	GAPS**	gsrc = glist? glist->glst: 0;
	sd = aggregate(sd, slist, gsrc, wtlst);
	if (alprm2.spb <= 0 || !gsrc) return (sd);
	sd->sigII = new SigII((const Seq**) slist, (const GAPS**) gsrc, wtlst);
	if (sd->sigII->pfqnum) {
	    sd->sigII->resetend(sd->len);
	} else {
	    delete sd->sigII;
	    sd->sigII = 0;
	}
	return (sd);
}

/* return a subset of glst */ 
GAPS* gatherseq(mSeq* sd, mSeq** slst, GapsList* glst, 
	GapsList* gbuf, FTYPE* wlst, int* group)
{
	GAPS*	hh = 0;

	if (group) {
	    GapsList*	grsv = gbuf;
	    if (!gbuf) gbuf = new GapsList(*glst, true);
	    glst->copy_to(gbuf, group);
	    hh = gbuf->delcommongap();
	    if (gbuf->num == 1) {
		mSeq*&	src = slst[*group];
		if (src->many == 1) src->aliaseq(sd);
		else {
		    src->copyseq(sd);
#if USE_WEIGHT
		    sd->clean();
		    for (int j = 0; j < sd->many; ++j)
			sd->weight[j] *= wlst[*group];
		    sd->mkthick();
#endif
		}
	    } else {
		mSeq**	sbuf = new mSeq*[gbuf->num + 1];
		mSeq**	sq = sbuf;
		FTYPE*	wbuf =  wlst? new FTYPE[gbuf->num + 1]: 0;
		FTYPE*	ww = wbuf;
		for (int* gr = group; *gr >= 0; ++gr) {
		    *sq++ = slst[*gr];
		    if (wlst) *ww++ = wlst[*gr];
		}
		*sq = 0;
		aggregate_sb(sd, sbuf, gbuf, wbuf);
		delete[] sbuf;
		delete[] wbuf;
	    }
	    if (!grsv) delete gbuf;
	} else {
	    aggregate_sb(sd, slst, glst, wlst);
	}
	sd->exg_seq(sd->inex.exgl, sd->inex.exgr);
	return (hh);
}

void incrgap(int* gg, CHAR* ss, int n)
{
	if (ss) {
	    for ( ; n--; ++ss) {
		if (IsGap(*ss)) ++(*gg++);
		else    *gg++ = 0;
	    }
	} else
	    while (n--) ++(*gg++);
}

void elongap(int* gl, int* pr, CHAR* s, int n)
{
	if (s) {
	    for ( ; n--; pr++, ++s)
		*gl++ = IsGap(*s)? *pr + 1: 0;
	} else {
	    while (n--)
		*gl++ = *pr++ + 1;
	}
}
