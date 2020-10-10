/*****************************************************************************
*
*	Consistency of two sequence alignments
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
#include "css.h"

inline	int	height(SKL* s) {return (s->m + s->n);}
inline	bool	neoskl(SKL* s) {return (s->m != EOS);}
inline	bool	eofskl(SKL* s) {return (s->m == EOS);}
inline	bool	isdelim(int c) {return (isspace(c) || c == ',' || c == '-');}

static	int	lowerskl(SKL* wsk[]);
static	bool	isedge(SKL* x, SKL* y, SKL* y0);
static	int	commonloop(SKL* lwrb, SKL* here[], SKL* base[], int i, RANGE* wkr[]);

/*	Is *x an edge of y? when    */
/*	height(x) <= height(y)	    */

void fprintrng(FILE* fd, RANGE* rng)
{
	for ( ; neorng(rng); ++rng) {
	    fprintf(fd, "%2d %2d; ", rng->left, rng->right);
	}
	putc('\n', fd);
}


RANGE* getrng(const char* ps)
{
	char	str[MAXL];
	RANGE	rng = {0, 0};
	int	p = 1;
	Mfile	mfd(sizeof(RANGE));
	FILE*	fd = 0;

	if (ps) {
	    if (*ps) {
		while (*ps && isdelim(*ps)) ++ps;
		if (!isdigit(*ps)) fd = fopen(ps, "r");
	    } else	fd = stdin;
	} else		fd = stdin;
		
	mfd.write(&rng);
	for (;;) {
	    if (fd) {
		prompt("~ ");
		ps = fgets(str, MAXL, fd);
	    }
	    if (!ps || !*ps || *ps == '\n') break;
	    while (*ps) {
		while (*ps && isdelim(*ps)) ++ps;
		if (*ps == '\0') break;
		if (p) {
		    rng.left = atoi(ps) - 1;
		    if (rng.left < 0) goto ool;
		} else {
		    rng.right = atoi(ps);
		    if (rng.right < rng.left) goto ool;
		    mfd.write(&rng);
		}
		p = 1 - p;
		while (*ps && !isdelim(*ps)) ++ps;
	    }
	}
ool:
	mfd.write(&endrng);
	p = mfd.size();
	RANGE*	prg  = (RANGE*) mfd.flush();
	if (fd && fd != stdin) fclose(fd);
	prg->left = prg->right = p;
	return (prg);
}

int sumrng(RANGE* rng)
{
	int	sum = 0;

	while (neorng(++rng))
	    sum += rng->right - rng->left;
	return (sum);
}

static int lowerskl(SKL* wsk[])
{
	if (height(wsk[0]) < height(wsk[1])) return (0);
	if (height(wsk[0]) > height(wsk[1])) return (1);
	return (wsk[0]->m <= wsk[1]->m)? 0: 1;
}

static bool isedge(SKL* x, SKL* y, SKL* y0)
{
	while (y > y0 && height(y) > height(x)) --y;
	SKL*	z = y + 1;
	if (eofskl(z)) return false;
	return ((x->m - y->m) * (z->n - y->n)
	     == (x->n - y->n) * (z->m - y->m));
}

static int commonloop(SKL* lwrb, SKL* here[], SKL* base[], int i, RANGE* wkr[])
{
	SKL*	wsk[2];
	int	j = 1 -i;
	int	c = 0;
	int	hlwrb = height(lwrb);
	RANGE	grg[2][2];

	for (wsk[0] = here[0]; wsk[0] > base[0] && height(wsk[0]) > hlwrb; --wsk[0]) ;
	if (height(wsk[0]) < hlwrb) ++wsk[0];
	for (wsk[1] = here[1]; wsk[1] > base[1] && height(wsk[1]) > hlwrb; --wsk[1]) ;
	if (height(wsk[1]) < hlwrb) ++wsk[1];
	++wsk[0]; ++wsk[1];
	grg[0][0] = grg[0][1] = grg[1][0] = grg[1][1] = zerorng;
	hlwrb = height(here[i]);
	for (;;) {
	    i = lowerskl(wsk);
	    j = 1 - i;
	    if (eofskl(wsk[i]) || height(wsk[i]) > hlwrb) break;
	    int	dm = wsk[i]->m - wsk[i][-1].m;
	    int	dn = wsk[i]->n - wsk[i][-1].n;
	    if (dm > 0 && dn == 0) {
		grg[i][0].left = wsk[i][-1].m;
		grg[i][0].right = wsk[i]->m;
		wkr[0]->right = min(grg[i][0].right, grg[j][0].right);
		wkr[0]->left = max(grg[i][0].left, grg[j][0].left);
		dm = wkr[0]->right - wkr[0]->left;
		if (dm > 0 && wkr[0][-1].right < wkr[0]->left) {
		    c += dm;
		    wkr[0]++;
		}
	    } else if (dm == 0 && dn > 0) {
		grg[i][1].left = wsk[i][-1].n;
		grg[i][1].right = wsk[i]->n;
		wkr[1]->right = min(grg[i][1].right, grg[j][1].right);
		wkr[1]->left = max(grg[i][1].left, grg[j][1].left);
		dn = wkr[1]->right - wkr[1]->left;
		if (dn > 0 && wkr[1][-1].right < wkr[1]->left) {
		    c += dn;
		    wkr[1]++;
		}
	    }
	    ++wsk[i];
	}
	return (c);
}

int css(RANGE* rngs[2], SKL* skl[2])
{
	bool	f = false;
	int 	c = 0;
	RANGE*	wkr[2];
	int	nsk[2];
	SKL*	wsk[2];

	nsk[0] = skl[0]->n;
	nsk[1] = skl[1]->n;
	SKL*	prv = ++skl[0];
	wsk[0] = prv;
	wsk[1] = ++skl[1];
	int	i = nsk[0] + nsk[1];
	rngs[0] = new RANGE[i];
	rngs[1] = new RANGE[i];
	wkr[0] = rngs[0];
	wkr[1] = rngs[1];
	(wkr[0]++)->right = wsk[0]->m;
	(wkr[1]++)->right = wsk[0]->n;

	while (neoskl(wsk[0]) && neoskl(wsk[1])) {
	    i = lowerskl(wsk);
	    int	j = 1 - i;
	    if (isedge(wsk[i], wsk[j], skl[j])) {
		if (!f) {
		    c += commonloop(prv, wsk, skl, i, wkr);
		    wkr[0]->left = wsk[i]->m;
		    wkr[1]->left = wsk[i]->n;
		    f = true;
		}
		prv = wsk[i];
	    } else if (f) {
		c += prv->m - wkr[0]->left + prv->n - wkr[1]->left;
		(wkr[0]++)->right = prv->m;
		(wkr[1]++)->right = prv->n;
		f = false;
	    }
	    wsk[i]++;
	}
	if (f) {
	    c += prv->m - wkr[0]->left + prv->n - wkr[1]->left;
	    (wkr[0]++)->right = prv->m;
	    (wkr[1]++)->right = prv->n;
	}
	*wkr[0] = *wkr[1] = endrng;
	rngs[0]->left = rngs[0]->right = ++wkr[0] - rngs[0];
	rngs[1]->left = rngs[1]->right = ++wkr[1] - rngs[1];
	--skl[0]; --skl[1];
	return (c);
}

void unfoldrng(RANGE* rng, GAPS* gg)
{
	int	len = 0;

	GAPS*	g = gg + 1;
	while (neorng(++rng)) {
	    while (gaps_intr(g) && g->gps < rng->left)
		len += (g++)->gln;
	    rng->left += len;
	    while (gaps_intr(g) && g->gps <= rng->right)
		len += (g++)->gln;
	    rng->right += len;
	}
}

void foldrng(RANGE* rng, GAPS* gg)
{
	int	len = 0;
	RANGE*	org = rng;
	RANGE*	fd = rng + 1;

	GAPS*	g = gg + 1;
	while (neorng(++org)) {
	    while (gaps_intr(g) && g->gps + g->gln < org->left)
		len += (g++)->gln;
	    fd->left = min(org->left, g->gps) - len;
	    while (gaps_intr(g) && g->gps + g->gln <= org->right)
		len += (g++)->gln;
	    (fd++)->right = min(org->right, g->gps) - len;
	}
	*fd = endrng;
	rng->left = rng->right = ++fd - rng;
}

RANGE*	copyrng(RANGE* s)
{
	RANGE*	d = new RANGE[sizerng(s)];
	RANGE*	w = d;

	do {
	    *w++ = *s;
	} while (neorng(s++));
	return (d);
}

RANGE*	cmnrng(RANGE* a, RANGE* b, int leave0)
{
	RANGE*	rc = new RANGE[sizerng(a++) + sizerng(b++) - NSENT];
	RANGE*	rw = rc + 1;
	RANGE** ll;
	RANGE** uu;

	while (neorng(a) && neorng(b)) {
	    if (a->right <= b->right)
		 {ll = &a; uu = &b;}
	    else {ll = &b; uu = &a;}
	    if ((*ll)->right >= (*uu)->left) {
		rw->left = max(a->left, b->left);
		rw->right = (*ll)->right;
		if (rw->right - rw->left > 0 || leave0) rw++;
	    }
	    (*ll)++;
	}
	*rw = endrng;
	rc->left = rc->right = ++rw - rc;
	return (rc);
}

RANGE*	uniterng(RANGE* a, RANGE* b, int leave0)
{
	RANGE*	rc = new RANGE[sizerng(a++) + sizerng(b++) - NSENT];
	RANGE*	rw = rc + 1;
	RANGE** ll;
	RANGE** uu;
	int	link = 0;

	while (neorng(a) || neorng(b)) {
	    if (!neorng(b) || (neorng(a) && a->right <= b->right))
		 {ll = &a; uu = &b;}
	    else {ll = &b; uu = &a;}
	    if (!link)
		rw->left = (neorng(*uu) && (*uu)->left < (*ll)->left)?
		    (*uu)->left: (*ll)->left;
	    if (neorng(*uu) && (*ll)->right >= (*uu)->left)
		link = 1;
	    else {
		link = 0;
		rw->right = (*ll)->right;
		if (rw->right - rw->left > 0 || leave0) rw++;
	    }
	    (*ll)++;
	}
	*rw = endrng;
	rc->left = rc->right = ++rw - rc;
	return (rc);
}

RANGE*	complerng(RANGE* c, RANGE* a)
{
	if (emptyrng(a)) return (copyrng(c));
	int	n = sizerng(a) + 1;
	RANGE*	b = lastrng(a);
	RANGE*	bb;

	if ((++c)->right == b->right) --n;
	if (c->left == (++a)->left) --n;
	bb = b = new RANGE[n];
	(++b)->left = c->left;
	if (n > NSENT) {	// not empty
	    for ( ; neorng(a); ++a) {
		b->right = a->left;
		if (b->left < b->right) ++b;
		b->left = a->right;
	    }
	    b->right = c->right;
	    if (b->left < b->right) ++b;
	}
	*b = endrng;
	bb->left = bb->right = ++b - bb;
	return (bb);
}

RANGE*	singlerng(RANGE* rng, int left, int right)
{
	if (!rng) rng = new RANGE[3];
	rng[0].left = rng[0].right = 3;
	rng[1].left = left;
	rng[1].right = right;
	rng[2] = endrng;
	return (rng);
}

#if MAIN
#include "css.inc"
#endif
