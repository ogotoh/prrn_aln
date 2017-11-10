/*****************************************************************************
*
*	pairs.cc
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

#include "seq.h"
#include "mseq.h"
#include "phyl.h"
#include "pairs.h"
#include <math.h>

static	void	fminidx(int& ii, int& jj, const FTYPE* array, const int n);

static void fminidx(int& ii, int& jj, const FTYPE* array, const int n)
{
	FTYPE	fmin = *array++;

	ii = 0;
	jj = 1;
	for (int j = 2; j < n; j++) 
	    for (int i = 0; i < j; i++, array++)
		if (*array < fmin) {
		    fmin = *array;
		    ii = i;
		    jj = j;
		}
}

void Rnode::extreme(int*& leaf)
{
	if (id >= 0) {	// leaf
	    *leaf = id;
	    rm = lm = leaf++;
	} else {	// internal node
	    left->extreme(leaf);
	    lm = left->lm;
	    right->extreme(leaf);
	    rm = right->rm;
	}
}

Bitree::Bitree(FTYPE* dist, int members)
{
	pnode = new Rnode*[members];
	nodes = new Rnode[2 * members];
	leaves = new int[members];

	root = nodes;
	for (int m = 0; m < members; ++m, ++root) {
	    root->ndesc = 1;
	    root->id = m;
	    root->left = root->right = 0;
	    root->height = 0.;
	    root->lm = root->rm = 0;
	    pnode[m] = root;
	}

	int	inode = members;
	for (int m = members; m > 1; ) {
	    int	i, j;
	    fminidx(i, j, dist, m--);
	    root->left = pnode[i];
	    root->right = pnode[j];
	    root->id = -inode++;
	    root->ndesc = pnode[i]->ndesc + pnode[j]->ndesc;
	    root->height = dist[elem(i, j)] / 2;
	    for (int k = 0; k <= m; k++) {
		if (k == i || k == j) continue;
		FTYPE	dd = pnode[i]->ndesc * dist[elem(k, i)]
			   + pnode[j]->ndesc * dist[elem(k, j)];
		dist[elem(k, i)] = dd / root->ndesc;
	    }
	    for (int k = 0; k < m; k++)
		if (k != j) dist[elem(k, j)] = dist[elem(k, m)];
	    pnode[i] = root++;
	    pnode[j] = pnode[m];
	}
	int*	leaf = leaves;
	root->extreme(leaf);
}

Bitree::~Bitree()
{
	delete[] pnode;
	delete[] nodes;
	delete[] leaves;
}

bool PairPrm::pair_prm()
{
	double	uc = upcutoff;
	double	lc = lwcutoff;

	promptin("Cutoff (%6.3lf > %6.2lf) : ", &upcutoff, &lwcutoff);
	Among	c = species;
	prompt("Species: ");
	promptin("Intra[1]/Inter[2]/All[3]/Give[4]/Tree[5] (%d) : ", &c);
	promptin("Show pairs [0-2] (%d) : ", &showmode);
	bool	update = (c != species || uc != upcutoff || lc != lwcutoff);
	if (1 <= c && c <= 5) species = c;
	return (update);
}

void Pairs::givenpair(Seq* sd)
{
	int	n = sd->many*(sd->many-1)/2;
	PAIR*	wk = pairs;

	for (npair = 0; npair < n; npair++, wk++) {
	    promptin("Seq #\'s (%d, %d): ", &(wk->x), &(wk->y));
	    if (wk->x < 0 || wk->y < 0 || wk->x == wk->y)
		break;
	    wk->w = 1.;
	}
}

void Pairs::nodepair(Rnode* node, const PairPrm* pp)
{
	double	hi = node->height * 2;

	if (!node->left || !node->right) return;
	nodepair(node->left, pp);
	nodepair(node->right, pp);
	if (pp->lwcutoff > hi || hi >= pp->upcutoff) return;
	double	wt = 1. / ((node->left->rm - node->left->lm + 1) *
		(node->right->rm - node->right->lm + 1));
	for (int* pi = node->left->lm; pi <= node->left->rm; pi++)
	    for (int* pj = node->right->lm; pj <= node->right->rm; pj++) {
		pairs[npair].x = *pi;
		pairs[npair].y = *pj;
		pairs[npair++].w = wt;
	    }
}

void Pairs::treepair(Seq* sd, const PairPrm* pp)
{
	FTYPE*	dist = calcdist(sd, 0);
	Bitree	bitree(dist, sd->many);
	npair = 0;
	nodepair(bitree.root, pp);
	delete[] dist;
}

void Pairs::combpair(Seq* sd, const PairPrm* pp)
{

	FTYPE*	dstmat = calcdist(sd, 0);
	FTYPE*	wk = dstmat;
	for (int j = 1; j < sd->many; j++) {
	    for (int i = 0; i < j; i++, wk++) {
		bool	condition = pp->lwcutoff <= *wk && *wk < pp->upcutoff;
		if (pp->species != ALL) {
		    bool	samespec = !strncmp((*sd->sname)[i], (*sd->sname)[j], pp->Nspc);
		    bool	intraspc = pp->species == INTRA;
		    condition = condition && 
			((intraspc && samespec) || (!intraspc && !samespec));
		}
		if (condition) {
		    pairs[npair].x = i;
		    pairs[npair].y = j;
		    pairs[npair++].w = 1.;
		}
	    }
	}
	delete[] dstmat;
}

void Pairs::intergroup(Seq* sd)
{
	Subset* ss = new Subset(sd->many, stdin);
	if (!ss) return;
	for (int j = 1; j < ss->num; ++j) {
	    for (int i = 0; i < j; ++i) {
		for (int* gi = ss->group[i]; *gi >= 0; ++gi) {
		    for (int* gj = ss->group[j]; *gj >= 0; ++gj) {
			pairs[npair].x = *gi;
			pairs[npair].y = *gj;
			pairs[npair++].w = 1.;
		    }
		}
	    }
	}
	delete ss;
}

void Pairs::intragroup(Seq* sd)
{
	Subset* ss = new Subset(sd->many, stdin);
	if (!ss) return;
	for (int i = 0; i < ss->num; ++i) {
	    for (int* gj = ss->group[i]; *gj >= 0; ++gj) {
		for (int* gi = ss->group[i]; gi < gj; ++gi) {
		    pairs[npair].x = *gi;
		    pairs[npair].y = *gj;
		    pairs[npair++].w = 1.;
		}
	    }
	}
	delete ss;
}

void Pairs::showpairs(Seq* sd)
{
	PAIR*	pr = pairs;
	for (int i = 0; i < npair; i++, pr++) {
	    printf("%4d %4d %4d %8.4f   %-12s %-12s\n", 
		elem(pr->x, pr->y), pr->x, pr->y, pr->w, 
		(*sd->sname)[pr->x], (*sd->sname)[pr->y]);
	}
}

Pairs::Pairs(Seq* sd, const PairPrm* pp)
{
	pairs = new PAIR[sd->many*(sd->many-1)/2];
	npair = 0;
	switch (pp->species) {
	    case INTRA:		intragroup(sd); break;
	    case INTER:		intergroup(sd); break;
	    case GIVE:		givenpair(sd); break;
	    case BITREE:	treepair(sd, pp); break;
	    default:		combpair(sd, pp); break;
	}
	switch (pp->showmode) {
	    case 1: fprintf(stderr, "# of pairs = %d\n", npair); break;
	    case 2: showpairs(sd); break;
	    default: break;
	}
}
