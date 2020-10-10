/*****************************************************************************
*
*	Generate a series of psuedo random numbers each of which codes
*	for a binary expression of a two-way partition of a set
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
class	mSeq;
#include "mgaps.h"
#include "phyl.h"
#include "consreg.h"
#include "randiv.h"

#define	 EOL	(-1)

//	 Mixed congruence method

McRand::McRand(int p, int rn)
{
	if (rn == 0) {			// non-random
	    mrand = false;
	    mcoef = 1;
	} else {			// seed
	    mrand = true;
	    mcval = (MCTYPE) (rn == 1? rand(): rn);
	    srand(mcval);
	}
	if (p >= WordBits) {
	    prompt("MCRAND: Too large p value %d!\n", p);
	    p %= WordBits;
	}
	mcmod = 1L << p;
	if (mrand) {
	    mcoef = (rand() / 4 * 4 + 5) % mcmod;
	    mcval %= mcmod;
	} else	mcval = mcmod - 1;
}

int setrandiv(int dm)
{
	static	int	divmode = TREEDIV;
	static	int	divmodeB = TREEDIV;

	setintval("Div mode: No[0]/Single[1]/Tree[2]/Every[3]/Part[4] (%d) : ",
		&divmode, &divmodeB, dm);
	return (divmode);
}

void Randiv::bin2lst(int* wk, LRAND lnbr)
{
	MCTYPE	m = 0;

	for (int i = 0; i < rand_bit; i++, m >>= 1) {
	    if (i % WordBits == 0) m = *lnbr++;
	    if (m & 1) *wk++ = i;
	}
	*wk = EOL;
}

void Randiv::bin2lst2(int* lst[2], LRAND lnbr)
{
	MCTYPE	m = 0;
	int*	wk[2] = {lst[0], lst[1]};

	for (int i = 0; i < rand_bit; m >>= 1) {
	    if (i % WordBits == 0) m = *lnbr++;
	    *wk[m & 1]++ = i++;
	}
	*wk[0] = *wk[1] = EOL;
}

void Randiv::bin2lst2g(int* lst[2], LRAND lnbr, Subset* ss)
{
	MCTYPE	m = 0;
	int*	wk[2] = {lst[0], lst[1]};

	for (int i = 0; i < rand_bit; i++, m >>= 1) {
	    if (i % WordBits == 0) m = *lnbr++;
	    int	k = m & 1;
	    for (int* g = ss->group[i]; *g >= 0; )
		*wk[k]++ = *g++;
	}
	*wk[0] = *wk[1] = EOL;
}

void putlst(FILE* fd, int* lst)
{
	while (*lst != EOL) fprintf(fd, "%2d ", *lst++);
}

int Randiv::count1(LRAND lnbr)
{
	int	k = 0;
	MCTYPE	m = 0;

	if (!lnbr) lnbr = partab;
	for (int i = 0; i < rand_bit; i++, m >>= 1) {
	    if (i % WordBits == 0) m = *lnbr++;
	    if (m & 1) k++;
	}
	return (min(k, rand_bit - k));
}

MCTYPE* Randiv::fill_tree_tab(Knode* node)
{
	MCTYPE* dst = partab + node->tid * rand_dim;

	if (node->isleaf()) {
	    dst[node->tid / WordBits] = 1L << (node->tid % WordBits);
	} else {
	    MCTYPE* a = fill_tree_tab(node->left);
	    MCTYPE* b = fill_tree_tab(node->right);
	    for (int i = 0; i < rand_dim; i++)
		dst[i] = a[i] | b[i];
	}
	return (dst);
}

Randiv::~Randiv()
{
	delete[] partab;
	delete	dtree;
	delete	mcr;
}

void Randiv::one_tab(int nn)
{
	MCTYPE*	dst = partab = new MCTYPE[nn * rand_dim];
	for (int i = 0; i < nn * rand_dim; ++i) partab[i] = 0;
	for (int i = 0; i < nn; ++i, dst += rand_dim) {
	    dst[i / WordBits] = 1L << (i % WordBits);
	}
}

void Randiv::tree_tab(Knode* root, int nn)
{
	partab = new MCTYPE[2 * nn * rand_dim];
	for (int i = 0; i < 2 * nn * rand_dim; ++i) partab[i] = 0;
	fill_tree_tab(root);
}

Randiv::Randiv(Seq* sd, Ssrel* srl, DivMode dm, int rn) : divmode(dm)
{
	int	nn = srl->ss? srl->ss->num: sd->many;

	rand_bit = nn;
	rand_dim = (nn + WordBits - 1) / WordBits;
	partab = 0; mcr = 0; dtree = 0;

	if (divmode == ONE_DIV) {
	    cycle = nn;
	    one_tab(nn);
	} else {
	    cycle = 2 * nn - 3;
	    tree_tab(srl->ktree->root, nn);
	}
	nn = 0;
	for (MCTYPE x = 1L; x < cycle; nn++) x <<= 1;
	mcr = new McRand(nn, rn);
}

Randiv::Randiv(Seq* sd, Subset* ss, DivMode dm, int rn) : divmode(dm)
{
	INT	nn = ss? ss->num: sd->many;

	rand_bit = nn;
	rand_dim = (nn + WordBits - 1) / WordBits;
	partab = 0; mcr = 0; dtree = 0;

	switch (divmode) {
	    case NON_DIV:   return;
	    case ONE_DIV:
		cycle = nn;
		one_tab(nn);
		nn = 0;
		for (MCTYPE x = 1L; x < cycle; nn++) x <<= 1;
		break;
	    case TREEDIV:
		cycle = 2 * nn - 3;
		dtree = new DistTree(sd, ss);
		tree_tab(dtree->root, nn);
		nn = 0;
		for (MCTYPE x = 1L; x < cycle; nn++) x <<= 1;
		break;
	    case ALL_DIV:
		cycle = (1 << --nn) - 1; break;
	    case PARTDIV:
	    default:
		cycle = rand_bit * rand_bit / 2;
		partab = new MCTYPE[rand_dim];
		break;
	}
	mcr = new McRand(nn, rn);
}

LRAND Randiv::randiv(int n)
{
	return (partab + n * rand_dim);
}

LRAND Randiv::nextrandiv()
{
	switch (divmode) {
	    case NON_DIV: return (0);
	    case ONE_DIV:
	    case TREEDIV:
		while ((rnbr = mcr->mcrand()) >= cycle) ;
		return (partab + rnbr * rand_dim);
	    case ALL_DIV:
		while (!(rnbr = mcr->mcrand())) ;
		rnbr += cycle;
		return (&rnbr);
	    default: break;
	}
	int	bit = rand() % (rand_bit / 2 + divmode - PARTDIV) + 1;
	if (bit > rand_bit / 2) bit = 1;
	for (int i = 0; i < rand_dim; i++) partab[i] = 0L;
	for (int i = 0; i < bit; ++i) {
	    rnbr = rand() % rand_bit;
	    partab[rnbr / WordBits] |= (1L << (rnbr % WordBits));
	}
	return (partab);
}
