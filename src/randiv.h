/*****************************************************************************
*
*	Header to randiv.c
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

#ifndef _RANDIV_H_
#define _RANDIV_H_

typedef	LONG	MCTYPE;
typedef MCTYPE* LRAND;

static	const	int	WordBits = 8 * sizeof(MCTYPE);
static	const	MCTYPE	seed0 = 0x234;

enum DivMode {NON_DIV, ONE_DIV, TREEDIV, ALL_DIV, PARTDIV};

class McRand {	// Mixed congruence method
	MCTYPE	mcoef;
	MCTYPE	mcmod;
	MCTYPE	mcval;
	bool	mrand;
public:
	McRand(int p, int rn);
	~McRand() {}
	MCTYPE	mcrand() {return (mcval = (mcoef * mcval + 1) % mcmod);}
	MCTYPE	mcrand_now() {return (mcval);}
};

class Randiv {
	LRAND	partab;
	DistTree*	dtree;
	DivMode	divmode;
	int	rand_dim;
	int	rand_bit;
	MCTYPE	rnbr;
	MCTYPE* fill_tree_tab(Knode* node);
	void	one_tab(int nn);
	void	tree_tab(Knode* root, int nn);
public:
	McRand*	mcr;
	MCTYPE	cycle;
	Randiv(Seq* sd, Subset* ss, DivMode dm, int rn = 1);
	Randiv(Seq* sd, Ssrel* srl, DivMode dm, int rn = 1);
	~Randiv();
	void	bin2lst(int* wk, LRAND lnbr);
	void	bin2lst2(int* lst[2], LRAND lnbr);
	void	bin2lst2g(int* lst[2], LRAND lnbr, Subset* ss);
	int	count1(LRAND lnbr = 0);
	LRAND	randiv(int n);
	LRAND	nextrandiv();
};

extern	void	putlst(FILE* fd, int* lst);
extern	int	setrandiv(int dm);

#endif
