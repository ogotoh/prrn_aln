/*****************************************************************************
*
*	Header file for use of pair-comparison analysis
*
**	Osamu Gotoh, ph.D.	(-2001)
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

#ifndef _PAIRS_H_
#define _PAIRS_H_

enum Among {INTRA = 1, INTER, ALL, GIVE, BITREE, SPECIES};

static	const	double	def_lwcutoff = 0.;
static	const	double	def_upcutoff = 20.;
static	const	int	def_Nspc = 3;

struct	PAIR {int x, y; FTYPE w;};
class	Rnode;

struct PairPrm {
	double	upcutoff;
	double	lwcutoff;
	int	showmode;
	Among	species;
	int	Nspc;
	PairPrm() {
	    lwcutoff = def_lwcutoff;
	    upcutoff = def_upcutoff;
	    showmode = 1;
	    species = ALL;
	    Nspc = def_Nspc;
	}
	bool	pair_prm();
};

struct Pairs {
	PAIR*	pairs;
	int	npair;
	Pairs(Seq* sd, const PairPrm* pp);
	~Pairs() {delete[] pairs;}
	void	intragroup(Seq* sd);
	void	intergroup(Seq* sd);
	void	givenpair(Seq* sd);
	void	treepair(Seq* sd, const PairPrm* pp);
	void	combpair(Seq* sd, const PairPrm* pp);
	void	showpairs(Seq* sd);
	void	nodepair(Rnode* node, const PairPrm* pp);
};

class Bitree {
	Rnode**	pnode;
	Rnode*	nodes;
	int*	leaves;
public:
	Rnode*	root;
	Bitree(FTYPE* dst, int members);
	~Bitree();
};

class Rnode {
	Rnode*	left;
	Rnode*	right;
	int	id;
	int	ndesc;
	FTYPE	height;
	int*	lm;
	int*	rm;
public:
	void	extreme(int*& leaves);
friend	Bitree::Bitree(FTYPE* dst, int members);
friend	void Pairs::nodepair(Rnode* node, const PairPrm* pp);
};

#endif
