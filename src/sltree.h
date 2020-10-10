/****************************************************************************
*
*	Header to sltree (single linkage tree)
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
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japang
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#include "adjmat.h"

// Union Find

static	const	int	Grad = -2;

class FindUnion {
	std::vector<int>	dad;
	std::vector<int>	npr;	// number of progenies
	int	find_dad(int x) {
	    return dad[x] < 0 ? x : dad[x] = find_dad(dad[x]);
	}
public:
	FindUnion(int n = 0) : dad(n, -1), npr(n, 1) {}

	bool	merge(int& x, int& y);
	void	graduate(int x) {dad[find_dad(x)] = Grad;}
};

// single linkage node

class Slnode {
	Slnode*	unpacked();
	void	graft(bool swp);
public:
	Slnode*	left;
	Slnode*	right;
	Slnode*	root;
	int	tid;
	int	ndesc;
	int	nunit;
	FTYPE	dist;
	bool	isleaf() {return (!left && !right);}
	int	ndescend() {return (ndesc);}
	void	divsltree(std::vector<Slnode*>* subtrees, std::vector<Slnode*>* outlier = 0);
};

// single linkage forest 

class Slforest : public AdjacentMat {
	Slnode*	nodes;
	FILE*	fo;
	int	clmpos;
	char	str[MAXL];
	void	sltree();
	void	cruck();
	void	newick_i(Slnode* node, bool newpar);
public:
	std::vector<Slnode*>*	trees;
	Slforest(const char* fname, FTYPE distthr, int clmn, bool sim)
	    : AdjacentMat(fname, distthr, clmn, sim), 
		nodes(0), fo(0), clmpos(0), trees(0) {
	    if (n_edge) sltree();
	    else usage();
	}
	Slforest(int argc, const char** argv, int molc, 
	    const char* catalog = 0, int min_seqs = 0)
	    : AdjacentMat(argc, argv, molc, catalog, min_seqs),
		nodes(0), fo(0), clmpos(0), trees(0)
	    {if (n_edge) sltree();}
	Slforest(Seq* sd) 
	    : AdjacentMat(sd),
		nodes(0), fo(0), clmpos(0), trees(0)
	    {if (n_edge) sltree();}
	~Slforest() {
	    delete[] nodes; delete trees;
	}
	void	node_union(Slnode* xy, int rx, int ry);
	void	newick(Slnode* node);
	int	get_trees(bool crk = false);
	void	tree_output();
	const	int*	get_verdid() {return verdid;}
};

typedef std::vector<Slnode*>::iterator TreePtrItr;

extern	void	sltree_defparam(int molc = PROTEIN);
extern	void	sltree_getoption(int& argc, const char**& argv);
extern	void	set_max_memb(int m);
extern	void	reset_min_memb();
extern	void	set_min_memb(int m);
extern	void	reset_max_memb();
extern	void	set_max_height(float h);
extern	void	reset_max_height();
