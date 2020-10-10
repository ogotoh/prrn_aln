/****************************************************************************
*
*	sltree.h: Kruskal' algorithm to construct a single-linkage forest
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

#include "sltree.h"

static	int	treewd = 80;
static	int	min_memb = 2;
static	int	min_membB = 2;
static	int	max_memb = INT_MAX - 2;
static	int	max_membB = max_memb;
static	float	max_heightB = 100.;
static	bool	show_dist = false;

void set_max_memb(int maxm) {
	max_membB = maxm;
	std::swap(max_memb, max_membB);
}

void reset_max_memb() {
	std::swap(max_memb, max_membB);
}

void set_min_memb(int minm) {
	min_membB = minm;
	std::swap(min_memb, min_membB);
}

void reset_min_memb() {
	std::swap(min_memb, min_membB);
}

void set_max_height(float maxh) {
	max_heightB = maxh;
	std::swap(alprm.thr, max_heightB);
}

void reset_max_height() {
	std::swap(alprm.thr, max_heightB);
}

bool FindUnion::merge(int& x, int& y)
{
	if (dad[x = find_dad(x)] == Grad) return (false);
	if (dad[y = find_dad(y)] == Grad) return (false);
	if (x == y) return (false);
	if (npr[x] < npr[y]) { std::swap(x, y); }
	npr[x] += npr[y];
	if (npr[x] > max_memb && npr[y] >= min_memb) {
	    dad[x] = dad[y] = Grad;
	    return (false);
	}
	dad[y] = x;
	return (true);
}

void Slnode::graft(bool swp)
{
	Slnode*& ra = swp? right: left;
	Slnode*& rb = swp? left: right;
	Slnode*  a = swp? right: left;
	Slnode*  b = swp? left: right;
	int major = a->left->ndesc;
	int minor = a->right->ndesc;
	swp = minor > major;
	Slnode*  c = swp? a->right: a->left;
	Slnode*  d = swp? a->left: a->right;
	ra = c;
	rb = a;
	rb->left = d;
	rb->right = b;
	rb->ndesc = rb->nunit = d->ndesc + b->ndesc;
}

Slnode* Slnode::unpacked()
{
	if (ndesc <= max_memb) return (0);
	int major = left->ndesc;
	int minor = right->ndesc;
	bool	swp = minor > major;
	if (swp) std::swap(major, minor);
	if (minor > min_memb) return (this);
	if (major <= max_memb) return (0);
	graft(swp);
	return unpacked();
}

void Slnode::divsltree(std::vector<Slnode*>* subtrees,
	std::vector<Slnode*>* outlier)
{
	Slnode*	up = unpacked();
	if (up) {
	    up->left->divsltree(subtrees);
	    up->right->divsltree(subtrees);
	} else if (ndesc > min_memb) {
	    subtrees->push_back(this);
	} else if (outlier) {
	    outlier->push_back(this);
	}
}

void Slforest::newick_i(Slnode* node, bool newpar)
{
	int	l_gt_r = 0;

	if (node->isleaf()) {
	    if (show_dist) sprintf(str, "%s:%.2f", mname[node->tid], node->dist);
	    else	sprintf(str, "%s", mname[node->tid]);
	    fputs(str, fo);
	    clmpos += strlen(str);
	    return;
	}
	if (newpar) {
	    putc('(', fo);
	    clmpos++;
	}
	if (node->left) newick_i(node->left, l_gt_r >= 0);
	putc(',', fo);
	clmpos++;
	if (clmpos >= treewd) {
	    putc('\n', fo);
	    clmpos = 0;
	}
	if (node->right) newick_i(node->right, l_gt_r <= 0);
	if (newpar) {
	    putc(')', fo);
	    clmpos++;
	}
}

void Slforest::newick(Slnode* node)
{
	clmpos = 0;
	newick_i(node, true);
	putc(';', fo); putc('\n', fo);
}

void Slforest::sltree()
{
	fo = stdout;
	nodes = new Slnode[2 * n_vertex];
	vclear(nodes, 2 * n_vertex);
	trees = new std::vector<Slnode*>;
#if SPNRANK
	INT*	rank = new INT[OutPrm.MaxOut + 1];
	INT	max_rank = 0;
	vclear(rank, OutPrm.MaxOut + 1);
#endif
	Slnode*	i_node = nodes;
	for (int i = 0; i < n_vertex; ++i) {
	    i_node->root = i_node;
	    i_node->ndesc = i_node->nunit = 1;	// leaves
	    (i_node++)->tid = i;
	}
	Slnode*	l_node = i_node + n_vertex - 1;
	FindUnion fu(n_vertex);
	PrQueue_idx<Edge>	pq(edges, n_edge, n_edge, false);
	while (!pq.empty()) {
	    Edge*	ei = pq.shift_ptr();
	    if (ei->dist > alprm.thr) break;
	    int	x = ei->u;
	    int	y = ei->v;
	    if (fu.merge(x, y)) {
		node_union(i_node, x, y);
		i_node->dist = ei->dist;
#if SPNRANK
		if (ei->rank <= OutPrm.MaxOut) ++rank[ei->rank];
		if (ei->rank > max_rank) max_rank =  ei->rank;
#endif
		if (i_node->ndesc > max_memb) fu.graduate(x);
		if (++i_node >= l_node) break;	// connected
	    }
	}
#if SPNRANK
	for (INT i = 1; i <= min(max_rank, OutPrm.MaxOut); ++i)
	    printf("%d\t%d\n", i, rank[i]);
	delete[] rank;
#endif
}

void Slforest::node_union(Slnode* xy, int rx, int ry)
{
	xy->left = nodes[rx].root;
	xy->right = nodes[ry].root;
	xy->ndesc = xy->nunit = xy->left->ndesc + xy->right->ndesc;
	nodes[rx].root = xy;
	nodes[ry].root = 0;
}

inline bool compf(const Slnode* left, const Slnode* right)
{
	return (left->ndesc > right->ndesc);
}

int Slforest::get_trees(bool crk)
{
	if (!n_vertex) return(0);
	for (int i = 0; i < n_vertex; ++i)
	    if (nodes[i].root) trees->push_back(nodes[i].root); 
	if (crk) cruck();
	sort(trees->begin(), trees->end(), compf);
	return trees->size();
}

void Slforest::cruck() 
{
	std::vector<Slnode*>* subtree = new std::vector<Slnode*>;
	int	mm = 1;
	std::swap(min_memb, mm);
	for (TreePtrItr st = trees->begin(); st < trees->end(); ++st) {
	    if ((*st)->ndesc > max_memb) {
		(*st)->divsltree(subtree);
	    } else {
		subtree->push_back(*st);
	    }
	}
	std::swap(min_memb, mm);
	std::swap(trees, subtree);
	delete subtree;
}

void Slforest::tree_output() 
{
	for (TreePtrItr st = trees->begin(); st < trees->end(); ++st) {
	    if ((*st)->ndesc > max_memb) {
		std::vector<Slnode*> subtree;
		std::vector<Slnode*> outlier;
		(*st)->divsltree(&subtree, &outlier);
		for (TreePtrItr sst = subtree.begin(); sst < subtree.end(); ++sst) {
		   newick(*sst);
		}
	    } else {
		newick(*st);
	    }
	}
}

void sltree_defparam(int molc)
{
	algmode.lcl = 0;	// default global
	alprm.ls = 2;		// affine gap penalty
	algmode.lsg = 0;	// non-splicing alignment
	algmode.qck = 0;	// dp alignment
	algmode.blk = 1;	// block search ON
	algmode.mlt = 1;	// sigle alignment for each query
	algmode.mns = 1;	// singel strand and single direction
	algmode.thr = 1;	// filter out weak matches
	algmode.crs = 1;	// distant pairs
	setNpam(4, -6); 	// default int mismatch score
	setpam(DefPam, 0);	// change to this pam value
	setpam(WlpPam, 1);	// pam of secondary sim matrix
	OutPrm.MaxOut = 0;	// M-nearest 0: data-dependent
	OutPrm.SkipLongGap = 1;	// suppress display of long gaps
	OutPrm.supself = 1;	// suppress self comparison
	if (molc == PROTEIN) {
	    setQ4prm("k4");	// tuple
	    setQ4prm("A18");	// alphabet
	    setQ4prm("B11101,11011");
	    setQ4prm("C4");	// no bitpattern
	} else {
	    setQ4prm("k8");	// tuple
	    setQ4prm("B110011101,11101011");
	    setQ4prm("C4");	// no bitpattern
	}
}

static	int	defmolc = UNKNOWN;
static	const	char*	catalog = 0;
static	int	clmn = 1;
static	bool	sim = false;

void sltree_getoption(int& argc, const char**& argv)
{
	while (--argc > 0 && **++argv == OPTCHAR) {
	  const char*	pn = argv[0] + 2;
	  switch (argv[0][1]) {
	    case 'A':
	      if ((pn = getarg(argc, argv, true)))
		algmode.any = (int) atoi(pn);
	      break;
	    case 'D': show_dist= true; break;
	    case 'H':
	      if ((pn = getarg(argc, argv, true)))
		alprm.thr = (float) atof(pn);
	      break;
	    case 'K': 
	      if ((pn = getarg(argc, argv, false))) {
		switch (toupper(*pn)) {
		  case 'A': case 'P': defmolc = PROTEIN; break;
		  case 'D': case 'N': case 'R': defmolc = DNA; break;
		}
	      }
	    case 'M': 
	      if ((pn = getarg(argc, argv, true)))
		OutPrm.MaxOut = atoi(pn);
	      break;
	    case 'N': 
	      if ((pn = getarg(argc, argv, true)))
		min_memb = atoi(pn);
	      break;
	    case 'S':	sim = true; break;
	    case 'X': 
	      if ((pn = getarg(argc, argv, true)))
		max_memb = atoi(pn);
	      break;
	    case 'c':		// N-th column (1)
	      if ((pn = getarg(argc, argv, true)))
		clmn = atoi(pn);
	      break;
	    case 'h': usage();
	    case 'i':
	      if (*pn == ':') catalog = ++pn;
	      else if ((pn = getarg(argc, argv, false)))
		catalog = pn;
	      break;
            case 'p':
         	switch (*pn) {
                  case 'i': OutPrm.ColorEij = 1; break;
                  case 'J': case 'm': case 'u': case 'v':
			setexprm_z(argc, argv); break;
		  default: break;
		}
		break;
	    case 's':
	      if ((pn = getarg(argc, argv, false)))
		setdfn(pn);
	      break;
	    case 't': 
	      thread_num = ((pn = getarg(argc, argv, true)))?
		atoi(pn): -1;
	      break;
	    case 'x':		// x2: lower triangle only
	      if ((pn = getarg(argc, argv, true)))
		OutPrm.supself = atoi(pn);
	      break;
	    default: setQ4prm(argv[0] + 1); break;
	  }
	}
}

#if MAIN

void usage()
{
	fputs("Usage:\n\tsltree [-Hn -Nn -Xn -cn -S -D] distance file\n", stderr);
	fputs("  or\tsltree -K[A|D] [-Hn -Nn -Xn -D] fasta sequence files\n",
	   stderr);
	exit (1);
}

int main(int argc, const char** argv)
{
	sltree_defparam();
	sltree_getoption(argc, argv);
	if (!*argv) usage();
	Slforest*	slf = defmolc? new Slforest(argc, argv, defmolc, catalog):
			new Slforest(*argv, alprm.thr, clmn, sim);
	slf->get_trees();
	slf->tree_output();
	delete	slf;
	EraDbsDt();
	eraStrPhrases();
	return (0);
}

#endif	// MAIN
