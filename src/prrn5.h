/*****************************************************************************
*
*	prrn.h: header to prrn, Refinement of alignment by iteration
* 
*	See prrn.doc for details
* 
*	Osamu Gotoh, ph.D.      (-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.      (2001-)
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

enum {AS_INPUT, BY_TREE, SHORTER = 1, LONGER, AT_RANDOM};
enum {EIJ, MCH, FSTI, FSTD, MDLI, MDLD, LSTI, LSTD};

struct Outlier {int eij, match, ins_f, del_f, ins_m, del_m, ins_l, del_l;};

static	const	int	NGLIST = 3;
static	const	char*	OLR_EXT = ".olr";
static	const	char*	def_header = "MSA";
static	const	int	N_Stamp = 10;
static	const	int	artp_bit = INT_MIN;
static	const	int	bit_mask = INT_MAX;

//	construction of MSA with progressive methods

template <class node_t>
class ProgMsa {
	mSeq**	seqs;
	int	no_seqs;
	DbsDt*	wdbf;
	const	int*	verdid;
	int*	lst_odr;
	int	lst_idx;
	Subset*	ss;
	int*	pl;
	mSeq**	leaves;
	ALPRM*	alnprm;
public:
	int*	lstodr() {return lst_odr;}
	int	lstodr(int i) {return lst_odr[i];}
	mSeq*	getseq(node_t*);
	ProgMsa(DbsDt* db, const int* did, mSeq** subseqs, Subset* ss_ = 0, ALPRM* alp = 0) :
	    seqs(subseqs), no_seqs(0), wdbf(db), verdid(did), lst_odr(0), 
	    lst_idx(0), ss(ss_), pl(0), leaves(0), alnprm(alp) {
	    if (ss) {
		leaves = new mSeq*[ss->num];
		ss->num = ss->elms = 0;
		pl = ss->pool;
	    }
	}
	ProgMsa(mSeq** sqs, int nn, ALPRM* alp = 0) : seqs(sqs), no_seqs(nn), wdbf(0),
	    verdid(0), lst_idx(0), ss(0), pl(0), leaves(0), alnprm(alp)
	    {lst_odr = new int[nn];}
	ProgMsa(Btree<node_t>& ltree, const char* guidetree) :
	    seqs(0), no_seqs(0), wdbf(0), verdid(0), lst_odr(0), lst_idx(0),
	    ss(0), pl(0), leaves(0), alnprm(&alprm) {
	    if (!ltree.root) fatal(not_found, guidetree);
	    ltree.fill_tname();
	}
	~ProgMsa() {
	    if (leaves) {
		clearseq(leaves, ss->num);
		delete[] leaves;
	    }
	    delete[] lst_odr;
	    delete ss;
	}
	mSeq*	prog_up(node_t* node);
	Subset*	subset() {return ss;}
	mSeq**	get_leaves() {return leaves;}
};

template <class node_t>
mSeq* ProgMsa<node_t>::prog_up(node_t* node)
{
	mSeq*	msd = 0;

	if (node->isleaf()) {
	    msd = getseq(node);
	    if (!msd) fatal("Missing leaf !\n");
	} else {
	    mSeq*	sqs[3];
	    sqs[0] = prog_up(node->left);
	    sqs[1] = prog_up(node->right);
	    sqs[2] = 0;
	    sqs[0]->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
	    sqs[1]->exg_seq(algmode.lcl & 4, algmode.lcl & 8);
	    msd = align2(sqs, 0, alnprm);
	    if (!(node->left->isleaf()  && leaves)) delete sqs[0];
	    if (!(node->right->isleaf() && leaves)) delete sqs[1];
	}
	return (msd);
}

//	Iterative refinement of MSA

class IterMsa : public Ssrel {
	VTYPE	sp_score(FSTAT* fst);
	mSeq*	infc[4];
	mSeq*&  msd;
	mSeq**	seqs;
	int	no_seqs;
	mSeq*	rsv;
	ALPRM*	alnprm;
	int*	lst_odr;
	DivMode	dm;
public:
	IterMsa(mSeq* sd, mSeq** sqs, int nn, ALPRM* alp = 0, 
	    int* lstodr = 0, DivMode dm_ = TREEDIV);
	IterMsa(mSeq* sd, mSeq** sqs = 0, Subset *ss = 0, 
	    ALPRM* alp = 0, DivMode dm_ = TREEDIV);
	~IterMsa() {clearseq(infc, 3);}
	int	checkss();
	bool	preprrn(bool bestn = false);
	void	phyl_pwt();
	void	readgap(FILE* fd);
	void	prntgap(FILE* fd);
	Outlier*	findoutliers(float* sdv);
	mSeq*	msa(bool bestn = false);
};

struct ThreadArg {
	mSeq**	sqs;
	SKL**	skl;
	int**	lst;
	VTYPE	scr;
	MCTYPE	bch;
	ALPRM*	alp;
	VTYPE	pwt;
	GapsList*	gsub[2];
};

class Prrn {
	mSeq**  seqs;
	mSeq*&	msd;
	ALPRM*	alnprm;
	Ssrel*  srl;
	int	num;
	int	no_thread;
	INT	countaln;
	DivMode	dm;
	pthread_t*	handle;
	ThreadArg*	thargs;
	mSeq**	thseqs;
	SKL**	thskls;
	PwdM**	thpwds;
	int**	thlsts;
	int*	thlbuf;
	GapsList*  glists[NGLIST];
	mSeq**  slst;
	mSeq**  sbuf;
	VTYPE*  wlst;
	VTYPE*  wbuf;
#if USE_WEIGHT
	VTYPE*	wfact;
	void	childfact(Knode* node, double fact);
	VTYPE	calcfact(VTYPE* w, Knode* node);
#endif
	GapsList*	glst;   // glists[0]
	int	best_k;
	void	gather(mSeq* sd, int* gr = 0, GapsList* gbuf = 0);
	VTYPE	onecycle(SKL** skl, VTYPE pwt = 1.);
	VTYPE	onecycle(SKL** skl, int** lb, Randiv* rdiv);
	VTYPE	onecycle(SKL** skl, int** lb, int k = -1);
	VTYPE	best_of_n(SKL** skl, int** lb, Randiv* rdiv);
public:
	Prrn(mSeq** sqs, ALPRM* alp, Ssrel* trl, VTYPE& ref, 
	    bool bestn = false, DivMode dm_ = TREEDIV);
	~Prrn();
	SKL*    divideseq(Randiv* rdiv, int* lst[], GapsList* gsubs[], VTYPE* pwt, MCTYPE* bch = 0);
	SKL*    divideseq(MCTYPE mcran, int* lst[], GapsList* gsubs[], VTYPE* pwt);
	VTYPE   rir(VTYPE prv);
	int     totalNoSeq(int* lst);
};

class AddEijInfo
{
	FILE*	fd;
	StrHash<long>*	cdsmem;
	int	n_mem;
	int	n_eij;
	const	char*	coordinate;
	char	str[MAXL];
public:
	AddEijInfo(const char* fn);
	~AddEijInfo() {delete cdsmem; fclose(fd);}
	void	mksigii(Seq* sd);
	void	mksigii(Seq** seqs, int nn);
};

typedef int	(*FGET)(FILE* fd, const char* fn);

struct Sltree {
	Slnode*	root;
	int	sid;
	int	tid;
	int	vrtl;
	FGET	fget;
	bool	operator<(const Sltree& right) const {
	    return (root->ndesc > right.root->ndesc);	// descending order
	}
};

class SlfPrrn {
	int	no_trees;
	int	no_outsider;
	ALPRM	alnprm;
	const	int*	verdid;
	DbsDt*	wdbf;
	std::vector<int>	outsiders;
	std::vector<Slnode*>*	trees;
	int	no_seqs;
	void	findap(Slnode* node, std::vector<Slnode*>& subtrees);
	void	move_to_outsiders(Slnode*& node);
	int	restruct(Slnode* node);
	int	make_artps(std::vector<Slnode*>& subtrees);
public:
	mSeq**	seqs;
	mSeq*	make_msa(int molc);
	mSeq*	make_msa(Slnode* node, bool sglthrd);
	mSeq*	print_submsas();
	void	prune_subtrees();
	SlfPrrn(Slforest& slf) 
	    : no_trees(slf.get_trees(true)), no_outsider(slf.no_outsiders()),
		alnprm(alprm), verdid(slf.get_verdid()), wdbf(slf.dbf()), 
		trees(slf.trees), no_seqs(0), seqs(0)  {
		int*	osdr = slf.outsiders();
		outsiders.assign(&osdr[0], &osdr[no_outsider]);
	}
	~SlfPrrn() {}
};

class Msa {
	Ssrel*	srl;
	int	mtx_no;
	float	u0;
	int	lst_idx;
	int*	lst_odr;
	void	lstodr(Knode* node);
	mSeq*	msd;
	void	phylsort();
	FTYPE	nomal_pairwt();
public:
	Msa(mSeq* isd) : srl(0), mtx_no(1), u0(0.),
	    lst_idx(0), lst_odr(0), msd(isd) { }
	~Msa() {delete srl;}
	Outlier*	findoutliers(float* sdv);
	mSeq*	output();
};

class RunStat {
	FILE*	fmessg;
	int	timepoint;
	int	previoust;
	int	values[N_Stamp];
	time_t	timestamp[N_Stamp];
public:
	void	stamp(int val = 0);
	void	conclude();
	void	setfmessg(int& argc, const char**& argv) {
	    const	char* val = getarg(argc, argv);
	    if (val && *val) {
		if (!(fmessg = fopen(val, "w"))) fatal(no_file, val);
	    } else	fmessg = stderr;
	}
	RunStat() : fmessg(0), timepoint(0), previoust(0) {
	    vclear(values, N_Stamp);
	    vclear(timestamp, N_Stamp);
	}
	~RunStat() {if (fmessg && fmessg != stderr) fclose(fmessg);}
};

