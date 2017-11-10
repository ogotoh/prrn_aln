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

static	const	int	NGLIST = 7;
static	const	char*	OLR_EXT = ".olr";

//	construction of MSA with progressive methods

class Msa : public Ssrel {
	int	mtx_no;
	float	u0;
	int	lst_idx;
	int*	lst_odr;
	void	lstodr(Knode* node);
	VTYPE	nomal_pairwt();
	FTYPE	sp_score(FSTAT* fst);
public:
	mSeq*	seqs[4];
	mSeq*&	msd;
	Msa(AlnServer<mSeq>&svr);
	~Msa() {clearseq(seqs, 3);}
	int	checkss();
	int	do_job();
	void	swapprm();
	void	stuckupseq(mSeq** argsq);
	void	progressive();
	VTYPE	phyl_pwt(int dyn_aln);
	void	readgap(FILE* fd);
	void	prntgap(FILE* fd);
	void	individuallen(int* leng);
	void	prrn_main(AlnServer<mSeq>* svr);
	void	updatesrl(Ssrel* srl);
	void	makemsa(const char* gtree);
	mSeq*	prog_up(Knode* node, mSeq** argsq);
	Outlier*	findoutliers(float* sdv);
};

class Prrn {
	mSeq**	seqs;
	Ssrel*	srl;
	mSeq**	slst;
	mSeq**	sbuf;
	FTYPE*	wlst;
	FTYPE*	wbuf;
	GapsList*	glists[NGLIST];
	GapsList*&	glst;	// glists[0]
	GapsList*&	gorg;	// glists[3]
	GapsList*&	gopt;	// glists[4]
	GapsList*&	gmax;	// glists[5]
	GapsList*&	grsv;	// glists[6]
public:
	Prrn(mSeq** sqs, Ssrel* trl, VTYPE& ref);
	~Prrn() {
	    delete[] slst; delete[] sbuf;
	    delete[] wlst; delete[] wbuf;
	}
	GAPS*	gather(int k, int* gr, GapsList* gbuf);
	SKL* 	divideseq(Randiv* rdiv, int* lst[], GapsList* gsubs[]);
	VTYPE	rir(VTYPE prv);
	int	totalNoSeq(int* lst);
	VTYPE	profscore(FTYPE* pw);
};

template <class node_t>
class ProgMsa {
	mSeq*	msd;
	mSeq**	seqs;
	int	no_seqs;
	DbsDt*	dbf;
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
	    seqs(subseqs), no_seqs(0), dbf(db), verdid(did), lst_odr(0), 
	    lst_idx(0), ss(ss_), pl(0), leaves(0), alnprm(alp) {
	    if (ss) {
		leaves = new mSeq*[ss->num];
		ss->num = ss->elms = 0;
		pl = ss->pool;
	    }
	}
	ProgMsa(mSeq** sqs, int nn, ALPRM* alp = 0) : seqs(sqs), no_seqs(nn), dbf(0),
	    verdid(0), lst_idx(0), ss(0), pl(0), leaves(0), alnprm(alp)
	    {lst_odr = new int[nn];}
	ProgMsa(Btree<node_t>& ltree, const char* guidetree) :
	    seqs(0), no_seqs(0), dbf(0), verdid(0), lst_odr(0), lst_idx(0),
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
