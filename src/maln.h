/*****************************************************************************
*
*	Collection of headers commonly used for sequence comparison
*
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

#ifndef  _MALN_H_
#define  _MALN_H_

#include "aln.h"

struct	GFREQ;
struct	IDELTA;
class	PwdM;
class	Msap;

static	const	int	thr_gfq_11 = 15;
static	const	int	thr_gfq_12 = 30;
static	const	int	thr_gfq_21 = 8;
static	const	int	thr_gfq_22 = 14;
static	const	int	thr_gfq_23 = 60;
static  const   char warn_mess[] =
    "alg.mode (%d) is not supported!  Default mode will be used.\n";
static  const   char warn_mess2[] =
    "Sorry, this combination is not supported!\n";

struct	RVPDJ_nv : public RVPDJ {
	int	gla;
	int     glb;
	void	reseth(bool all = true) {
	    if (all) RVPDJ::reseth();
	    gla = glb = 0;
	}
	RVPDJ_nv&	operator=(const RVPDJ_nv& src) {
	    val = src.val;
	    ptr = src.ptr;
	    dir = src.dir;
	    jnc = src.jnc;
	    gla = src.gla;
	    glb = src.glb;
	    return (*this);
	}
};

struct	RVPDJ_hf : public RVPDJ {
	IDELTA* dla;
	int     glb;
	void	reseth(bool all = true) {
	    if (all) RVPDJ::reseth();
	    cleardelta(dla);
	    glb = 0;
	}
	RVPDJ_hf&	operator=(const RVPDJ_hf& src) {
	    val = src.val;
	    ptr = src.ptr;
	    dir = src.dir;
	    jnc = src.jnc;
	    copydelta(dla, src.dla);
	    glb = src.glb;
	    return (*this);
	}
};

struct	RVPDJ_pf : public RVPDJ {
	IDELTA* dla;
	IDELTA*	dlb;
	void	reseth(bool all = true) {
	    if (all) RVPDJ::reseth();
	    cleardelta(dla);
	    cleardelta(dlb);
	}
	RVPDJ_pf&	operator=(const RVPDJ_pf& src) {
	    val = src.val;
	    ptr = src.ptr;
	    dir = src.dir;
	    jnc = src.jnc;
	    copydelta(dla, src.dla);
	    copydelta(dlb, src.dlb);
	    return (*this);
	}
};

class Msap {
private:
	mSeq*	msd;
	int	an;
	INT	a_mode;
	int	base_code;
	int	codonk1;
	VTYPE	diffu;
	VTYPE	van;
	VTYPE	Basic_GOP;
const	Simmtx*	simmtx;
	ALPRM	alnprm;

	VTYPE	selfi(const mSeqItr&) const;
	VTYPE	crg1i(const int* gl, const mSeqItr&) const;
#if USE_WEIGHT
	FTYPE*	wta;
	FTYPE*	pairwt;
	VTYPE	selfw(const mSeqItr&) const;
	VTYPE	crg1w(const int* gl, const mSeqItr&) const;
#endif	// USE_WEIGHT
	int	eth_code;
	VTYPE	ab_u0;
	VTYPE	self_p(const mSeqItr&) const;
	VTYPE	selfw_p(const mSeqItr&) const;
public:
	int	calc_mode;
	VTYPE	Vab;
	VTYPE	(Msap::*sim1)(const mSeqItr&) const;
	VTYPE	(Msap::*crg1)(const int* gl, const mSeqItr&) const;

	Msap(mSeq* seq, int mode = 0, const ALPRM* ap = 0);
	~Msap() {}
	VTYPE	gop(VTYPE x) const {return (VTYPE(Basic_GOP * x));}
	VTYPE	newgap3(const GFREQ* acf, const GFREQ* bdf) const {
	    return (newgap(acf, bdf) * Basic_GOP);
	}
	VTYPE	newgap3(const GFREQ* acf, const IDELTA* dla,
		const GFREQ* bdf, const IDELTA* dlb) const {
	    return (newgap(acf, dla, bdf, dlb) * Basic_GOP);
	}
	void	selSpMode(mSeq* sd);
	VTYPE	ps_nml();
};

class PwdM : public PwdB {
	mSeq*&	a;
	mSeq*&	b;
	INT	abgfq: 1;
	INT	aprof: 1;
	INT	bprof: 1;
	INT	prof;
	int	base_code;
	VTYPE	van;
	VTYPE	vbn;
	VTYPE	Weighted_GOP;
	VTYPE	Basic_GOP;	// scale * v
	VTYPE	diff_u;		// scale * (u - u1)

	bool	advised_sim2(mSeq* seqs[]);
	VTYPE	sim00(const mSeqItr& asi, const mSeqItr& bsi) const {return (0);}
	VTYPE	sim11(const mSeqItr& asi, const mSeqItr& bsi) const {
	    return (simmtx->mtx[*asi.res][*bsi.res]);
	}
	VTYPE	sim12i(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim13(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim21i(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim22i(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim23i(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim31(const mSeqItr& asi, const mSeqItr& bsi) const {
	    return (asi.vss[asi.felm + *bsi.res]);
	}
	VTYPE	sim32i(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim33(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim33_n(const mSeqItr&, const mSeqItr&) const;

	void	stt11(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	void	stt12i(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	void	stt13(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	void	stt21i(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	void	stt22i(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	void	stt23i(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	void	stt31(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	void	stt32i(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	void	stt33(FSTAT*, const mSeqItr&, const mSeqItr&) const;

	VTYPE	unp1(const mSeqItr& asi, const mSeqItr& bsi) const {
	    return (asi.dns->cfq * bsi.dns->efq * -alnprm.u);
	}
	VTYPE	crg11(const int*, const int*, const mSeqItr&, const mSeqItr&, int) const;
	VTYPE	crg12i(const int*, const int*, const mSeqItr&, const mSeqItr&, int) const;
	VTYPE	crg21i(const int*, const int*, const mSeqItr&, const mSeqItr&, int) const;
	VTYPE	crg22i(const int*, const int*, const mSeqItr&, const mSeqItr&, int) const;
	VTYPE	frw11(const int* a_gla, const int* a_glb, 
		const int* b_gla, const int* b_glb) const {
	    return (((*a_gla > *a_glb && *b_gla <= *b_glb) ||
		 (*a_gla < *a_glb && *b_gla >= *b_glb))? -Basic_GOP: 0);
	}
	VTYPE	frw12i(const int*, const int*, const int*, const int*) const;
	VTYPE	frw21i(const int*, const int*, const int*, const int*) const;
	VTYPE	frw22i(const int*, const int*, const int*, const int*) const;

#if USE_WEIGHT
	VTYPE*	wta;
	VTYPE*	wtb;

	VTYPE	sim12w(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim21w(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim22w(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim23w(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim32w(const mSeqItr&, const mSeqItr&) const;
	void	stt12w(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	void	stt21w(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	void	stt22w(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	void	stt23w(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	void	stt32w(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	VTYPE	crg12w(const int*, const int*, const mSeqItr&, const mSeqItr&, int) const;
	VTYPE	crg21w(const int*, const int*, const mSeqItr&, const mSeqItr&, int) const;
	VTYPE	crg22w(const int*, const int*, const mSeqItr&, const mSeqItr&, int) const;
	VTYPE	frw12w(const int*, const int*, const int*, const int*) const;
	VTYPE	frw21w(const int*, const int*, const int*, const int*) const;
	VTYPE	frw22w(const int*, const int*, const int*, const int*) const;
#endif	// USE_WEIGHT

// >> Use_ether
	VTYPE	a_u0;
	VTYPE	b_u0;
	VTYPE	wa_u0;
	VTYPE	wb_u0;
	VTYPE	ab_u0;

	VTYPE	sim12i_p(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim13_p(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim21i_p(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim22i_p(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim23i_p(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim31_p(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim32i_p(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim33_p(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim33_np(const mSeqItr&, const mSeqItr&) const;
	VTYPE	unpa_p(const mSeqItr& asi, const mSeqItr& bsi) const {
	    return (asi.dns->cfq * bsi.dns->efq * -alnprm.u +
		asi.dns->dfq * a_u0 + bsi.dns->efq * b_u0);
	}
	VTYPE	unpb_p(const mSeqItr& bsi, const mSeqItr& asi) const {
	    return (bsi.dns->cfq * asi.dns->efq * -alnprm.u +
		bsi.dns->dfq * b_u0 + asi.dns->efq * a_u0);
	}

#if USE_WEIGHT
	VTYPE	sim11w_p(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim12w_p(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim21w_p(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim22w_p(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim23w_p(const mSeqItr&, const mSeqItr&) const;
	VTYPE	sim32w_p(const mSeqItr&, const mSeqItr&) const;
#endif	// USE_WEIGHT
// << Use_ether
public:
	INT	a_mode;
	INT	b_mode;
	ALPRM	alnprm;
	int	alnmode;
	int	an;
	int	bn;
	bool	swp;
	VTYPE	(PwdM::*Sim2)(const mSeqItr&, const mSeqItr&) const;
	VTYPE	(PwdM::*sim2)(const mSeqItr&, const mSeqItr&) const;
	VTYPE	(PwdM::*unpa)(const mSeqItr&, const mSeqItr&) const;
	VTYPE	(PwdM::*unpb)(const mSeqItr&, const mSeqItr&) const;
	void	(PwdM::*stt1)(FSTAT*, const mSeqItr&) const;
	void	(PwdM::*stt2)(FSTAT*, const mSeqItr&, const mSeqItr&) const;
	VTYPE	(PwdM::*crg2)(const int*, const int*, const mSeqItr&, const mSeqItr&, int) const;
	VTYPE	(PwdM::*frw2)(const int*, const int*, const int*, const int*) const;
#if SSHP
	VTYPE	sim2_sshp(const mSeqItr&, const mSeqItr&) const;
#endif
	PwdM(mSeq** seqs, const ALPRM* alp = 0, bool prf = true);
	~PwdM() {};
	int	selAlnMode(mSeq** seqs, bool renew = true);
	void	swapDvsP() {if (DvsP == 1 || DvsP == 2) DvsP = 3 - DvsP;}
	VTYPE	newgapc(const GFREQ* acf, const IDELTA* dla, int glb) const {
	    return ((dla->nins + acf->glen >= glb)?
		(VTYPE) (Weighted_GOP * acf->freq): 0);
	}
	VTYPE	newgapc(const GFREQ* acf, const IDELTA* dla, const IDELTA* dlb) const {
	    return ((dla->nins + acf->glen >= dlb->nins)?
		(VTYPE) (Weighted_GOP * acf->freq): 0);
	}
	VTYPE	newgapd(const GFREQ* adf, int glb, const IDELTA* dla) const {
	    return ((glb >= dla->nins + adf->glen)?
		(VTYPE) (Weighted_GOP * adf->freq): 0);
	}
	VTYPE	newgapd(const GFREQ* adf, const IDELTA* dlb, const IDELTA* dla) const {
	    return ((dlb->nins >= dla->nins + adf->glen)?
		(VTYPE) (Weighted_GOP * adf->freq): 0);
	}
	VTYPE	newgap1(const GFREQ* acf, const IDELTA* dla, int glb) const {
	    return (neogfq(acf)? (neogfq(acf+1)? Weighted_GOP * newgap(acf, dla, glb):
		newgapc(acf, dla, glb)): 0);
	}
	VTYPE	newgap1(const GFREQ* acf, IDELTA* dla, IDELTA* dlb) const {
	    return (neogfq(acf)? (neogfq(acf+1)? Weighted_GOP * newgap(acf, dla, dlb->nins):
		newgapc(acf, dla, dlb)): 0);
	}
	VTYPE	newgap2(const mSeqItr& asi, const CHAR* bs, const int* glb) const ;
	VTYPE	newgap2(const GFREQ* adf, int glb, const IDELTA* dla) const {
	    return (neogfq(adf)? (neogfq(adf+1)? Weighted_GOP * newgap(adf, glb, dla):
		newgapd(adf, glb, dla)): 0);
	}
	VTYPE	newgap2(const GFREQ* adf, const IDELTA* dlb, const IDELTA* dla) const {
	    return (neogfq(adf)? (neogfq(adf+1)? Weighted_GOP * newgap(adf, dlb->nins, dla):
		newgapd(adf, dlb, dla)): 0);
	}
	VTYPE	newgap3(const GFREQ* acf, const GFREQ* bdf) const {
	    return (newgap(acf, bdf) * Basic_GOP);
	}
	VTYPE	newgap3(const GFREQ* acf, const IDELTA* dla, 
		const GFREQ* bdf, const IDELTA* dlb) const {
	    return (newgap(acf, dla, bdf, dlb) * Basic_GOP);
	}
	VTYPE	vgop(VTYPE x) const {return (Basic_GOP * x);}
	VTYPE	wgop(VTYPE v, VTYPE lu = 0) const {
	    bool	half = alnmode == HLF_ALB || alnmode == HLF_ALN
			    || alnmode == RHF_ALB || alnmode == RHF_ALN;
	    return (v * (half? Weighted_GOP: Basic_GOP) + diff_u * lu);
	}
	void	reset_ab(mSeq** sqs) {a = sqs[0]; b = sqs[1];}
	void	resetuab();
	void	rescale(VTYPE& scr, VTYPE tgap, FSTAT* fst);
};

extern	int	getoutmode();

/*	Header to hsim2.cc	*/

extern	VTYPE	HomSelf(mSeq* sd);

/*	Header to fwd2.cc	*/

extern	void	convertseqs(mSeq** seqs, int nn, Simmtx* sm);
extern	int	selAlnMode(mSeq** seqs, PwdM* pwa);
extern	VTYPE	selfAlnScr(mSeq* sd, Simmtx* sm = 0, FSTAT* stt = 0);
extern	VTYPE	HomScore(mSeq* seqs[], PwdM* pwd, long rr[] = 0);
extern	VTYPE	HomScoreS_ng(mSeq* seqs[], PwdM* pwdm);
extern	SKL*	align2(mSeq* seqs[], PwdM* pwa, VTYPE* scr, Gsinfo* gsi = 0);
extern	mSeq*	align2(mSeq* seqs[], mSeq* msd = 0, ALPRM* alp = 0);
extern	mSeq*	align2(mSeq* a, mSeq* b, mSeq* msd = 0, ALPRM* alp = 0);
extern	mSeq*	syntheseq(mSeq* sd, mSeq** sqs, SKL* skl);
extern	FTYPE	alnscore2dist(mSeq* sqs[], PwdM* pwd, int* ends = 0, FTYPE denom = 0);
extern	SKL*	alignS_ng(mSeq* seqs[], PwdM* pwd, VTYPE* scr);
extern	SKL*	swg2ndH_ng(mSeq* seqs[], PwdM* pwa, VTYPE* scr, COLONY* clny);
extern	SKL*	swg2ndS_ng(mSeq* seqs[], PwdM* pwa, VTYPE* scr, COLONY* clny);
extern  SKL*	swg2nd(mSeq* seqs[], PwdM* pwd, Gsinfo* scr, COLONY* clny);
extern	Gsinfo*	localB_ng(mSeq* seqs[], PwdM* pwa);
extern  Colonies*	swg1st(mSeq* seqs[], PwdM* pwa);
extern  Colonies*	swg1stB_ng(mSeq* seqs[], PwdM* pwa, VTYPE* scr);
extern	Colonies*	swg1stH_ng(mSeq* seqs[], PwdM* pwa, VTYPE* scr);
extern	Colonies*	swg1stS_ng(mSeq* seqs[], PwdM* pwa, VTYPE* scr);
extern	VTYPE	skl_rngH_hf(mSeq* seqs[], Gsinfo* gsi, PwdM* pwa);
extern	VTYPE	skl_rngH_ng(mSeq* seqs[], Gsinfo* gsi, PwdM* pwd);
extern	VTYPE	skl_rngS_ng(mSeq* seqs[], Gsinfo* gsi, PwdM* pwd);
extern	void	freeCodeReg();
extern	void	report_ild(mSeq* seqs[], ISLAND* ild, INT n, int ud);
extern	VTYPE	pairsum(mSeq* sd);

#endif
