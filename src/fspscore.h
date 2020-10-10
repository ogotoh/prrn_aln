/*****************************************************************************
*
*	 Header to fspscore.c
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

#ifndef  _FSPSCORE_H_
#define  _FSPSCORE_H_

struct SPunit {
	int	an;
	int	bn;
	SPunit(mSeq* seqs[]) {
	    an = seqs[0]->many;
	    bn = seqs[1]->many;
	}
};

struct SPunit_nv : public SPunit {
	int*	gla;
	int*	glb;
	SPunit_nv(mSeq* seqs[]) : SPunit(seqs) {
	    gla = new int[an + bn];
	    glb = gla + an;
	    seqs[0]->pregap(gla);
	    seqs[1]->pregap(glb);
	}
	~SPunit_nv() {delete[] gla;}
};

struct SPunit_w11 : public SPunit {
	int	gla;
	int	glb;
	SPunit_w11(mSeq* seqs[]) : SPunit(seqs) {
	    gla = glb = 0;
	}
	~SPunit_w11() {}
};

#if USE_WEIGHT
struct SPunit_w21 : public SPunit {
	int*	gla;
	int	glb;
	FTYPE*	wta;
	SPunit_w21(mSeq* seqs[]) : SPunit(seqs) {
	    gla = new int[an];
	    wta = seqs[0]->weight;
	    seqs[0]->pregap(gla);
	    glb = 0;
	}
	~SPunit_w21() {delete[] gla;}
};

struct SPunit_w22 : public SPunit_nv {
	FTYPE*	wta;
	FTYPE*	wtb;
	SPunit_w22(mSeq* seqs[]) : SPunit_nv(seqs) {
	    wta = seqs[0]->weight;
	    wtb = seqs[1]->weight;
	}
	~SPunit_w22() {}
};

#endif	// USE_WEIGHT

struct SPunit_hf : public SPunit {
	IDELTA*	dla;
	int	glb;
	SPunit_hf(mSeq* seqs[]) : SPunit(seqs) {
	    dla = new IDELTA[seqs[0]->gfq->hetero + 2];
	    cleardelta(dla);
	    glb = 0;
	}
	~SPunit_hf() {delete[] dla;}
};

struct SPunit_pf : public SPunit {
	IDELTA*	dla;
	IDELTA*	dlb;
	SPunit_pf(mSeq* seqs[]) : SPunit(seqs) {
	    int	ah = seqs[0]->gfq->hetero + 2;
	    dla = new IDELTA[ah + seqs[1]->gfq->hetero + 2];
	    dlb = dla + ah;
	    cleardelta(dla);
	    cleardelta(dlb);
	}
	~SPunit_pf() {delete[] dla;}
};

template <class recd_t>
class SpScore : public recd_t {
	mSeq*	a;
	mSeq*	b;
	mSeq**	seqs;
	mSeqItr	asi;
	mSeqItr	bsi;
	mSeqItr	azsi;
	mSeqItr bzsi;
	int	an;
	int	bn;
	VTYPE	vbn;
	VTYPE	scr;	// alignment score
	VTYPE	tgap;	// total number of gaps
	VTYPE	lunp;	// number of unpaired residues in long gaps
	Gep1st*	agep;
	Gep1st*	bgep;
	PwdM*	pwd;
	Gsinfo*	gsi;
	SKL*	skl;
	FSTAT*	fst;
	void	initialize();
	void	calscr(int mi, int ni);
	void	calcstat(int d3);
public:
	SpScore(mSeq* seqs[], PwdM* pwd_, FSTAT* fst_);
	SpScore(mSeq* seqs[], PwdM* pwd_, SKL* skl_);
	SpScore(mSeq* seqs[], PwdM* pwd_, Gsinfo* gsi_);
	~SpScore()  {delete agep; delete bgep;}
	VTYPE	calcSkl();
	VTYPE	calcJxt();
};

template <class recd_t>
void SpScore<recd_t>::initialize()
{
	a = seqs[0];
	b = seqs[1];
	an = a->many;
	bn = b->many;
	vbn = b->sumwt;
	agep = alprm.ls > 2? new Gep1st(a, pwd->codonk1): 0;
	bgep = alprm.ls > 2? new Gep1st(b, pwd->codonk1): 0;
	if (fst) vclear(fst);
}

template <class recd_t>
SpScore<recd_t>::SpScore(mSeq* seqs_[], PwdM* pwd_, FSTAT* fst_) : 
	recd_t(seqs_), seqs(seqs_), 
	scr(0), tgap(0), lunp(0), 
	pwd(pwd_), gsi(0), skl(0), fst(fst_)
{
	initialize();
}

template <class recd_t>
SpScore<recd_t>::SpScore(mSeq* seqs_[], PwdM* pwd_, SKL* skl_) : 
	recd_t(seqs_), seqs(seqs_), 
	scr(0), tgap(0), lunp(0), 
	pwd(pwd_), gsi(0), skl(skl_), fst(0)
{
	initialize();
}

template <class recd_t>
SpScore<recd_t>::SpScore(mSeq* seqs_[], PwdM* pwd_, Gsinfo* gsi_) : 
	recd_t(seqs_), seqs(seqs_), 
	scr(0), tgap(0), lunp(0), 
	pwd(pwd_), gsi(gsi_), skl(gsi->skl), fst(&gsi->fstat)
{
	initialize();
}

template <class recd_t>
VTYPE SpScore<recd_t>::calcJxt()
{
	int	l = a->right - a->left;

	asi.reset(a->left - 1, a);
	bsi.reset(b->left - 1, b);
	SeqThk	athk = {a->sumwt, 0, a->sumwt};
	SeqThk	bthk = {b->sumwt, 0, b->sumwt};
	azsi.dns = &athk;
	bzsi.dns = &bthk;
	Iiinfo*	iif = (SpbFact > 0 && ((a->sigII && a->sigII->pfqnum) || (b->sigII && b->sigII->pfqnum)))?
	    new Iiinfo((const Seq**) seqs, a->left, b->left): 0;
	calscr(l, l);
	if (iif) scr += iif->StoreIIinfo(a->left + l, b->left + l);
	if (fst) {
	    fst->val = scr;
	    fst->gap = tgap;
	}
	delete iif;
	return (scr - tgap * alprm.v);
}

template <class recd_t>
VTYPE SpScore<recd_t>::calcSkl()
{
	int	num = (skl++)->n;
	int	m = skl->m;
	int	n = skl->n;
	asi.reset(m - 1, a);
	bsi.reset(n - 1, b);
	SeqThk	athk = {a->sumwt, 0, a->sumwt};
	SeqThk	bthk = {b->sumwt, 0, b->sumwt};
	azsi.dns = &athk;
	bzsi.dns = &bthk;
	int	span = 0;
	Iiinfo*	iif = (SpbFact > 0 && ((a->sigII && a->sigII->pfqnum) || (b->sigII && b->sigII->pfqnum)))?
	    new Iiinfo((const Seq**) seqs, a->left, b->left, true): 0;
	while (--num) {
	    ++skl;
	    int	mi = skl->m - m;
	    int	ni = skl->n - n;
	    int	i = mi - ni;
	    int	d = (i >= 0)? ni: mi;
	    span += max(mi, ni);
	    if (!i || !mi || !ni) {
		calscr(mi, ni);
	    } else if (i > 0) {
		calscr(ni, ni);
		calscr(i, 0);
	    } else {
		calscr(mi, mi);
		calscr(0, -i);
	    }
	    if (d) {
		m += d;
		n += d;
		if (iif) scr += iif->StoreIIinfo(m, n);
	    }
	    if (i < 0)	{d = -i; n -= i;}
	    else if (i > 0)	{d = i;  m += i;}
	    else		d = 0;
	    if (d && iif) {
		d *= iif->step;
		if (i > 0)	iif->bgap += d;
		else	iif->agap += d;
		scr += iif->StoreIIinfo(m, n);
	    }
	    m = skl->m;
	    n = skl->n;
	}
	scr += pwd->wgop(tgap, lunp);
	pwd->rescale(scr, tgap, fst);
	if (gsi) gsi->SaveGsInfo(iif, span);
	delete iif;
	return scr;
}

template <class recd_t>
VTYPE grp_sttB(mSeq** seqs, PwdM* pwd, FSTAT* fst)
{
	SpScore<recd_t> sps(seqs, pwd, fst);
	return sps.calcJxt();
}

template <class recd_t>
VTYPE skl_sttB(mSeq** seqs, PwdM* pwd, SKL* skl)
{
	SpScore<recd_t> sps(seqs, pwd, skl);
	return sps.calcSkl();
}

template <class recd_t>
VTYPE skl_sttB(mSeq* seqs[], PwdM* pwd, Gsinfo* gsi)
{
	SpScore<recd_t> sps(seqs, pwd, gsi);
	return sps.calcSkl();
}

class PreSpScore {
	mSeq**	seqs;
	PwdM*	pwd;
	bool	newpwd;
	SKL*	skl;
public:
	PreSpScore(mSeq* seqs_[], PwdM* pwd_)
	    : seqs(seqs_), pwd(pwd_), newpwd(false), skl(0) {}
	PreSpScore(mSeq* seqs_[], SKL* skl_)
	    : seqs(seqs_), newpwd(true), skl(skl_) {
		pwd = new PwdM(seqs, &alprm, true);
		seqs[0]->convseq(seqs[0]->inex.prof);
		seqs[1]->convseq(seqs[1]->inex.prof);
		if (skl && pwd->swp) swapskl(skl);
	}
	~PreSpScore() {
	    if (newpwd) {
		if (pwd->swp) swapseq(seqs, seqs + 1);
		delete pwd;
	    }
	}
	VTYPE	calcSpScore(SKL* wsk = 0);
	VTYPE	calcSpScore(Gsinfo* gsi);
};

class Sptree {
	int*	group;
	int*	wkgr;
	FTYPE*	slwt;
	FTYPE*	wkwt;
	mSeq*	cur_sd;
	mSeq**	slst;
	mSeq**	slct;
	mSeq**	wksl;
	Ktree*	ktree;
	void	addleaf(Knode* node);
	void	collectleaf(mSeq* sub, mSeq* sd, Knode* node);
	void	addleaf_ss(Knode* node);
	int 	collect_ss(mSeq* sub, mSeq** slct, Knode* node);
public:
	Sptree(mSeq* sd, Ssrel* srl);
	Sptree(mSeq* sd, int num, Ktree* kt);
	~Sptree();
	VTYPE	sptree(mSeq* sd, Knode* node, bool use_pw = true);
	VTYPE	sptree_ss(mSeq* sd, Knode* node, bool use_pw = true);
};

extern	VTYPE	calcscore_grp(mSeq* seqs[], FSTAT* stt = 0);
#if USE_WEIGHT
extern	VTYPE	pairsum_ss(mSeq* sd, Ssrel* srl, bool use_pw = true);
#else
extern	VTYPE	pairsum_ss(mSeq* sd, Ssrel* srl, bool use_pw = false);
#endif
extern	VTYPE	groupsum(mSeq* sd, int usg);
#endif
