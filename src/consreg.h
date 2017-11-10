/*****************************************************************************
*
*	Header to consreg.c
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

#ifndef _CONSREG_H_
#define _CONSREG_H

#define	DISSIM	1

class Randiv;
class Msap;
class PwdM;

class Conserved1 {
	mSeq*	sd;
	VTYPE	vthr;
	VTYPE	fthr;
	int	length;
	float*	results;
public:
	Msap*	msp;
	Conserved1(mSeq* sd_);
	~Conserved1();
	RANGE*	cons1_nv();
	RANGE*	cons1_pf();
};

class Conserved2 {
	mSeq**	seqs;
	mSeq*	a;
	mSeq*	b;
	VTYPE	fthr;
	int	length;
	float*	results;
public:
	PwdM*	pwd;
	Conserved2(mSeq** seqs_);
	~Conserved2();
	RANGE*	cons2_ng();
	RANGE*	cons2_nv();
	RANGE*	cons2_hf();
	RANGE*	cons2_pf();
	bool	may_fuse_ng();
	bool	may_fuse_nv();
	bool	may_fuse_hf();
	bool	may_fuse_pf();
};

class ConservedT {
	int	ntid;
	Knode*	lead;
	mSeq**	seqs;
	int**	gg_son;
	int**	gg_par;
	int**	gg;
	int*	pp;
	int*	leaf;
	bool	mayfuse(mSeq** seqs);
public:
	ConservedT(mSeq** sqs, Subset* ss, Subset* tt, bool lf = false);
	~ConservedT() {delete[] leaf;}
	int	grcopy(int i);
	void	setlead(Knode* ld) {lead = ld;}
	int*	testcons(Knode* node);
	int	testssrel(Knode* node);
	Knode*	rearrange(Knode* node);
	int	grnum() {return gg - gg_son;}
	void	setntid(int id) {ntid = id;}
	void	endstamp() {*gg = 0;}
};

class Ssrel {
	void	dividseq(mSeq* sons[], mSeq* fath, 
	    Subset* ss, Randiv& rdv, int* lstbuf[]);
	void	shiftrng(RANGE* rng, int by);
	RANGE*	constwo(mSeq** seqs);
public:
	Subset*	ss;
	Ktree*	ktree;
	VTYPE*	spscr;
	FTYPE*	pairwt;
	Ssrel() : ss(0), ktree(0), spscr(0), pairwt(0) {}
	Ssrel(const Ssrel& srl) {*this = srl;}
	Ssrel(Subset* ss_, int n = 0) : ktree(0), spscr(0), pairwt(0) {
	    ss = ss_? new Subset(*ss_): (n? new Subset(n): 0);
	}
	~Ssrel() {
	    delete ss; delete ktree; delete[] spscr; delete[] pairwt;
	}
	void	clear();
	RANGE*	consreg(mSeq* sd, int dissim);
	Ssrel*	consregss(mSeq** sqs);
	mSeq**	takeapart(GapsList* glst, mSeq* sd, bool save);
#if USE_WEIGHT
	VTYPE	pairsum_ss(mSeq* sd, bool use_pw = true);
	FTYPE	sumofpwt();
	void	calcpw(Seq* sd);
	void	recalcpw() { pairwt = ktree->recalcpw(); }
	Ssrel(Seq* sd) : ss(new Subset(sd->many)), spscr(0), pairwt(0) {calcpw(sd);}
#else
	VTYPE	pairsum_ss(mSeq* sd, bool use_pw = false);
	Ssrel(Seq* sd) : ss(new Subset(sd->many)), ktree(new Ktree(sd)), spscr(0), pairwt(0)
	    {}
#endif	// USE_WEIGHT
	void	fprint(FILE* fd);
};

extern	RANGE*	bluntreg(Seq* sd, RANGE* rng);
extern	void	set_localprm(float u = -1., float v = -1., float t = -1.);

#endif
