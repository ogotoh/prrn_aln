/*****************************************************************************
*
*	Header for sequence analysis
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
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#ifndef  _MSEQ_H_
#define  _MSEQ_H_

struct	mSeq;

#include "seq.h"
#include "gfreq.h"
#include "simmtx.h"

#if !FVAL
#undef	SSHP
#endif
#if SSHP
#include "ssp.h"
class	mSeqItr;
#endif  

// attribute characters

static	const	int	decompact[6] = {nil_code, gap_code, A, C, G, T};
static	const	int	nbits[16] = {
	0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
//	0, 1, 2, 3, 4, 5, 6, 7, 8, 9, a, b, c, d, e, f
static	const	int	amblist[] = {
   A,C,A,C,G,A,G,C,G,A,C,G,T,A,T,C,T,A,C,T,G,T,A,G,T,C,G,T,A,C,G,T};
// A,C,M,M,G,R,R,S,S,V,V,V,T,W,W,Y,Y,H,H,H,K,K,D,D,D,B,B,B,N,N,N,N
// 0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1
static	const	int	ambaddr[] = {
   0, 0, 0, 1, 2, 4, 5, 7, 9, 12, 13, 15, 17, 20, 22, 25, 28};
// -  -  A  C  M  G  R  S  V  T   W   Y   H   K   D   B   N

/*******************************************************************
	structure of profile vector
	{nil, unp, [selm], 0, gap_prof, [max_code - 2], ether}
	selm = 4 | 20 | 23
	nelm = 2 + selm + max_code + ether? 1: 0
	felm = selm + 2
********************************************************************/

struct MsaState {
	INT	thk:	1;
	INT	gfq:	1;
	INT	psq:	1;
	INT	cmp:	1;
	INT	snd:	1;
};

struct SeqThk {
	VTYPE	cfq;
	VTYPE	dfq;
	VTYPE	efq;
};

struct DiThk {
	VTYPE	cc;
	VTYPE	cd;
	VTYPE	dc;
	VTYPE	dd;
	VTYPE	cx;
	VTYPE	dx;
};

static	const	SeqThk	unit_dns = {1, 0, 1};

class mSeq : public Seq {
	VTYPE*	pseq;	// frequency & profile vector
	Simmtx*	simmtx;
	MsaState	wasmst;
	CHAR**	internalres;
#if SSHP
	VTYPE	evenwt;
	bool	is_tron;
	int	sshp_step;
	int	sshp_bias;
	int	sshp_size;
	VTYPE*	sshp_prof;
	CHAR*	seq_end;
	void	ssprof(mSeqItr& asi);
	void	hyprof(mSeqItr& asi);
	void	hmprof(mSeqItr& asi, int aid);
#endif
public:
	VTYPE	sumwt;
	int	nelm;	// vector dimension
	int	felm;	// offset to profile vector
	VTYPE	height;
	VTYPE	resist;
	int	eth_code;	// background del prof
#if SSHP
	VTYPE*	sat(int i) {return sshp_prof? sshp_prof + sshpprm->sshpelems * (i - sshp_bias): 0;}
	VTYPE*	hat(int i) {return sshp_prof? sat(i) + sshpprm->sndstates: 0;}
	void	makesshpprof(bool fill = true);
	void	sshp_normalize(int i);
#endif
	mSeq**	anti_;
	Gfq*	gfq;	// gap profile
	SeqThk*	thk;	// density of residues
	DiThk*	dithk;	// di-density of residues
	CHAR*	at(int i) {return (seq_ + many * i);}
	VTYPE*	fat(int i) {return (pseq + nelm * ++i);}
	VTYPE*	pat(int i) {return (pseq + nelm * ++i + felm);}
	CHAR*	cfat(int i) {return (inex.prof? (CHAR*) fat(i): at(i));}
	CHAR*	cpat(int i) {return (CHAR*) pat(i);}
	mSeq*	refresh(const int& m = 0, const int& l = 0);
	int	unit_mode()	{return inex.prof? 2: (many == 1? 0: 1);}
	VTYPE	thickness() {return (sumwt);}
	void	mkthick();
	mSeq*	clean();
	void	back2seq();
	mSeq*	aliaseq(mSeq* dest, bool this_is_alias = false);
	mSeq*	getseq(const char* str, DbsDt* dbf = 0);
	mSeq*	copyseq(mSeq* dst, int snl = CPY_ALL);
	mSeq*	catseq(mSeq* head);
	mSeq*	cutseq(mSeq* dst, int snl);
	mSeq*	extseq(mSeq* dst, int* grp, int snl, FTYPE nfact = 1);
	mSeq*	extract(Seq* src, int* grp, int snl);
	mSeq*	fgetseq(FILE* fd, const char* attr = 0, const char* attr2 = 0);
	mSeq*	read_dbseq(DbsDt* dbf, long pos);
	mSeq*	postseq(CHAR* ss) {Seq::postseq(ss); return (this);}
	mSeq&	operator=(mSeq& src) {
	    mSeq& dest = *src.copyseq(this);
	    dest.sid = src.sid;
	    return (dest);
	}
	VTYPE	gapdensity(const CHAR* s, int i) {
	    if (*s > gap_code) return (0);
	    if (*s == gap_code || !internalres) return (1);
	    if (s < internalres[i]) return (inex.exgl? 0: alprm.tgapf);
	    else	return (inex.exgr? 0: alprm.tgapf);
	}
	VTYPE	postgapdensity(const CHAR* s, int i) {
	    if (*s == nil_code && s < internalres[i]) return (inex.exgl? 0: alprm.tgapf);
	    if (s[many] == nil_code && s >= internalres[i]) return (inex.exgr? 0: alprm.tgapf);
	    return (1);
	}
	void	setsimmtx(Simmtx* sm) {simmtx = sm;}
	void	ntor(VTYPE* v, int k, VTYPE w);
	void 	seq2vec(VTYPE* v, CHAR* s, FTYPE* w);
	void	nuc2cvec(VTYPE* v, CHAR* s, FTYPE* w);
	void	aas2cvec(VTYPE* v, CHAR* s, FTYPE* w);
	void	res2prof(VTYPE v[], CHAR r, VTYPE w);
	void	profile_n(VTYPE v[], VTYPE w[]);
	void	profile_p(VTYPE v[], VTYPE w[]);
	void	profile(VTYPE v[], VTYPE w[]);
	VTYPE	unp_prof(VTYPE *w);
	void	prepseq(bool mksp);
	mSeq*	convseq(INT vect, int step = 1);
	Seq*	consenseq(Seq* cs = 0);
	void	clear(bool snd = true) {
	    pseq = 0; gfq = 0; thk = 0; dithk = 0; simmtx = 0; 
	    anti_ = 0; sumwt = height = resist = 0; internalres = 0;
	    inex.vect = inex.prof = inex.gfrq = 0; 
	    wasmst.thk = wasmst.gfq = wasmst.psq = wasmst.cmp = wasmst.snd = 0;
#if SSHP
	    is_tron = false; sshp_step = 0;
	    sshp_prof = 0; sshp_bias = sshp_size = 0; seq_end = 0;
#endif
	}
	mSeq() : Seq(0) {clear();}	// consruct an empty object;
	mSeq(const int many, int len = 0) : Seq(many, len) {clear();}
	mSeq(const char* fname) : Seq(fname) {clear();}
	mSeq(const mSeq& msd) : Seq((Seq&) msd, 0, CPY_ALL) {clear();}
	mSeq(const mSeq* msd) : Seq((Seq*) msd, 0, CPY_ALL) {clear();}
	mSeq(Seq& sd, int* which = 0, int snl = CPY_ALL) 
	    : Seq(sd, which, snl) {clear();}
	mSeq(Seq* sd, int* which = 0, int snl = CPY_ALL) 
	    : Seq(sd, which, snl) {clear();}
	~mSeq();
};

/*******************************************************************
	Sequence iterator
*******************************************************************/

class mSeqItr : public SeqItr {
	mSeq*	msd;
	void	add_res_to_prof(mSeqItr& rsi, VTYPE wt);
public:
	int	nelm;
	int	felm;
	int	eth_code;
	VTYPE*	vss;
	int	thk_mode;
#if SSHP
	VTYPE*	sshp;
#endif
	GFREQ**	sfq;
	GFREQ**	tfq;
	GFREQ**	rfq;
const	SeqThk*	dns;
	DiThk*	didns;

	void repos()
	{
	    int	d = 0;
	    if (pos < 0) d = -1; else
	    if (pos >= msd->len - 1) d = 1;
	    dns = msd->thk? msd->thk + d: 0;
	    didns = msd->dithk? msd->dithk + d: 0;
	    if (sfq) {
		sfq = msd->gfq->sfrq + d;
		tfq = msd->gfq->tfrq + d;
		sfq = msd->gfq->rfrq + d;
	    }
	}
	mSeqItr& operator++()
	{
	    ++pos; res += many;
	    if (vss) vss += nelm;
	    if (bb) ++bb;
#if SSHP
	    if (sshp) sshp += sshpprm->sshpelems;
#endif
	    if (thk_mode == 2 || (thk_mode == 1 && (pos == 0 || pos == msd->len - 1))) {
		++dns;
		if (didns) ++didns;
		if (sfq) {
		    ++sfq; ++tfq; ++rfq;
		}
	    }
	    return (*this);
	}
	mSeqItr& operator--()
	{
	    --pos; res -= many;
	    if (bb) --bb;
	    if (vss) vss -= nelm;
#if SSHP
	    if (sshp) sshp -= sshpprm->sshpelems;
#endif
	    if (thk_mode == 2 || (thk_mode == 1 && (pos < 0 || pos == msd->len - 2))) {
		--dns;
		if (didns) --didns;
		if (sfq) {
		    --sfq; --tfq; --rfq;
		}
	    }
	    return (*this);
	}
	mSeqItr& operator+=(int n)
	{
	    pos += n; res += many * n;
	    if (bb) bb += n;
	    if (vss) vss += n * nelm;
#if SSHP
	    if (sshp) sshp += n * sshpprm->sshpelems;
#endif
	    if (thk_mode == 2) {
		dns += n;
		if (didns) didns += n;
		if (sfq) {
		    sfq += n; tfq += n; rfq += n;
		}
	    } else if (thk_mode == 1) repos();
	    return (*this);
	}
	mSeqItr& operator-=(int n)
	{
	    pos -= n; res -= many * n;
	    if (bb) bb -= n;
	    if (vss) vss -= n * nelm;
#if SSHP
	    if (sshp) sshp -= n * sshpprm->sshpelems;
#endif
	    if (thk_mode == 2) {
		dns -= n;
		if (didns) didns -= n;
		if (sfq) {
		    sfq -= n; tfq -= n; rfq -= n;
		}
	    } else if (thk_mode == 1) repos();
	    return (*this);
	}
	mSeqItr operator+(int n)
	{
	    mSeqItr	temp(*this);
	    temp.pos += n; temp.res += many * n;
	    if (bb) bb += n;
	    if (vss) temp.vss += n * nelm;
#if SSHP
	    if (sshp) temp.sshp += n * sshpprm->sshpelems;
#endif
	    if (thk_mode == 2) {
		temp.dns += n;
		if (temp.didns) temp.didns += n;
		if (sfq) {
		    temp.sfq += n; temp.tfq += n; temp.rfq += n;
		}
	    } else if (thk_mode == 1) temp.repos();
	    return (temp);
	}
	mSeqItr operator-(int n)
	{
	    mSeqItr	temp(*this);
	    temp.pos -= n;  temp.res -= many * n;
	    if (bb) bb -= n;
	    if (vss) temp.vss -= n * nelm;
#if SSHP
	    if (sshp) temp.sshp -= n * sshpprm->sshpelems;
#endif
	    if (thk_mode == 2) {
		temp.dns -= n;
		if (temp.didns) temp.didns -= n;
		if (sfq) {
		    temp.sfq -= n; temp.tfq -= n; temp.rfq -= n;
		}
	    } else if (thk_mode == 1) temp.repos();
	    return (temp);
	}
	bool	dullend();
	void	reset(int n = 0, mSeq* sd = 0);
	void	sqset(CHAR* ss, int n = 0) {SeqItr::sqset(ss, n);}
	void	insertres(mSeqItr& rsi, int bias, VTYPE wt);
	void	insertgap(mSeqItr& rsi, int bias);
	mSeqItr(mSeq* sd = 0, int n = 0) : SeqItr(sd, n), msd(sd)
		{res = 0; reset(n, sd);}
	mSeqItr(mSeqItr& src) {*this = src;}
};

// number or sum of weights of unpaired residues in gaps 
// longer than the critical gap length = k1 in an MSA column

class Gep1st {
	int	many;
#if USE_WEIGHT
	FTYPE*	weight;
	FTYPE*	pairwt;
#endif
	Queue<int>**	gep1;
public:
	Gep1st(mSeq* sd, int k1, int initialvalue = 0);
	~Gep1st();
	void	shift(const int i, const int n) {gep1[i]->shift(n);}
	void	shift(const CHAR* s, const int n);
	bool	longup(const int i, const int n, const int tally_gap_len)
		{return tally_gap_len > n - gep1[i]->oldest();}
	VTYPE	longup(const CHAR* s, const int n, const int tally_gap_len, bool sht = true);
	VTYPE	longup(const CHAR* s, const int n, const int* gap_lens);
	VTYPE	longup(GFREQ* df, IDELTA* dld, const int n);
	VTYPE	longup(GFREQ* df, IDELTA* dld, SeqItr& si);
};

/*	Headers to sqio.c	*/

inline	void	swapseq(mSeq** a, mSeq** b)
	{if (a && b) {mSeq* t = *a; *a = *b; *b = t;}}
extern	void	initseq(mSeq* seqs[], int n);
extern	void	clearseq(mSeq* seqs[], int n);
extern	void	antiseq(mSeq**);
extern	mSeq*	inputseq(mSeq* seqs[], char* ps);

#endif	// _MSEQ_H_
