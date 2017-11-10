/*****************************************************************************
*
*	Quickly search for similar blocks
*	Using oligomer compositions
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

#ifndef	_QDIV_H_
#define	_QDIV_H_

#include "bitpat.h"

struct WCPRM {
	INT	Nalpha;
	INT	Ktuple;
	DistMes	distm;
	INT	nBitPat;
const	char*	sBitPat;
const	char*	reducap;
};

class Kcomp : public Dhash<INT, int> {
	INT*	ConvTab;
	void	makeHash(Seq* sd, int i, Bitpat_wq* bpp);
public:
	int	sid;
	bool	vrtl;
	int	many;
	int	TtlKmers;
	int	UnqKmers;
	KVpair<INT, int>*	kmers;
	bool	empty() {return UnqKmers == 0;}
	Kcomp() : Dhash<INT, int>(0), ConvTab(0), kmers(0)
	   {sid = many = TtlKmers = UnqKmers = 0;}
	Kcomp(Kcomp& src) {
	    *this = src;
	    kmers = new KVpair<INT, int>[UnqKmers];
	    vcopy(kmers, src.kmers, UnqKmers);
	}
	Kcomp(Seq* sd, INT sz, INT* ct, int i, int m = 1);
	~Kcomp() {}
	void	aliaseq(Kcomp* kd, bool tia = false) {
	    *kd = *this;
	    if (tia)	vrtl = true;
	    else	kd->vrtl = true;
	}
	Kcomp* copyseq(Kcomp* kd, int mode = 0) {
	    delete kd;
	    return (kd = new Kcomp(*this));
	}
	Kcomp*	rndseq() {return (this);}
	int	fget(FILE* fd, const char* fn = 0) {return 1;}
};

class Kcomps : public ReducWord {
public:
	int	num;
	int	inodr;
	int	input_mode;
	FTYPE	vthr;
	Kcomp**	kcomp;
	Strlist*	sname;
	FTYPE*	dist;
	int*	nmidx;
	FTYPE*	getdist() { FTYPE* rv = dist; dist = 0; return (rv); }
	void	fout(int i, int j, int k);
	void	dvxout();
	Kcomps(Seq* sd, int argc, const char** argv, int nn, AlnServer<Seq>* ssvr);
	Kcomps(Seq* sd, int argc, const char** argv, const char* cat, int nn, int inm);
	Kcomps(Seq* msd, int im = IM_EVRY, Subset* ss = 0);
	Kcomps(Seq* sqs[], int nn, int im = IM_EVRY);
	~Kcomps();
};

extern	FTYPE*	calcdist_kmer(Seq* sd, int inm, Subset* ss);
extern	FTYPE*	calcdist_kmer(Seq* sd, int argc, const char** argv,
		    const char* cat, int num, int inm);
extern	FTYPE*	calcdist_kmer(Seq** seqs, int nn, int inm, Strlist*& snm);
extern	WCPRM*	setqdiv(INT k, const char* ap, const char* bp, DistMes dm);
extern	double	qdiv(Kcomp* a, Kcomp* b);
extern	WCPRM*	resetqdiv(Seq* sd);

#endif

