/*****************************************************************************
*
*	Management of 'static' and 'dynamic' gap states
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

#ifndef	_GFREQ_H_
#define _GFREQ_H_

struct	GFREQ	{int glen; VTYPE freq; int nres;};
struct	IDELTA	{int glen, nins;};

extern	int	warnrecovf;

inline	bool	neogfq(const GFREQ* gb) {return (gb->glen >= 0);}
inline	bool	neodelta(const IDELTA* dl) {return (dl->glen < INT_MAX);}
inline	void	resetwarnovf()	{warnrecovf = 0;}
inline	void	warnovf(int w)	{warnrecovf |= w;}

static	const	int	DGSRECOVF = 1;
static	const	int	CANRECOVF = 2;
static	const	IDELTA	ZeroDelta = {0, 0};
static	const	IDELTA	LastDelta = {INT_MAX, 0};
static	const	GFREQ	zerogfq = {0, 0, 0};
static	const	GFREQ	delmgfq = {-1, 0, 0};
static	const	GFREQ	unitgfq[2] = {{0, 1, 1}, {-1, 0, 0}};

class	Gfq {
	GFREQ*	sbuf;
	GFREQ*	tbuf;
	GFREQ*	rbuf;
	mSeq*&	msd;
	int*	lbuf;
	int	seq2gfq(int kk, CHAR* ps);
	void	store(const int& kk, const int* wk);
public:
	int	hetero;
	int	grain;
	GFREQ**	sfrq;
	GFREQ**	tfrq;
	GFREQ**	rfrq;
	GFREQ*	gbuf;
	Gfq(mSeq* sd, int gr = 1);
	Gfq(mSeq* seqs[], FTYPE* wtlst, GAPS** glist = 0);
	~Gfq() {delete[] gbuf; if (sfrq) delete[] (sfrq - 1);}
};

//	Headers to gfreq.c

inline int GapLenSD(const GFREQ* gf, const IDELTA* dl)
{
	while (gf->glen >= dl[1].glen) ++dl;
	return (gf->glen + dl->nins);
}

extern	void	printseqgfq(FILE* fd, mSeq* sd, int fb);
extern	void	copygfq(GFREQ* nf, const GFREQ* df);
extern	FTYPE	meanhetero(mSeq* sd);

//	Headers to idelta.c

extern	void	printdlt(FILE* fd, const IDELTA* dlt);
extern	void	warnovfmsg();
extern	void	copydelta(IDELTA* dlt, const IDELTA* dln);
extern	void	cleardelta(IDELTA* dlt);
extern	VTYPE	newgap(const GFREQ* cf, const GFREQ* df);
extern	VTYPE	newgap(const GFREQ* cf, const IDELTA* dlc, const GFREQ* df, const IDELTA* dld);
extern	VTYPE	newgap(const GFREQ* cf, const IDELTA* dlc, int j);
extern	VTYPE	newgap(const GFREQ* df, int i, const IDELTA* dld);
extern	void	newdelta(IDELTA* dlt, const GFREQ* df, const IDELTA* dln, int n = 1);
extern	void	incdelta(IDELTA* dlt, int n = 1);
extern	void	incdelta(IDELTA* dlt, IDELTA* dln, int n = 1);
extern	void	print_min_afq();
#endif
