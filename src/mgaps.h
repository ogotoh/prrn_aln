/*****************************************************************************
*
*	Header for abstract structures of alingmnet
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
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

#ifndef _MGAPS_H_
#define _MGAPS_H

/*	Header to mgaps.c	*/

class GapsList {
	int	size;		// capacity of records
	int	span;		// span of a record
	bool	work;		// provide work place
	GAPS*	tmp;		// temporary buffer
	GAPS*	res;		// output buffer
public:
	int	num;		// current # of records
	GAPS**	glst;		// records
	GapsList(int n, int spn = 0, bool wk = false);
	~GapsList();
	GapsList(GapsList& src, bool wk = false);
	GapsList(FILE* fd);
	GAPS*&	operator[](int i) {return (glst[i]);}
	GapsList*	resize(int spn);
	void	copy_to(GapsList* dst);
	void	copy_to(GapsList* dst, int* gr);
	void	merge_from(GapsList* src, int *gr);
	void	print(FILE* fd);
	void	insertgap(GAPS* hh);
	GAPS*	delcommongap();
};

extern	void	incrgap(int* gg, CHAR* ss, int n);
extern	void	elongap(int* gl, int* pr, CHAR* s, int n);
extern	bool	synthgap(GapsList* glists[], SKL* skl, int* lst[]);
extern	mSeq*	aggregate(mSeq* sd, mSeq** slist, GAPS** glst, FTYPE* wtlst);
extern	mSeq*	aggregate_sb(mSeq* sd, mSeq** slist, GapsList* gsrc, FTYPE* wtlst);
extern	GAPS*	gatherseq(mSeq* sd, mSeq** slst, GapsList* glst, 
		GapsList* gbuf = 0, VTYPE* wlst = 0, int* group = 0);

#endif
