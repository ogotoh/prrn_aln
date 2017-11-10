/*****************************************************************************
*
*	Collection of headers commonly used for sequence comparison
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

#ifndef _CSS_H_
#define _CSS_H

#if USE_CSRTAB
typedef struct {
	int	dim;
	RANGE**	dat;
	RANGE*** tab;
} CSRTAB;
#endif

extern	void	fprintrng(FILE* fd, RANGE* rng);
extern	RANGE*	getrng(const char* ps);
extern	int	sumrng(RANGE* rng);
extern	int	css(RANGE* rngs[2], SKL* skl[2]);
extern	void	unfoldrng(RANGE* rng, GAPS* gg);
extern	void	foldrng(RANGE* rng, GAPS* gg);
extern	RANGE*	copyrng(RANGE* s);
extern	RANGE*	cmnrng(RANGE* a, RANGE* b, int leave0);
extern	RANGE*	uniterng(RANGE* a, RANGE* b, int leave0);
extern	RANGE*	complerng(RANGE* c, RANGE* a);
extern	RANGE*	singlerng(RANGE* rng, int left, int right);

#if USE_CSRTAB
extern	void	freecsrtab (CSRTAB* ct);
extern	CSRTAB* makecsrtab (int nn);
extern	void	calccsrtab (CSRTAB* ct, GAPS** gg, GAPS** hh);
extern	CSRTAB* cmncsrtab (CSRTAB* ct1, CSRTAB* ct2, int* ga, int* gb);
extern	RANGE**	cmntorow (RANGE** clm, CSRTAB* ct);
#endif

#endif
