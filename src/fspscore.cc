/*****************************************************************************
*
*	Calculate Sum-of-Pairs score --- Floating point version
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

#include "aln.h"
#include "mseq.h"
#include "maln.h"
#include "mgaps.h"
#include "gfreq.h"
#include "phyl.h"
#include "consreg.h"
#include "fspscore.h"

#ifdef SSPH
#include "ssp.h"
#endif

#define	ST_PS_ALG = 0;

int	ndesc_thr = thr_gfq_23;

#ifdef DEBUG
static	void	printstat(FSTAT* st);
static	void	prval(PVTYPE x);
#endif
static	VTYPE	ps_lst(mSeq** slst, FTYPE* pw = 0);
static	mSeq*	fuseseq(mSeq* sd, mSeq** sqs, FTYPE wt = 0);

#ifdef DEBUG
static void printstat(FSTAT* st)
{
    printf("%6.1f %5.1f %5.1f %5.1f %5.1f\n", st->val, st->mch, st->mmc, 
	st->gap, st->unp);
}

static void prval(PVTYPE x)
{
    if (NEG_INT + 20 < x && x < POS_INT - 20) printf(" %5.1lf", (double) x);
    else fputs("   *  ", stdout);
}
#endif

template<>
void SpScore<SPunit_nv>::calcstat(int d3)
{
	int	mch = 0;
	int	mmc = 0;
	VTYPE	unp = 0;
	int*	ag = gla;
	int*	bg = glb;
	CHAR*	as = asi.res;
	CHAR*	bs = bsi.res;

	if (d3 == 0) {
	    scr += (pwd->*pwd->sim2)(asi, bsi);
	    for (int i = 0; i < an; ++i, ++as, ++ag) {
		bool	ar = IsntGap(*as);
		VTYPE	au = a->gapdensity(as, i);
		int*	cg = bg;
		CHAR*	cs = bs;
		for (int j = 0; j < bn; ++j, ++cs, ++cg) {
		    bool	br = IsntGap(*cs);
		    VTYPE	bu = b->gapdensity(cs, j);
		    if (ar && br) {		// * *
			if (*as == *cs)	++mch;
			else		++mmc;
		    } else if (ar && bu > 0) {	//* -
			unp += bu;
			if (*cg <= *ag) tgap += bu; else 
			if (agep && agep->longup(i, asi.pos, *cg + 1)) lunp += bu;
		    } else if (au > 0 && br) {	// - *
			unp += au;
			if (*ag <= *cg) tgap += au; else 
			if (bgep && bgep->longup(j, bsi.pos, *ag + 1)) lunp += au;
		    }
		}
		if (ar && agep) agep->shift(i, asi.pos);
	    }
	    if (bgep) bgep->shift(bsi.res, bsi.pos);
	} else if (d3 > 0) {	// * -
	    scr += (pwd->*pwd->unpa)(asi, bsi);
	    unp += asi.dns->cfq * bsi.dns->efq;
	    for (int i = 0; i < an; ++i, ++as, ++ag) {
		if (IsntGap(*as)) {
		    int*	cg = bg;
		    CHAR*	cs = bs;
		    for (int j = 0; j < bn; ++j, ++cs, ++cg) {
			VTYPE	bu = b->postgapdensity(cs, j);
			if (bu > 0) {
			    if (*cg <= *ag) tgap += bu; else 
			    if (agep && agep->longup(i, asi.pos, *cg + 1))
				lunp += bu;
			}
		    }
		    if (agep) agep->shift(i, asi.pos);
		}
	    }
	} else {		// - *
	    scr += (pwd->*pwd->unpb)(bsi, asi);
	    unp += asi.dns->efq * bsi.dns->cfq;
	    for (int j = 0; j < bn; ++j, ++bs, ++bg) {
		if (IsntGap(*bs)) {
		    int*	cg = ag;
		    CHAR*	cs = as;
		    for (int i = 0; i < an; ++i, ++cs, ++cg) {
			VTYPE	au = a->postgapdensity(cs, i);
			if (au > 0) {
			    if (*cg <= *bg) tgap += au; else
			    if (bgep && bgep->longup(j, bsi.pos, *cg + 1))
				lunp += au;
			}
		    }
		    if (bgep) bgep->shift(j, bsi.pos);
		}
	    }
	}
	if (fst) {
	    fst->mch += mch;
	    fst->mmc += mmc;
	    fst->unp += unp;
	}
}

template<>
void SpScore<SPunit_w11>::calcstat(int d3)
{
	CHAR*	as = asi.res;
	CHAR*	bs = bsi.res;

	if (d3 == 0) {
	    scr += (pwd->*pwd->sim2)(asi, bsi);
	    int	ar = IsntGap(*as);
	    int	br = IsntGap(*bs);
	    VTYPE	au = a->gapdensity(as, 0);
	    VTYPE	bu = b->gapdensity(bs, 0);
	    if (ar && br) {
		if (fst) {
		    if (*as == *bs)	++fst->mch;
		    else		++fst->mmc;
		}
	    } else if (ar && bu > 0) {
		if (fst) fst->unp += bu;
		if (glb <= gla) tgap += bu; else
		if (agep) lunp += agep->longup(as, asi.pos, glb + 1) * bu;
	    } else if (au > 0 && br) {
		if (fst) fst->unp += au;
		if (gla <= glb) tgap += au; else
		if (bgep) lunp += bgep->longup(bs, bsi.pos, gla + 1) * au;
	    }
	} else if (d3 > 0) {
	    scr += (pwd->*pwd->unpa)(asi, bsi);
	    if (IsntGap(*as)) {
		const	VTYPE&	bu = bsi.dns->efq;
		if (bu > 0) {
		    if (fst) fst->unp += bu;
		    if (glb <= gla) tgap += bu; else
		    if (agep) lunp += agep->longup(as, asi.pos, glb + 1) * bu;
		}
	    }
	} else {
	    scr += (pwd->*pwd->unpb)(bsi, asi);
	    if (IsntGap(*bs)) {
		const	VTYPE&	au = asi.dns->efq;
		if (au > 0) {
		    if (fst) fst->unp += au;
		    if (gla <= glb) tgap += au; else 
		    if (bgep) lunp += bgep->longup(bs, bsi.pos, gla + 1) * au;
		}
	    }
	}
}

#if USE_WEIGHT

template<>
void SpScore<SPunit_w21>::calcstat(int d3)
{
	int*	ag = gla;
	CHAR*	as = asi.res;
	CHAR*	bs = bsi.res;
	int	br = IsntGap(*bs);
	VTYPE	bu = bsi.dns->efq;

	if (d3 == 0) {
	    scr += (pwd->*pwd->sim2)(asi, bsi);
	    FTYPE* wa = wta;
	    for (int i = 0; i < an; ++i, ++as, ++ag, ++wa) {
		int	ar = IsntGap(*as);
		VTYPE	au = a->gapdensity(as, i);
		if (ar && br) {
		    if (fst) {
			if (*as == *bs)	fst->mch += *wa;
			else		fst->mmc += *wa;
		    }
		} else if (ar && bu > 0) {
		    if (fst) fst->unp += *wa * bu;
		    if (glb <= *ag) tgap += *wa * bu; else
		    if (agep && agep->longup(i, asi.pos, glb + 1)) lunp += *wa * bu;
		} else if (au > 0 && br) {
		    if (fst) fst->unp += *wa * au;
		    if (*ag <= glb) tgap += *wa * au; else
		    if (bgep && bgep->longup(0, bsi.pos, *ag + 1)) lunp += *wa * au;
		}
		if (ar && agep) agep->shift(i, asi.pos);
	    }
	    if (br && bgep) bgep->shift(0, bsi.pos);
	} else if (d3 > 0) {
	    scr += (pwd->*pwd->unpa)(asi, bsi);
	    FTYPE*	wa = wta;
	    for (int i = 0; i < an; ++i, ++as, ++ag, ++wa) {
		if (IsntGap(*as)) {
		    if (bu > 0) {
		        if (fst) fst->unp += *wa * bu;
		        if (glb <= *ag) tgap += *wa * bu; else
		        if (agep && agep->longup(i, asi.pos, glb + 1)) lunp += *wa * bu;
		    }
		    if (agep) agep->shift(i, asi.pos);
		}
	    }
	} else {
	    scr += (pwd->*pwd->unpb)(bsi, asi);
	    if (IsntGap(*bs)) {
		FTYPE* wa = wta;
		for (int i = 0; i < an; ++i, ++as, ++ag, ++wa) {
		    VTYPE	au = a->postgapdensity(as, i);
		    if (au > 0) {
			if (fst) fst->unp += *wa * au;
		        if (*ag <= glb) tgap += *wa * au; else
		        if (bgep && bgep->longup(0, bsi.pos, *ag + 1)) lunp += *wa * au;
		    }
		}
		if (bgep) bgep->shift(0, bsi.pos);
	    }
	}
}

template<>
void SpScore<SPunit_w22>::calcstat(int d3)
{
	int*	ag = gla;
	int*	bg = glb;
	CHAR*	as = asi.res;
	CHAR*	bs = bsi.res;

	if (d3 == 0) {
	    scr += (pwd->*pwd->sim2)(asi, bsi);
	    FTYPE* wa = wta;
	    for (int i = 0; i < an; ++i, ++as, ++ag, ++wa) {
		int	ar = IsntGap(*as);
		VTYPE	au = a->gapdensity(as, i);
		FTYPE*	wb = wtb;
		int*	cg = bg;
		CHAR*	cs = bs;
		VTYPE	s = 0, g = 0, u = 0, l = 0;
		for (int j = 0; j < bn; ++j, ++cs, ++cg, ++wb) {
		    int	br = IsntGap(*cs);
		    VTYPE	bu = b->gapdensity(cs, j);
		    if (ar && br) {
			if (*as == *cs)	s += *wb;
		    } else if (ar && bu > 0) {
			u += *wb * bu;
			if (*cg <= *ag) g += *wb * bu; else
			if (agep && agep->longup(i, asi.pos, *cg + 1)) l += *wb * bu;
		    } else if (au > 0 && br) {
			u += *wb * au;
			if (*ag <= *cg) g += *wb * au; else
			if (bgep && bgep->longup(j, bsi.pos, *ag + 1)) l += *wb * au;
		    }
		}
		tgap += g * *wa;
		lunp += l * *wa;
		if (fst) {
		    fst->mch += s * *wa;
		    fst->mmc += (vbn - s - u) * *wa;
		    fst->unp += u * *wa;
		}
		if (ar && agep) agep->shift(i, asi.pos);
	    }
	    if (bgep) bgep->shift(bsi.res, bsi.pos);
	} else if (d3 > 0) {
	    scr += (pwd->*pwd->unpa)(asi, bsi);
	    if (fst) fst->unp += bsi.dns->efq * asi.dns->cfq;
	    FTYPE* 	wa = wta;
	    for (int i = 0; i < an; ++i, ++as, ++ag, ++wa) {
		if (IsntGap(*as)) {
		    VTYPE	g = 0, l = 0;
		    int*	cg = bg;
		    FTYPE*	wb = wtb;
		    CHAR*	cs = bs;
		    for (int j = 0; j < bn; ++j, ++cs, ++cg, ++wb) {
			VTYPE	bu = b->postgapdensity(cs, j);
			if (bu > 0) {
			    if (*cg <= *ag) g += *wb * bu; else
			    if (agep && agep->longup(i, asi.pos, *cg + 1)) l += *wb * bu;
			}
		    }
		    tgap += g * *wa;
		    lunp += l * *wa;
		    if (agep) agep->shift(i, asi.pos);
		}
	    }
	} else {
	    scr += (pwd->*pwd->unpb)(bsi, asi);
	    if (fst) fst->unp += asi.dns->efq * bsi.dns->cfq;
	    FTYPE* wb = wtb;
	    for (int j = 0; j < bn; ++j, ++bs, ++bg, ++wb) {
		if (IsntGap(*bs)) {
		    VTYPE	g = 0, l = 0;
		    int*	cg = ag;
		    FTYPE*	wa = wta;
		    CHAR*	cs = as;
		    for (int i = 0; i < an; ++i, ++cs, ++cg, ++wa) {
			VTYPE	au = a->postgapdensity(cs, i);
			if (au > 0) {
			    if (*cg <= *bg) g += *wa * au; else
			    if (bgep && bgep->longup(j, bsi.pos, *cg + 1)) l += *wa * au;
			}
		    }
		    tgap += g * *wb;
		    lunp += l * *wb;
		    if (bgep) bgep->shift(j, bsi.pos);
		}
	    }
	}
}

#endif	// USE_WEIGHT

template<>
void SpScore<SPunit>::calscr(int mi, int ni)
{
	if (mi == ni) {
	    while (mi--) {
		scr += (pwd->*pwd->sim2)(++asi, ++bsi);
		if (fst) (pwd->*pwd->stt2)(fst, asi, bsi);
	    }
	} else if (mi) {
	    if (fst) fst->unp += mi;
	    tgap += pwd->Vab * bsi.dns->efq;
	    scr += pwd->UnpPenalty(mi) * bsi.dns->efq;
	    asi += mi;
	} else if (ni) {
	    if (fst) fst->unp += ni;
	    tgap += pwd->Vab * asi.dns->efq;
	    scr += pwd->UnpPenalty(ni) * asi.dns->efq;
	    bsi += ni;
	}
}

template<>
void SpScore<SPunit_nv>::calscr(int mi, int ni)
{
	if (mi == ni) {
	    while (mi--) {
		++asi; ++bsi;
		calcstat(0);
		incrgap(gla, asi.res, an);
		incrgap(glb, bsi.res, bn);
	    }
	} else if (mi) {
	    while (mi--) {
		++asi;
		calcstat(1);
		incrgap(gla, asi.res, an);
		incrgap(glb, 0, bn);
	    }
	} else if (ni) {
	    while (ni--) {
		++bsi;
		calcstat(-1);
		incrgap(gla, 0, an);
		incrgap(glb, bsi.res, bn);
	    }
	}
}

template<>
void SpScore<SPunit_w11>::calscr(int mi, int ni)
{
	if (mi == ni) {
	    while (mi--) {
		++asi; ++bsi;
		calcstat(0);
		gla = glb = 0;
	    }
	} else if (mi) {
	    while (mi--) {
		++asi;
		calcstat(1);
		gla = 0; ++glb;
	    }
	} else if (ni) {
	    while (ni--) {
		++bsi;
		calcstat(-1);
		++gla; glb = 0;
	    }
	}
}

#if USE_WEIGHT

template<>
void SpScore<SPunit_w21>::calscr(int mi, int ni)
{
	if (mi == ni) {
	    while (mi--) {
		++asi; ++bsi;
		calcstat(0);
		incrgap(gla, asi.res, an);
		glb = 0;
	    }
	} else if (mi) {
	    while (mi--) {
		++asi;
		calcstat(1);
		incrgap(gla, asi.res, an);
		++glb;
	    }
	} else if (ni) {
	    while (ni--) {
		++bsi;
		calcstat(-1);
		incrgap(gla, 0, an);
		glb = 0;
	    }
	}
}

template<>
void SpScore<SPunit_w22>::calscr(int mi, int ni)
{
	if (mi == ni) {
	    while (mi--) {
		++asi; ++bsi;
		calcstat(0);
		incrgap(gla, asi.res, an);
		incrgap(glb, bsi.res, bn);
	    }
	} else if (mi) {
	    while (mi--) {
		++asi;
		calcstat(1);
		incrgap(gla, asi.res, an);
		incrgap(glb, 0, bn);
	    }
	} else if (ni) {
	    while (ni--) {
		++bsi;
		calcstat(-1);
		incrgap(gla, 0, an);
		incrgap(glb, bsi.res, bn);
	    }
	}
}

#endif	// USE_WEIGHT

template<>
void SpScore<SPunit_hf>::calscr(int mi, int ni)
{
	if (mi == ni) {
	    while (mi--) {
		scr += (pwd->*pwd->sim2)(++asi, ++bsi);
		tgap += newgap(*asi.tfq, glb, dla);
		if (fst) (pwd->*pwd->stt2)(fst, asi, bsi);
		newdelta(dla, *asi.tfq, dla);
		if (bgep) lunp += bgep->longup(*asi.tfq, dla, bsi.pos);
		if (agep) agep->shift(asi.res, asi.pos);
		glb = 0;
	    }
	} else if (mi) {
	    while (mi--) {
		scr += (pwd->*pwd->unpa)(++asi, bsi);
		tgap += newgap(*asi.sfq, dla, glb);
		if (fst) (pwd->*pwd->stt2)(fst, asi, bzsi);
		newdelta(dla, *asi.tfq, dla);
		++glb;
		if (agep) lunp += agep->longup(asi.res, asi.pos, glb);
	    }
	} else if (ni) {
	    while (ni--) {
		scr += (pwd->*pwd->unpb)(++bsi, asi);
		tgap += newgap(*asi.rfq, glb, dla);
		if (fst) (pwd->*pwd->stt2)(fst, azsi, bsi);
		incdelta(dla, dla);
		if (bgep) lunp += bgep->longup(*asi.rfq, dla, bsi.pos);
	    }
	}
}

template<>
void SpScore<SPunit_pf>::calscr(int mi, int ni)
{
	if (mi == ni) {
	    while (mi--) {
		scr += (pwd->*pwd->sim2)(++asi, ++bsi);
		tgap += newgap(*asi.sfq, dla, *bsi.tfq, dlb) +
		     newgap(*bsi.sfq, dlb, *asi.tfq, dla);
		if (fst) (pwd->*pwd->stt2)(fst, asi, bsi);
		newdelta(dla, *asi.tfq, dla);
		newdelta(dlb, *bsi.tfq, dlb);
		if (agep) lunp += agep->longup(*bsi.tfq, dlb, asi);
		if (bgep) lunp += bgep->longup(*asi.tfq, dla, bsi);
	    }
	} else if (mi) {
	    while (mi--) {
		scr += (pwd->*pwd->unpa)(++asi, bsi);
		tgap += newgap(*asi.sfq, dla, *bsi.rfq, dlb);
		if (fst) (pwd->*pwd->stt2)(fst, asi, bzsi);
		newdelta(dla, *asi.tfq, dla);
		incdelta(dlb, dlb);
		if (agep) lunp += agep->longup(*bsi.rfq, dlb, asi);
	    }
	} else if (ni) {
	    while (ni--) {
		scr += (pwd->*pwd->unpb)(++bsi, asi);
		tgap += newgap(*bsi.sfq, dlb, *asi.rfq, dla);
		if (fst) (pwd->*pwd->stt2)(fst, azsi, bsi);
		newdelta(dlb, *bsi.tfq, dlb);
		incdelta(dla, dla);
		if (bgep) lunp += bgep->longup(*asi.rfq, dla, bsi);
	    }
	}
}

VTYPE PreSpScore::calcSpScore(SKL* wsk)
{
	if (wsk) skl = wsk;
	int	lastp = skl->n;
	SKL	rsv = skl[lastp];
	VTYPE	scr = 0;
	mSeq*&	a = seqs[0];
#if USE_WEIGHT
	mSeq*&	b = seqs[1];
#endif	// USE_WEIGHT

	switch (pwd->alnmode) {
	  case NGP_ALB:
	  case NGP_ALN: scr = skl_sttB<SPunit>(seqs, pwd, skl); break;
	  case NTV_ALB:
	  case NTV_ALN:
#if USE_WEIGHT
	   if (a->weight && b->weight) {
	    if (a->many == 1) scr = skl_sttB<SPunit_w11>(seqs, pwd, skl); else
	    if (b->many == 1) scr = skl_sttB<SPunit_w21>(seqs, pwd, skl);
	    else	scr = skl_sttB<SPunit_w22>(seqs, pwd, skl);
	   }
#else
	    if (a->many == 1) scr = skl_sttB<SPunit_w11>(seqs, pwd, skl);
#endif	// USE_WEIGHT
	    else	scr = skl_sttB<SPunit_nv>(seqs, pwd, skl);
	   break;
	  case RHF_ALB:
	  case RHF_ALN:
	  case HLF_ALB:
	  case HLF_ALN: scr = skl_sttB<SPunit_hf>(seqs, pwd, skl); break;
	  case GPF_ALB: scr = skl_sttB<SPunit_pf>(seqs, pwd, skl); break;
	  default:
	    fprintf(stderr, "Score Mode %d is not supported!\n", pwd->alnmode);
	    return (0);
	}
	skl[lastp] = rsv;
	return (scr);
}

VTYPE PreSpScore::calcSpScore(Gsinfo* gsi)
{
	SKL*	skl = gsi->skl;
	int	lastp = skl->n;
	SKL	rsv = skl[lastp];
	VTYPE	scr = 0;
	mSeq*&	a = seqs[0];
#if USE_WEIGHT
	mSeq*&	b = seqs[1];
#endif	// USE_WEIGHT

	switch (pwd->alnmode) {
	  case NGP_ALB:
	  case NGP_ALN: scr = skl_sttB<SPunit>(seqs, pwd, gsi); break;
	  case NTV_ALB:
	  case NTV_ALN:
#if USE_WEIGHT
	   if (a->weight && b->weight) {
	    if (a->many == 1) scr = skl_sttB<SPunit_w11>(seqs, pwd, gsi); else
	    if (b->many == 1) scr = skl_sttB<SPunit_w21>(seqs, pwd, gsi);
	    else	scr = skl_sttB<SPunit_w22>(seqs, pwd, gsi);
	   }
#else
	    if (a->many == 1) scr = skl_sttB<SPunit_w11>(seqs, pwd, gsi);
#endif	// USE_WEIGHT
	    else	scr = skl_sttB<SPunit_nv>(seqs, pwd, gsi);
	   break;
	  case RHF_ALB:
	  case RHF_ALN:
	  case HLF_ALB:
	  case HLF_ALN: scr = skl_sttB<SPunit_hf>(seqs, pwd, gsi); break;
	  case GPF_ALB: scr = skl_sttB<SPunit_pf>(seqs, pwd, gsi); break;
	  default:
	    fprintf(stderr, "Score Mode %d is not supported!\n", pwd->alnmode);
	    return (0);
	}
	skl[lastp] = rsv;
	return (scr);
}

VTYPE calcscore_grp(mSeq* seqs[], FSTAT* stt)
{
	VTYPE	scr = 0;
	PwdM	pwd(seqs, 0, !stt);
	mSeq*&	a = seqs[0];
#if USE_WEIGHT
	mSeq*&	b = seqs[1];
#endif	// USE_WEIGHT

	switch (pwd.alnmode) {
	  case NGP_ALB:
	  case NGP_ALN: scr = grp_sttB<SPunit>(seqs, &pwd, stt); break;
	  case NTV_ALB:
	  case NTV_ALN: 
#if USE_WEIGHT
	   if (a->weight && b->weight) {
	    if (a->many == 1) scr = grp_sttB<SPunit_w11>(seqs, &pwd, stt); else
	    if (b->many == 1) scr = grp_sttB<SPunit_w21>(seqs, &pwd, stt);
	    else	scr = grp_sttB<SPunit_w22>(seqs, &pwd, stt);
	   }
#else
	    if (a->many == 1) scr = grp_sttB<SPunit_w11>(seqs, &pwd, stt);
#endif	// USE_WEIGHT
	    else	scr = grp_sttB<SPunit_nv>(seqs, &pwd, stt);
	    break;
	  case RHF_ALB:
	  case RHF_ALN:
	  case HLF_ALB:
	  case HLF_ALN: scr = grp_sttB<SPunit_hf>(seqs, &pwd, stt); break;
	  case GPF_ALB:
	  default: scr = grp_sttB<SPunit_pf>(seqs, &pwd, stt); break;
	}
	if (stt) stt->val = scr;
	if (pwd.swp) swapseq(seqs, seqs + 1);
	return (scr);
}

static VTYPE ps_lst(mSeq** slst, FTYPE* pw)
{
	mSeq*	sqs[2];
	VTYPE	scr = 0;

#ifdef DEBUG
	printf("fscore\n");
#endif

	for (mSeq** s2 = slst + 1; *s2; ++s2) {
	    sqs[1] = *s2;
	    for (mSeq** s1 = slst; s1 < s2; ++s1) {
		sqs[0] = *s1;
		FTYPE	pwv = pw? *pw++: 1;
#ifdef DEBUG
		VTYPE	spg = calcscore_grp(sqs);
		printf("%15.7e %10.2f %10.2f %3d %3d\n",
		pwv, spg, pwv * spg, sqs[0]->many, sqs[1]->many);
		scr += pwv * spg;
#else
		scr += pwv * calcscore_grp(sqs);
#endif
	    }
	}
	return ((VTYPE) scr);
}

static mSeq* fuseseq(mSeq* sd, mSeq** sqs, FTYPE wt)
{
	mSeq*&	a = sqs[0];
	mSeq*&	b = sqs[1];

	sd->refresh(a->many + b->many, a->len);
	sd->inex.molc = a->inex.molc;
#if USE_WEIGHT
	if (a->weight && b->weight) {
	    sd->weight = new FTYPE[sd->many];
	    FTYPE*	w = sd->weight;
	    for (int i = 0; i < a->many; ++i)
		*w++ = wt * a->weight[i];
	    for (int i = 0; i < b->many; ++i)
		*w++ = wt * b->weight[i];
	}
#endif	// USE_WEIGHT
	CHAR*	ss = sd->at(0);
	CHAR*	as = a->at(0);
	CHAR*	bs = b->at(0);
	CHAR*	ts = a->at(a->len);
	while (as < ts) {
	    memcpy(ss, as, a->many);
	    memcpy(ss += a->many, bs, b->many);
	    ss += b->many;
	    as += a->many;
	    bs += b->many;
	}
	return (sd->postseq(ss));
}

void Sptree::addleaf(Knode* node)
{
	if (node->isleaf()) {
	    *wkgr++ = node->tid;
	    if (wkwt) *wkwt++ = (FTYPE) node->vol;
	} else {
	    addleaf(node->left);
	    addleaf(node->right);
	}
}

void Sptree::addleaf_ss(Knode* node)
{
	if (node->isleaf()) {
	    *wksl++ = slst[node->tid];
	    if (wkwt) *wkwt++ = (FTYPE) node->vol;
	} else {
	    addleaf_ss(node->left);
	    addleaf_ss(node->right);
	}
}

void Sptree::collectleaf(mSeq* sub, mSeq* sd, Knode* node)
{
#if USE_WEIGHT
        FTYPE*  weight  = new FTYPE[node->ndesc + 1];
        wkwt = weight;
        wkgr = group;
        addleaf(node);
        *wkwt = *wkgr = -1;
	for (wkwt = weight; *wkwt >= 0; )
	    *wkwt++ *= node->cur / node->vol;
        sub = sd->extseq(sub, group, CPY_SEQ);
        delete[] sub->weight;
        sub->weight = weight;
#else
	wkwt = 0;
        wkgr = group;
        addleaf(node);
        *wkgr = -1;
        sub = sd->extseq(sub, group, CPY_SEQ);
#endif	// USE_WEIGHT
}
 
int Sptree::collect_ss(mSeq* sub, mSeq** slct, Knode* node)
{
	wkwt = slwt;
	wksl = slct;
	addleaf_ss(node);
	*wksl = 0;
#if USE_WEIGHT
	int	mm = 0;
	for (wksl = slct, wkwt = slwt; *wksl; ) {
	    mm += (*wksl++)->many;
	    *wkwt++ *= node->cur / node->vol;
	}
	delete[] sub->weight;
	sub->weight = new FTYPE[mm];
#endif
	int	nn = wksl - slct;
	aggregate_sb(sub, slct, 0, slwt);
	return (nn);
}

VTYPE Sptree::sptree(mSeq* sprf, Knode* node, bool use_pw)
{
	VTYPE	scr = 0;
	int	grp[2] = {0, -1};

	if (node->isleaf()) {
	    grp[0] = node->tid;
	    sprf = cur_sd->extseq(sprf, grp, CPY_SEQ);
#if USE_WEIGHT
	    if (!sprf->weight) sprf->weight = new FTYPE[1];
	    sprf->weight[0] = node->vol / node->parent->vol;
#endif
	} else if (node->ndesc <= ndesc_thr) {
	    collectleaf(sprf, cur_sd, node);
#if USE_WEIGHT
	    sprf->pairwt = use_pw? ktree->recalcpw(node): 0;
	    Msap	msd(sprf, NTV_ALN);
	    scr = msd.ps_nml();
	    delete[] sprf->pairwt;
	    sprf->pairwt = 0;
#else
	    Msap	msd(sprf, NTV_ALN);
	    scr = msd.ps_nml();
#endif
	} else {
	    mSeq*	sqs[3];
	    initseq(sqs, 2); sqs[2] = 0;
	    sprf->copyattr(sqs[0]);
	    sprf->copyattr(sqs[1]);
	    scr  = sptree(sqs[0], node->left, use_pw);
	    scr += sptree(sqs[1], node->right, use_pw);
	    sqs[0]->convseq(VECTOR);
	    sqs[1]->convseq(VECTOR);
	    sprf = fuseseq(sprf, sqs, (FTYPE) node->cur);
	    scr += calcscore_grp(sqs);
	    clearseq(sqs, 2);
	}
	return (scr);
}

VTYPE Sptree::sptree_ss(mSeq* sprf, Knode* node, bool use_pw)
{
	VTYPE	scr = 0;

	if (node->isleaf()) {
	    sprf = slst[node->tid]->copyseq(sprf, CPY_SEQ);
#if USE_WEIGHT
	    if (!sprf->weight) {
		sprf->weight = new FTYPE[1];
		sprf->weight[0] = node->vol / node->parent->vol;
	    } else
		for (int i = 0; i < sprf->many; ++i)
		    sprf->weight[i] *= node->vol / node->parent->vol;
#endif
	} else if (node->ndesc <= ndesc_thr) {
#if USE_WEIGHT
	    int	i = collect_ss(sprf, slct, node);
	    FTYPE*	pw = ktree->recalcpw(node, i);
	    scr = ps_lst(slct, (FTYPE*) pw);
	    delete[] pw;
#else
	    collect_ss(sprf, slct, node);
	    scr = ps_lst(slct);
#endif
	} else {
	    mSeq*	sqs[3];
	    initseq(sqs, 2); sqs[2] = 0;
	    sprf->copyattr(sqs[0]); sprf->copyattr(sqs[1]);
	    scr  = sptree_ss(sqs[0], node->left, use_pw);
	    scr += sptree_ss(sqs[1], node->right, use_pw);
	    sqs[0]->convseq(VECTOR);
	    sqs[1]->convseq(VECTOR);
	    sprf = fuseseq(sprf, sqs, (FTYPE) node->cur);
	    scr += calcscore_grp(sqs);
	    clearseq(sqs, 2);
	}
	return (scr);
}

Sptree::Sptree(mSeq* sd, Ssrel* srl) : cur_sd(sd), ktree(srl->ktree)
{
	if (srl->ss->num < srl->ss->elms) {
	    slst = srl->takeapart(0, sd, true);
	    slct = new mSeq*[srl->ss->num + 1];
#if USE_WEIGHT
	    slwt = new FTYPE[srl->ss->num + 1];
#else
	    slwt = 0;
#endif
	    group = 0;
	} else {
	    slst = slct = 0; slwt = 0;
	    group = new int[srl->ss->num + 1];
	}
}

Sptree::Sptree(mSeq* sd, int num, Ktree* kt) : cur_sd(sd), ktree(kt)
{
	slst = slct = 0; slwt = 0;
	group = new int[num + 1];
}

Sptree::~Sptree()
{
	delete[] group;
	if (slst) for (mSeq** sq = slst; *sq; ++sq) delete *sq;
	delete[] slst; delete[] slct; delete[] slwt;
}

/******************************************************************
*	pairsum --	Calculate sum-of-pairs measure with weight
******************************************************************/

VTYPE Ssrel::pairsum_ss(mSeq* sd, bool use_pw)
{
	if (sd->many < 2) return (0);
	VTYPE	ps;

	if (ss) {
	    Sptree	spt(sd, this);
	    mSeq	tmp;
	    sd->copyattr(&tmp);
	    if (ss->num < ss->elms)
		ps = spt.sptree_ss(&tmp, ktree->root, use_pw);
	    else
		ps = spt.sptree(&tmp, ktree->root, use_pw);
	} else {
#if USE_WEIGHT
	    swap(sd->pairwt, pairwt);
	    ps = pairsum(sd);
	    swap(sd->pairwt, pairwt);
#else
	    ps = pairsum(sd);
#endif
	}
	if (alprm.u0 != 0)
	    ps -= (VTYPE) (alprm.u0 * sd->countunps());
	if (sd->sigII) ps += spSigII(sd);
	return ((VTYPE) ps);
}

#if TST_PS_ALG

#include <sys/types.h>
#include <sys/times.h>

/*   pairsum -- to test performance of each method */

VTYPE test_pairsum(mSeq* sd, int usg)
{
	if (sd->many == 1) return (0);

	mSeq*	temp = 0;
	Ktree*	ktree = 0;
	VTYPE	ps = 0;
struct	tms t0, t1, t2;
	long	cput1;
	long	cput2;

	times(&t0);
	times(&t1);
	Msap*	msd = 0;
	switch (usg) {
	    case 0:
		msd = new Msap(sd, 1);
		ps = msd->ps_nml(sd);
		break;
	    case 1:
	    case 2:
	    case 3:
		msd = new Msap(sd, 2);
		ps = msd->ps_gfq(sd);
		break;
	    case 4:
		ktree = new Ktree(sd);
#if USE_WEIGHT
		sd->pairwt = ktree->calcpw();
		times(&t1);
		msd = new Msap(sd, 1);
		ps = msd->ps_nml(sd);
		delete[] sd->pairwt;
		sd->pairwt = 0;
#else
		msd = new Msap(sd, 1);
		ps = msd->ps_nml(sd);
#endif
		break;
	    case 5:
	    case 6:
		ktree = new Ktree(sd);
		if (usg & 1) ktree->root->findcenter();
		ktree->root->ros = fInfinit;
#if USE_WEIGHT
		sd->pairwt = ktree->recalcpw();
		times(&t1);
#endif
		Sptree	spt(sd, sd->many, ktree);
		ps = spt.sptree(temp, ktree->root);
		break;
	}
	delete msd;
	times(&t2);
	cput1 = t1.tms_utime + t1.tms_stime - t0.tms_utime - t0.tms_stime;
	cput2 = t2.tms_utime + t2.tms_stime - t1.tms_utime - t1.tms_stime;
	printf("%6ld %6ld %2d", cput1, cput2, ndesc_thr); 
	return (ps);
}

#endif


