/*****************************************************************************
*
*	dpunit.h
*
*	Dada structure used in fwd2c and fwd2d
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-4-7 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-2012)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#ifndef	_DPUNIT_H_
#define	_DPUNIT_H_

#include "aln.h"
#include "mseq.h"
#include "gfreq.h"

struct DPunit {
	VTYPE   val;
	int     dir;
	long    ptr;
};

struct DPunit_nv : public DPunit {
	int     dtsize;
	int*    gla;
	int*    glb;
};

struct DPunit_hf : public DPunit {
	IDELTA* dla;
	int     glb;
};

struct DPunit_pf : public DPunit {
	IDELTA* dla;
	IDELTA* dlb;
};

struct SwgDPunit {
	VTYPE   val;
	int     dir;
	int     lwr;
	int     upr;
	int     mlb;
	int     nlb;
	int     mrb;
	int     nrb;
	COLONY* clny;
};

struct SwgDPunit_nv : public SwgDPunit {
	int     dtsize;
	int*    gla;
	int*    glb;
};

struct SwgDPunit_hf : public SwgDPunit {
	IDELTA* dla;
	int     glb;
};

struct SwgDPunit_pf : public SwgDPunit {
	IDELTA* dla;
	IDELTA* dlb;
};

static const DPunit blank_dpunit = {0, 0, 0};
static const DPunit black_dpunit = {NEVSEL, 0, 0};
static const SwgDPunit blank_swgdpunit = {0, 0, POS_INT, NEG_INT};
static const SwgDPunit black_swgdpunit = {NEVSEL, 0, POS_INT, NEG_INT};

template <class recd_t>
void    clear(recd_t* rcd);

template <class recd_t>
void    reset(recd_t* rcd);

template <class recd_t>
void    copy(recd_t* dst, const recd_t* src);

#endif	// _DPUNIT_H_

