/*****************************************************************************
*
*	Alignment of two protein/nucleotide sequences 
*	5' and 3' splice site signals and coding potential
*	calculated from ditron frequnecies are considered.
*	Assumes there are internal gap(s) in the reference protein sequence.
*	The maximum number of candidates is limited to three.
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

#include "fwd2h.h"

//	RVPDJ_nv

template<>
void Fwd2h<RVPDJ_nv>::initializeH(int recsize)
{
	RVPDJ_nv*	r = recbuf;
	RVPDJ_nv*	t = r + recsize;
	for ( ; r < t; ++r) r->reseth();
}

template<>
VTYPE Fwd2h<RVPDJ_nv>::gapopen(const RVPDJ_nv* rcd, const mSeqItr& asi, int d3)
{
	return (((rcd->gla >= rcd->glb && d3 > 0) ||
		 (rcd->gla <= rcd->glb && d3 < 0))? pwd->BasicGOP: 0);
}

template<>
void Fwd2h<RVPDJ_nv>::update(RVPDJ_nv* dst, const RVPDJ_nv* src, const mSeqItr& asi, 
	VTYPE gop, int d3)
{
	if (d3 == 0) {
	    dst->gla = dst->glb = 0;
	} else if (d3 > 0) {
	    dst->gla = 0;
	    dst->glb = src->glb + d3;
	} else {
	    dst->gla = src->gla - d3;
	    dst->glb = 0;
	}
	dst->ptr = src->ptr;
	dst->jnc = src->jnc;
	dst->val = src->val + gop;
}

template <>
void Fwd2h<RVPDJ_nv>::update(RVPDJ_nv* dst, const mSeqItr& asi, int d3)
{
	if (d3 == 0) {
	    dst->gla = dst->glb = 0;
	} else if (d3 > 0) {
	    dst->gla = 0;
	    dst->glb += d3;
	} else {
	    dst->gla -= d3;
	    dst->glb = 0;
	}
}

//	RVPDJ_hf

template<>
void Fwd2h<RVPDJ_hf>::initializeH(int recsize)
{
	int	apb = a->gfq->hetero + 1;
	gfab = new IDELTA[apb * recsize];
	IDELTA*	w = gfab - apb;
	RVPDJ_hf* r = recbuf;
	for (int n = 0; n < recsize; ++n, ++r) {
	    r->val = NEVSEL; r->ptr = r->dir = 0;
	    cleardelta(r->dla = w += apb);
	    r->glb = 0;
	}
}

template<>
VTYPE Fwd2h<RVPDJ_hf>::gapopen(const RVPDJ_hf* rcd, const mSeqItr& asi, int d3)
{
	if (d3 == 0)
	    return pwd->newgap2(*asi.tfq, rcd->glb, rcd->dla);
	else if (d3 > 0)
	    return pwd->newgap1(*asi.sfq, rcd->dla, rcd->glb);
	else
	    return pwd->newgap2(*asi.rfq, rcd->glb, rcd->dla);
}

template<>
void Fwd2h<RVPDJ_hf>::update(RVPDJ_hf* dst, const RVPDJ_hf* src,
	const mSeqItr& asi, VTYPE gop, int d3)
{
	if (d3 == 0) {
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    dst->glb = 0;
	} else if (d3 > 0) {
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    dst->glb = src->glb + d3;
	} else {
	    incdelta(dst->dla, src->dla, -d3);
	    dst->glb = 0;
	}
	dst->ptr = src->ptr;
	dst->jnc = src->jnc;
	dst->val = src->val + gop;
}

template<>
void Fwd2h<RVPDJ_hf>::update(RVPDJ_hf* rcd, const mSeqItr& asi, int d3)
{
	if (d3 == 0) {
	    newdelta(rcd->dla, *asi.tfq, rcd->dla);
	    rcd->glb = 0;
	} else if (d3 > 0) {
	    newdelta(rcd->dla, *asi.tfq, rcd->dla);
	    rcd->glb = rcd->glb + d3;
	} else {
	    incdelta(rcd->dla, rcd->dla, -d3);
	    rcd->glb = 0;
	}
}

