/*****************************************************************************
*
*	Alignment of two protein or nucleotide sequences.
*	Splicing is NOT considered.
*	Assumes there is no internal gap in either sequence.
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-2012)
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

#include "fwd2b.h"

// Non-internal gap algorithm

template <>
void Fwd2b<DPunit>::initializeB()
{
	an = a->many;
	bn = b->many;
	glab = 0;
	gfab = 0;
	vset(buf, black_dpunit, bufsize);
}

template <>
void clear(DPunit* rcd)
{
	*rcd = blank_dpunit;
}

template <>
void reset(DPunit* rcd)
{
	*rcd = black_dpunit;
}
template <>
void copy(DPunit* dst, const DPunit* src)
{
	*dst = *src;
}

template<>
VTYPE Fwd2b<DPunit>::gapopen(const DPunit* rcd, const mSeqItr& asi, const mSeqItr& bsi, const int d3)
{
	return ((!isvert(rcd) && d3 > 0) || (!ishori(rcd) && d3 < 0))? pwd->BasicGOP: 0;
}

template<>
void Fwd2b<DPunit>::update(DPunit* dst, const DPunit* src, 
	const mSeqItr& asi, const mSeqItr& bsi, const VTYPE gpn, const int d3)
{
	if (d3 > 0) 	dst->dir = VERT; else
	if (d3 < 0)	dst->dir = HORI;
	else		dst->dir = isdiag(src)? DIAG: NEWD;
	dst->val = src->val + gpn;
	dst->ptr = src->ptr;
}

// Naive algorithm

template <>
void clear(DPunit_nv* rcd)
{
	clear((DPunit*) rcd);
	vclear(rcd->gla, rcd->dtsize);
}

template <>
void reset(DPunit_nv* rcd)
{
	reset((DPunit*) rcd);
	vclear(rcd->gla, rcd->dtsize);
}

template <>
void copy(DPunit_nv* dst, const DPunit_nv* src)
{
	copy((DPunit*) dst, (const DPunit*) src);
	vcopy(dst->gla, src->gla, src->dtsize);
}

template<>
VTYPE Fwd2b<DPunit_nv>::gapopen(const DPunit_nv* rcd, 
	const mSeqItr& asi, const mSeqItr& bsi, const int d3)
{
	return (pwd->*pwd->crg2)(rcd->gla, rcd->glb, asi, bsi, d3);
}

template<>
void Fwd2b<DPunit_nv>::update(DPunit_nv* dst, const DPunit_nv* src, 
	const mSeqItr& asi, const mSeqItr& bsi, const VTYPE gpn, const int d3)
{
	if (d3 == 0)
	    dst->dir = isdiag(src)? DIAG: NEWD;
	else if (d3 > 0)
	    dst->dir = VERT;
	else
	    dst->dir = HORI;

	elongap(dst->gla, src->gla, d3 >= 0? asi.res: 0, an);
	elongap(dst->glb, src->glb, d3 <= 0? bsi.res: 0, bn);
	dst->val = src->val + gpn;
	dst->ptr = src->ptr;
}

template<>
void Fwd2b<DPunit_nv>::initializeB()
{
	an = a->many;
	bn = b->many;
	int	apb = an + bn;
	glab = new int[apb * bufsize];
	vclear(glab, apb * bufsize);
	DPunit_nv* r = buf;
	int*	w = glab - bn;
	for (INT n = 0; n < bufsize; ++n, ++r) {
	    r->val = NEVSEL; r->ptr = r->dir = 0;
	    r->dtsize = apb;
	    r->gla = w += bn;
	    r->glb = w += an;
	}
	gfab = 0;
}

// Half-general profile

template<>
void clear(DPunit_hf* rcd)
{
	clear((DPunit*) rcd);
	rcd->glb = 0;
	cleardelta(rcd->dla);
}

template<>
void reset(DPunit_hf* rcd)
{
	reset((DPunit*) rcd);
	rcd->glb = 0;
	cleardelta(rcd->dla);
}

template<>
void copy(DPunit_hf* dst, const DPunit_hf* src)
{
	copy((DPunit*) dst, (const DPunit*) src);
	copydelta(dst->dla, src->dla);
	dst->glb = src->glb;
}

template<>
VTYPE Fwd2b<DPunit_hf>::gapopen(const DPunit_hf* rcd, const mSeqItr& asi, const mSeqItr& bsi, const int d3)
{
	if (d3 == 0)
	    return pwd->newgap2(*asi.tfq, rcd->glb, rcd->dla);
	else if (d3 > 0)
	    return pwd->newgap1(*asi.sfq, rcd->dla, rcd->glb);
	else
	    return pwd->newgap2(*asi.rfq, rcd->glb, rcd->dla);
}

template<>
void Fwd2b<DPunit_hf>::update(DPunit_hf* dst, const DPunit_hf* src, 
	const mSeqItr& asi, const mSeqItr& bsi, const VTYPE gpn, const int d3)
{
	if (d3 == 0) {
	    dst->dir = isdiag(src)? DIAG: NEWD;
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    dst->glb = 0;
	} else if (d3 > 0) {
	    dst->dir = VERT;
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    dst->glb = src->glb + 1;
	} else {
	    dst->dir = HORI;
	    incdelta(dst->dla, src->dla);
	    dst->glb = 0;
	}
	dst->val = src->val + gpn;
	dst->ptr = src->ptr;
}

template<>
void Fwd2b<DPunit_hf>::initializeB()
{
	an = a->gfq->hetero;
	bn = b->many;
	glab = 0;
	int	apb = an + 1;
	gfab = new IDELTA[apb * bufsize];
	IDELTA*	w = gfab - apb;
	DPunit_hf* r = buf;
	for (INT n = 0; n < bufsize; ++n, ++r) {
	    r->val = NEVSEL; r->ptr = r->dir = 0;
	    cleardelta(r->dla = w += apb);
	    r->glb = 0;
	}
}

// Full-general profile

template<>
void clear(DPunit_pf* rcd)
{
	clear((DPunit*) rcd);
	cleardelta(rcd->dla);
	cleardelta(rcd->dlb);
}

template<>
void reset(DPunit_pf* rcd)
{
	reset((DPunit*) rcd);
	cleardelta(rcd->dla);
	cleardelta(rcd->dlb);
}

template<>
void copy(DPunit_pf* dst, const DPunit_pf* src)
{
	copy((DPunit*) dst, (const DPunit*) src);
	copydelta(dst->dla, src->dla);
	copydelta(dst->dlb, src->dlb);
}

template<>
VTYPE Fwd2b<DPunit_pf>::gapopen(const DPunit_pf* rcd, const mSeqItr& asi, const mSeqItr& bsi, const int d3)
{
	if (d3 == 0)		// diagonal
	    return pwd->newgap3(*asi.sfq, rcd->dla, *bsi.tfq, rcd->dlb)
		 + pwd->newgap3(*bsi.sfq, rcd->dlb, *asi.tfq, rcd->dla);
	else if (d3 > 0)	// vertical
	    return pwd->newgap3(*asi.sfq, rcd->dla, *bsi.rfq, rcd->dlb);
	else 			// horizontal
	    return pwd->newgap3(*bsi.sfq, rcd->dlb, *asi.rfq, rcd->dla);
}

template<>
void Fwd2b<DPunit_pf>::update(DPunit_pf* dst, const DPunit_pf* src, 
	const mSeqItr& asi, const mSeqItr& bsi, const VTYPE gpn, const int d3)
{
	if (d3 == 0) {		// diagonal
	    dst->dir = isdiag(src)? DIAG: NEWD;
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    newdelta(dst->dlb, *bsi.tfq, src->dlb);
	} else if (d3 > 0) {	// vertical
	    dst->dir = VERT;
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    incdelta(dst->dlb, src->dlb);
	} else {		// horizontal
	    dst->dir = HORI;
	    newdelta(dst->dlb, *bsi.tfq, src->dlb);
	    incdelta(dst->dla, src->dla);
	}
	dst->val = src->val + gpn;
	dst->ptr = src->ptr;
}

template<>
void Fwd2b<DPunit_pf>::initializeB()
{
	an = a->gfq->hetero;
	bn = b->gfq->hetero;
	glab =  0;
	int	apb = an + bn + 2;
	gfab = new IDELTA[apb * bufsize];
	IDELTA*	w = gfab;
	DPunit_pf* r = buf;
	w = gfab - bn - 1;
	for (INT n = 0; n < bufsize; ++n, ++r) {
	    r->val = NEVSEL; r->ptr = r->dir = 0;
	    cleardelta(r->dla = w += bn + 1);
	    cleardelta(r->dlb = w += an + 1);
	}
}
