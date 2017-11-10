/*****************************************************************************
*
*	fwd2c.c:
*
*	Alignment of two protein or nucleotide sequence
*	 or groups of sequences with Algorithm <<C>>
*
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

#include "fwd2c.h"

// Non-internal gap algorithm

template <>
void Fwd2c<DPunit>::initializeC()
{
	an = a->many;
	bn = b->many;
	glab = 0;
	gfab = 0;
	vset(buf, black_dpunit, bufsize);
}

/*******************************************************
	0	1	2	3	4	5
	**	-*	*-	--	*?	-?
0  **	0	0	1	0
1  -*	0	0	1	g
2  *-	1	1	0	0
3  --	0	g	0	0
*******************************************************/

template<>
VTYPE Fwd2c<DPunit>::gapopen(const DPunit* rcd, const mSeqItr& asi, const mSeqItr& bsi, const int d3)
{
	VTYPE	axb = 0;
	if (d3 > 0) {
	    if (!isvert(rcd)) {		// \| || -|	//
		axb = (asi.didns && bsi.didns)? 
		    asi.dns->cfq * bsi.didns->cx + asi.didns->dc * bsi.didns->dx * alprm.gamma :
		    asi.dns->cfq * bsi.dns->efq;
	    } else if (asi.didns) {	// ||		//
		axb = asi.didns->dc * bsi.dns->efq * alprm.gamma;
	    }
	} else if (d3 < 0) {
	    if (!ishori(rcd)) {		// \- || |-	//
		axb = (asi.didns && bsi.didns)? 
		    bsi.dns->cfq * asi.didns->cx + bsi.didns->dc * asi.didns->dx * alprm.gamma :
		    bsi.dns->cfq * asi.dns->efq;
	    } else if (bsi.didns) {	// --		//
		axb = bsi.didns->dc * asi.dns->efq * alprm.gamma;
	    }
	} else if (!asi.didns && !bsi.didns) {
	    return (0);
	} else if (isdiag(rcd)) {	// \\	//
	  if (asi.didns) {
	    if (bsi.didns) {
		axb = asi.didns->cd * (bsi.didns->cc + bsi.didns->dc)
		    + bsi.didns->cd * (asi.didns->cc + asi.didns->dc)
		    + (asi.didns->dc * bsi.didns->dd + asi.didns->dd * bsi.didns->dc) * alprm.gamma;
	    } else {
		axb = asi.didns->cd * bsi.dns->cfq;
	    }
	  } else {
	    axb = bsi.didns->cd * asi.dns->cfq;
	  }
	} else if (isvert(rcd) && asi.didns) {	// |\ 	//
	    axb = bsi.dns->cfq * (asi.didns->cd + asi.didns->dd * alprm.gamma);
	} else if (ishori(rcd) && bsi.didns) {	// -\	// 
	    axb = asi.dns->cfq * (bsi.didns->cd + bsi.didns->dd * alprm.gamma);
	}
	return (pwd->vgop(axb));
}

template<>
void Fwd2c<DPunit>::update(DPunit* dst, const DPunit* src, 
	const mSeqItr& asi, const mSeqItr& bsi, const VTYPE gpn, const int d3)
{
	if (d3 > 0) 	dst->dir = ishori(src)? NEWV: VERT; else
	if (d3 < 0)	dst->dir = isvert(src)? NEWH: HORI;
	else		dst->dir = isdiag(src)? DIAG: NEWD;
	dst->val = src->val + gpn;
	dst->ptr = src->ptr;
}

// Naive algorithm

template<>
VTYPE Fwd2c<DPunit_nv>::gapopen(const DPunit_nv* rcd, 
	const mSeqItr& asi, const mSeqItr& bsi, const int d3)
{
	return (pwd->*pwd->crg2)(rcd->gla, rcd->glb, asi, bsi, d3);
}

template<>
void Fwd2c<DPunit_nv>::update(DPunit_nv* dst, const DPunit_nv* src, 
	const mSeqItr& asi, const mSeqItr& bsi, const VTYPE gpn, const int d3)
{
	if (d3 == 0)
	    dst->dir = isdiag(src)? DIAG: NEWD;
	else if (d3 > 0)
	    dst->dir = ishori(src)? NEWV: VERT;
	else
	    dst->dir = isvert(src)? NEWH: HORI;

	elongap(dst->gla, src->gla, d3 >= 0? asi.res: 0, an);
	elongap(dst->glb, src->glb, d3 <= 0? bsi.res: 0, bn);
	dst->val = src->val + gpn;
	dst->ptr = src->ptr;
}

template<>
void Fwd2c<DPunit_nv>::initializeC()
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
VTYPE Fwd2c<DPunit_hf>::gapopen(const DPunit_hf* rcd, const mSeqItr& asi, const mSeqItr& bsi, const int d3)
{
	if (d3 == 0)
	    return pwd->newgap2(*asi.tfq, rcd->glb, rcd->dla);
	else if (d3 > 0)
	    return pwd->newgap1(*asi.sfq, rcd->dla, rcd->glb);
	else
	    return pwd->newgap2(*asi.rfq, rcd->glb, rcd->dla);
}

template<>
void Fwd2c<DPunit_hf>::update(DPunit_hf* dst, const DPunit_hf* src, 
	const mSeqItr& asi, const mSeqItr& bsi, const VTYPE gpn, const int d3)
{
	if (d3 == 0) {
	    dst->dir = isdiag(src)? DIAG: NEWD;
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    dst->glb = 0;
	} else if (d3 > 0) {
	    dst->dir = ishori(src)? NEWV: VERT;
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    dst->glb = src->glb + 1;
	} else {
	    dst->dir = isvert(src)? NEWH: HORI;
	    incdelta(dst->dla, src->dla);
	    dst->glb = 0;
	}
	dst->val = src->val + gpn;
	dst->ptr = src->ptr;
}

template<>
void Fwd2c<DPunit_hf>::initializeC()
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
VTYPE Fwd2c<DPunit_pf>::gapopen(const DPunit_pf* rcd, const mSeqItr& asi, const mSeqItr& bsi, const int d3)
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
void Fwd2c<DPunit_pf>::update(DPunit_pf* dst, const DPunit_pf* src, 
	const mSeqItr& asi, const mSeqItr& bsi, const VTYPE gpn, const int d3)
{
	if (d3 == 0) {		// diagonal
	    dst->dir = isdiag(src)? DIAG: NEWD;
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    newdelta(dst->dlb, *bsi.tfq, src->dlb);
	} else if (d3 > 0) {	// vertical
	    dst->dir = ishori(src)? NEWV: VERT;
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    incdelta(dst->dlb, src->dlb);
	} else {		// horizontal
	    dst->dir = isvert(src)? NEWH: HORI;
	    newdelta(dst->dlb, *bsi.tfq, src->dlb);
	    incdelta(dst->dla, src->dla);
	}
	dst->val = src->val + gpn;
	dst->ptr = src->ptr;
}

template<>
void Fwd2c<DPunit_pf>::initializeC()
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

// Swg local alignment

// Non-internal gap algorithm

template <>
void Fwd2c<SwgDPunit>::initializeC()
{
	an = a->many;
	bn = b->many;
	glab = 0;
	gfab = 0;
	vset(buf, black_swgdpunit, bufsize);
}

template<>
VTYPE Fwd2c<SwgDPunit>::gapopen(const SwgDPunit* rcd, const mSeqItr& asi, const mSeqItr& bsi, const int d3)
{
	return (isdiag(rcd) && d3)? pwd->BasicGOP: 0;
}

void swg_gdpunit_update(SwgDPunit* dst, const SwgDPunit* src,
	const mSeqItr& asi, const mSeqItr& bsi, const VTYPE gpn, const int d3)
{
	*dst = *src;
	dst->val += gpn;
	int     r = bsi.pos - asi.pos;
	if (d3 > 0) {
	    dst->dir = VERT;
	    if (r < dst->lwr) dst->lwr = r;
	} else if (d3 < 0) {
	    dst->dir = HORI;
	    if (r > dst->upr) dst->upr = r;
	} else  {
	    dst->dir = isdiag(src)? DIAG: NEWD;
	}
}

template<>
void Fwd2c<SwgDPunit>::update(SwgDPunit* dst, const SwgDPunit* src, 
	const mSeqItr& asi, const mSeqItr& bsi, const VTYPE gpn, const int d3)
{
	swg_gdpunit_update(dst, src, asi, bsi, gpn, d3);
}

// Naive algorithm

template<>
VTYPE Fwd2c<SwgDPunit_nv>::gapopen(const SwgDPunit_nv* rcd, 
	const mSeqItr& asi, const mSeqItr& bsi, const int d3)
{
	return (pwd->*pwd->crg2)(rcd->gla, rcd->glb, asi, bsi, d3);
}

template<>
void Fwd2c<SwgDPunit_nv>::update(SwgDPunit_nv* dst, const SwgDPunit_nv* src, 
	const mSeqItr& asi, const mSeqItr& bsi, const VTYPE gpn, const int d3)
{
	swg_gdpunit_update((SwgDPunit*) dst, (const SwgDPunit*) src, asi, bsi, gpn, d3);
	elongap(dst->gla, src->gla, d3 >= 0? asi.res: 0, an);
	elongap(dst->glb, src->glb, d3 <= 0? bsi.res: 0, bn);
}

template<>
void Fwd2c<SwgDPunit_nv>::initializeC()
{
	an = a->many;
	bn = b->many;
	int	apb = an + bn;
	glab = new int[apb * bufsize];
	vclear(glab, apb * bufsize);
	SwgDPunit_nv* r = buf;
	int*	w = glab - bn;
	for (INT n = 0; n < bufsize; ++n, ++r) {
	    r->dtsize = apb;
	    r->gla = w += bn;
	    r->glb = w += an;
	}
	gfab = 0;
}

// Half-general profile

template<>
VTYPE Fwd2c<SwgDPunit_hf>::gapopen(const SwgDPunit_hf* rcd, const mSeqItr& asi, const mSeqItr& bsi, const int d3)
{
	if (d3 == 0)
	    return pwd->newgap2(*asi.tfq, rcd->glb, rcd->dla);
	else if (d3 > 0)
	    return pwd->newgap1(*asi.sfq, rcd->dla, rcd->glb);
	else
	    return pwd->newgap2(*asi.rfq, rcd->glb, rcd->dla);
}

template<>
void Fwd2c<SwgDPunit_hf>::update(SwgDPunit_hf* dst, const SwgDPunit_hf* src, 
	const mSeqItr& asi, const mSeqItr& bsi, const VTYPE gpn, const int d3)
{
	swg_gdpunit_update((SwgDPunit*) dst, (const SwgDPunit*) src, asi, bsi, gpn, d3);
	if (d3 == 0) {
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    dst->glb = 0;
	} else if (d3 > 0) {
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    dst->glb = src->glb + 1;
	} else {
	    incdelta(dst->dla, src->dla);
	    dst->glb = 0;
	}
}

template<>
void Fwd2c<SwgDPunit_hf>::initializeC()
{
	an = a->gfq->hetero;
	bn = b->many;
	glab = 0;
	int	apb = an + 1;
	gfab = new IDELTA[apb * bufsize];
	IDELTA*	w = gfab - apb;
	SwgDPunit_hf* r = buf;
	for (INT n = 0; n < bufsize; ++n, ++r) {
	    cleardelta(r->dla = w += apb);
	    r->glb = 0;
	}
}

// Full-general profile

template<>
VTYPE Fwd2c<SwgDPunit_pf>::gapopen(const SwgDPunit_pf* rcd, const mSeqItr& asi, const mSeqItr& bsi, const int d3)
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
void Fwd2c<SwgDPunit_pf>::update(SwgDPunit_pf* dst, const SwgDPunit_pf* src, 
	const mSeqItr& asi, const mSeqItr& bsi, const VTYPE gpn, const int d3)
{
	swg_gdpunit_update((SwgDPunit*) dst, (const SwgDPunit*) src, asi, bsi, gpn, d3);
	if (d3 == 0) {		// diagonal
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    newdelta(dst->dlb, *bsi.tfq, src->dlb);
	} else if (d3 > 0) {	// vertical
	    newdelta(dst->dla, *asi.tfq, src->dla);
	    incdelta(dst->dlb, src->dlb);
	} else {		// horizontal
	    newdelta(dst->dlb, *bsi.tfq, src->dlb);
	    incdelta(dst->dla, src->dla);
	}
}

template<>
void Fwd2c<SwgDPunit_pf>::initializeC()
{
	an = a->gfq->hetero;
	bn = b->gfq->hetero;
	glab =  0;
	int	apb = an + bn + 2;
	gfab = new IDELTA[apb * bufsize];
	IDELTA*	w = gfab;
	SwgDPunit_pf* r = buf;
	w = gfab - bn - 1;
	for (INT n = 0; n < bufsize; ++n, ++r) {
	    cleardelta(r->dla = w += bn + 1);
	    cleardelta(r->dlb = w += an + 1);
	}
}

