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
*	Copyright(c) Osamu Gotoh <<o.gotoh@i.kyoto-u.ac.jp>>
*****************************************************************************/

#include "fwd2h.h"

//	DPUH

template <>
void reseth(DPUH* rcd, bool all) {if (all) *rcd = blackDPUH;}

template <>
void copyh(DPUH* dst, DPUH* src, bool all) {if (all) *dst = *src;}

template<>
void Fwd2h<DPUH>::initializeH(int recsize)
{
	DPUH*	r = recbuf;
	for (int n = 0; n < recsize; ++n, ++r) *r = blackDPUH;
}

template<>
VTYPE Fwd2h<DPUH>::gapopen(DPUH* rcd, mSeqItr& asi, int d3)
{
	return (d3 == 0? 0: pwd->BasicGOP);
}

template<>
VTYPE Fwd2h<DPUH>::ins_penalty(DPUH* rcd, mSeqItr& asi, int d3)
{
	return pwd->UnpPenalty3(d3);
}

template<>
VTYPE Fwd2h<DPUH>::del_penalty(DPUH* rcd, mSeqItr& asi, int d3, FSTAT* fst)
{
	fst->unp += pwd->Vab * d3;
	fst->gap += gapopen(rcd, asi, d3);
	asi += d3 / 3;
	return pwd->UnpPenalty3(d3);
}

template<>
void Fwd2h<DPUH>::update(DPUH* dst, DPUH* src, mSeqItr& asi, int d3)
{
	dst->ptr = src->ptr;
	dst->jnc = src->jnc;
}

template <>
void Fwd2h<DPUH>::update(DPUH* rcd, mSeqItr& asi, int d3)
{
}

//	DPUH_hf

template <>
void reseth(DPUH_hf* rcd, bool all)
{
	if (all) reseth((DPUH*) rcd);
	rcd->glb = 0;
	cleardelta(rcd->dla);
}

template <>
void copyh(DPUH_hf* dst, DPUH_hf* src, bool all)
{
	if (all) copyh((DPUH*) dst, (DPUH*) src);
	copydelta(dst->dla, src->dla);
	dst->glb = src->glb;
}

template<>
void Fwd2h<DPUH_hf>::initializeH(int recsize)
{
	int	apb = a->gfq->hetero + 1;
	gfab = new IDELTA[apb * recsize];
	IDELTA*	w = gfab - apb;
	DPUH_hf* r = recbuf;
	for (int n = 0; n < recsize; ++n, ++r) {
	    r->val = NEVSEL; r->ptr = r->dir = 0;
	    cleardelta(r->dla = w += apb);
	    r->glb = 0;
	}
}

template<>
VTYPE Fwd2h<DPUH_hf>::gapopen(DPUH_hf* rcd, mSeqItr& asi, int d3)
{
	if (d3 == 0)
	    return pwd->newgap2(*asi.tfq, rcd->glb, rcd->dla);
	else if (d3 > 0)
	    return pwd->newgap1(*asi.sfq, rcd->dla, rcd->glb);
	else
	    return pwd->newgap2(*asi.rfq, rcd->glb, rcd->dla);
}

template<>
VTYPE Fwd2h<DPUH_hf>::ins_penalty(DPUH_hf* rcd, mSeqItr& asi, int d3)
{
	VTYPE	uscr = 0;
	for (GFREQ* cf = *asi.rfq; neogfq(cf); ++cf) {
	    int	i = GapLenSD(cf, rcd->dla);
	    uscr += cf->freq * pwd->UnpPenalty(i);
	}
	return uscr;
}

template<>
void Fwd2h<DPUH_hf>::update(DPUH_hf* dst, DPUH_hf* src,
	mSeqItr& asi, int d3)
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
}

template<>
void Fwd2h<DPUH_hf>::update(DPUH_hf* rcd, mSeqItr& asi, int d3)
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

template<>
VTYPE Fwd2h<DPUH_hf>::del_penalty(DPUH_hf* rcd, mSeqItr& asi, int d3, FSTAT* fst)
{
	VTYPE	uscr = 0;
	for (int i = 0; ; ++asi) {
	    (pwd->*pwd->stt2)(fst, asi, bzsi);
	    fst->gap += gapopen(rcd, asi, i += 3);
	    GFREQ* cf = *asi.sfq;
	    while ((d3 - cf->glen) > pwd->codonk1) ++cf;
	    uscr += cf->freq * pwd->BasicGEP;
	    if (cf > *asi.sfq)
		uscr += ((*asi.sfq)->freq - cf->freq) * pwd->LongGEP;
	    update(rcd, asi, 3);
	    if (i >= d3) break;
	}
	return uscr;
}

VTYPE SumCodePot(EXIN* bb, int i, CHAR* cs, PwdM* pwd)
{
	VTYPE	y = 0;

	for (++bb; i > 0; i -= 3, bb += 3) {
	    VTYPE	hvl;
	    if (cs) {
		hvl = pwd->pmt->prematT(cs);
		cs = 0;
	    } else
		hvl = bb->sigE;
	    y += hvl;
	}
	return (y);
}

