/*****************************************************************************
*
*	dpunit.c
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

#include "dpunit.h"

// Non-internal gap algorithm

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

// Non-internal gap algorithm for local alignment

template <>
void clear(SwgDPunit* rcd)
{
	*rcd = blank_swgdpunit;
}

template <>
void reset(SwgDPunit* rcd)
{
	*rcd = black_swgdpunit;
}
template <>
void copy(SwgDPunit* dst, const SwgDPunit* src)
{
	*dst = *src;
}

// Naive algorithm for local alignment

template <>
void clear(SwgDPunit_nv* rcd)
{
	clear((SwgDPunit*) rcd);
	vclear(rcd->gla, rcd->dtsize);
}

template <>
void reset(SwgDPunit_nv* rcd)
{
	reset((SwgDPunit*) rcd);
	vclear(rcd->gla, rcd->dtsize);
}

template <>
void copy(SwgDPunit_nv* dst, const SwgDPunit_nv* src)
{
	copy((SwgDPunit*) dst, (const SwgDPunit*) src);
	vcopy(dst->gla, src->gla, src->dtsize);
}

// Half-general profile for local alignment

template<>
void clear(SwgDPunit_hf* rcd)
{
	clear((SwgDPunit*) rcd);
	rcd->glb = 0;
	cleardelta(rcd->dla);
}

template<>
void reset(SwgDPunit_hf* rcd)
{
	reset((SwgDPunit*) rcd);
	rcd->glb = 0;
	cleardelta(rcd->dla);
}

template<>
void copy(SwgDPunit_hf* dst, const SwgDPunit_hf* src)
{
	copy((SwgDPunit*) dst, (const SwgDPunit*) src);
	copydelta(dst->dla, src->dla);
	dst->glb = src->glb;
}

// Full-general profile for local alignment

template<>
void clear(SwgDPunit_pf* rcd)
{
	clear((SwgDPunit*) rcd);
	cleardelta(rcd->dla);
	cleardelta(rcd->dlb);
}

template<>
void reset(SwgDPunit_pf* rcd)
{
	reset((SwgDPunit*) rcd);
	cleardelta(rcd->dla);
	cleardelta(rcd->dlb);
}

template<>
void copy(SwgDPunit_pf* dst, const SwgDPunit_pf* src)
{
	copy((SwgDPunit*) dst, (const SwgDPunit*) src);
	copydelta(dst->dla, src->dla);
	copydelta(dst->dlb, src->dlb);
}

