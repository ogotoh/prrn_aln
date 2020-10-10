/*****************************************************************************
*
*	Alignment of two protein/nucleotide sequences 
*	Selecter of various algorithms
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
*	Osamu Gotoh, Ph.D.      (2003-2012)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#include "aln.h"
#include "wln.h"
#include "mseq.h"
#include "maln.h"
#include "fwd2c.h"
#include "fwd2h.h"
#include "fwd2s.h"
#include "phyl.h"
#include "consreg.h"
#include "fspscore.h"
#include "autocomp.h"

#define	DEBUG	0

#define NotGap(x) ((x) && IsntGap(*(x)))
#define	TrueGap(x) (!(x) || IsTrueGap(*(x)))

static	int	by_upr(ISLAND* ald, ISLAND* bld);
static	int	by_val(ISLAND* ald, ISLAND* bld);

bool PwdM::advised_sim2(mSeq* seqs[])
{
	int	i = seqs[0]->many < seqs[1]->many;
	int	j = 1 - i;
	mSeq*	a = seqs[i];
	mSeq*	b = seqs[j];
	int	ni = a->many;
	int	nj = b->many;
	int	nt = 2 * nj + ni;

	if (nt < thr_gfq_21) return (false);
	bool	apf = nt >= thr_gfq_22 || nj == 1;
	bool	bpf = nj > b->code->max_code;
	if (i) swap(apf, bpf);
	aprof = apf;
	bprof = bpf;
	return (true);
}

static int convert_main(CalcServer<mSeq>* svr, mSeq* sqs[], ThQueue<mSeq>* q)
{
	Simmtx*	simmtx = (Simmtx*) svr->prm;
	mSeq*&	sd = sqs[0];
	if (sd->many > 1) {
	    sd->exg_seq(algmode.lcl, algmode.lcl);
	    sd->setsimmtx(simmtx);
	    sd->convseq(VECTOR);
	    sd->convseq(VECPRO);
	}
	return (OK);
}

void convertseqs(mSeq** seqs, int nos, Simmtx* sm)
{
	CalcServer<mSeq> svr1(IM_SNGL, sm, &convert_main, 0, 0, seqs, nos);
	svr1.autocomp();
}

int PwdM::selAlnMode(mSeq** seqs, bool renew)
{
	mSeq*&	a = seqs[0];
	mSeq*&	b = seqs[1];
	int	htr = (DvsP == 1 || DvsP == 2)? 5: 0;
	abgfq = aprof = bprof = 0;

	if ((a->inex.molc == GENOME && b->istron()) ||
	    (b->inex.molc == GENOME && a->istron()))
	    return (BAD_ALN);
	if (renew) {
	    algmode.lsg = a->inex.intr || b->inex.intr;	// splice
	    abgfq = advised_sim2(seqs);
	    bool	agfq = a->inex.dels;
	    bool	bgfq = b->inex.dels;
	    if (!agfq && !bgfq)	alnmode = NGP_ALN;	// No internal gap
	    else if (!abgfq)	alnmode = (htr? NGP_ALN: NTV_ALN);	// Without profile
	    else if (!agfq) alnmode = RHF_ALN;	// 1 profile & swap
	    else if (!bgfq) alnmode = HLF_ALN;	// 1 profile
	    else	alnmode = GPF_ALN; 				// 2 profiles
	    prof = aprof || bprof;
	    if (algmode.bnd) alnmode += NTV_ALN;
	    if (algmode.lsg) alnmode += 5 + htr;	// splice
	    switch (alnmode) {
		case HLF_ALN: case HLF_ALB: swp = false; break;
		case RHF_ALN: case RHF_ALB: swp = true; break;
		case GPF_ALN: case GPF_ALB:
		    swp = !aprof && bprof; break;
		case NTV_ALN:
		    swp = (a->right - a->left) < (b->right - b->left);
		    break;
		default: swp = a->inex.intr; break;
	    }
	    if (swp) {
		swap(seqs[0], seqs[1]);
		swapDvsP();
		swap(an, bn);
		INT	t = aprof; aprof = bprof; bprof = t;
	    }
	}

	if (a->len) {
	    if (algmode.lcl & 16)  // local
		a->exg_seq(1, 1);
	    else		// global
		a->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
	    if (DvsP == TxP && a->at(0) && !a->istron()) a->nuc2tron();
	    if (aprof) {
		if (!a->inex.vect) a->convseq(VECTOR);
		if (a->inex.vect != VECPRO) {
		    a->setsimmtx(simmtx);
		    a->convseq(VECPRO);
		}
	    } else	a->convseq(RAWSEQ);

	}
	if (b->len) {
	    if (algmode.lcl & 16)  // local
		b->exg_seq(1, 1);
	    else		// global
		b->exg_seq(algmode.lcl & 4, algmode.lcl & 8);
	    if (DvsP == PxT && b->at(0) && !b->istron()) b->nuc2tron();
	    if (bprof) {
		if (!b->inex.vect) b->convseq(VECTOR);
		if (!aprof && b->inex.vect != VECPRO) {
		    b->setsimmtx(simmtx);
		    b->convseq(VECPRO);
		}
	    } else	b->convseq(RAWSEQ);
	}
	a->mkthick();
	b->mkthick();
	return (alnmode);
}

void Msap::selSpMode(mSeq* sd)
{
	if (!calc_mode) {
	    if (sd->many == 1) calc_mode = BAD_ALN;
	    if (sd->many < thr_gfq_11) calc_mode = NTV_ALN; else
	    if (sd->many < thr_gfq_12) calc_mode = HLF_ALN; else
	    calc_mode = GPF_ALN;
	}
	switch (calc_mode) {
	    case NTV_ALN: sd->inex.gfrq = sd->inex.prof = 0; break;
	    case HLF_ALN: sd->inex.gfrq = 1; sd->inex.prof = 0; break;
	    default:	sd->inex.gfrq = sd->inex.prof = 1; break;
	}
}

VTYPE selfAlnScr(mSeq* sd, Simmtx* sm, FSTAT* stt)
{
	if (!sm) sm = getSimmtx(0);
	mSeq*	sqs[2] = {sd, sd};
	return calcscore_grp(sqs, stt);
}

Msap::Msap(mSeq* sd, int mode, const ALPRM* ap) 
	: msd(sd), alnprm(ap? *ap: alprm), calc_mode(mode)
{
	selSpMode(sd);
	an = sd->many;
	codonk1 = alnprm.ls > 2? alnprm.k1: LARGEN;
	diffu = alnprm.ls > 2? (VTYPE) (alnprm.scale * (alnprm.u - alnprm.u1)): 0;
	if (sd->inex.gfrq && !sd->gfq) sd->gfq = new Gfq(sd);
	a_mode = sd->unit_mode();
	simmtx = getSimmtx(alnprm.mtx_no);
	if (sd->inex.prof) sd->setsimmtx(simmtx);
	sd->convseq(sd->inex.prof);
	base_code = sd->code->base_code;
	van = sd->sumwt;
	Basic_GOP = (VTYPE) (alnprm.scale * -alnprm.v);
	eth_code = sd->eth_code;
#if USE_WEIGHT
	wta = sd->weight;
	pairwt = sd->pairwt;
	if (pairwt) {
	    FTYPE*	pw = pairwt;
	    VTYPE	sum = 0.;
	    for (int i = 0; i < an * (an - 1) / 2; ++i)
		sum += *pw++;
	    Vab = (VTYPE) sum;
	    crg1 = &Msap::crg1w;
	} else if (wta) {
	    FTYPE*	wt = wta;
	    FTYPE	w2 = 0;
	    for (int i = 0; i < sd->many; ++i, ++wt)
		w2 += *wt * *wt;
	    Vab = (VTYPE) ((sd->sumwt * sd->sumwt - w2) / 2);
	    crg1 = &Msap::crg1i;
	} else {
	    Vab = an * (an - 1) / 2;
	    crg1 = &Msap::crg1i;
	}
	if (pairwt) {
	    sim1 = (alnprm.u0 > 0)? &Msap::selfw_p: &Msap::selfw;
	} else {
	    sim1 = (alnprm.u0 > 0)? &Msap::self_p: &Msap::selfi;
	}
#else	// !USE_WEIGHT
	Vab = an * (an - 1) / 2;
	crg1 = &Msap::crg1i;
	sim1 = &Msap::selfi;
#endif	// USE_WEIGHT
}

void PwdM::resetuab()
{
	VTYPE	wa = a_mode? van: 1;
	VTYPE	wb = b_mode? vbn: 1;

	Basic_GOP = (VTYPE) (-alnprm.scale * alnprm.v);
	diff_u = (VTYPE) (alnprm.scale * (alnprm.u - alnprm.u1));
	Vab = (VTYPE) (alnprm.scale * wa * wb);
	Vthr = Vab * alnprm.thr;
	if (!a_mode) wb = wa;
	Weighted_GOP = (VTYPE) -alnprm.v;
	a_u0 = b_u0 = ab_u0 = (VTYPE) (alnprm.scale * -alnprm.u0);
	if (!a_mode) {a_u0 /= van; b_u0 /= van;}
	if (!b_mode) {b_u0 /= vbn; a_u0 /= vbn;}
	wa_u0 = van * a_u0;
	wb_u0 = vbn * b_u0;
}

void PwdM::rescale(VTYPE& scr, VTYPE tgap, FSTAT* fst)
{
	if (fst) {
	    fst->val = scr / Vab;
	    fst->gap = tgap / Vab; fst->unp /= Vab;
	    fst->mch /= Vab; fst->mmc /= Vab;
	}
}

PwdM::PwdM(mSeq** seqs, const ALPRM* alp, bool prf) 
	: PwdB((const Seq**) seqs, alp), a(seqs[0]), b(seqs[1]), prof(prf), 
	  alnprm(alp? *alp: alprm), an(a->many), bn(b->many)
{
	if (seqs == 0) {
	    Sim2 = &PwdM::sim11;
	    sim2 = &PwdM::sim11;
	    unpa = &PwdM::unp1;
	    unpb = &PwdM::unp1;
	    return;
	}
	selAlnMode(seqs);
#if USE_WEIGHT
	bool	wwt = false;
	wta = a->weight;
	wtb = b->weight;
	if (wta && !wtb && !bprof) {
	    wtb = b->weight = new FTYPE[bn];
	    for (int i = 0; i < bn; ++i) wtb[i] = 1;
	}
	if (wtb && !wta && !aprof) {
	    wta = a->weight = new FTYPE[an];
	    for (int i = 0; i < an; ++i) wta[i] = 1;
	}
	wwt = wta && wtb;
#endif	// USE_WEIGHT
	base_code = a->code->base_code;
	a_mode = aprof? 2: a->many > 1;
	b_mode = bprof? 2: b->many > 1;
	van = a->sumwt;
	vbn = b->sumwt;
	resetuab();

#if USE_WEIGHT
	if (alnprm.u0 > 0) {
	  unpa = &PwdM::unpa_p;
	  unpb = &PwdM::unpb_p;
	  switch (3 * a_mode + b_mode) {
	    case 0:
		Sim2 = wwt? &PwdM::sim11w_p: &PwdM::sim11;
		stt2 = &PwdM::stt11;
		crg2 = &PwdM::crg11;
		frw2 = &PwdM::frw11; break;
	    case 1: if (wtb)
		{Sim2 = &PwdM::sim12w_p;
		 stt2 = &PwdM::stt12w;
		 crg2 = &PwdM::crg12w;
		 frw2 = &PwdM::frw12w;} else
		{Sim2 = &PwdM::sim12i_p;
		 stt2 = &PwdM::stt12i;
		 crg2 = &PwdM::crg12i;
		 frw2 = &PwdM::frw12i;} break;
	    case 2:
		Sim2 = &PwdM::sim13_p;
		stt2 = &PwdM::stt13; break;
	    case 3: if (wta)
		{Sim2 = &PwdM::sim21w_p;
		 stt2 = &PwdM::stt21w;
		 crg2 = &PwdM::crg21w;
		 frw2 = &PwdM::frw21w;} else
		{Sim2 = &PwdM::sim21i_p;
		 stt2 = &PwdM::stt21i;
		 crg2 = &PwdM::crg21i;
		 frw2 = &PwdM::frw21i;} break;
	    case 4: if (wwt)
		{Sim2 = &PwdM::sim22w_p;
		 stt2 = &PwdM::stt22w;
		 crg2 = &PwdM::crg22w;
		 frw2 = &PwdM::frw22w;} else
		{Sim2 = &PwdM::sim22i_p;
		 stt2 = &PwdM::stt22i;
		 crg2 = &PwdM::crg22i;
		 frw2 = &PwdM::frw22i;} break;
	    case 5: if (wta)
		{Sim2 = &PwdM::sim23w_p;
		 stt2 = &PwdM::stt23w;} else
		{Sim2 = &PwdM::sim23i_p;
		 stt2 = &PwdM::stt23i;} break;
	    case 6:
		Sim2 = &PwdM::sim31_p;
		stt2 = &PwdM::stt31; break;
	    case 7: if (wtb)
		{Sim2 = &PwdM::sim32w_p;
		 stt2 = &PwdM::stt32w;} else
		{Sim2 = &PwdM::sim32i_p;
		 stt2 = &PwdM::stt32i;} break;
	    case 8:
		Sim2 = DvsP == DxD? &PwdM::sim33_np: &PwdM::sim33_p;
		stt2 = &PwdM::stt33; break;
	  }
	} else {	// don't use ether
	  unpa = &PwdM::unp1;
	  unpb = &PwdM::unp1;
	  switch (3 * a_mode + b_mode) {
	    case 0:
		Sim2 = &PwdM::sim11;
		stt2 = &PwdM::stt11; 
		crg2 = &PwdM::crg11;
		frw2 = &PwdM::frw11; break;
	    case 1: if (wtb)
		{Sim2 = &PwdM::sim12w;
		 stt2 = &PwdM::stt12w;
		 crg2 = &PwdM::crg12w;
		 frw2 = &PwdM::frw12w;} else
		{Sim2 = &PwdM::sim12i;
		 stt2 = &PwdM::stt12i;
		 crg2 = &PwdM::crg12i;
		 frw2 = &PwdM::frw12i;} break;
	    case 2:
		Sim2 = &PwdM::sim13;
		stt2 = &PwdM::stt13; break;
	    case 3: if (wta)
		{Sim2 = &PwdM::sim21w;
		 stt2 = &PwdM::stt21w;
		 crg2 = &PwdM::crg21w;
		 frw2 = &PwdM::frw21w;} else
		{Sim2 = &PwdM::sim21i;
		 stt2 = &PwdM::stt21i;
		 crg2 = &PwdM::crg21i;
		 frw2 = &PwdM::frw21i;} break;
	    case 4: if (wwt)
		{Sim2 = &PwdM::sim22w;
		 stt2 = &PwdM::stt22w;
		 crg2 = &PwdM::crg22w;
		 frw2 = &PwdM::frw22w;} else
		{Sim2 = &PwdM::sim22i;
		 stt2 = &PwdM::stt22i;
		 crg2 = &PwdM::crg22i;
		 frw2 = &PwdM::frw22i;} break;
	    case 5: if (wta)
		{Sim2 = &PwdM::sim23w;
		 stt2 = &PwdM::stt23w;} else
		{Sim2 = &PwdM::sim23i;
		 stt2 = &PwdM::stt23i;} break;
	    case 6:
		Sim2 = &PwdM::sim31;
		stt2 = &PwdM::stt31; break;
	    case 7: if (wtb)
		{Sim2 = &PwdM::sim32w;
		 stt2 = &PwdM::stt32w;} else
		{Sim2 = &PwdM::sim32i;
		 stt2 = &PwdM::stt32i;} break;
	    case 8:
		Sim2 = DvsP == DxD? &PwdM::sim33_n: &PwdM::sim33;
		stt2 = &PwdM::stt33; break;
	  }
	}
#else	// USE_WEIGHT
	if (alnprm.u0 > 0) {
	  unpa = &PwdM::unpa_p;
	  unpb = &PwdM::unpb_p;
	  switch (3 * a_mode + b_mode) {
	    case 0:
		Sim2 = &PwdM::sim11;
		stt2 = &PwdM::stt11;
		crg2 = &PwdM::crg11;
		frw2 = &PwdM::frw11; break;
	    case 1: 
		Sim2 = &PwdM::sim12i_p;
		stt2 = &PwdM::stt12i;
		crg2 = &PwdM::crg12i;
		frw2 = &PwdM::frw12i; break;
	    case 2:
		Sim2 = &PwdM::sim13_p;
		stt2 = &PwdM::stt13; break;
	    case 3:
		Sim2 = &PwdM::sim21i_p;
		stt2 = &PwdM::stt21i;
		crg2 = &PwdM::crg21i;
		frw2 = &PwdM::frw21i; break;
	    case 4:
		Sim2 = &PwdM::sim22i_p;
		stt2 = &PwdM::stt22i;
		crg2 = &PwdM::crg22i;
		frw2 = &PwdM::frw22i; break;
	    case 5:
		Sim2 = &PwdM::sim23i_p;
		stt2 = &PwdM::stt23i; break;
	    case 6:
		Sim2 = &PwdM::sim31_p;
		stt2 = &PwdM::stt31; break;
	    case 7:
		Sim2 = &PwdM::sim32i_p;
		stt2 = &PwdM::stt32i; break;
	    case 8:
		Sim2 = DvsP == DxD? &PwdM::sim33_np: &PwdM::sim33_p;
		stt2 = &PwdM::stt33; break;
	  }
	} else {
	  unpa = &PwdM::unp1;
	  unpb = &PwdM::unp1;
	  switch (3 * a_mode + b_mode) {
	    case 0:
		Sim2 = &PwdM::sim11;
		stt2 = &PwdM::stt11; 
		crg2 = &PwdM::crg11;
		frw2 = &PwdM::frw11; break;
	    case 1: 
		Sim2 = &PwdM::sim12i;
		stt2 = &PwdM::stt12i;
		crg2 = &PwdM::crg12i;
		frw2 = &PwdM::frw12i; break;
	    case 2:
		Sim2 = &PwdM::sim13;
		stt2 = &PwdM::stt13; break;
	    case 3:
		Sim2 = &PwdM::sim21i;
		stt2 = &PwdM::stt21i;
		crg2 = &PwdM::crg21i;
		frw2 = &PwdM::frw21i; break;
	    case 4:
		Sim2 = &PwdM::sim22i;
		stt2 = &PwdM::stt22i;
		crg2 = &PwdM::crg22i;
		frw2 = &PwdM::frw22i; break;
	    case 5:
		Sim2 = &PwdM::sim23i;
		stt2 = &PwdM::stt23i; break;
	    case 6:
		Sim2 = &PwdM::sim31;
		stt2 = &PwdM::stt31; break;
	    case 7:
		Sim2 = &PwdM::sim32i;
		stt2 = &PwdM::stt32i; break;
	    case 8:
		Sim2 = DvsP == DxD? &PwdM::sim33_n: &PwdM::sim33;
		stt2 = &PwdM::stt33; break;
	  }
	}
#endif	// USE_WEIGHT
	if (!prof && (aprof || bprof))
	    Sim2 = &PwdM::sim00;	// always return 0
#if SSHP
	sim2 = (sshpprm)? &PwdM::sim2_sshp: Sim2;
#else
	sim2 = Sim2;
#endif
}

// similariy for self comparison

VTYPE Msap::selfi(const mSeqItr& asi) const
{
	VTYPE 	s = 0;
const 	CHAR*	ca = asi.res;
const 	CHAR*	cb = ca;
const 	CHAR*	ct = ca + an;

	while (++cb < ct) {
	    VTYPE*	mtxas = simmtx->mtx[*cb];
	    for (const CHAR* cs = ca; cs < cb; )
		s += mtxas[*cs++];
	}
	return (s);
}

VTYPE Msap::crg1i(const int* gla, const mSeqItr& asi) const
{
const 	CHAR*	as = asi.res;
	VTYPE	g = 0;

	for (int i = 1; i < asi.many; ++i) {
	    if (IsntGap(as[i])) {
		for (int j = 0; j < i; ++j) {
		    VTYPE	au = msd->gapdensity(as + j, j);
		    if (au > 0 && gla[i] >= gla[j]) g += au;
		}
	    } else {
		VTYPE	au = msd->gapdensity(as + i, i);
		if (au > 0) {
		    for (int j = 0; j < i; ++j) 
			if (IsntGap(as[j]) && gla[i] <= gla[j]) g += au;
		}
	    }
	}
	return (g * Basic_GOP);
}

// similarity and distance functions

VTYPE PwdM::sim12i(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	cb = bsi.res;
const 	CHAR*	ct = cb + bsi.many;
	VTYPE	s = 0;
const 	VTYPE*	mtxas = simmtx->mtx[*asi.res];

	while (cb < ct) s += mtxas[*cb++];
	return (s);
}

VTYPE PwdM::sim13(const mSeqItr& asi, const mSeqItr& bsi) const
{
	return (bsi.vss[bsi.felm + *asi.res]);
}

VTYPE PwdM::sim21i(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	ct = ca + asi.many;
	VTYPE	s = 0;
	VTYPE*	mtxbs = simmtx->mtx[*bsi.res];

	while (ca < ct) s += mtxbs[*ca++];
	return (s);
}

VTYPE PwdM::sim22i(const mSeqItr& asi, const mSeqItr& bsi) const
{
	VTYPE 	s = 0;
const 	CHAR*	ca = asi.res;
const 	CHAR*	cb = bsi.res;
const 	CHAR*	ta = ca + asi.many;
const 	CHAR*	tb = cb + bsi.many;

	while (ca < ta) {
	    VTYPE*	mtxas = simmtx->mtx[*ca++];
	    for (const CHAR* cs = cb; cs < tb; )
		s += mtxas[*cs++];
	}
	return (s);
}

VTYPE PwdM::sim23i(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	ta = ca + asi.many;
const 	VTYPE*	vb = bsi.vss + bsi.felm;
	VTYPE	s = 0;

	while (ca < ta) s += vb[*ca++];
	return (s);
}

VTYPE PwdM::sim32i(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	VTYPE*	va = asi.vss + asi.felm;
const 	CHAR*	cb = bsi.res;
const 	CHAR*	tb = cb + bsi.many;
	VTYPE	s = 0;

	while (cb < tb) s += va[*cb++];
	return (s);
}

//	Now 'va' is a profile vector and 'vb' is a frequency vector.
//	If 'vb' is real type it is usually normalized.
//	If 'vb' is interger type it is not normalized.

VTYPE PwdM::sim33(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	VTYPE*	va = asi.vss + asi.felm;
const 	VTYPE*	vb = bsi.vss;
const 	VTYPE*	tb = vb + bsi.felm;
	VTYPE	s = 0;

	while (vb < tb) s += *va++ * *vb++;
	return (s);
}

VTYPE PwdM::sim33_n(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	VTYPE*	va = asi.vss + asi.felm;
const 	VTYPE*	vb = bsi.vss;
	VTYPE	s = 0;

	for (int i = 0; i < bsi.felm; ++i)
	    s += va[decompact[i]] * *vb++;
	return (s);
}

//	calculate statistics : # of matches, mismatches, gaps, unps

void PwdM::stt11(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	cb = bsi.res;
	if (NotGap(ca)) {
	    if (NotGap(cb)) {
		if (same(*ca, *cb))	++stt->mch;
		else		++stt->mmc;
	    } else if (TrueGap(ca))
		++stt->unp;
	} else if (NotGap(cb) && TrueGap(ca))
		++stt->unp;
}

void PwdM::stt12i(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	cb = bsi.res;
const 	CHAR*	ct = cb + bsi.many;
	int	s = 0;
	int	g = 0;
	int	u = 0;

	if (NotGap(ca)) {
	    if (cb) {
		for ( ; cb < ct; ++cb) {
		    if (same(*cb, *ca)) ++s;
		    else {
			if (IsGap(*cb)) ++u;
			if (IsTrueGap(*cb)) ++g;
		    }
		}
		stt->mch += s;
		stt->mmc += bsi.many - s - u;
		stt->unp += g;
	    } else
		stt->unp += bsi.dns->efq;
	} else if (cb && TrueGap(ca)) {
	    for ( ; cb < ct; ++cb)
		if (IsntGap(*cb)) ++g;
	    stt->unp += g;
	}
}

void PwdM::stt13(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	VTYPE*	vb = bsi.vss;
	if (NotGap(ca)) {
	    if (vb) {
		CHAR cca = *ca;
		if (DvsP == DxD) cca = nccmpctab[cca];
		else if (DvsP == TxP && cca == SER2) cca = SER;
		stt->mch += vb[cca];
		stt->mmc += vbn - vb[gap_code] - vb[cca] - vb[nil_code];
		stt->unp += vb[gap_code];
	    } else
		stt->unp += bsi.dns->efq;
	} else if (vb && TrueGap(ca)) {
	    stt->unp += vbn - vb[gap_code] - vb[nil_code];
	}
}

void PwdM::stt21i(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	ct = ca + asi.many;
const 	CHAR*	cb = bsi.res;
	int	s = 0;
	int	g = 0;
	int	u = 0;

	if (NotGap(cb)) {
	    if (ca) {
		for ( ; ca < ct; ++ca) {
		    if (same(*ca, *cb)) ++s; else
		    if (IsGap(*ca)) ++u;
		    if (IsTrueGap(*ca)) ++g;
		}
		stt->mch += s;
		stt->mmc += asi.many - s - u;
		stt->unp += g;
	    } else
		stt->unp += asi.dns->efq;
	} else if (ca && TrueGap(cb)) {
	    for ( ; ca < ct; ++ca)
		if (IsntGap(*ca)) ++g;
	    stt->unp += g;
	}
}

void PwdM::stt22i(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	cb = bsi.res;
const 	CHAR*	ta = ca + asi.many;
const 	CHAR*	tb = cb + bsi.many;
	int	s = 0;
	int	m = 0;
	int	g = 0;

	if (ca && cb) {
	  for ( ; cb < tb; ++cb) {
	    if (IsntGap(*cb)) {
	      for (const CHAR* cs = ca; cs < ta; ++cs) {
		if (*cs == *cb) ++s; else
		if (IsntGap(*cs)) ++m;
		if (IsTrueGap(*cs)) ++g;
	      }
	    } else if (IsTrueGap(*cb)) {
	      for (const CHAR* cs = ca; cs < ta; ++cs)
		if (IsntGap(*cs)) ++g;
	    }
	  }
	} else if (ca) {
	  for ( ; ca < ta; ++ca)
	    if (IsntGap(*ca)) g += (int) bsi.dns->efq;
	} else if (cb) {
	  for ( ; cb < tb; ++cb)
	    if (IsntGap(*cb)) g += (int) asi.dns->efq;
	}
	stt->mch += s;
	stt->mmc += m;
	stt->unp += g;
}

void PwdM::stt23i(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	VTYPE*	vb = bsi.vss;

	if (ca) {
	  const CHAR*	ta = ca + asi.many;
	  for ( ; ca < ta; ++ca) {
	    int	cca = *ca;
	    if (DvsP == DxD) cca = nccmpctab[cca];
	    else if (DvsP == TxP && cca == SER2) cca = SER;
	    if (IsntGap(cca)) {
	      if (vb) {
		stt->mch += vb[cca];
		stt->mmc += vbn - vb[gap_code] - vb[cca] - vb[nil_code];
		stt->unp += vb[gap_code];
	      } else if (IsTrueGap(cca))
		stt->unp += bsi.dns->efq;
	    } else if (vb && IsTrueGap(cca)) {
		stt->unp += vbn - vb[gap_code] - vb[nil_code];
	    }
	  }
	} else if (vb)
	    stt->unp += asi.dns->efq * (vbn - vb[gap_code] - vb[nil_code]);
}

void PwdM::stt31(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	VTYPE*	va = asi.vss;
const 	CHAR*	cb = bsi.res;
	if (NotGap(cb)) {
	    if (va) {
		CHAR	ccb = *cb;
		if (DvsP == DxD) ccb = nccmpctab[ccb];
		else if (DvsP == PxT && ccb == SER2) ccb = SER;
		stt->mch += va[ccb];
		stt->mmc += van - va[gap_code] - va[ccb] - va[nil_code];
		stt->unp += va[gap_code];
	    } else
		stt->unp += asi.dns->efq;
	} else if (va && TrueGap(cb)) {
	    stt->unp += van - va[gap_code] - va[nil_code];
	}
}

void PwdM::stt32i(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	VTYPE*	va = asi.vss;
const 	CHAR*	cb = bsi.res;

	if (cb) {
	  const CHAR*	tb = cb + bsi.many;
	  for ( ; cb < tb; ++cb) {
	    int	ccb = *cb;
	    if (DvsP == DxD) ccb = nccmpctab[ccb];
	    else if (DvsP == PxT && ccb == SER2) ccb = SER;
	    if (IsntGap(ccb)) {
	      if (va) {
		stt->mch += va[ccb];
		stt->mmc += van - va[gap_code] - va[ccb] - va[nil_code];
		stt->unp += va[gap_code];
	      } else if (IsTrueGap(ccb))
		stt->unp += asi.dns->efq;
	    } else if (va && IsTrueGap(ccb)) {
		stt->unp += van - va[gap_code] - va[nil_code];
	    }
	  }
	} else if (va)
	    stt->unp += bsi.dns->efq * (van - va[gap_code] - va[nil_code]);
}

void PwdM::stt33(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	VTYPE*	va = asi.vss;
const 	VTYPE*	vb = bsi.vss;

	if (va && vb) {
	    stt->mmc += (van - va[gap_code] - va[nil_code])
		      * (vbn - vb[gap_code] - vb[nil_code]);
	    stt->unp += va[gap_code] * (vbn - vb[gap_code] - vb[nil_code]) + 
		    (van - va[gap_code] - va[nil_code]) * vb[gap_code];
	    const VTYPE*	ta = va + asi.felm - 2;
	    VTYPE	s = 0;
	    va += base_code;
	    vb += base_code;
	    while (va < ta) s += *va++ * *vb++;
	    stt->mch += s;
	    stt->mmc -= s;
	} else if (va) {
	    stt->unp += (van - va[gap_code] - va[nil_code]) * bsi.dns->efq;
	} else if (vb) {
	    stt->unp += (vbn - vb[gap_code] - vb[nil_code]) * asi.dns->efq;
	}
}

// gap open penalty

VTYPE  PwdM::newgap2(const mSeqItr& asi, const CHAR* bs, const int* glb) const
{
	VTYPE	g = 0;
const 	CHAR*	ts = bs + bn;
#if USE_WEIGHT
const 	FTYPE*	wt = wtb;
#endif
	for ( ; bs < ts; ++bs, ++glb) {
#if USE_WEIGHT
	    FTYPE	w = wt? *wt++: 1;
#else
	    FTYPE	w = 1;
#endif
	    if (IsntGap(*bs)) {
		for (GFREQ* adf = *asi.tfq; neogfq(adf); ++adf) {
		    if (*glb < adf->glen) break;
		    g += adf->freq * w;
		}
	    } else {
		VTYPE	bu = b->gapdensity(bs, 0);
		if (bu > 0) {
		    for (GFREQ* acf = *asi.sfq; neogfq(acf); ++acf) {
			if (acf->glen >= *glb) g += acf->freq * w * bu;
			else	break;
		    }
		}
	    }
	}
	return (Weighted_GOP * g);
}

VTYPE PwdM::crg11(const int* gla, const int* glb, 
	const mSeqItr& asi, const mSeqItr& bsi, int d3) const
{
const 	CHAR*	as = asi.res;
const 	CHAR*	bs = bsi.res;
	if (d3 == 0) {
	    VTYPE	bu = b->gapdensity(bs, 0);
	    if (IsntGap(*as) && bu > 0 && *gla >= *glb) return (bu * Basic_GOP);
	    VTYPE	au = a->gapdensity(as, 0);
	    if (IsntGap(*bs) && au > 0 && *glb >= *gla) return (au * Basic_GOP);
	} else if (d3 > 0) {
	    VTYPE	bu = b->postgapdensity(bs, 0);
	    if (bu > 0 && *gla >= *glb) return (bu * Basic_GOP);
	} else {
	    VTYPE	au = a->postgapdensity(as, 0);
	    if (au > 0 && *glb >= *gla) return (au * Basic_GOP);
	}
	return (0);
}

VTYPE PwdM::crg12i(const int* gla, const int* glb, 
	const mSeqItr& asi, const mSeqItr& bsi, int d3) const
{
const 	CHAR*	as = asi.res;
const 	CHAR*	bs = bsi.res;
const 	int*	tl = glb + bsi.many;
	VTYPE	g = 0;

	if (d3 == 0) {
	    if (IsntGap(*as)) {
		for ( ; glb < tl; ++glb, ++bs) {
		    VTYPE	bu = b->gapdensity(bs, bs - bsi.res);
		    if (bu > 0 && *gla >= *glb) g += bu;
		}
	    } else {
		VTYPE	au = a->gapdensity(as, 0);
		if (au > 0) {
		    for ( ; glb < tl; ++glb, ++bs)
			if (IsntGap(*bs) && *glb >= *gla) g += au;
		}
	    }
	} else if (d3 > 0) {		/* a x --- */
	    if (IsntGap(*as)) {
		for ( ; glb < tl; ++glb, ++bs) {
		    VTYPE	bu = b->postgapdensity(bs, bs - bsi.res);
		    if (bu > 0 && *gla >= *glb) g += bu;
		}
	    }
	} else {	/* - x bbb */
	    VTYPE	au = a->postgapdensity(as, 0);
	    if (au > 0) {
		for ( ; glb < tl; ++glb, ++bs)
		    if (IsntGap(*bs) && *glb >= *gla) g += au;
	    }
	}
	return (g * Basic_GOP);
}

VTYPE PwdM::crg21i(const int* gla, const int* glb, 
	const mSeqItr& asi, const mSeqItr& bsi, int d3) const
{
const 	CHAR*	as = asi.res;
const 	CHAR*	bs = bsi.res;
const 	int*	tl = gla + asi.many;
	VTYPE	g = 0;

	if (d3 == 0) {
	    if (IsntGap(*bs)) {
		for ( ; gla < tl; ++gla, ++as) {
		    VTYPE	au = a->gapdensity(as, as - asi.res);
		    if (au > 0 && *glb >= *gla) g += au;
		}
	    } else {
		VTYPE	bu = b->gapdensity(bs, 0);
		if (bu > 0) {
		    for ( ; gla < tl; ++gla, ++as)
			if (IsntGap(*as) && *gla >= *glb) g += bu;
		}
	    }
	} else if (d3 < 0) {		/* --- x b */
	    if (IsntGap(*bs)) {
		for ( ; gla < tl; ++gla, ++as) {
		    VTYPE	au = a->postgapdensity(as, as - asi.res);
		    if (au > 0 && *glb >= *gla) g += au;
		}
	    }
	} else {	/* aaa x - */
	    VTYPE	bu = b->postgapdensity(bs, 0);
	    if (bu > 0) {
		for ( ; gla < tl; ++gla, ++as)
		    if (IsntGap(*as) && *gla >= *glb) g += bu;
	    }
	}
	return (g * Basic_GOP);
}

VTYPE PwdM::crg22i(const int* gla, const int* glb, 
	const mSeqItr& asi, const mSeqItr& bsi, int d3) const
{
const 	CHAR*	as = asi.res;
const 	CHAR*	bs = bsi.res;
const 	int*	gta = gla + asi.many;
const 	int*	gtb = glb + bsi.many;
	VTYPE	g = 0;

	if (d3 == 0) {
	    for ( ; gla < gta; ++gla, ++as) {
		const int*	glc = glb;
		if (IsntGap(*as)) {
		    for (const CHAR* cs = bs; glc < gtb; ++glc, ++cs) {
			VTYPE	bu = b->gapdensity(cs, cs - bs);
			if (bu > 0 && *gla >= *glc) g += bu;
		    }
		} else {
		    VTYPE	au = a->gapdensity(as, as - asi.res);
		    if (au > 0) {
			for (const CHAR* cs = bs ; glc < gtb; ++glc, ++cs)
			    if (IsntGap(*cs) && *glc >= *gla) g += au;
		    }
		}
	    }	
	} else if (d3 > 0) {	/* aaa x --- */
	    for ( ; gla < gta; ++gla, ++as) {
		if (IsntGap(*as)) {
		    const int*	glc = glb;
		    for (const CHAR* cs = bs; glc < gtb; ++glc, ++cs) {
			VTYPE	bu = b->postgapdensity(cs, cs - bs);
			if (bu > 0 && *gla >= *glc) g += bu;
		    }
		}
	    }
	} else {		/* --- x bbb */
	    for ( ; glb < gtb; ++glb, ++bs) {
		if (IsntGap(*bs)) {
		    const int*	glc = gla;
		    for (const CHAR* cs = as; glc < gta; ++glc, ++cs) {
			VTYPE   au = a->postgapdensity(cs, cs - as);
			if (au > 0 && *glb >= *glc) g += au;
		    }
		}
	    }
	}
	return (g * Basic_GOP);
}

// 

VTYPE PwdM::frw12i(const int* a_gla, const int* a_glb, 
	const int* b_gla, const int* b_glb) const
{
const 	int*	tb = b_glb + bn;
	int	nn = 0;

	for ( ; b_glb < tb; ++a_glb, ++b_glb) {
	    if ((*a_gla > *a_glb && *b_gla <= *b_glb) ||
		(*a_gla < *a_glb && *b_gla >= *b_glb)) ++nn;
	}
	return (-Basic_GOP * nn);
}

VTYPE PwdM::frw21i(const int* a_gla, const int* a_glb, 
	const int* b_gla, const int* b_glb) const
{
const 	int*	ta = a_gla + an;
	int	nn = 0;

	for ( ; a_gla < ta; ++a_gla, ++b_gla) {
	    if ((*a_gla > *a_glb && *b_gla <= *b_glb) ||
		(*a_gla < *a_glb && *b_gla >= *b_glb)) ++nn;
	}
	return (-Basic_GOP * nn);
}

VTYPE PwdM::frw22i(const int* a_gla, const int* a_glb, 
	const int* b_gla, const int* b_glb) const
{
const 	int*	ta = a_gla + an;
const 	int*	tb = b_glb + bn;
	int	nn = 0;

	for ( ; a_gla < ta; ++a_gla, ++b_gla) {
	    for (const int *ab = a_glb, *bb = b_glb; bb < tb; ++ab, ++bb) {
		if ((*a_gla > *ab && *b_gla <= *bb) ||
		    (*a_gla < *ab && *b_gla >= *bb)) ++nn;
	    }
	}
	return (-Basic_GOP * nn);
}

VTYPE Msap::self_p(const mSeqItr& asi) const
{
const 	CHAR*	ca = asi.res;
	VTYPE 	s = 0;
const 	CHAR*	cb = ca;
const 	CHAR*	tb = cb + asi.many;
	VTYPE	g = msd->gapdensity(cb, 0);

	while (++cb < tb) {
	    g += msd->gapdensity(cb, cb - ca);
	    VTYPE*	mtxas = simmtx->mtx[*cb];
	    for (const CHAR* cs = ca; cs < cb; )
		s += mtxas[*cs++];
	}
	return (s + ab_u0 * g);
}

VTYPE PwdM::sim12i_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	cb = bsi.res;
const 	CHAR*	tb = cb + bsi.many;
	VTYPE	s = (*ca) * a_u0;
	VTYPE	g = 0;
const 	VTYPE*	mtxas = simmtx->mtx[*ca];

	while (cb < tb) {
	    g += b->gapdensity(cb, cb - bsi.res);
	    s += mtxas[*cb++];
	}
	return (s + b_u0 * g);
}

VTYPE PwdM::sim13_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	VTYPE*	vb = bsi.vss + bsi.felm;
	VTYPE	g = vb[bsi.eth_code] * b_u0;

	g += wa_u0 * a->gapdensity(ca, 0);
	return (vb[*ca] + g);
}

VTYPE PwdM::sim21i_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	ta = ca + asi.many;
	VTYPE	s = b->gapdensity(bsi.res, 0) * b_u0;
const 	VTYPE*	mtxbs = simmtx->mtx[*bsi.res];
	VTYPE	g = 0;

	while (ca < ta) {
	    g += a->gapdensity(ca, ca - asi.res);
	    s += mtxbs[*ca++];
	}
	return (s + a_u0 * g);
}

VTYPE PwdM::sim22i_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	cb = bsi.res;
const 	CHAR*	ta = ca + asi.many;
const 	CHAR*	tb = cb + bsi.many;
	VTYPE 	s = 0;
	VTYPE	g = 0;

	while (ca < ta) {
	    g += a->gapdensity(ca, ca - asi.res);
	    VTYPE*	mtxas = simmtx->mtx[*ca++];
	    for (const CHAR* cs = cb; cs < tb; )
		s += mtxas[*cs++];
	}
	for ( ; cb < tb; ++cb)
	    g += b->gapdensity(cb, cb - bsi.res);
	return (s + ab_u0 * g);
}

VTYPE PwdM::sim23i_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	ta = ca + asi.many;
const 	VTYPE*	vb = bsi.vss + bsi.felm;
	VTYPE	s = 0;
	VTYPE	g = 0;

	while (ca < ta) {
	    g += a->gapdensity(ca, ca - asi.res);
	    s += vb[*ca++];
	}
	return (s + ab_u0 * (bsi.vss[bsi.eth_code] + g));
}

VTYPE PwdM::sim31_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
	VTYPE	g = asi.vss[asi.eth_code] * a_u0;

	g += b->gapdensity(bsi.res, 0) * wb_u0;
	return (asi.vss[asi.felm + *bsi.res] + g);
}

VTYPE PwdM::sim32i_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	VTYPE*	va = asi.vss + asi.felm;
const 	CHAR*	cb = bsi.res;
const 	CHAR*	tb = cb + bsi.many;
	VTYPE	s = 0;
	VTYPE	g = 0;

	while (cb < tb) {
	    g += b->gapdensity(cb, cb - bsi.res);
	    s += va[*cb++];
	}
	return (s + ab_u0 * (asi.vss[asi.eth_code] + g));
}

VTYPE PwdM::sim33_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	VTYPE*	va = asi.vss;
const 	VTYPE*	vb = bsi.vss;
const 	VTYPE*	tb = vb + bsi.felm;
	VTYPE	s = ab_u0 * (va[asi.eth_code] + vb[bsi.eth_code]);

	va += asi.felm;
	while (vb < tb) s += *va++ * *vb++;
	return (s);
}

VTYPE PwdM::sim33_np(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	VTYPE*	va = asi.vss;
const 	VTYPE*	vb = bsi.vss;
	VTYPE	s = ab_u0 * (va[asi.eth_code] + vb[bsi.eth_code]);

	va += asi.felm;
	for (int i = 0; i < asi.felm; ++i)
	    s += va[decompact[i]] * *vb++;
	return (s);
}

#if USE_WEIGHT

// similarity and distance

VTYPE Msap::selfw(const mSeqItr& asi) const
{
	VTYPE 	s = 0;
const 	CHAR*	ca = asi.res;
const 	CHAR*	cb = ca;
const 	CHAR*	ct = ca + asi.many;
const 	FTYPE*	pw = pairwt;

	while (++cb < ct) {
	    VTYPE*	mtxas = simmtx->mtx[*cb];
	    for (const CHAR* cs = ca; cs < cb; )
		s += mtxas[*cs++] * *pw++;
	}
	return ((VTYPE) s);
}

VTYPE PwdM::sim12w(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	cb = bsi.res;
const 	CHAR*	ct = cb + bsi.many;
	VTYPE	s = 0;
const 	VTYPE*	mtxas = simmtx->mtx[*asi.res];
const 	FTYPE*	wb = wtb;

	while (cb < ct) s += mtxas[*cb++] * *wb++;
	return (VTYPE) (s);
}

VTYPE PwdM::sim21w(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	ct = ca + asi.many;
	VTYPE	s = 0;
const 	VTYPE*	mtxbs = simmtx->mtx[*bsi.res];
const 	FTYPE*	wa = wta;

	while (ca < ct) s += mtxbs[*ca++] * *wa++;
	return (VTYPE) (s);
}

VTYPE PwdM::sim22w(const mSeqItr& asi, const mSeqItr& bsi) const
{
	VTYPE 	s = 0;
const 	CHAR*	ca = asi.res;
const 	CHAR*	cb = bsi.res;
const 	CHAR*	ta = ca + asi.many;
const 	CHAR*	tb = cb + bsi.many;
const 	FTYPE*	wa = wta;

	while (ca < ta) {
	    VTYPE*	mtxas = simmtx->mtx[*ca++];
	    FTYPE*	wb = wtb;
	    VTYPE	st = 0;
	    for (const CHAR* cs = cb; cs < tb; )
		st += mtxas[*cs++] * *wb++;
	    s += st * *wa++;
	}
	return ((VTYPE) s);
}

VTYPE PwdM::sim23w(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	ta = ca + asi.many;
const 	VTYPE*	vb = bsi.vss + bsi.felm;
	VTYPE	s = 0;
const 	FTYPE*	wa = wta;

	while (ca < ta) s += vb[*ca++] * *wa++;
	return ((VTYPE) s);
}

VTYPE PwdM::sim32w(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	VTYPE*	va = asi.vss + asi.felm;
const 	CHAR*	cb = bsi.res;
const 	CHAR*	tb = cb + bsi.many;
	VTYPE	s = 0;
const 	FTYPE*	wb = wtb;

	while (cb < tb) s += va[*cb++] * *wb++;
	return ((VTYPE) s);
}

// calculate statistics

void PwdM::stt12w(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	cb = bsi.res;
const 	CHAR*	ct = cb + bsi.many;
	VTYPE	s = 0;
	VTYPE	g = 0;
	VTYPE	u = 0;
const 	FTYPE*	wb = wtb;

	if (NotGap(ca)) {
	    if (cb) {
		for ( ; cb < ct; ++cb, ++wb) {
		    if (*cb == *ca) s += *wb;
		    else {
			if (IsGap(*cb)) u += *wb;
			if (IsTrueGap(*cb)) g += *wb;
		    }
		}
		stt->mch += s;
		stt->mmc += vbn - s - u;
		stt->unp += g;
	    } else
		stt->unp += bsi.dns->efq;
	} else if (cb && TrueGap(ca)) {
	    for ( ; cb < ct; ++wb, ++cb)
		if (IsntGap(*cb))
		    stt->unp += *wb;
	}
}

void PwdM::stt21w(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	ct = ca + asi.many;
const 	CHAR*	cb = bsi.res;
	VTYPE	s = 0;
	VTYPE	g = 0;
	VTYPE	u = 0;
const 	FTYPE*	wa = wta;

	if (NotGap(cb)) {
	    if (ca) {
		for ( ; ca < ct; ++ca, ++wa) {
		    if (*ca == *cb) s += *wa; else
		    if (IsGap(*ca)) u += *wa;
		    if (IsTrueGap(*ca)) g += *wa;
		}
		stt->mch += s;
		stt->mmc += van - s - u;
		stt->unp += g;
	    } else
		stt->unp += asi.dns->efq;
	} else if (ca && TrueGap(cb)) {
	    for ( ; ca < ct; ++wa, ++ca)
		if (IsntGap(*ca))
		    stt->unp += *wa;
	}
}

void PwdM::stt22w(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	cb = bsi.res;
const 	CHAR*	ta = ca + asi.many;
const 	CHAR*	tb = cb + bsi.many;
const 	FTYPE*	wb = wtb;

	if (ca && cb) {
	  for ( ; cb < tb; ++cb, ++wb) {
	    VTYPE	s = 0, g = 0, u = 0;
	    if (IsntGap(*cb)) {
	      const FTYPE*	wa  = wta;
	      for (const CHAR* cs = ca; cs < ta; ++cs, ++wa) {
		if (*cs == *cb) s += *wa;
		else {
		  if (IsGap(*cs)) u += *wa;
		  if (IsTrueGap(*cs)) g += *wa;
		}
	      }
	      stt->mmc += (van - s - u) * *wb;
	      stt->mch += s * *wb;
	      stt->unp += g * *wb;
	    } else if (IsTrueGap(*cb)) {
	      const FTYPE*	wa  = wta;
	      for (const CHAR* cs = ca; cs < ta; ++wa, ++cs)
		if (IsntGap(*cs)) g += *wa;
	      stt->unp += g * *wb;
	    }
	  }
	} else if (ca) {
	  for (const FTYPE* wa = wta; ca < ta; ++ca, ++wa)
	    if (IsntGap(*ca))  stt->unp += bsi.dns->efq * *wa;
	} else if (cb) {
	  for ( ; cb < tb; ++cb, ++wb)
	    if (IsntGap(*cb))  stt->unp += asi.dns->efq * *wb;
	}
}

void PwdM::stt23w(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	VTYPE*	vb = bsi.vss;

	if (ca) {
	  const CHAR*	ta = ca + asi.many;
	  for (FTYPE* wa = wta; ca < ta; ++ca, ++wa) {
	    int	cca = *ca;
	    if (DvsP == DxD) cca = nccmpctab[cca];
	    else if (DvsP == PxT && cca == SER2) cca = SER;
	    if (IsntGap(cca)) {
	      if (vb) {
		stt->mch += vb[cca] * *wa;
		stt->mmc += (vbn - vb[gap_code] - vb[cca] - vb[nil_code]) * *wa;
		stt->unp += vb[gap_code] * *wa;
	      } else if (IsTrueGap(*ca))
		stt->unp += bsi.dns->efq;
	    } else if (vb && IsTrueGap(*ca)) {
		stt->unp += (vbn - vb[gap_code] - vb[nil_code]) * *wa;
	    }
	  }
	} else if (vb)
	    stt->unp += asi.dns->efq * (vbn - vb[gap_code] - vb[nil_code]);
}

void PwdM::stt32w(FSTAT* stt, const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	VTYPE*	va = asi.vss;
const 	CHAR*	cb = bsi.res;

	if (cb) {
	  const CHAR*	tb = cb + bsi.many;
	  FTYPE*	wb = wtb;
	  for ( ; cb < tb; ++cb, ++wb) {
	    int	ccb = *cb;
	    if (DvsP == DxD) ccb = nccmpctab[ccb];
	    else if (DvsP == PxT && ccb == SER2) ccb = SER;
	    if (IsntGap(ccb)) {
	      if (va) {
		stt->mch += va[ccb] * *wb;
		stt->mmc += (van - va[gap_code] - va[ccb] - va[nil_code]) * *wb;
		stt->unp += va[gap_code] * *wb;
	      } else if (IsTrueGap(*cb))
		stt->unp += asi.dns->efq;
	    } else if (va && IsTrueGap(*cb)) {
		stt->unp += (van - va[gap_code] - va[nil_code]) * *wb;
	    }
	  }
	} else if (va)
	    stt->unp += bsi.dns->efq * (van - va[gap_code] - va[nil_code]);
}

// gap open penalty

VTYPE PwdM::crg12w(const int* gla, const int* glb, 
	const mSeqItr& asi, const mSeqItr& bsi, int d3) const
{
const 	CHAR*	as = asi.res;
const 	CHAR*	bs = bsi.res;
const 	int*	tl = glb + bsi.many;
const 	FTYPE*	wb = wtb;
	VTYPE	g = 0;

	if (d3 == 0) {
	    if (IsntGap(*as)) {
		for ( ; glb < tl; ++glb, ++bs, ++wb) {
		    VTYPE	bu = b->gapdensity(bs, bs - bsi.res);
		    if (bu > 0 && *gla >= *glb) g += *wb * bu;
		}
	    } else {
		VTYPE	au = a->gapdensity(as, 0);
		if (au > 0) {
		    for ( ; glb < tl; ++glb, ++bs, ++wb)
			if (IsntGap(*bs) && *glb >= *gla) g += *wb * au;
		}
	    }
	} else if (d3 > 0) {		/* a x --- */
	    if (IsntGap(*as)) {
		for ( ; glb < tl; ++glb, ++wb, ++bs) {
		    VTYPE	bu = b->postgapdensity(bs, bs - bsi.res);
		    if (bu > 0 && *gla >= *glb) g += *wb * bu;
		}
	    }
	} else {			/* - x bbb */
	    VTYPE	au = a->postgapdensity(as, 0);
	    if (au > 0) {
		for ( ; glb < tl; ++glb, ++bs, ++wb)
		    if (IsntGap(*bs) && *glb >= *gla) g += *wb * au;
	    }
	}
	return ((VTYPE) (g * Basic_GOP));
}

VTYPE PwdM::crg21w(const int* gla, const int* glb, 
	const mSeqItr& asi, const mSeqItr& bsi, int d3) const
{
const 	CHAR*	as = asi.res;
const 	CHAR*	bs = bsi.res;
const	int*	tl = gla + asi.many;
const 	FTYPE*	wa = wta;
	VTYPE	g = 0;

	if (d3 == 0) {
	    if (IsntGap(*bs)) {
		for ( ; gla < tl; ++gla, ++as, ++wa) {
		    VTYPE	au = a->gapdensity(as, as - asi.res);
		    if (au > 0 && *glb >= *gla) g += *wa * au;
		}
	    } else {
		VTYPE	bu = b->gapdensity(bs, 0);
		if (bu > 0) {
		    for ( ; gla < tl; ++gla, ++as, ++wa)
			if (IsntGap(*as) && *gla >= *glb) g += *wa * bu;
		}
	    }
	} else if (d3 < 0) {		/* --- x b */
	    if (IsntGap(*bs)) {
		for ( ; gla < tl; ++gla, ++wa, ++as) {
		    VTYPE	au = a->postgapdensity(as, as - asi.res);
		    if (au > 0 && *glb >= *gla) g += *wa * au;
		}
	    }
	} else if (*bs && bs[1]) {	/* aaa x - */
	    VTYPE	bu = b->postgapdensity(bs, 0);
	    if (bu > 0) {
		for ( ; gla < tl; ++gla, ++as, ++wa)
		    if (IsntGap(*as) && *gla >= *glb) g += *wa * bu;
	    }
	}
	return ((VTYPE) (g * Basic_GOP));
}

VTYPE Msap::crg1w(const int* gla, const mSeqItr& asi) const
{
const	CHAR*	as = asi.res;
const 	FTYPE*	pw = pairwt;
	VTYPE	g = 0;

	for (int i = 1; i < asi.many; ++i) {
	    if (IsntGap(as[i])) {
		for (int j = 0; j < i; ++j, ++pw) {
		    VTYPE	au = msd->gapdensity(as + j, j);
		    if (au > 0 && gla[i] >= gla[j]) g += *pw * au;
		}
	    } else {
		VTYPE	au = msd->gapdensity(as + i, i);
		if (au > 0) {
		    for (int j = 0; j < i; ++j, ++pw)
			if (IsntGap(as[j]) && gla[i] <= gla[j]) g += *pw * au;
		} else pw += i;
	    } 
	}
	return ((VTYPE) (g * Basic_GOP));
}

VTYPE PwdM::crg22w(const int* gla, const int* glb, 
	const mSeqItr& asi, const mSeqItr& bsi, int d3) const
{
const 	CHAR*	as = asi.res;
const 	CHAR*	bs = bsi.res;
const 	int*	gta = gla + asi.many;
const 	int*	gtb = glb + bsi.many;
const 	FTYPE*	wa = wta;
const 	FTYPE*	wb = wtb;
	VTYPE	g = 0;

	if (d3 == 0) {
	    for ( ; gla < gta; ++gla, ++as, ++wa) {
		VTYPE	s = 0;
		const int*	glc = glb;
		wb = wtb;
		if (IsntGap(*as)) {
		    for (const CHAR* cs = bs; glc < gtb; ++glc, ++cs, ++wb) {
			VTYPE	bu = b->gapdensity(cs, cs - bs);
			if (bu > 0 && *gla >= *glc) s += *wb * bu;
		    }
		    g += s * *wa;
		} else {
		    VTYPE	au = a->gapdensity(as, as - asi.res);
		    if (au > 0) {
			for (const CHAR* cs = bs; glc < gtb; ++glc, ++cs, ++wb)
			    if (IsntGap(*cs) && *glc >= *gla) s += *wb * au;
		    }
		    g += s * *wa;
		}
	    }
	} else if (d3 > 0) {	/* aaa x --- */
	    for ( ; gla < gta; ++gla, ++as, ++wa) {
		if (IsntGap(*as)) {
		    wb = wtb;
		    const int*	glc = glb;
		    VTYPE	s = 0;
		    for (const CHAR* cs = bs; glc < gtb; ++glc, ++wb, ++cs) {
			VTYPE   bu = b->postgapdensity(cs, cs - bs);
			if (bu > 0 && *gla >= *glc) s += *wb * bu;
		    }
		    g += s * *wa;
		}
	    }
	} else {		/* --- x bbb */
	    for ( ; glb < gtb; ++glb, ++bs, ++wb) {
		if (IsntGap(*bs)) {
		    wa = wta;
		    const int*	glc = gla;
		    VTYPE	s = 0;
		    for (const CHAR* cs = as; glc < gta; ++glc, ++wa, ++cs) {
			VTYPE	au = a->postgapdensity(cs, cs - as);
			if (au > 0 && *glb >= *glc) s += *wa * au;
		    }
		    g += s * *wb;
		}
	    }
	}
	return ((VTYPE) (g * Basic_GOP));
}

VTYPE PwdM::frw12w(const int* a_gla, const int* a_glb, 
	const int* b_gla, const int* b_glb) const
{
const 	int*	tb = b_glb + bn;
	VTYPE	t = 0;

	for (FTYPE* wb = wtb; b_glb < tb; ++a_glb, ++b_glb, ++wb) {
	    if ((*a_gla > *a_glb && *b_gla <= *b_glb) ||
		(*a_gla < *a_glb && *b_gla >= *b_glb)) t += *wb;
	}
	return ((VTYPE) (-Basic_GOP * t));
}

VTYPE PwdM::frw21w(const int* a_gla, const int* a_glb, 
	const int* b_gla, const int* b_glb) const
{
const 	int*	ta = a_gla + an;
	VTYPE	t = 0;

	for (FTYPE* wa = wta; a_gla < ta; ++a_gla, ++b_gla, ++wa) {
	    if ((*a_gla > *a_glb && *b_gla <= *b_glb) ||
		(*a_gla < *a_glb && *b_gla >= *b_glb)) t +=  *wa;
	}
	return ((VTYPE) (-Basic_GOP * t));
}

VTYPE PwdM::frw22w(const int* a_gla, const int* a_glb, 
	const int* b_gla, const int* b_glb) const
{
const 	int*	ta = a_gla + an;
const 	int*	tb = b_glb + bn;
	VTYPE	g = 0;

	for (FTYPE* wa = wta; a_gla < ta; ++a_gla, ++b_gla, ++wa) {
	    const VTYPE*	wb = wtb;
	    VTYPE	t = 0;
	    for (const int *ab = a_glb, *bb = b_glb; bb < tb; ++ab, ++bb, ++wb) {
		if ((*a_gla > *ab && *b_gla <= *bb) ||
		    (*a_gla < *ab && *b_gla >= *bb)) t += *wb;
	    }
	    g += *wa * t;
	}
	return ((VTYPE) (-Basic_GOP * g));
}

VTYPE Msap::selfw_p(const mSeqItr& asi) const
{
const 	CHAR*	ca = asi.res;
	VTYPE 	s = 0;
const 	CHAR*	cb = ca;
const 	CHAR*	tb = cb + asi.many;
const 	FTYPE*	wt = wta;
const 	FTYPE*	pw = pairwt;
	VTYPE	g = msd->gapdensity(cb, 0) * *wt;

	while (++cb < tb) {
	    ++wt;
	    g += msd->gapdensity(cb, cb - ca) * *wt;
	    VTYPE*	mtxas = simmtx->mtx[*cb];
	    for (const CHAR* cs = ca; cs < cb; )
		s += mtxas[*cs++] * *pw++;
	}
	return ((VTYPE) (s + ab_u0 * g));
}

VTYPE PwdM::sim11w_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	cb = bsi.res;
	VTYPE	d = a->gapdensity(ca, 0) * wa_u0;

	d += a->gapdensity(cb, 0) * wb_u0;
	return (simmtx->mtx[*ca][*cb] + d);
}

VTYPE PwdM::sim12w_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	cb = bsi.res;
const 	CHAR*	tb = cb + bsi.many;
	VTYPE	s = a->gapdensity(asi.res, 0) * wa_u0;
	VTYPE	g = 0;
const 	VTYPE*	mtxas = simmtx->mtx[*asi.res];

	for (FTYPE* wb = wtb; cb < tb; ) {
	    g += b->gapdensity(cb, cb - bsi.res) * *wb;
	    s += mtxas[*cb++] * *wb++;
	}
	return ((VTYPE) (s + b_u0 * g));
}

VTYPE PwdM::sim21w_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	ta = ca + asi.many;
	VTYPE	s = b->gapdensity(bsi.res, 0) * wb_u0;
const 	VTYPE*	mtxbs = simmtx->mtx[*bsi.res];
	VTYPE	g = 0;

	for (FTYPE* wa = wta; ca < ta; ) {
	    g += a->gapdensity(ca, ca - asi.res) * *wa;
	    s += mtxbs[*ca++] * *wa++;
	}
	return ((VTYPE) (s + a_u0 * g));
}

VTYPE PwdM::sim22w_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	cb = bsi.res;
const 	CHAR*	ta = ca + asi.many;
const 	CHAR*	tb = cb + bsi.many;
	VTYPE 	s = 0;
	VTYPE	g = 0;

	for (const FTYPE* wa = wta; ca < ta; ) {
	    g += a->gapdensity(ca, ca - asi.res) * *wa;
	    const VTYPE*	mtxas = simmtx->mtx[*ca++];
	    const FTYPE*	wb = wtb;
	    VTYPE	st = 0;
	    for (const CHAR* cs = cb; cs < tb; )
		st += mtxas[*cs++] * *wb++;
	    s += st * *wa++;
	}
	for (const FTYPE* wb = wtb; cb < tb; ++cb, ++wb)
	    g += b->gapdensity(cb, cb - bsi.res) * *wb;
	return ((VTYPE) (s + ab_u0 * g));
}

VTYPE PwdM::sim23w_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	ca = asi.res;
const 	CHAR*	ta = ca + asi.many;
const 	FTYPE*	vb = bsi.vss + bsi.felm;
	VTYPE	s = 0;
	VTYPE	g = bsi.vss[bsi.eth_code];

	for (FTYPE* wa = wta; ca < ta; ) {
	    g += a->gapdensity(ca, ca - asi.res) * *wa;
	    s += vb[*ca++] * *wa++;
	}
	return ((VTYPE) (s + ab_u0 * g));
}

VTYPE PwdM::sim32w_p(const mSeqItr& asi, const mSeqItr& bsi) const
{
const 	CHAR*	cb = bsi.res;
const 	CHAR*	tb = cb + bsi.many;
	VTYPE	s = 0;
	VTYPE	g = asi.vss[asi.eth_code];
const 	VTYPE*	va = asi.vss + asi.felm;

	for (FTYPE* wb = wtb; cb < tb; ) {
	    g += b->gapdensity(cb, cb - bsi.res) * *wb;
	    s += va[*cb++] * *wb++;
	}
	return ((VTYPE) (s + ab_u0 * g));
}

#endif	// USE_WEIGHT

#if SSHP

VTYPE	PwdM::sim2_sshp(const mSeqItr& asi, const mSeqItr& bsi) const
{
	VTYPE	scr = 0;
const 	VTYPE*	ash = asi.sshp;
const 	VTYPE*	bsh = bsi.sshp;
	if (ash && bsh) {
	    for (int i = 0; i < sshpprm->sndstates; ++i)
		scr += (*ash++) * (*bsh++) * alprm3.scnd;
	    for (int i = 0; i < sshpprm->hphstates; ++i)
		scr += (*ash++) * (*bsh++) * alprm3.hydr;
	    for (int i = 0; i < sshpprm->hmtstates; ++i)
		scr += (*ash++) * (*bsh++) * alprm3.hpmt;
	}
	return ((this->*this->Sim2)(asi, bsi) + scr);
}

#endif	// SSHP

/******************************************************
	pairsum --	Calculate pair-of-sum measure 
	by the most efficient method
******************************************************/

VTYPE Msap::ps_nml()
{
	VTYPE	scr = 0;
	int*	gla = new int[msd->many];
	Gep1st*	gep = alnprm.ls > 2? new Gep1st(msd, codonk1): 0;
	mSeqItr	asi(msd, msd->left);
	mSeqItr	tsi(msd, msd->right);
	VTYPE	lunp = 0;

	msd->pregap(gla);
	for ( ; asi < tsi; ++asi) {
	    scr += (this->*this->crg1)(gla, asi)
		+  (this->*this->sim1)(asi);
	    incrgap(gla, asi.res, msd->many);
	    if (gep) {
		lunp += gep->longup(asi.res, asi.pos, gla);
		gep->shift(asi.res, asi.pos);
	    }
	}
	if (gep) scr += lunp * diffu;
	delete[] gla;
	delete gep;
	return (scr);
}

VTYPE pairsum(mSeq* sd)
{
	if (sd->many == 1) return (0);

	Msap	msd(sd);
	VTYPE	ps = msd.ps_nml();
	if (alprm.u0 != 0)
	    ps -= (VTYPE) (alprm.u0 * sd->countunps());
	return (ps);
}

VTYPE HomScore(mSeq* seqs[], PwdM* pwdm, long rr[])
{
	mSeq*&	a = seqs[0];
	mSeq*&	b = seqs[1];
	VTYPE	scr = 0;

	if (a->left == a->right || b->left == b->right) {
	    if (rr) {
		rr[0] = b->left - a->left;
		rr[1] = b->right - a->right;
	    }
	    return (0);
	}

	switch (pwdm->alnmode) {
	    case BAD_ALN: fatal(warn_mess2);
	    case NGP_ALB: scr = HomScoreC<DPunit>(seqs, pwdm, rr); break;
	    case HLF_ALB:
	    case RHF_ALB: scr = HomScoreC<DPunit_hf>(seqs, pwdm, rr); break;
	    case GPF_ALB: scr = HomScoreC<DPunit_pf>(seqs, pwdm, rr); break;
	    case NTV_ALB: scr = HomScoreC<DPunit_nv>(seqs, pwdm, rr); break;
	    case NGP_ALN: scr = HomScoreC<DPunit>(seqs, pwdm, rr, true); break;
	    case NTV_ALN: scr = HomScoreC<DPunit_nv>(seqs, pwdm, rr, true); break;
	    case HLF_ALN:
	    case RHF_ALN: scr = HomScoreC<DPunit_hf>(seqs, pwdm, rr, true); break;
	    case GPF_ALN: scr = HomScoreC<DPunit_pf>(seqs, pwdm, rr, true); break;
	    case NGP_ALH: scr = HomScoreH<RVPDJ_nv>(seqs, pwdm); break;
	    case HLF_ALH:
	    case RHF_ALH: scr = HomScoreH<RVPDJ_hf>(seqs, pwdm); break;
	    case NGP_ALS: scr = HomScoreS<RVPDJ_nv>(seqs, pwdm); break;
	    case HLF_ALS:
	    case RHF_ALS: scr = HomScoreH<RVPDJ_hf>(seqs, pwdm); break;
	    default:
		fatal("Mode %d is not supported !\n", pwdm->alnmode);
	}
	return (scr);
}

SKL* align2(mSeq* seqs[], PwdM* pwdm, VTYPE* scr, Gsinfo* GsI)
{
	mSeq*	a = seqs[0];
	mSeq*	b = seqs[1];
	SKL*	skl = 0;

retry:
	if (a->left == a->right || b->left == b->right) {
	    skl = nogap_skl(a, b);
	} else if (algmode.qck & 1) {
	  switch (pwdm->alnmode) {
	    case BAD_ALN: prompt(warn_mess2); break; 
	    case NGP_ALB: case HLF_ALB: case RHF_ALB: case GPF_ALB: case NTV_ALB:
		skl = alignC<DPunit>(seqs, pwdm, scr); break;
	    case NGP_ALN: case NTV_ALN: case HLF_ALN: case RHF_ALN: case GPF_ALN:
		skl = alignC<DPunit_nv>(seqs, pwdm, scr, true); break;
	    case NGP_ALH: case HLF_ALH: case RHF_ALH: 
		skl = alignH<RVPDJ_nv>(seqs, pwdm, scr); break;
	    case NGP_ALS: case HLF_ALS: case RHF_ALS:
		skl = alignS<RVPDJ_nv>(seqs, pwdm, scr); break;
	    default: 
		fatal("Mode %d is not supported !\n", pwdm->alnmode);
	  }
	} else {
	  switch (pwdm->alnmode) {
	    case BAD_ALN: prompt(warn_mess2); break; 
	    case NGP_ALB: skl = alignC<DPunit>(seqs, pwdm, scr); break;
	    case HLF_ALB:
	    case RHF_ALB: skl = alignC<DPunit_hf>(seqs, pwdm, scr); break;
	    case GPF_ALB: skl = alignC<DPunit_pf>(seqs, pwdm, scr); break;
	    case NTV_ALB: skl = alignC<DPunit_nv>(seqs, pwdm, scr); break;
	    case NGP_ALN: skl = alignC<DPunit>(seqs, pwdm, scr, true); break;
	    case NTV_ALN: skl = alignC<DPunit_nv>(seqs, pwdm, scr, true); break;
	    case HLF_ALN:
	    case RHF_ALN: skl = alignC<DPunit_hf>(seqs, pwdm, scr, true); break;
	    case GPF_ALN: skl = alignC<DPunit_pf>(seqs, pwdm, scr, true); break;
	    case NGP_ALH: skl = alignH<RVPDJ_nv>(seqs, pwdm, scr); break;
	    case HLF_ALH:
	    case RHF_ALH: skl = alignH<RVPDJ_hf>(seqs, pwdm, scr); break;
	    case NGP_ALS: skl = alignS<RVPDJ_nv>(seqs, pwdm, scr); break;
	    case HLF_ALS:
	    case RHF_ALS: skl = alignS<RVPDJ_hf>(seqs, pwdm, scr); break;
	    default: 
		fatal("Mode %d is not supported !\n", pwdm->alnmode);
	  }
	}

	if (!skl || !skl->n || !GsI) {
	    delete[] skl; return (0);
	} else {			// do traceback 
	    if (pwdm->alnmode >= NGP_ALS) {
		GsI->skl = skl;
		if (skl->m) {
		  switch (pwdm->alnmode) {
		    case NGP_ALH: skl_rngH<RVPDJ_nv>(seqs, GsI, pwdm); break;
		    case RHF_ALH:
		    case HLF_ALH: skl_rngH<RVPDJ_hf>(seqs, GsI, pwdm); break;
		    case NGP_ALS: skl_rngS<RVPDJ_nv>(seqs, GsI, pwdm); break;
		    case RHF_ALS:
		    case HLF_ALS: skl_rngS<RVPDJ_hf>(seqs, GsI, pwdm); break;
		    default:
			fprintf(stderr, "Score Mode %d is not supported!\n",
			pwdm->alnmode);
			break;
		  }
		  GsI->eiscr2rng();	// add bonus to introns
		  GsI->scr -= pwdm->GapPenalty(Ip_equ_k) * (GsI->noeij - 1);
		}
	    } else {
		GsI->skl = stdskl(&skl);
		int	num = skl->n;
		if (skl[1].m != a->left || skl[num].m != a->right || 
		    skl[1].n != b->left || skl[num].n != b->right) {
		    delete[] skl;
//		    if (pwdm->alnprm.sh == -100) fatal("align2 error !\n");
		    pwdm->alnprm.sh = -100;
		    goto retry;
		}
#if DEBUG
		if (badskl(GsI->skl, (Seq**) seqs)) {
		    delete[] skl;
		    FILE*   fd = fopen("a", "w");
		    a->typeseq(fd);
		    fclose(fd);
		    fd = fopen("b", "w");
		    b->typeseq(fd);
		    fclose(fd);
		    skl = alignC<DPunit_pf>(seqs, pwdm, scr);
		    fatal("Alignment failed !\n");
		}
#endif
		PreSpScore	pss(seqs, pwdm);
		pss.calcSpScore(GsI);
		if (OutPrm.trimend) skl = trimskl((const Seq**) seqs, skl);
		GsI->skl = 0;
	    }
	}
	return (skl);
}

SKL* swg2nd(mSeq* seqs[], PwdM* pwd, Gsinfo* gsi, COLONY* clny)
{
	if (algmode.qck & 1) {
	  switch (pwd->alnmode) {
	    case BAD_ALN: prompt(warn_mess2); break; 
	    case NGP_ALB: case HLF_ALB: case RHF_ALB: case GPF_ALB: case NTV_ALB:
		return swg2ndC<DPunit>(seqs, pwd, gsi, clny);
	    default: 
		fatal("Mode %d is not supported !\n", pwd->alnmode);
	  }
	} else {  
	  switch (pwd->alnmode) {
	    case BAD_ALN: prompt(warn_mess2); break; 
	    case NGP_ALB: return swg2ndC<DPunit>(seqs, pwd, gsi, clny);
	    case HLF_ALB:
	    case RHF_ALB: return swg2ndC<DPunit_hf>(seqs, pwd, gsi, clny);
	    case GPF_ALB: return swg2ndC<DPunit_pf>(seqs, pwd, gsi, clny);
	    case NTV_ALB: return swg2ndC<DPunit_nv>(seqs, pwd, gsi, clny);
	    default: 
		fatal("Mode %d is not supported !\n", pwd->alnmode);
	  }
	}
	return (0);
}

Colonies* swg1st(mSeq** seqs, PwdM* pwd)
{
	if (seqs[0]->left == seqs[0]->right || seqs[1]->left == seqs[1]->right)
	    return (0);
	if (algmode.qck & 1) {
	  switch (pwd->alnmode) {
	    case BAD_ALN: prompt(warn_mess2); break; 
	    case NGP_ALB: case HLF_ALB: case RHF_ALB: case GPF_ALB: case NTV_ALB:
		return swg1stC<SwgDPunit>(seqs, pwd);
	    default: 
		fatal("Mode %d is not supported !\n", pwd->alnmode);
	  }
	} else {
	  switch (pwd->alnmode) {
	    case BAD_ALN: prompt(warn_mess2); break; 
	    case NGP_ALB: return swg1stC<SwgDPunit>(seqs, pwd);
	    case HLF_ALB:
	    case RHF_ALB: return swg1stC<SwgDPunit_hf>(seqs, pwd);
	    case GPF_ALB: return swg1stC<SwgDPunit_pf>(seqs, pwd);
	    case NTV_ALB: return swg1stC<SwgDPunit_nv>(seqs, pwd);
	    default: 
		fatal("Mode %d is not supported !\n", pwd->alnmode);
	  }
	}
	return (0);
}

mSeq* syntheseq(mSeq* sd, mSeq** sqs, SKL* skl)
{
	GAPS*   gps[2];

	skl2gaps(gps, skl, true);
	unfoldgap(gps[0], 1, true);
	unfoldgap(gps[1], 1, true);
	sd = aggregate(sd, sqs, gps, 0);
	int     n = 0;
	for ( ; n < sqs[0]->many; ++n) {
	    sd->nbr[n] = sqs[0]->nbr[n];
	    sd->sname->push((*sqs[0]->sname)[n]);
	}
	for (int i = 0; i < sqs[1]->many; ++i, ++n) {
	    sd->nbr[n] = sqs[1]->nbr[i];
	    sd->sname->push((*sqs[1]->sname)[i]);
	}
	delete[] gps[0]; delete[] gps[1];
	return (sd);
}

mSeq* align2(mSeq** sqs, mSeq* msd, ALPRM* alp)
{
	Gsinfo	alninf;
	VTYPE	dst;
	PwdM	pwd(sqs, alp);
	SKL*	skl = align2(sqs, &pwd, &dst, &alninf);
	SKL*	nsk = skl;
	msd = syntheseq(msd, sqs, nsk);
	if (nsk != skl) delete[] nsk;
	delete[] skl;
	swap(msd->sigII, alninf.sigII);
	if (pwd.swp) swap(sqs[0], sqs[1]);
	return (msd);
}

mSeq* align2(mSeq* a, mSeq* b, mSeq* msd, ALPRM* alp)
{
	mSeq*	sqs[3] = {a, b, 0};
	return align2(sqs, msd, alp);
}

FTYPE alnscore2dist(mSeq* sqs[], PwdM* pwd, int* ends, FTYPE denome)
{
	mSeq*&	a = sqs[0];
	mSeq*&	b = sqs[1];
	VTYPE	scr;
	Gsinfo	gsi;
	SKL*	skl = align2(sqs, pwd, &scr, &gsi);
	if (!skl) return (fInfinit);
	if (algmode.lcl) trimskl((const Seq**) sqs, stdskl(&skl));
	int	al = skl[1].m;
	int	bl = skl[1].n;
	int	ar = skl[skl->n].m;
	int	br = skl[skl->n].n;
	delete[] skl;
	if (ends) {
	    ends[0] = bl - al;
	    ends[1] = br - ar;
	}
	scr /= pwd->Vab;
	if (algmode.lcl || algmode.nsa | 2) return (scr);
	if (!denome)	denome = sqrt((selfAlnScr(a, pwd->simmtx) / (a->sumwt * a->sumwt)
		* selfAlnScr(b, pwd->simmtx) / (b->sumwt * b->sumwt)));
	FTYPE	dist = (FTYPE) scr + alprm.u * abs(ar - al - br + bl);
	return (denome < fepsilon)? fInfinit: (1 - dist / denome);
}
	
static int by_upr(ISLAND* ald, ISLAND* bld)
{
	return (ald->upr - bld->upr);
}

static int by_val(ISLAND* ald, ISLAND* bld)
{
	if (ald->val < bld->val) return (1);
	return ((ald->val > bld->val)? -1: 0);
}

void report_ild(mSeq* seqs[], ISLAND* ild, INT n, int ud)
{
static	const char	fmt[] = "%.1f\t%s %d %d %c\t%s %d %d %c\n";
static	const char	smark[] = ">-^<";
	mSeq*	a = seqs[0];
	mSeq*	b = seqs[1];
	ISLAND*	tld = ild + n;

	if (n == 0) return;
	qsort((UPTR) ild, n, sizeof(ISLAND), (CMPF) by_upr);
	for (ISLAND* ald = ild; ald < tld - 1; ++ald) {// remove overlap
	    if (ald->val == NEVSEL) continue;
	    for (ISLAND* bld = ald; ++bld < tld; ) {
		if (bld->lwr > ald->upr) break;
		if (ald->val > bld->val) bld->val = NEVSEL;
		else		ald->val = NEVSEL;
	    }
	}
	qsort((UPTR) ild, n, sizeof(ISLAND), (CMPF) by_val);
	for (ISLAND* ald = ild; ald < tld; ++ald) {
	    if (ald->val == NEVSEL) break;
	    int	ml = a->left;
	    int	nl = b->left;
	    int	r = nl - ml;
	    if (ald->lwr > r)   nl = a->left + ald->lwr;
	    else		ml = b->left - ald->lwr;
	    int	mr = a->right;
	    int	nr = b->right;
	    r = nr - mr;
	    if (ald->upr > r)   mr = b->right - ald->upr;
	    else		nr = a->right + ald->upr;
	    if (ud) {
		fprintf(stderr, fmt, float(ald->val),
		b->sqname(), b->SiteNo(nl), b->SiteNo(nr - 1), smark[b->inex.sens],
		a->sqname(), a->SiteNo(ml), a->SiteNo(mr - 1), smark[a->inex.sens]);
	    } else {
		fprintf(stderr, fmt, float(ald->val),
		a->sqname(), a->SiteNo(ml), a->SiteNo(mr - 1), smark[a->inex.sens],
		b->sqname(), b->SiteNo(nl), b->SiteNo(nr - 1), smark[b->inex.sens]);
	    }
	}
}

Wilip::Wilip(mSeq* seqs[], const PwdB* pwd, INT level)
	: top(0), nwlu(0), wlu(0)
{
	Seq*	temp[2];
	temp[0] = seqs[0]->many > 1? seqs[0]->consenseq(): seqs[0];
	temp[1] = seqs[1];
	Wlp	wln((const Seq**) temp, pwd, level);
	JUXT*	jxt = wln.run_dmsnno(nwlu);
	if (!jxt) return;
	wlu = wln.willip(&top, nwlu, jxt);
	if (seqs[0]->many > 1) delete temp[0];
}
