/*****************************************************************************
*
*	Quickly calculate distance between sequences 
*	Using oligomer compositions
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
#include "mseq.h"
#include "autocomp.h"
#include "qdiv.h"

extern	int	main(int argc, const char** argv);
inline	int	numcomp(const KVpair<INT, int>* a, const KVpair<INT, int>* b)
	{return (a->key - b->key);}

static	const	WCPRM	def_wcp_nc = {4, 10, QJuCa, 0, 0, 0};
static	const	WCPRM	def_wcp_aa = {20, 5, Qpamd, 0, 0, 0};
static	WCPRM	wcp = def_wcp_nc;

WCPRM* setqdiv(INT k, const char* ap, const char* bp, DistMes dm)
{
	if (k > 0 && k < MaxBitPat) wcp.Ktuple = k;
	if (ap) wcp.reducap = ap;
	if (bp) wcp.sBitPat = bp;
	wcp.distm = dm;
	return (&wcp);
}

WCPRM* resetqdiv(Seq* sd)
{
	if (sd->isprotein()) {
	    if (wcp.Nalpha == def_wcp_nc.Nalpha) wcp.Nalpha = def_wcp_aa.Nalpha;
	    if (wcp.Ktuple == def_wcp_nc.Ktuple) wcp.Ktuple = def_wcp_aa.Ktuple;
	    if (wcp.distm == def_wcp_nc.distm) wcp.distm = def_wcp_aa.distm;
	} else {
	    if (wcp.Nalpha == def_wcp_aa.Nalpha) wcp.Nalpha = def_wcp_nc.Nalpha;
	    if (wcp.Ktuple == def_wcp_aa.Ktuple) wcp.Ktuple = def_wcp_nc.Ktuple;
	    if (wcp.distm == def_wcp_aa.distm) wcp.distm = def_wcp_nc.distm;
	}
	return (&wcp);
}

//	construct word hash table

void Kcomp::makeHash(Seq* sd, int i, Bitpat_wq* bpp)
{
	CHAR*	ps = sd->at(sd->left);
	CHAR*	ts = sd->at(sd->right - bpp->width);

	bpp->clear();
	if (i < sd->many) {
	    for (ps += i; ps < ts; ps += sd->many) {
		INT	c = ConvTab[*ps];
		if (bpp->good(c)) {
		    INT	w = bpp->word(c);
		    if (bpp->flawless()) incr(w);
		} else
		    bpp->flaw();
	    }
	} else {
	    for (int k = 0; k < sd->many; ++k) 
		makeHash(sd, k, bpp);
	}
}

Kcomp::Kcomp(Seq* sd, INT sz, CHAR* ct, int i, int m) 
	: Dhash<INT, int>(sz), ConvTab(ct), sid(i)
{
	many = (m < sd->many)? 1: sd->many;
	Bitpat_wq	bpp(wcp.Nalpha, 1, false, bitmask(wcp.Ktuple), wcp.sBitPat);
	makeHash(sd, m, &bpp);
	kmers = press(&UnqKmers, &TtlKmers);
	qsort((UPTR) hash, UnqKmers, sizeof(KVpair<INT, int>), (CMPF) numcomp);
}

static int qdiv_main(CalcServer<Seq>* svr, Seq* seqs[], ThQueue<Seq>* q = 0)
{
	Seq*&	sd = seqs[0];
	Kcomps*	kcp = (Kcomps*) svr->prm;
	INT	sz = (sd->right - sd->left) * sd->many;
	Kcomp*	kc = new Kcomp(sd, sz, kcp->iConvTab, sd->sid, sd->many);
	m_thread_Lock(q);
	kcp->sname->push(sd->sqname());
	m_thread_Unlock(q);
	kcp->nmidx[sd->sid] = kcp->inodr++;
	kcp->kcomp[sd->sid] = kc;
	return (OK);
}

Kcomps::Kcomps(Seq* msd, int im, Subset* ss)
	: ReducWord(msd, wcp.Nalpha, wcp.reducap), input_mode(im)
{
	num = ss? ss->num: msd->many;
	kcomp = new Kcomp*[num];
	nmidx = new int[num];
	inodr = 0;
	sname = new Strlist;
	int	nn = im == IM_EVRY? ncomb(num): num;
	dist = new FTYPE[nn];
	vclear(dist, nn);
	vthr = (FTYPE) (alprm.thr > 0.? alprm.thr: FLT_MAX);
	AlnServer<Seq> svr(msd, IM_SNGL, ss, 0, this, &qdiv_main);
	svr.auto_comp();
}

Kcomps::Kcomps(Seq* sd, int argc, const char** argv, int nn, AlnServer<Seq>* ssvr)
	: ReducWord(sd, wcp.Nalpha, wcp.reducap), num(nn)
{
	if (!num) fatal("Specify data or %s", no_group_2);
	kcomp = new Kcomp*[num];
	nmidx = new int[num];
	inodr = 0;
	sname = new Strlist;
	vthr = (FTYPE) (alprm.thr > 0.? alprm.thr: FLT_MAX);
	AlnServer<Seq> ksvr(argc, argv, IM_SNGL, ssvr->catalog, this, &qdiv_main);
	ksvr.auto_comp();
	nn = ssvr->calcsize(num, ssvr->calc_mode);
	dist = new FTYPE[nn];
	vclear(dist, nn);
}

Kcomps::Kcomps(Seq* sd, int argc, const char** argv, 
	const char* cat, int nn, int imn)
	: ReducWord(sd, wcp.Nalpha, wcp.reducap), num(nn)
{
	if (!num) fatal("Specify data or %s", no_group_2);
	kcomp = new Kcomp*[num];
	nmidx = new int[num];
	inodr = 0;
	sname = new Strlist;
	vthr = (FTYPE) (alprm.thr > 0.? alprm.thr: FLT_MAX);
	AlnServer<Seq> ksvr(argc, argv, IM_SNGL, cat, this, &qdiv_main);
	ksvr.auto_comp();
	nn = ksvr.calcsize(num, imn);
	dist = new FTYPE[nn];
	vclear(dist, nn);
}

Kcomps::Kcomps(Seq* sqs[], int nn, int im) : 
	ReducWord(sqs[0], wcp.Nalpha, wcp.reducap), num(nn), input_mode(im)
{
	kcomp = new Kcomp*[num];
	nmidx = new int[num];
	inodr = 0;
	sname = new Strlist;
	nn = im == IM_EVRY? ncomb(num): num;
	dist = new FTYPE[nn + 1];
	vclear(dist, nn + 1);
	vthr = (FTYPE) (alprm.thr > 0.? alprm.thr: FLT_MAX);
	AlnServer<Seq> svr(sqs, num, IM_SNGL, this, &qdiv_main);
	svr.auto_comp(sqs);
}

Kcomps::~Kcomps()
{
	for (int n = 0; n < num; ++n) delete kcomp[n];
	delete sname;
	delete[] kcomp;
	delete[] nmidx;
	delete[] dist;
}

double qdiv(Kcomp* a, Kcomp* b)
{
static	double	param[4][2] = {
{0.92042, 0.18677},	// d314
{0.34576, 0.07108},	// d414
{0.22333, 0.03164},	// d514
{0.18704, 0.00501}	// d420
};
	double	f = 0;
	if (a->kmers && b->kmers) { 
	    int	s = 0;
	    KVpair<INT, int>*	ka = a->kmers;
	    KVpair<INT, int>*	kb = b->kmers;
	    KVpair<INT, int>*	ta = a->kmers + a->UnqKmers;
	    KVpair<INT, int>*	tb = b->kmers + b->UnqKmers;

	    while (ka < ta && kb < tb) {
		if (ka->key == kb->key) {
		    s += min(ka->val * b->many, kb->val * a->many);
		    ++ka; ++kb;
		} else if (ka->key < kb->key) ++ka;
		else	++kb;
	    }
	    f = (double) a->TtlKmers / a->many;
	    double e = (double) b->TtlKmers / b->many;
	    f = min(f, e) * a->many * b->many;
	    f = (double) s / f;
	} 
	double d = 1. - f;
	switch (wcp.distm) {
	    case QFdiv: return (d);
	    case QFidn: return (f);
	    case NJuCa:	return jukcan(d);
	    default: break;
	}
	int	s = 3;
	if (wcp.Nalpha == 14) {
	    switch (wcp.Ktuple) {
		case 3: s = 0; break;
		case 4: s = 1; break;
		case 5: s = 2; break;
		default: s = 3; break;
	    }
	}
	f = param[s][0] * log((param[s][1] + f)/(param[s][1] + 1.)) + 1.;
	d = 1. - f;
	switch (wcp.distm) {
	    case QCidn: return(f);
	    case QJuCa: return jukcanpro(d);
	    case Qpamd: return pamcorrect(d) / 100.;
	    default: break;
	}
	return (d);
}

int kcmp_main(CalcServer<Kcomp>* svr, Kcomp* kcmp[], ThQueue<Kcomp>* q = 0)
{
	Kcomp*&	a = kcmp[0];
	Kcomp*&	b = kcmp[1];
	Kcomps*	kcp = (Kcomps*) svr->prm;
	FTYPE	dst = 100. * qdiv(a, b);
	int	nn = svr->calcnbr(a->sid, b->sid);
	kcp->dist[nn] = dst;
	return (OK);
}

static FTYPE* calcdist_kmer(int inm, Kcomps* kcomps)
{
	CalcServer<Kcomp> svr(inm, kcomps, &kcmp_main, 
		0, 0, kcomps->kcomp, kcomps->num);
	svr.auto_comp();
	FTYPE*	dist = kcomps->getdist();
	return (dist);
}

FTYPE* calcdist_kmer(Seq* sd, int inm, Subset* ss)
{
	resetqdiv(sd);
	Kcomps*	kcomps = new Kcomps(sd, inm, ss);
	FTYPE*	dist = calcdist_kmer(inm, kcomps);
	delete kcomps;
	return (dist);
}

FTYPE* calcdist_kmer(Seq* sd, int argc, const char** argv, 
	const char* catalog, int num, int inm)
{
	resetqdiv(sd);
	Kcomps*	kcomps = new Kcomps(sd, argc, argv, catalog, num, inm);
	FTYPE*	dist = calcdist_kmer(inm, kcomps);
	delete kcomps;
	return (dist);
}

FTYPE* calcdist_kmer(Seq** seqs, int nn, int inm, Strlist*& snm)
{
	resetqdiv(seqs[0]);
	Kcomps*	kcomps = new Kcomps(seqs, nn, inm);
	FTYPE*	dist = calcdist_kmer(inm, kcomps);
	snm = kcomps->sname;
	kcomps->sname = 0;
	delete kcomps;
	return (dist);
}

/*************************************************************************
*
*	main function of qdiv
*
*************************************************************************/

#ifdef MAIN

void usage()
{
	fputs("Usage: qdiv [-i:e|j catalog] [Options] [Seq1..SeqN]\n", stderr);
	fputs("\t[-A AaClass] [-B BitPat] [-k tuple] [-K D|P] [-s dir]\n", stderr);
	fputs("\t[-ON: 0:rawD; 1:rawI; 2:fitD; 3:fitI; 4:JuCa; 5:Pam]\n", stderr);
}

int kcmp_out(CalcServer<Kcomp>* svr, Kcomp* kcmp[], ThQueue<Kcomp>* q = 0)
{
	Kcomp*&	a = kcmp[0];
	Kcomp*&	b = kcmp[1];
	Kcomps*	kcp = (Kcomps*) svr->prm;
	FTYPE	dst = kcp->dist[svr->calcnbr(a->sid, b->sid)];
	int	i = kcp->nmidx[a->sid];
	int	j = kcp->nmidx[b->sid];
	if (dst < kcp->vthr) {
	    fprintf(out_fd, "%7.2f\t%d\t%d\t%s\t%s\n", (float) dst, 
		kcp->kcomp[a->sid]->TtlKmers / kcp->kcomp[a->sid]->many,
		kcp->kcomp[b->sid]->TtlKmers / kcp->kcomp[b->sid]->many,
		(*kcp->sname)[i], (*kcp->sname)[j]);
	}
	return (OK);
}

template <class Kcomp>
int AlnServer<Kcomp>::localoption(int& argc, const char**& argv)
{
	int	n = 1;

	switch (argv[0][1]) {
	  case 'A':	// aa classification pattern
	    wcp.reducap = getarg(argc, argv); break;
	  case 'B':	// matching bit pattern
	    wcp.sBitPat = getarg(argc, argv); break;
	  case 'O':	// output mode
	    wcp.distm = (DistMes) atoi(getarg(argc, argv)); break;
	  case 'R':
	    wcp.Nalpha = (INT) atoi(getarg(argc, argv)); break;
	  case 'k':	// tuple size
	    wcp.Ktuple = atoi(getarg(argc, argv)); break;
	  default: n = 0; break;
	}
	return (n);
}

int	main(int argc, const char** argv)
{
	alprm.thr = 0;
	InputSeqTest	tv = zero_InputSeqTest;
	AlnServer<Seq> ssvr(argc, argv, IM_SNGL, IM_ALTR, (void*) &tv, &testInputSeq);
	ssvr.autocomp(false);
	if (tv.bad) fatal("%d seqs were not found !\n", tv.bad);
	if (tv.num < 2) fatal("not enough data !\n");
	if (tv.molc == PROTEIN && wcp.Nalpha == 4) resetqdiv(ssvr.in_face[0]);
	if (ssvr.calc_mode == IM_SNGL) ssvr.calc_mode = IM_ALTR;
	Kcomps	kcomps(ssvr.in_face[0], argc, argv, tv.num, &ssvr);
	CalcServer<Kcomp> ksvr(ssvr.calc_mode, &kcomps, 
	    &kcmp_main, 0, 0, kcomps.kcomp, kcomps.num, 
	    ssvr.calc_mode == IM_GRUP? ssvr.getgrp2(): 0);
	ksvr.autocomp(true);
	ksvr.change_job(&kcmp_out);
	ksvr.autocomp(false);
	eraStrPhrases();
	return (0);
}

#endif // MAIN
