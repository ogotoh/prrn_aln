/*****************************************************************************
*
*	adjmat.c: calculate adjacent matrix from a set of sequences
*
*	Mapping and alignment of protein/cDNA sequences onto
*	genomic sequence
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

#include "adjmat.h"
#include "divseq.h"
#include "qdiv.h"
#include <math.h>

typedef	KVpair<INT, int>	kv_t;

static	int	MinQueryLen = 0;
static	FTYPE	max_dist = FLT_MAX;
static	FTYPE	max_distB = FLT_MAX;

#if M_THREAD
static	void*	worker_func(void* arg);
static	void*	master_func(void* targ);

class AdjMatThQueue {
	int     rp, wp;
	int     remain;
	int     done;
	pthread_cond_t  not_full;
	pthread_cond_t  not_empty;
public:
	Seq**   sinp;
	Seq**   sque;
	pthread_mutex_t mutex;
	AdjMatThQueue(Seq** sqs);
	~AdjMatThQueue() {delete[] sque;}
	void    enqueue(Seq** fsd);
	void    dequeue(Seq** fsd);
};

struct adj_thread_arg_t {
	AdjMatThQueue*  q;
	int     cpuid;
	Seq**   seqs;
	SrchBlk*   sbk;
	AdjacentMat*    ajm;
};

struct adj_mast_arg_t {
	AdjMatThQueue*  q;
	DbsDt*	  dbf;
};
#endif	// M_THREAD

void	set_max_dist(FTYPE val) {
	std::swap(max_dist, max_distB);
	std::swap(max_dist, val);
}

void	reset_max_dist() {
	std::swap(max_dist, max_distB);
}

FTYPE AdjacentMat::dist_align2(Seq* sqs[], const PwdB* pwd)
{
	Gsinfo	GsI;
	sqs[0]->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
	sqs[1]->exg_seq(algmode.lcl & 4, algmode.lcl & 8);
	GsI.skl = alignB_ng((const Seq**) sqs, pwd, &GsI.scr);
	if (!GsI.skl || GsI.skl->n == 0) return (NEVSEL);
	GsI.scr = skl_rngB_ng((const Seq**) sqs, &GsI, pwd);
	return (degdiv(&GsI));
}

FTYPE AdjacentMat::dist_kmer(Seq* sqs[])
{
	Seq*&   a = sqs[0];
	Seq*&	b = sqs[1];
	Kcomp	ka(a, a->right - a->left, rdword->iConvTab, 0);
	Kcomp	kb(b, b->right - b->left, rdword->iConvTab, 1);
	return ((FTYPE) qdiv(&ka, &kb));
}

#if SPNRANK
void AdjacentMat::fillmat(Seq* sqs[], FTYPE dist, int rnk)
{
	Edge	e = {sqs[0]->did, sqs[1]->did, dist, rnk};
	if (e.dist < max_dist) {
	    if (e.u > e.v) std::swap(e.u, e.v);
	    vedges.push_back(e);
	}
}
#else
void AdjacentMat::fillmat(Seq* sqs[], FTYPE dist)
{
	Edge	e = {(INT) sqs[0]->did, (INT) sqs[1]->did, dist};
	if (e.dist < max_dist) {
	    if (e.u > e.v) std::swap(e.u, e.v);
	    vedges.push_back(e);
	}
}
#endif	// SPNRANK

void AdjacentMat::spaln_job(Seq* sqs[], SrchBlk* bks, AdjMatThQueue* q)
{
	int	nparalog = bks->finds(sqs);
	if (nparalog == ERROR) return;		// no alignment
	int	n_out = nparalog > 1? OutPrm.supself: 0;
	for (int n = 1; n <= nparalog && n_out < nparalog; ++n) {
	    if (n > 1) std::swap(sqs[1], sqs[n]);
	    if (sqs[0]->did == sqs[1]->did || 
		(OutPrm.supself > 1 && sqs[0]->did > sqs[1]->did)) continue;
	    FTYPE	dist = 0;
	    switch (distcal_mode) {
		case Composition: dist = dist_kmer(sqs); break;
		case DynAln:	dist = dist_align2(sqs, bks->pwd); break;
		case DynScr: default:
			dist = alnscore2dist(sqs, bks->pwd); break;
	    }
	    if ((dist *= 100.) < alprm.thr) {
#if M_THREAD
		if (q) {
		    pthread_mutex_lock(&q->mutex);
#if SPNRANK
		    fillmat(sqs, dist, n);
#else
		    fillmat(sqs, dist);
#endif
		    pthread_mutex_unlock(&q->mutex);
		} else
#endif	// M_THREAD
#if SPNRANK
		    fillmat(sqs, dist, n);
#else
		    fillmat(sqs, dist);
#endif	// SPNRANK
		 ++n_out;
	    }
	}
	if (nparalog > 1) std::swap(sqs[1], sqs[2]);
}

SrchBlk* AdjacentMat::getblkinf(Seq* sqs[], MakeBlk* mb)
{
	Seq*&	a = sqs[0];	// query
	Seq*&	b = sqs[1];	// database

	a->inex.intr = b->inex.intr = 0;
	SrchBlk* bks = new SrchBlk(sqs, mb);
	bks->setseqs(sqs);
	if (MinQueryLen == 0) MinQueryLen = bks->MinQuery();
	return bks;
}

void AdjacentMat::all_in_func(Seq** sqs, SrchBlk* sbk)
{
	sbk->dbf = target_dbf;
	Seq*&	a = sqs[0];
	for (INT i = 0; i < target_dbf->numidx; ) {
	    if (a->right - a->left < MinQueryLen) {
		prompt("%s (%d - %d) %d is too short!\n",
                a->sqname(), a->left, a->right, a->right - a->left);
	    } else
// perform the job
	        spaln_job(sqs, sbk);
// read query sequence
	    if (a->read_dbseq(target_dbf, ++i) == 0) break;
	}
}

#if M_THREAD
AdjMatThQueue::AdjMatThQueue(Seq** sqs) : sinp(sqs)
{
	sque = new Seq*[max_queue_num];
	initseq(sque, max_queue_num);
	rp = wp = remain = 0;
	pthread_mutex_init(&mutex, 0);
	pthread_cond_init(&not_full, 0);
	pthread_cond_init(&not_empty, 0);
}

void AdjMatThQueue::enqueue(Seq** sqs)
{
	pthread_mutex_lock(&mutex);
	while (remain == max_queue_num)
	    pthread_cond_wait(&not_full, &mutex);
	if (sqs) {
#if QDEBUG
	    fprintf(stderr, "e%d: %s %d %d %d\n", n, sqs[0]->sqname(),
	    sqs[0]->sid, sqs[0]->len, sqs[0]->many);
#endif
	    if (sqs[0]->sid > 0) swapseq(sque + wp, sqs);
	} else	sque[wp]->refresh(0);
	++wp; ++remain;
	if (wp == max_queue_num) wp = 0;
	pthread_cond_signal(&not_empty);
	pthread_mutex_unlock(&mutex);
}

void AdjMatThQueue::dequeue(Seq** sqs)
{
	pthread_mutex_lock(&mutex);
	while (remain == 0)
	     pthread_cond_wait(&not_empty, &mutex);
#if QDEBUG
	if (sque[rp]->many)	{swapseq(sqs, sque + rp);
	fprintf(stderr, "d%d: %s %d %d %d\n", 0, sqs[0]->sqname(),
	sqs[0]->sid, sqs[0]->len, sqs[0]->many);}
#else
	if (sque[rp]->many)	swapseq(sqs, sque + rp);
#endif
	else	sqs[0]->refresh(0);
	++rp; --remain;
	if (rp == max_queue_num) rp = 0;
	pthread_cond_signal(&not_full);
	pthread_mutex_unlock(&mutex);
}

static	void* master_func(void* arg)
{
	adj_mast_arg_t*	targ = (adj_mast_arg_t*) arg;
	AdjMatThQueue*	q = targ->q;
	DbsDt*	dbf = targ->dbf;
	Seq**	sqs = q->sinp;
	Seq*&	a = sqs[0];

	for (INT i = 0; i < dbf->numidx; ++i) {
// read query and push it in the queue
	    if (a->read_dbseq(dbf, i) == 0) break;
	    if (a->right - a->left < MinQueryLen) {
		prompt("%s (%d - %d) %d is too short!\n",
                a->sqname(), a->left, a->right, a->right - a->left);
	    } else	q->enqueue(sqs);
	}
	for (int n = 0; n < thread_num; ++n)
	    q->enqueue(0);
	return (void*) 0;
}

static	void* worker_func(void* arg)
{
	adj_thread_arg_t*	targ = (adj_thread_arg_t*) arg;

#ifdef __CPU_SET
	cpu_set_t	mask;
	__CPU_ZERO(&mask);
	__CPU_SET(targ->cpuid, &mask);
	if (sched_setaffinity(0, sizeof(mask), &mask) == -1)
	    prompt("Warning: faild to set CPU affinity !\n");
#endif

	while (true) {
	    targ->q->dequeue(targ->seqs);
	    if (targ->seqs[0]->many == 0) break;
	    targ->ajm->spaln_job(targ->seqs, targ->sbk, targ->q);
	}
	return (void*) 0;
}

void AdjacentMat::MasterWorker(Seq** sqs, SrchBlk* primaty)
{
	adj_mast_arg_t	maarg;
	pthread_t	master;
	int	cpu_num = sysconf(_SC_NPROCESSORS_CONF);

	if (thread_num <= 0) thread_num = cpu_num;
	if (max_queue_num == 0) max_queue_num = int(FACT_QUEUE * thread_num);
	max_queue_num = max_queue_num;
	adj_thread_arg_t*	targ = new adj_thread_arg_t[thread_num];
	pthread_t*	worker = new pthread_t[thread_num];

	AdjMatThQueue	q(sqs);
	maarg.q = &q;
	maarg.dbf = target_dbf;
	targ[0].seqs = new Seq*[no_seqs * thread_num];
	initseq(targ[0].seqs, no_seqs * thread_num);
	for (int n = 0; n < thread_num; ++n) {
	    targ[n].q = &q;
	    targ[n].cpuid = n % cpu_num;
	    if (n > 0) targ[n].seqs = targ[n - 1].seqs + no_seqs;
	    targ[n].seqs[0]->inex = sqs[0]->inex;
	    targ[n].seqs[1]->inex = sqs[1]->inex;
	    targ[n].ajm = this;
	    if (n) {
		targ[n].sbk = new SrchBlk(primaty, target_dbf);
	    } else {
		primaty->reset(target_dbf);
		targ[n].sbk = primaty;
	    }
	}
	q.enqueue(sqs);
	pthread_create(&master, 0, master_func, (void*) &maarg);
	for (int n = 0; n < thread_num; ++n)
	    pthread_create(worker + n, 0, worker_func, (void*) (targ + n));
	for (int n = 0; n < thread_num; ++n)
	    pthread_join(worker[n], 0);
	clearseq(targ[0].seqs, no_seqs * thread_num);
	clearseq(q.sque, max_queue_num);
	for (int n = 1; n < thread_num; ++n) {
	    targ[n].sbk->dbf = 0;
	    delete targ[n].sbk;
	}
	delete[] targ[0].seqs;
	delete[] targ;
	delete[] worker;
}
#endif	// M_THREAD

AdjacentMat::AdjacentMat(const char* fn, FTYPE distthr, int clmn, bool sim)
	: no_seqs(0), namebuf(0), singleton(0), rdword(0), 
	  distcal_mode((DistCal) algmode.any), target_dbf(0), 
	  n_vertex(0), n_edge(0), mname(0), verdid(0), edges(0)
{
	if (distPrm.contbt_spj > 0) distcal_mode = DynAln;
	char    str[MAXL];

	FILE*   fd = fopen(fn, "r");
	if (!fd) fatal("%s not found !n", fn);
	n_edge = n_vertex = 0;
	while (fgets(str, MAXL, fd)) {
	    if (*str == '>') {	// FASTA format
		n_edge = 0;
		return;
	    }
	    ++n_edge;  // number of edges
	}
	rewind(fd);
	int	maxver = 2 * n_edge;
	StrHash<int>	verhash(maxver);
	int     nerr = 0;
	float	maxdist = 0;
	while (fgets(str, MAXL, fd)) {	  // hash table
	    Strlist	stl(str, stddelim);
	    int		nf = stl.size();
	    if (nf - clmn < 3) continue;
	    FTYPE       dst = atof(stl[clmn - 1]);
	    if ((dst > distthr) ^ sim) continue;
	    if (dst > maxdist) maxdist = dst;
	    char*	as = stl[nf - 2];
	    char*	bs = stl[nf - 1];
	    if (strcmp(as, bs)) {
		verhash.pile(as);
		verhash.pile(bs);
	    } else if (dst != 0) {
		prompt("%s\n", as);
		++nerr;
	    }
	}
	if (nerr) fatal("%d same id but non-zero distance !\n", nerr);
	n_vertex = verhash.size() - 1;
	mname = new const char*[n_vertex];
	kv_t*	kv = verhash.begin();
	kv_t*	kz = verhash.end();
	for ( ; kv < kz; ++kv) {
	    if (kv->val) {	      // str to number
		mname[kv->val - 1] = verhash.strkey(kv->key);
	    }
	}
	rewind(fd);
	edges = new Edge[n_edge];
	Edge*   e = edges;	      // make adjacent matrix
	Dhash<INT, int> edghash(n_edge);
	while (fgets(str, MAXL, fd)) {
	    Strlist	stl(str, stddelim);
	    int		nf = stl.size();
	    if (nf - clmn < 3) continue;
	    FTYPE       dst = atof(stl[clmn - 1]);
	    if ((dst > distthr) ^ sim) continue;
	    e->dist = sim? maxdist - dst: dst;
	    char*	as = stl[nf - 2];
	    char*	bs = stl[nf - 1];
	    kv = verhash.find(as);
	    e->u = kv->val - 1;
	    kv = verhash.find(bs);
	    e->v = kv->val - 1;
	    if (e->u > e->v) std::swap(e->u, e->v);
	    KVpair<INT, int>*   edgkv = edghash.incr(n_vertex * e->u + e->v);
	    if (edgkv->val == 1) ++e;
	}
	n_edge = e - edges;
	namebuf = verhash.squeeze();
	fclose(fd);
}

AdjacentMat::AdjacentMat(int argc, const char** argv, 
	int molc, const char* catalog, int min_seqs)
	: no_seqs(0), namebuf(0), singleton(0), rdword(0), 
	  distcal_mode((DistCal) algmode.any), target_dbf(0), 
	  n_vertex(0), n_edge(0), mname(0), verdid(0), edges(0)
{
	if (distPrm.contbt_spj > 0) distcal_mode = DynAln;
	if (distcal_mode == Composition) {
	    Seq	sd(1);
	    setSeqCode(&sd, molc);
	    WCPRM*	wcp = resetqdiv(&sd);
	    rdword = new ReducWord(&sd, wcp->Nalpha, wcp->reducap);
	}
	MakeBlk*	mb = 0;
	if (alprm2.spb > 0.) {
	    SeqServer	svr(argc, argv, IM_SNGL, catalog, molc);
	    mb = makeblock(&svr);
	} else
	    mb = makeblock(argc, argv, molc);
	prePwd(molc);
	if (mb->no_entry() >= min_seqs) calcAdjMat(mb);
	delete mb;
}

AdjacentMat::AdjacentMat(Seq* sd)
	: no_seqs(0), namebuf(0), singleton(0), rdword(0), 
	  distcal_mode((DistCal) algmode.any), target_dbf(0), 
	  n_vertex(0), n_edge(0), mname(0), verdid(0), edges(0)
{
	if (distPrm.contbt_spj > 0) distcal_mode = DynAln;
	if (distcal_mode == Composition) {
	    WCPRM*	wcp = resetqdiv(sd);
	    rdword = new ReducWord(sd, wcp->Nalpha, wcp->reducap);
	}
	MakeBlk*	mb = makeblock(sd);
	prePwd(sd->inex.molc);
	calcAdjMat(mb);
	delete mb;
}

void AdjacentMat::calcAdjMat(MakeBlk* mb)
{
	target_dbf = mb->wdbf;
	if ((OutPrm.MaxOut > target_dbf->numidx - 1) || (OutPrm.MaxOut == 0))
	    OutPrm.MaxOut = int(log(double(target_dbf->numidx) + 1));
	no_seqs = OutPrm.MaxOut + 2;
	Seq**	seqs = new Seq*[no_seqs];
	initseq(seqs, no_seqs);	// 0: query, 1: database, 2: reverse
	if (seqs[0]->read_dbseq(target_dbf, 0) == 0)
	    fatal("Can't open query !\n");
	SrchBlk*	bprm = getblkinf(seqs, mb);
	INT	maxver = target_dbf->numidx;

#if M_THREAD
	if (thread_num) MasterWorker(seqs, bprm); else
#endif
	all_in_func(seqs, bprm);
	bprm->dbf = 0;
	delete bprm;

	n_vertex = 0;
	n_edge = vedges.size();		// tentative
	if (n_edge) {
	  verhash.resize(2 * n_edge);
	  Dhash<LONG, int> edghash(n_edge);
	  n_edge = 0;
	  for (EdgeItr ei = vedges.begin(); ei < vedges.end(); ++ei) {
	    KVpair<LONG, int>*   edgkv = edghash.incr(ei->u * maxver + ei->v);
	    if (edgkv->val == 1) {
		++n_edge;
		verhash.incr(ei->u);
		verhash.incr(ei->v);
	    }
	  }
	  for (KVpair<INT, int>* verkv = verhash.begin(); verkv < verhash.end(); ++verkv)
	    if (verkv->val) ++n_vertex;
	  mname = new const char*[n_vertex];
	  verdid = new int[n_vertex];

	  n_vertex = 0;
	  for (KVpair<INT, int>* verkv = verhash.begin(); verkv < verhash.end(); ++verkv) {
	    if (verkv->val) {
		mname[n_vertex] = target_dbf->entname(verkv->key);
		verdid[n_vertex] = verkv->key;
		verhash.assign(verkv->key, n_vertex++ + 1);
	    }
	  }
// remove redundancy
	  edges = new Edge[n_edge];
	  Edge*	eg = edges;
	  edghash.clear();
	  for (EdgeItr ei = vedges.begin(); ei < vedges.end(); ++ei) {
	    KVpair<LONG, int>*   edgkv = edghash.incr(ei->u * maxver + ei->v);
	    if (edgkv->val == 1) {
		KVpair<INT, int>* verkv = verhash.find(ei->u);
		*eg = *ei;
		eg->u = verkv->val - 1;
		verkv = verhash.find(ei->v);
		eg->v = verkv->val - 1;
		if (eg->u > eg->v) std::swap(eg->u, eg->v);
		++eg;
	    }
	  }
	  n_edge = eg - edges;
	} else {
	    mname = 0; verdid = 0; edges = 0;
	}
	clearseq(seqs, no_seqs);
	delete[] seqs;
}

int* AdjacentMat::outsiders()
{
	if (singleton) return singleton;
	singleton = new int[target_dbf->numidx - n_vertex + 1];
	int*	wls = singleton;
	KVpair<INT, int>	zerkv = {0, 0};
	for (INT i = 0; i < target_dbf->numidx; ++i) {
	    KVpair<INT, int>* verkv = n_vertex? verhash.find(i): &zerkv;
	    if (!verkv || verkv->val == 0) *wls++ = i;
	}
	*wls = -1;
	return (singleton);
}

void AdjacentMat::print_singletons(const char* fn)
{
	FILE*	fd = fn? fopen(fn, "w"): 0;
	if  (!fd) {
	    prompt("fail to write to %s !\n", fn);
	    return;
	}
	int*	wls = outsiders();
	while (*wls >= 0) 
	    fprintf(fd, "%s\n", target_dbf->entname(*wls++));
	fclose(fd);
}


