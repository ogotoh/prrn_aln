/****************************************************************************
*
*	adjmat.h: Header to adjacent matrix calclation
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
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japang
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#define	SPNRANK	0

#include "aln.h"
#include "mseq.h"
#include "autocomp.h"
#include "wln.h"
#include "blksrc.h"
#include <vector>
#include <algorithm>
#include <functional>

#ifdef M_THREAD
#include <pthread.h>
#include <unistd.h>
#include <sched.h>
#endif

extern	int	thread_num;
static	const	int	DefPam = 150;
static	const	int	WlpPam = 50;
static	const	char*	singleton_fn = "singleton.lst";

class	AdjMatThQueue;

struct Edge {
	INT	u;
	INT	v;
	FTYPE	dist;
#if SPNRANK
	INT	rank;
#endif
	bool	operator<(const Edge& right) const {
	    return (dist < right.dist);
	}
	bool	operator<=(const Edge& right) const {
	    return (dist <= right.dist);
	}
};

class AdjacentMat {
private:
	int	no_seqs;
	char*	namebuf;
	int*	singleton;
	ReducWord*	rdword;
	DistCal	distcal_mode;
#if M_THREAD
	void	MasterWorker(Seq** sqs, SrchBlk* prm);
#endif
	FTYPE	dist_align2(Seq* sqs[], PwdB* pwd);
	FTYPE	dist_kmer(Seq* sqs[]);
#if SPNRANK
	void	fillmat(Seq* sqs[], FTYPE dist, int rnk);
#else
	void	fillmat(Seq* sqs[], FTYPE dist);
#endif
	SrchBlk*	getblkinf(Seq* sqs[], MakeBlk* mb);
	void	all_in_func(Seq** sqs, SrchBlk* sbk);
	void	calcAdjMat(MakeBlk* mb);
protected:
	DbsDt*	target_dbf;
	int	n_vertex;
	int	n_edge;
	const	char**	mname;
	int*	verdid;
	Edge*	edges;
	std::vector<Edge>	vedges;
	Dhash<INT, int> verhash;
public:
	void	spaln_job(Seq* sqs[], SrchBlk* bks, AdjMatThQueue* q = 0);
	AdjacentMat(const char* fn, FTYPE distthr = 100.);
	AdjacentMat(int argc, const char** argv, int molc, 
	    const char* catalog = 0, int min_seqs = 0);
	AdjacentMat(Seq* sd);
	~AdjacentMat() {
	    delete[] mname; delete[] edges; delete[] singleton;
	    delete[] verdid;  delete[] namebuf; delete rdword;
	    delete target_dbf;
	}
	int*	outsiders();
	DbsDt*	dbf() {return (target_dbf);}
	int	no_outsiders() {return (target_dbf->numidx - n_vertex);}
	int	no_edges() {return (n_edge);}
	int	no_verteces() {return (n_vertex);}
	void	print_singletons(const char* fn = singleton_fn);
};

typedef std::vector<Edge>::iterator EdgeItr;

extern	void	set_max_dist(FTYPE val);
extern	void	reset_max_dist();

