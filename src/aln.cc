/*****************************************************************************
*
*	Alignment of two protein/nucleotide sequences or
*	pre-aligned groups of sequences 
*
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

#include "aln.h"
#include "vmf.h"
#include "wln.h"
#include "mseq.h"
#include "maln.h"
#include "mgaps.h"
#include "autocomp.h"
#include "consreg.h"
#include "fspscore.h"

static	const	VTYPE	JobError = INT_MIN;
static	const	int	WlpPam = 50;
static	const	char*	TREE = "tree";
static	const	char*	Def_out_file = DEF_OUT_FN;
static	const	char*	genomedb = 0;
static	const	char*	guidetree = 0;

void usage()
{
	fputs("Specific Options:\n", stderr);
	fputs("\t-a catalog; pileupt msa\n", stderr);
	fputs("\t-b tree; progressive msa\n", stderr);
	fputs("\t-ix catalog; input mode\n", stderr);
	fputs("\t\tx? 'a': alternate; 'e'; every pair; 'f' first vs others;\n", stderr);
	fputs("\t\tx? 'g': between groups; 'i': self; 'p' parallel;\n", stderr);
	fputs("\t-L \t;Semi-local\n", stderr);
	fputs("\t-M \t;Both directions\n", stderr);
	fputs("\t-r \t;Don\'t remove temp files\n", stderr);
	fputs("\t-u#\t;Gap ext. weight\n", stderr); 
	fputs("\t-v#\t;Gap open weight\n", stderr);
	fputs("\t-w#\t;Band width around main diagonals\n", stderr);
	fputs("\t-yx#: splicing options \n", stderr);
}

template <class seq_t>
void AlnServer<seq_t>::setparam(int level)
{
	if (level & 1) {
	    setalprm();
	    setstrip(QUERY);
	}
	if (level & 2) {
	    setlpw(QUERY);
	    setlsegs(QUERY);
	    setprmode(QUERY, QUERY, QUERY);
	    setalgmode(QUERY, QUERY);
	    setminus(QUERY);
	    setshuffle(SILENT, QUERY, QUERY);
	}
}

template <class seq_t>
int AlnServer<seq_t>::localoption(int& argc, const char**& argv)
{
	const	char*	opt = argv[0] + 1;
	const	char*	val = argv[0] + 2;
	int	n = 1;

	switch (*opt) {
	  case 'a': 
	    this->calc_mode = IM_ADON;
	    if ((val = getarg(argc, argv)))
		this->catalog = val;
	    else	this->catalog = CATALOG;
	    break;
	  case 'b': 
	    this->calc_mode = IM_TREE;
	    if ((val = getarg(argc, argv)))
		guidetree = val;
	    else	guidetree = TREE;
	    break;
	  case 'G':
	    genomedb = *val? val: "";
	    algmode.lcl = 15;
	    algmode.lsg = 1;
	    break;
	  case 'L': 
	    val = getarg(argc, argv, true);
	    if (val && isdigit((*val)))
		algmode.lcl = atoi(val);
	    else {
		switch (tolower(opt[1])) {
		  case 'l': algmode.lcl = 5; break;
		  case 'r': algmode.lcl = 10; break;
		  case 'd': algmode.lcl = 3; break;
		  case 'u': algmode.lcl = 12; break;
		  case 's': algmode.lcl = 16; break;	// SWG
		  default:  algmode.lcl = 15; break;
		}
	    }
	    break;
	  case 'M': 
	    if ((val = getarg(argc, argv, true)))
			OutPrm.NoOut = atoi(val);
	    else	OutPrm.NoOut = Def_COLONY;
	    if (OutPrm.NoOut < 4) {
		algmode.mlt = OutPrm.NoOut;
		if (OutPrm.NoOut != 1) OutPrm.NoOut = MAX_COLONY;
	    } else	algmode.mlt = 2;
	    break;
	  case 'r':
	    OutPrm.RemoveTmp = false;
	    break;
	  case 'U': algmode.aut = 0; break;
	  case 'V':
	    if ((val = getarg(argc, argv, true))) {
		long	m = long(ktof(val));
		setVmfSpace(m);
	    }
	    break;
	  default: n = 0; break;
	}
	return (n);
}

/********************************************************************
*
*	Construct multiple sequence alignment by a progressive method
*
********************************************************************/

template <class seq_t>
int AlnServer<seq_t>::makemsa()
{
	Btree<Tnode>*	ltree = 0;
	int	rv = this->interactive? 1: 2;

	if (logfile && !(logfd = fopen(logfile, "w"))) fatal(erropen, logfile);
	switch (this->calc_mode) {
	  case IM_ADON:
	    MakeMsa<seq_t>(this);
	    rv = 0;
	    break;
	  case IM_TREE:
	    if (!guidetree) fatal("Give guide tree !\n");
	    ltree = new Btree<Tnode>(guidetree);
	    if (!ltree->root) fatal(not_found, guidetree);
	    ltree->fill_tname();
	    MakeMsa<seq_t>(this, ltree->root);
	    delete ltree;
	    rv = 0;
	    break;
	}
	if (logfd) fclose(logfd);
	return (rv);
}

template <class seq_t>
VTYPE valonly(CalcServer<seq_t>* svr, seq_t** sqs, PwdM* pwd, ThQueue<seq_t>* q)
{
	long	rr[2];
	seq_t*&	a = sqs[0];
	seq_t*&	b = sqs[1];
	seq_t*&	shorter = (a->right - a->left) < (b->right - b->left) ? a: b;
	VTYPE	scr = JobError;

	int 	sz = shorter->right - shorter->left;
	if (sz) {
	    scr = HomScore(sqs, pwd, rr);
	    VTYPE	ref = selfAlnScr(shorter, pwd->simmtx);
	    m_thread_Lock(q);
	    fprintf(out_fd, "%6.1f\t%6.1f\t%6.1f\t", 
		(float) scr, (float) ref, 100. * (scr - ref) / sz);
	    fprint_seq_mem((const Seq**) sqs, 2, out_fd);
	    m_thread_Unlock(q);
	}
	return (scr);
}

template <class seq_t>
VTYPE alnoutput(CalcServer<seq_t>* svr, seq_t** sqs, PwdM* pwd, Gsinfo* GsI, 
	ThQueue<seq_t>* q)
{
	seq_t*&	a = sqs[0];
	seq_t*&	b = sqs[1];
	if (pwd->swp) {
	    if (GsI && GsI->sigII && GsI->sigII->lst)
		GsI->sigII->swaplst(a->many, b->many);
	    swapskl(GsI->skl);
	    swap(a, b);
	    pwd->swapDvsP();
	    pwd->swp = false;
	}
	m_thread_Lock(q);
	int	omode = algmode.nsa & 7;
	Seq*	gene = 0;
	int	print2Skip = 0;
	GAPS*	gaps[2] = {0, 0};
	if (a->inex.intr)	gene = a;
	else if (b->inex.intr)	gene = b;
	if (omode == 0 || omode == 4 || (omode == 5 && gene)) setlpw(0);

	if (pwd->DvsP == 1 || pwd->DvsP == 2) {	// DNA vs protein
	    if (omode == 1) {
		if (gene) {
		    fphseqs((const Seq**) sqs, 2, out_fd);
		    GBcdsForm(GsI->CDSrng, gene, out_fd);
		    print2Skip = 1;
		}
		if (getlpw()) {		/* print alignment */
		    skl2gaps3(gaps, GsI->skl, pwd->DvsP);
		    if (gene) gene->exons = GsI->eiscrunfold(gaps[2 - pwd->DvsP]);
		    unfoldgap(gaps[0], pwd->DvsP == 1? 3: 1);
		    unfoldgap(gaps[1], pwd->DvsP == 2? 3: 1);
		    GsI->print2((Seq**) sqs, (const GAPS**) gaps, 1, 1, print2Skip, out_fd);
		} else
		    GsI->repalninf((Seq**) sqs, omode, out_fd);
	    } else if (gene) {
		GsI->printgene((Seq**) sqs, omode, out_fd);
	    }
	} else {			/* DNA vs DNA or Pro vs Pro */
	    if (omode == 1 || (omode == 5 && !gene)) {
		if (getlpw()) {
		    skl2gaps(gaps, GsI->skl);
		    if (gene) gene->exons = GsI->eiscrunfold(gaps[0]);
		    toimage(gaps, 2);
		}
		if (gene) {
		    fphseqs((const Seq**) sqs, 2, out_fd);
		    GBcdsForm(GsI->CDSrng, gene, out_fd);
		    print2Skip = 1;
		}
		if (getlpw())
		    GsI->print2((Seq**) sqs, (const GAPS**) gaps, 1, 1, print2Skip, out_fd);
		else
		    GsI->repalninf((Seq**) sqs, omode, out_fd);
	    } else if (gene) {
		GsI->printgene((Seq**) sqs, omode, out_fd);
	    } else if (omode == 0 && !OutPrm.ColorEij) {
		svr->fpavsd(GsI->scr);
		put_stat(out_fd, GsI);
		fprint_seq_mem((const Seq**) sqs, 2, out_fd);
	    } else {
		GsI->repalninf((Seq**) sqs, omode, out_fd);
	    }
	}
	if (gene) {delete[] gene->exons; gene->exons = 0;}
	delete[] gaps[0];
	delete[] gaps[1];
	m_thread_Unlock(q);
	return (GsI->scr);
}

template <class seq_t>
void swgscore(CalcServer<seq_t>* svr, seq_t** sqs, VTYPE val, ThQueue<seq_t>* q = 0)
{
	m_thread_Lock(q);
	svr->fpavsd(val);
	fprintf(out_fd, "%6.2lf ", (double) val / sqs[0]->many / sqs[1]->many);
	fprint_seq_mem((const Seq**) sqs, 2, out_fd);
	m_thread_Unlock(q);
}

template <class seq_t>
VTYPE match_2(CalcServer<seq_t>* svr, seq_t* sqs[], PwdM* pwd, ThQueue<seq_t>* q)
{
	seq_t*&	a = sqs[0];
	seq_t*&	b = sqs[1];
	VTYPE	scr = 0;
	RANGE	rng[2];
static	const	char smark[] = ">-^<";

	if (pwd->DvsP == 0 && !a->inex.intr && !b->inex.intr &&
	    algmode.nsa == 7) return (valonly(svr, sqs, pwd, q));
	if ((algmode.lcl & 16) && !algmode.qck && !algmode.lsg) {
	    Colonies*	clns = 0;			// swg DP mode
	    a->saverange(rng);
	    b->saverange(rng + 1);
	    if (algmode.mlt == 0) {
		scr = HomScore(sqs, pwd, 0);
		swgscore(svr, sqs, scr, q);
		goto match2end;
	    }
	    clns = swg1st(sqs, pwd);
	    if (!clns)  scr = JobError;
	    else if (algmode.mlt == 1 || (algmode.mlt > 2 && clns->size() == 0)) {
		Gsinfo AlnI;
		scr = (swg2nd(sqs, pwd, &AlnI, clns->at()) &&
		   (!algmode.thr || AlnI.fstat.val >= pwd->alnprm.thr))?
		    alnoutput(svr, sqs, pwd, &AlnI, q): JobError;
	    } else {
		for (int n = 0; n < clns->size(); ++n) {
		    Gsinfo AlnI;
		    scr = (swg2nd(sqs, pwd, &AlnI, clns->at(n)))?
			alnoutput(svr, sqs, pwd, &AlnI, q): JobError;
		    if (scr == JobError) break;
		}
	    }
	    delete clns;
match2end:
	    a->restrange(rng);
	    b->restrange(rng + 1);
	} else if (algmode.mlt == 0) {			// score & range
	    long	rr[2];
	    scr = HomScore(sqs, pwd, rr);
	    int	ml = a->left;
	    int	nl = b->left;
	    int	r = nl - ml;
	    if (rr[0] > r)  nl = a->left + rr[0];
	    else	    ml = b->left - rr[0];
	    int	mr = a->right;
	    int	nr = b->right;
	    r = nr - mr;
	    if (rr[1] > r)  mr = b->right - rr[1];
	    else	    nr = a->right + rr[1];
	    m_thread_Lock(q);
	    fprintf(out_fd, "%.1f\t%s %d %d %c\t%s %d %d %c\n", (double) scr,
		a->sqname(), a->SiteNo(ml), a->SiteNo(mr - 1), smark[a->inex.sens],
		b->sqname(), b->SiteNo(nl), b->SiteNo(nr - 1), smark[b->inex.sens]);
	    m_thread_Unlock(q);
	} else {					// (semi) global
	    Gsinfo GsI[2];	// abondan returned skl
	    int	dir = ((algmode.mns & 1) && 
		(GsI->skl = align2(sqs, pwd, &GsI->scr, GsI)))? 0: ERROR;
	    if (algmode.mns & 2) {
		a->comrev();
		antiseq(sqs + 1);
		if (b->jxt)
		    b->revjxt();
		if ((GsI[1].skl = align2(sqs, pwd, &GsI[1].scr, GsI + 1)))
		    dir = (dir == 0)? GsI[1].scr > GsI[0].scr: 1;
		if (dir == 0) {
		    a->comrev();
		    antiseq(sqs + 1);
		}
	    }
	    scr = (dir != ERROR && (!algmode.thr || GsI[dir].fstat.val >= pwd->alnprm.thr))?
		alnoutput(svr, sqs, pwd, GsI + dir, q): JobError;
	    if (dir == 1) {
		a->comrev();
		antiseq(sqs + 1);
	    }
	}
	return (scr);
}

void aln_setup(CalcServer<mSeq>* svr)
{
	setup_output(algmode.nsa);
#if SSHP
	initSsHpPrm();
#endif
	mSeq*&	a = svr->in_face[0];
	mSeq*&	b = svr->in_face[1];
	int	dvsp = prePwd((const Seq**) svr->in_face);
	if (algmode.qck || dvsp == 1 || dvsp == 2 || b->inex.intr)
	    makeWlprms(dvsp);
	dvsp = a->isprotein() + 2 * b->isprotein();
	if (dvsp == 1) b->inex.intr = algmode.lsg = algmode.aut;
	if (dvsp == 2) a->inex.intr = algmode.lsg = algmode.aut;
	if (dvsp == 1 || dvsp == 2) initcodon(genomedb);
	if (dvsp == 0) a->inex.intr = genomedb? 1: 0;
	if (a->inex.intr || b->inex.intr) makeStdSig53();
	mSeq*	sqs[2] = {a, b};
	bool	mlt[2] = {a->many > 1, b->many > 1};
	if (mlt[0]) {
	    sqs[0] = new mSeq(1, 1);
	    a->copyattr(sqs[0]);
	    sqs[0]->inex.vect = sqs[0]->inex.prof = sqs[0]->inex.algn = 0;
	}
	if (mlt[1]) {
	    sqs[1] = new mSeq(1, 1);
	    b->copyattr(sqs[1]);
	    sqs[1]->inex.vect = sqs[1]->inex.prof = sqs[1]->inex.algn = 0;
	}
	PwdM*	pwd = new PwdM(sqs);
	svr->prm = (void*) pwd;
	if (pwd->swp) swap(mlt[0], mlt[1]);
	if (mlt[0]) delete sqs[0];
	if (mlt[1]) delete sqs[1];
	spb_fact();
}

static PwdM* aln_reset(PwdM* pwd, mSeq** sqs)
{
	mSeq*&	a = sqs[0];
	mSeq*&	b = sqs[1];
	int	dvsp = a->isprotein() + 2 * b->isprotein();
	if (dvsp == 1) b->inex.intr = algmode.lsg;
	if (dvsp == 2) a->inex.intr = algmode.lsg;
	if (dvsp == 0) a->inex.intr = genomedb? 1: 0;
	if (!pwd) {
	    prePwd((const Seq**) sqs);
	    pwd = new PwdM(sqs);
	} else if (a->many > 1 || b->many > 1) {
	    int	swp = pwd->swp;
	    pwd = new PwdM(sqs);
	    if (dvsp == 1 || dvsp == 2) pwd->swp = swp;
	} else {
	    if (pwd->swp) swap(a, b);
	    a->mkthick();
	    b->mkthick();
	}
	return pwd;
}

void aln_cleanup(CalcServer<mSeq>* svr)
{
#if SSHP
	eraseSsHpPrm();
#endif
	delete (PwdM*) svr->prm;
}

static void genomicseq(Seq** seqs, PwdB* pwd, int ori)
{
	if (pwd->DvsP == 1) seqs[0]->nuc2tron();
	Intron53(seqs[0], pwd, ori == 3);
	if (ori == 1) return;
	seqs[0]->comrev(seqs + 1);
	seqs[1]->setanti(seqs);
	if (pwd->DvsP == 1) seqs[1]->nuc2tron();
	Intron53(seqs[1], pwd, ori == 3);
}

int aln_main(CalcServer<mSeq>* svr, mSeq** sqs, ThQueue<mSeq>* q = 0)
{
	mSeq*&	a = sqs[0];
	polyA.rmpolyA(a, 1);
	mSeq*&	b = sqs[1];
	polyA.rmpolyA(b, 1);
	if (algmode.lcl & 16) {	// local
	    a->exg_seq(1, 1);
	    b->exg_seq(1, 1);
	} else {		// global
	    a->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
	    b->exg_seq(algmode.lcl & 4, algmode.lcl & 8);
	}
	PwdM*	pwd = aln_reset((PwdM*) svr->prm, sqs);
	switch (pwd->DvsP) {
	    case 1: a->deamb(2); break;
	    case 2: b->deamb(2); break;
	    default: break;
	}
	if (algmode.lsg) genomicseq((Seq**) (sqs + 1), pwd, 
	    pwd->DvsP == 0? algmode.mns: 1);
	VTYPE	rv = match_2(svr, sqs, pwd, q);
	if (pwd->swp) {
	    swap(a, b);
	    pwd->swapDvsP();
	}
	if (pwd != (PwdM*) svr->prm) delete pwd;
	return (rv == JobError? ERROR: OK);
}

int aln_sub(CalcServer<mSeq>* svr, mSeq** sqs, ThQueue<mSeq>* q = 0)
{
	VTYPE	scr = HomScore(sqs, (PwdM*) svr->prm);
	m_thread_Lock(q);
	svr->avsd[0] += scr;
	svr->avsd[1] += scr * scr;
	m_thread_Unlock(q);
	return (OK);
}

/********************************************************************
*
*	Construct multiple sequence alignment by a progressive method
*
********************************************************************/

inline bool isntleaf(Tnode* node) {return (node->left || node->right);}

static	const	char* prg_frmt = "%10s <=== %10s + %10s\n";

template<>
char* MakeMsa<mSeq>::seqnam(Tnode* root)
{
	if (!root->tname) {
	    char*	sn1 = seqnam(root->left);
	    char*	sn2 = seqnam(root->right);
	    sprintf(wrkvar, "%d", root->tid);
	    prompt(prg_frmt, wrkstr, sn1, sn2);
	    mSeq*&	a = svr->in_face[0];
	    mSeq*&	b = svr->in_face[1];
	    if (svr->logfd) fprintf(svr->logfd, prg_frmt, wrkstr, sn1, sn2);
	    if (!a->getseq(sn1)) fatal(not_found, sn1);
	    if (!b->getseq(sn2)) fatal(not_found, sn2);
	    root->tname = strrealloc(0, wrkstr);
	    if ((out_fd = fopen(wrkstr, "w")) == 0)
		fatal(erropen, wrkstr);
	    if (!svr->input_ns++) prePwd(a);
	    aln_main(svr, svr->in_face);
	    fclose(out_fd); out_fd = stdout;
	    if (OutPrm.RemoveTmp) {
		if (isntleaf(root->left)) remove(sn1);
		if (isntleaf(root->right)) remove(sn2);
	    }
	    if (isntleaf(root->left)) delete[] sn1;
	    if (isntleaf(root->right)) delete[] sn2;
	}
	return (root->tname);
}

template<>
MakeMsa<mSeq>::MakeMsa(AlnServer<mSeq>* sv, Tnode* root) : svr(sv)
{
	sprintf(wrkstr, "%smsa%d_", svr->tmpdir, getpid());
	int	wrklen = strlen(wrkstr);
	wrkvar = wrkstr + wrklen;
	seqnam(root);
	if (rename(root->tname, OutPrm.out_file && *OutPrm.out_file? OutPrm.out_file: Def_out_file))
	    perror("Error renaming file ");
	if (isntleaf(root)) delete[] root->tname;
}

template<>
MakeMsa<mSeq>::MakeMsa(AlnServer<mSeq>* sv): svr(sv)
{
	char	tname[LINE_MAX];
	int	tid = 0;

	strcpy(wrkstr, svr->tmpdir);
	strcat(wrkstr, tmpchr);
	int	wrklen = strlen(wrkstr);
	wrkvar = wrkstr + wrklen;
	int	i = 0;
	char*	ss = 0;
	mSeq*&	ins_a = svr->in_face[0];
	mSeq*&	ins_b = svr->in_face[1];
	if (!ins_a->getseq(svr->memb->name[i++]))
	    fatal("Nothing to do !\n");
	while (i < svr->memb->membersize &&
	    ins_b->getseq(svr->memb->name[i++])) {
	    sprintf(wrkvar, "%d", tid++);
	    prompt(prg_frmt, wrkstr, ins_a->sqname(), ins_b->sqname());
	    if (svr->logfd) fprintf(svr->logfd, prg_frmt, wrkstr, 
		ins_a->sqname(), ins_b->sqname());
	    if ((out_fd = fopen(wrkstr, "w")) == 0)
		fatal(erropen, wrkstr);
	    if (!svr->input_ns++) prePwd(ins_a);
	    svr->main_job(svr, svr->in_face, 0);
	    fclose(out_fd); out_fd = stdout;
	    if (ss && OutPrm.RemoveTmp) remove(ss);
	    ss = strcpy(tname, wrkstr);
	    if (!ins_a->getseq(ss)) {
		prompt("Incomplete !\n");
		ss = 0;
		break;
	    }
	}
	if (OutPrm.RemoveTmp && ss && rename(ss, OutPrm.out_file && *OutPrm.out_file?
		OutPrm.out_file: Def_out_file))
	    perror("Error renaming file ");
}

static void setdefparam()
{
	optimize(GLOBAL, MAXIMUM);
	alprm.sh = -50;			// band shoulder = 1/2 of shorter
	algmode.thr = 0;		// disable threshold
	algmode.mns = 1;
	algmode.aut = 1;
	algmode.crs = 1;
	OutPrm.trimend = true;		// trim dangling ends with -L
//	setdefNprm(-2., 2., 4.);	// n, u, v
	setdefPprm(250, 2., 9.);	// pam, u, v
	setpam(WlpPam, WlnPamNo);	// pam for HSP search
	setalgmode(1, 0);
}

int main(int argc, const char** argv)
{
	setdefparam();
	AlnServer<mSeq>	svr(argc, argv, IM_NONE, IM_ALTR, (void*) 0, 
	&aln_main, &aln_setup, &aln_cleanup, &aln_sub);
	int	rv = (svr.calc_mode == IM_ADON || svr.calc_mode == IM_TREE)?
		svr.makemsa(): svr.autocomp();
	if (rv == 1) svr.menucomp();
	EraDbsDt();
	EraStdSig53();
	eraWlprms();
	eraStrPhrases();
	resetSimmtxes(true);
	return (rv);
}
