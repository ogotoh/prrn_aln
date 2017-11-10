/*****************************************************************************
*
*	Automatic performance of program "perform" for a given set
*	of protein or nucleotide sequences.
*
*	Subcommand: 'e' all; 		'f' single scan; 	'g' group
*	Subcommand: 'i' internal; 	'j' juxtapose		'k' addon
*	Subcommand: '1' onesequence;
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

#include	"autocomp.h"
#include	"seq.h"

InputSeqTest	zero_InputSeqTest = {0, 0, 0, 0, 0, 0, 0};

int testInputSeq(CalcServer<Seq>* svr, Seq* sqs[], ThQueue<Seq>* q)
{
	InputSeqTest*	tv = (InputSeqTest*) svr->prm;
	Seq*&	sd = sqs[0];
	if (tv->molc && sd->inex.molc != tv->molc) {
	    prompt("%s is mixed type !\n", (*sd->spath)[0]);
	    ++tv->bad;
	} else if (!sd->many) {
	    prompt("%s is empty !\n", (*sd->spath)[0]);
	    ++tv->bad;
	} else {
	    tv->molc = sd->inex.molc;
	    tv->many += sd->many;
	    if (sd->many > tv->maxmany) tv->maxmany = sd->many;
	    tv->space += sd->many * sd->len;
	    if (tv->sname) {
		const char*	name = (*sd->sname)[0];
		if (sd->many > 1) {
		    const char* path = (*sd->spath)[0];
		    if (path) {
			const char*	slash = strrchr(path, '/');
			name = slash? slash + 1: path;
		    }
		}
		tv->sname->push(name);
	    }
	    ++tv->num;
	}
	return (OK);
}

