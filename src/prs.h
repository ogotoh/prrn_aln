/*****************************************************************************
*
*	Header file for use of prosite data
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

#ifndef _PRS_H_
#define _PRS_H_

enum {DAT = 1, DOC = 2, PAT = 4, ALL = 7};

struct DOC_IDX {long doc; long doc_ptr;};
struct PAT_IDX {long pat; long pat_ptr, dat_ptr, doc_ptr;};

class ProSite {
	const	char*	path;
	DOC_IDX* doc_idx;
	PAT_IDX* pat_idx;
	INT	doc_nel;
	INT	pat_nel;
	long	max_pat;
	long	max_doc;
	FILE*	fdat;
	FILE*	fdoc;
	FILE*	fpat;
	void	read_doc_idx();
	void	read_pat_idx();
	int	propos(int loc[], int maxloc, Seq* sd, long psn);
public:
	ProSite(int wh, const char* dir = 0);
	~ProSite();
	DOC_IDX*	prs_doc_idx(long ptr);
	PAT_IDX*	prs_pat_idx(long ptr);
	void	put_prs_doc(FILE* fd, long ptr);
	void	put_prs_dat(FILE* fd, long ptr);
	char*	get_prs_pat(char* str, long ptr);
	void	prosites(FILE* fd, Seq* sd, char* pat);
	void	allprs(FILE* fd, Seq* sd);
	void 	set_prs_dir(const char* dir);
};

#endif
