/*****************************************************************************
*
*	Use prosite.dat and prosite.doc files through indeces
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

#include "seq.h"
#include "prs.h"
#include "pattern.h"

static 	int	cmpdocidx(DOC_IDX* a, DOC_IDX* b);
static 	int	cmppatidx(PAT_IDX* a, PAT_IDX* b);

static	const	char*	FN_DAT = "prosite.dat";
static	const	char*	FN_DOC = "prosite.doc";
static	const	char*	FN_PAT = "prosite.pat";
static	const	char*	FN_DOC_IDX = "prosite.doc_idx";
static	const	char*	FN_PAT_IDX = "prosite.pat_idx";
static	const	char*	DOC_END = "{END}";
static	const	char*	DAT_END = "//";

static int cmpdocidx(DOC_IDX* a, DOC_IDX* b)
{
	return (a->doc - b->doc);
}

static int cmppatidx(PAT_IDX* a, PAT_IDX* b)
{
	return (a->pat - b->pat);
}

void ProSite::read_doc_idx()
{
	FILE* fd = fopenpbe(path, FN_DOC_IDX, 0, "r", 1);

	if (fseek(fd, 0L, 2)) fatal("Seek error!\n");
	long  tail = ftell(fd);
	doc_nel = tail / sizeof(DOC_IDX);
	doc_idx = new DOC_IDX[doc_nel];
	rewind(fd);
	if (fread(doc_idx, sizeof(DOC_IDX), doc_nel, fd) != doc_nel)
	    fatal("%s: read error!\n", FN_DOC_IDX);
	fclose(fd);
	max_doc = doc_idx[doc_nel - 1].doc;
}

void ProSite::read_pat_idx()
{
	FILE* fd = fopenpbe(path, FN_PAT_IDX, 0, "r", 1);

	if (fseek(fd, 0L, 2)) fatal("Seek error!\n");
	long  tail = ftell(fd);
	pat_nel = tail / sizeof(PAT_IDX);
	pat_idx = new PAT_IDX[pat_nel];
	rewind(fd);
	if (fread(pat_idx, sizeof(PAT_IDX), pat_nel, fd) != pat_nel)
	    fatal("%s: read error!", FN_PAT_IDX);
	fclose(fd);
	max_pat = pat_idx[pat_nel - 1].pat;
}

DOC_IDX* ProSite::prs_doc_idx(long ptr)
{
	DOC_IDX	    key;

	key.doc = ptr;
	return (DOC_IDX*)
	    bsearch(&key, doc_idx, doc_nel, sizeof(DOC_IDX), (CMPF) cmpdocidx);
}

PAT_IDX* ProSite::prs_pat_idx(long ptr)
{
	PAT_IDX	    key;

	key.pat = ptr;
	return (PAT_IDX*)
	    bsearch(&key, pat_idx, pat_nel, sizeof(PAT_IDX), (CMPF) cmppatidx);
}

void ProSite::put_prs_doc(FILE* fd, long ptr)
{
	char	str[MAXL];

	fseek(fdoc, ptr, 0);
	while (fgets(str, MAXL, fdoc)) {
	    fputs(str, fd);
	    if (!wordcmp(str, DOC_END)) break;
	}
}

void ProSite::put_prs_dat(FILE* fd, long ptr)
{
	char	str[MAXL];

	fseek(fdat, ptr, 0);
	while (fgets(str, MAXL, fdat)) {
	    fputs(str, fd);
	    if (!wordcmp(str, DAT_END)) break;
	}
}

char* ProSite::get_prs_pat(char* str, long ptr)
{
	fseek(fpat, ptr, 0);
	chop(fgets(str, MAXL, fpat));
	return (str);
}

ProSite::ProSite(int wh, const char* dir)
{
	path = dir? dir: DBS_DIR;
	read_doc_idx();
	read_pat_idx();
	if ((wh & DAT) && !fdat) fdat = fopenpbe(path, FN_DAT, 0, "r", 1);
	if ((wh & DOC) && !fdoc) fdoc = fopenpbe(path, FN_DOC, 0, "r", 1);
	if ((wh & PAT) && !fpat) fpat = fopenpbe(path, FN_PAT, 0, "r", 1);
}

ProSite::~ProSite()
{
	delete[] doc_idx;
	delete[] pat_idx;
	if (fdat) fclose(fdat);
	if (fdoc) fclose(fdoc);
	if (fpat) fclose(fpat);
}

int ProSite::propos(int loc[], int maxloc, Seq* sd, long psn)
{
	PAT_IDX*    pat_rec;
	char	str[MAXL];
	int	nn = 0;

	if (psn > max_pat) return (-1);
	pat_rec = prs_pat_idx(psn);
	if (pat_rec) {
	    get_prs_pat(str, pat_rec->pat_ptr);
	    nn = findseq(loc, maxloc, sd, str);
	}
	return (nn);
}

void ProSite::prosites(FILE* fd, Seq* sd, char* pat)
{
	char*	pb;
	int	loc[MAXLOC];

	while (*(pb = pat)) {
	    while (*pat && *pat != ' ') pat++;
	    if (*pat) *pat++ = '\0';
	    int	num = propos(loc, MAXLOC, sd, atol(pb));
	    if (num > 0) putloc(fd, sd, loc, num);
	}
}

void ProSite::allprs(FILE* fd, Seq* sd)
{
	long	i = 0L;
	int	loc[MAXLOC];
	int	num;

	while ((num = propos(loc, MAXLOC, sd, ++i)) >= 0) {
	    if (num) {
		fprintf(fd, "PS%05ld: ", i);
		putloc(fd, sd, loc, num);
	    }
	}
}

