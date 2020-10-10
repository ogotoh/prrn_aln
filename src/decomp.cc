/******************************************************************************		*
*	Decompose bundled sequence data files
*		GenBank, EMBL, SeqDB, NBRF(new  and old)
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
#include <time.h>

static	const	int	DATE_COL = 62;
static	const	int	GenSpcLen = 8;

static	void	usage();
static	int	newer(const char* entry);
static	void	setpath(const char* str, const char* gs);
static	char*	fname(char* str, SeqDb* dbf);
static	void	decompdb(FILE* fin, const char* genspc = 0);

static	char	path[LINE_MAX];
static	char*	fn;
static	const	char*	gdate = 0;
static	struct	tm	given_day;
static	int	AfterEnd = 0;
static	int	field_sps = 0;
static	int	field_bar = INT_MAX;
static	bool	quiet = false;

static void usage()
{
	fputs("Usage:\n", stderr);
	fputs("\tdecomp [-pOutPut_Path] [-n\"Day-MON-Year\"] [-fieldN] [-hGenspc] [Input_File]\n", stderr);
}

static int newer(const char* entry)
{
	struct tm	entry_day;

	strptime(entry, "%d-%b-%Y", &entry_day);
	return (entry_day.tm_year > given_day.tm_year ||
		(entry_day.tm_year == given_day.tm_year && 
		entry_day.tm_mon > given_day.tm_mon) ||
		(entry_day.tm_mon == given_day.tm_mon &&
		entry_day.tm_mday >= given_day.tm_mday));
}
	
static void setpath(const char* str, const char* gs)
{
	topath(path, str);
	fn = strrchr(path, PATHDELM);
	if (fn == 0)	fn = path;
	else		++fn;
	if (gs) {
	    strncpy(fn, gs, GenSpcLen);
	    int	r = GenSpcLen - strlen(fn);
	    if (r) memset(fn + GenSpcLen - r, '_', r);
	    fn += GenSpcLen;
	}
}

static char* fname(char* str, SeqDb* dbf)
{
	char	locus[LINE_MAX + 1];
static	const	char*	delm = " \t\n";

	if (dbf->FormID == ProDB) {
	    strcpy(locus, "P");
	    strncat(locus, str, LINE_MAX - 3);
	} else
	    strncpy(locus, str, LINE_MAX);
	char*	ps = strtok(locus, delm);
	char*	ss = 0;
	for (int i = 0; i < field_sps; ++i)
	    ps = strtok(0, delm);
	if (field_bar >= 0) {
	    ps = strtok(ps, "|");
	    for (int i = 0; ps && i < field_bar; ++i) {
		ss = ps;
		ps = strtok(0, "|");
	    }
	}
	if (!ps) ps = ss;
	if (!ps) return (0);
	for (ss = fn; (*ss = *ps); ++ps)
	    if (isalnum(*ss) || *ss == '.' || *ss == '_') ++ss;
	return (path);
}

//	>sid || SDI sid
//	path/[genspc]fn

static void decompdb(FILE* fin, const char* genspc)
{
	char	str[MAXL];

	*str = '\0';
	while (!feof(fin)) {
loophead:
	    SeqDb*	dbf = whichdb(str, fin);
	    if (dbf == 0)
		fatal("%s\nNot a DataBase File!", str);
	    bool	fasta = dbf->FormID == FASTA;
	    char*	sid = fasta? str + 1: cdr(str);
	    char*	fnm = (!gdate || newer(str + DATE_COL))?
			fname(sid, dbf): 0;
	    if (genspc)	{
		int	sidlen = strlen(sid);
		memmove(sid + GenSpcLen, sid, sidlen);
		sid[sidlen + GenSpcLen] = '\0';
		memcpy(sid, fn - GenSpcLen, GenSpcLen);
	    }
	    FILE*	fout = fnm? fopen(fnm, "w"): 0;
	    if (fnm && !quiet) printf("%s: ", fnm);
	    if (!quiet) fputs(sid, stdout);
	    for (;;) {
		if (fout && fputs(str, fout) == EOF)
		    fatal("Write Error %s!", path);
		if (!fgets(str, MAXL, fin)) {
		    fclose(fout);
		    return;
		}
		if (fasta && *str == '>') {
		    fclose(fout);
		    goto loophead;
		}
		if (!AfterEnd && dbf->EndLabel && !wordcmp(str, dbf->EndLabel))
		    break;
	    }
	    if (fout) {
		fputs(str, fout);
		fclose(fout);
	    }
	    while (fgets(str, MAXL, fin) && isBlankLine(str)) ;
	}
}

int main(int argc, const char** argv)
{
	time_t	gmt;
	struct	tm *gd;
	bool	defheader = false;
	const	char*	dir = 0;
	const	char*	genspc = 0;

	while (--argc && (*++argv)[0] == '-') {
	    const char*	opt = argv[0] + 1;
	    if (isdigit(*opt)) {
		field_sps = atoi(opt);
		if ((opt = strchr(opt, '.')))
		    field_bar = atoi(opt + 1);
	    } else {
	      switch (*opt) {
		case 'p':
		case 'P':
		    if (opt[1])	dir = opt + 1;
		    else if (argc > 1 && argv[1][0] != '-')
			{argc--; dir = *++argv;}
		    else	dir = "";;
		    break;
		case 'h':
		    if (opt[1])	genspc = opt + 1;
		    else if (argc > 1 && argv[1][0] != '-')
			{argc--; genspc = *++argv;}
		    else	defheader = true;;
		    break;
		case 'n':
		    if (opt[1])	gdate = opt + 1;
		    else if (argc > 1 && argv[1][0] != '-')
			{argc--; gdate = *++argv;}
		    break;
		case 'e': AfterEnd = 1; break;
		case 'q': quiet = true; break;
		default:
		    usage();
		    exit (0);
	      }
	    }
	}
	if (gdate) {
		gmt = time(0);
		gd = gmtime(&gmt);
		given_day = *gd;
	 	strptime(gdate, "%d-%b-%Y", &given_day);
	}
	setpath(dir, genspc);
	if (argc > 0) {
	    for ( ; argc--; ++argv) {
		FILE*	fin = fopen(*argv, "r");
		if (fin == 0) fatal("%s not found", argv);
		if (defheader) {
		    genspc = strrchr(*argv, PATHDELM);
		    setpath(dir, genspc? genspc + 1: *argv);
		}
		decompdb(fin, genspc);
		fclose(fin);
	    }
	} else
		decompdb(stdin);

	exit (0);
}
