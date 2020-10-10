/*****************************************************************************
*
*       Structure-associated properties of protein sequence
*
*       Osamu Gotoh, ph.D.      (-2001)
*       Saitama Cancer Center Research Institute
*       818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*       Osamu Gotoh, Ph.D.      (2001-)
*       National Institute of Advanced Industrial Science and Technology
*       Computational Biology Research Center (CBRC)
*       2-41-6 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*       Osamu Gotoh, Ph.D.      (2003-2012)
*       Department of Intelligence Science and Technology
*       Graduate School of Informatics, Kyoto University
*       Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*       Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#include "seq.h"
#include "simmtx.h"
#include "ssp.h"
#include <math.h>

static  const   char    HphPrmFile[] = "hydro";
static  const   char    SndPrmFile[] = "gor3";
static	const	char*	sshplbl[] = {
	"helix", "sheet", "coil", "hydpho", "hm100", "hm180"
};

SsHpPrm* sshpprm = 0;		   // global 
static  float   angle[] = {100., 180};  // hydrophobicity moment

SsHpPrm::SsHpPrm()
{
	FILE*	fd = ftable.fopen(sshp_data, "r");
	if (!fd) fatal("Can't open %s!\n", sshp_data);
	if ((fread(params, sizeof(params), 1, fd) != 1) ||
	    (fread(phptbl, sizeof(phptbl), 1, fd) != 1) ||
	    (fread(psstbl, sizeof(psstbl), 1, fd) != 1) ||
	    (fread(sshpav, sizeof(sshpav), 1, fd) != 1) ||
	    (fread(sshpsd, sizeof(sshpsd), 1, fd) != 1) ||
	    (fread(sincrv, sizeof(sincrv), 1, fd) != 1) ||
	    (fread(coscrv, sizeof(coscrv), 1, fd) != 1))
		fatal("Error in read sshp.dat !\n");
	fclose(fd);
	hphtbl = phptbl[params[0]];		// hydrophobicity
	hmttbl = phptbl[params[1]];		// hydrophobicity moment
	sndstates = (alprm3.scnd > 0)? NOSS: 0;
	hphstates = (alprm3.hydr > 0)? 1: 0;
	hmtstates = alprm3.no_angle;
	sshpelems = sndstates + hphstates + hmtstates;
	hpwidth = 2 * alprm3.hpwing + 1;
	if (sndstates) sshpsd[NOSS] /= sqrt(double(hpwidth));
}

SsHpPrm::SsHpPrm(HPSCALE hps, HPSCALE hms)
{
	read_prm((PROFL*) psstbl, SndPrmFile, NOSS * SSWIDTH);  // 2nd structure
	read_prm(phptbl, HphPrmFile, NOHP);
	for (int th = 0; th < 2; ++th) {	// hydrophobicity moment
	    double   period = PI * angle[th] / 360.;
	    double   ang = 0;
	    for (int i = 0; i <= HMWING; ++i, ang += period) {
		sincrv[th][i] = sin(ang);
		coscrv[th][i] = cos(ang);
	    }
	}
	params[0] = hps;
	params[1] = hms;
	hphtbl = phptbl[(int) hps];
	hmttbl = phptbl[(int) hms];
	hpwidth = sndstates = hphstates = hmtstates = sshpelems = 0;
}

void SsHpPrm::read_prm(PROFL* tbl, const char* fname, int rec)
{
	char    str[MAXL];
	FILE*   fd = ftable.fopen(fname, "r");
	if (!fd) {
	    strcpy(str, fname);
	    strcat(str, ".txt");
	    fd = ftable.fopen(str, "r");
	}
	if (!fd) fatal("%s cannot open !\n", fname);

	for (int i = 0; i < rec && fgets(str, MAXL, fd); ) {
	    if (*str == '\n') continue;
	    if (*str == _LCOMM) {
		flush_line(fd);
		continue;
	    }
	    float*      pf = tbl[i++];
	    char*       ps = str;
	    for (int j = 0; j < 20; ++j) {
		ps = cdr(ps);
		if (!ps || !*ps)
		    fatal("Insufficient data in %s\n", fname);
		*pf++ = atof(ps);
	    }
	}
	fclose(fd);
}

void SsHpPrm::avsd()
{
	double*	pcmp = get_mdmcmp();
	vclear(sshpav, NOSS + HPDIM);
	vclear(sshpsd, NOSS + HPDIM);
	int	s = 0;
	for ( ; s < NOSS; ++s) {
	    for (int j = 0; j < SSWIDTH; ++j) {
		double* pc = pcmp;
		for (int a = 0; a < 20; ++a, ++pc) {
		    float&	v = psstbl[s][j][a];
		    sshpav[s] += *pc * v;
		    sshpsd[s] += *pc * v * v;
		}
	    }
	    sshpsd[s] = sqrt(sshpsd[s] - sshpav[s] * sshpav[s]);
	}
	double*	pc = pcmp;
	for (int a = 0; a < 20; ++a, ++pc) {
	    float&      v = hphtbl[a];
	    sshpav[s] += *pc * v;
	    sshpsd[s] += *pc * v * v;
	}
	sshpsd[s] = sqrt(sshpsd[s] - sshpav[s] * sshpav[s]);
	float	hm = 0;
	float	hm2 = 0;
	pc = pcmp;
	int	w = 2 * HMWING + 1;
	for (int a = 0; a < 20; ++a, ++pc) {
	    float&	v = hmttbl[a];
	    hm += *pc * v;
	    hm2 += *pc * v * v;
	}
	for (int k = 0; k < 2; ++k) {
	    ++s;
	    float	ks = 0;
	    float	period = PI * angle[k] / 360.;
	    float	th = period;
	    for (int d = 1; d < w; ++d, th += period) 
		ks += (w - d) * cos(th);
	    sshpsd[s] =  w * hm2 + 2 * ks * hm * hm;
	    sshpav[s] =  0;
	    for (int a = 0; a < 20; ++a) {
		float	t = 0;
		for (int b = 0; b < 20; ++b)
		    t += pcmp[b] * sqrt(w * hm2 + ks * hmttbl[a] * hmttbl[b]);
		sshpav[s] += pcmp[a] * t;
	    }
	    sshpsd[s] = sqrt(sshpsd[s] - sshpav[s] * sshpav[s]);
	}
}

void SsHpPrm::simul_avsd(int nrep)
{
	if (nrep < 2) return;
	vclear(sshpav, NOSS + HPDIM);
	vclear(sshpsd, NOSS + HPDIM);
	RandNumGen	rn(get_mdmcmp(), 20);
	CHAR	que[SSWIDTH];
	for (int q = 0; q < SSWIDTH; ++q)
	    que[q] = rn.get();
	float	ss[NOSS + HPDIM];
	for (int n = 0, qp = 0; n < nrep; ++n) {
	    vclear(ss, NOSS + HPDIM);
	    for (int j = 0, q = qp; j < SSWIDTH; ++j) {
		for (int s = 0; s < NOSS; ++s)
		    ss[s] += psstbl[s][j][que[q]];
		if (++q == SSWIDTH) q = 0;
	    }
	    ss[NOSS] += hphtbl[que[qp]];
	    float	hpp[2][2] = {{0., 0.,}, {0., 0.,}};
	    for (int j = 0, q = qp; j <= 2 * HMWING; ++j) {
		int	jj = j - HMWING;
		for (int k = 0; k < 2; ++k) {
		    float sv = jj? sincrv[k][abs(jj)] * hmttbl[que[q]]: 0;
		    hpp[k][0] += ((jj > 0)? sv: -sv);
		    hpp[k][1] += coscrv[k][abs(jj)] * hmttbl[que[q]];
		}
		if (++q == SSWIDTH) q = 0;
	    }
	    for (int k = 0; k < 2; ++k)
		ss[NOSS + 1 + k] = sqrt(hpp[k][0] * hpp[k][0]
				      + hpp[k][1] * hpp[k][1]);
	    for (int s = 0; s < NOSS + HPDIM; ++s) {
		sshpav[s] += ss[s];
		sshpsd[s] += ss[s] * ss[s];
	    }
	    que[qp] = rn.get();
	    if (++qp == SSWIDTH) qp = 0;
	}
	for (int s = 0; s < NOSS + HPDIM; ++s) {
	    float	sm = sshpav[s];
	    sshpav[s] /= nrep;
	    sshpsd[s] = sqrt((sshpsd[s] - sm * sshpav[s]) / (nrep - 1));
	}
}

void SsHpPrm::write_data(FILE* fd)
{
	if ((fwrite(params, sizeof(params), 1, fd) != 1) ||
	    (fwrite(phptbl, sizeof(phptbl), 1, fd) != 1) ||
	    (fwrite(psstbl, sizeof(psstbl), 1, fd) != 1) ||
	    (fwrite(sshpav, sizeof(sshpav), 1, fd) != 1) ||
	    (fwrite(sshpsd, sizeof(sshpsd), 1, fd) != 1) ||
	    (fwrite(sincrv, sizeof(sincrv), 1, fd) != 1) ||
	    (fwrite(coscrv, sizeof(coscrv), 1, fd) != 1))
		fatal("Error in writing sshp.dat !\n");
}

void SsHpPrm::report()
{
	for (int s = 0; s < NOSS + HPDIM; ++s)
	    printf("%s\t%7.2f\t%7.2f\n", sshplbl[s], sshpav[s], sshpsd[s]);
}

SsHpPrm* initSsHpPrm()
{
	if (alprm3.scnd == 0.&& alprm3.hydr == 0. && alprm3.hpmt == 0) return (0);
	if (alprm3.hpmt > 0. && !alprm3.no_angle)
	    alprm3.no_angle = 1;
	if (alprm3.no_angle && alprm3.hpmt == 0)
	    alprm3.hpmt = alprm3.hydr;
	return (sshpprm = new SsHpPrm());
}

void eraseSsHpPrm() {delete sshpprm;}

#if MAIN

int main(int argc, const char** argv)
{
	int	nrepeat = 0;
	while (--argc && (*++argv)[0] == '-') {
	    switch (argv[0][1]) {
		case 'n': nrepeat = atoi(argv[0] + 2); break;
		default: break;
	    }
	}
	FILE*	fd = argc? fopen(*argv, "w"): ftable.fopen(sshp_data, "w");
	if (!fd) fatal("Can't write to %s\n", argc? argv[1]: sshp_data);
	alprm3.scnd = alprm3.hydr = alprm3.hpmt = 1;
	alprm3.no_angle = 2;
	sshpprm = new SsHpPrm(KYTDO, PRIFT);
	if (nrepeat) sshpprm->simul_avsd(nrepeat);
	else	sshpprm->avsd();
	sshpprm->report();
	sshpprm->write_data(fd);
	fclose(fd);
	delete sshpprm;
	return (0);
}

#endif	// MAIN
