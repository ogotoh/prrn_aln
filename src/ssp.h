/*****************************************************************************
*
*       Use predicted secondary structure and hydrophobicity profiles
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

#ifndef  _SSP_H_
#define  _SSP_H_

enum HPSCALE {PRIFT, KYTDO, HOPOW, CHARGE};

static  const   double  PI = 3.14159265358979323846;
static  const   int     NOSS = 3;
static  const   int     NOHP = 4;
static	const	int	HPDIM = 3;
static  const   int     SSWING = 8;
static  const   int     SSWIDTH = 17;
static  const   int     HMWING = 4;
static	const	int	ssidx = 0;
static	const	int	hpidx = NOSS;
static	const	int	hmidx = hpidx + 1;
static	const	char	sshp_data[] = "sshp.data";

typedef float   PROFL[20];
typedef float   TRIGONO[HMWING + 1];
typedef	PROFL	(*SCNDRY)[SSWIDTH];

class SsHpPrm {
	void    read_prm(PROFL* tbl, const char* fname, int rec);
public:
	int	params[2];
	PROFL	phptbl[4];
	PROFL	psstbl[NOSS][SSWIDTH];
	float	sshpav[NOSS + HPDIM];
	float	sshpsd[NOSS + HPDIM];
	TRIGONO	sincrv[2];
	TRIGONO	coscrv[2];
	float*	hphtbl;
	float*	hmttbl;
	int	hpwidth;
	int	sndstates;
	int	hphstates;
	int	hmtstates;
	int	sshpelems;
	SsHpPrm();
	SsHpPrm(HPSCALE hps, HPSCALE hms);
	~SsHpPrm() {}
	void	avsd();
	void	simul_avsd(int nrep);
	void	write_data(FILE* fd);
	void	report();
};

extern  SsHpPrm* sshpprm;		       // global 
extern  SsHpPrm* initSsHpPrm();
extern  void eraseSsHpPrm();

#endif  // _SSP_H

