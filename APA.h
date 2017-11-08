#ifndef __APA__H__
#define __APA__H__

#include"helperfunc.h"
#include<random>

class SUS
{
public:
	SUS(int par_nCend, double *par_u01, double *par_sig01, double par_t_C, int par_randseed, int par_nSimiter);
	~SUS();
	std::mt19937 e;
	vector<double> APA_probEstList;
	vector<double> APA_probCILowList;
	vector<double> APA_probCIUpList;
	vector<int> APA_simtotalList;
	int nCend;
	double u01[2];
	double sig01[2];
	double t_C;
	int randseed;
	int nSimiter;

	bool MC(vector<int> &src,vector<double> &dst,int simtype,int nsim=1E9);	
	bool sus_delta_sim(vector<vector<double> >&XSeed, vector<vector<double> >&Yseed, vector<vector<double> >&XsaSeed, vector<vector<double> >&ylimSeed,int simtype);
	bool genX(vector<vector<double> >&XSeed, vector<vector<double> >&YSeed, vector<double>&Y1Seed, vector<double>&Y2Seed, vector<vector<double> >&XsaSeed, vector<vector<double> >&ylimSeed, int &nSimTotal, int  nSimiter, vector<double> perfSpec, vector<vector<vector<double> > >&Res_X, vector<vector<vector<double> > > &Res_Y, vector<vector<vector<double> > >&Res_Y1, vector<vector<vector<double> > >&Res_Y2, vector<vector<vector<double> > >&Res_Xsa, vector<vector<vector<double> > >&Res_ylim,int simtype);
private:
	bool sim_SA_fake(vector<vector<double> >&Xsa, vector<double>&Res_ylimit,int simtype);
	bool simout(vector<vector<double> >&X, vector<vector<double> >&Y, int &nSimTotal, bool epo,int simtype);
	bool simout2(vector<vector<double> >&X, vector<vector<double> >&Y, int &nSimTotal, bool epo,int simtype);
	bool simX(vector<vector<double> > &src, vector<vector<double> >&dst, int index,int simtype);
	bool output(vector<double>probList, vector<double>sigList,int &nSimTotal);
	bool getspec(vector<double> Y1, vector<double> Y2, double probtarget, vector<double> &perdeltaCur);
	bool filterX(vector<vector<vector<double> > > &X, vector<vector<vector<double> > >&Y, vector<vector<vector<double> > >&Y1, vector<vector<vector<double> > >&Y2, vector<vector<vector<double> > >&Xsa, vector<vector<vector<double> > >&ylim, vector<double>&probList, vector<double>&sigList, vector<vector<double> >&perfDeltaList, int &nSimTotal, int nSimiter, vector<vector<double> > &XSeed, vector<vector<double> > &YSeed, vector<double>&Y1Seed, vector<double>&Y2Seed, vector<vector<double> >&XsaSeed, vector<vector<double> >&ylimSeed, int &eop,int simtype);
	bool expandSeed(vector<vector<double> >&XSeed, vector<vector<double> >&YSeed, vector<vector<double> >&XsaSeed, vector<vector<double> >&ylimSeed, int & nSimTotal, int nSimiter, vector<vector<double> >&Res_XSeed, vector<vector<double> >&Res_YSeed, vector<vector<double> >&Res_Xsa, vector<vector<double> >&Res_ylim, vector<vector<vector<double> > >&Xcell, vector<vector<vector<double> > >&Ycell,int simtype);
	bool metropolis(vector<double>&X, int nVar, vector<vector<double> > &XNext);
	double normpdf(double src, double mu = 0, double sigma = 1);
	bool transpose(vector<double>&src, vector<vector<double> > &dst);
	bool transpose(vector<vector<double> >&src, vector<double> &dst);
	bool gen_corr(vector<double> &bias,int simtype);
	vector<vector<double> > transpose(vector<vector<double> >&src);
};



#endif // !__APA__h__

