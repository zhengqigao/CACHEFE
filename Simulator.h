#ifndef __SIMULATOR__H__
#define __SIMULATOR__H__

#include"helperfunc.h"
#include"Eigen"
#include <random>

//sram circuit simulation class
class simulator
{
public:
	simulator();
	~simulator();
	int M_ND_[NumTrans] = { 2,1,2 }, M_NG_[NumTrans] = { 1,1,1 }, M_NS_[NumTrans] = { -1,2,1 };
	int M_TYPE_[2 * NumTrans] = { NMOS_,NMOS_,PMOS_,NMOS_,NMOS_,PMOS_ };
	double M_VT_[2 * NumTrans] = { 0.1750,0.1750,-0.2640, 0.1750,0.1750,-0.2640 };
	double M_W_[2*NumTrans] = { 0.16E-6,0.07E-6,0.08E-6,0.16E-6,0.07E-6,0.08E-6 };
	double M_L_[2 * NumTrans] = { 0.6E-7 ,0.6E-7,0.6E-7,0.6E-7 ,0.6E-7,0.6E-7 };
	double M_MU_[2 * NumTrans] = { 0.0204,0.0204,.0090 ,0.0204,0.0204,.0090 };
	double M_COX_[2 * NumTrans] = { 0.013595,0.013595,0.01298170,0.013595,0.013595,0.01298170 };
	double M_LAMBDA_[2*NumTrans] = { 0,0,0,0,0,0 };
	double M_COV_[2 * NumTrans] = {1.79454E-4,1.79454E-4, 1.6097308E-4,1.79454E-4,1.79454E-4,1.6097308E-4 };
	//double M_COV_[2 * NumTrans] = { 1.29454E-4,1.29454E-4, 1.1097308E-4,1.29454E-4,1.29454E-4,1.1097308E-4 };
	//double M_COV_[2 * NumTrans] = { 0 };
	bool dcsim(vector<vector<double> > varValue, vector<double> &Res_iPG, vector<double>&Res_V2, bool dis = 0);
	bool readsim(vector<vector<double> > &varValue,vector<double> &delta_V);
	bool writesim(vector<vector<double> > &varValue, vector<double> &delta_t);
	bool MC(vector<int> &src,vector<double> &dst,int simtype);
	bool displayall();
	bool addcorelation(int numsample,vector<double> &core_V);
private:
	bool display(int *array, string name);
	bool display(double *array, string name);
	bool stamp_MOSFET(vector<vector<double> > &M, vector<double>&J, vector<double> &old_V, int MOS_Num);
	bool alter(vector<double> &paramMean, vector<vector<double> >&varValue, vector<vector<double> >&Res_paramValue,int simtype);
	bool set_M_VT_(vector<double> src);
	bool set_M_VT_all(vector<double> src);
	bool reset_M_VT_all();
	bool simple_solve(vector<vector<double> >cur_M, vector<double>cur_J, vector<double> &V); // M[3,3] * V[3,1] = J[3,1] solve V
	double Icurrent(int MOS_num, int MOS_type, double Vg, double Vd, double Vs);
	double Capitance(int MOS_num, int MOS_type, double Vg, double Vd, double Vs, int C_type);
	double cal_integral(double lower, double upper, int intervals=1000);
	double helper_integral(double cur_v);
};
#endif 

