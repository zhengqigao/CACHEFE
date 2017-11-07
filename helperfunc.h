#ifndef __HELPERFUNC__H__
#define __HELPERFUNC__H__

#include"const.h"

namespace helperfunc{
	int combntns(int maxx, int minn);
	bool choose_ref(vector<double>&src, vector<double> &dst, vector<bool>ref, int axis);
	bool choose_ref(vector<vector<double> >&src, vector<vector<double> > &dst, vector<bool>ref, int axis);
	bool display_matrix_vector(vector<vector<double> >src, string name, int precision = 6);
	bool display_matrix_vector(vector<double>src, string name, int precision = 6);
	bool delete_ref(vector<double>&src, vector<bool>ref, int axis = 1);
	bool delete_ref(vector<vector<double> >&src, vector<bool>ref, int axis = 1);
	bool expand(vector<vector<double> >&src1, vector<vector<double> >&src2, vector<vector<double> >&dst, int axis);
	double mean(vector<double>src);
	double var(vector<double> src);
	void normcdf(double *src, int size, double mu, double sigma, double *dst);
	void normpdf(double *src, int size, double mu, double sigma, double *dst);
	double erf(double x);
	bool write_to_txt(vector<double> &Fp, const char *file);
	bool write_to_txt(vector<vector<double> >&Fp, const char *file);
	bool Mysolve(double A, double B, double C, double *dst);
	bool write_res(vector<int> &order,vector<double> &res,const char *file);
}

#endif 
