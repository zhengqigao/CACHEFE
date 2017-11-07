#include"demo.h"

bool demo::safail(){

	system("mkdir -p ./CACHEFE_res");
	// get parameters ready
	double par_t_C = 1019.6867251526046;
	int par_nCend = 8, par_randseed = 7, par_nsimiter = 10000;
	double *par_u01 = new double[2], *par_sig01 = new double[2];
	par_u01[0] = 0; par_u01[1] = INFINITY;
	par_sig01[0] = 2.451733398437500e-05; par_sig01[1] = 2.451733398437500e-05;

	// initialzation
	SUS SUS_instance(par_nCend, par_u01, par_sig01, par_t_C, par_randseed, par_nsimiter);

	// let's do simulation !
	vector<vector<double> >XSeed, Yseed, XsaSeed, ylimSeed;
	SUS_instance.sus_delta_sim(XSeed, Yseed, XsaSeed, ylimSeed,SAFAIL_);
	vector<int>cells={1,2,3,4,5,6,7,8,10,13,16,32,48,64};
	// get results
	double prob = 1;
	vector<double>probAnd(SUS_instance.nCend);
	vector<vector<double> >Fp(par_nCend, vector<double>(cells.size()));
	for (int j = 1; j <= par_nCend; j++) {
		prob *= SUS_instance.APA_probEstList[j - 1];
		//cout<<SUS_instance.APA_probEstList[j-1]<<endl;
		probAnd[j - 1] = prob;
		// APA result
		for (int m = 0; m < cells.size(); m++) {
			Fp[j - 1][m] = (cells[m] - j + 1)*pow(-1, j - 1)*probAnd[j - 1];
		}
		for (int k = 1; k <= j - 1; k++) {
			for (int m = 0; m < Fp[0].size(); m++) {
				Fp[j - 1][m] = Fp[j - 1][m] + pow(-1, k - 1)*probAnd[k - 1] * (helperfunc::combntns(j, k) + (cells[m] - j)*helperfunc::combntns(j - 1, k - 1));
			
			}
		}
	}
	helperfunc::write_to_txt(Fp, "./CACHEFE_res/APAdemo_safail__res.txt");

	delete []par_u01;
	delete []par_sig01;
	return 1;
}



bool demo::readfail() {
	
	system("mkdir -p ./CACHEFE_res");
	// get parameters ready
	double par_t_C = 1019.6867251526046;
	int par_nCend = 8, par_randseed = 7, par_nsimiter = 10000;
	double *par_u01 = new double[2], *par_sig01 = new double[2];
	par_u01[0] = 0; par_u01[1] = INFINITY;
	par_sig01[0] = 2.451733398437500e-05; par_sig01[1] = 2.451733398437500e-05;

	// initialzation
	SUS SUS_instance(par_nCend, par_u01, par_sig01, par_t_C, par_randseed, par_nsimiter);

	// let's do simulation !
	vector<vector<double> >XSeed, Yseed, XsaSeed, ylimSeed;
	SUS_instance.sus_delta_sim(XSeed, Yseed, XsaSeed, ylimSeed,READFAIL_);
	vector<int>cells={1,2,3,4,5,6,7,8,10,13,16,32,48,64};
	// get results
	double prob = 1;
	vector<double>probAnd(SUS_instance.nCend);
	vector<vector<double> >Fp(par_nCend, vector<double>(cells.size()));
	for (int j = 1; j <= par_nCend; j++) {
		prob *= SUS_instance.APA_probEstList[j - 1];
		probAnd[j - 1] = prob;
		// APA result
		for (int m = 0; m < cells.size(); m++) {
			Fp[j - 1][m] = (cells[m] - j + 1)*pow(-1, j - 1)*probAnd[j - 1];
		}
		for (int k = 1; k <= j - 1; k++) {
			for (int m = 0; m < Fp[0].size(); m++) {
				Fp[j - 1][m] = Fp[j - 1][m] + pow(-1, k - 1)*probAnd[k - 1] * (helperfunc::combntns(j, k) + (cells[m] - j)*helperfunc::combntns(j - 1, k - 1));
			
			}
		}
	}
	delete []par_u01;
	delete []par_sig01;
	helperfunc::write_to_txt(Fp, "./CACHEFE_res/APAdemo_readfail__res.txt");
	return 1;
	
}

bool demo::writefail(){
	system("mkdir -p ./CACHEFE_res");
	// get parameters ready
	double par_t_C = 1019.6867251526046;
	int par_nCend = 8, par_randseed = 7, par_nsimiter = 10000;
	double *par_u01 = new double[2], *par_sig01 = new double[2];
	par_u01[0] = 0; par_u01[1] = INFINITY;
	par_sig01[0] = 2.451733398437500e-05; par_sig01[1] = 2.451733398437500e-05;

	// initialzation
	SUS SUS_instance(par_nCend, par_u01, par_sig01, par_t_C, par_randseed, par_nsimiter);

	// let's do simulation !
	vector<vector<double> >XSeed, Yseed, XsaSeed, ylimSeed;
	SUS_instance.sus_delta_sim(XSeed, Yseed, XsaSeed, ylimSeed,WRITEFAIL_);
	vector<int>cells={1,2,3,4,5,6,7,8,10,13,16,32,48,64};
	// get results
	double prob = 1;
	vector<double>probAnd(SUS_instance.nCend);
	vector<vector<double> >Fp(par_nCend, vector<double>(cells.size()));
	for (int j = 1; j <= par_nCend; j++) {
		prob *= SUS_instance.APA_probEstList[j - 1];
		probAnd[j - 1] = prob;
		// APA result
		for (int m = 0; m < cells.size(); m++) {
			Fp[j - 1][m] = (cells[m] - j + 1)*pow(-1, j - 1)*probAnd[j - 1];
		}
		for (int k = 1; k <= j - 1; k++) {
			for (int m = 0; m < Fp[0].size(); m++) {
				Fp[j - 1][m] = Fp[j - 1][m] + pow(-1, k - 1)*probAnd[k - 1] * (helperfunc::combntns(j, k) + (cells[m] - j)*helperfunc::combntns(j - 1, k - 1));
			
			}
		}
	}
	delete []par_u01;
	delete []par_sig01;
	helperfunc::write_to_txt(Fp, "./CACHEFE_res/APAdemo_writefail__res.txt");
	return 1;
}


bool demo::cleardemo(){
	system("rm -r ./CACHEFE");
	return 1;
}