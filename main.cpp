#include"APA.h"
#include"cmdline.h"
int APAdemo();

int main(int argc,char *argv[]) {
	cmdline::parser Myparser;
	Myparser.add<string>("Algo",'A',"Algorithm used by CACHEFE",false,"APA",cmdline::oneof<string>("APA","APE"));
	Myparser.add<string>("Type",'T',"Failure Type simulated by CACHEFE",false,"All",cmdline::oneof<string>("all","read","write","sa"));
	Myparser.add<int>("Cell",'C',"cell number of SRAM simulated by CACHEFE",true,0,cmdline::range(1, 65));
	Myparser.parse_check(argc, argv);
	string Algo = Myparser.get<string>("Algo");
	string simtype = Myparser.get<string>("Type");
	int Cellnum = Myparser.get<int>("Cell");
	std::vector<int> Simtypelist;
	

	transform(simtype.begin(), simtype.end(), simtype.begin(), ::tolower);  
	if (simtype=="all" || simtype=="read"){
		Simtypelist.push_back(READFAIL_);
	}
	if (simtype=="all" || simtype=="write"){
		Simtypelist.push_back(WRITEFAIL_);
	}
	if (simtype=="all" || simtype=="sa"){
		Simtypelist.push_back(SAFAIL_);
	}

	cout << "\n----------CACHEFE for SRAM Failure Rate Estimation ----------\n\n\n";

	for (int i=0;i<Simtypelist.size();i++){
		// get parameters ready
		double par_t_C = 1019.6867251526046;
		int par_nCend = 8, par_randseed = 7, par_nsimiter = 1000;
		double *par_u01 = new double[2], *par_sig01 = new double[2];
		par_u01[0] = 0; par_u01[1] = INFINITY;
		par_sig01[0] = 2.451733398437500e-05; par_sig01[1] = 2.451733398437500e-05;

		// initialzation
		SUS SUS_instance(par_nCend, par_u01, par_sig01, par_t_C, par_randseed, par_nsimiter);

		// let's do simulation !
		vector<vector<double> >XSeed, Yseed, XsaSeed, ylimSeed;

		string cur;

		switch (Simtypelist[i]){
			case READFAIL_:
				cur="Read";
				break;
			case WRITEFAIL_:
				cur="write";
				break;
			case SAFAIL_:
				cur="sense amplifier";
				break;
			default:
				cout<<"[Error]:wrong simulatinon type for CACHFE\n";
				exit(1);
		}

		string wrk1_pre="\nSRAM 【";
		string wrk1_post=" Failure Rate】 Estimation...... \n\n";
		wrk1_pre.append(cur);
		wrk1_pre.append(wrk1_post);
		cout << wrk1_pre;

		SUS_instance.sus_delta_sim(XSeed, Yseed, XsaSeed, ylimSeed,Simtypelist[i]);
		// get results
		double prob = 1;
		//string tmp; helperfunc::display_matrix_vector(SUS_instance.APA_probEstList, tmp);
		vector<double>probAnd(SUS_instance.nCend);
		vector<vector<double> >Fp(par_nCend, vector<double>(1));
		for (int j = 1; j <= par_nCend; j++) {
			prob *= SUS_instance.APA_probEstList[j - 1];
			probAnd[j - 1] = prob;
			// APA result
			for (int m = 0; m < 1; m++) {
				Fp[j - 1][m] = (Cellnum - j + 1)*pow(-1, j - 1)*probAnd[j - 1];
			}
			for (int k = 1; k <= j - 1; k++) {
				for (int m = 0; m < Fp[0].size(); m++) {
					Fp[j - 1][m] = Fp[j - 1][m] + pow(-1, k - 1)*probAnd[k - 1] * (helperfunc::combntns(j, k) + (Cellnum - j)*helperfunc::combntns(j - 1, k - 1));
				}
			}
		}

		std::vector<double> res;
		cout<<endl;
		for (int i = 0, cur = 0; i<Fp.size(); i++) {
			if (Fp[i][0]>0) {
				res.push_back(Fp[i][0]);
				cur++;
				cout << "APA order " << cur << "th Failure Rate Estimation Result: " << Fp[i][0] << endl;
			}
		}
		cout << endl;
		string wrk2_pre="./APA_";
		string wrk2_post="fail_res.txt";
		wrk2_pre.append(cur);
		wrk2_pre.append(wrk2_post);
   		const char* file = wrk2_pre.c_str();
		helperfunc::write_to_txt(res, file);
	}
	cout << "CACHEFE finish Failure Rate Estimation...Press any key to return...\n";
	system("read");
	return 1;
}


int APAdemo() {
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
	helperfunc::write_to_txt(Fp, "./APAdemo__res.txt");
	return 1;
}