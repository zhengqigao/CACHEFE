#include"APA.h"
#include"cmdline.h"
#include"demo.h"

int main(int argc,char *argv[]) {
	
	// Yu can uncomment the lines below to see demo results. The following lines will create a new directory ./CACHEFE_res/ 
	// And it will put the demo results in the directory.

	//demo::writefail();
	//cout<<"finish write"<<endl;
	//demo::readfail();
	//cout<<"finish read"<<endl;
	//demo::safail();
	//cout<<"finish sa"<<endl<<endl;


	// parse Input 
	cmdline::parser Myparser;
	Myparser.add<string>("Algo",'A',"Algorithm used by CACHEFE",false,"APA",cmdline::oneof<string>("APA","APE"));
	Myparser.add<string>("Type",'T',"Failure Type simulated by CACHEFE",false,"all");
	Myparser.add<int>("Cell",'C',"cell number of SRAM simulated by CACHEFE",true,0,cmdline::range(1, 65));
	Myparser.parse_check(argc, argv);
	
	string Algo = Myparser.get<string>("Algo");
	string simtype = Myparser.get<string>("Type");
	transform(simtype.begin(), simtype.end(), simtype.begin(), ::tolower); 
	int Cellnum = Myparser.get<int>("Cell");
	std::vector<int> Simtypelist;
	if (simtype != "all" && simtype !="read" && simtype != "write" && simtype !="sa"){
		cout<<"[Error]: Simulation type is wrong, must bt one of ['all','read','write','sa']"<<endl;
		exit(1);
	}

	if (simtype=="all" || simtype=="read"){
		Simtypelist.push_back(READFAIL_);
	}
	if (simtype=="all" || simtype=="write"){
		Simtypelist.push_back(WRITEFAIL_);
	}
	if (simtype=="all" || simtype=="sa"){
		Simtypelist.push_back(SAFAIL_);
	}

	//begin the main part of the tool

	cout << "\n----------CACHEFE for SRAM Failure Rate Estimation ----------\n\n\n";

	//begin simulating very type of failure
	for (int i=0;i<Simtypelist.size();i++){
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
		string wrk1_post=" Failure Rate】 Estimation :  \n\n";
		wrk1_pre.append(cur);
		wrk1_pre.append(wrk1_post);
		cout << wrk1_pre;

		// get parameters ready
		double par_t_C = 1019.6867251526046;
		int par_nCend = 8, par_randseed = 7, par_nsimiter = 1000;
		double *par_u01 = new double[2], *par_sig01 = new double[2];
		par_u01[0] = 0; par_u01[1] = INFINITY;
		par_sig01[0] = 2.451733398437500e-05; par_sig01[1] = 2.451733398437500e-05;

		// initialzation
		SUS SUS_instance(par_nCend, par_u01, par_sig01, par_t_C, par_randseed, par_nsimiter);


		SUS_instance.sus_delta_sim(XSeed, Yseed, XsaSeed, ylimSeed,Simtypelist[i]);
		// get results
		double prob = 1;
		vector<double>probAnd(SUS_instance.nCend);
		vector<double> Fp(par_nCend);// very order estimation results, may be negative
		for (int j = 1; j <= par_nCend; j++) {
			prob *= SUS_instance.APA_probEstList[j - 1];
			probAnd[j - 1] = prob;
			// APA result
			Fp[j - 1] = (Cellnum - j + 1)*pow(-1, j - 1)*probAnd[j - 1];
			for (int k = 1; k <= j - 1; k++) {
				Fp[j - 1] = Fp[j - 1] + pow(-1, k - 1)*probAnd[k - 1] * (helperfunc::combntns(j, k) + (Cellnum - j)*helperfunc::combntns(j - 1, k - 1));
			}
		}
		std::vector<double> res;
		double lowerbound = MYMAX,upperbound = MYMIN; //the range of the final result
		for (int i=0,cur=0;i<Fp.size();i++){
			if (Fp[i]>0 && Fp[i]<1){
				lowerbound = (lowerbound>Fp[i]?Fp[i]:lowerbound);
				upperbound = (upperbound>Fp[i]?upperbound:Fp[i]);
				res.push_back(Fp[i]);
				cur++;
				cout<<cur<<"th Estimation Result: "<<Fp[i]<<endl;
			}
		}
		//cout<<endl<<"\t"<<cur<<" Failure Rate Result: "<<"[ "<<lowerbound<<" , "<<upperbound<<"]"<<endl;
		cout << endl;
	
	}
	cout << "CACHEFE finish Failure Rate Estimation...Press any key to return...\n";
	system("read");
	return 1;
}