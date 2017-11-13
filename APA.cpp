#include"APA.h"
#include"Simulator.h"

SUS::SUS(int par_nCend, double *par_u01, double *par_sig01, double par_t_C, int par_randseed, int par_nSimiter) :e(par_randseed) {
	nCend = par_nCend;
	u01[0] = par_u01[0]; u01[1] = par_u01[1];
	sig01[0] = par_sig01[0]; sig01[1] = par_sig01[1];
	t_C = par_t_C;
	randseed = par_randseed;
	nSimiter = par_nSimiter;
	/*
	cout<<par_nCend<<endl;
	cout<<par_u01[0]<<" "<<par_u01[1]<<endl;
	cout<<par_sig01[0]<<" "<<par_sig01[1]<<endl;
	cout<<par_t_C<<endl;
	cout<<par_randseed<<endl;
	cout<<par_nSmimiter<<endl;
	*/
}

SUS::~SUS() {};

bool SUS::sus_delta_sim(vector<vector<double> >&XSeed, vector<vector<double> >&YSeed, vector<vector<double> >&XsaSeed, vector<vector<double> >&ylimSeed, int simtype) {
	// setup
	int nSimTotal = 0, nIter = 0;
	vector<double>probList, sigList;
	vector<vector<double> >perfDeltaList(2, vector<double>());
	vector<double> perfDelta(2); perfDelta[0] = 0; perfDelta[1] = 0;
	int nC = YSeed.size(), nS;
	if (nC != 0) { nS = YSeed[0].size(); } // the first time Yseed will be an empty matrix,thus nC==0,nS==0
	vector<vector<double> > ylim(2, vector<double>(nSimiter));
	vector <vector<double> >Xsa(2, vector<double>(nSimiter));
	//default_random_engine e(randseed);
	//normal_distribution<double> n(0, 1);
	if (nC == 0) {
		nS = nSimiter;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < nSimiter; j++) {
				Xsa[i][j] = n(e);
			}
		}
		vector<double>ylimit_V(nSimiter);
		sim_SA_fake(Xsa, ylimit_V, simtype);
		for (int i = 0; i < ylim[0].size(); i++) {
			ylim[0][i] = ylimit_V[i] / t_C + u01[0];
			ylim[1][i] = sig01[1] * n(e) + u01[1];
			//cout << ylim[0][i] << "\t" << ylim[1][i] << endl;
		}
	}
	else {
		Xsa = XsaSeed;
		ylim = ylimSeed;
		if (YSeed[0].size()<nSimiter) {
			vector<vector<vector<double> > >Xcell, Ycell;
			expandSeed(XSeed, YSeed, Xsa, ylimSeed, nSimTotal, nSimiter, XSeed, YSeed, Xsa, ylim, Xcell, Ycell, simtype);
		}
		nC = YSeed.size();
		nS = YSeed[0].size();
	}
	vector<vector<double> >Xnew(2 * NumTrans, vector<double>(nS));
	for (int i = 0; i < Xnew.size(); i++) {
		for (int j = 0; j < Xnew[0].size(); j++) {
			Xnew[i][j] = n(e);
		}
	}
	//cout<<"IN sus_delta_sim "<<Xnew.size()<<"\t"<<Xnew[0].size()<<endl;
	//system("read");
	vector<vector<double> >Ynew(Xnew.size() / (2 * NumTrans), vector<double>(nS));// Ynew is the simulation result of Xnew, 6 transistors will be a cell
	simout(Xnew, Ynew, nSimTotal, 0, simtype);// here because of epo==0, Ynew will be a vector
	vector<vector<double> >X(XSeed.size() + Xnew.size(), vector<double>(Xnew[0].size()));
	vector<vector<double> >Y(YSeed.size() + Ynew.size(), vector<double>(Ynew[0].size()));
	// obtain X,Y
	helperfunc::expand(XSeed, Xnew, X, 2);
	helperfunc::expand(YSeed, Ynew, Y, 2);

	nIter = nIter + 1;
	vector<bool>indSimFail(nS); // here beacuse Ynew is actually a vector...
	double probSimFail = 0;
	for (int i = 0; i < nS; i++) {
		indSimFail[i] = isnan(Ynew[0][i]);
		probSimFail += indSimFail[i];
	}
	probSimFail /= (double(nS));

	helperfunc::delete_ref(Ynew, indSimFail, 2);
	helperfunc::delete_ref(X, indSimFail, 2);
	helperfunc::delete_ref(Y, indSimFail, 2);
	helperfunc::delete_ref(Xsa, indSimFail, 2);
	helperfunc::delete_ref(ylim, indSimFail, 2);

	vector<double>Y1(Ynew[0].size()), Y2(Ynew[0].size());
	for (int i = 0; i < Y1.size(); i++) {
		Y1[i] = Ynew[0][i] - ylim[0][i];
		Y2[i] = Ynew[0][i] - ylim[1][i];
	}

	vector<bool>indParamFail(Ynew[0].size()), wrk_ref(Ynew[0].size());
	double probParamFail = 0;
	for (int i = 0; i < indParamFail.size(); i++) {
		indParamFail[i] = (Y1[i]<perfDelta[0] || Y2[i]>perfDelta[1]);
		wrk_ref[i] = !indParamFail[i];// !  ~  CAUTION
		probParamFail += indParamFail[i];
	}
	probParamFail /= (double(indParamFail.size()));
	double probCur = probParamFail*(1 - probSimFail) + probSimFail;
	if (probCur >= PROBSUB) {
		double sigCur = sqrt((1 - probCur)*probCur / nS);
		vector<double>perfDeltaCur = perfDelta;
		perfDeltaList.push_back(perfDeltaCur);
		probList.push_back(probCur);
		sigList.push_back(sigCur);
		cout << "\tFinish " << nC + 1 << "th order APA......current order simulation times : " << nSimTotal << endl;
		output(probList, sigList, nSimTotal);
		if (nC + 1 < nCend) {
			//next APA order
			helperfunc::delete_ref(X, wrk_ref, 2);
			helperfunc::delete_ref(Y, wrk_ref, 2);
			helperfunc::delete_ref(Xsa, wrk_ref, 2);
			helperfunc::delete_ref(ylim, wrk_ref, 2);
			sus_delta_sim(X, Y, Xsa, ylim, simtype);
		}
	}
	else {
		double probTarget = (PROBSUB - probSimFail) / (1 - probSimFail);

		vector<double>perfDeltaCur(2);
		getspec(Y1, Y2, probTarget, perfDeltaCur);

		vector<bool>indParamFail(Ynew[0].size());
		double probParamFail = 0;
		for (int i = 0; i < indParamFail.size(); i++) {
			indParamFail[i] = (Y1[i]<perfDeltaCur[0] || Y2[i]>perfDeltaCur[1]);
			probParamFail += indParamFail[i];
		}
		probParamFail /= (double(indParamFail.size()));
		double probCur = probParamFail*(1 - probSimFail) + probSimFail;
		double sigCur = sqrt((1 - probCur)*probCur / nS);
		vector<double>Y1Seed, Y2Seed;
		helperfunc::choose_ref(X, XSeed, indParamFail, 2);
		helperfunc::choose_ref(Y, YSeed, indParamFail, 2);
		helperfunc::choose_ref(Y1, Y1Seed, indParamFail, 2);
		helperfunc::choose_ref(Y2, Y2Seed, indParamFail, 2);
		helperfunc::choose_ref(Xsa, XsaSeed, indParamFail, 2);
		helperfunc::choose_ref(ylim, ylimSeed, indParamFail, 2);
		probList.push_back(probCur);
		sigList.push_back(sigCur);
		perfDeltaList.push_back(perfDeltaCur);
		while (nIter <= NITERMAX) {
			vector<vector<vector<double> > >X, Y, Y1, Y2, Xsa, ylim;
			genX(XSeed, YSeed, Y1Seed, Y2Seed, XsaSeed, ylimSeed, nSimTotal, nSimiter, perfDeltaList[perfDeltaList.size() - 1], X, Y, Y1, Y2, Xsa, ylim, simtype);
			nIter = nIter + 1;
			int eop = 0;
			filterX(X, Y, Y1, Y2, Xsa, ylim, probList, sigList, perfDeltaList, nSimTotal, nSimiter, XSeed, YSeed, Y1Seed, Y2Seed, XsaSeed, ylimSeed, eop, simtype);
			if (eop == 1) {
				break;
			}
		}
	}
	return 1;
}

bool SUS::genX(vector<vector<double> >&XSeed, vector<vector<double> >&YSeed, vector<double>&Y1Seed, vector<double>&Y2Seed, vector<vector<double> >&XsaSeed, vector<vector<double> >&ylimSeed, int &nSimTotal, int  nSimiter, vector<double> perfSpec, vector<vector<vector<double> > >&Res_X, vector<vector<vector<double> > > &Res_Y, vector<vector<vector<double> > >&Res_Y1, vector<vector<vector<double> > >&Res_Y2, vector<vector<vector<double> > >&Res_Xsa, vector<vector<vector<double> > >&Res_ylim, int simtype) {
	int nC = YSeed.size(), nSeed = YSeed[0].size();
	vector<vector<double> >XsapSeed(XsaSeed);
	vector<vector<double> >ylimpSeed(ylimSeed);
	vector<vector<vector<double> > >Xcell, Ycell;
	if (nC >= 2) {
		vector<vector<double> >XpSeed(XSeed.size() - 6, vector<double>(XSeed[0].size()));
		vector<bool> ref(XSeed.size());
		for (int i = 0; i < XSeed.size() - 6; i++) {
			ref[i] = 1;
		}
		helperfunc::choose_ref(XSeed, XpSeed, ref, 1);
		helperfunc::delete_ref(XSeed, ref, 1);
		vector<vector<double> >YpSeed(YSeed.size() - 1, vector<double>(YSeed[0].size()));
		vector<bool>ref2(YSeed.size());
		for (int i = 0; i < YSeed.size() - 1; i++) {
			ref2[i] = 1;
		}
		helperfunc::choose_ref(YSeed, YpSeed, ref2, 1);
		helperfunc::delete_ref(YSeed, ref2, 1);
		expandSeed(XpSeed, YpSeed, XsapSeed, ylimpSeed, nSimTotal, nSimiter, XpSeed, YpSeed, XsapSeed, ylimpSeed, Xcell, Ycell, simtype);
	}
	vector<vector<vector<double> > >X(nSeed), Y(nSeed), Y1(nSeed), Y2(nSeed), Xsa(nSeed), ylim(nSeed);
	int nSim = nSeed, nAccept = 0;
	for (int i = 0; i < nSeed; i++) {
		// store the ith chain in X[i]
		vector<double>tmp1(XSeed.size());
		for (int j = 0; j < tmp1.size(); j++) {
			tmp1[j] = XSeed[j][i];
		}
		X[i].push_back(tmp1);

		vector<double>tmp2(YSeed.size());
		for (int j = 0; j < tmp2.size(); j++) {
			tmp2[j] = YSeed[j][i];
		}
		Y[i].push_back(tmp2);

		vector<double>tmp3(XsaSeed.size());
		for (int j = 0; j < tmp3.size(); j++) {
			tmp3[j] = XsaSeed[j][i];
		}
		Xsa[i].push_back(tmp3);

		vector<double>tmp4(ylimSeed.size());
		for (int j = 0; j < tmp4.size(); j++) {
			tmp4[j] = ylimSeed[j][0];
		}
		ylim[i].push_back(tmp4);

		vector<double>tmp5(1);
		for (int j = 0; j < tmp5.size(); j++) {
			tmp5[j] = Y1Seed[i];
		}
		Y1[i].push_back(tmp5);

		vector<double>tmp6(1);
		for (int j = 0; j < tmp6.size(); j++) {
			tmp6[j] = Y2Seed[i];
		}
		Y2[i].push_back(tmp6);
	}
	while (nSim < nSimiter) {
		for (int i = 0; i < nSeed; i++) {
			vector<double> wrk(X[i][X[i].size() - 1]);
			vector<vector<double> >XTemp((X[i][X[i].size() - 1]).size(), vector<double>(1));
			metropolis(wrk, wrk.size(), XTemp);
			vector<vector<double> >YTemp(XTemp.size() / (2 * NumTrans), vector<double>(1));
			simout(XTemp, YTemp, nSimTotal, 0, simtype);
			vector<vector<double> >XsaTemp((Xsa[i][Xsa[i].size() - 1]).size(), vector<double>(1));
			vector<double> ylimTemp(2);
			if (nC < 2) {
				vector<double> wrk2(Xsa[i][Xsa[i].size() - 1]);
				metropolis(wrk2, 2, XsaTemp);
				vector<double>ylimit_V(1);
				sim_SA_fake(XsaTemp, ylimit_V, simtype);
				ylimTemp[0] = ylimit_V[0] / t_C + u01[0];
				ylimTemp[1] = ylimit_V[0] / t_C + u01[1];
			}
			else {
				for (int ii = 0; ii < XsaSeed.size(); ii++) {
					XsaTemp[ii][0] = XsapSeed[ii][nSim];
				}
				for (int ii = 0; ii < ylimpSeed.size(); ii++) {
					ylimTemp[ii] = ylimpSeed[ii][nSim];
				}
			}
			vector<vector<double> >Y1Temp(YTemp), Y2Temp(YTemp);
			double max1 = MYMIN, min2 = MYMAX;
			int hasnan = 0;
			for (int i = 0; i < YTemp.size(); i++) {
				hasnan += isnan(Y1Temp[i][0]);
				Y1Temp[i][0] -= ylimTemp[0];
				Y2Temp[i][0] -= ylimTemp[1];
				max1 = max1 > Y1Temp[i][0] ? max1 : Y1Temp[i][0];
				min2 = min2 < Y2Temp[i][0] ? min2 : Y2Temp[i][0];
			}
			if (hasnan > 0 || max1 < perfSpec[0] || min2>perfSpec[1]) {
				vector<double>wrk_XTemp(XTemp.size());
				transpose(XTemp, wrk_XTemp);
				X[i].push_back(wrk_XTemp);

				vector<double>wrk_YTemp(YTemp.size());
				transpose(YTemp, wrk_YTemp);
				Y[i].push_back(wrk_YTemp);

				vector<double>wrk_XsaTemp(XsaTemp.size());
				transpose(XsaTemp, wrk_XsaTemp);
				Xsa[i].push_back(wrk_XsaTemp);

				vector<double>wrk_Y1Temp(Y1Temp.size());
				transpose(Y1Temp, wrk_Y1Temp);
				Y1[i].push_back(wrk_Y1Temp);

				vector<double>wrk_Y2Temp(Y2Temp.size());
				transpose(Y2Temp, wrk_Y2Temp);
				Y2[i].push_back(wrk_Y2Temp);

				ylim[i].push_back(ylimTemp);
				nAccept += 1;
			}
			else {
				X[i].push_back(X[i][X[i].size() - 1]);
				Y[i].push_back(Y[i][Y[i].size() - 1]);
				if (nC > 2) {
					Xcell[i][X[i].size() - 1] = Xcell[i][X[i].size() - 2];
					Ycell[i][Y[i].size() - 1] = Ycell[i][Y[i].size() - 2];
				}
				Xsa[i].push_back(Xsa[i][Xsa[i].size() - 1]);
				ylim[i].push_back(ylim[i][ylim[i].size() - 1]);
				Y1[i].push_back(Y1[i][Y1[i].size() - 1]);
				Y2[i].push_back(Y2[i][Y2[i].size() - 1]);
			}
			nSim += 1;
			if (nSim >= nSimiter) {
				break;
			}
		}
	}
	vector<vector<vector<double> > >tmp_Res_X(nSeed), tmp_Res_Y(nSeed);
	if (nC >= 2) {
		for (int i = 0; i < nSeed; i++) {
			for (int j = 0; j < Xcell[i].size(); j++) {
				vector<double>wrk(Xcell[i][j].size() + X[i][j].size());
				for (int m = 0; m < Xcell[i][j].size() + X[i][j].size(); m++) {
					wrk[m] = m < Xcell[i][j].size() ? Xcell[i][j][m] : X[i][j][m - Xcell[i][j].size()];
				}
				tmp_Res_X[i].push_back(wrk);
			}
		}
		for (int i = 0; i < nSeed; i++) {
			for (int j = 0; j < Ycell[i].size(); j++) {
				vector<double>wrk(Ycell[i][j].size() + Y[i][j].size());
				for (int m = 0; m < Ycell[i][j].size() + Y[i][j].size(); m++) {
					wrk[m] = m < Ycell[i][j].size() ? Ycell[i][j][m] : Y[i][j][m - Ycell[i][j].size()];
				}
				tmp_Res_Y[i].push_back(wrk);
			}
		}
		Res_X = tmp_Res_X;
		Res_Y = tmp_Res_Y;
	}
	else {
		Res_X = X;
		Res_Y = Y;
	}
	Res_Y1 = Y1;
	Res_Y2 = Y2;
	Res_Xsa = Xsa;
	Res_ylim = ylim;
	return 1;
}


bool SUS::filterX(vector<vector<vector<double> > > &X, vector<vector<vector<double> > >&Y, vector<vector<vector<double> > >&Y1, vector<vector<vector<double> > >&Y2, vector<vector<vector<double> > >&Xsa, vector<vector<vector<double> > >&ylim, vector<double>&probList, vector<double>&sigList, vector<vector<double> >&perfDeltaList, int &nSimTotal, int nSimiter, vector<vector<double> > &Res_XSeed, vector<vector<double> > &Res_YSeed, vector<double>&Res_Y1Seed, vector<double>&Res_Y2Seed, vector<vector<double> >&Res_XsaSeed, vector<vector<double> >&Res_ylimSeed, int &Res_eop, int simtype) {
	Res_eop = 0;
	int nSeed = Y.size();
	vector<vector<vector<double> > >indSimFail(nSeed), indParamFail(nSeed), indFail(nSeed);
	vector<double> nSim(nSeed);
	vector<double>perfDelta(2);
	vector<vector<double> >XAll, YAll, XsaAll, ylimAll;
	vector<double> Y1All, Y2All;
	vector<double>indFailAll, indSimFailAll;
	for (int i = 0; i < nSeed; i++) {
		//cout << i << endl;
		nSim[i] = Y1[i].size();
		vector<double>wrk(Y[i].size());
		for (int p = 0; p < wrk.size(); p++) {
			for (int j = 0; j < Y[i][p].size(); j++) {
				if (isnan(Y[i][p][j])) {
					wrk[p] = true;
					break;
				}
			}
		}
		indSimFail[i].push_back(wrk);
		vector<double>wrk2(Y1[i].size()), wrk3(Y1[i].size());
		for (int j = 0; j < wrk2.size(); j++) {
			wrk2[j] = Y1[i][j][0] < perfDelta[0] || Y2[i][j][0] > perfDelta[1];
			wrk3[j] = wrk2[j] || indSimFail[i][0][j];
		}
		indParamFail[i].push_back(wrk2);
		indFail[i].push_back(wrk3);

		for (int j = 0; j < wrk.size(); j++) {
			indSimFailAll.push_back(wrk[j]);
		}
		for (int j = 0; j < wrk3.size(); j++) {
			indFailAll.push_back(wrk3[j]);
		}
		for (int j = 0; j < X[i].size(); j++) {
			XAll.push_back(X[i][j]);
		}
		for (int j = 0; j < Y[i].size(); j++) {
			YAll.push_back(Y[i][j]);
		}
		for (int j = 0; j < Y1[i].size(); j++) {
			Y1All.push_back(Y1[i][j][0]);
		}
		for (int j = 0; j < Y2[i].size(); j++) {
			Y2All.push_back(Y2[i][j][0]);
		}
		for (int j = 0; j < Xsa[i].size(); j++) {
			XsaAll.push_back(Xsa[i][j]);
		}
		for (int j = 0; j < ylim[i].size(); j++) {
			ylimAll.push_back(ylim[i][j]);
		}
	}
	XAll = transpose(XAll);
	YAll = transpose(YAll);
	XsaAll = transpose(XsaAll);
	ylimAll = transpose(ylimAll);
	double probCur = 0;
	for (int i = 0; i < indFailAll.size(); i++) {
		probCur += indFailAll[i];
	}
	probCur /= (double(indFailAll.size()));
	vector<double>perfDeltaCur(2);
	if (probCur >= PROBSUB) {
		Res_eop = 1;
		vector<double>perfDeltCur(perfDelta);
		// first find indFailAll==1
		vector<bool>wrk_ref(indFailAll.size());
		for (int i = 0; i < indFailAll.size(); i++) {
			bool tmp = (indFailAll[i] == 1 ? 1 : 0);
			wrk_ref[i] = tmp;
		}
		helperfunc::choose_ref(XAll, Res_XSeed, wrk_ref, 2);
		helperfunc::choose_ref(YAll, Res_YSeed, wrk_ref, 2);
		helperfunc::choose_ref(Y1All, Res_Y1Seed, wrk_ref, 2);
		helperfunc::choose_ref(Y2All, Res_Y2Seed, wrk_ref, 2);
		helperfunc::choose_ref(XsaAll, Res_XsaSeed, wrk_ref, 2);
		helperfunc::choose_ref(ylimAll, Res_ylimSeed, wrk_ref, 2);
	}
	else {
		vector<bool>wrk_ref;
		for (int i = 0; i < indSimFailAll.size(); i++)
		{
			bool tmp = (indSimFailAll[i] > 0 ? true : false);
			wrk_ref.push_back(tmp);
		}
		helperfunc::delete_ref(XAll, wrk_ref, 2);
		helperfunc::delete_ref(YAll, wrk_ref, 2);
		helperfunc::delete_ref(Y1All, wrk_ref, 2);
		helperfunc::delete_ref(Y2All, wrk_ref, 2);
		helperfunc::delete_ref(XsaAll, wrk_ref, 2);
		helperfunc::delete_ref(ylimAll, wrk_ref, 2);
		double probSimFail = 0;
		for (int i = 0; i < indSimFailAll.size(); i++) {
			probSimFail += indSimFailAll[i];
		}
		probSimFail /= (double(indSimFailAll.size()));
		double   probTarget = (PROBSUB - probSimFail) / (1 - probSimFail);
		getspec(Y1All, Y2All, probTarget, perfDeltaCur);

		indFailAll.clear();
		for (int i = 0; i < nSeed; i++) {
			nSim[i] = Y[i].size();
			for (int j = 0; j < Y1[i].size(); j++) {
				indParamFail[i][0][j] = Y1[i][j][0] < perfDeltaCur[0] || Y2[i][j][0] > perfDeltaCur[1];
				indSimFail[i][0][j] = isnan(Y[i][j][0]);
				indFail[i][0][j] = indParamFail[i][0][j] || indSimFail[i][0][j];
			}
			for (int j = 0; j < indFail[i][0].size(); j++) {
				indFailAll.push_back(indFail[i][0][j]);
			}
		}
		probCur = helperfunc::mean(indFailAll);
		vector<bool>wrk_ref2;
		for (int i = 0; i < indSimFailAll.size(); i++) {
			bool tmp = (indSimFailAll[i] > 0 ? true : false);
			wrk_ref2.push_back(tmp);
		}
		helperfunc::delete_ref(indFailAll, wrk_ref2);
		vector<bool>wrk_ref3;
		for (int i = 0; i < indFailAll.size(); i++) {
			bool tmp = (indFailAll[i] > 0 ? true : false);
			wrk_ref3.push_back(tmp);
		}
		helperfunc::choose_ref(XAll, Res_XSeed, wrk_ref3, 2);
		helperfunc::choose_ref(YAll, Res_YSeed, wrk_ref3, 2);
		helperfunc::choose_ref(Y1All, Res_Y1Seed, wrk_ref3, 2);
		helperfunc::choose_ref(Y2All, Res_Y2Seed, wrk_ref3, 2);
		helperfunc::choose_ref(XsaAll, Res_XsaSeed, wrk_ref3, 2);
		helperfunc::choose_ref(ylimAll, Res_ylimSeed, wrk_ref3, 2);
	}

	vector<double>Z(nSeed);
	double sigCur;
	if (probCur == 1) {
		for (int i = 0; i < nSeed; i++) {
			Z[i] = helperfunc::mean(indFail[i][0]);
		}
		sigCur = 0;
	}
	else {
		vector<double>Z(nSeed);
		for (int i = 0; i < nSeed; i++) {
			Z[i] = helperfunc::mean(indFail[i][0]);
		}
		sigCur = sqrt(helperfunc::var(Z) / nSeed);
	}
	probList.push_back(probCur);
	sigList.push_back(sigCur);
	perfDeltaList.push_back(perfDeltaCur);
	if (Res_eop == 1) {
		cout << "\tFinish " << Res_YSeed.size() << "th order APA......current order simulation times : " << nSimTotal << endl;
		output(probList, sigList, nSimTotal);
		if (Res_YSeed.size() < nCend) {
			sus_delta_sim(Res_XSeed, Res_YSeed, Res_XsaSeed, Res_ylimSeed, simtype);
		}
	}
	return 1;
};


bool SUS::expandSeed(vector<vector<double> >&XSeed, vector<vector<double> >&YSeed, vector<vector<double> >&XsaSeed, vector<vector<double> >&ylimSeed, int & nSimTotal, int nSimiter, vector<vector<double> >&Res_XSeed, vector<vector<double> >&Res_YSeed, vector<vector<double> >&Res_Xsa, vector<vector<double> >&Res_ylim, vector<vector<vector<double> > >&Xcell, vector<vector<vector<double> > >&Ycell, int simtype) {
	int nVar = XSeed.size(), nSeed = XSeed[0].size();
	vector<vector<vector<double> > >X(nSeed), Y(nSeed), Xsa(nSeed), ylim(nSeed);
	int nSim = nSeed, nAccept = 0;
	for (int i = 0; i < nSeed; i++) {
		// store the ith chain in X[i]
		vector<double>tmp1(XSeed.size());
		for (int j = 0; j < tmp1.size(); j++) {
			tmp1[j] = XSeed[j][i];
		}
		X[i].push_back(tmp1);

		vector<double>tmp2(YSeed.size());
		for (int j = 0; j < tmp2.size(); j++) {
			tmp2[j] = YSeed[j][i];
		}
		Y[i].push_back(tmp2);
		vector<double>tmp3(XsaSeed.size());
		for (int j = 0; j < tmp3.size(); j++) {
			tmp3[j] = XsaSeed[j][i];
		}
		Xsa[i].push_back(tmp3);
		vector<double>tmp4(ylimSeed.size());
		for (int j = 0; j < tmp4.size(); j++) {
			tmp4[j] = ylimSeed[j][0];
		}
		ylim[i].push_back(tmp4);
	}
	while (nSim < nSimiter) {
		for (int i = 0; i < nSeed; i++) {
			vector<double> wrk(X[i][X[i].size() - 1]);
			vector<vector<double> >XTemp((X[i][X[i].size() - 1]).size(), vector<double>(1));
			metropolis(wrk, nVar, XTemp);
			vector<vector<double> >YTemp(XTemp.size() / (2 * NumTrans), vector<double>(1));
			simout(XTemp, YTemp, nSimTotal, 1, simtype);

			vector<double> wrk2(Xsa[i][Xsa[i].size() - 1]);
			vector<vector<double> >XsaTemp((Xsa[i][Xsa[i].size() - 1]).size(), vector<double>(1));
			metropolis(wrk2, 2, XsaTemp);
			vector<double>ylimit_V(1);
			sim_SA_fake(XsaTemp, ylimit_V, simtype);

			vector<double> ylimTemp(2);
			ylimTemp[0] = ylimit_V[0] / t_C + u01[0];
			ylimTemp[1] = ylimit_V[0] / t_C + u01[1];

			vector<vector<double> >Y1Temp(YTemp), Y2Temp(YTemp);
			double max1 = MYMIN, min2 = MYMAX;
			int hasnan = 0;
			for (int i = 0; i < YTemp.size(); i++) {
				hasnan += isnan(Y1Temp[i][0]);
				Y1Temp[i][0] -= ylimTemp[0];
				Y2Temp[i][0] -= ylimTemp[1];
				max1 = max1 > Y1Temp[i][0] ? max1 : Y1Temp[i][0];
				min2 = min2 < Y2Temp[i][0] ? min2 : Y2Temp[i][0];
			}
			if (hasnan > 0 || max1 < 0 || min2>0) {
				vector<double>wrk_XTemp(XTemp.size());
				transpose(XTemp, wrk_XTemp);
				X[i].push_back(wrk_XTemp);

				vector<double>wrk_YTemp(YTemp.size());
				transpose(YTemp, wrk_YTemp);
				Y[i].push_back(wrk_YTemp);

				vector<double>wrk_XsaTemp(XsaTemp.size());
				transpose(XsaTemp, wrk_XsaTemp);
				Xsa[i].push_back(wrk_XsaTemp);

				ylim[i].push_back(ylimTemp);
				nAccept += 1;
			}
			else {
				X[i].push_back(X[i][X[i].size() - 1]);
				Y[i].push_back(Y[i][Y[i].size() - 1]);
				Xsa[i].push_back(Xsa[i][Xsa[i].size() - 1]);
				ylim[i].push_back(ylim[i][ylim[i].size() - 1]);
			}
			nSim += 1;
			if (nSim >= nSimiter) {
				break;
			}
		}
	}
	nSeed = X.size();
	vector<vector<double> > tmp_XSeed, tmp_Yseed, tmp_Xsa, tmp_ylim;
	for (int i = 0; i < nSeed; i++) {
		for (int j = 0; j < X[i].size(); j++) {
			tmp_XSeed.push_back(X[i][j]);
		}
		for (int j = 0; j < Y[i].size(); j++) {
			tmp_Yseed.push_back(Y[i][j]);
		}
		for (int j = 0; j < Xsa[i].size(); j++) {
			tmp_Xsa.push_back(Xsa[i][j]);
		}
		for (int j = 0; j < ylim[i].size(); j++) {
			tmp_ylim.push_back(ylim[i][j]);
		}

	}
	Res_XSeed = transpose(tmp_XSeed);
	Res_YSeed = transpose(tmp_Yseed);
	Res_Xsa = transpose(tmp_Xsa);
	Res_ylim = transpose(tmp_ylim);
	Xcell = X; Ycell = Y;
	return 1;
}

bool SUS::transpose(vector<double>&src, vector<vector<double> > &dst) {
	if (dst[0].size() != 1 || dst.size() != src.size()) {
		cout << "enter function transpose, dimension is wrong,error!\n";
		return 0;
	}
	for (int i = 0; i < src.size(); i++) {
		dst[i][0] = src[i];
	}
	return 1;
}

vector<vector<double> > SUS::transpose(vector<vector<double> >&src) {
	vector<vector<double> > dst(src[0].size(), vector<double>(src.size()));
	for (int i = 0; i < dst.size(); i++) {
		for (int j = 0; j < dst[0].size(); j++) {
			dst[i][j] = src[j][i];
		}
	}
	return dst;
}

bool SUS::transpose(vector<vector<double> >&src, vector<double> &dst) {
	if (src[0].size() != 1 || src.size() != dst.size()) {
		cout << "enter function transpose, dimension is wrong,error!\n";
		return 0;
	}
	for (int i = 0; i < src.size(); i++) {
		dst[i] = src[i][0];
	}
	return 1;
}

bool SUS::metropolis(vector<double>&X, int nVar, vector<vector<double> > &XNext) {
	for (int i = 0; i < X.size(); i++) {
		XNext[i][0] = X[i];
	}
	double step = 1;
	//std::uniform_real_distribution<double>randengine(0, 1);
	//std::normal_distribution<double>n(0, 1);
	vector<double>Xtemp(nVar);
	for (int i = 0; i < nVar; i++) {
		Xtemp[i] = step*n(e) + X[i];
	}
	//vector<double>r1(nVar),r2(nVar);
	double r1, r2;
	for (int i = 0; i < nVar; i++) {
		r1 = normpdf(Xtemp[i]) / normpdf(X[i]);
		r2 = randengine(e);
		if (r2 <= r1) {
			XNext[i][0] = Xtemp[i];
		}
	}
	return 1;
}

double SUS::normpdf(double src, double mu, double sigma) {
	return exp(-0.5*(src - mu) / sigma*(src - mu) / sigma) / (sqrt(2 * PI)*sigma);
}

bool SUS::getspec(vector<double> Y1, vector<double> Y2, double probTarget, vector<double> &perfDeltaSub) {
	sort(Y1.begin(), Y1.end());
	sort(Y2.begin(), Y2.end());
	int nY = Y1.size();
	if (nY != Y2.size()) {
		cout << "enter function getspec,dimension is wrong,error!\n";
		return 0;
	}
	//double nFailTarget = round(probTarget*nY);
	double nFailTarget = ceil(probTarget*nY);
	if (u01[0] == -INFINITY) {
		perfDeltaSub[0] = 0;
		perfDeltaSub[1] = (Y2[nY - 1 - nFailTarget] + Y2[nY - 1 - nFailTarget + 1]) / 2;
		//perfDeltaSub[1] = (1+MYEPSILON)*(Y2[nY - 1 - nFailTarget] + Y2[nY - 1 - nFailTarget + 1]) / 2;
	}
	else {
		if (u01[1] == INFINITY) {
			perfDeltaSub[0] = (Y1[nFailTarget - 1] + Y1[nFailTarget]) / 2;
			//perfDeltaSub[0] = (1+MYEPSILON)*(Y1[nFailTarget - 1] + Y1[nFailTarget]) / 2;
			perfDeltaSub[1] = 0;
		}
		else {
			int nFailTarget1 = round(nFailTarget / 2.0);
			int nFailTarget2 = nFailTarget - nFailTarget1;
			perfDeltaSub[1] = (Y2[nY - 1 - nFailTarget2] + Y2[nY - nFailTarget2]) / 2;
			perfDeltaSub[0] = (Y1[nFailTarget1 - 1] + Y1[nFailTarget1]) / 2;
			//perfDeltaSub[1] = (1+MYEPSILON)*(Y2[nY - 1 - nFailTarget2] + Y2[nY - nFailTarget2]) / 2;
			//perfDeltaSub[0] = (1+MYEPSILON)*(Y1[nFailTarget1 - 1] + Y1[nFailTarget1]) / 2;
		}
	}
	return 1;
}


bool SUS::output(vector<double>probList, vector<double>sigList, int &nSimTotal) {
	double probEst = 1;
	//string tmp; tmp = "problist"; display_matrix_vector(probList, tmp);
	int length = probList.size();
	if (length != sigList.size()) {
		cout << "enter function output,dimension is wrong,error!\n";
		return 0;
	}
	string tmp = "subsets probability"; helperfunc::display_matrix_vector(probList, tmp);
	for (int i = 0; i < length; i++) {
		probEst *= probList[i];
	}

	double muLogP = log10(probEst);
	vector<double>sigLogPList(length);
	double sigLogP = 0;
	for (int i = 0; i<length; i++) {
		sigLogPList[i] = sigList[i] / probList[i] / log(10);
		sigLogP += sigLogPList[i] * sigLogPList[i];
		if (i > 0) {
			sigLogP += 2 * sigLogPList[i] * sigLogPList[i - 1];
		}
	}
	sigLogP = sqrt(sigLogP);
	double probCILow = pow(10, muLogP - 1.96*sigLogP);
	double probCIUp = pow(10, muLogP + 1.96*sigLogP);

	//store values
	APA_probCILowList.push_back(probCILow);
	APA_probCIUpList.push_back(probCIUp);
	APA_probEstList.push_back(probEst);
	APA_simtotalList.push_back(nSimTotal);
	return 1;
}

bool SUS::sim_SA_fake(vector<vector<double> >&Xsa, vector<double>&Res_ylimit, int simtype) {
	switch (simtype)
	{
	case SAFAIL_: {
		if (Xsa.size() != 2 || Xsa[0].size() != Res_ylimit.size()) {
			cout << "Enter function sim_SA_fake, dimensiont is wronr,error!\n";
			return 0;
		}
		double	rate = 0.1;
		double vthn0 = 0.175;
		for (int i = 0; i < Res_ylimit.size(); i++) {
			double vthn1 = vthn0 + rate*vthn0*Xsa[0][i];
			double vthn2 = vthn0 + rate*vthn0*Xsa[1][i];
			Res_ylimit[i] = vthn1 - vthn2;
		}
		break;
	}
	case READFAIL_: {
		if (Xsa.size() != 2 || Xsa[0].size() != Res_ylimit.size()) {
			cout << "Enter function sim_SA_fake, dimensiont is wronr,error!\n";
			return 0;
		}
		for (int i = 0; i < Res_ylimit.size(); i++) {
			//Res_ylimit[i] = 0;
			Res_ylimit[i] = n(e)*REALCOR;
			//Res_ylimit[i] = n(e) * 0;
			//cout << Res_ylimit[i] << endl;
		}
		break;
	}
	case WRITEFAIL_: {
		if (Xsa.size() != 2 || Xsa[0].size() != Res_ylimit.size()) {
			cout << "Enter function sim_SA_fake, dimensiont is wronr,error!\n";
			return 0;
		}
		for (int i = 0; i < Res_ylimit.size(); i++) {
			//Res_ylimit[i] = 0;
			Res_ylimit[i] = REALCOR2*n(e)*TWL;
		}
		break;
	}
	default: {cout << "in function sim_AS_fake, wrong simulation type, error!\n"; break; };
	}
	return 1;

}
/*
if (simtype == SAFAIL_) {
if (Xsa.size() != 2 || Xsa[0].size() != Res_ylimit.size()) {
cout << "Enter function sim_SA_fake, dimensiont is wronr,error!\n";
return 0;
}
double	rate = 0.1;
double vthn0 = 0.175;
for (int i = 0; i < Res_ylimit.size(); i++) {
double vthn1 = vthn0 + rate*vthn0*Xsa[0][i];
double vthn2 = vthn0 + rate*vthn0*Xsa[1][i];
Res_ylimit[i] = vthn1 - vthn2;
}
return 1;
}
else {
if (simtype == READFAIL_) {
if (Xsa.size() != 2 || Xsa[0].size() != Res_ylimit.size()) {
cout << "Enter function sim_SA_fake, dimensiont is wronr,error!\n";
return 0;
}
for (int i = 0; i < Res_ylimit.size(); i++) {
Res_ylimit[i] = 0;
}
}
}
}*/

bool SUS::simout(vector<vector<double> >&X, vector<vector<double> >&Y, int &nSimTotal, bool epo, int simtype) {
	//  when epo==1, it will return a matrix in Y;
	//cout<<"In simout "<<X.size()<<"\t"<<X[0].size()<<endl;
	//system("read");
	if (epo == 1) {
		if (Y[0].size() != X[0].size() || Y.size() * 2 * NumTrans != X.size()) {
			cout << "enter function simout, dimension is wrong for epo==1,error!\n";
			return 0;
		}
		int ncells = X.size(), nSim = X[0].size();
		vector<double> bias(nSim);
		gen_corr(bias, simtype);
		for (int i = 0; i < ncells / (2 * NumTrans); i++)
		{
			vector<vector<double> > tmp(2 * NumTrans, vector<double>(nSim));
			for (int m = 0; m < 2 * NumTrans; m++) {
				for (int n = 0; n < nSim; n++) {
					tmp[m][n] = X[2 * NumTrans*i + m][n] + bias[n];
				}
			}
			simX(tmp, Y, i, simtype);
		}
		nSimTotal += nSim*ncells / 6;
		return 1;
	}
	else {
		if (Y[0].size() != X[0].size() || Y.size() != 1) {
			cout << "enter function simout, dimension is wrong for epo==0,error!\n";
			return 0;
		}
		int ncells = X.size(), nSim = X[0].size();
		vector<double> bias(nSim);
		gen_corr(bias, simtype);
		vector<vector<double> >Yall(ncells / (2 * NumTrans), vector<double>(nSim));
		for (int i = 0; i < ncells / (2 * NumTrans); i++)
		{
			vector<vector<double> > tmp(2 * NumTrans, vector<double>(nSim));
			for (int m = 0; m < 2 * NumTrans; m++) {
				for (int n = 0; n < nSim; n++) {
					tmp[m][n] = X[2 * NumTrans*i + m][n] + bias[n];
				}
			}
			simX(tmp, Yall, i, simtype);
		}
		for (int i = 0; i < nSim; i++) {
			double tmp = MYMIN;
			for (int j = 0; j < ncells / (2 * NumTrans); j++) {
				tmp = tmp > Yall[j][i] ? tmp : Yall[j][i];
			}
			Y[0][i] = tmp;
		}
		nSimTotal += nSim*ncells / 6;
		return 1;
	}
}

bool SUS::simX(vector<vector<double> > &src, vector<vector<double> >&dst, int index, int simtype) {
	//cout<<"IN simX : "<<src.size()<<"\t"<<src[0].size()<<endl;
	//system("read");
	switch (simtype)
	{
	case SAFAIL_: {
		if (index < 0 || index >= dst.size()) {
			cout << "Entering function simX. dimension is wrong,error!\n";
			return 0;
		}
		simulator Simulator;
		int length = dst[0].size();
		vector<double>iPG(length);
		vector<double>V2(length);
		vector<vector<double> >wrk(NumTrans, vector<double>(length));
		for (int i = 0; i < NumTrans; i++) {
			for (int j = 0; j < length; j++) {
				wrk[i][j] = src[i][j];
			}
		}
		Simulator.dcsim(wrk, iPG, V2);
		for (int i = 0; i < length; i++) {
			dst[index][i] = iPG[i];
		}
		break;
	}
	case READFAIL_: {
		if (index < 0 || index >= dst.size()) {
			cout << "Entering function simX. dimension is wrong,error!\n";
			return 0;
		}
		int length = dst[0].size();
		simulator Simulator;
		vector<double>delta_V(length);
		Simulator.readsim(src, delta_V);
		//cout<<"in simx "<<src.size()<<"\t"<<src[0].size()<<endl;

		for (int i = 0; i < length; i++) {
			dst[index][i] = delta_V[i];
		}
		break;
	}
	case WRITEFAIL_: {
		if (index < 0 || index >= dst.size()) {
			cout << "Entering function simX. dimension is wrong,error!\n";
			return 0;
		}
		int length = dst[0].size();
		simulator Simulator;
		vector<double>delta_t(length);
		Simulator.writesim(src, delta_t);
		for (int i = 0; i < length; i++) {
			dst[index][i] = delta_t[i];
		}
		break;
	}
	default: {cout << "in function simX, Wrong SImulation type!\n"; }
	}
	return 1;
}

bool SUS::simout2(vector<vector<double> >&X, vector<vector<double> >&Y, int &nSimTotal, bool epo, int simtype) {
	//  when epo==1, it will return a matrix in Y;
	//cout<<"In simout "<<X.size()<<"\t"<<X[0].size()<<endl;
	//system("read");
	if (epo == 1) {
		if (Y[0].size() != X[0].size() || Y.size() * 2 * NumTrans != X.size()) {
			cout << "enter function simout, dimension is wrong for epo==1,error!\n";
			return 0;
		}
		int ncells = X.size(), nSim = X[0].size();
		for (int i = 0; i < ncells / (2 * NumTrans); i++)
		{
			vector<vector<double> > tmp(2 * NumTrans, vector<double>(nSim));
			for (int m = 0; m < 2 * NumTrans; m++) {
				for (int n = 0; n < nSim; n++) {
					tmp[m][n] = X[2 * NumTrans*i + m][n];
				}
			}
			simX(tmp, Y, i, simtype);
		}
		nSimTotal += nSim*ncells / 6;
		return 1;
	}
	else {
		if (Y[0].size() != X[0].size() || Y.size() != 1) {
			cout << "enter function simout, dimension is wrong for epo==0,error!\n";
			return 0;
		}
		int ncells = X.size(), nSim = X[0].size();
		vector<vector<double> >Yall(ncells / (2 * NumTrans), vector<double>(nSim));
		for (int i = 0; i < ncells / (2 * NumTrans); i++)
		{
			vector<vector<double> > tmp(2 * NumTrans, vector<double>(nSim));
			for (int m = 0; m < 2 * NumTrans; m++) {
				for (int n = 0; n < nSim; n++) {
					tmp[m][n] = X[2 * NumTrans*i + m][n];
				}
			}
			simX(tmp, Yall, i, simtype);
		}
		for (int i = 0; i < nSim; i++) {
			double tmp = MYMAX;
			for (int j = 0; j < ncells / (2 * NumTrans); j++) {
				tmp = tmp < Yall[j][i] ? tmp : Yall[j][i];
			}
			Y[0][i] = tmp;
		}
		nSimTotal += nSim*ncells / 6;
		return 1;
	}
}



bool SUS::MC(vector<int> &src, vector<double> &dst, int simtype, int nsim) {
	if (src.size() != dst.size()) {
		cout << "[Warning]:In function APA::MC,dimension is wrong!\n";
	}
	for (int j = 0; j<src.size(); j++) {
		int cur_cell = src[j];
		int failnum = 0;
		for (int i = 0; i<nsim; i++) {
			vector<vector<double> > cur_x(cur_cell * 2 * NumTrans, vector<double>(1));
			for (int p = 0; p<cur_x.size(); p++) {
				cur_x[p][0] = n(e);
			}
			vector<vector<double> > wrk(cur_cell, vector<double>(1));
			int simtimes;
			simout2(cur_x, wrk, simtimes, 1, simtype);
			vector<vector<double> > ylim(2, vector<double>(1));
			vector <vector<double> > Xsa(2, vector<double>(1));
			vector<double>ylimit_V(1);
			Xsa[0][0] = n(e);
			Xsa[1][0] = n(e);
			sim_SA_fake(Xsa, ylimit_V, simtype);
			ylim[0][0] = ylimit_V[0] / t_C + u01[0];
			ylim[1][0] = sig01[1] * n(e) + u01[1];
			for (int p = 0; p<wrk.size(); p++) {
				if (isnan(wrk[p][0]) || wrk[p][0]<ylim[0][0] || wrk[p][0]>ylim[1][0]) {
					failnum += 1;
					break;
				}
			}
			cout << "Finish " << i << "/" << nsim << " of the " << j << "th cell" << endl;
		}
		dst[j] = double(failnum) / nsim;
	}
	return true;
}



bool SUS::gen_corr(vector<double> &bias, int simtype) {
	for (int i = 0; i < bias.size(); i++)
	{
		if (simtype == READFAIL_) { bias[i] = n(e)*FAKECOR; }
		if (simtype == SAFAIL_) { bias[i] = 0; }//bias[i]=0;
		if (simtype == WRITEFAIL_) { bias[i] = n(e)*FAKECOR2; }

	}
	return true;
}