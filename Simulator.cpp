#include"Simulator.h"

simulator::~simulator() {};

simulator::simulator() {
	M_ND_[0] = 2; M_ND_[1] = 1; M_ND_[2] = 2;
	M_NG_[0] = M_NG_[1] = M_NG_[2] = 1;
	M_NS_[0] = -1; M_NS_[1] = 2; M_NS_[2] = 1;
	M_TYPE_[0] = M_TYPE_[1] = M_TYPE_[3] = M_TYPE_[4] = NMOS_; M_TYPE_[2] = M_TYPE_[5] = PMOS_;
	M_VT_[0] = M_VT_[1] = M_VT_[3] = M_VT_[4] = 0.1750; M_VT_[2] = M_VT_[5] = -0.2640;
	M_W_[0] = M_W_[3] = 0.16E-6; M_W_[1] = M_W_[4] = 0.07E-6; M_W_[2] = M_W_[5] = 0.08E-6;
	M_LAMBDA_[0] = M_LAMBDA_[1] = M_LAMBDA_[2] = M_LAMBDA_[3] = M_LAMBDA_[4] = M_LAMBDA_[5] = 0;
	M_L_[0] = M_L_[1] = M_L_[2] = M_L_[3] = M_L_[4] = M_L_[5] = 0.6E-7;
	M_MU_[0] = M_MU_[1] = M_MU_[3] = M_MU_[4] = 0.0204; M_MU_[2] = M_MU_[5] = 0.0090;
	M_COX_[0] = M_COX_[1] = M_COX_[3] = M_COX_[4] = 0.013595; M_COX_[2] = M_COX_[5] = 0.01298170;
	M_COV_[0] = M_COV_[1] = M_COV_[3] = M_COV_[4] = 1.79454E-4; M_COV_[2] = M_COV_[5] = 1.6097308E-4;
	/*M_ND_[NumTrans] = { 2,1,2 }, M_NG_[NumTrans] = { 1,1,1 }, M_NS_[NumTrans] = { -1,2,1 };
	M_TYPE_[2 * NumTrans] = { NMOS_,NMOS_,PMOS_,NMOS_,NMOS_,PMOS_ };
	M_VT_[2 * NumTrans] = { 0.1750,0.1750,-0.2640, 0.1750,0.1750,-0.2640 };
	M_W_[2 * NumTrans] = { 0.16E-6,0.07E-6,0.08E-6,0.16E-6,0.07E-6,0.08E-6 };
	M_L_[2 * NumTrans] = { 0.6E-7 ,0.6E-7,0.6E-7,0.6E-7 ,0.6E-7,0.6E-7 };
	M_MU_[2 * NumTrans] = { 0.0204,0.0204,.0090 ,0.0204,0.0204,.0090 };
	M_COX_[2 * NumTrans] = { 0.013595,0.013595,0.01298170,0.013595,0.013595,0.01298170 };
	M_LAMBDA_[2 * NumTrans] = { 0,0,0,0,0,0 };
	M_COV_[2 * NumTrans] = { 1.79454E-4,1.79454E-4, 1.6097308E-4,1.79454E-4,1.79454E-4,1.6097308E-4 };*/
};


bool simulator::set_M_VT_all(vector<double>src) {
	if (src.size() != 2 * NumTrans) {
		cout << "Cannot set value for VT beacuse of dimension difference!\n";
		return 0;
	}
	else {
		for (int i = 0; i < 2 * NumTrans; i++) {
			M_VT_[i] = src[i];
		}
	}
	return 1;
};


bool simulator::set_M_VT_(vector<double>src) {
	if (src.size() != NumTrans) {
		cout << "Cannot set value for VT beacuse of dimension difference!\n";
		return 0;
	}
	else {
		for (int i = 0; i < NumTrans; i++) {
			M_VT_[i] = src[i];
		}
	}
	return 1;
};

bool simulator::dcsim(vector<vector<double> > varValue, vector<double> &Res_iPG, vector<double>&Res_V2, bool dis) {
	if (Res_iPG.size() != Res_V2.size() || varValue[0].size() != Res_V2.size() || Res_iPG.size() != varValue[0].size()) {
		cout << "entering function dcsim,dimension is wrong,error!\n";
		return 0;
	}
	simulator::reset_M_VT_all();
	int numSample = varValue[0].size();
	vector<double> paramMean(NumTrans); for (int i = 0; i < NumTrans; i++) { paramMean[i] = M_VT_[i]; }
	int num_NR_max = NewTonIer;
	double v_err_tol = VerrTol;
	double i_err_tol = IerrTol;
	int num_node = NumNodes;
	int num_nlelem = NumTrans;
	int num_vsource = 1;
	//set up for Newton Method
	vector<vector<double> > M(NumTrans, vector<double>(NumTrans)); M[0][NumTrans - 1] = 1; M[NumTrans - 1][0] = 1;
	vector<double>J(NumTrans); J[NumTrans - 1] = 1.2;
	vector<double> PerfSpec(num_node + num_vsource);
	for (int i = 0; i < num_node + num_vsource; i++)
	{
		PerfSpec[i] = i < num_node ? v_err_tol : i_err_tol;
	}
	vector<vector<double> >paramValue(NumTrans, vector<double>(numSample));
	alter(paramMean, varValue, paramValue, SAFAIL_);

	for (int sampleInd = 0; sampleInd < numSample; sampleInd++) {
		vector<double>paramValueCur(NumTrans); for (int i = 0; i < NumTrans; i++) { paramValueCur[i] = paramValue[i][sampleInd]; };
		//tmp_pri = "paramValueCur"; display_matrix_vector(paramValueCur, tmp_pri);
		vector <double>temp = paramValueCur;
		temp[2] = -temp[2];
		int equalzero = 0;
		for (int i = 0; i < NumTrans; i++) { equalzero = temp[i] > 0 ? (equalzero + 1) : equalzero; };
		int indSimFail = equalzero == 3 ? 0 : 1;
		//cout << "indSimFail :" << indSimFail << endl;
		if (indSimFail == 1) {
			Res_iPG[sampleInd] = NAN;
			Res_V2[sampleInd] = NAN;
		}
		else {
			simulator::set_M_VT_(paramValueCur);
			//tmp_pri = "M_VT_"; display_matrix_vector(M_VT_, tmp_pri);
			vector<double> V(NumTrans);
			int num_NR = 0;
			while (true) {
				vector<vector<double> >cur_M = M;
				vector<double>cur_J = J;
				for (int ind = 0; ind < num_nlelem; ind++) {
					stamp_MOSFET(cur_M, cur_J, V, ind);
					//tmp_pri = "cur_M"; display_matrix_vector(cur_M, tmp_pri, 10);
					//tmp_pri = "cur_J"; display_matrix_vector(cur_J, tmp_pri, 10);
				}
				vector<bool> NZ_M(cur_M[0].size()), wrk_NZ_M(cur_M[0].size());
				for (int i = 0; i < cur_M[0].size(); i++) {
					double cur = 0;
					for (int j = 0; j < cur_M.size(); j++) {
						cur += abs(cur_M[j][i]);
					}
					NZ_M[i] = cur > 0 ? 1 : 0;
					wrk_NZ_M[i] = cur > 0 ? 0 : 1;
				};
				helperfunc::delete_ref(cur_M, wrk_NZ_M, 2);
				//tmp_pri = "cur_M"; display_matrix_vector(cur_M, tmp_pri,10);
				//tmp_pri = "cur_J"; display_matrix_vector(cur_J, tmp_pri,10);
				vector<double>NZ_V(cur_M[0].size());
				if (simple_solve(cur_M, cur_J, NZ_V) == 0) {
					cout << "Solve matrix wrong !\n";
					return 0;
				}
				num_NR = num_NR + 1;
				bool hasnan = false;
				for (int i = 0; i < NZ_V.size(); i++) {
					if (isnan(NZ_V[i]) == true) {
						hasnan = true;
						break;
					}
				}
				if (hasnan) {
					cout << "Newton Method fails because of NAN!" << endl;
				}
				vector<double>	cur_V(V.size());
				for (int i = 0, cur = 0; i < cur_V.size(); i++) {
					if (NZ_M[i] == 1) {
						cur_V[i] = NZ_V[cur];
						cur++;
					}
					else {
						cur_V[i] = 0;
					}
				}
				bool allsmall = true;
				for (int i = 0; i < PerfSpec.size(); i++) {
					if (abs(cur_V[i] - V[i]) >= PerfSpec[i]) {
						allsmall = false;
						break;
					}
				}
				if (allsmall == true) {
					V = cur_V;
					break;
				}
				else {
					if (num_NR > num_NR_max) {
						cout << "Increase the # of iterations!\n";
					}
					else
						V = cur_V;
				}
				//display_matrix_vector(V, tmp_pri,10);
			}

			double	Beta = M_MU_[1] * M_COX_[1] * M_W_[1] / M_L_[1];
			double Vgs = V[0] - V[1];
			double Vds = V[0] - V[1];
			double delta_Vgs = Vgs - M_VT_[1];
			double Ids = 0;
			if (delta_Vgs > 0)
			{
				Ids = Beta * delta_Vgs *delta_Vgs * (1 + M_LAMBDA_[1] * Vds) / 2;
			}
			Res_iPG[sampleInd] = Ids;
			Res_V2[sampleInd] = V[1];
		}
	}
	if (dis) {
		string tmp;
		tmp = "Res_iPG"; helperfunc::display_matrix_vector(Res_iPG, tmp, 8);
		tmp = "Res_V2"; helperfunc::display_matrix_vector(Res_V2, tmp, 8);
	}
	return 1;
};

bool simulator::stamp_MOSFET(vector<vector<double> > &M, vector<double>&J, vector<double> &old_v, int index) {
	//int index = MOS_Num - 1;  MOs_num is 1,2,3 ; index is 0,1,2
	int Nd = M_ND_[index], Ng = M_NG_[index], Ns = M_NS_[index];
	double Beta = M_MU_[index] * M_COX_[index] * M_W_[index] / M_L_[index];

	double Vd = Nd > 0 ? old_v[Nd - 1] : 0;
	double Vs = Ns > 0 ? old_v[Ns - 1] : 0;
	double Vg = Ng > 0 ? old_v[Ng - 1] : 0;

	double gm = 0, Vgs = 0, gds = 0, Vds = 0, Ids = 0, Vgd = 0, delta_Vgs = 0, delta_Vgd = 0;

	if (M_TYPE_[index] == PMOS_) {
		if (Vd > Vs) {
			int tempN = Nd;
			double tempV = Vd;
			Nd = Ns;
			Vd = Vs;
			Ns = tempN;
			Vs = tempV;
		}
		Vgs = Vg - Vs;
		Vds = Vd - Vs;
		Vgd = Vg - Vd;
		delta_Vgs = Vgs - M_VT_[index];
		delta_Vgd = Vgd - M_VT_[index];
		if (delta_Vgs >= 0) {
			Ids = 0;
			gm = 0;
			gds = 0;
		}
		else
		{
			if (delta_Vgd <= 0) {
				Ids = -Beta * (2 * delta_Vgs*Vds - Vds *Vds) * (1 - M_LAMBDA_[index] * Vds) / 2;
				gm = -Beta * Vds * (1 - M_LAMBDA_[index] * Vds);
				gds = -Beta * (delta_Vgs * (1 - 2 * M_LAMBDA_[index] * Vds) - Vds * (1 - 1.5*M_LAMBDA_[index] * Vds));
			}
			else {
				if (delta_Vgd > 0) {
					Ids = -Beta * delta_Vgs *delta_Vgs * (1 - M_LAMBDA_[index] * Vds) / 2;
					gm = 2 * Ids / delta_Vgs;
					gds = -Ids * M_LAMBDA_[index] / (1 - M_LAMBDA_[index] * Vds);
				}
			}
		}
	}
	else {
		if (M_TYPE_[index] == NMOS_) {
			if (Vd < Vs) {
				int tempN = Nd;
				double tempV = Vd;
				Nd = Ns;
				Vd = Vs;
				Ns = tempN;
				Vs = tempV;
			}
			Vgs = Vg - Vs;
			Vds = Vd - Vs;
			Vgd = Vg - Vd;
			delta_Vgs = Vgs - M_VT_[index];
			delta_Vgd = Vgd - M_VT_[index];
			if (delta_Vgs <= 0) {
				Ids = 0;
				gm = 0;
				gds = 0;
			}
			else {
				if (delta_Vgd >= 0) {
					Ids = Beta * (2 * delta_Vgs*Vds - Vds*Vds) * (1 + M_LAMBDA_[index] * Vds) / 2;
					gm = Beta * Vds * (1 + M_LAMBDA_[index] * Vds);
					gds = Beta * (delta_Vgs * (1 + 2 * M_LAMBDA_[index] * Vds) - Vds * (1 + 1.5*M_LAMBDA_[index] * Vds));
				}
				else {
					if (delta_Vgd < 0) {
						Ids = Beta * delta_Vgs *delta_Vgs * (1 + M_LAMBDA_[index] * Vds) / 2;
						gm = 2 * Ids / delta_Vgs;
						gds = Ids * M_LAMBDA_[index] / (1 + M_LAMBDA_[index] * Vds);
					}
				}
			}
		}
		else { cout << "Simulating " << index << "th(index) MOSFET, Not NOMS and not PMOS,type is wrong!" << endl; }
	}

	double	IdsIeq = Ids - gm * Vgs - gds * Vds;
	if (Nd > 0 && Ns > 0 && Ng > 0) {
		J[Nd - 1] = J[Nd - 1] - IdsIeq;
		J[Ns - 1] = J[Ns - 1] + IdsIeq;
		M[Nd - 1][Nd - 1] = M[Nd - 1][Nd - 1] + gds;
		M[Nd - 1][Ns - 1] = M[Nd - 1][Ns - 1] - gds - gm;
		M[Ns - 1][Nd - 1] = M[Ns - 1][Nd - 1] - gds;
		M[Ns - 1][Ns - 1] = M[Ns - 1][Ns - 1] + gds + gm;
		M[Nd - 1][Ng - 1] = M[Nd - 1][Ng - 1] + gm;
		M[Ns - 1][Ng - 1] = M[Ns - 1][Ng - 1] - gm;
	}
	else {
		if (Nd > 0 && Ns > 0 && Ng < 0) {
			J[Nd - 1] = J[Nd - 1] - IdsIeq;
			J[Ns - 1] = J[Ns - 1] + IdsIeq;
			M[Nd - 1][Nd - 1] = M[Nd - 1][Nd - 1] + gds;
			M[Nd - 1][Ns - 1] = M[Nd - 1][Ns - 1] - gds - gm;
			M[Ns - 1][Nd - 1] = M[Ns - 1][Nd - 1] - gds;
			M[Ns - 1][Ns - 1] = M[Ns - 1][Ns - 1] + gds + gm;
		}
		else {
			if (Nd > 0 && Ns < 0 && Ng > 0) {
				J[Nd - 1] = J[Nd - 1] - IdsIeq;
				M[Nd - 1][Nd - 1] = M[Nd - 1][Nd - 1] + gds;
				M[Nd - 1][Ng - 1] = M[Nd - 1][Ng - 1] + gm;
			}
			else {
				if (Nd > 0 && Ns < 0 && Ng < 0) {
					J[Nd - 1] = J[Nd - 1] - IdsIeq;
					M[Nd - 1][Nd - 1] = M[Nd - 1][Nd - 1] + gds;
				}
				else {
					if (Nd < 0 && Ns > 0 && Ng > 0) {
						J[Ns - 1] = J[Ns - 1] + IdsIeq;
						M[Ns - 1][Ns - 1] = M[Ns - 1][Ns - 1] + gds + gm;
						M[Ns - 1][Ng - 1] = M[Ns - 1][Ng - 1] - gm;
					}
					else {
						if (Nd < 0 && Ns > 0 && Ng < 0) {
							J[Ns - 1] = J[Ns - 1] + IdsIeq;
							M[Ns - 1][Ns - 1] = M[Ns - 1][Ns - 1] + gds + gm;
						}
						else {
							cout << "Drain and Source are shorted!/n";
						}
					}
				}
			}
		}
	}
	return 1;
};

bool simulator::alter(vector<double> &paramMean, vector<vector<double> >&varValue, vector<vector<double> >&Res_paramValue, int simtype) {
	int numSample = varValue[0].size();
	//cout<<"Entering alter\n";
	vector<double> cor_v(numSample);
	simulator::addcorelation(numSample, cor_v);
	switch (simtype)
	{
	case SAFAIL_: {
		vector<double>delta(3); delta[0] = 0.111*0.175 / sqrt(2); delta[1] = 0.111*0.175; delta[2] = 0.111*0.264;
		for (int i = 0; i < varValue.size(); i++)
		{
			for (int j = 0; j < numSample; j++)
			{
				Res_paramValue[i][j] = paramMean[i] + varValue[i][j] * delta[i] + cor_v[j];
			}
		}
		break;
	}
	case READFAIL_: {
		vector<double>delta{ 0.111*0.175 / sqrt(2),  0.111*0.175,  0.111*0.264,0.111*0.175 / sqrt(2),  0.111*0.175, 0.111*0.264 };
		for (int i = 0; i < varValue.size(); i++)
		{
			for (int j = 0; j < numSample; j++)
			{
				Res_paramValue[i][j] = paramMean[i] + varValue[i][j] * delta[i] + cor_v[j];
			}
		}
		break;
	}
	case WRITEFAIL_: {
		vector<double>delta{ 0.111*0.175 / sqrt(2),  0.111*0.175,  0.111*0.264,0.111*0.175 / sqrt(2),  0.111*0.175, 0.111*0.264 };
		for (int i = 0; i < varValue.size(); i++)
		{
			for (int j = 0; j < numSample; j++)
			{
				Res_paramValue[i][j] = paramMean[i] + varValue[i][j] * delta[i] + cor_v[j];
			}
		}
		break;
	}
	default: {
		cout << "in function alter,wrong simulation type,error!\n";
		break;
	}
	}
	return 1;
}
/*
int numSample = varValue[0].size();
if (simtype == SAFAIL_) {
vector<double>delta(3); delta[0] = 0.111*0.175 / sqrt(2); delta[1] = 0.111*0.175; delta[2] = 0.111*0.264;
for (int i = 0; i < varValue.size(); i++)
{
for (int j = 0; j < numSample; j++)
{
Res_paramValue[i][j] = paramMean[i] + varValue[i][j] * delta[i];
}
}
return 1;
}
else {
if (simtype == READFAIL_) {
vector<double>delta{ 0.111*0.175 / sqrt(2),  0.111*0.175,  0.111*0.264,0.111*0.175 / sqrt(2),  0.111*0.175, 0.111*0.264 };
for (int i = 0; i < varValue.size(); i++)
{
for (int j = 0; j < numSample; j++)
{
Res_paramValue[i][j] = paramMean[i] + varValue[i][j] * delta[i];
}
}
return 1;
}
}
};*/

bool simulator::simple_solve(vector<vector<double> >cur_M, vector<double>cur_J, vector<double> &V) {
	// solve M (M1*M2) * V (V1,1) = J(J1,1)

	/*
	string tmp;
	tmp = "cur_M"; display_matrix_vector(cur_M, tmp);
	tmp = "cur_J"; display_matrix_vector(cur_J, tmp);
	tmp = "V"; display_matrix_vector(V, tmp);
	*/

	int M1 = cur_M.size(), M2 = cur_M[0].size();
	int J1 = cur_J.size();
	int V1 = V.size();
	//cout << "M1 is " << M1 << "J1 is " << J1 << "V1 is " << V1 << endl;
	if (M2 != V1 || J1 != M1) {
		cout << "Dimension error! cannot solve !";
		return 0;
	}
	Eigen::MatrixXd  M(M1, M2);
	Eigen::MatrixXd J(J1, 1);
	for (int i = 0; i < M1; i++) {
		for (int j = 0; j < M2; j++) {
			M(i, j) = cur_M[i][j];
		}
	}
	for (int i = 0; i < J1; i++) {
		J(i) = cur_J[i];
	}
	Eigen::MatrixXd x(V1, 1);
	x = M.colPivHouseholderQr().solve(J);
	//x = M.jacobiSvd(ComputeThinU | ComputeThinV).solve(J);

	for (int i = 0; i < x.rows(); i++) {
		V[i] = x(i);
	}

	return 1;
};


bool simulator::display(int *array, string name) {
	cout << name << "\t";
	for (int i = 0; i < NumTrans; i++) {
		cout << array[i] << "\t";
	}
	cout << endl;
	return 1;
};

bool simulator::display(double *array, string name) {
	cout << name << "\t";
	for (int i = 0; i < NumTrans; i++) {
		cout << array[i] << "\t";
	}
	cout << endl;
	return 1;
};

bool simulator::displayall() {
	string s;
	s = "M_ND_"; simulator::display(M_ND_, s);
	s = "M_NG_"; simulator::display(M_NG_, s);
	s = "M_NS_"; simulator::display(M_NS_, s);
	s = "M_TYPE_"; simulator::display(M_TYPE_, s);
	s = "M_VT_"; simulator::display(M_VT_, s);
	s = "M_W_"; simulator::display(M_W_, s);
	s = "M_L_"; simulator::display(M_L_, s);
	s = "M_MU_"; simulator::display(M_MU_, s);
	s = "M_COX_"; simulator::display(M_COX_, s);
	s = "M_LAMBDA_"; simulator::display(M_LAMBDA_, s);
	return 1;
};

bool simulator::readsim(vector<vector<double> > &varValue, vector<double> &delta_V) {
	//cout<<"in readsim : "<<varValue.size()<<"\t"<<varValue[0].size()<<endl;
	if (delta_V.size() != varValue[0].size()) {
		cout << "entering function readsim,dimension is wrong,error!\n";
		return 0;
	}
	simulator::reset_M_VT_all();
	int numSample = varValue[0].size();
	vector<double> paramMean(2 * NumTrans); for (int i = 0; i < 2 * NumTrans; i++) { paramMean[i] = M_VT_[i]; }
	vector<vector<double> >paramValue(2 * NumTrans, vector<double>(numSample));
	//cout<<varValue.size()<<varValue[0].size()<<endl;
	alter(paramMean, varValue, paramValue, READFAIL_);
	for (int sampleInd = 0; sampleInd < numSample; sampleInd++) {
		vector<double>paramValueCur(2 * NumTrans); for (int i = 0; i < 2 * NumTrans; i++) { paramValueCur[i] = paramValue[i][sampleInd]; };
		//tmp_pri = "paramValueCur"; display_matrix_vector(paramValueCur, tmp_pri);
		vector <double>temp = paramValueCur;
		temp[PL_] = -temp[PL_];
		temp[PR_] = -temp[PR_];
		int equalzero = 0;
		for (int i = 0; i <2 * NumTrans; i++) { equalzero = (temp[i] > 0 ? (equalzero + 1) : equalzero); };
		int indSimFail = (equalzero == 2 * NumTrans ? 0 : 1);
		//cout << "indSimFail :" << indSimFail << endl;
		if (indSimFail == 1) {
			delta_V[sampleInd] = NAN;
			continue;
		}
		else {
			double vtripL = VDD / 2.0, vread = 0;
			simulator::set_M_VT_all(paramValueCur);
			double beta_NL = M_MU_[NL_] * M_COX_[NL_] * M_W_[NL_] / M_L_[NL_];
			double beta_AXL = M_MU_[AXL_] * M_COX_[AXL_] * M_W_[AXL_] / M_L_[AXL_];
			double beta_PL = M_MU_[PL_] * M_COX_[PL_] * M_W_[PL_] / M_L_[PL_];
			double beta_NR = M_MU_[NR_] * M_COX_[NR_] * M_W_[NR_] / M_L_[NR_];
			double beta_AXR = M_MU_[AXR_] * M_COX_[AXR_] * M_W_[AXR_] / M_L_[AXR_];
			double beta_PR = M_MU_[PR_] * M_COX_[PR_] * M_W_[PR_] / M_L_[PR_];
			double vth_NL = M_VT_[NL_];
			double vth_AXL = M_VT_[AXL_];
			double vth_PL = M_VT_[PL_];
			double vth_AXR = M_VT_[AXR_];
			double vth_NR = M_VT_[NR_];
			//map to equations parameters
			double beta1 = 0.5*beta_NL, beta2 = 0.5*beta_AXL, beta3 = 0.5*beta_PL;
			double tmp_a1 = vth_NL, tmp_a2 = VDD - vth_AXL, tmp_a3 = VDD + vth_PL;
			double tmp_A = beta1 - beta2 - beta3;
			double tmp_B = -(2 * tmp_a1*beta1 - 2 * tmp_a2*beta2 - 2 * tmp_a3*beta3);
			double tmp_C = beta1*tmp_a1*tmp_a1 - beta2*tmp_a2*tmp_a2 - beta3*tmp_a3*tmp_a3;
			//solve Ax^2+Bx+C=0
			double candi[2] = { VDD / 2.0 ,VDD / 2.0 };
			if (helperfunc::Mysolve(tmp_A, tmp_B, tmp_C, candi) == false) {
				delta_V[sampleInd] = NAN;
				continue;
			}
			//check VtripL candidats' validity
			if (candi[0] != candi[1]) {
				bool vali1 = (candi[0] > vth_NL && VDD - candi[0] > vth_AXL && candi[0] - VDD < vth_PL);
				bool vali2 = (candi[1] > vth_NL && VDD - candi[1] > vth_AXL && candi[1] - VDD < vth_PL);
				if (vali1 == true && vali2 == false) {
					vtripL = candi[0];
				}
				if (vali1 == false && vali2 == true) {
					vtripL = candi[1];
				}
				if (vali1 == false && vali2 == false) {
					delta_V[sampleInd] = NAN;
					continue;
				}
				if (vali1 == true && vali2 == true) {
					cout << "[Error]: two possible solutions for VtripL\n";
					delta_V[sampleInd] = NAN;
					continue;
				}
			}
			else {
				bool vali = (candi[0] > vth_NL && VDD - candi[0] > vth_AXL && candi[0] - VDD < vth_PL);
				if (vali == true) {
					vtripL = candi[0];
				}
				else {
					delta_V[sampleInd] = NAN;
					continue;
				}
			}
			//vtripL = 0.14;
			vtripL *= MYRATE;
			//solve vread,the following parameters are based on : IdsatAXR=IdlinNR,IdlinAXL=IdlinPL
			double a1 = 0.5*beta_AXR, a2 = VDD - vth_AXR, a3 = beta_NR, a4 = VDD - vth_NR;
			tmp_A = a1 + 0.5*a3, tmp_B = -(2 * a1*a2 + a3*a4), tmp_C = a1*a2*a2;
			double candi_R[2] = { 0,0 };
			if (helperfunc::Mysolve(tmp_A, tmp_B, tmp_C, candi_R) == false) {
				delta_V[sampleInd] = NAN;
				continue;
			}
			//check Vread validity
			if (candi_R[0] != candi_R[1]) {
				bool vali1 = (VDD - candi_R[0] > vth_AXR && VDD > vth_NR && vth_NR < VDD - candi_R[0]);
				bool vali2 = (VDD - candi_R[1] > vth_AXR && VDD > vth_NR && vth_NR < VDD - candi_R[1]);
				if (vali1 == true && vali2 == true) {
					cout << "[Error]: two possible solutions for Vread\n";
					delta_V[sampleInd] = NAN;
					continue;
				}
				if (vali1 == false && vali2 == false) {
					delta_V[sampleInd] = NAN;
					continue;
				}
				if (vali1 == true && vali2 == false) {
					vread = candi_R[0];
				}
				if (vali1 == false && vali2 == true) {
					vread = candi_R[1];
				}
			}
			else {
				bool vali = (VDD - candi_R[0] > vth_AXR && VDD > vth_NR && vth_NR < VDD - candi_R[0]);
				if (vali == true)
				{
					vread = candi_R[0];
				}
				else {
					delta_V[sampleInd] = NAN;
					continue;
				}
			}
			//cout << "vtripL : " << vtripL << " vread : " << vread << endl;
			delta_V[sampleInd] = vtripL - vread;
		};
	}
	return 1;
};

bool simulator::reset_M_VT_all() {
	M_VT_[0] = 0.1750, M_VT_[1] = 0.1750, M_VT_[2] = -0.2640;
	M_VT_[3] = 0.1750, M_VT_[4] = 0.1750, M_VT_[5] = -0.2640;
	return true;
};

bool simulator::writesim(vector<vector<double> > &varValue, vector<double> &delta_t) {
	if (delta_t.size() != varValue[0].size()) {
		cout << "entering function writesim,dimension is wrong,error!\n";
		return 0;
	}
	simulator::reset_M_VT_all();
	int numSample = varValue[0].size();
	vector<double> paramMean(2 * NumTrans); for (int i = 0; i < 2 * NumTrans; i++) { paramMean[i] = M_VT_[i]; }
	vector<vector<double> >paramValue(2 * NumTrans, vector<double>(numSample));
	alter(paramMean, varValue, paramValue, WRITEFAIL_);
	for (int sampleInd = 0; sampleInd < numSample; sampleInd++) {
		vector<double>paramValueCur(2 * NumTrans); for (int i = 0; i < 2 * NumTrans; i++) { paramValueCur[i] = paramValue[i][sampleInd]; };
		//tmp_pri = "paramValueCur"; display_matrix_vector(paramValueCur, tmp_pri);
		vector <double>temp = paramValueCur;
		temp[PL_] = -temp[PL_];
		temp[PR_] = -temp[PR_];
		int equalzero = 0;
		for (int i = 0; i < 2 * NumTrans; i++) { equalzero = (temp[i] > 0 ? (equalzero + 1) : equalzero); };
		int indSimFail = (equalzero == 2 * NumTrans ? 0 : 1);
		//cout << "indSimFail :" << indSimFail << endl;
		if (indSimFail == 1) {
			delta_t[sampleInd] = NAN;
			continue;
		}
		else
		{
			simulator::set_M_VT_all(paramValueCur);
			double vtripR = VDD / 2.0;
			// solve vtripR
			double beta_NL = M_MU_[NL_] * M_COX_[NL_] * M_W_[NL_] / M_L_[NL_];
			double beta_AXL = M_MU_[AXL_] * M_COX_[AXL_] * M_W_[AXL_] / M_L_[AXL_];
			double beta_PL = M_MU_[PL_] * M_COX_[PL_] * M_W_[PL_] / M_L_[PL_];
			double beta_NR = M_MU_[NR_] * M_COX_[NR_] * M_W_[NR_] / M_L_[NR_];
			double beta_AXR = M_MU_[AXR_] * M_COX_[AXR_] * M_W_[AXR_] / M_L_[AXR_];
			double beta_PR = M_MU_[PR_] * M_COX_[PR_] * M_W_[PR_] / M_L_[PR_];
			double vth_NL = M_VT_[NL_];
			double vth_AXL = M_VT_[AXL_];
			double vth_PL = M_VT_[PL_];
			double vth_AXR = M_VT_[AXR_];
			double vth_NR = M_VT_[NR_];
			double vth_PR = M_VT_[PR_];
			//map to equations parameters
			double beta1 = 0.5*beta_NR, beta2 = 0.5*beta_AXR, beta3 = 0.5*beta_PR;
			double tmp_a1 = vth_NR, tmp_a2 = VDD - vth_AXR, tmp_a3 = VDD + vth_PR;
			double tmp_A = beta1 - beta2 - beta3;
			double tmp_B = -(2 * tmp_a1*beta1 - 2 * tmp_a2*beta2 - 2 * tmp_a3*beta3);
			double tmp_C = beta1*tmp_a1*tmp_a1 - beta2*tmp_a2*tmp_a2 - beta3*tmp_a3*tmp_a3;
			//solve Ax^2+Bx+C=0
			double candi[2] = { VDD / 2.0 ,VDD / 2.0 };
			if (helperfunc::Mysolve(tmp_A, tmp_B, tmp_C, candi) == false) {
				delta_t[sampleInd] = NAN;
				continue;
			}
			//check VtripR candidats' validity
			if (candi[0] != candi[1]) {
				bool vali1 = (candi[0] > vth_NR && VDD - candi[0] > vth_AXR && candi[0] - VDD < vth_PR);
				bool vali2 = (candi[1] > vth_NR && VDD - candi[1] > vth_AXR && candi[1] - VDD < vth_PR);
				if (vali1 == true && vali2 == false) {
					vtripR = candi[0];
				}
				if (vali1 == false && vali2 == true) {
					vtripR = candi[1];
				}
				if (vali1 == false && vali2 == false) {
					delta_t[sampleInd] = NAN;
					continue;
				}
				if (vali1 == true && vali2 == true) {
					cout << "[Error]: two possible solutions for VtripL\n";
					delta_t[sampleInd] = NAN;
					continue;
				}
			}
			else {
				bool vali = (candi[0] > vth_NR && VDD - candi[0] > vth_AXR && candi[0] - VDD < vth_PR);
				if (vali == true) {
					vtripR = candi[0];
				}
				else {
					delta_t[sampleInd] = NAN;
					continue;
				}
			}
			//vtripR = 0.14;
			//vtripR *= MYRATE;
			double T = cal_integral(VDD, vtripR, 1000);
			delta_t[sampleInd] = TWL - T;
			//cout << "TWL : " << TWL << " T : " << T << endl;
		}
	}
	return 1;
}

//[function]Icurrent : calculate Ids of MOSFET given vg,vd,vs.
double simulator::Icurrent(int MOS_num, int MOS_type, double Vg, double Vd, double Vs) {
	double beta = M_MU_[MOS_num] * M_COX_[MOS_num] * M_W_[MOS_num] / M_L_[MOS_num];
	if (MOS_type == NMOS_ && Vs > Vd) {
		cout << "[Warning]:NMOS drain source exchange! Vd " << Vd << " Vs " << Vs << endl;
		double tmp = Vs;
		Vs = Vd;
		Vd = tmp;
	}
	if (MOS_type == PMOS_ && Vs < Vd) {
		cout << "[Warning]:PMOS drain source exchange! Vd " << Vd << " Vs " << Vs << endl;
		double tmp = Vs;
		Vs = Vd;
		Vd = tmp;
	}
	double vgs = Vg - Vs, vds = Vd - Vs, vth = M_VT_[MOS_num];
	if ((vgs < vth && MOS_type == NMOS_) || (vgs > vth && MOS_type == PMOS_)) {
		//cout << " off";
		return 0; // off 
	}
	if ((vgs > vth && vds <= (vgs - vth) && MOS_type == NMOS_) || (vgs < vth && vds >= (vgs - vth) && MOS_type == PMOS_)) {
		//cout << " linear";
		return beta*((vgs - vth)*vds - 0.5*vds*vds);// linear region
	}
	else {
		//cout << " saturation";
		return 0.5*beta*(vgs - vth)*(vgs - vth);// saturation region
	}
};

double simulator::Capitance(int MOS_num, int MOS_type, double Vg, double Vd, double Vs, int C_type) {
	double vgs = Vg - Vs, vds = Vd - Vs, vth = M_VT_[MOS_num];
	if ((vgs < vth && MOS_type == NMOS_) || (vgs > vth && MOS_type == PMOS_)) {
		return M_W_[MOS_num] * M_COV_[MOS_num];// off state
	}
	if ((vgs > vth && vds <= (vgs - vth) && MOS_type == NMOS_) || (vgs < vth && vds >= (vgs - vth) && MOS_type == PMOS_)) {
		return 0.5*M_W_[MOS_num] * M_L_[MOS_num] * M_COX_[MOS_num] + M_W_[MOS_num] * M_COV_[MOS_num];// linear region
	}
	else {
		if (C_type == GD) {
			return M_W_[MOS_num] * M_COV_[MOS_num];
		}
		else {
			return M_W_[MOS_num] * M_COV_[MOS_num] + 2 / 3 * M_W_[MOS_num] * M_L_[MOS_num] * M_COX_[MOS_num];
		}
	}
};

double simulator::cal_integral(double lower, double upper, int intervals) {
	double step = (upper - lower) / intervals;
	double res = 0;
	for (int i = 0; i < intervals; i++) {
		double cur_V = lower + i*step;
		res += helper_integral(cur_V)*step;
	}
	return abs(res);
};

double simulator::helper_integral(double cur_v) {
	//cout << "PMOS region:";  double I1 = Icurrent(PL_, PMOS_, 0, cur_v, VDD); cout << "current : "<<I1 << endl;
	//cout << "NMOS region:"; ; double I2 = Icurrent(AXL_, NMOS_, VDD, cur_v, 0); cout << "current : "<<I2 << endl;
	double I = Icurrent(PL_, PMOS_, 0, cur_v, VDD) - Icurrent(AXL_, NMOS_, VDD, cur_v, 0);
	//double I = I1 - I2;
	//double C = Capitance(PL_, PMOS_, 0, cur_v, VDD, GD);// + Capitance(AXL_, NMOS_, VDD, cur_v, 0, GD);// +Capitance(NL_, NMOS_, 0, cur_v, 0, GD);
	double C = Capitance(AXL_, NMOS_, VDD, cur_v, 0, GD);// +Capitance(NL_, NMOS_, 0, cur_v, 0, GD);
	return C / I;
};


bool simulator::addcorelation(int numSample, vector<double> &cor_v) {
	//cout<<"entering addcorelation,numSample is : "<<numSample<<endl;
	for (int i = 0; i < numSample; i++) {
		//cor_v[i]=n(e)*0.111*0.2;
		cor_v[i] = 0;
		//cout<<cor_v[i];
	}
	return 1;
}
