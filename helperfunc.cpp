#include"helperfunc.h"

bool helperfunc::write_to_txt(vector<vector<double> >&Fp, const char *file) {
	FILE *fp = NULL;
	fp = fopen(file, "w+");
	for (int i = 0; i < Fp.size(); i++) {
		for (int j = 0; j < Fp[0].size(); j++) {
			fprintf(fp, "%.12f", Fp[i][j]);
			fprintf(fp, " ");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return 1;
}

bool helperfunc::write_to_txt(vector<double> &Fp, const char *file) {
	FILE *fp = NULL;
	fp = fopen(file, "w+");
	for (int i = 0; i < Fp.size(); i++) {
		{
			fprintf(fp, "%E", Fp[i]);
			fprintf(fp, " ");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return 1;
}

int helperfunc::combntns(int maxx, int minn) {
	minn = minn > maxx - minn ? maxx - minn : minn;
	int res = 1;
	for (int i = maxx; i >= maxx - minn + 1; i--) {
		res = res*i;
	}
	for (int i = minn; i >= 1; i--) {
		res /= i;
	}
	return res;
}

void helperfunc::normpdf(double *src, int size, double mu, double sigma, double *dst) {
	for (int i = 0; i < size; i++) {
		dst[i] = exp(-0.5*(src[i] - mu) / sigma*(src[i] - mu) / sigma) / (sqrt(2 * PI)*sigma);
	}
}

double helperfunc::erf(double x)
{
	double y = 1.0 / (1.0 + 0.3275911 * x);
	return 1 - (((((
		+1.061405429  * y
		- 1.453152027) * y
		+ 1.421413741) * y
		- 0.284496736) * y
		+ 0.254829592) * y)
		* exp(-x * x);
};

void helperfunc::normcdf(double *src, int size, double mu, double sigma, double *dst)
{
	// constants
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p = 0.3275911;
	for (int i = 0; i < size; i++) {
		src[i] = (src[i] - mu) / fabs(sigma);
		int sign = 1;
		if (src[i] < 0)
			sign = -1;
		src[i] = fabs(src[i]) / sqrt(2.0);
		double t = 1.0 / (1.0 + p*src[i]);
		double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-src[i] * src[i]);
		dst[i] = 0.5*(1.0 + sign*y);
	}
}

bool helperfunc::choose_ref(vector<double>&src, vector<double> &dst, vector<bool>ref, int axis) {
	// choose elements in vector src,according to vector ref, i.e. dst=src(ref) in matlab
	if (src.size() != ref.size()) {
		cout << "In function choose_ref, the dimension of the vector is not same, error !\n";
		return 0;
	}
	dst.clear();
	for (int i = 0; i < ref.size(); i++) {
		if (ref[i] == 1) {
			dst.push_back(src[i]);
		}
	}
	return 1;
}

bool helperfunc::choose_ref(vector<vector<double> >&src, vector<vector<double> >&dst, vector<bool>ref, int axis) {
	// choose elements in matrix src,according to vector ref, i.e. if axis=1,dst=src(ref,:) ; if axis=2,dst=src(:,ref) in matlab
	if (axis == 1) {
		if (ref.size() != src.size()) {
			cout << "In function choose_ref, the dimension of the matrix and vector is not same, error !\n";
			return 0;
		}
		else {
			dst.clear();
			for (int i = 0; i < ref.size(); i++) {
				if (ref[i] == 1) {
					dst.push_back(src[i]);
				}
			}
			return 1;
		}
	}
	else {
		if (axis == 2) {
			if (ref.size() != src[0].size()) {
				cout << "In function choose_ref, the dimension of the matrix and vector is not same, error !\n";
				return 0;
			}
			else {
				dst.clear();
				int newlength = 0;
				for (int i = 0; i < ref.size(); i++) {
					if (ref[i] == 1) {
						newlength += 1;
					}
				}
				vector<vector<double> > tmp(src.size(), vector<double>(newlength));
				for (int i = 0, cur = 0; i < newlength; i++) {
					while (ref[cur] != 1) { cur++; }
					for (int j = 0; j < src.size(); j++) {
						tmp[j][i] = src[j][cur];
					}
					cur++;
				}
				dst = tmp;
				return 1;
			}
		}
		else {
			cout << "In function choose_ref, parameter axis should be 1 or 2, error !\n";
		}
	}
	return 1;
};

bool helperfunc::expand(vector<vector<double> >&src1, vector<vector<double> >&src2, vector<vector<double> >&dst, int axis) {
	//axis==1: dst=[src1 src2] ; axis==2 ; dst=[src1 ; src2];

	// src1==NULL src2==NULL
	if (src1.size() == 0 && src2.size() == 0) {
		cout << "in function expand, both src are NULL\n";
		return 0;
	}
	//src1==NULL src2!=NULL
	if (src1.size() == 0 && src2.size() != 0) {
		int dst1 = dst.size(), dst2 = dst[0].size();
		int src2_1 = src2.size(), src2_2 = src2[0].size();
		if (dst1 != src2_1 || dst2 != src2_2) {
			cout << "in function expand,dimension error\n";
		}
		for (int i = 0; i < dst1; i++) {
			for (int j = 0; j < dst2; j++) {
				dst[i][j] = src2[i][j];
			}
		}
		return 1;
	}
	//src1!=NULL src2==NULL
	if (src1.size() != 0 && src2.size() == 0) {
		int dst1 = dst.size(), dst2 = dst[0].size();
		int src1_1 = src1.size(), src1_2 = src1[0].size();
		if (dst1 != src1_1 || dst2 != src1_2) {
			cout << "in function expand,dimension error\n";
		}
		for (int i = 0; i < src1_1; i++) {
			for (int j = 0; j < src1_2; j++) {
				dst[i][j] = src1[i][j];
			}
		}
		return 1;
	}
	int src1_1 = src1.size(), src1_2 = src1[0].size(), src2_1 = src2.size(), src2_2 = src2[0].size(), dst_1 = dst.size(), dst_2 = dst[0].size();
	if (axis == 1) {
		if (src1_1 != src2_1 || src1_1 != dst_1 || dst_1 != src2_1 || src1_2 + src2_2 != dst_2) {
			cout << "enter function expand,dimension is wrong,error!\n";
			return 0;
		}
		for (int i = 0; i < dst_1; i++) {
			for (int j = 0; j < dst_2; j++) {
				dst[i][j] = j < src1_2 ? src1[i][j] : src2[i][j - src1_2];
			}
		}
	}
	else {
		if (axis == 2) {
			if (src1_2 != src2_2 || src1_2 != dst_2 || dst_2 != src2_2 || src1_1 + src2_1 != dst_1) {
				cout << "enter function expand,dimension is wrong,error!\n";
				return 0;
			}
			for (int i = 0; i < dst_1; i++) {
				for (int j = 0; j < dst_2; j++) {
					dst[i][j] = i < src1_1 ? src1[i][j] : src2[i - src1_1][j];
				}
			}
		}
		else {
			cout << "enter function expand,axis is wrong,error!\n";
			return 0;
		}
	}
	return 1;
};

bool helperfunc::delete_ref(vector<double>&src, vector<bool>ref, int axis) {
	// delete elements in vector src,according to vector ref, i.e. src(ref)=[] in matlab
	if (src.size() != ref.size()) {
		cout << "In function delete_ref, the dimension of the vector is not same, error !\n";
		return 0;
	}
	vector<double>dst;
	for (int i = 0; i < ref.size(); i++) {
		if (ref[i] == 0) {
			dst.push_back(src[i]);
		}
	}
	src = dst;
	return 1;
}

bool helperfunc::delete_ref(vector<vector<double> >&src, vector<bool>ref, int axis) {
	// delete elements in matrix src,according to vector ref, i.e. if axis=1,src(ref,:)=[] ; if axis=2,src(:,ref)=[] in matlab
	if (axis == 1) {
		if (ref.size() != src.size()) {
			cout << "In function delete_ref, the dimension of the matrix and vector is not same, error !\n";
			return 0;
		}
		else {
			vector<vector<double> > dst;
			for (int i = 0; i < ref.size(); i++) {
				if (ref[i] == 0) {
					dst.push_back(src[i]);
				}
			}
			src = dst;
			return 1;
		}
	}
	else {
		if (axis == 2) {
			if (ref.size() != src[0].size()) {
				cout << "In function delete_ref, the dimension of the matrix and vector is not same, error !\n";
				return 0;
			}
			else {
				int newlength = 0;
				for (int i = 0; i < ref.size(); i++) {
					if (ref[i] == 0) {
						newlength += 1;
					}
				}
				vector<vector<double> > dst(src.size(), vector<double>(newlength));
				for (int i = 0, cur = 0; i < newlength; i++) {
					while (ref[cur] != 0) { cur++; }
					for (int j = 0; j < src.size(); j++) {
						dst[j][i] = src[j][cur];
					}
					cur++;
				}
				src = dst;
				return 1;
			}
		}
		else {
			cout << "In function delete_ref, parameter axis should be 1 or 2, error !\n";
		}
	}
	return 1;
};

bool helperfunc::display_matrix_vector(vector<vector<double> >src, string name, int precision) {
	cout << name << endl;
	for (int i = 0; i < src.size(); i++) {
		for (int j = 0; j < src[0].size(); j++) {
			//cout << fixed << setprecision(precision)<< src[i][j] << "\t";
			cout << src[i][j] << "\t";
		}
		cout << endl << endl;
	}
	return 1;
}

bool helperfunc::display_matrix_vector(vector<double>src, string name, int precision) {
	cout << name << endl;
	for (int i = 0; i < src.size(); i++) {
		//cout << fixed << setprecision(precision) << src[i] << "\t";
		cout << src[i] << "\t";
	}
	cout << endl << endl;
	return 1;
}

double helperfunc::mean(vector<double>src) {
	double ans = 0;
	for (int i = 0; i < src.size(); i++) {
		ans += src[i];
	}
	ans /= (double(src.size()));
	return ans;
}


double helperfunc::var(vector<double> src) {
	int n = src.size();
	double ans = 0;
	for (int i = 0; i < n; i++) {
		ans += src[i];
	}
	return ans /= n;
}

bool helperfunc::Mysolve(double A, double B, double C, double *dst) {
	double deta = sqrt(B*B - 4 * A*C);
	if (A == 0) {
		if (B == 0) {
			if (C != 0) {
				cout << "[Error]: Equations Ax^2+Bx+C=0 has values A=0,B=0,C!=0\n";
				return false;
			}
			else {
				cout << "[Warning]: Equations Ax^2+Bx+C=0 has values A=0,B=0,C=0\n";
				return false;
			}
		}
		else {
			cout << "[Warning]: Equations Ax^2+Bx+C=0 has values A=0,B!=0\n";
			dst[0] = dst[1] = -C / B;
			return false;
		}
	}
	else {
		if (deta < 0) {
			cout << "[Error]: Equations Ax^2+Bx+C=0 has values A!=0, deta<0!\n";
			return false;
		}
		else {
			double tmp1 = (-B - deta) / (2 * A);
			double tmp2 = (-B + deta) / (2 * A);
			dst[0] = tmp1 < tmp2 ? tmp1 : tmp2;
			dst[1] = tmp1 < tmp2 ? tmp2 : tmp1;
			return true;
		}
	}
};