#include <iostream>
#include <vector>
#include <cmath>

#define N_STEPS 20.0
#define N_PERIODS 3.0

#define U0 1.0
#define V0 0.0

#define R 0.3
#define L 2.0
#define G 4.0

using namespace std;

vector<double> operator+(vector<double> const &v1, vector<double> const &v2)
{
	vector<double> result(v1.size());
	for (size_t i = 0; i < result.size(); ++i)
		result[i] = v1[i] + v2[i];
	return result;
}


vector<double> operator*(double a, vector<double> const &v)
{
	vector<double> result(v.size());
	for (size_t i = 0; i < result.size(); ++i)
		result[i] = a*v[i];
	return result;
}

class Experiment {
public:
	//friend Experiment operator+()
	vector<double> FuncCalc(const vector<double> &cur_v) {
		vector<double> new_v (2);
		new_v[0] = cur_v[1];
		new_v[1] = -2*R/L*cur_v[1] - G/L*sin(cur_v[0]);
		return new_v;
	}	


	vector<double> EulerMethod(const vector<double> &cur_vec)
      	{
		vector<double> next_vec (2);
		const double t = CalcStep(); 

		next_vec[0] = cur_vec[0] + t*cur_vec[1];
		next_vec[1] = cur_vec[1] - t*(2*R*cur_vec[1] + G*sin(cur_vec[0]))/L;
		return next_vec;
	}



	vector<double> RK_Method(const vector<double> &cur_vec)
	{
		vector<double> next_vec (2);

		vector<double> k1 (2);
		vector<double> k2 (2);
		vector<double> k3 (2);
		vector<double> k4 (2);

		k1 = FuncCalc(cur_vec);

		const double t = CalcStep();

		k2 = (2+t)/2*k1;
		k3 = (4+2*t+t*t)/4*k1;
		k4 = (4+4*t+2*t*t+t*t*t)/4*k1;

		next_vec = cur_vec + t/6*(k1 + 2*k2 + 2*k3 + k4);

		return next_vec;
	}

	double SolutionEq(double t)
       	{
		double q = R/L;
		double w = G/L - q*q;
		w = sqrt(w);
		return exp(-q*t)*cos(w*t);
	}

	double CalcStep() {
		return (2*M_PI*sqrt(L/G))/N_STEPS;
	}
	
};


int main() 
{
	Experiment exper;

	vector<double> em_cur (2);

	em_cur[0] = U0;
	em_cur[1] = V0;

	vector<double> rk_cur (2);

	rk_cur[0] = U0;
	rk_cur[1] = V0;

	double t = exper.CalcStep();

	for(int i = 0; i < N_STEPS*N_PERIODS ; i++) {
		vector<double> em_new = exper.EulerMethod(em_cur);
		em_cur = em_new;

		vector<double> rk_new = exper.RK_Method(rk_cur);
		rk_cur = rk_new;

		double sol = exper.SolutionEq(t*i);
		cout << "time: " << t*i <<" em: " << em_cur[0] << " rk: "
	       	     << rk_cur[0] << " sol: " << sol << endl;
	}

	return 0;
}
