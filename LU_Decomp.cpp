#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

typedef std::vector<double> dvector;

void lu_decomp(dvector const &mat, dvector &lu_mat, int dim)
{
	for (int i = 0; i < dim; ++i) {
		for (int j = i; j < dim; ++j) {
			lu_mat[i * dim + j] = mat[i * dim + j];
			for (int k = 0; k < i; ++k)
				lu_mat[i * dim + j] -= lu_mat[i * dim + j] * lu_mat[k * dim + j];
			if (j == i) continue;
			lu_mat[j* dim +i] = mat[j* dim +i];
			for (int k = 0; k < i; ++k) {
				if (j == k)
					lu_mat[j * dim + i] -= lu_mat[k* dim +i];
				else
					lu_mat[j* dim +i] -= lu_mat[j* dim +k] * lu_mat[k* dim +i];
			}
			lu_mat[j* dim + i] /= lu_mat[i* dim +i];
		}
	}
}

void solve_lu(dvector &vec_x, dvector &mat, dvector &vec_b,
	dvector &tmp_mat, dvector &tmp_vec, int dim)
{
	lu_decomp(mat, tmp_mat, dim);
	for (long i = 0; i < dim; ++i) {
		double val = vec_b[i];
		for (long j = 0; j < i; ++j)
			val -= tmp_mat[i * dim + j] * tmp_vec[j];
		tmp_vec[i] = val;
	}
	for (long i = dim - 1; i >= 0; --i) {
		double val = tmp_vec[i];
		for (long j = dim - 1; j > i; --j)
			val -= tmp_mat[i * dim + j] * vec_x[j];
		vec_x[i] = val / tmp_mat[i * dim +i];
	}
}

void mult(dvector &res, dvector const &mat, dvector const &vec)
{
	int size = vec.size();
	for (int i = 0; i < size; ++i) {
		res[i] = 0;
		for (int j = 0; j < size; ++j)
			res[i] += mat[i * size + j] * vec[j];
	}
}

void add(dvector &res, dvector const &v1, dvector const &v2)
{
	for (int i = 0; i < v1.size(); ++i)
		res[i] = v1[i] + v2[i];
}

void sub(dvector &res, dvector const &v1, dvector const &v2)
{
	for (int i = 0; i < v1.size(); ++i)
		res[i] = v1[i] - v2[i];
}

void copy(dvector &res, dvector const &v1)
{
	for (int i = 0; i < v1.size(); ++i)
		res[i] = v1[i];
}

double euclid_norm(dvector const &vec)
{
	double val = 0;
	for (int i = 0; i < vec.size(); i++)
		val += vec[i] * vec[i];
	return std::sqrt(val);
}

void my_func(dvector &res, dvector const &arg, dvector const &mat)
{
	mult(res, mat, arg);
	for (int i = 0; i < res.size(); ++i)
		res[i] -= std::exp(-arg[i]);
}

void my_jakobian(dvector &jk, dvector const &arg, dvector const &mat, int dim)
{
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			jk[i * dim + j] = mat[i * dim + j];
			if (i == j) jk[i * dim + j] += std::exp(-arg[i]);
		}
	}
}

void newton_iter(dvector &next, dvector const &prev, dvector const &mat,
		dvector &tmp_mat_jk, dvector &tmp_mat_lu,
		dvector &tmp_vec_b, dvector &tmp_vec_lu)
{
	my_func(tmp_vec_b, prev, mat);
	my_jakobian(tmp_mat_jk, prev, mat, prev.size());
	solve_lu(next, tmp_mat_jk, tmp_vec_b, tmp_mat_lu, tmp_vec_lu, prev.size());
	sub(next, prev, next);
}

double residual = 1e-06;

dvector import(char const *path)
{
	std::ifstream file;
	file.open(path);

	int size;
	file >> size;
	dvector mat(size * size);

	for (int i = 0; i < mat.size(); ++i) {
		double val;
		file >> val;
		mat[i] = val;
	}	
	file.close();
	return mat;
}

int main(int argc, char *argv[])
{
	if (argc != 2)
		return 1;

	dvector mat = import(argv[1]);
	int size = 2;

	dvector vec_x(size);
	dvector next(size);
	for (int i = 0; i < size; ++i)
		vec_x[i] = 1;

	dvector tmp_mat_jk(size * size);
	dvector tmp_mat_lu(size * size);
	dvector tmp_vec_b (size);
	dvector tmp_vec_lu(size);

	double res;
	do {
		newton_iter(next, vec_x, mat, tmp_mat_jk, tmp_mat_lu,
			tmp_vec_b, tmp_vec_lu);
		std::cout << "Vector: ";
		for (int i = 0; i < next.size(); ++i)
			std::cout << next[i] << " ";
		std::cout << std::endl;
		sub(vec_x, next, vec_x);
		res = euclid_norm(vec_x);
		std::cout << "Residual: " << res << std::endl;
		copy(vec_x, next);
	} while (res > residual);
	return 0;
}
