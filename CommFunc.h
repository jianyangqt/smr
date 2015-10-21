/*
 * Interface to the commonly-used functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#ifndef _COMMFUNC_H
#define _COMMFUNC_H

#define MAX_BUF_SIZE 0x40000000
#ifdef _WIN64
#define MAX_LINE_SIZE 655360
#else
#define MAX_LINE_SIZE 0x00753000
#endif

#define MAX_SNP_NAME 64
#define DENSE_FILE_TYPE_1 0  // float + floats
#define SPARSE_FILE_TYPE_1 1 // float + float + floats + floats +floats
#define SPARSE_FILE_TYPE_2 2 // float + uint64_t + uint32_ts + uint32_ts + floats
#define DENSE_FILE_TYPE_2 3

#include <limits>
#include <vector>
#include <ctime>
#include <fstream>
#include "StrFunc.h"
#include <string.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/Eigenvalues>

typedef unsigned long long         uint64_t;
typedef unsigned int         uint32_t;


using namespace Eigen;
using namespace std;

namespace CommFunc
{
	const double FloatErr=numeric_limits<double>::epsilon();
	double Abs(const double &x);
	double sum(const vector<double> &x);
    double mean(const vector<double> &x);
    double median(const vector<double> &x);
    double var(const vector<double> &x);
    double cov(const vector<double> &x, const vector<double> &y);
	bool FloatEqual(double lhs, double rhs);
	bool FloatNotEqual(double lhs, double rhs);
	const double Sqr(const double &a);
	const double Max(const double &a, const double &b);
	const double Min(const double &a, const double &b);
	const double Sign(const double &a, const double &b);
	int rand_seed(); //positive value, return random seed using the system time
    void FileExist(string filename);  
    int max_abs_id(VectorXd &zsxz);
    int max_abs_id(vector<double> &zsxz);
    int fopen_checked(FILE** in_file, const char* filename, const char* flag);
    void getRank(vector<int> &a, vector<int> &b);
}

#endif
