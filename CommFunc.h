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
#define MAX_LINE_SIZE 0x00753000
#define MAX_SNP_NAME 64


#include <limits>
#include <complex>
#include <vector>
#include <algorithm>
#include <ctime>
#include <fstream>
#include "StrFunc.h"

#include <string.h>

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
    int max_abs_id(vector<double> &zsxz);
    int fopen_checked(FILE** in_file, const char* filename, const char* flag);
    void getRank(vector<int> &a, vector<int> &b);
}

#endif
