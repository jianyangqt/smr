/*
 * Implementations of the commonly-used functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "CommFunc.h"

double CommFunc::Abs(const double &x)
{
	complex<double> cld(x);
	double ldAbs = abs(cld);
	return(ldAbs);
}

double CommFunc::sum(const vector<double> &x)
{
    int size = x.size();
    int i=0;
    double d_buf=0.0;
    for(i=0; i<size; i++) d_buf+=x[i];
    return (double)d_buf;
}

double CommFunc::mean(const vector<double> &x)
{
    int size = x.size();
    int i=0;
    double d_buf=0.0;
    for(i=0; i<size; i++) d_buf+=x[i];
    d_buf/=(double)size;
    return (double)d_buf;
}

double CommFunc::median(const vector<double> &x)
{
    vector<double> b(x);
    int size = b.size();
    if(size==1) return b[0];
    stable_sort(b.begin(), b.end());
    if(size%2==1) return b[(size-1)/2];
    else return (b[size/2]+b[size/2-1])/2;
}

double CommFunc::var(const vector<double> &x)
{
    int size = x.size();
    if(size<=1) return(0.0);
    int i=0;
    double mu=0.0, s2=0.0;
    for(i=0; i<size; i++) mu+=x[i];
    mu/=(double)size;
    for(i=0; i<size; i++) s2+=(x[i]-mu)*(x[i]-mu);
    s2/=(double)(size-1);
    return (double)s2;
}

double CommFunc::cov(const vector<double> &x, const vector<double> &y)
{
    int size = x.size();
    int i=0;
    double mu1=0.0, mu2=0.0, c=0.0;
    for(i=0; i<size; i++){
        mu1+=x[i];
        mu2+=y[i];
    }
    mu1/=(double)size;
    mu2/=(double)size;

    for(i=0; i<size; i++) c+=(x[i]-mu1)*(y[i]-mu2);
    c/=(double)(size-1);
    return c;
}

bool CommFunc::FloatEqual(double lhs, double rhs)
{
	if (Abs(lhs - rhs) < FloatErr) return true;
	return false;
}

bool CommFunc::FloatNotEqual(double lhs, double rhs)
{
	if (Abs(lhs - rhs) >= FloatErr) return true;
	return false;
}

const double CommFunc::Sqr(const double &a)
{
	return a*a;
}

const double CommFunc::Max(const double &a, const double &b)
{
	return b > a ? (b) : (a);
}

const double CommFunc::Min(const double &a, const double &b)
{
	return b < a ? (b) : (a);
}

const double CommFunc::Sign(const double &a, const double &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

int CommFunc::rand_seed()
{
  	stringstream str_strm;
	str_strm<<time(NULL);
	string seed_str=str_strm.str();
    reverse(seed_str.begin(), seed_str.end());
    seed_str.erase(seed_str.begin()+7, seed_str.end());
	return(abs(atoi(seed_str.c_str())));
}

void CommFunc::FileExist(string filename)
{
    ifstream ifile(filename.c_str());
    if(!ifile) throw("Error: can not open the file ["+filename+"] to read.");
}
int CommFunc::max_abs_id(vector<double> &zsxz)
{
    int id=0;
    double tmpVal, cmpVal=abs(zsxz[0]);
    for( int i=1;i<zsxz.size();i++)
    {
        tmpVal=abs(zsxz[i]);
        if( cmpVal-tmpVal < 1e-6)
        {
            cmpVal=tmpVal;
            id=i;
        }
    }
    return(id);
}

int CommFunc::max_abs_id(VectorXd &zsxz)
{
    int id=0;
    double tmpVal, cmpVal=abs(zsxz[0]);
    for( int i=1;i<zsxz.size();i++)
    {        
        tmpVal=abs(zsxz[i]);        
        if( cmpVal-tmpVal < 1e-6)
        {
            cmpVal=tmpVal;
            id=i;
         }
    }
    return(id);
}

int CommFunc::fopen_checked(FILE** in_file, const char* filename, const char* flag)
{
    *in_file=fopen(filename,flag);
    if(!*in_file)  return 0;
    return 1;
}

void CommFunc::getRank(vector<int> &a, vector<int> &b)
{
    b.resize(a.size());
    for (int i = a.size()-1; i >= 0; i--)
    {
        int count = 0;
        for (int j = 0; j < a.size(); j++) if (a[j] < a[i]) count++;
        b[i] = count;
    }
}
void CommFunc::getRank_norep(vector<int> &a, vector<int> &b)
{
    b.resize(a.size());
    map<int,int> rep_chck;
    long mapsize=0;
    for (int i = (int)a.size()-1; i >= 0; i--)
    {
        int count = 0;
        for (int j = 0; j < a.size(); j++) if (a[j] < a[i]) count++;
        rep_chck.insert(pair<int,int>(count,i));
        while(rep_chck.size()==mapsize)
        {
            count++;
            rep_chck.insert(pair<int,int>(count,i));
        }
        mapsize=rep_chck.size();
        b[i] = count;       
    }
}
void CommFunc::getUnique(vector<uint32_t> &a)
{
    sort(a.begin(),a.end());
    vector<uint32_t> ::iterator it=unique(a.begin(),a.end());
    a.erase(it,a.end());
}
void CommFunc::match(const vector<uint32_t> &VecA, const vector<uint32_t> &VecB, vector<int> &VecC)
{
    int i=0;
    map<uint32_t, int> id_map;
    map<uint32_t, int>::iterator iter;
    VecC.clear();
    for(i=0; i<VecB.size(); i++) id_map.insert(pair<uint32_t,int>(VecB[i], i));
    for(i=0; i<VecA.size(); i++){
        iter=id_map.find(VecA[i]);
        if(iter==id_map.end()) VecC.push_back(-9);
        else VecC.push_back(iter->second);
    }    
}
