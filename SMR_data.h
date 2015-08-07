//
//  SMR_data.h
//  SRM_CPP
//
//  Created by Futao Zhang on 29/06/15.
//  Copyright (c) 2015 Futao Zhang. All rights reserved.
//

#ifndef __SRM_CPP__SMR_data__
#define __SRM_CPP__SMR_data__


#include "CommFunc.h"
#include "StatFunc.h"
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <errno.h>
#include <zlib.h>
#include <bitset>

//#include <omp.h>

typedef unsigned int         uint32_t;
typedef MatrixXf eigenMatrix;
typedef VectorXf eigenVector;

using namespace std;
using namespace StatFunc;
using namespace CommFunc;

namespace SMRDATA
{
    //implementation of SoA(Struct of Array)
    typedef struct{
        // bim file
        int _autosome_num;
        vector<int> _chr;
        vector<string> _snp_name;
        map<string, int> _snp_name_map;
        vector<double> _genet_dst;
        vector<int> _bp;
        vector<char> _allele1;
        vector<char> _allele2;
        vector<char> _ref_A; // reference allele
        vector<char> _other_A; // the other allele
        int _snp_num;
        vector<double> _rc_rate;
        vector<int> _include; // initialized in the read_bimfile()
        eigenVector _maf;
        
        // fam file
        vector<string> _fid;
        vector<string> _pid;
        map<string, int> _id_map;
        vector<string> _fa_id;
        vector<string> _mo_id;
        vector<int> _sex;
        vector<double> _pheno;
        int _indi_num;
        vector<int> _keep; // initialized in the read_famfile()
        eigenMatrix _varcmp_Py; // BLUP solution to the total genetic effects of individuals
        
        // bed file
        vector< vector<bool> > _snp_1;
        vector< vector<bool> > _snp_2;
        
        // imputed data
        bool _dosage_flag;
        vector< vector<float> > _geno_dose;
        vector<double> _impRsq;
        
        // genotypes
        MatrixXf _geno;
        
        
        vector<double> _mu;
    } bInfo;
    

    typedef struct{
        long snpNum;
        vector<string> snpName;
        char* allele_1;
        char* allele_2;
        double* freq;
        double* byz;
        double* seyz;
        double* pvalue;
        uint32_t* splSize;
    } gwasData;
    
    typedef struct{
        vector<uint32_t> _esi_chr;
        vector<string> _esi_rs;
        vector<uint32_t> _esi_gd;
        vector<uint32_t> _esi_bp;
        vector<char> _esi_allele1;
        vector<char> _esi_allele2;
		vector<int> _esi_include; // initialized in the readesi

        vector<uint32_t> _epi_chr;
        vector<string> _epi_prbID;
        vector<uint32_t> _epi_gd;
        vector<uint32_t> _epi_bp;
        vector<string> _epi_gene;
        vector<char> _epi_orien;
        vector<int> _include; // initialized in the readepi
        //for sparse
        vector<uint32_t> _cols;
        vector<uint32_t> _rowid;
        vector<double> _val;
        // for dense
        vector< vector<float> > _bxz; // first dimension is probe, second is snp
        vector< vector<float> > _sexz;
        
        uint32_t _probNum;
        uint32_t _snpNum;
        uint32_t _valNum;
        
    } eqtlInfo;
    
    
    void read_bimfile(bInfo* bdata,string bimfile);
    void read_famfile(bInfo* bdata,string famfile);
    void read_bedfile(bInfo* bdata,string bedfile);
    void read_gwas_data(gwasData* gdata, char* gwasFileName);
    void read_esifile(eqtlInfo* eqtlinfo, string esifile);
    void read_epifile(eqtlInfo* eqtlinfo, string epifile);
    void read_besdfile(eqtlInfo* eqtlinfo, string besdfile);
    void read_esdfile(eqtlInfo* eqtlinfo, string esdfile);
   
   
    
    bool has_suffix(const std::string &str, const std::string &suffix);
    void ld_calcualte(double* ref_ld,float* ref_snpData,int rsize,int csize);
    void get_square_idxes(vector<int> &sn_ids,vector<double> &zsxz,double threshold);
    void est_cov_bxy(double* covbxy, vector<double> zsxz1,vector<double> bxy1,vector<double> seyz1,vector<double> bxz1,double* ref_ld1);
    double bxy_hetero3(vector<double> byz1, vector<double> bxz1,vector<double> seyz1,vector<double> sexz1,vector<double> zsxz1,double* ref_ld1, long* nsnp);
    
    bool make_XMat(bInfo* bdata, MatrixXf &X);
    void make_XMat_SNPs(bInfo* bdata, vector< vector<float> > &X, bool miss_with_mu);
    void LD_Blocks(bInfo* bdata, double wind_size, double alpha, bool IncldQ,vector<string> &_ld_target_snp);
    void makex_eigenVector(bInfo* bdata,int j, eigenVector &x, bool resize, bool minus_2p);
    void calcu_mu(bInfo* bdata, bool ssq_flag=false);
    // inline functions
    template<typename ElemType>
    void makex(bInfo* bdata,int j, vector<ElemType> &x, bool minus_2p = false) {
        int i = 0;
        x.resize(bdata->_keep.size());
        for (i = 0; i < bdata->_keep.size(); i++) {
            if (!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]) {
                if (bdata->_allele1[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) x[i] = (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                else x[i] = 2.0 - (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
            } else x[i] = bdata->_mu[bdata->_include[j]];
            if (minus_2p) x[i] -= bdata->_mu[bdata->_include[j]];
        }
    }
    void read_msglist(string msglistfile, vector<string> &msglist, string msg);
    void makeptrx(bInfo* bdata,int bsnpid,int cursnpid, float* x, bool minus_2p = false);
    void keep_indi(bInfo* bdata,string indi_list_file);
    void extract_snp(bInfo* bdata,string snplistfile);
    void filter_snp_maf(bInfo* bdata,double maf);
    void extract_prob(eqtlInfo* eqtlinfo,string problstName);
	void extract_eqtl_snp(eqtlInfo* eqtlinfo, string snplstName);

    string dtos(double value);
    string dtosf(double value);
    string itos(int value);
    
    void free_gwas_data(gwasData* gdata);  
  
    int file_read_check(ifstream* in_file, const char* filename);
    
    void smr(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf, char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero , char* indilst2remove, char* snplst2exclde, char* problst2exclde);
    void make_esd_file(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf, char* indilstName, char* snplstName,char* problstName,bool bFlag,bool make_besd_flag,bool make_esd_flag, char* indilst2remove, char* snplst2exclde, char* problst2exclde);
}
#endif /* defined(__SRM_CPP__SMR_data__) */
