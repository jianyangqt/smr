//
//  SMR_data_p1.h
//  SMR_CPP
//
//  Created by Futao Zhang on 10/03/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __SMR_CPP__SMR_data_p1__
#define __SMR_CPP__SMR_data_p1__

#include "SMR_data.h"

namespace SMRDATA
{
    typedef struct{
        long lineNum;
        vector<string> _Expo_id;
        vector<int> _Expo_chr;
        vector<string> _Expo_gene;
        vector<int> _Expo_bp;
        vector<string> _Outco_id;
        vector<int> _Outco_chr;
        vector<string> _Outco_gene;
        vector<int> _Outco_bp;
        vector<string> _snp_rs;
        vector<int> _snp_chr;
        vector<int> _snp_bp;
        vector<char> _snp_a1;
        vector<char> _snp_a2;
        vector<float> _snp_frq;
        vector<double> _b;
        vector<double> _se;
        vector<double> _p_smr;
        vector<double> _p_heidi;
         vector<int> _nsnp;
        vector<uint32_t> _include;
    } eSMRrlt;
    

    
    void est_effect_splsize(char* eqtlsmaslstName, char* eqtlFileName, char* snplstName,char* problstName,char* snplst2exclde, char* problst2exclde,float thres);
    void iternal_test(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,int cis_itvl,char* smrFileName);
    void make_cojo(char* outFileName,char* eqtlFileName, char* snplstName,char* snplst2exclde, char* problstName, char* problst2exclde, char* genelistName, bool bFlag);
}

#endif /* defined(__SMR_CPP__SMR_data_p1__) */
