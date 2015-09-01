//
//  eData.h
//  SMR_CPP
//
//  Created by Futao Zhang on 7/08/2015.
//  Copyright (c) 2015 Futao Zhang. All rights reserved.
//

#ifndef __SMR_CPP__eData__
#define __SMR_CPP__eData__


#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "CommFunc.h"

using namespace CommFunc;
using namespace StrFunc;

typedef struct{
    
    uint32_t prbNum;
    uint32_t indiNum;
    
    //eii
    vector< string > _eii_fid;
    vector< string > _eii_iid;
    vector<bool> _eii_fain;
    vector<bool> _eii_moin;
    vector<char> _eii_sex;
    vector<float> _eii_pheno;
    vector<uint32_t> _keep;
    
    //epi
    vector<uint32_t> _epi_chr;
    vector<string> _epi_prbid;
    vector<uint32_t> _epi_gd;
    vector<uint32_t> _epi_bp;
    vector<string> _epi_gene;
    vector<char> _epi_orien;
    vector<uint32_t> _include;
    
    //eed
    vector<double> _eprVal;// probe-major
   
} eData;

void read_efile(eData* edata, string efile);

#endif /* defined(__SMR_CPP__eData__) */

