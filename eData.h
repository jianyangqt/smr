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
    
    vector< string > FID;
    vector< string > IID;
    vector< string > PrbID;
  
    vector<uint32_t> include;
    vector<uint32_t> keep;
    
    vector<double> eprVal;
   
} eData;

void read_efile(eData* edata, string efile);

#endif /* defined(__SMR_CPP__eData__) */

