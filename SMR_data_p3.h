//
//  SMR_data_p3.h
//  SMR_CPP
//
//  Created by Futao Zhang on 15/06/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __SMR_CPP__SMR_data_p3__
#define __SMR_CPP__SMR_data_p3__

#include "SMR_data.h"

namespace SMRDATA
{
    void combineBesd(char* eqtlsmaslstName, char* outFileName,bool makecisflag, int cis_itvl, int trans_itvl, float transThres, float restThres);
    void make_sparse_besd(char* eqtlFileName, char* outFileName, int cis_itvl, int trans_itvl, float transThres, float restThres);
}

#endif /* defined(__SMR_CPP__SMR_data_p3__) */
