//
//  SMR_data_p2.h
//  SMR_CPP
//
//  Created by Futao Zhang on 10/06/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#ifndef __SMR_CPP__SMR_data_p2__
#define __SMR_CPP__SMR_data_p2__

#include "SMR_data.h"

namespace SMRDATA
{
    void make_besd(char*outFileName, char* syllabusName, bool gctaflag,bool plinkflag,bool gemmaflag,bool merlinflag,bool save_dense_flag,int cis_itvl, int trans_itvl, float transThres, float restThres);
    void make_besd_byQfile(char* qfile, char* outFileName,bool save_dense_flag,int cis_itvl, int trans_itvl, float transThres, float restThres);
}
#endif /* defined(__SMR_CPP__SMR_data_p2__) */
