//
//  eData.cpp
//  SMR_CPP
//
//  Created by Futao Zhang on 7/08/2015.
//  Copyright (c) 2015 Futao Zhang. All rights reserved.
//

#include "eData.h"
// main function declared in gcta.h
void read_efile(eData* edata, string efile)
{
    int numCol2Skip=2;
    int lineNum=0;
    FILE* eflptr;
    if(!fopen_checked(&eflptr,efile.c_str(),"r"))
    {
        fprintf (stderr, "%s: Couldn't open file %s\n",
                 efile.c_str(), strerror (errno));
        exit (EXIT_FAILURE);
    }
    cout << "Reading gene expression / methylation data from [" + efile + "] ..." << endl;
    
    char buf[MAX_LINE_SIZE];
    fgets(buf, MAX_LINE_SIZE, eflptr);//head
    edata->prbNum=split_string_skip(buf, edata->_epi_prbid, "\x20\t\n",numCol2Skip);
    for(int i=0;i<edata->prbNum;i++) edata->_include.push_back(i);
   
    while (fgets(buf, MAX_LINE_SIZE, eflptr)) lineNum++;
    if(buf[0]=='\0') lineNum--;
    edata->indiNum=lineNum;
    rewind(eflptr);
    fgets(buf, MAX_LINE_SIZE, eflptr); //head
    edata->_eprVal.resize(edata->prbNum*lineNum);
    edata->_keep.resize(lineNum);
    edata->_eii_fid.resize(lineNum);
    edata->_eii_iid.resize(lineNum);
    vector<string> tmp;
    for(int i=0;i<lineNum;i++)
    {
        tmp.clear();
        fgets(buf, MAX_LINE_SIZE, eflptr);
        split_string(buf, tmp, "\x20\t\n");
        edata->_keep[i]=i;
        edata->_eii_fid[i]=tmp[0];
        edata->_eii_iid[i]=tmp[1];
        for(int j=2;j<tmp.size();j++)
            if(tmp[j]=="NA") edata->_eprVal[(j-2)*lineNum+i]=1e10;
            else edata->_eprVal[(j-2)*lineNum+i]=atof(tmp[j].c_str());
        
    }
    fclose(eflptr);
    cout<<"Expression data for "<<edata->prbNum<<" probes of "<<edata->indiNum<<" individuals have been included from the file [" + efile + "]."<<endl;
    
    }
