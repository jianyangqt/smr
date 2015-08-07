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
    /*
    fseek(eflptr,0, SEEK_END);
    int fszie=ftell(eflptr);
    rewind(eflptr);
     */
    char buf[MAX_LINE_SIZE];
    fgets(buf, MAX_LINE_SIZE, eflptr);
    edata->prbNum=split_string_skip(buf, edata->PrbID, "\x20\t\n",numCol2Skip);
    for(int i=0;i<edata->prbNum;i++) edata->include.push_back(i);
    vector<string> tmp;
    while(fgets(buf, MAX_LINE_SIZE, eflptr))
    {
        tmp.clear();
        split_string(buf, tmp, "\x20\t\n");
        edata->keep.push_back(lineNum++);
        edata->FID.push_back(tmp[0]);
        edata->IID.push_back(tmp[1]);
        for(int i=2;i<tmp.size();i++)
            if(tmp[i]=="NA") edata->eprVal.push_back(1e10);
            else edata->eprVal.push_back(atof(tmp[i].c_str()));
        
    }
    edata->indiNum=lineNum;
    if(edata->eprVal.size()!=lineNum*edata->prbNum)
    {
        cout << "Error: Unmatched data from [" + efile + "]" << endl;
        exit(EXIT_FAILURE);
    }
    fclose(eflptr);
    cout<<"Expression data for "<<edata->prbNum<<" probes of "<<edata->indiNum<<" individuals have been included from the file [" + efile + "]."<<endl;
    
    }
