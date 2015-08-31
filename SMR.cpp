//
//  SMR.cpp
//  SRM_CPP
//
//  Created by Futao Zhang on 29/06/15.
//  Copyright (c) 2015 Futao Zhang. All rights reserved.
//


#include "StrFunc.h"
#include "SMR_data.h"
#include "eData.h"

using namespace SMRDATA;
using namespace StatFunc;

void option(int option_num, char* option_str[]);
static inline bool not_in_flags(vector<string> &flags, string str)
{
    return find(flags.begin(),flags.end(),str) == flags.end();
}
static inline void FLAG_VALID_CK(int option_num, char* option_str[])
{
    const char *flgs[] = { "--bfile","--gwas-summary","--beqtl-summary","--maf","--keep","--remove","--extract-snp","--exclude-snp","--extract-probe",
        "--exclude-probe","--eqtl-summary","--ld-pruning","--peqtl-hetero","--m-hetero","--make-besd","--make-esd","--out", "--peqtl-smr","--smr",
    "--cis-itvl","--trans-thres","--rest-thres","--efile"};
    
    vector<string> flags(flgs, flgs + sizeof(flgs)/sizeof(flgs[0]));
    
    if(option_num<3)
    {
        cout<<"flags include:"<<endl;
        int cur_mark=0;
        for(int i=0;i<flags.size();i++)
        {
            int tmp=i>>2;
            if(tmp>cur_mark)
            {
                cout<<endl;
                cur_mark=tmp;
            }
            cout<<flags[i]<<",";
        }
        exit (EXIT_FAILURE);
    }
    for(int i=0;i<option_num;i++)
    {
        if(has_prefix(option_str[i],"--"))
            if(not_in_flags(flags, option_str[i]))
            {
                
                fprintf (stderr, "%s: Invalid option\n",
                         option_str[i]);
                exit (EXIT_FAILURE);
            }
    }

}
int main(int argc, char** argv) {
    
    cout << "*******************************************************************" << endl;
    cout << "* Summary-data-based Mendelian Randomization (SMR)" << endl;
    cout << "* version 0.3" << endl;
    cout << "* (C) 2015 Futao Zhang, Zhihong Zhu and Jian Yang" << endl;
    cout << "* The University of Queensland" << endl;
    cout << "* MIT License" << endl;
    cout << "*******************************************************************" << endl;
    
    FLAG_VALID_CK(argc,argv);
    long int time_used = 0, start = time(NULL);
    time_t curr = time(0);
    cout << "Analysis started: " << ctime(&curr) << endl;
    cout << "Options:" << endl;
    try {
        option(argc, argv);
    } catch (const string &err_msg) {
        cerr << "\n" << err_msg << endl;
    } catch (const char *err_msg) {
        cerr << "\n" << err_msg << endl;
    }    
    curr = time(0);
    cout << "\nAnalysis finished: " << ctime(&curr);
    time_used = time(NULL) - start;
    cout << "Computational time: " << time_used / 3600 << ":" << (time_used % 3600) / 60 << ":" << time_used % 60 << endl;
    
    return 0;
}

void option(int option_num, char* option_str[])
{    
    char* bFileName=NULL;    
    char* gwasFileName=NULL;
    char* eqtlFileName=NULL;
    char* outFileName=NULL;
    char* indilstName=NULL;
    char* snplstName=NULL;
    char* indilst2remove=NULL;
    char* snplst2exclde=NULL;
    char* problst2exclde=NULL;
    bool bFlag=false;// for binary file
    double maf=0.0;
    double p_hetero=1.5654e-3;
    double p_smr=5.0e-8;
    double ld_prune=0.9;
    unsigned int m_hetero=10;
    
    // for data management
    bool make_besd_flag=false;
    bool make_esd_flag=false;
    char* problstName=NULL;
    bool smr_flag=false;
	
    char* eFileName=NULL;
    bool eremlFlag=false;
    
    // for filtering
    bool cis_flag=false;
    int cis=1;
    float transThres=5.0e-8;
    float restThres=1.0e-5;
    
    for(int i=0;i<option_num;i++)
    {
        // Plink files for LD from a reference sample
        if(0==strcmp(option_str[i],"--bfile")){
            bFileName=option_str[++i];
            printf("--bfile %s\n",bFileName);
        }
        // gwas data file as cojo format
        else if(0==strcmp(option_str[i],"--gwas-summary")){
            gwasFileName=option_str[++i];
            printf("--gwas-summary %s\n",gwasFileName);
            CommFunc::FileExist(gwasFileName);
        }
        // eQTL files
        else if(0==strcmp(option_str[i],"--eqtl-summary")){
            eqtlFileName=option_str[++i];
            printf("--eqtl-summary %s\n",eqtlFileName);
            
        }
        else if(0==strcmp(option_str[i],"--beqtl-summary")){
            eqtlFileName=option_str[++i];
            bFlag=true;
            printf("--beqtl-summary %s\n",eqtlFileName);
            
        }
        else if(strcmp(option_str[i],"--keep")==0){
            indilstName=option_str[++i];
            cout<<"--keep "<<indilstName<<endl;
            CommFunc::FileExist(indilstName);
        }
        else if(strcmp(option_str[i],"--remove")==0){
            indilst2remove=option_str[++i];
            cout<<"--remove "<<indilst2remove<<endl;
            CommFunc::FileExist(indilst2remove);
        }
        else if(strcmp(option_str[i],"--extract-snp")==0){
            snplstName=option_str[++i];
            cout<<"--extract-snp "<<snplstName<<endl;
            CommFunc::FileExist(snplstName);
        }
        else if(strcmp(option_str[i],"--extract-probe")==0){
            problstName=option_str[++i];
            cout<<"--extract-probe "<<problstName<<endl;
            CommFunc::FileExist(problstName);
        }
        else if(strcmp(option_str[i],"--exclude-snp")==0){
            snplst2exclde=option_str[++i];
            cout<<"--exclude-snp "<<snplst2exclde<<endl;
            CommFunc::FileExist(snplst2exclde);
        }
        else if(strcmp(option_str[i],"--exclude-probe")==0){
            problst2exclde=option_str[++i];
            cout<<"--exclude-probe "<<problst2exclde<<endl;
            CommFunc::FileExist(problst2exclde);
        }
        else if(strcmp(option_str[i],"--maf")==0){
            maf=atof(option_str[++i]);
            cout<<"--maf "<<maf<<endl;
            if(maf<0 || maf>0.5)
            {
                fprintf (stderr, "Error: --maf should be within the range from 0 to 0.5.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--make-besd")==0){
            make_besd_flag=true;
            cout<<"--make-besd "<<endl;
        }
        else if(strcmp(option_str[i],"--make-esd")==0){
            make_esd_flag=true;
            cout<<"--make-esd "<<endl;
        }
        else if (0 == strcmp(option_str[i], "--out")){
            outFileName = option_str[++i];
            printf("--out %s\n", outFileName);
        }
        else if (0 == strcmp(option_str[i], "--peqtl-smr")){
            p_smr = atof(option_str[++i]);
            printf("--peqtl-smr %10.2e\n", p_smr);
            if(p_smr<0 || p_smr>1)
            {
                fprintf (stderr, "Error: --peqtl-smr should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--peqtl-hetero")){
            p_hetero = atof(option_str[++i]);
            printf("--peqtl-hetero %10.2e\n", p_hetero);
            if(p_hetero<0 || p_hetero>1)
            {
                fprintf (stderr, "Error: --peqtl-hetero should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
            
        }
        else if (0 == strcmp(option_str[i], "--m-hetero")){
            m_hetero = atoi(option_str[++i]);
            printf("--m-hetero %d\n", m_hetero);
        }
        else if (0 == strcmp(option_str[i], "--ld-pruning")){
            ld_prune = atof(option_str[++i]);
            printf("--ld-pruning %f\n", ld_prune);
            if(ld_prune<0 || ld_prune>1)
            {
                fprintf (stderr, "Error: --ld-pruning should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--smr")){
            smr_flag=true;
            printf("--smr \n");
        }
        else if(strcmp(option_str[i],"--cis-itvl")==0){
            cis_flag=true;
            cis=atoi(option_str[++i]);
            printf("--cis-itvl %d Mb\n", cis);
            if(cis<0 )
            {
                fprintf (stderr, "Error: --cis-itvl should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--trans-thres")==0){
            transThres=atof(option_str[++i]);
            printf("--trans-thres %10.2e\n", transThres);
            if(transThres<0 || transThres>1)
            {
                fprintf (stderr, "Error: --trans-thres should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--rest-thres")==0){
            restThres=atof(option_str[++i]);
            printf("--rest-thres %10.2e\n", restThres);
            if(restThres<0 || restThres>1)
            {
                fprintf (stderr, "Error: --rest-thres should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--efile")){
            eremlFlag=true;
            eFileName = option_str[++i];
            printf("--efile %s\n", eFileName);
        }


    }
    cout<<endl;
    eData edata;
    char tmpch[4]="smr";
    if(outFileName == NULL) outFileName=tmpch;
    if(make_besd_flag || make_esd_flag || cis_flag) make_esd_file(outFileName, bFileName,gwasFileName, eqtlFileName, maf,indilstName, snplstName,problstName,bFlag,make_besd_flag,make_esd_flag, indilst2remove, snplst2exclde, problst2exclde, cis_flag, cis, transThres, restThres);
    else if(smr_flag)  smr(outFileName, bFileName,gwasFileName, eqtlFileName, maf,indilstName, snplstName,problstName,bFlag,p_hetero,ld_prune,m_hetero, indilst2remove, snplst2exclde, problst2exclde,p_smr);
    else if(eremlFlag) read_efile(&edata, eFileName);
    
   }
