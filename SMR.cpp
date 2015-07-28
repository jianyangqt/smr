//
//  SMR.cpp
//  SRM_CPP
//
//  Created by Futao Zhang on 29/06/15.
//  Copyright (c) 2015 Futao Zhang. All rights reserved.
//

#include "SMR.h"

#include "StrFunc.h"
#include "SMR_data.h"


using namespace SMRDATA;
using namespace StatFunc;

void option(int option_num, char* option_str[]);
int main(int argc, char** argv) {    
    
    cout << "*******************************************************************" << endl;
    cout << "* Summary-data-based Mendelian Randomization (SMR)" << endl;
    cout << "* version 0.3" << endl;
    cout << "* (C) 2015 Futao Zhang, Zhihong Zhu and Jian Yang" << endl;
    cout << "* The University of Queensland" << endl;
    cout << "* MIT License" << endl;
    cout << "*******************************************************************" << endl;
    
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
    bool bFlag=false;// for binary file
    double maf=0.0;
    
    // for data management
    bool make_besd_flag=false;
    bool make_esd_flag=false;
    char* problstName=NULL;
    
	
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
        else if(strcmp(option_str[i],"--maf")==0){
            maf=atof(option_str[++i]);
            cout<<"--maf "<<maf<<endl;
            if(maf<0 || maf>0.5) throw("\nError: --maf should be within the range from 0 to 0.5.\n");
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
    }
    cout<<endl;
    
    if(make_besd_flag || make_esd_flag) make_esd_file(outFileName, bFileName,gwasFileName, eqtlFileName, maf,indilstName, snplstName,problstName,bFlag,make_besd_flag,make_esd_flag);
    else   smr(outFileName, bFileName,gwasFileName, eqtlFileName, maf,indilstName, snplstName,problstName,bFlag);
    
   }
