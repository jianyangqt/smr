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
#include "SMR.h"

using namespace SMRDATA;
using namespace StatFunc;
using namespace E_DATA;
int thread_num;

int main(int argc, char** argv) {
        
    cout << "*******************************************************************" << endl;
    cout << "* Summary-data-based Mendelian Randomization (SMR)" << endl;
    cout << "* version 0.6" << endl;
    cout << "* (C) 2015 Futao Zhang, Zhihong Zhu and Jian Yang" << endl;
    cout << "* The University of Queensland" << endl;
    cout << "* MIT License" << endl;
    cout << "*******************************************************************" << endl;
   
    FLAGS_VALID_CK(argc,argv);
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
    char* eproblst2exclde=NULL;
     char* mproblst2exclde=NULL;
    bool bFlag=false;// for binary file
    double maf=0.0;
    double p_hetero=1.5654e-3;
    double p_smr=5.0e-8;
    double ld_prune=0.9;
    unsigned int m_hetero=10;
    bool smr_flag=true;
    bool smr_trans_flag=false;
    
    // for data management
    bool make_besd_flag=false;
    bool make_esd_flag=false;
    char* problstName=NULL;
    char* eproblstName=NULL;
    char* mproblstName=NULL;
	
    char* eFileName=NULL;
    bool eremlFlag=false;
    
    // for filtering
    bool cis_flag=false;
    int cis_itvl=2000;
    int trans_itvl=1000;
    float transThres=5.0e-8;
    float restThres=1.0e-5;
    
    // for lookup
    float plookup=5e-8;
    bool lookup_flag=false;
    char* genelistName=NULL;
    
    
    char* refSNP=NULL;
    bool heidioffFlag = false;
    
     thread_num = 1;
    
    bool combineFlg = false;
    char* eqtlsmaslstName = NULL;
    //
    char* gwasFileName2=NULL;
    char* eqtlFileName2=NULL;
    char* traitlstName=NULL;
    //
    bool plotflg=false;
    //
    char* syllabusName=NULL;
    bool gctaflag=false;
    bool plinkflag=false;
    bool gemmaflag=false;
    bool merlinflag=false;
    
    bool tosbesdflag=false;
    
    char* freqName=NULL;
    bool esdstd=false;
    
    bool metaflg=false;
    
    
    
    for(int i=0;i<option_num;i++)
    {
        // Plink files for LD from a reference sample
        if(0==strcmp(option_str[i],"--bfile")){
            bFileName=option_str[++i];
            FLAG_VALID_CK("--bfile", bFileName);
            printf("--bfile %s\n",bFileName);
        }
        // gwas data file as cojo format
        else if(0==strcmp(option_str[i],"--gwas-summary")){
            if(gwasFileName ==NULL){
                gwasFileName=option_str[++i];
                FLAG_VALID_CK("--gwas-summary", gwasFileName);
                printf("--gwas-summary %s\n",gwasFileName);
                CommFunc::FileExist(gwasFileName);
            }else
            {
                gwasFileName2=option_str[++i];
                FLAG_VALID_CK("--gwas-summary", gwasFileName2);
                printf("--gwas-summary %s\n",gwasFileName2);
                CommFunc::FileExist(gwasFileName2);
            }
        }
        // eQTL files
        else if(0==strcmp(option_str[i],"--eqtl-summary")){
            if(eqtlFileName==NULL)
            {
                eqtlFileName=option_str[++i];
                FLAG_VALID_CK("--eqtl-summary", eqtlFileName);
                printf("--eqtl-summary %s\n",eqtlFileName);
            }else{
                eqtlFileName2=option_str[++i];
                FLAG_VALID_CK("--eqtl-summary", eqtlFileName2);
                printf("--eqtl-summary %s\n",eqtlFileName2);
            }
           
            
        }
        else if(0==strcmp(option_str[i],"--beqtl-summary")){
             bFlag=true;
            if(eqtlFileName==NULL)
            {
                eqtlFileName=option_str[++i];
                FLAG_VALID_CK("--beqtl-summary", eqtlFileName);
                printf("--beqtl-summary %s\n",eqtlFileName);
            }else{
                eqtlFileName2=option_str[++i];
                FLAG_VALID_CK("--beqtl-summary", eqtlFileName2);
                printf("--beqtl-summary %s\n",eqtlFileName2);
            }
        }
        else if(strcmp(option_str[i],"--keep")==0){
            indilstName=option_str[++i];
            FLAG_VALID_CK("--keep", indilstName);
            cout<<"--keep "<<indilstName<<endl;
            CommFunc::FileExist(indilstName);
        }
        else if(strcmp(option_str[i],"--remove")==0){
            indilst2remove=option_str[++i];
            FLAG_VALID_CK("--remove", indilst2remove);
            cout<<"--remove "<<indilst2remove<<endl;
            CommFunc::FileExist(indilst2remove);
        }
        else if(strcmp(option_str[i],"--extract-snp")==0){
            snplstName=option_str[++i];
            FLAG_VALID_CK("--extract-snp", snplstName);
            cout<<"--extract-snp "<<snplstName<<endl;
            CommFunc::FileExist(snplstName);
        }
        else if(strcmp(option_str[i],"--extract-probe")==0){
            problstName=option_str[++i];
            FLAG_VALID_CK("--extract-probe", problstName);
            cout<<"--extract-probe "<<problstName<<endl;
            CommFunc::FileExist(problstName);
        }
        else if(strcmp(option_str[i],"--extract-eprobe")==0){
            eproblstName=option_str[++i];
            FLAG_VALID_CK("--extract-eprobe", eproblstName);
            cout<<"--extract-eprobe "<<eproblstName<<endl;
            CommFunc::FileExist(eproblstName);
        }
        else if(strcmp(option_str[i],"--extract-mprobe")==0){
            mproblstName=option_str[++i];
            FLAG_VALID_CK("--extract-mprobe", mproblstName);
            cout<<"--extract-mprobe "<<mproblstName<<endl;
            CommFunc::FileExist(mproblstName);
        }
        else if(strcmp(option_str[i],"--exclude-snp")==0){
            snplst2exclde=option_str[++i];
            FLAG_VALID_CK("--exclude-snp", snplst2exclde);
            cout<<"--exclude-snp "<<snplst2exclde<<endl;
            CommFunc::FileExist(snplst2exclde);
        }
        else if(strcmp(option_str[i],"--exclude-probe")==0){
            problst2exclde=option_str[++i];
            FLAG_VALID_CK("--exclude-probe", problst2exclde);
            cout<<"--exclude-probe "<<problst2exclde<<endl;
            CommFunc::FileExist(problst2exclde);
        }
        else if(strcmp(option_str[i],"--exclude-eprobe")==0){
            eproblst2exclde=option_str[++i];
            FLAG_VALID_CK("--exclude-eprobe", eproblst2exclde);
            cout<<"--exclude-eprobe "<<eproblst2exclde<<endl;
            CommFunc::FileExist(eproblst2exclde);
        }
        else if(strcmp(option_str[i],"--exclude-mprobe")==0){
            mproblst2exclde=option_str[++i];
            FLAG_VALID_CK("--exclude-mprobe", mproblst2exclde);
            cout<<"--exclude-mprobe "<<mproblst2exclde<<endl;
            CommFunc::FileExist(mproblst2exclde);
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
            if(outFileName !=NULL && has_prefix(outFileName, "--"))
            {
                outFileName=NULL;
                i--;
            }else   printf("--out %s\n", outFileName);
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
        else if (0 == strcmp(option_str[i], "--peqtl-heidi")){
            p_hetero = atof(option_str[++i]);
            printf("--peqtl-heidi %10.2e\n", p_hetero);
            if(p_hetero<0 || p_hetero>1)
            {
                fprintf (stderr, "Error: --peqtl-heidi should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
            
        }
        else if (0 == strcmp(option_str[i], "--heidi-m")){
            m_hetero = atoi(option_str[++i]);
            printf("--heidi-m %d\n", m_hetero);
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
        else if(strcmp(option_str[i],"--cis-wind")==0){
            cis_flag=true;
            cis_itvl=atoi(option_str[++i]);
            printf("--cis-wind %d Kb\n", cis_itvl);
            if(cis_itvl<0 )
            {
                fprintf (stderr, "Error: --cis-wind should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--trans-wind")==0){
            cis_flag=true;
            trans_itvl=atoi(option_str[++i]);
            printf("--trans-wind %d Kb\n", trans_itvl);
            if(trans_itvl<0 )
            {
                fprintf (stderr, "Error: --trans-wind should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--peqtl-trans")==0){
            transThres=atof(option_str[++i]);
            printf("--peqtl-trans %10.2e\n", transThres);
            if(transThres<0 || transThres>1)
            {
                fprintf (stderr, "Error: --peqtl-trans should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--peqtl-other")==0){
            restThres=atof(option_str[++i]);
            printf("--peqtl-other %10.2e\n", restThres);
            if(restThres<0 || restThres>1)
            {
                fprintf (stderr, "Error: --peqtl-other should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if (0 == strcmp(option_str[i], "--efile")){
            eremlFlag=true;
            eFileName = option_str[++i];
            FLAG_VALID_CK("--efile", eFileName);
            printf("--efile %s\n", eFileName);
        }
        else if (0 == strcmp(option_str[i], "--lookup")){
            lookup_flag=true;
            if(i+1==option_num || has_prefix(option_str[i+1],"--"))  plookup = 5.0e-8;
            else plookup = atof (option_str[++i]);
            printf("--lookup %10.2e\n", plookup);
            if(plookup<0 || plookup>1)
            {
                fprintf (stderr, "Error: --lookup should be within the range from 0 to 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        else if(strcmp(option_str[i],"--gene-list")==0){
            genelistName=option_str[++i];
            FLAG_VALID_CK("--gene-list", genelistName);
            cout<<"--gene-list "<<genelistName<<endl;
            CommFunc::FileExist(genelistName);
        }
        else if (0 == strcmp(option_str[i], "--target-snp")){
            refSNP = option_str[++i];
            FLAG_VALID_CK("--target-snp", refSNP);
            printf("--target-snp %s\n", refSNP);
        }
        else if (0 == strcmp(option_str[i], "--heidi-off")){
            heidioffFlag = true;            
            printf("--heidi-off \n");
        }
        else if(0 == strcmp(option_str[i],"--thread-num")){
            thread_num=atoi(option_str[++i]);
            printf("--thread-num %d\n", thread_num);
        }
        else if (0 == strcmp(option_str[i], "--combine-besd")){
            combineFlg = true;
            printf("--combine-besd \n");
        }
        else if (0 == strcmp(option_str[i], "--beqtl-summaries")){
            eqtlsmaslstName = option_str[++i];
            FLAG_VALID_CK("--beqtl-summaries", eqtlsmaslstName);
            printf("--beqtl-summaries %s\n", eqtlsmaslstName);
        }
        else if (0 == strcmp(option_str[i], "--plot")){
            plotflg = true;
            printf("--plot \n" );
        }
        else if (0 == strcmp(option_str[i], "--trans")){
            smr_trans_flag = true;
            printf("--trans \n" );
        }
        else if (0 == strcmp(option_str[i], "--eqtl-flist")){
            syllabusName = option_str[++i];
             gctaflag = true;
            if(syllabusName !=NULL && has_prefix(syllabusName, "--"))
            {
                syllabusName=NULL;
                i--;
            }else   printf("--eqtl-flist %s\n", syllabusName);
        }
        else if (0 == strcmp(option_str[i], "--smr-format")){
            gctaflag = true;
            printf("--smr-format \n" );
        }
        else if (0 == strcmp(option_str[i], "--plink-qassoc-format")){
            plinkflag = true;
             gctaflag = false;
            printf("--plink-qassoc-format \n" );
        }
        else if (0 == strcmp(option_str[i], "--gemma-format")){
            gemmaflag = true;
            gctaflag = false;
            printf("--gemma-format \n" );
        }
        else if (0 == strcmp(option_str[i], "--merlin-fastassoc-format")){
            merlinflag = true;
            gctaflag = false;
            printf("----merlin-fastassoc-format \n" );
        }
        else if (0 == strcmp(option_str[i], "--make-sbesd")){
            tosbesdflag = true;
            printf("--make-sbesd \n" );
        }
        else if(strcmp(option_str[i],"--extract-trait")==0){
            traitlstName=option_str[++i];
            FLAG_VALID_CK("--extract-trait", traitlstName);
            cout<<"--extract-trait "<<traitlstName<<endl;
            CommFunc::FileExist(traitlstName);
        }
        else if(strcmp(option_str[i],"--esd-std")==0){
            esdstd=true;
            cout<<"--esd-std "<<endl;
        }
        else if(strcmp(option_str[i],"--freq")==0){
            freqName=option_str[++i];
            FLAG_VALID_CK("--freq", freqName);
            cout<<"--freq "<<freqName<<endl;
            CommFunc::FileExist(freqName);
        }
        else if(strcmp(option_str[i],"--meta")==0){
            metaflg=true;
            cout<<"--meta "<<endl;
        }

    }
    
#ifndef __APPLE__
    stringstream ss;
    ss << thread_num;
    setenv("OMP_NUM_THREADS", ss.str().c_str(), 1);
    omp_set_num_threads(thread_num);
#endif
    
    cout<<endl;
    eData edata;
    char tmpch[4]="smr";
    if(outFileName == NULL) outFileName=tmpch;
    if(make_besd_flag && (gctaflag || plinkflag || gemmaflag || merlinflag) ) make_besd(outFileName, syllabusName, gctaflag, plinkflag, gemmaflag,merlinflag); // from text to besd
    else if (esdstd) standardization(outFileName, eqtlFileName,bFlag,freqName);
    else if (metaflg) meta(outFileName,eqtlFileName, eqtlFileName2);
    else if (combineFlg) combineCis(eqtlsmaslstName, outFileName, cis_flag, cis_itvl,trans_itvl, transThres, restThres);
    else if(tosbesdflag)  esd2sbesd(outFileName, eqtlFileName );
    else if(make_besd_flag || make_esd_flag ) make_esd_file(outFileName, bFileName,gwasFileName, eqtlFileName, maf,indilstName, snplstName,problstName,bFlag,make_besd_flag,make_esd_flag, indilst2remove, snplst2exclde, problst2exclde,  cis_flag, cis_itvl,trans_itvl, transThres, restThres);
    else if(gwasFileName2 != NULL) smr_g2g(gwasFileName,gwasFileName2,snplstName,snplst2exclde); // gwas summary by gwas summary , not finished.
    else if(eqtlFileName2 != NULL) smr_e2e(outFileName, bFileName,eqtlFileName, eqtlFileName2, maf,indilstName, snplstName,eproblstName,mproblstName,bFlag,p_hetero,ld_prune,m_hetero, indilst2remove, snplst2exclde, eproblst2exclde,mproblst2exclde,p_smr,refSNP, heidioffFlag,cis_itvl,traitlstName,plotflg);
    else if(eremlFlag) read_efile(&edata, eFileName);
    else if(lookup_flag) lookup(outFileName,eqtlFileName, snplstName, problstName, genelistName, plookup, bFlag);
    else if(smr_flag && !smr_trans_flag) smr(outFileName, bFileName,gwasFileName, eqtlFileName, maf,indilstName, snplstName,problstName,bFlag,p_hetero,ld_prune,m_hetero, indilst2remove, snplst2exclde, problst2exclde,p_smr,refSNP, heidioffFlag,cis_itvl, plotflg);
    else if(smr_flag && smr_trans_flag) smr_trans_wholeInOne(outFileName, bFileName,gwasFileName, eqtlFileName, maf,indilstName, snplstName,problstName, bFlag,p_hetero,ld_prune,m_hetero, indilst2remove, snplst2exclde, problst2exclde,transThres,refSNP, heidioffFlag,trans_itvl, plotflg);
    
   }
