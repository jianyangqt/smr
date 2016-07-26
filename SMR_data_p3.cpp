//
//  SMR_data_p3.cpp
//  SMR_CPP
//
//  Created by Futao Zhang on 15/06/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "SMR_data_p3.h"

namespace SMRDATA
{
    void combine_esi(vector<snpinfolst> &snpinfo, vector<string> &smasNames, bool genouni)
    {
        long counter = 0;
        map<string, int> rs_map;
        map<string, int> rsbp_map;
        long f2r=smasNames.size();
        if(genouni) f2r=1;
     
        for (int i = 0; i < f2r; i++)
        {
            eqtlInfo etmp;
            string esifile = smasNames[i]+".esi";
            read_esifile(&etmp, esifile);
         
            for (int j = 0; j<etmp._snpNum; j++)
            {
                string crsbpstr=etmp._esi_rs[j]+":"+atos(etmp._esi_bp[j]);
                rs_map.insert(pair<string, int>(etmp._esi_rs[j].c_str(), counter));
                rsbp_map.insert(pair<string, int>(crsbpstr.c_str(), counter));
                if(rs_map.size() != rsbp_map.size())
                {
                    printf("ERROR: SNP %s has multiple BPs, please make sure that the genomic coordinates for the SNPs are on the same page.\n", etmp._esi_rs[j].c_str()) ;
                    exit(EXIT_FAILURE);
                }
                if (counter < rs_map.size())
                {
                    snpinfolst snpinfotmp;
                    counter = rs_map.size();
                    snpinfotmp.snpchr=etmp._esi_chr[j];
                    snpinfotmp.snprs=etmp._esi_rs[j];
                    snpinfotmp.bp=etmp._esi_bp[j];
                    snpinfotmp.gd=etmp._esi_gd[j];
                    snpinfotmp.a1=etmp._esi_allele1[j];
                    snpinfotmp.a2=etmp._esi_allele2[j];
                    snpinfotmp.freq=etmp._esi_freq[j];
                    snpinfo.push_back(snpinfotmp);
                }
            }
        }
        printf("Total %ld SNPs to be included.\n",snpinfo.size());
    }
    void combine_epi(vector<probeinfolst2> &probeinfo, vector<string> &smasNames)
    {
        long counter = 0;
        map<string, int> prb_map;
        map<string, int> prbbp_map;
        map<string, int>::iterator iter;
        for (int i = 0; i < smasNames.size(); i++)
        {
            eqtlInfo etmp;
            string epifile = smasNames[i]+".epi";
            read_epifile(&etmp, epifile);
            for (int j = 0; j<etmp._probNum; j++)
            {
                string crsbpstr=etmp._epi_prbID[j]+":"+atos(etmp._epi_bp[j]);
                prb_map.insert(pair<string, int>(etmp._epi_prbID[j].c_str(), counter));
                prbbp_map.insert(pair<string, int>(crsbpstr.c_str(), counter));
                if(prb_map.size() != prbbp_map.size())
                {
                    printf("ERROR: Probe %s has multiple BPs, please make sure that the genomic coordinates for the probes are on the same page.\n", etmp._epi_prbID[j].c_str()) ;
                    exit(EXIT_FAILURE);
                }
               
                if (counter < prb_map.size())
                {
                    probeinfolst2 probinfotmp;
                    counter=prb_map.size();
                    probinfotmp.probechr=etmp._epi_chr[j];
                     probinfotmp.probeId=etmp._epi_prbID[j];
                     probinfotmp.bp=etmp._epi_bp[j];
                    probinfotmp.gd=etmp._epi_gd[j];
                     probinfotmp.genename=etmp._epi_gene[j];
                     probinfotmp.orien=etmp._epi_orien[j];
                    probinfotmp.vnum=0;
                    probinfotmp.rowid=NULL;
                    probinfotmp.beta_se=NULL;
                    probinfotmp.besdpath.push_back(smasNames[i]);
                    probeinfo.push_back(probinfotmp);
                    
                } else {
                    iter=prb_map.find(etmp._epi_prbID[j]);
                    if(iter!=prb_map.end())
                    {
                        probeinfo[iter->second].besdpath.push_back(smasNames[i]);
                    }
                    else
                    {
                        printf("ERROR: This would never happen. please help to report this bug.\n") ;
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
        printf("Total %ld probes to be included.\n",probeinfo.size());
    }
    
    float readfloat(FILE *f) {
        float v;
        fread((void*)(&v), sizeof(v), 1, f);
        return v;
    }
    uint64_t readuint64(FILE *f) {
        uint64_t v;
        fread((void*)(&v), sizeof(v), 1, f);
        return v;
    }
    uint64_t countNotNullNum(vector<string> &smasNames)
    {
        uint64_t count=0;
        for (int i = 0; i < smasNames.size(); i++)
        {
            eqtlInfo etmp;
            string besdfile = smasNames[i]+".besd";
            FILE *fptr=fopen(besdfile.c_str(), "rb");
            if(!fptr)
            {
                printf ( "ERROR: Couldn't open file %s\n", besdfile.c_str());
                exit (EXIT_FAILURE);
            }
            float filetype=readfloat(fptr);
            uint64_t valnum=0;
            if((int)filetype==SPARSE_FILE_TYPE_3 || (int)filetype==SPARSE_FILE_TYPE_2) valnum=readuint64(fptr);
            if((int)filetype==SPARSE_FILE_TYPE_1) valnum=(uint64_t)readfloat(fptr);
            if((int)filetype==DENSE_FILE_TYPE_1) {
                while(!feof(fptr))
                {
                    float tmpfloat=readfloat(fptr);
                    if(abs(tmpfloat+9)>1e-6) valnum++;
                }
            }
           fclose(fptr);
           count+=(valnum>>1);
        }
        
        return count;
    }
    void save_besds_dbesd(char* outFileName, vector<snpinfolst> &snpinfo, vector<probeinfolst2> &probeinfo, vector<string> &esi_rs,vector<string> &esi_a1,vector<string> &esi_a2, bool genouni)
    {
        map<string, int> esi_map;
        for(int j=0;j<esi_rs.size();j++)
        {
            esi_map.insert(pair<string,int>(esi_rs[j],j));
        }
        
        // get esd info
        long esiNum=snpinfo.size();
        long epiNum=probeinfo.size();
        string esdfile=string(outFileName)+string(".besd");
        FILE * smr1;
        smr1 = fopen (esdfile.c_str(), "wb");
        if (!(smr1)) {
            printf("ERROR: Failed to open file %s.\n",esdfile.c_str());
            exit(EXIT_FAILURE);
        }
        float filetype=DENSE_FILE_TYPE_1;
        fwrite (&filetype,sizeof(float), 1, smr1);
        
        uint64_t bsize=(uint64_t)esiNum<<1;
        float* buffer=(float*)malloc (sizeof(float)*bsize);
        if (NULL == buffer) {
            fprintf(stderr, "Malloc failed\n");
            exit(-1);
        }
        bool prtscr=false;
        map<string, int>::iterator iter;
        eqtlInfo eqtlinfo;
        if(genouni)
        {
            read_esifile(&eqtlinfo, string(probeinfo[0].besdpath[0])+".esi", prtscr);
        }
        for(int j=0;j<epiNum;j++)
        {
            printf("Saving... %3.0f%%\r", 100.0*j/epiNum);
            fflush(stdout);
            string prbname=probeinfo[j].probeId;
            for(int k=0;k<bsize;k++) buffer[k]=-9; //init
            for(int k=0;k<probeinfo[j].besdpath.size();k++)
            {
                string besdfilepre = probeinfo[j].besdpath[k];
                vector<string> _rs;
                vector<float> _beta;
                vector<float> _se;
                vector<int> _chr;
                vector<string> _a1;
                vector<string> _a2;
                vector<int> _bp;
                vector<float> _freq;
                if(genouni)
                {
                    read_epifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".epi", prtscr);
                    extract_eqtl_single_probe(&eqtlinfo, prbname, prtscr);
                    read_besdfile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".besd", prtscr);
                    
                } else {
                    read_epifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".epi", prtscr);
                    extract_eqtl_single_probe(&eqtlinfo, prbname, prtscr);
                    read_esifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".esi", prtscr);
                    read_besdfile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".besd", prtscr);
                }
                
                if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
                {
                    //printf("No data included of probe %s from %s.\n",prbname.c_str(),probeinfo[j].besdpath[k].c_str());
                    continue;
                    
                }
                if(eqtlinfo._valNum==0)
                {
                    for(uint32_t jj=0;jj<eqtlinfo._snpNum;jj++)
                    {
                        float beta=eqtlinfo._bxz[0][jj];
                        float se=eqtlinfo._sexz[0][jj];
                        if(ABS(se+9)<1e-6) continue;
                        _beta.push_back(beta);
                        _se.push_back(se);
                        _rs.push_back(eqtlinfo._esi_rs[jj]);
                        _chr.push_back(eqtlinfo._esi_chr[jj]);
                        _a1.push_back(eqtlinfo._esi_allele1[jj]);
                        _a2.push_back(eqtlinfo._esi_allele2[jj]);
                        _bp.push_back(eqtlinfo._esi_bp[jj]);
                        _freq.push_back(eqtlinfo._esi_freq[jj]);
                    }
                }
                else
                {
                    if(eqtlinfo._val.size()==0)
                    {
                        throw ("Error: No data extracted from the input, please check.\n");
                    }
                    
                    for(uint32_t ii=0;ii<eqtlinfo._probNum;ii++)
                    {
                        uint64_t proid=eqtlinfo._include[ii];
                        uint64_t pos=eqtlinfo._cols[proid<<1];
                        uint64_t pos1=eqtlinfo._cols[(proid<<1)+1];
                        uint64_t num=pos1-pos;
                        for(int jj=0;jj<num;jj++)
                        {
                            double beta=eqtlinfo._val[pos+jj];
                            double se=eqtlinfo._val[pos+jj+num];
                            int rowid=eqtlinfo._rowid[pos+jj];
                            _beta.push_back(beta);
                            _se.push_back(se);
                            _rs.push_back(eqtlinfo._esi_rs[rowid]);
                            _chr.push_back(eqtlinfo._esi_chr[rowid]);
                            _a1.push_back(eqtlinfo._esi_allele1[rowid]);
                            _a2.push_back(eqtlinfo._esi_allele2[rowid]);
                            _bp.push_back(eqtlinfo._esi_bp[rowid]);
                            _freq.push_back(eqtlinfo._esi_freq[rowid]);
                            
                        }
                    }
                }
                
                
                vector<int> rsid(_rs.size());
#pragma omp parallel for private(iter)
                for (int l = 0; l<_rs.size(); l++){
                    iter = esi_map.find(_rs[l]);
                    if (iter != esi_map.end()) rsid[l]=iter->second;
                    else {
                        printf("ERROR: SNP is not in SNP map. if you are using --geno-uni, please disable it then try again, else please report this bug.");
                        exit(EXIT_FAILURE);
                    }
                }
                for(int l=0;l<rsid.size();l++)
                {
                    if(esi_a1[rsid[l]]==_a1[l] && esi_a2[rsid[l]]==_a2[l]){
                        buffer[rsid[l]]=_beta[l];
                        buffer[esiNum+rsid[l]]=_se[l];
                    } else if(esi_a1[rsid[l]]==_a2[l] && esi_a2[rsid[l]]==_a1[l]){
                        buffer[rsid[l]]=-1.0*_beta[l];
                        buffer[esiNum+rsid[l]]=_se[l];
                    } else {
                        int did=-9;
                        float sig=1.0;
                        for(int m=0;m<esi_rs.size();m++)
                        {
                            if(esi_rs[m]==_rs[l])
                            {
                                if(esi_a1[m]==_a1[l] && esi_a2[m]==_a2[l])
                                {
                                    did=m;
                                    break;
                                }
                                if(esi_a1[m]==_a2[l] && esi_a2[m]==_a1[l])
                                {
                                    did=m;
                                    sig=-1.0;
                                    break;
                                }
                            }
                        }
                        if(did==-9)
                        {
                            printf("ERROR: This would not go to happen. Please report this bug.");
                            exit(EXIT_FAILURE);
                        }
                        buffer[did]=sig*_beta[l];
                        buffer[esiNum+did]=_se[l];
                    }
                }
            }
            fwrite (buffer,sizeof(float), bsize, smr1);
        }
        fclose (smr1);
        free(buffer);
        cout<<"Beta values and SE values for "<<epiNum<<" Probes and "<<esiNum<<" SNPs have been saved in the dense binary file [" + esdfile + "]." <<endl;
        
    }

    void save_besds_dbesd(char* outFileName, vector<snpinfolst> &snpinfo, vector<probeinfolst2> &probeinfo)
    {
        //because read only once .esi file. so the alleles of each SNP in every .esi should be the same and in the same alignment.
         // get esd info
        long esiNum=snpinfo.size();
        long epiNum=probeinfo.size();
        string esdfile=string(outFileName)+string(".besd");
        FILE * smr1;
        smr1 = fopen (esdfile.c_str(), "wb");
        if (!(smr1)) {
            printf("ERROR: Failed to open file %s.\n",esdfile.c_str());
            exit(EXIT_FAILURE);
        }
        float filetype=DENSE_FILE_TYPE_1;
        fwrite (&filetype,sizeof(float), 1, smr1);
        
        uint64_t sizeperprb=sizeof(float)*esiNum*2;
        uint64_t bsize=0x7FFFFFC0;
        float* buffer=(float*)malloc (bsize);
        while(NULL == buffer && bsize>sizeperprb) {
            buffer=(float*)malloc (bsize>>=1);
        }
        if(NULL == buffer)
        {
            printf("ERROR: Can't malloc the reading buffer for %llu MB.\n",(bsize>>20));
            exit(EXIT_FAILURE);
        }
        long prbperloop=bsize/sizeperprb;
        bool prtscr=false;
        map<string, int>::iterator iter;
        eqtlInfo eqtlinfo;
        read_esifile(&eqtlinfo, string(probeinfo[0].besdpath[0])+".esi", prtscr);
       if(esiNum != eqtlinfo._snpNum)
       {
           printf("ERROR: The SNPs in your .esi files are not in consistency, please disable --geno-uni, then have a try.\n");
           exit(EXIT_FAILURE);
       }
        int loops=ceil(1.0*epiNum/prbperloop);
        for(int j=0;j<loops;j++)
        {
            printf("Saving... %3.0f%%\r", 100.0*j/loops);
            fflush(stdout);
            uint64_t numprbcurloop=prbperloop;
            if(j==loops-1) numprbcurloop=epiNum-j*numprbcurloop;
            uint64_t vnum=numprbcurloop*esiNum*2;
            for(int k=0;k<vnum;k++) buffer[k]=-9; //init
            map<string,int> fcurloop;
            long fnum=0;
            vector<F2Prb> f2prb;
            for(int k=0;k<numprbcurloop;k++)
            {
                long curPrid=j*prbperloop+k;
                if(curPrid>=probeinfo.size()) break;
                string prbname=probeinfo[curPrid].probeId;
                if(probeinfo[curPrid].besdpath.size()>1)
                {
                    printf("ERROR: Probe %s is found in %ld BESD files. please check wether the SNPs are in consistency among the .esi files, if yes, please remove one of the duplicate probe then have a try, if not, please disable --geno-uni, then have a try.\n",prbname.c_str(),probeinfo[curPrid].besdpath.size());
                    exit(EXIT_FAILURE);
                }
                iter=fcurloop.find(probeinfo[curPrid].besdpath[0]);
                if(iter!=fcurloop.end())
                {
                    f2prb[iter->second].pid.push_back(k);
                } else {
                    fcurloop.insert(pair<string,int>(probeinfo[curPrid].besdpath[0],fnum));
                    F2Prb fptmp;
                    fptmp.besdpath=probeinfo[curPrid].besdpath[0];
                    fptmp.pid.push_back(k);
                    f2prb.push_back(fptmp);
                    fnum++;
                }
               
            }
            
            for(int l=0;l<f2prb.size();l++)
            {
                
                string besdfilepre = f2prb[l].besdpath;
                read_epifile(&eqtlinfo, besdfilepre+".epi", prtscr);
                
                string besdfile = string(besdfilepre)+".besd";
                FILE *fptr=fopen(besdfile.c_str(), "rb");
                if(!fptr)
                {
                    printf ( "ERROR: Couldn't open file %s\n", besdfile.c_str());
                    exit (EXIT_FAILURE);
                }
                
                float filetype=readfloat(fptr);
               
                uint64_t valNum=0;
                uint64_t* ptr=NULL;
                uint64_t rowSTART=0;
                uint64_t valSTART=0;
                if((int)filetype==SPARSE_FILE_TYPE_3 ){
                    uint64_t colNum=(eqtlinfo._probNum<<1)+1;
                    fseek(fptr, 0L, SEEK_END);
                    uint64_t lSize = ftell(fptr);
                    fseek(fptr, 0L, SEEK_SET);
                    readfloat(fptr);
                    valNum=readuint64(fptr);
                    if( lSize - (sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float)) != 0) {fputs ("wrong element number",stderr); exit (3);}
                    
                    
                    uint64_t colsize=colNum*sizeof(uint64_t);
                    uint64_t* colbuf=(uint64_t*)malloc(colsize);
                    if(NULL == colbuf)
                    {
                        printf("ERROR: Can't malloc the reading cols buffer for %llu MB.\n",(colsize>>20));
                        exit(EXIT_FAILURE);
                    }
                    fread(colbuf,colNum,sizeof(uint64_t),fptr);
                    
                    ptr=colbuf;
                    rowSTART=sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t);
                    valSTART=sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t)+valNum*sizeof(uint32_t);
                    
                    for(int m=0;m<f2prb[l].pid.size();m++)
                    {
                        long curPrid=j*prbperloop+f2prb[l].pid[m];
                        string curprb=probeinfo[curPrid].probeId;
                        iter=eqtlinfo._probe_name_map.find(curprb);
                        if(iter!=eqtlinfo._probe_name_map.end())
                        {
                            uint64_t pid=iter->second;
                            uint64_t pos=*(ptr+(pid<<1)); //BETA START
                            uint64_t pos1=*(ptr+(pid<<1)+1); //SE START
                            uint64_t num=pos1-pos;
                            
                            char* row_char_ptr;
                            row_char_ptr = (char*) malloc (sizeof(char)*2*num*sizeof(uint32_t));
                            if (row_char_ptr == NULL) {fputs ("Memory error",stderr); exit (1);}
                            char* val_char_ptr;
                            val_char_ptr = (char*) malloc (sizeof(char)*2*num*sizeof(float));
                            if (val_char_ptr == NULL) {fputs ("Memory error",stderr); exit (1);}
                            memset(row_char_ptr,0,sizeof(char)*2*num*sizeof(uint32_t));
                            memset(val_char_ptr,0,sizeof(char)*2*num*sizeof(float));
                            fseek(fptr, rowSTART+pos*sizeof(uint32_t), SEEK_SET);
                            fread(row_char_ptr, sizeof(uint32_t),2*num,fptr);
                            uint32_t* row_ptr=(uint32_t *)row_char_ptr;
                            fseek(fptr,valSTART+pos*sizeof(float),SEEK_SET);
                            fread(val_char_ptr,sizeof(float), 2*num,fptr);
                            float* val_ptr=(float*)val_char_ptr;
                            for(int j=0;j<num;j++)
                            {
                                uint32_t rid=*(row_ptr+j);
                                *(buffer+f2prb[l].pid[m]*esiNum*2+rid)=*(val_ptr+j);
                                *(buffer+f2prb[l].pid[m]*esiNum*2+rid+esiNum)=*(val_ptr+j+num);
                            }

                            free(row_char_ptr);
                            free(val_char_ptr);
                        }
                    }
                    free(colbuf);
                }
                else if((int)filetype==DENSE_FILE_TYPE_1){
                     for(int m=0;m<f2prb[l].pid.size();m++)
                     {
                         long curPrid=j*prbperloop+f2prb[l].pid[m];
                         string curprb=probeinfo[curPrid].probeId;
                         iter=eqtlinfo._probe_name_map.find(curprb);
                         if(iter!=eqtlinfo._probe_name_map.end())
                         {
                             uint64_t pid=iter->second;
                             fseek(fptr,((pid<<1)*eqtlinfo._snpNum+1)<<2, SEEK_SET);
                             float* wptr=buffer+f2prb[l].pid[m]*esiNum*2;
                             fread(wptr, sizeof(char),eqtlinfo._snpNum<<3,fptr);
                         }
                     }
                }
                else {
                    printf("Your file is in the old sparse format. please first re-make it by --beqtl-summary and --make-besd.\n");
                    exit(EXIT_FAILURE);
                }
                fclose(fptr);
            }
            fwrite (buffer,sizeof(float), vnum, smr1);
        }
        fclose (smr1);
        free(buffer);
        cout<<"Beta values and SE values for "<<epiNum<<" Probes and "<<esiNum<<" SNPs have been saved in the dense binary file [" + esdfile + "]." <<endl;
        
    }
    
    void save_besds_sbesd(char* outFileName, vector<snpinfolst> &snpinfo, vector<probeinfolst2> &probeinfo, vector<string> &esi_rs,vector<string> &esi_a1,vector<string> &esi_a2, bool genouni)
    {
        map<string, int> esi_map;
        for(int j=0;j<esi_rs.size();j++)
        {
            esi_map.insert(pair<string,int>(esi_rs[j],j));
        }
        // get esd info
        long esiNum=snpinfo.size();
        long epiNum=probeinfo.size();
        string esdfile=string(outFileName)+string(".besd");
        FILE * smr1;
        smr1 = fopen (esdfile.c_str(), "wb");
        if (!(smr1)) {
            printf("ERROR: Failed to open file %s.\n",esdfile.c_str());
            exit(EXIT_FAILURE);
        }
        float filetype=SPARSE_FILE_TYPE_3;
        fwrite (&filetype,sizeof(float), 1, smr1);
        
        vector<uint64_t> cols((epiNum<<1)+1);;
        vector<uint32_t> rowids;
        vector<float> val;
        cols[0]=0;
        bool prtscr=false;
        map<string, int>::iterator iter;
        eqtlInfo eqtlinfo;
        if(genouni)
        {
            read_esifile(&eqtlinfo, string(probeinfo[0].besdpath[0])+".esi", prtscr);
        }
        for(int j=0;j<epiNum;j++)
        {
            printf("Saving... %3.0f%%\r", 100.0*j/epiNum);
            fflush(stdout);
            string prbname=probeinfo[j].probeId;
            map<string, int> rsa_map;
            long rsNum=0;
            
            vector<uint32_t> tmprid;
            vector<float> tmpse;
            for(int k=0;k<probeinfo[j].besdpath.size();k++)
            {
                string besdfilepre = probeinfo[j].besdpath[k];
                vector<string> _rs;
                vector<float> _beta;
                vector<float> _se;
                vector<int> _chr;
                vector<string> _a1;
                vector<string> _a2;
                vector<int> _bp;
                if(genouni)
                {
                    read_epifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".epi", prtscr);
                    extract_eqtl_single_probe(&eqtlinfo, prbname, prtscr);
                    read_besdfile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".besd", prtscr);
                    
                } else {
                    read_epifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".epi", prtscr);
                    extract_eqtl_single_probe(&eqtlinfo, prbname, prtscr);
                    read_esifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".esi", prtscr);
                    read_besdfile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".besd", prtscr);
                }
                if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
                {
                    //printf("No data included of probe %s from %s.\n",prbname.c_str(),probeinfo[j].besdpath[k].c_str());
                    continue;
                    
                }
                if(eqtlinfo._valNum==0)
                {
                    for(uint32_t jj=0;jj<eqtlinfo._snpNum;jj++)
                    {
                        float beta=eqtlinfo._bxz[0][jj];
                        float se=eqtlinfo._sexz[0][jj];
                        if(ABS(se+9)<1e-6) continue;
                        _beta.push_back(beta);
                        _se.push_back(se);
                        _rs.push_back(eqtlinfo._esi_rs[jj]);
                        _chr.push_back(eqtlinfo._esi_chr[jj]);
                        _a1.push_back(eqtlinfo._esi_allele1[jj]);
                        _a2.push_back(eqtlinfo._esi_allele2[jj]);
                        _bp.push_back(eqtlinfo._esi_bp[jj]);
                    }
                }
                else
                {
                    if(eqtlinfo._val.size()==0)
                    {
                        throw ("Error: No data extracted from the input, please check.\n");
                    }
                    
                    for(uint32_t ii=0;ii<eqtlinfo._probNum;ii++)
                    {
                        uint64_t proid=eqtlinfo._include[ii];
                        uint64_t pos=eqtlinfo._cols[proid<<1];
                        uint64_t pos1=eqtlinfo._cols[(proid<<1)+1];
                        uint64_t num=pos1-pos;
                        for(int jj=0;jj<num;jj++)
                        {
                            double beta=eqtlinfo._val[pos+jj];
                            double se=eqtlinfo._val[pos+jj+num];
                            int rowid=eqtlinfo._rowid[pos+jj];
                            _beta.push_back(beta);
                            _se.push_back(se);
                            _rs.push_back(eqtlinfo._esi_rs[rowid]);
                            _chr.push_back(eqtlinfo._esi_chr[rowid]);
                            _a1.push_back(eqtlinfo._esi_allele1[rowid]);
                            _a2.push_back(eqtlinfo._esi_allele2[rowid]);
                            _bp.push_back(eqtlinfo._esi_bp[rowid]);
                            
                        }
                    }
                }
                vector<int> rsid(_rs.size());
                #pragma omp parallel for private(iter)
                for (int l = 0; l<_rs.size(); l++){
                    iter = esi_map.find(_rs[l]);
                    if (iter != esi_map.end()) rsid[l]=iter->second;
                    else {
                        printf("ERROR: SNP is not in SNP map. if you are using --geno-uni, please disable it then try again, else please report this bug.");
                        exit(EXIT_FAILURE);
                    }
                }
                for(int l=0;l<rsid.size();l++)
                {
                    if(abs(_se[l]+9)>1e-6)
                    {
                        string chckstr=_rs[l]+":"+_a1[l]+":"+_a2[l];
                        rsa_map.insert(pair<string,int>(chckstr,l));
                        if(rsNum<rsa_map.size())
                        {
                            if(esi_a1[rsid[l]]==_a1[l] && esi_a2[rsid[l]]==_a2[l] ){
                                val.push_back(_beta[l]);
                                rowids.push_back(rsid[l]);
                                tmpse.push_back(_se[l]);
                                tmprid.push_back(rsid[l]);
                            } else if(esi_a1[rsid[l]]==_a2[l] && esi_a2[rsid[l]]==_a1[l] ){
                                val.push_back(-1.0*_beta[l]);
                                rowids.push_back(rsid[l]);
                                tmpse.push_back(_se[l]);
                                tmprid.push_back(rsid[l]);
                            } else {
                                int did=-9;
                                float sig=1.0;
                                for(int m=0;m<esi_rs.size();m++)
                                {
                                    if(esi_rs[m]==_rs[l])
                                    {
                                        if(esi_a1[m]==_a1[l] && esi_a2[m]==_a2[l])
                                        {
                                            did=m;
                                            break;
                                        }
                                        if(esi_a1[m]==_a2[l] && esi_a2[m]==_a1[l])
                                        {
                                            did=m;
                                            sig=-1.0;
                                            break;
                                        }
                                    }
                                }
                                if(did==-9)
                                {
                                    printf("ERROR: This would not go to happen. Please report this bug.");
                                    exit(EXIT_FAILURE);
                                }
                                val.push_back(sig*_beta[l]);
                                rowids.push_back(did);
                                tmpse.push_back(_se[l]);
                                tmprid.push_back(did);
                            }
                            rsNum=rsa_map.size();
                        } else {
                            printf("WARNING: duplicate SNP %s with the same alleles belonging to the same probe %s found in current summary data file\"%s\" and the esd file(s) read before. \n",_rs[l].c_str(),probeinfo[j].probeId.c_str(),probeinfo[j].besdpath[k].c_str());
                        }
                    } else {
                        printf("WARNING: SNP %s  with \"NA\" value found and skipped.\n",_rs[l].c_str());
                    }
                    
                    
                }
            }
            for(int k=0;k<tmpse.size();k++)
            {
                val.push_back(tmpse[k]);
                rowids.push_back(tmprid[k]);
            }
            uint64_t real_num=tmpse.size();
            cols[(j<<1)+1]=real_num+cols[j<<1];
            cols[j+1<<1]=(real_num<<1)+cols[j<<1];
        }
        uint64_t valNum=val.size();
        fwrite (&valNum,sizeof(uint64_t), 1, smr1);
        fwrite (&cols[0],sizeof(uint64_t), cols.size(), smr1);
        fwrite (&rowids[0],sizeof(uint32_t), rowids.size(), smr1);
        fwrite (&val[0],sizeof(float), val.size(), smr1);
        fclose (smr1);
        
        cout<<"Beta values and SE values for "<<epiNum<<" Probes and "<<esiNum<<" SNPs have been saved in the sparse binary file [" + esdfile + "]." <<endl;
    }
    void save_besds_sbesd(char* outFileName, vector<snpinfolst> &snpinfo, vector<probeinfolst2> &probeinfo ,vector<string> &smasNames)
    {
        map<string, int> epi_map;
        map<string,int>::iterator iter;
        for(int j=0;j<probeinfo.size();j++)
        {
            epi_map.insert(pair<string,int>(probeinfo[j].probeId,j));
        }
        // get esd info
        long esiNum=snpinfo.size();
        long epiNum=probeinfo.size();
        string esdfile=string(outFileName)+string(".besd");
        FILE * smr1;
        smr1 = fopen (esdfile.c_str(), "wb");
        if (!(smr1)) {
            printf("ERROR: Failed to open file %s.\n",esdfile.c_str());
            exit(EXIT_FAILURE);
        }
        float filetype=SPARSE_FILE_TYPE_3;
        fwrite (&filetype,sizeof(float), 1, smr1);
        
         bool prtscr=false;
       
        eqtlInfo etmp;
        //read_esifile(&etmp, string(probeinfo[0].besdpath[0])+".esi", prtscr);
        //if(esiNum != etmp._snpNum)
        //{
        //    printf("ERROR: The SNPs in your .esi files are not in consistency, please disable --geno-uni, then have a try.\n");
        //    exit(EXIT_FAILURE);
        //}
        printf("Reading besd files....\n");

        for (int i = 0; i < smasNames.size(); i++)
        {
            printf("Reading... %3.0f%%\r", 100.0*i/(smasNames.size()));
            fflush(stdout);
            
            string epifile = smasNames[i]+".epi";
            read_epifile(&etmp, epifile, prtscr);

            string besdfile = smasNames[i]+".besd";
            FILE *fptr=fopen(besdfile.c_str(), "rb");
            if(!fptr)
            {
                printf ( "ERROR: Couldn't open file %s\n", besdfile.c_str());
                exit (EXIT_FAILURE);
            }
            float filetype=readfloat(fptr);
            uint64_t valNum=0;
            uint64_t* ptr=NULL;
            uint64_t rowSTART=0;
            uint64_t valSTART=0;
            if((int)filetype==SPARSE_FILE_TYPE_3 ){
                uint64_t colNum=(etmp._probNum<<1)+1;
                fseek(fptr, 0L, SEEK_END);
                uint64_t lSize = ftell(fptr);
                fseek(fptr, 0L, SEEK_SET);
                readfloat(fptr);
                valNum=readuint64(fptr);
                if( lSize - (sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float)) != 0) {fputs ("wrong element number",stderr); exit (3);}
                
                
                uint64_t colsize=colNum*sizeof(uint64_t);
                uint64_t* colbuf=(uint64_t*)malloc(colsize);
                if(NULL == colbuf)
                {
                    printf("ERROR: Can't malloc the reading cols buffer for %llu MB.\n",(colsize>>20));
                    exit(EXIT_FAILURE);
                }
                fread(colbuf,colNum,sizeof(uint64_t),fptr);
                
                ptr=colbuf;
                rowSTART=sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t);
                valSTART=sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t)+valNum*sizeof(uint32_t);
                
                for( int j=0;j<etmp._probNum;j++)
                {
                    string curprb=etmp._epi_prbID[j];
                    iter=epi_map.find(curprb);
                    int prbindx=iter->second;
                    if(probeinfo[prbindx].besdpath.size()>1) // or if(probeinfo[prbindx].beta_se.size()>0)
                    {
                        printf("ERROR: Probe %s is found in %ld BESD files. please check wether the SNPs are in consistency among .esi files, if yes, please remove one of the duplicate probe then have a try, if no, please disable --geno-uni, then have a try.\n",curprb.c_str(),probeinfo[prbindx].besdpath.size());
                        exit(EXIT_FAILURE);
                    }
                    
                    uint64_t pos=*(ptr+(j<<1)); //BETA START
                    uint64_t pos1=*(ptr+(j<<1)+1); //SE START
                    uint64_t num=pos1-pos;
                    if(num>0)
                    {
                        uint32_t* ridbuff=(uint32_t*)malloc(num*2*sizeof(uint32_t));
                        if(NULL == ridbuff)
                        {
                            printf("ERROR: Memory allocation error when when reading probe %s in the file %s.\n",curprb.c_str(),besdfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                        memset(ridbuff,0,num*2*sizeof(uint32_t));
                        float* betasebuff=(float*)malloc(num*2*sizeof(float));
                        if(NULL == betasebuff)
                        {
                            printf("ERROR: Memory allocation error when when reading probe %s in the file %s.\n",curprb.c_str(),besdfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                         memset(betasebuff,0,num*2*sizeof(float));
                        
                        fseek(fptr, rowSTART+pos*sizeof(uint32_t), SEEK_SET);
                        fread(ridbuff, sizeof(uint32_t),2*num,fptr);
                        probeinfo[prbindx].rowid=ridbuff;
                        fseek(fptr,valSTART+pos*sizeof(float),SEEK_SET);
                        fread(betasebuff,sizeof(float), 2*num,fptr);
                        probeinfo[prbindx].beta_se=betasebuff;
                        probeinfo[prbindx].vnum=num;
                        
                        // no memory free here. later.
                    } else {
                        //printf("Probe %s has no values.\n",curprb.c_str());
                    }
                    
                }
                
                free(colbuf);
            }
            else if((int)filetype==DENSE_FILE_TYPE_1){
                
                float* tmpbetase=(float*)malloc(sizeof(float)*etmp._snpNum<<1);
                if(NULL == tmpbetase)
                {
                    printf("ERROR: Can't malloc the reading cols buffer for %u MB.\n",(etmp._snpNum>>17));
                    exit(EXIT_FAILURE);
                }

                for( int j=0;j<etmp._probNum;j++)
                {
                    string curprb=etmp._epi_prbID[j];
                    iter=epi_map.find(curprb);
                    int prbindx=iter->second;
                    if(probeinfo[prbindx].besdpath.size()>1)
                    {
                        printf("ERROR: Probe %s is found in %ld BESD files. please check wether the SNPs are in consistency among .esi files, if yes, please remove one of the duplicate probe then have a try, if no, please disable --geno-uni, then have a try.\n",curprb.c_str(),probeinfo[prbindx].besdpath.size());
                        exit(EXIT_FAILURE);
                    }
                    if(probeinfo[prbindx].vnum>0)
                    {
                        printf("ERROR: Probe %s is found in multiple BESD files. please check wether the SNPs are in consistency among .esi files, if yes, please remove one of the duplicate probe then have a try, if no, please disable --geno-uni, then have a try.\n",curprb.c_str());
                        exit(EXIT_FAILURE);
                    }
                    memset(tmpbetase,0,sizeof(float)*etmp._snpNum<<1);
                    fseek(fptr,((j<<1)*etmp._snpNum+1)<<2, SEEK_SET);
                    fread(tmpbetase, sizeof(float),etmp._snpNum<<1,fptr);
                    uint64_t realnum=0;
                    for(int k=0;k<etmp._snpNum;k++) if(tmpbetase[etmp._snpNum+k]+9>1e-6) realnum++;
                    probeinfo[prbindx].vnum=realnum;
                    uint32_t* ridbuff=(uint32_t*)malloc(realnum*2*sizeof(uint32_t));
                    if(NULL == ridbuff)
                    {
                        printf("ERROR: Memory allocation error when when reading probe %s in the file %s.\n",curprb.c_str(),besdfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    memset(ridbuff,0,realnum*2*sizeof(uint32_t));
                    float* betasebuff=(float*)malloc(realnum*2*sizeof(float));
                    if(NULL == betasebuff)
                    {
                        printf("ERROR: Memory allocation error when when reading probe %s in the file %s.\n",curprb.c_str(),besdfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    memset(betasebuff,0,realnum*2*sizeof(uint32_t));
                    realnum=0;
                    for(int k=0;k<etmp._snpNum;k++)
                    {
                        if(tmpbetase[etmp._snpNum+k]+9>1e-6){
                            ridbuff[realnum]=k;
                            betasebuff[realnum]=tmpbetase[k];
                            realnum++;
                        }
                    }
                    for(int k=0;k<etmp._snpNum;k++)
                    {
                        if(tmpbetase[etmp._snpNum+k]+9>1e-6){
                            ridbuff[realnum]=k;
                            betasebuff[realnum]=tmpbetase[etmp._snpNum+k];
                            realnum++;
                        }
                    }
                     probeinfo[prbindx].rowid=ridbuff;
                     probeinfo[prbindx].beta_se=betasebuff;

                }
                  free(tmpbetase);
            }
            else {
                printf("Your file is in the old sparse format. please first re-make it by --beqtl-summary and --make-besd.\n");
                exit(EXIT_FAILURE);
            }
            fclose(fptr);
            
        }

        vector<uint64_t> cols((epiNum<<1)+1);
        uint64_t valNum=0;
        cols[0]=0;
       
        for(int j=0;j<epiNum;j++)
        {
            uint64_t real_num=probeinfo[j].vnum;
            cols[(j<<1)+1]=real_num+cols[j<<1];
            cols[j+1<<1]=(real_num<<1)+cols[j<<1];
            valNum+=real_num*2;
        }
        
        fwrite (&valNum,sizeof(uint64_t), 1, smr1);
        fwrite (&cols[0],sizeof(uint64_t), cols.size(), smr1);
        
        for(int j=0;j<epiNum;j++)
        {
            printf("Saving...  %3.0f%%\r", 100.0*j/(2*epiNum));
            fflush(stdout);
            fwrite (probeinfo[j].rowid,sizeof(uint32_t), probeinfo[j].vnum*2, smr1);
        }

        for(int j=0;j<epiNum;j++)
        {
            printf("Saving...  %3.0f%%\r", (100.0*j/(2*epiNum)+50));
            fflush(stdout);
            fwrite (probeinfo[j].beta_se,sizeof(float), probeinfo[j].vnum*2, smr1);
        }
        fclose (smr1);
        for(int j=0;j<epiNum;j++) {
            free(probeinfo[j].beta_se);
            probeinfo[j].beta_se=NULL;
            free(probeinfo[j].rowid);
            probeinfo[j].rowid=NULL;
        }
        cout<<"Beta values and SE values for "<<epiNum<<" Probes and "<<esiNum<<" SNPs have been saved in the sparse binary file [" + esdfile + "]." <<endl;
    }
    
    void save_slct_besds_sbesd(char* outFileName,  long esiNum, vector<probeinfolst2> &probeinfo,vector<string> &esi_rs,vector<string> &esi_a1,vector<string> &esi_a2,int cis_itvl,int trans_itvl,float transThres,float restThres, bool genouni)
    {
        map<string, int> esi_map;
        for(int j=0;j<esi_rs.size();j++)
        {
            esi_map.insert(pair<string,int>(esi_rs[j],j));
        }
        // get esd info
        long epiNum=probeinfo.size();
        string esdfile=string(outFileName)+string(".besd");
        FILE * smr1;
        smr1 = fopen (esdfile.c_str(), "wb");
        if (!(smr1)) {
            printf("ERROR: Failed to open file %s.\n",esdfile.c_str());
            exit(EXIT_FAILURE);
        }
        float filetype=SPARSE_FILE_TYPE_3;
        fwrite (&filetype,sizeof(float), 1, smr1);
        
        vector<uint64_t> cols((epiNum<<1)+1);;
        vector<uint32_t> rowids;
        vector<float> val;
        cols[0]=0;
        bool prtscr=false;
        map<string, int>::iterator iter;
        
        
        // log file
        FILE* logfile=NULL;
        string logfname = string(outFileName)+".summary";
        logfile=fopen(logfname.c_str(), "w");
        if (!(logfile)) {
            printf("Error: Failed to open log file.\n");
            exit(1);
        }
        string logstr="cis-window:\t"+itos(cis_itvl)+"Kb\ntrans-window:\t"+itos(trans_itvl)+"Kb\np-value threshold of trans:\t"+dtos(transThres)+"\np-value threshold of others:\t"+dtos(restThres)+"\n";
        logstr+="\ncis region is indicated by [Chr, Start bp, End bp, nsnp];\ntrans region is indicated by <Chr, Start bp, End bp, nsnp>;\nthe number of other SNPs selected is indicated by (NumSNPs beyond cis and trans).\n";
        
        logstr+="\n{ProbeID, ProbeChr, ProbeBP}\t[Chr,cis_startBP,cis_endBP,NumSNPs]\t<Chr,trans_startBP,trans_endBP,NumSNPs>\t(NumSNPs)\n";
        fputs(logstr.c_str(),logfile);
        fflush(logfile);
        cis_itvl=cis_itvl*1000;
        trans_itvl=trans_itvl*1000;
        vector<snpinfolst> snpinfo;
        eqtlInfo eqtlinfo;
        if(genouni)
        {
            read_esifile(&eqtlinfo, string(probeinfo[0].besdpath[0])+".esi", prtscr);
        }
        for(int j=0;j<epiNum;j++)
        {
            printf("Saving... %3.0f%%\r", 100.0*j/epiNum);
            fflush(stdout);
            string prbname=probeinfo[j].probeId;
            vector<uint32_t> tmprid;
            vector<float> tmpse;
            
            snpinfo.clear();
            map<string,int> snp_map;
            long snpmapsize=0;
            for(int k=0;k<probeinfo[j].besdpath.size();k++)
            {
                if(genouni)
                {
                    read_epifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".epi", prtscr);
                    extract_eqtl_single_probe(&eqtlinfo, prbname, prtscr);
                    read_besdfile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".besd", prtscr);
                    
                } else {
                    read_epifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".epi", prtscr);
                    extract_eqtl_single_probe(&eqtlinfo, prbname, prtscr);
                    read_esifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".esi", prtscr);
                    read_besdfile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".besd", prtscr);
                }
                if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
                {
                    //No data included of this file of the current probe does not mean no data in other files of the current probe;
                    continue;
                    
                }
                if(eqtlinfo._valNum==0)
                {
                    for(uint32_t jj=0;jj<eqtlinfo._snpNum;jj++)
                    {
                        float beta=eqtlinfo._bxz[0][jj];
                        float se=eqtlinfo._sexz[0][jj];
                        if(ABS(se+9)<1e-6) continue;
                        snp_map.insert(pair<string,int>(eqtlinfo._esi_rs[jj],snpmapsize));
                        if(snpmapsize<snp_map.size())
                        {
                            snpinfolst tmpinfo;
                            tmpinfo.beta=beta;
                            tmpinfo.se=se;
                            tmpinfo.snpchr=eqtlinfo._esi_chr[jj];
                            tmpinfo.snprs=eqtlinfo._esi_rs[jj];
                            tmpinfo.a1=eqtlinfo._esi_allele1[jj];
                            tmpinfo.a2=eqtlinfo._esi_allele2[jj];
                            tmpinfo.bp=eqtlinfo._esi_bp[jj];
                            snpinfo.push_back(tmpinfo);
                            snpmapsize=snp_map.size();
                        }
                        else{
                            printf("WARNING: duplicate SNP %s of probe %s found in different summary data files.\n",eqtlinfo._esi_rs[jj].c_str(),prbname.c_str());
                        }
                       
                    }
                }
                else
                {
                    if(eqtlinfo._val.size()==0)
                    {
                        throw ("Error: No data extracted from the input, please check.\n");
                    }
                    
                    for(uint32_t ii=0;ii<eqtlinfo._probNum;ii++) //eqtlinfo._probNum should be 1 here
                    {
                        uint64_t proid=eqtlinfo._include[ii];
                        uint64_t pos=eqtlinfo._cols[proid<<1];
                        uint64_t pos1=eqtlinfo._cols[(proid<<1)+1];
                        uint64_t num=pos1-pos;
                        for(int jj=0;jj<num;jj++)
                        {
                            double beta=eqtlinfo._val[pos+jj];
                            double se=eqtlinfo._val[pos+jj+num];
                            if(ABS(se+9)<1e-6) continue;
                            int rowid=eqtlinfo._rowid[pos+jj];
                            snp_map.insert(pair<string,int>(eqtlinfo._esi_rs[jj],snpmapsize));
                            if(snpmapsize<snp_map.size())
                            {
                                snpinfolst tmpinfo;
                                tmpinfo.beta=beta;
                                tmpinfo.se=se;
                                tmpinfo.snpchr=eqtlinfo._esi_chr[rowid];
                                tmpinfo.snprs=eqtlinfo._esi_rs[rowid];
                                tmpinfo.a1=eqtlinfo._esi_allele1[rowid];
                                tmpinfo.a2=eqtlinfo._esi_allele2[rowid];
                                tmpinfo.bp=eqtlinfo._esi_bp[rowid];
                                snpinfo.push_back(tmpinfo);
                                snpmapsize=snp_map.size();
                            }
                            else{
                                printf("WARNING: duplicate SNP %s of probe %s found in different summary data files.\n",eqtlinfo._esi_rs[jj].c_str(),prbname.c_str());
                            }
                        }
                    }
                }

            }
            if(snpinfo.size()>0)
            {
                snpinfolst* sortptr=&snpinfo[0];
                qsort(sortptr,snpinfo.size(),sizeof(snpinfolst),comp_esi); // when reading sparse, the rowids in file is not in order, so snpinfo is not in order, so sort should be implement here.
                
                probeinfolst prbifo;
                prbifo.bp=probeinfo[j].bp;
                prbifo.probechr=probeinfo[j].probechr;
                prbifo.probeId=probeinfo[j].probeId;
                
                vector<int> slct_idx;
                slct_sparse_per_prb(slct_idx, &prbifo, snpinfo,  cis_itvl,  trans_itvl, transThres, restThres,logfile); //slct_idx with no order if there are trans-rgeions
                stable_sort(slct_idx.begin(),slct_idx.end());
                vector<string> _rs(slct_idx.size()), _a1(slct_idx.size()),_a2(slct_idx.size());
                vector<float> _beta(slct_idx.size()), _se(slct_idx.size());
                
                for(int l=0;l<slct_idx.size();l++) {
                    _rs[l]=snpinfo[slct_idx[l]].snprs;
                    _a1[l]=snpinfo[slct_idx[l]].a1;
                    _a2[l]=snpinfo[slct_idx[l]].a2;
                    _beta[l]=snpinfo[slct_idx[l]].beta;
                    _se[l]=snpinfo[slct_idx[l]].se;
                }
                vector<int> rsid(_rs.size());
                #pragma omp parallel for private(iter)
                for (int l = 0; l<_rs.size(); l++){
                    iter = esi_map.find(_rs[l]);
                    if (iter != esi_map.end()) rsid[l]=iter->second;
                    else {
                        printf("ERROR: SNP is not in SNP map. if you are using --geno-uni, please disable it then try again, else please report this bug.");
                        exit(EXIT_FAILURE);
                    }
                }
               // map<string, int> rsa_map;
              //  long rsNum=0;
                for(int l=0;l<rsid.size();l++) // here rsid is not in order, so the rowids in the file is not in order
                {
                   // if(abs(_se[l]+9)>1e-6) // can move this. the NA is controled in slct_sparse_per_prb
                   // {
                    //    string chckstr=_rs[l]+":"+_a1[l]+":"+_a2[l];
                     //   rsa_map.insert(pair<string,int>(chckstr,l)); // in slct_sparse_per_prb, ras_map can privent selecting duplicate SNPs and double-slelecting SNPs. so we can move rsa_map here.
                     //   if(rsNum<rsa_map.size())
                      //  {
                            
                            if(esi_a1[rsid[l]]==_a1[l] && esi_a2[rsid[l]]==_a2[l] ){
                                val.push_back(_beta[l]);
                                rowids.push_back(rsid[l]);
                                tmpse.push_back(_se[l]);
                                tmprid.push_back(rsid[l]);
                            } else if(esi_a1[rsid[l]]==_a2[l] && esi_a2[rsid[l]]==_a1[l] ){
                                val.push_back(-1.0*_beta[l]);
                                rowids.push_back(rsid[l]);
                                tmpse.push_back(_se[l]);
                                tmprid.push_back(rsid[l]);
                            } else {
                                int did=-9;
                                float sig=1.0;
                                for(int m=0;m<esi_rs.size();m++)
                                {
                                    if(esi_rs[m]==_rs[l])
                                    {
                                        if(esi_a1[m]==_a1[l] && esi_a2[m]==_a2[l])
                                        {
                                            did=m;
                                            break;
                                        }
                                        if(esi_a1[m]==_a2[l] && esi_a2[m]==_a1[l])
                                        {
                                            did=m;
                                            sig=-1.0;
                                            break;
                                        }
                                    }
                                }
                                if(did==-9)
                                {
                                    printf("ERROR: This would not go to happen. Please report this bug.");
                                    exit(EXIT_FAILURE);
                                }
                                val.push_back(sig*_beta[l]);
                                rowids.push_back(did);
                                tmpse.push_back(_se[l]);
                                tmprid.push_back(did);
                            }
                            
                       //     rsNum=rsa_map.size();
                      //  } else {
                       //     printf("WARNING: duplicate SNP %s with the same alleles belonging to the same probe %s. \n",_rs[l].c_str(),prbname.c_str());
                       // }
                  //  } else {
                  //      printf("WARNING: SNP %s  with \"NA\" value found and skipped.\n",_rs[l].c_str());
                  //  }
                }
                for(int k=0;k<tmpse.size();k++)
                {
                    val.push_back(tmpse[k]);
                    rowids.push_back(tmprid[k]);
                }
                uint64_t real_num=tmpse.size();
                cols[(j<<1)+1]=real_num+cols[j<<1];
                cols[j+1<<1]=(real_num<<1)+cols[j<<1];
            } else {
                cols[(j<<1)+1]=cols[j<<1];
                cols[j+1<<1]=cols[j<<1];
            }            
        }
        uint64_t valNum=val.size();
        fwrite (&valNum,sizeof(uint64_t), 1, smr1);
        fwrite (&cols[0],sizeof(uint64_t), cols.size(), smr1);
        fwrite (&rowids[0],sizeof(uint32_t), rowids.size(), smr1);
        fwrite (&val[0],sizeof(float), val.size(), smr1);
        fclose (smr1);
        
        printf("The summary information of selection has been saved in %s.\n", logfname.c_str());
        cout<<"\nBeta values and SE values for "<<epiNum<<" Probes have been saved in the sparse binary file [" + esdfile + "]." <<endl;
        fclose(logfile);

    }
    void save_slct_besds_sbesd(char* outFileName, vector<probeinfolst2> &probeinfo,int cis_itvl,int trans_itvl,float transThres,float restThres,vector<string> &smasNames)
    {
       
        map<string, int> epi_map;
        map<string,int>::iterator iter;
        for(int j=0;j<probeinfo.size();j++)
        {
            epi_map.insert(pair<string,int>(probeinfo[j].probeId,j));
        }
        // get esd info
        long epiNum=probeinfo.size();
        string esdfile=string(outFileName)+string(".besd");
        FILE * smr1;
        smr1 = fopen (esdfile.c_str(), "wb");
        if (!(smr1)) {
            printf("ERROR: Failed to open file %s.\n",esdfile.c_str());
            exit(EXIT_FAILURE);
        }
        float filetype=SPARSE_FILE_TYPE_3;
        fwrite (&filetype,sizeof(float), 1, smr1);
       
        bool prtscr=false;
        
        // log file
        FILE* logfile=NULL;
        string logfname = string(outFileName)+".summary";
        logfile=fopen(logfname.c_str(), "w");
        if (!(logfile)) {
            printf("Error: Failed to open log file.\n");
            exit(1);
        }
        string logstr="cis-window:\t"+itos(cis_itvl)+"Kb\ntrans-window:\t"+itos(trans_itvl)+"Kb\np-value threshold of trans:\t"+dtos(transThres)+"\np-value threshold of others:\t"+dtos(restThres)+"\n";
        logstr+="\ncis region is indicated by [Chr, Start bp, End bp, nsnp];\ntrans region is indicated by <Chr, Start bp, End bp, nsnp>;\nthe number of other SNPs selected is indicated by (NumSNPs beyond cis and trans).\n";
        
        logstr+="\n{ProbeID, ProbeChr, ProbeBP}\t[Chr,cis_startBP,cis_endBP,NumSNPs]\t<Chr,trans_startBP,trans_endBP,NumSNPs>\t(NumSNPs)\n";
        fputs(logstr.c_str(),logfile);
        fflush(logfile);
        cis_itvl=cis_itvl*1000;
        trans_itvl=trans_itvl*1000;
        vector<snpinfolst> snpinfo;
        eqtlInfo eqtlinfo;
       
        read_esifile(&eqtlinfo, string(probeinfo[0].besdpath[0])+".esi", prtscr);
        
        for (int i = 0; i < smasNames.size(); i++)
        {
            printf("Reading... %3.0f%%\r", 100.0*i/(smasNames.size()));
            fflush(stdout);
            
            string epifile = smasNames[i]+".epi";
            read_epifile(&eqtlinfo, epifile, prtscr);
            
            string besdfile = smasNames[i]+".besd";
            FILE *fptr=fopen(besdfile.c_str(), "rb");
            if(!fptr)
            {
                printf ( "ERROR: Couldn't open file %s\n", besdfile.c_str());
                exit (EXIT_FAILURE);
            }
            float filetype=readfloat(fptr);
            uint64_t valNum=0;
            uint64_t* ptr=NULL;
            uint64_t rowSTART=0;
            uint64_t valSTART=0;
            if((int)filetype==SPARSE_FILE_TYPE_3 ){
                uint64_t colNum=(eqtlinfo._probNum<<1)+1;
                fseek(fptr, 0L, SEEK_END);
                uint64_t lSize = ftell(fptr);
                fseek(fptr, 0L, SEEK_SET);
                readfloat(fptr);
                valNum=readuint64(fptr);
                if( lSize - (sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float)) != 0) {fputs ("wrong element number",stderr); exit (3);}
                
                
                uint64_t colsize=colNum*sizeof(uint64_t);
                uint64_t* colbuf=(uint64_t*)malloc(colsize);
                if(NULL == colbuf)
                {
                    printf("ERROR: Can't malloc the reading cols buffer for %llu MB.\n",(colsize>>20));
                    exit(EXIT_FAILURE);
                }
                fread(colbuf,colNum,sizeof(uint64_t),fptr);
                
                ptr=colbuf;
                rowSTART=sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t);
                valSTART=sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t)+valNum*sizeof(uint32_t);
                
                for( int j=0;j<eqtlinfo._probNum;j++)
                {
                    snpinfo.clear();
                    string curprb=eqtlinfo._epi_prbID[j];
                    iter=epi_map.find(curprb);
                    int prbindx=iter->second;
                    if(probeinfo[prbindx].besdpath.size()>1) // or if(probeinfo[prbindx].beta_se.size()>0)
                    {
                        printf("ERROR: Probe %s is found in %ld BESD files. please check wether the SNPs are in consistency among .esi files, if yes, please remove one of the duplicate probe then have a try, if no, please disable --geno-uni, then have a try.\n",curprb.c_str(),probeinfo[prbindx].besdpath.size());
                        exit(EXIT_FAILURE);
                    }
                    vector<uint32_t> rowidx;
                    uint64_t pos=*(ptr+(j<<1)); //BETA START
                    uint64_t pos1=*(ptr+(j<<1)+1); //SE START
                    uint64_t num=pos1-pos;
                    if(num>0)
                    {
                        uint32_t* ridbuff=(uint32_t*)malloc(num*2*sizeof(uint32_t));
                        if(NULL == ridbuff)
                        {
                            printf("ERROR: Memory allocation error when when reading probe %s in the file %s.\n",curprb.c_str(),besdfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                        memset(ridbuff,0,num*2*sizeof(uint32_t));
                        float* betasebuff=(float*)malloc(num*2*sizeof(float));
                        if(NULL == betasebuff)
                        {
                            printf("ERROR: Memory allocation error when when reading probe %s in the file %s.\n",curprb.c_str(),besdfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                        memset(betasebuff,0,num*2*sizeof(float));
                        
                        fseek(fptr, rowSTART+pos*sizeof(uint32_t), SEEK_SET);
                        fread(ridbuff, sizeof(uint32_t),2*num,fptr);
                        fseek(fptr,valSTART+pos*sizeof(float),SEEK_SET);
                        fread(betasebuff,sizeof(float), 2*num,fptr);
                        for(int jj=0;jj<num;jj++)
                        {
                            double beta=betasebuff[jj];
                            double se=betasebuff[jj+num];
                            if(ABS(se+9)<1e-6) continue;
                            int rowid=ridbuff[jj];
                            snpinfolst tmpinfo;
                            tmpinfo.beta=beta;
                            tmpinfo.se=se;
                            tmpinfo.snpchr=eqtlinfo._esi_chr[rowid];
                            tmpinfo.snprs=eqtlinfo._esi_rs[rowid];
                            tmpinfo.a1=eqtlinfo._esi_allele1[rowid];
                            tmpinfo.a2=eqtlinfo._esi_allele2[rowid];
                            tmpinfo.bp=eqtlinfo._esi_bp[rowid];
                            snpinfo.push_back(tmpinfo);
                            rowidx.push_back(rowid);
                        }
                        free(ridbuff);
                        free(betasebuff);
                        
                    } else {
                        //printf("Probe %s has no values.\n",curprb.c_str());
                    }
                    
                    if(snpinfo.size()>0)
                    {
                        snpinfolst* sortptr=&snpinfo[0];
                        qsort(sortptr,snpinfo.size(),sizeof(snpinfolst),comp_esi); // when reading sparse, the rowids in file is not in order, so snpinfo is not in order, so sort should be implement here.
                        
                        probeinfolst prbifo;
                        prbifo.bp=eqtlinfo._epi_bp[j];// probeinfo[prbindx].bp;
                        prbifo.probechr=eqtlinfo._epi_chr[j];//probeinfo[prbindx].probechr;
                        prbifo.probeId=curprb;//probeinfo[prbindx].probeId;
                        
                        vector<int> slct_idx;
                        slct_sparse_per_prb(slct_idx, &prbifo, snpinfo,  cis_itvl,  trans_itvl, transThres, restThres,logfile); //slct_idx with no order if there are trans-rgeions
                        stable_sort(slct_idx.begin(),slct_idx.end());
                        
                        uint32_t* ridbuff=(uint32_t*)malloc(slct_idx.size()*2*sizeof(uint32_t));
                        if(NULL == ridbuff)
                        {
                            printf("ERROR: Memory allocation error when when dealing with probe %s in the file %s.\n",curprb.c_str(),besdfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                        memset(ridbuff,0,slct_idx.size()*2*sizeof(uint32_t));
                        float* betasebuff=(float*)malloc(slct_idx.size()*2*sizeof(float));
                        if(NULL == betasebuff)
                        {
                            printf("ERROR: Memory allocation error when when dealing with probe %s in the file %s.\n",curprb.c_str(),besdfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                        memset(betasebuff,0,slct_idx.size()*2*sizeof(float));
                        /*
                         // have relaced this with rowidx
                        vector<string> _rs(slct_idx.size());
                        for(int l=0;l<slct_idx.size();l++) _rs[l]=snpinfo[slct_idx[l]].snprs;
                        vector<uint32_t> rsid(_rs.size());
                        #pragma omp parallel for private(iter)
                        for (int l = 0; l<_rs.size(); l++){
                            iter = esi_map.find(_rs[l]);
                            if (iter != esi_map.end()) rsid[l]=iter->second;
                            else {
                                printf("ERROR: SNP is not in SNP map. if you are using --geno-uni, please disable it then try again, else please report this bug.");
                                exit(EXIT_FAILURE);
                            }
                        }
                         */

                        for(int l=0;l<slct_idx.size();l++)
                        {
                            //ridbuff[l]=rsid[l];
                            ridbuff[l]=rowidx[slct_idx[l]];
                            betasebuff[l]=snpinfo[slct_idx[l]].beta;
                        }
                        
                        for(int l=0;l<slct_idx.size();l++)
                        {
                            //ridbuff[l+slct_idx.size()]=rsid[l];
                            ridbuff[l+slct_idx.size()]=rowidx[slct_idx[l]];
                            betasebuff[l+slct_idx.size()]=snpinfo[slct_idx[l]].se;
                        }
                        probeinfo[prbindx].vnum=slct_idx.size();
                        probeinfo[prbindx].rowid=ridbuff;
                        probeinfo[prbindx].beta_se=betasebuff;
                    }
                    
                }
                
                free(colbuf);
            }
            else if((int)filetype==DENSE_FILE_TYPE_1) {
                
                float* tmpbetase=(float*)malloc(sizeof(float)*eqtlinfo._snpNum<<1);
                if(NULL == tmpbetase)
                {
                    printf("ERROR: Can't malloc the reading cols buffer for %u MB.\n",(eqtlinfo._snpNum>>17));
                    exit(EXIT_FAILURE);
                }
                
                for( int j=0;j<eqtlinfo._probNum;j++)
                {
                    snpinfo.clear();
                    vector<uint32_t> rowidx;
                    string curprb=eqtlinfo._epi_prbID[j];
                    iter=epi_map.find(curprb);
                    int prbindx=iter->second;
                    if(probeinfo[prbindx].besdpath.size()>1)
                    {
                        printf("ERROR: Probe %s is found in %ld BESD files. please check wether the SNPs are in consistency among .esi files, if yes, please remove one of the duplicate probe then have a try, if no, please disable --geno-uni, then have a try.\n",curprb.c_str(),probeinfo[prbindx].besdpath.size());
                        exit(EXIT_FAILURE);
                    }
                    memset(tmpbetase,0,sizeof(float)*eqtlinfo._snpNum<<1);
                    fseek(fptr,((j<<1)*eqtlinfo._snpNum+1)<<2, SEEK_SET);
                    fread(tmpbetase, sizeof(float),eqtlinfo._snpNum<<1,fptr);
                    for(uint32_t jj=0;jj<eqtlinfo._snpNum;jj++)
                    {
                        float beta=tmpbetase[jj];
                        float se=tmpbetase[jj+eqtlinfo._snpNum];
                        if(ABS(se+9)<1e-6) continue;
                       
                        snpinfolst tmpinfo;
                        tmpinfo.beta=beta;
                        tmpinfo.se=se;
                        tmpinfo.snpchr=eqtlinfo._esi_chr[jj];
                        tmpinfo.snprs=eqtlinfo._esi_rs[jj];
                        tmpinfo.a1=eqtlinfo._esi_allele1[jj];
                        tmpinfo.a2=eqtlinfo._esi_allele2[jj];
                        tmpinfo.bp=eqtlinfo._esi_bp[jj];
                        snpinfo.push_back(tmpinfo);
                        rowidx.push_back(jj);
                    }
                    
                    if(snpinfo.size()>0)
                    {
                        snpinfolst* sortptr=&snpinfo[0];
                        qsort(sortptr,snpinfo.size(),sizeof(snpinfolst),comp_esi); // when reading sparse, the rowids in file is not in order, so snpinfo is not in order, so sort should be implement here.
                        
                        probeinfolst prbifo;
                        prbifo.bp=eqtlinfo._epi_bp[j];// probeinfo[prbindx].bp;
                        prbifo.probechr=eqtlinfo._epi_chr[j];//probeinfo[prbindx].probechr;
                        prbifo.probeId=curprb;//probeinfo[prbindx].probeId;
                        
                        vector<int> slct_idx;
                        slct_sparse_per_prb(slct_idx, &prbifo, snpinfo,  cis_itvl,  trans_itvl, transThres, restThres,logfile); //slct_idx with no order if there are trans-rgeions
                        stable_sort(slct_idx.begin(),slct_idx.end());
                        
                        uint32_t* ridbuff=(uint32_t*)malloc(slct_idx.size()*2*sizeof(uint32_t));
                        if(NULL == ridbuff)
                        {
                            printf("ERROR: Memory allocation error when when dealing with probe %s in the file %s.\n",curprb.c_str(),besdfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                        memset(ridbuff,0,slct_idx.size()*2*sizeof(uint32_t));
                        float* betasebuff=(float*)malloc(slct_idx.size()*2*sizeof(float));
                        if(NULL == betasebuff)
                        {
                            printf("ERROR: Memory allocation error when when dealing with probe %s in the file %s.\n",curprb.c_str(),besdfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                        memset(betasebuff,0,slct_idx.size()*2*sizeof(float));
                        
                        for(int l=0;l<slct_idx.size();l++)
                        {
                            //ridbuff[l]=rsid[l];
                            ridbuff[l]=rowidx[slct_idx[l]];
                            betasebuff[l]=snpinfo[slct_idx[l]].beta;
                        }
                        
                        for(int l=0;l<slct_idx.size();l++)
                        {
                            //ridbuff[l+slct_idx.size()]=rsid[l];
                            ridbuff[l+slct_idx.size()]=rowidx[slct_idx[l]];
                            betasebuff[l+slct_idx.size()]=snpinfo[slct_idx[l]].se;
                        }
                        probeinfo[prbindx].vnum=slct_idx.size();
                        probeinfo[prbindx].rowid=ridbuff;
                        probeinfo[prbindx].beta_se=betasebuff;
                    }
                }
                free(tmpbetase);
            }
            else {
                printf("Your file is in the old sparse format. please first re-make it by --beqtl-summary and --make-besd.\n");
                exit(EXIT_FAILURE);
            }
            fclose(fptr);
            
        }
        
        vector<uint64_t> cols((epiNum<<1)+1);
        uint64_t valNum=0;
        cols[0]=0;
        
        for(int j=0;j<epiNum;j++)
        {
            uint64_t real_num=probeinfo[j].vnum;
            cols[(j<<1)+1]=real_num+cols[j<<1];
            cols[j+1<<1]=(real_num<<1)+cols[j<<1];
            valNum+=real_num*2;
        }
        
        fwrite (&valNum,sizeof(uint64_t), 1, smr1);
        fwrite (&cols[0],sizeof(uint64_t), cols.size(), smr1);
        
        for(int j=0;j<epiNum;j++)
        {
            printf("Saving...  %3.0f%%\r", 100.0*j/(2*epiNum));
            fflush(stdout);
            fwrite (probeinfo[j].rowid,sizeof(uint32_t), probeinfo[j].vnum*2, smr1);
        }
        
        for(int j=0;j<epiNum;j++)
        {
            printf("Saving...  %3.0f%%\r", (100.0*j/(2*epiNum)+50));
            fflush(stdout);
            fwrite (probeinfo[j].beta_se,sizeof(float), probeinfo[j].vnum*2, smr1);
        }
        fclose (smr1);
        for(int j=0;j<epiNum;j++) {
            free(probeinfo[j].beta_se);
            probeinfo[j].beta_se=NULL;
            free(probeinfo[j].rowid);
            probeinfo[j].rowid=NULL;
        }
        printf("The summary information of selection has been saved in %s.\n", logfname.c_str());
        cout<<"\nBeta values and SE values for "<<epiNum<<" Probes have been saved in the sparse binary file [" + esdfile + "]." <<endl;
        fclose(logfile);
        
    }
    void combineBesd(char* eqtlsmaslstName, char* outFileName,bool save_dense_flag, int cis_itvl, int trans_itvl, float transThres, float restThres, bool genouni)
    {
        vector<string> smasNames;
        vector<snpinfolst> snpinfo;
        vector<probeinfolst2> probeinfo;
        if(genouni) {
            printf("WARNING: --geno-uni is actived. Please ensure the SNPs in all the .esi files are the same and the allele pairs <the effect allele, the other allele> of each SNP in all the .esi files are the same!\n");
        }
        read_smaslist(smasNames, string(eqtlsmaslstName));
        if(smasNames.size()==0) throw("No eqtl summary file list in [ "+ string(eqtlsmaslstName)  +" ]");
        combine_esi(snpinfo, smasNames, genouni);
        if(snpinfo.size()==0)
        {
            printf("ERROR: No SNPs to be included!");
            exit(EXIT_FAILURE);
        }
        snpinfolst* esiptr=&snpinfo[0];
        qsort(esiptr,snpinfo.size(),sizeof(snpinfolst),comp_esi);
        combine_epi(probeinfo, smasNames);
        if(probeinfo.size()==0)
        {
            printf("ERROR: No probes to be included!");
            exit(EXIT_FAILURE);
        }
        probeinfolst2* epiptr=&probeinfo[0];
        qsort(epiptr,probeinfo.size(),sizeof(probeinfolst2),comp2);
        
        printf("\nGenerating epi file...\n");       
        string epifile = string(outFileName)+string(".epi");
        ofstream epi(epifile.c_str());
        if (!epi) throw ("Error: can not open the epi file " + epifile + " to save!");
                for (int j = 0;j <probeinfo.size(); j++) {
            epi<<probeinfo[j].probechr<<'\t'<<probeinfo[j].probeId<<'\t'<<0<<'\t'<<probeinfo[+j].bp<<'\t'<<probeinfo[j].genename<<'\t'<<probeinfo[j].orien<<'\n';
        }
        epi.close();
        printf("%ld probes have been saved in the file %s.\n",probeinfo.size(),epifile.c_str());
        
        printf("\nGenerating esi file...\n");
        vector<string> esi_rs,esi_a1,esi_a2;
        if(!genouni)
        {
            esi_rs.resize(snpinfo.size());
            esi_a1.resize(snpinfo.size());
            esi_a2.resize(snpinfo.size());
        }
        string esifile =  string(outFileName)+string(".esi");
        ofstream esi(esifile.c_str());
        if (!esi) throw ("Error: can not open the esi file to save!");
        for (long j = 0;j <snpinfo.size(); j++) {
            esi<<snpinfo[j].snpchr<<'\t'<<snpinfo[j].snprs<<'\t'<<snpinfo[j].gd<<'\t'<<snpinfo[j].bp<<'\t'<<snpinfo[j].a1<<'\t'<<snpinfo[j].a2<<'\t'<<snpinfo[j].freq<<'\n';
            if(!genouni)
            {
                esi_rs[j]=snpinfo[j].snprs;
                esi_a1[j]=snpinfo[j].a1;
                esi_a2[j]=snpinfo[j].a2;
            }
            
        }
        esi.close();
        printf("%ld SNPs have been saved in the file %s.\n",snpinfo.size(),esifile.c_str());
        
        uint64_t valnum=countNotNullNum(smasNames);
        double sparsity=1.0*valnum/(probeinfo.size()*snpinfo.size());
        if(save_dense_flag)
        {
            if(sparsity>=0.4)
            {
                printf("The sparsity of your data matrix is %f. We are going to save it in dense format!\n", sparsity);
                if(genouni) save_besds_dbesd(outFileName, snpinfo, probeinfo);
                else save_besds_dbesd(outFileName, snpinfo, probeinfo,esi_rs,esi_a1,esi_a2,genouni);
                
            } else {
                printf("The sparsity of your data matrix is %f. We are going to save it in sparse format!\n", sparsity);
                if(genouni) save_besds_sbesd( outFileName, snpinfo, probeinfo,smasNames);
                else save_besds_sbesd( outFileName, snpinfo, probeinfo,esi_rs,esi_a1,esi_a2,genouni);
            }
        } else
        {
             long esiNum=snpinfo.size();
            if(genouni) save_slct_besds_sbesd(outFileName, probeinfo, cis_itvl,  trans_itvl,  transThres,  restThres,smasNames);
            else save_slct_besds_sbesd(outFileName, esiNum, probeinfo,esi_rs,esi_a1,esi_a2, cis_itvl,  trans_itvl,  transThres,  restThres,genouni);
        }
    }
    float est_sample_size(float freq, float beta, float se)
    {
        return (1.0-2*freq*(1-freq)*beta*beta)/(2*freq*(1-freq)*se*se);
    }
    void get_snpinfo_cur_prb_sparse(vector<snpinfolst> &snpinfo,FILE* fptr,  uint64_t pid, uint64_t* ptr, uint64_t rowSTART,uint64_t valSTART,eqtlInfo* etmp,map<int, int> &_incld_id_map, bool qcflag){
        
        uint64_t pos=*(ptr+(pid<<1)); //BETA START
        uint64_t pos1=*(ptr+(pid<<1)+1); //SE START
        uint64_t num=pos1-pos;
        uint64_t real_num=0;
        bool nufreqwarnflg=false;
        char* row_char_ptr;
        row_char_ptr = (char*) malloc (sizeof(char)*2*num*sizeof(uint32_t));
        if (row_char_ptr == NULL) {fputs ("Memory error",stderr); exit (1);}
        char* val_char_ptr;
        val_char_ptr = (char*) malloc (sizeof(char)*2*num*sizeof(float));
        if (val_char_ptr == NULL) {fputs ("Memory error",stderr); exit (1);}
        memset(row_char_ptr,0,sizeof(char)*2*num*sizeof(uint32_t));
        memset(val_char_ptr,0,sizeof(char)*2*num*sizeof(float));
        fseek(fptr, rowSTART+pos*sizeof(uint32_t), SEEK_SET);
        fread(row_char_ptr, sizeof(uint32_t),2*num,fptr);
        
        uint32_t* row_ptr=(uint32_t *)row_char_ptr;
        fseek(fptr,valSTART+pos*sizeof(float),SEEK_SET);
        fread(val_char_ptr,sizeof(float), 2*num,fptr);
        float* val_ptr=(float*)val_char_ptr;
        for(int j=0;j<num;j++)
        {
            snpinfolst snpinfotmp;
            uint32_t rid=*(row_ptr+j);
            
            map<int, int>::iterator iter;
            iter=_incld_id_map.find(rid);
            if(iter!=_incld_id_map.end())
            {
                snpinfotmp.snprs=etmp->_esi_rs[rid];
                snpinfotmp.snpchr=etmp->_esi_chr[rid];
                snpinfotmp.bp=etmp->_esi_bp[rid];
                snpinfotmp.gd=etmp->_esi_gd[rid];
                snpinfotmp.a1=etmp->_esi_allele1[rid];
                snpinfotmp.a2=etmp->_esi_allele2[rid];
                snpinfotmp.freq=etmp->_esi_freq[rid];
                snpinfotmp.beta=*(val_ptr+j);
                snpinfotmp.se=*(val_ptr+j+num);
                if(qcflag){
                    if(abs(snpinfotmp.freq+9)<1e-6 && !nufreqwarnflg)
                    {
                        printf("WARNING: one or more NA freqencies found. This SNP would be excluded.\n");
                        nufreqwarnflg=true;
                        continue;
                    }
                    if(snpinfotmp.freq<1e-8)
                    {
                        printf("WARNING: %s freqency is 0. This SNP would be excluded.\n",snpinfotmp.snprs.c_str());
                        continue;
                    }
                    snpinfotmp.estn=est_sample_size(snpinfotmp.freq, snpinfotmp.beta, snpinfotmp.se);
                    if(snpinfotmp.estn<0)
                    {
                        printf("ERROR: Negative estimated sample size found of SNP %s.\n",snpinfotmp.snprs.c_str());
                        exit(EXIT_FAILURE);
                    }
                    snpinfo.push_back(snpinfotmp);
                } else {
                    snpinfotmp.estn=-9;
                    snpinfo.push_back(snpinfotmp);
                }
                
                //int sid=iter->second;
                //cout<<rid<<":"<<etmp._esi_include[sid]<<endl; // test passed
                real_num++;
            }
            
        }
        free(row_char_ptr);
        free(val_char_ptr);

    }
    void get_snpinfo_cur_prb_dense(vector<snpinfolst> &snpinfo,FILE* fptr,  uint64_t pid, char** buffer, eqtlInfo* etmp)
    {
        fseek(fptr,((pid<<1)*etmp->_snpNum+1)<<2, SEEK_SET);
        
        memset(*buffer,0,sizeof(char)*etmp->_snpNum<<3);
        fread(*buffer, sizeof(char),etmp->_snpNum<<3,fptr);
        float* ft=(float *)*buffer;
        float* se_ptr = ft + etmp->_snpNum;
        for (int j = 0; j<etmp->_esi_include.size(); j++) {
            float se=*(se_ptr + etmp->_esi_include[j]);
            if(abs(se+9)>1e-6)
            {
                snpinfolst snpinfotmp;
                snpinfotmp.snprs=etmp->_esi_rs[etmp->_esi_include[j]];
                snpinfotmp.snpchr=etmp->_esi_chr[etmp->_esi_include[j]];
                snpinfotmp.bp=etmp->_esi_bp[etmp->_esi_include[j]];
                snpinfotmp.gd=etmp->_esi_gd[etmp->_esi_include[j]];
                snpinfotmp.a1=etmp->_esi_allele1[etmp->_esi_include[j]];
                snpinfotmp.a2=etmp->_esi_allele2[etmp->_esi_include[j]];
                snpinfotmp.freq=etmp->_esi_freq[etmp->_esi_include[j]];
                snpinfotmp.beta=*(ft + etmp->_esi_include[j]);
                snpinfotmp.se=se;
                snpinfo.push_back(snpinfotmp);
            }
        }
        

    }
    
    void qc(vector<snpinfolst> &snpinfo,probeinfolst* prbifo,int qc_mtd,int z_thresh, vector<float> &suminfo,FILE* outlierfptr)
    {
        snpinfolst* sortptr=&snpinfo[0];
        qsort(sortptr,snpinfo.size(),sizeof(snpinfolst),comp_estn);
        long ttlnum=snpinfo.size();
        float minn=snpinfo[0].estn;
        float maxn=snpinfo[ttlnum-1].estn;
        float mediann=snpinfo[ttlnum/2].estn;
        float sumn=0;
        for(int i=0;i<ttlnum;i++)
        {
            sumn+=snpinfo[i].estn;
        }
        float avgn=sumn/ttlnum;
        float sdn=0;
        for(int i=0;i<ttlnum;i++)
        {
            sdn+=(snpinfo[i].estn-avgn)*(snpinfo[i].estn-avgn);
        }
        sdn=sqrt(sdn/(ttlnum-1));
        suminfo.push_back(maxn);
        suminfo.push_back(minn);
        suminfo.push_back(avgn);
        suminfo.push_back(mediann);
        suminfo.push_back(sdn);
        suminfo.push_back(ttlnum);
        long lowidx=-1;
        long upidx=ttlnum;
        if(qc_mtd==0)
        {
            float cutoff_ratio=0.0005;
            long cutnum=ceil(cutoff_ratio*ttlnum);
            lowidx=cutnum;
            upidx=ttlnum-cutnum+1;
           
        }
        else if(qc_mtd==1)
        {
            for(long i=0;i<ttlnum;i++)
            {
                float z=(snpinfo[i].estn-avgn)/sdn;
                if(abs(z)>3) lowidx=i;
                else break;
            }
            for(long i=ttlnum-1;i>=0;i--)
            {
                float z=(snpinfo[i].estn-avgn)/sdn;
                if(abs(z)>z_thresh) upidx=i;
                else break;
            }
            
        }
        else if(qc_mtd==2)
        {
            printf("ERROR: MAD method is comming soon.\n");
            exit(EXIT_FAILURE);
        }
        if(lowidx>-1)
        {
            for(long i=0;i<=lowidx;i++)
            {
                string logstr=snpinfo[i].snprs+'\t'+atos(snpinfo[i].snpchr)+'\t'+atos(snpinfo[i].bp)+'\t'+snpinfo[i].a1+'\t'+snpinfo[i].a2+'\t'+atos(snpinfo[i].freq)+'\t'+'\t'+prbifo->probeId+'\t'+atos(prbifo->probechr)+'\t'+atos(prbifo->bp)+'\t'+prbifo->genename+'\t'+prbifo->orien+'\t'+atos(snpinfo[i].beta)+'\t'+atos(snpinfo[i].se)+'\t'+atos(snpinfo[i].estn)+'\n';
                fputs(logstr.c_str(),outlierfptr);
                fflush(outlierfptr);
            }
        }
        if(upidx<ttlnum)
        {
            for(long i=ttlnum-1;i>=upidx;i--)
            {
                string logstr=snpinfo[i].snprs+'\t'+atos(snpinfo[i].snpchr)+'\t'+atos(snpinfo[i].bp)+'\t'+snpinfo[i].a1+'\t'+snpinfo[i].a2+'\t'+atos(snpinfo[i].freq)+'\t'+'\t'+prbifo->probeId+'\t'+atos(prbifo->probechr)+'\t'+atos(prbifo->bp)+'\t'+prbifo->genename+'\t'+prbifo->orien+'\t'+atos(snpinfo[i].beta)+'\t'+atos(snpinfo[i].se)+'\t'+atos(snpinfo[i].estn)+'\n';
                fputs(logstr.c_str(),outlierfptr);
                fflush(outlierfptr);
            }
        }
        long rmnum=0;
        if(lowidx>-1){
            rmnum=lowidx+1;
            snpinfo.erase(snpinfo.begin(), snpinfo.begin()+rmnum);
        }
        if(upidx<ttlnum){
            rmnum+=ttlnum-upidx;
            snpinfo.erase(snpinfo.end()-ttlnum+upidx, snpinfo.end());
        }
        suminfo.push_back(rmnum);
    }
    
    void make_sparse_besd(char* eqtlFileName, char* outFileName, int cis_itvl, int trans_itvl, float transThres, float restThres,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* snplstName,char* problstName, char* snplst2exclde, char* problst2exclde, bool qcflag, int qc_mtd, int z_thresh)
    {
        
        eqtlInfo etmp;
        read_esifile(&etmp, string(eqtlFileName)+".esi");
        esi_man(&etmp, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,  cis_itvl, prbname);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&etmp, snplst2exclde);
        read_epifile(&etmp, string(eqtlFileName)+".epi");
        epi_man(&etmp, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
        if(problst2exclde != NULL) exclude_prob(&etmp, problst2exclde);
        

        printf("\nGenerating epi file...\n");
        string epifile = string(outFileName)+string(".epi");
        ofstream epi(epifile.c_str());
        if (!epi) throw ("Error: can not open the epi file " + epifile + " to save!");
        for (int j = 0;j <etmp._include.size(); j++) {
            epi<<etmp._epi_chr[etmp._include[j]]<<'\t'<<etmp._epi_prbID[etmp._include[j]]<<'\t'<<etmp._epi_gd[etmp._include[j]]<<'\t'<<etmp._epi_bp[etmp._include[j]]<<'\t'<<etmp._epi_gene[etmp._include[j]]<<'\t'<<etmp._epi_orien[etmp._include[j]]<<'\n';
        }
        epi.close();
        printf("%ld probes have been saved in the file %s.\n",etmp._include.size(),epifile.c_str());
        
        printf("\nGenerating esi file...\n");
        map<string, int> esi_map;
        vector<string> esi_rs(etmp._esi_include.size());
        vector<string> esi_a1(etmp._esi_include.size());
        vector<string> esi_a2(etmp._esi_include.size());
        string esifile =  string(outFileName)+string(".esi");
        ofstream esi(esifile.c_str());
        if (!esi) throw ("Error: can not open the esi file to save!");
        for (long j = 0;j <etmp._esi_include.size(); j++) {
            esi<< etmp._esi_chr[etmp._esi_include[j]]<<'\t'<<etmp._esi_rs[etmp._esi_include[j]]<<'\t'<<etmp._esi_gd[etmp._esi_include[j]]<<'\t'<<etmp._esi_bp[etmp._esi_include[j]]<<'\t'<<etmp._esi_allele1[etmp._esi_include[j]]<<'\t'<<etmp._esi_allele2[etmp._esi_include[j]]<<'\t'<<etmp._esi_freq[etmp._esi_include[j]]<<'\n';
            esi_map.insert(pair<string,int>(etmp._esi_rs[etmp._esi_include[j]],j));
            esi_rs[j]=etmp._esi_rs[etmp._esi_include[j]];
            esi_a1[j]=etmp._esi_allele1[etmp._esi_include[j]];
            esi_a2[j]=etmp._esi_allele2[etmp._esi_include[j]];
        }
        esi.close();
        printf("%ld SNPs have been saved in the file %s.\n",etmp._esi_include.size(),esifile.c_str());
        
        printf("\nGenerating besd file...\n");
        string esdfile=string(outFileName)+string(".besd");
        FILE * smr1;
        smr1 = fopen (esdfile.c_str(), "wb");
        if (!(smr1)) {
            printf("ERROR: Failed to open file %s.\n",esdfile.c_str());
            exit(EXIT_FAILURE);
        }
        float filetype2write=SPARSE_FILE_TYPE_3;
        fwrite (&filetype2write,sizeof(float), 1, smr1);
        
        vector<uint64_t> cols((etmp._include.size()<<1)+1);;
        vector<uint32_t> rowids;
        vector<float> val;
        cols[0]=0;

        
        string besdfile = string(eqtlFileName)+".besd";
        FILE *fptr=fopen(besdfile.c_str(), "rb");
        if(!fptr)
        {
            printf ( "ERROR: Couldn't open file %s\n", besdfile.c_str());
            exit (EXIT_FAILURE);
        }
        
        map<int, int > _incld_id_map;
        long size = 0;
        for (int i = 0; i<etmp._esi_include.size(); i++)
        {
            _incld_id_map.insert(pair<int, int>(etmp._esi_include[i], i));
            if (size == _incld_id_map.size()) throw ("Error: Duplicated SNP IDs found: \"" + etmp._esi_rs[etmp._esi_include[i]] + "\".");
            size = _incld_id_map.size();
        }

        float filetype=readfloat(fptr);
        char* buffer=NULL;
        uint64_t valNum=0;
        uint64_t* ptr=NULL;
        uint64_t rowSTART=0;
        uint64_t valSTART=0;
        if((int)filetype==SPARSE_FILE_TYPE_3 ){
            uint64_t colNum=(etmp._probNum<<1)+1;
            fseek(fptr, 0L, SEEK_END);
            uint64_t lSize = ftell(fptr);
            fseek(fptr, 0L, SEEK_SET);
            readfloat(fptr);
            valNum=readuint64(fptr);
            if( lSize - (sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float)) != 0) {fputs ("wrong element number",stderr); exit (3);}
            uint64_t colsize=colNum*sizeof(uint64_t);
            buffer = (char*) malloc (sizeof(char)*(colsize));
            if (buffer == NULL) {fputs ("Memory error when reading sparse BESD file.",stderr); exit (1);}
            fread(buffer,colsize,sizeof(char),fptr);
            
            ptr=(uint64_t *)buffer;
            rowSTART=sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t);
            valSTART=sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t)+valNum*sizeof(uint32_t);
        } else if((int)filetype==DENSE_FILE_TYPE_1){
            buffer = (char*) malloc (sizeof(char)*etmp._snpNum<<3);
            if (buffer == NULL) {fputs ("Memory error when reading dense BESD file.",stderr); exit (1);}
        } else {
            printf("Your file is in the old sparse format. please first re-make it by --beqtl-summary and --make-besd.\n");
            exit(EXIT_FAILURE);
        }
        
        // log file
        FILE* logfile=NULL;
        string logfname = string(outFileName)+".summary";
        logfile=fopen(logfname.c_str(), "w");
        if (!(logfile)) {
            printf("Error: Failed to open log file.\n");
            exit(1);
        }
        string logstr="cis-window:\t"+itos(cis_itvl)+"Kb\ntrans-window:\t"+itos(trans_itvl)+"Kb\np-value threshold of trans:\t"+dtos(transThres)+"\np-value threshold of others:\t"+dtos(restThres)+"\n";
        logstr+="\ncis region is indicated by [Chr, Start bp, End bp, nsnp];\ntrans region is indicated by <Chr, Start bp, End bp, nsnp>;\nthe number of other SNPs selected is indicated by (NumSNPs beyond cis and trans).\n";
        
        logstr+="\n{ProbeID, ProbeChr, ProbeBP}\t[Chr,cis_startBP,cis_endBP,NumSNPs]\t<Chr,trans_startBP,trans_endBP,NumSNPs>\t(NumSNPs)\n";
        fputs(logstr.c_str(),logfile);
        fflush(logfile);
        cis_itvl=cis_itvl*1000;
        trans_itvl=trans_itvl*1000;
        vector<snpinfolst> snpinfo;
        map<string, int>::iterator iter;
        FILE* outlierfptr=NULL;
        FILE* qcsryfptr=NULL;
        if(qcflag){
            // outlier file
            string outlierfname = string(outFileName)+".outlier";
            outlierfptr=fopen(outlierfname.c_str(), "w");
            if (!(outlierfptr)) {
                printf("Error: Failed to open log file.\n");
                exit(1);
            }
            logstr="SNP\tChr\tBP\tA1\tA2\tFreq\tProbe\tProbe_Chr\tProbe_bp\tGene\tOrientation\tb\tse\test_n\n";
            fputs(logstr.c_str(),outlierfptr);
            fflush(outlierfptr);
            // qc summary file
            string qcsryfname = string(outFileName)+".qc.summary";
            qcsryfptr=fopen(qcsryfname.c_str(), "w");
            if (!(qcsryfptr)) {
                printf("Error: Failed to open log file.\n");
                exit(1);
            }
            logstr="ProbeID\tmaxn\tminn\tavgn\tmediann\tsdn\tttlnum\trmnum\n";
            fputs(logstr.c_str(),qcsryfptr);
            fflush(qcsryfptr);
        }
        
        for(int i=0;i<etmp._include.size();i++)
        {
            printf("Saving... %3.0f%%\r", 100.0*i/etmp._include.size());
            fflush(stdout);
            bool nufreqwarnflg=false;
            string prbname=etmp._epi_prbID[etmp._include[i]];
            vector<uint32_t> tmprid;
            vector<float> tmpse;
            snpinfo.clear();
            uint64_t pid=etmp._include[i];
            if((int)filetype==SPARSE_FILE_TYPE_3) get_snpinfo_cur_prb_sparse(snpinfo,fptr, pid, ptr,  rowSTART, valSTART, &etmp,_incld_id_map, qcflag);
            else {
                //get_snpinfo_cur_prb_dense(snpinfo,fptr, pid, &buffer ,&etmp);
                fseek(fptr,((pid<<1)*etmp._snpNum+1)<<2, SEEK_SET);
                memset(buffer,0,sizeof(char)*etmp._snpNum<<3);
                fread(buffer, sizeof(char),etmp._snpNum<<3,fptr);
                float* ft=(float *)buffer;
                float* se_ptr = ft + etmp._snpNum;
                for (int j = 0; j<etmp._esi_include.size(); j++) {
                    float se=*(se_ptr + etmp._esi_include[j]);
                    if(abs(se+9)>1e-6)
                    {
                        snpinfolst snpinfotmp;
                        snpinfotmp.snprs=etmp._esi_rs[etmp._esi_include[j]];
                        snpinfotmp.snpchr=etmp._esi_chr[etmp._esi_include[j]];
                        snpinfotmp.bp=etmp._esi_bp[etmp._esi_include[j]];
                        snpinfotmp.gd=etmp._esi_gd[etmp._esi_include[j]];
                        snpinfotmp.a1=etmp._esi_allele1[etmp._esi_include[j]];
                        snpinfotmp.a2=etmp._esi_allele2[etmp._esi_include[j]];
                        snpinfotmp.freq=etmp._esi_freq[etmp._esi_include[j]];
                        snpinfotmp.beta=*(ft + etmp._esi_include[j]);
                        snpinfotmp.se=se;
                        if(qcflag) {
                            if(abs(snpinfotmp.freq+9)<1e-6 && !nufreqwarnflg)
                            {
                                printf("WARNING: one or more NA freqencies found. This SNP would be excluded.\n");
                                nufreqwarnflg=true;
                                continue;
                            }
                            if(abs(snpinfotmp.freq)<1e-8)
                            {
                                printf("WARNING: %s freqency is 0. This SNP would be excluded.\n",snpinfotmp.snprs.c_str());
                                continue;
                            }
                            snpinfotmp.estn=est_sample_size(snpinfotmp.freq, snpinfotmp.beta, se);
                            if(snpinfotmp.estn<0)
                            {
                                printf("ERROR: Negative estimated sample size found of SNP %s.\n",snpinfotmp.snprs.c_str());
                                exit(EXIT_FAILURE);
                            }
                            snpinfo.push_back(snpinfotmp);
                        } else {
                            snpinfotmp.estn=-9;
                            snpinfo.push_back(snpinfotmp);
                        }
                    }
                }
            }
            
            probeinfolst prbifo;
            prbifo.bp=etmp._epi_bp[etmp._include[i]];
            prbifo.probechr=etmp._epi_chr[etmp._include[i]];
            prbifo.probeId=etmp._epi_prbID[etmp._include[i]];
            prbifo.genename=etmp._epi_gene[etmp._include[i]];
            prbifo.orien=etmp._epi_orien[etmp._include[i]];
            
            //  QC
            if(qcflag && snpinfo.size()>0) {
                vector<float> suminfo; // maxn,minn,avgn,mediann,sdn,ttnum,rmnum
                qc(snpinfo,&prbifo,qc_mtd,z_thresh,suminfo,outlierfptr);
                logstr=prbname;
                for(int j=0;j<suminfo.size();j++)
                {
                    logstr+='\t'+atos(suminfo[j]);
                }
                logstr+='\n';
                fputs(logstr.c_str(),qcsryfptr);
                fflush(qcsryfptr);
            }
            // end of QC
            
            snpinfolst* sortptr=&snpinfo[0];
            qsort(sortptr,snpinfo.size(),sizeof(snpinfolst),comp_esi);
            
           
            
            vector<int> slct_idx;
            slct_sparse_per_prb(slct_idx, &prbifo, snpinfo,  cis_itvl,  trans_itvl, transThres, restThres,logfile); //slct_idx with no order if there are trans-rgeions
            stable_sort(slct_idx.begin(),slct_idx.end());
            vector<string> _rs(slct_idx.size()), _a1(slct_idx.size()),_a2(slct_idx.size());
            vector<float> _beta(slct_idx.size()), _se(slct_idx.size());
            
            for(int l=0;l<slct_idx.size();l++) {
                _rs[l]=snpinfo[slct_idx[l]].snprs;
                _a1[l]=snpinfo[slct_idx[l]].a1;
                _a2[l]=snpinfo[slct_idx[l]].a2;
                _beta[l]=snpinfo[slct_idx[l]].beta;
                _se[l]=snpinfo[slct_idx[l]].se;
            }
            vector<int> rsid(_rs.size());
            #pragma omp parallel for private(iter)
            for (int l = 0; l<_rs.size(); l++){
                iter = esi_map.find(_rs[l]);
                if (iter != esi_map.end()) rsid[l]=iter->second;
                else {
                    printf("ERROR: SNP is not in SNP map. Please report this bug.");
                    exit(EXIT_FAILURE);
                }
            }
            map<string, int> rsa_map;
            long rsNum=0;
            for(int l=0;l<rsid.size();l++)
            {
                if(abs(_se[l]+9)>1e-6) // can move this. the NA is controled in slct_sparse_per_prb
                {
                   // string chckstr=_rs[l]+":"+_a1[l]+":"+_a2[l];
                  //  rsa_map.insert(pair<string,int>(chckstr,l)); // in slct_sparse_per_prb, ras_map can privent selecting duplicate SNPs and double-slelecting SNPs. so we can move rsa_map here.
                  //  if(rsNum<rsa_map.size())
                  //  {
                        
                        if(esi_a1[rsid[l]]==_a1[l] && esi_a2[rsid[l]]==_a2[l] ){
                            val.push_back(_beta[l]);
                            rowids.push_back(rsid[l]);
                            tmpse.push_back(_se[l]);
                            tmprid.push_back(rsid[l]);
                        } else if(esi_a1[rsid[l]]==_a2[l] && esi_a2[rsid[l]]==_a1[l] ){
                            val.push_back(-1.0*_beta[l]);
                            rowids.push_back(rsid[l]);
                            tmpse.push_back(_se[l]);
                            tmprid.push_back(rsid[l]);
                        } else {
                            int did=-9;
                            float sig=1.0;
                            for(int m=0;m<esi_rs.size();m++)
                            {
                                if(esi_rs[m]==_rs[l])
                                {
                                    if(esi_a1[m]==_a1[l] && esi_a2[m]==_a2[l])
                                    {
                                        did=m;
                                        break;
                                    }
                                    if(esi_a1[m]==_a2[l] && esi_a2[m]==_a1[l])
                                    {
                                        did=m;
                                        sig=-1.0;
                                        break;
                                    }
                                }
                            }
                            if(did==-9)
                            {
                                printf("ERROR: This would not go to happen. Please report this bug.");
                                exit(EXIT_FAILURE);
                            }
                            val.push_back(sig*_beta[l]);
                            rowids.push_back(did);
                            tmpse.push_back(_se[l]);
                            tmprid.push_back(did);
                        }
                        
                        rsNum=rsa_map.size();
                 //   } else {
                  //      printf("WARNING: duplicate SNP %s with the same alleles belonging to the same probe %s. \n",_rs[l].c_str(),prbname.c_str());
                 //   }
                } else {
                    printf("WARNING: SNP %s  with \"NA\" value found and skipped.\n",_rs[l].c_str());
                }
            }
            for(int k=0;k<tmpse.size();k++)
            {
                val.push_back(tmpse[k]);
                rowids.push_back(tmprid[k]);
            }
            uint64_t real_num=tmpse.size();
            cols[(i<<1)+1]=real_num+cols[i<<1];
            cols[i+1<<1]=(real_num<<1)+cols[i<<1];
        }
        
        uint64_t valNum2write=val.size();
        fwrite (&valNum2write,sizeof(uint64_t), 1, smr1);
        fwrite (&cols[0],sizeof(uint64_t), cols.size(), smr1);
        fwrite (&rowids[0],sizeof(uint32_t), rowids.size(), smr1);
        fwrite (&val[0],sizeof(float), val.size(), smr1);
        fclose (smr1);
        
        printf("The summary information of selection has been saved in %s.\n", logfname.c_str());
        cout<<"\nBeta values and SE values for "<<etmp._include.size()<<" Probes have been saved in the sparse binary file [" + esdfile + "]." <<endl;
       

        if(buffer) free(buffer);
        fclose(fptr);
        fclose(logfile);
        if(qcflag)
        {
            fclose(outlierfptr);
            fclose(qcsryfptr);
        }
        
        
    }
    
    void diff(char* eqtlFileName1,char* eqtlFileName2)
    {
        if(eqtlFileName1==NULL) throw("Error: please input eQTL summary data for diff by the flag --eqtl-summary.");
        if(eqtlFileName2==NULL) throw("Error: please input eQTL summary data for diff by the flag --eqtl-summary.");
        bool diff=false;
        eqtlInfo edata1;
        eqtlInfo edata2;
        read_esifile(&edata1, string(eqtlFileName1)+".esi");
        read_esifile(&edata2, string(eqtlFileName2)+".esi");
        read_epifile(&edata1, string(eqtlFileName1)+".epi");
        read_epifile(&edata2, string(eqtlFileName2)+".epi");
        // log file
        FILE* logfile=NULL;
        string logfname = "diff.out";
        logfile=fopen(logfname.c_str(), "w");
        if (!(logfile)) {
            printf("Error: Failed to open log file.\n");
            exit(1);
        }
        string logstr="";
        if(edata1._snpNum!=edata2._snpNum  || edata1._probNum!=edata2._probNum )
        {
            diff=true;
            logstr="Failed: no consistence in SNP or probe information.\n";
            fputs(logstr.c_str(),logfile);
            fflush(logfile);
            printf("%s",logstr.c_str());
            return;
        }
        for(int i=0;i<edata1._snpNum;i++)
        {
            if(edata1._esi_rs[i]!=edata2._esi_rs[i]  || edata1._esi_bp[i]!=edata2._esi_bp[i] || edata1._esi_chr[i]!=edata2._esi_chr[i])
            {
                diff=true;
                logstr="Failed: no consistence in SNP information.\n";
                fputs(logstr.c_str(),logfile);
                fflush(logfile);
                printf("%s",logstr.c_str());
                return;
            }
        }
        for(int i=0;i<edata1._probNum;i++)
        {
            if(edata1._epi_chr[i]!=edata2._epi_chr[i]  || edata1._epi_prbID[i]!=edata2._epi_prbID[i] || edata1._epi_bp[i]!=edata2._epi_bp[i])
            {
                diff=true;
                logstr="Failed: no consistence in probe information.\n";
                fputs(logstr.c_str(),logfile);
                fflush(logfile);
                printf("%s",logstr.c_str());
                return;
            }
        }
     
        
        string besdfile1 = string(eqtlFileName1)+".besd";
        FILE *fptr1=fopen(besdfile1.c_str(), "rb");
        if(!fptr1)
        {
            printf ( "ERROR: Couldn't open file %s\n", besdfile1.c_str());
            exit (EXIT_FAILURE);
        }
        float filetype1=readfloat(fptr1);
        
        string besdfile2 = string(eqtlFileName1)+".besd";
        FILE *fptr2=fopen(besdfile1.c_str(), "rb");
        if(!fptr2)
        {
            printf ( "ERROR: Couldn't open file %s\n", besdfile2.c_str());
            exit (EXIT_FAILURE);
        }
        float filetype2=readfloat(fptr2);
        if(filetype1!=filetype2) {
            printf("different besd file format.\n");
            return;
        }
        if((int)filetype1==SPARSE_FILE_TYPE_3 ){
            
            uint64_t colNum=(edata1._probNum<<1)+1;
            fseek(fptr1, 0L, SEEK_END);
            uint64_t lSize = ftell(fptr1);
            fseek(fptr1, 0L, SEEK_SET);
            readfloat(fptr1);
            uint64_t valNum=readuint64(fptr1);
            if( lSize - (sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float)) != 0) {fputs ("wrong element number",stderr); exit (3);}
            fseek(fptr2, 0L, SEEK_END);
            uint64_t lSize2 = ftell(fptr2);
            fseek(fptr2, 0L, SEEK_SET);
            readfloat(fptr2);
            uint64_t valNum2=readuint64(fptr2);
            if( lSize2 - (sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum2*sizeof(uint32_t) + valNum2*sizeof(float)) != 0) {fputs ("wrong element number",stderr); exit (3);}
            

            
            uint64_t colsize=colNum*sizeof(uint64_t);
            uint64_t* colbuf1=(uint64_t*)malloc(colsize);
            if(NULL == colbuf1)
            {
                printf("ERROR: Can't malloc the reading cols buffer for %llu MB.\n",(colsize>>20));
                exit(EXIT_FAILURE);
            }
            fread(colbuf1,colNum,sizeof(uint64_t),fptr1);
            uint64_t* colbuf2=(uint64_t*)malloc(colsize);
            if(NULL == colbuf2)
            {
                printf("ERROR: Can't malloc the reading cols buffer for %llu MB.\n",(colsize>>20));
                exit(EXIT_FAILURE);
            }
            fread(colbuf2,colNum,sizeof(uint64_t),fptr2);
            for(int i=0;i<colNum;i++)
            {
                if(colbuf1[i]!=colbuf2[i]) {
                    diff=true;
                    logstr="Failed: no consistence in eQTL number in probe "+ edata1._epi_prbID[i>>1]+".\n";
                    fputs(logstr.c_str(),logfile);
                    fflush(logfile);
                    printf("%s",logstr.c_str());
                    return;
                }
            }
            free(colbuf2);
            uint64_t* ptr=NULL;
            uint64_t rowSTART=0;
            uint64_t valSTART=0;
            ptr=colbuf1;
            rowSTART=sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t);
            valSTART=sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t)+valNum*sizeof(uint32_t);
            
            for( int j=0;j<edata1._probNum;j++){
                uint64_t pos=*(ptr+(j<<1)); //BETA START
                uint64_t pos1=*(ptr+(j<<1)+1); //SE START
                uint64_t num=pos1-pos;
                if(num>0)
                {
                    uint32_t* ridbuff1=(uint32_t*)malloc(num*2*sizeof(uint32_t));
                    if(NULL == ridbuff1)
                    {
                        printf("ERROR: Memory allocation error .\n");
                        fclose(fptr1);
                        fclose(fptr2);
                        exit(EXIT_FAILURE);
                    }
                    memset(ridbuff1,0,num*2*sizeof(uint32_t));
                    float* betasebuff1=(float*)malloc(num*2*sizeof(float));
                    if(NULL == betasebuff1)
                    {
                        printf("ERROR: Memory allocation error.\n");
                        fclose(fptr1);
                        fclose(fptr2);
                        exit(EXIT_FAILURE);
                    }
                    memset(betasebuff1,0,num*2*sizeof(float));
                    
                    fseek(fptr1, rowSTART+pos*sizeof(uint32_t), SEEK_SET);
                    fread(ridbuff1, sizeof(uint32_t),2*num,fptr1);
                    fseek(fptr1,valSTART+pos*sizeof(float),SEEK_SET);
                    fread(betasebuff1,sizeof(float), 2*num,fptr1);
                    
                    uint32_t* ridbuff2=(uint32_t*)malloc(num*2*sizeof(uint32_t));
                    if(NULL == ridbuff2)
                    {
                        printf("ERROR: Memory allocation error .\n");
                        fclose(fptr1);
                        fclose(fptr2);
                        exit(EXIT_FAILURE);
                    }
                    memset(ridbuff2,0,num*2*sizeof(uint32_t));
                    float* betasebuff2=(float*)malloc(num*2*sizeof(float));
                    if(NULL == betasebuff2)
                    {
                        printf("ERROR: Memory allocation error.\n");
                        fclose(fptr1);
                        fclose(fptr2);
                        exit(EXIT_FAILURE);
                    }
                    memset(betasebuff2,0,num*2*sizeof(float));
                    
                    fseek(fptr2, rowSTART+pos*sizeof(uint32_t), SEEK_SET);
                    fread(ridbuff2, sizeof(uint32_t),2*num,fptr2);
                    fseek(fptr2,valSTART+pos*sizeof(float),SEEK_SET);
                    fread(betasebuff2,sizeof(float), 2*num,fptr2);
                    vector<uint32_t> rowid1;
                    vector<float> beta1;
                    vector<float> se1;
                    for(int jj=0;jj<num;jj++)
                    {
                        double beta=betasebuff1[jj];
                        double se=betasebuff1[jj+num];
                        if(ABS(se+9)<1e-6) continue;
                        rowid1.push_back(ridbuff1[jj]);
                        beta1.push_back(beta);
                        se1.push_back(se);
                        
                    }
                    vector<uint32_t> rowid2;
                    vector<float> beta2;
                    vector<float> se2;
                    for(int jj=0;jj<num;jj++)
                    {
                        double beta=betasebuff2[jj];
                        double se=betasebuff2[jj+num];
                        if(ABS(se+9)<1e-6) continue;
                        rowid2.push_back(ridbuff2[jj]);
                        beta2.push_back(beta);
                        se2.push_back(se);
                    }
                    vector<int> indx;
                    match(rowid1, rowid2, indx);
                    for(int jj=0;jj<indx.size();jj++)
                    {
                        if(indx[jj]==-9){
                            diff=true;
                            logstr="Failed: no consistence in probe "+ edata1._epi_prbID[j]+".\n";
                            fputs(logstr.c_str(),logfile);
                            fflush(logfile);
                            printf("%s",logstr.c_str());
                            fclose(fptr1);
                            fclose(fptr2);
                            exit(EXIT_FAILURE);
                        }
                        if(abs(beta1[jj]-beta2[indx[jj]])>1e-6 || abs(se1[jj]-se2[indx[jj]])>1e-6)
                        {
                            diff=true;
                            logstr=edata1._epi_prbID[j]+'\t'+edata1._esi_rs[rowid1[jj]]+'\t'+atos(beta1[jj])+'\t'+atos(beta2[indx[jj]])+'\t'+atos(se1[jj])+'\t'+atos(se2[indx[jj]])+'\n';
                            fputs(logstr.c_str(),logfile);
                            fflush(logfile);
                            printf("%s",logstr.c_str());
                        }
                    }
                   

                    free(ridbuff1);
                    free(betasebuff1);
                    free(ridbuff2);
                    free(betasebuff2);
                    
                }
            }
            if(!diff){
                logstr="PASSED: the files identify with each other.\n";
                fputs(logstr.c_str(),logfile);
                fflush(logfile);
                printf("%s",logstr.c_str());
                
            }
        }
        else if((int)filetype1==DENSE_FILE_TYPE_1)
        {
            fseek(fptr1, 0L, SEEK_END);
            uint64_t lSize = ftell(fptr1);
            fseek(fptr1, 0L, SEEK_SET);
            fseek(fptr2, 0L, SEEK_END);
            uint64_t lSize2 = ftell(fptr1);
            fseek(fptr2, 0L, SEEK_SET);
            if(lSize!=lSize2)
            {
                diff=true;
                logstr="Failed: no consistence in file size.\n";
                fputs(logstr.c_str(),logfile);
                fflush(logfile);
                printf("%s",logstr.c_str());
                fclose(fptr1);
                fclose(fptr2);
                exit(EXIT_FAILURE);

            }
            while(!feof(fptr1) && !feof(fptr2))
            {
                if(abs(readfloat(fptr1)-readfloat(fptr2))>1e-6)
                {
                    diff=true;
                    logstr="Failed: no consistence in file content.\n";
                    fputs(logstr.c_str(),logfile);
                    fflush(logfile);
                    printf("%s",logstr.c_str());
                    fclose(fptr1);
                    fclose(fptr2);
                    exit(EXIT_FAILURE);
                }
            }
            if(!diff){
                logstr="PASSED: the files identify with each other.\n";
                fputs(logstr.c_str(),logfile);
                fflush(logfile);
                printf("%s",logstr.c_str());
                
            }
        }
        
        fclose(fptr1);
        fclose(fptr2);
        fclose(logfile);
        
        
    }
    
    void beqtl_qc_se(eqtlInfo* eqtlinfo, int qc_mtd, int z_thresh, char* outFileName)
    {
        FILE* outlierfptr=NULL;
        FILE* qcsryfptr=NULL;
        string logstr="";
        // outlier file
        string outlierfname = string(outFileName)+".outlier";
        outlierfptr=fopen(outlierfname.c_str(), "w");
        if (!(outlierfptr)) {
            printf("Error: Failed to open log file.\n");
            exit(1);
        }
        logstr="SNP\tChr\tBP\tA1\tA2\tFreq\tProbe\tProbe_Chr\tProbe_bp\tGene\tOrientation\tb\tse\test_n\n";
        fputs(logstr.c_str(),outlierfptr);
        fflush(outlierfptr);
        // qc summary file
        string qcsryfname = string(outFileName)+".qc.summary";
        qcsryfptr=fopen(qcsryfname.c_str(), "w");
        if (!(qcsryfptr)) {
            printf("Error: Failed to open log file.\n");
            exit(1);
        }
        logstr="ProbeID\tmaxn\tminn\tavgn\tmediann\tsdn\tttlnum\trmnum\n";
        fputs(logstr.c_str(),qcsryfptr);
        fflush(qcsryfptr);
        
        
        bool warnnullfrqflag=false;
        vector<snpinfolst> snpinfo;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            printf("Saving... %3.0f%%\r", 100.0*i/eqtlinfo->_include.size());
            fflush(stdout);
            int prbidx=eqtlinfo->_include[i];
            string prbid=eqtlinfo->_epi_prbID[prbidx];
            
            if(eqtlinfo->_rowid.empty())
            {
                for (int j = 0; j<eqtlinfo->_esi_include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
                {
                    int snpidx=eqtlinfo->_esi_include[j];
                    if (abs(eqtlinfo->_bxz[prbidx][snpidx] + 9) > 1e-6)
                    {
                        snpinfolst snptmp;
                        float se=eqtlinfo->_sexz[prbidx][snpidx];
                        if(se<=0)
                        {
                            printf("ERROR: se >  0 . \n");
                            exit(EXIT_FAILURE);
                        }
                        float freq=eqtlinfo->_esi_freq[snpidx];
                        if(freq>=1 || freq<=0)
                        {
                            printf("ERROR: freqency should be between 0 and 1. \n");
                            exit(EXIT_FAILURE);
                        }
                        if(abs(freq+9)<1e-6 && !warnnullfrqflag)
                        {
                            printf("WANRNING: one or more freqencies are NA, these SNPs would be removed. \n");
                            warnnullfrqflag=true;
                            continue;
                        }
                        snptmp.snprs=eqtlinfo->_esi_rs[snpidx];
                        snptmp.snpchr=eqtlinfo->_esi_chr[snpidx];
                        snptmp.bp=eqtlinfo->_esi_bp[snpidx];
                        snptmp.a1=eqtlinfo->_esi_bp[snpidx];
                        snptmp.a2=eqtlinfo->_esi_bp[snpidx];
                        snptmp.freq=freq;
                        snptmp.beta=eqtlinfo->_bxz[prbidx][snpidx];
                        snptmp.se=se;
                        snptmp.estn=est_sample_size(freq, snptmp.beta, se);
                        snptmp.gd=snpidx;
                        snpinfo.push_back(snptmp);
                    }
                }
            }
            else{
                uint64_t beta_start=eqtlinfo->_cols[prbidx<<1];
                uint64_t se_start=eqtlinfo->_cols[1+(prbidx<<1)];
                uint64_t numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    int snpidx=eqtlinfo->_rowid[beta_start+j];
                    snpinfolst snptmp;
                    float se=eqtlinfo->_val[se_start+j];
                    float freq=eqtlinfo->_esi_freq[snpidx];
                    if(se<=0)
                    {
                        printf("ERROR: se >0. \n");
                        exit(EXIT_FAILURE);
                    }
                    if(freq>=1 || freq<=0)
                    {
                        printf("ERROR: freqency should be between 0 and 1. \n");
                        exit(EXIT_FAILURE);
                    }
                    if(abs(freq+9)<1e-6 && !warnnullfrqflag)
                    {
                        printf("WANRNING: one or more freqencies are NA, these SNPs would be removed. \n");
                        warnnullfrqflag=true;
                        continue;
                    }
                    snptmp.beta=eqtlinfo->_val[beta_start+j];
                    snptmp.se=se;
                    snptmp.estn= est_sample_size(freq, snptmp.beta, se);
                    snptmp.gd=snpidx;
                    snptmp.snprs=eqtlinfo->_esi_rs[snpidx];
                    snptmp.snpchr=eqtlinfo->_esi_chr[snpidx];
                    snptmp.a1=eqtlinfo->_esi_allele1[snpidx];
                    snptmp.a2=eqtlinfo->_esi_allele2[snpidx];
                    snptmp.bp=eqtlinfo->_esi_bp[snpidx];
                    snptmp.freq=eqtlinfo->_esi_freq[snpidx];
                    snpinfo.push_back(snptmp);
                }
            }
            
            probeinfolst prbifo;
            prbifo.bp=eqtlinfo->_epi_bp[prbidx];
            prbifo.probechr=eqtlinfo->_epi_chr[prbidx];
            prbifo.probeId=eqtlinfo->_epi_prbID[prbidx];
            prbifo.genename=eqtlinfo->_epi_gene[prbidx];
            prbifo.orien=eqtlinfo->_epi_orien[prbidx];
            
            //  QC
            if( snpinfo.size()>0) {
                vector<float> suminfo; // maxn,minn,avgn,mediann,sdn,ttnum,rmnum
                qc(snpinfo,&prbifo,qc_mtd,z_thresh, suminfo,outlierfptr);
                logstr=prbifo.probeId;
                for(int j=0;j<suminfo.size();j++)
                {
                    logstr+='\t'+atos(suminfo[j]);
                }
                logstr+='\n';
                fputs(logstr.c_str(),qcsryfptr);
                fflush(qcsryfptr);
            }
            // end of QC
            
            snpinfolst* sortptr=&snpinfo[0];
            qsort(sortptr,snpinfo.size(),sizeof(snpinfolst),comp_esi);

            if(eqtlinfo->_rowid.empty())
            {
                
            }
            else {
                
            }

        }
            
        
    }
}