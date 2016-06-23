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
    void combine_esi(vector<snpinfolst> &snpinfo, vector<string> &smasNames)
    {
        long counter = 0;
        map<string, int> rs_map;
        map<string, int> rsbp_map;
     
        for (int i = 0; i < smasNames.size(); i++)
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
    void save_besds_dbesd(char* outFileName, vector<snpinfolst> &snpinfo, vector<probeinfolst2> &probeinfo, vector<string> &esi_rs,vector<string> &esi_a1,vector<string> &esi_a2)
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
                eqtlInfo eqtlinfo;
                read_epifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".epi", prtscr);
                extract_eqtl_single_probe(&eqtlinfo, prbname, prtscr);
                read_esifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".esi", prtscr);
                read_besdfile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".besd", prtscr);
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
                        printf("ERROR: SNP is not in SNP map. Please report this bug.");
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
    
    void save_besds_sbesd(char* outFileName, vector<snpinfolst> &snpinfo, vector<probeinfolst2> &probeinfo, vector<string> &esi_rs,vector<string> &esi_a1,vector<string> &esi_a2)
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
                eqtlInfo eqtlinfo;
                read_epifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".epi", prtscr);
                extract_eqtl_single_probe(&eqtlinfo, prbname, prtscr);
                read_esifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".esi", prtscr);
                read_besdfile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".besd", prtscr);
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
                        printf("ERROR: SNP is not in SNP map. Please report this bug.");
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

    void save_slct_besds_sbesd(char* outFileName,  long esiNum, vector<probeinfolst2> &probeinfo,vector<string> &esi_rs,vector<string> &esi_a1,vector<string> &esi_a2,int cis_itvl,int trans_itvl,float transThres,float restThres)
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
        logstr+="\ncis region is represent as [Chr, Start bp, End bp, nsnp];\ntrans region is represent as <Chr, Start bp, End bp, nsnp>;\nthe number of other SNPs selected is represent as (NumSNPs beyond cis and trans).\n";
        
        logstr+="\n{ProbeID, ProbeChr, ProbeBP}\t[Chr,cis_startBP,cis_endBP,NumSNPs]\t<Chr,trans_startBP,trans_endBP,NumSNPs>\t(NumSNPs)\n";
        fputs(logstr.c_str(),logfile);
        fflush(logfile);
        cis_itvl=cis_itvl*1000;
        trans_itvl=trans_itvl*1000;
        vector<snpinfolst> snpinfo;
        for(int j=0;j<epiNum;j++)
        {
            printf("Saving... %3.0f%%\r", 100.0*j/epiNum);
            fflush(stdout);
            string prbname=probeinfo[j].probeId;
            vector<uint32_t> tmprid;
            vector<float> tmpse;
            
            snpinfo.clear();
            map<string,int> snp_map;
            int snpmapsize=0;
            for(int k=0;k<probeinfo[j].besdpath.size();k++)
            {
                eqtlInfo eqtlinfo;
                read_epifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".epi",prtscr);
                extract_eqtl_single_probe(&eqtlinfo, prbname,prtscr);
                read_esifile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".esi",prtscr);
                read_besdfile(&eqtlinfo, string(probeinfo[j].besdpath[k])+".besd",prtscr);
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
                qsort(sortptr,snpinfo.size(),sizeof(snpinfolst),comp_esi);
                
                probeinfolst prbifo;
                prbifo.bp=probeinfo[j].bp;
                prbifo.probechr=probeinfo[j].probechr;
                prbifo.probeId=probeinfo[j].probeId;
                
                vector<int> slct_idx;
                slct_sparse_per_prb(slct_idx, &prbifo, snpinfo,  cis_itvl,  trans_itvl, transThres, restThres,logfile); //slct_idx with no order if there are trans-rgeions
                
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
                        string chckstr=_rs[l]+":"+_a1[l]+":"+_a2[l];
                        rsa_map.insert(pair<string,int>(chckstr,l)); // in slct_sparse_per_prb, ras_map can privent selecting duplicate SNPs and double-slelecting SNPs. so we can move rsa_map here.
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
                            printf("WARNING: duplicate SNP %s with the same alleles belonging to the same probe %s. \n",_rs[l].c_str(),prbname.c_str());
                        }
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
    void combineBesd(char* eqtlsmaslstName, char* outFileName,bool makecisflag, int cis_itvl, int trans_itvl, float transThres, float restThres)
    {
        vector<string> smasNames;
        vector<snpinfolst> snpinfo;
        vector<probeinfolst2> probeinfo;
        
        read_smaslist(smasNames, string(eqtlsmaslstName));
        if(smasNames.size()==0) throw("No eqtl summary file list in [ "+ string(eqtlsmaslstName)  +" ]");
        combine_esi(snpinfo, smasNames);
        snpinfolst* esiptr=&snpinfo[0];
        qsort(esiptr,snpinfo.size(),sizeof(snpinfolst),comp_esi);
        combine_epi(probeinfo, smasNames);
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
        vector<string> esi_rs(snpinfo.size());
        vector<string> esi_a1(snpinfo.size());
        vector<string> esi_a2(snpinfo.size());
        string esifile =  string(outFileName)+string(".esi");
        ofstream esi(esifile.c_str());
        if (!esi) throw ("Error: can not open the esi file to save!");
        for (long j = 0;j <snpinfo.size(); j++) {
            esi<<snpinfo[j].snpchr<<'\t'<<snpinfo[j].snprs<<'\t'<<snpinfo[j].gd<<'\t'<<snpinfo[j].bp<<'\t'<<snpinfo[j].a1<<'\t'<<snpinfo[j].a2<<'\n';
            esi_rs[j]=snpinfo[j].snprs;
            esi_a1[j]=snpinfo[j].a1;
            esi_a2[j]=snpinfo[j].a2;
        }
        esi.close();
        printf("%ld SNPs have been saved in the file %s.\n",snpinfo.size(),esifile.c_str());
        
        uint64_t valnum=countNotNullNum(smasNames);
        double sparsity=1.0*valnum/(probeinfo.size()*snpinfo.size());
        if(!makecisflag)
        {
            if(sparsity>=0.4)
            {
                printf("The sparsity of your data matrix is %f. We are going to save it in dense format!\n", sparsity);
                save_besds_dbesd(outFileName, snpinfo, probeinfo,esi_rs,esi_a1,esi_a2);
            } else {
                printf("The sparsity of your data matrix is %f. We are going to save it in sparse format!\n", sparsity);
                save_besds_sbesd( outFileName, snpinfo, probeinfo,esi_rs,esi_a1,esi_a2);
            }
        } else
        {
             long esiNum=snpinfo.size();
            save_slct_besds_sbesd(outFileName, esiNum, probeinfo,esi_rs,esi_a1,esi_a2, cis_itvl,  trans_itvl,  transThres,  restThres);
        }
    }
    void get_snpinfo_cur_prb_sparse(vector<snpinfolst> &snpinfo,FILE* fptr,  uint32_t pid, uint64_t* ptr, uint64_t rowSTART,uint64_t valSTART,eqtlInfo* etmp,map<int, int> &_incld_id_map){
        
        uint64_t pos=*(ptr+(pid<<1)); //BETA START
        uint64_t pos1=*(ptr+(pid<<1)+1); //SE START
        uint64_t num=pos1-pos;
        uint64_t real_num=0;
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
                snpinfotmp.beta=*(val_ptr+j);
                snpinfotmp.se=*(val_ptr+j+num);
                snpinfo.push_back(snpinfotmp);
                //int sid=iter->second;
                //cout<<rid<<":"<<etmp._esi_include[sid]<<endl; // test passed
                real_num++;
            }
            
        }
        free(row_char_ptr);
        free(val_char_ptr);

    }
    void get_snpinfo_cur_prb_dense(vector<snpinfolst> &snpinfo,FILE* fptr,  uint32_t pid, char* buffer, eqtlInfo* etmp)
    {
        fseek(fptr,((pid<<1)*etmp->_snpNum+1)<<2, SEEK_SET);
        
        memset(buffer,0,sizeof(char)*etmp->_snpNum<<3);
        fread(buffer, sizeof(char),etmp->_snpNum<<3,fptr);
        float* ft=(float *)buffer;
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
                snpinfotmp.beta=*(ft + etmp->_esi_include[j]);
                snpinfotmp.se=se;
                snpinfo.push_back(snpinfotmp);
            }
        }
        

    }
    void make_sparse_besd(char* eqtlFileName, char* outFileName, int cis_itvl, int trans_itvl, float transThres, float restThres)
    {
        
        eqtlInfo etmp;
        read_esifile(&etmp, string(eqtlFileName)+".esi");
        // epi_man(&eqtlinfo, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
       // if(problst2exclde != NULL) exclude_prob(&eqtlinfo, problst2exclde);
        read_epifile(&etmp, string(eqtlFileName)+".epi");
      //   esi_man(&eqtlinfo, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,  cis_itvl, prbname);
      //  if(snplst2exclde != NULL) exclude_eqtl_snp(&eqtlinfo, snplst2exclde);

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
            esi<< etmp._esi_chr[etmp._esi_include[j]]<<'\t'<<etmp._esi_rs[etmp._esi_include[j]]<<'\t'<<etmp._esi_gd[etmp._esi_include[j]]<<'\t'<<etmp._esi_bp[etmp._esi_include[j]]<<'\t'<<etmp._esi_allele1[etmp._esi_include[j]]<<'\t'<<etmp._esi_allele2[etmp._esi_include[j]]<<'\n';
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
        logstr+="\ncis region is represent as [Chr, Start bp, End bp, nsnp];\ntrans region is represent as <Chr, Start bp, End bp, nsnp>;\nthe number of other SNPs selected is represent as (NumSNPs beyond cis and trans).\n";
        
        logstr+="\n{ProbeID, ProbeChr, ProbeBP}\t[Chr,cis_startBP,cis_endBP,NumSNPs]\t<Chr,trans_startBP,trans_endBP,NumSNPs>\t(NumSNPs)\n";
        fputs(logstr.c_str(),logfile);
        fflush(logfile);
        cis_itvl=cis_itvl*1000;
        trans_itvl=trans_itvl*1000;
        vector<snpinfolst> snpinfo;
        map<string, int>::iterator iter;
        for(int i=0;i<etmp._include.size();i++)
        {
            printf("Saving... %3.0f%%\r", 100.0*i/etmp._include.size());
            fflush(stdout);
            string prbname=etmp._epi_prbID[etmp._include[i]];
            vector<uint32_t> tmprid;
            vector<float> tmpse;
            snpinfo.clear();
            uint32_t pid=etmp._include[i];
            if((int)filetype==SPARSE_FILE_TYPE_3) get_snpinfo_cur_prb_sparse(snpinfo,fptr, pid, ptr,  rowSTART, valSTART, &etmp,_incld_id_map);
            else get_snpinfo_cur_prb_dense(snpinfo,fptr, pid, buffer ,&etmp);
            
            snpinfolst* sortptr=&snpinfo[0];
            qsort(sortptr,snpinfo.size(),sizeof(snpinfolst),comp_esi);
            
            probeinfolst prbifo;
            prbifo.bp=etmp._epi_bp[etmp._include[i]];
            prbifo.probechr=etmp._epi_chr[etmp._include[i]];
            prbifo.probeId=etmp._epi_prbID[etmp._include[i]];
            
            vector<int> slct_idx;
            slct_sparse_per_prb(slct_idx, &prbifo, snpinfo,  cis_itvl,  trans_itvl, transThres, restThres,logfile); //slct_idx with no order if there are trans-rgeions
            
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
        
        
    }
}