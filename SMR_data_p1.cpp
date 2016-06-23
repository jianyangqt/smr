//
//  SMR_data_p1.cpp
//  SMR_CPP
//
//  Created by Futao Zhang on 10/03/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "SMR_data_p1.h"
namespace SMRDATA
{
    void get_top_sets(eqtlInfo* eqtlinfo, vector<string> &prbIds, vector<float> &beta, vector<float> &se, vector<string> &rs, float thres)
    {
        if(eqtlinfo->_rowid.empty())
        {
            for(int i=0;i<eqtlinfo->_include.size();i++)
            {
                string probeId=eqtlinfo->_epi_prbID[eqtlinfo->_include[i]];
                int probechr=eqtlinfo->_epi_chr[eqtlinfo->_include[i]];
                double top_zsqr=0, top_beta=0, top_se=0;
                string snprs="";
                for (int j = 0; j<eqtlinfo->_esi_include.size(); j++)
                {
                    if (abs(eqtlinfo->_bxz[eqtlinfo->_include[i]][eqtlinfo->_esi_include[j]] + 9) > 1e-6)
                    {
                        int snpchr=eqtlinfo->_esi_chr[eqtlinfo->_esi_include[j]];
                        string rstmp=eqtlinfo->_esi_rs[eqtlinfo->_esi_include[j]];
                        if(snpchr==probechr)
                        {
                            float bxz=eqtlinfo->_bxz[eqtlinfo->_include[i]][eqtlinfo->_esi_include[j]];
                            float sexz=eqtlinfo->_sexz[eqtlinfo->_include[i]][eqtlinfo->_esi_include[j]];
                            float zxz=bxz/sexz;
                            zxz*=zxz;
                            if(zxz-top_zsqr>1e-8) {
                                top_zsqr=zxz;
                                top_beta=bxz;
                                top_se=sexz;
                                snprs=rstmp;
                            }
                        }
                    }
                }
                double pval=pchisq(top_zsqr, 1);
                if(pval<thres)
                {
                    prbIds.push_back(probeId);
                    beta.push_back(top_beta);
                    se.push_back(top_se);
                    rs.push_back(snprs);
                }
            }
            
        }
        else
        {
            for(int i=0;i<eqtlinfo->_include.size();i++)
            {
                string probeId=eqtlinfo->_epi_prbID[eqtlinfo->_include[i]];
                int probechr=eqtlinfo->_epi_chr[eqtlinfo->_include[i]];
                double top_zsqr=0, top_beta=0, top_se=0;
                string snprs="";
                
                uint64_t beta_start=eqtlinfo->_cols[eqtlinfo->_include[i]<<1];
                uint64_t se_start=eqtlinfo->_cols[1+(eqtlinfo->_include[i]<<1)];
                uint64_t numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=eqtlinfo->_rowid[beta_start+j];
                    int snpchr=eqtlinfo->_esi_chr[ge_rowid];
                    string snptmp=eqtlinfo->_esi_rs[ge_rowid];
                    if(snpchr==probechr )
                    {
                        float bxz=eqtlinfo->_val[beta_start+j];
                        float sexz=eqtlinfo->_val[se_start+j];
                        float zxz=bxz/sexz;
                        zxz*=zxz;
                        if(zxz-top_zsqr>1e-8) {
                            top_zsqr=zxz;
                            top_beta=bxz;
                            top_se=sexz;
                            snprs=snptmp;
                        }
                    }
                }
                
                double pval=pchisq(top_zsqr, 1);
                if(pval<thres)
                {
                    prbIds.push_back(probeId);
                    beta.push_back(top_beta);
                    se.push_back(top_se);
                    rs.push_back(snprs);
                }
                
            }
        }

    }
    
    void get_thres_sets(eqtlInfo* eqtlinfo, vector<string> &prbIds, vector<float> &beta, vector<float> &se, vector<string> &rs, float thres)
    {
        if(eqtlinfo->_rowid.empty())
        {
            for(int i=0;i<eqtlinfo->_include.size();i++)
            {
                string probeId=eqtlinfo->_epi_prbID[eqtlinfo->_include[i]];
                int probechr=eqtlinfo->_epi_chr[eqtlinfo->_include[i]];
                for (int j = 0; j<eqtlinfo->_esi_include.size(); j++)
                {
                    if (abs(eqtlinfo->_bxz[eqtlinfo->_include[i]][eqtlinfo->_esi_include[j]] + 9) > 1e-6)
                    {
                        int snpchr=eqtlinfo->_esi_chr[eqtlinfo->_esi_include[j]];
                        string rstmp=eqtlinfo->_esi_rs[eqtlinfo->_esi_include[j]];
                        if(snpchr==probechr)
                        {
                            float bxz=eqtlinfo->_bxz[eqtlinfo->_include[i]][eqtlinfo->_esi_include[j]];
                            float sexz=eqtlinfo->_sexz[eqtlinfo->_include[i]][eqtlinfo->_esi_include[j]];
                            float zxz=bxz/sexz;
                            zxz*=zxz;
                            double pval=pchisq(zxz, 1);
                            if(pval<thres) {
                                prbIds.push_back(probeId);
                                beta.push_back(bxz);
                                se.push_back(sexz);
                                rs.push_back(rstmp);
                            }
                        }
                    }
                }
            }
            
        }
        else
        {
            for(int i=0;i<eqtlinfo->_include.size();i++)
            {
                string probeId=eqtlinfo->_epi_prbID[eqtlinfo->_include[i]];
                int probechr=eqtlinfo->_epi_chr[eqtlinfo->_include[i]];
                
                uint64_t beta_start=eqtlinfo->_cols[eqtlinfo->_include[i]<<1];
                uint64_t se_start=eqtlinfo->_cols[1+(eqtlinfo->_include[i]<<1)];
                uint64_t numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=eqtlinfo->_rowid[beta_start+j];
                    int snpchr=eqtlinfo->_esi_chr[ge_rowid];
                    string snptmp=eqtlinfo->_esi_rs[ge_rowid];
                    if(snpchr==probechr )
                    {
                        float bxz=eqtlinfo->_val[beta_start+j];
                        float sexz=eqtlinfo->_val[se_start+j];
                        float zxz=bxz/sexz;
                        zxz*=zxz;
                        double pval=pchisq(zxz, 1);
                        if(pval<thres) {
                            prbIds.push_back(probeId);
                            beta.push_back(bxz);
                            se.push_back(sexz);
                            rs.push_back(snptmp);
                        }
                    }
                }
            }
        }
        
    }

    
    long est_n(vector<float> &beta,vector<float> &se)
    {
        long n;
        VectorXd Y(se.size());
        MatrixXd X(beta.size(),2);
        for(int i=0;i<se.size();i++)
        {
            Y(i)=se[i]*se[i];
            X(i,0)=1.0;
            X(i,1)=beta[i]*beta[i]-Y(i);
        }
        
        MatrixXd XtX_i=(X.transpose()*X).inverse();
        VectorXd w_hat=XtX_i*X.transpose()*Y;
        n=-1/w_hat[1];
        return n;
    }
    void est_effect_splsize(char* eqtlsmaslstName, char* eqtlFileName, char* snplstName,char* problstName,char* snplst2exclde, char* problst2exclde,float thres)
    {
        vector<string> prbIds;
        vector<float> beta;
        vector<float> se;
        vector<string> rs;
        
        vector<string> smasNames;
        if(eqtlsmaslstName!=NULL)
        {
            read_smaslist(smasNames, string(eqtlsmaslstName));
            if(smasNames.size()==0) throw("No eqtl summary file list in [ "+ string(eqtlsmaslstName)  +" ]");
            for(int ii=0;ii<smasNames.size();ii++)
            {
                eqtlInfo eqtlinfo;
                read_esifile(&eqtlinfo, smasNames[ii]+".esi");
                if (snplstName != NULL) extract_eqtl_snp(&eqtlinfo, snplstName);
                if(snplst2exclde != NULL) exclude_eqtl_snp(&eqtlinfo, snplst2exclde);
                read_epifile(&eqtlinfo, smasNames[ii]+".epi");
                if(problstName != NULL) extract_prob(&eqtlinfo, problstName);
                if(problst2exclde != NULL) exclude_prob(&eqtlinfo, problst2exclde);
                read_besdfile(&eqtlinfo, smasNames[ii]+".besd");
                if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
                {
                    printf("No data included from %s under current condition.\n",smasNames[ii].c_str());
                    exit(EXIT_FAILURE);
                }
                //get_top_sets(&eqtlinfo,prbIds,beta,se,rs,thres);
                get_thres_sets(&eqtlinfo,prbIds,beta,se,rs,thres);
            
            }

        } else if(eqtlFileName!=NULL)
        {
            eqtlInfo eqtlinfo;
            read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&eqtlinfo, snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&eqtlinfo, snplst2exclde);
            read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
            if(problstName != NULL) extract_prob(&eqtlinfo, problstName);
            if(problst2exclde != NULL) exclude_prob(&eqtlinfo, problst2exclde);
            read_besdfile(&eqtlinfo, string(eqtlFileName)+".besd");
            if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
            {
                printf("No data included from %s under current condition.\n",eqtlFileName);
                exit(EXIT_FAILURE);
            }
            //get_top_sets(&eqtlinfo,prbIds,beta,se,rs,thres);
            get_thres_sets(&eqtlinfo,prbIds,beta,se,rs,thres);
           
        }else {
            throw("Error: please input eQTL summary data files list or eQTL summary data by the flag --beqtl-summaries or --beqtl-summary.");
        }
        
        cout<<prbIds.size()<<" eQTLs are inluded to estimate the effective population size."<<endl;
        ///////
        
        string tstfile = "tst.top.bse.txt";
        ofstream tst(tstfile.c_str());
        if (!tst) throw ("Error: can not open the fam file " + tstfile + " to save!");
        
        tst << "Probe" <<'\t'<<"SNP"<<'\t'<< "b" <<'\t' << "se"<<'\n';
        
        for (int i = 0;i <beta.size(); i++) {
            tst<<prbIds[i]<<'\t'<<rs[i]<<'\t'<<beta[i]<<'\t'<<se[i]<< '\n';
        }
        
        tst.close();
        
        //////
        
        long n=est_n(beta,se);
        cout<<"The estimated effective population size is: "<<n<<endl;
    }
    
    void read_esmr(eSMRrlt* erlt, string smrFileName)
    {
        cout << "Reading SMR result information from [" + smrFileName + "]." << endl;
       
        ifstream flptr;
        flptr.open(smrFileName.c_str());
        if (!flptr) throw ("Error: can not open the file [" + smrFileName + "] to read.");
        
       
        char tbuf[MAX_LINE_SIZE];
        int lineNum(0);
       
        flptr.getline(tbuf,MAX_LINE_SIZE);// the header
        while(!flptr.eof())
        {
            string tmpStr;
            flptr.getline(tbuf,MAX_LINE_SIZE);
            if(tbuf[0]!='\0'){
                istringstream iss(tbuf);
                iss>>tmpStr;
                erlt->_Expo_id.push_back(tmpStr);
                iss>>tmpStr;
                erlt->_Expo_chr.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_Expo_gene.push_back(tmpStr);
                iss>>tmpStr;
                erlt->_Expo_bp.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_Outco_id.push_back(tmpStr);
                iss>>tmpStr;
                erlt->_Outco_chr.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_Outco_gene.push_back(tmpStr);
                iss>>tmpStr;
                erlt->_Outco_bp.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_snp_rs.push_back(tmpStr);
                iss>>tmpStr;
                erlt->_snp_chr.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_snp_bp.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr;
                to_upper(tmpStr);
                erlt->_snp_a1.push_back(tmpStr.c_str());
                iss>>tmpStr;
                to_upper(tmpStr);
                erlt->_snp_a2.push_back(tmpStr.c_str());
                iss>>tmpStr;
                erlt->_snp_frq.push_back(atof(tmpStr.c_str()));
                iss>>tmpStr; // b_outcome
                iss>>tmpStr; // se_outcome
                iss>>tmpStr; // p_outcome
                iss>>tmpStr; // b_exposure
                iss>>tmpStr; // se_exposure
                iss>>tmpStr; // p_exposure
                iss>>tmpStr;
                erlt->_b.push_back(atof(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_se.push_back(atof(tmpStr.c_str()));
                iss>>tmpStr;
                erlt->_p_smr.push_back(atof(tmpStr.c_str()));
                iss>>tmpStr;
                if(!tmpStr.compare("NA")) erlt->_p_heidi.push_back(-9);
                else  erlt->_p_heidi.push_back(atof(tmpStr.c_str()));
                iss>>tmpStr;
                if(!tmpStr.compare("NA")) erlt->_nsnp.push_back(-9);
                else  erlt->_nsnp.push_back(atoi(tmpStr.c_str()));
                erlt->_include.push_back(lineNum);
                lineNum++;
            }
        }
        erlt->lineNum=lineNum;
        flptr.close();
        cout << lineNum << " SNPs summary info to be included from [" + smrFileName + "]." << endl;
        
    }
    
    
    void read_probevarfile(eqtlInfo* eqtlinfo, string vpFileName)
    {
        cout << "Reading variance information from [" + vpFileName + "]." << endl;
        
        ifstream flptr;
        flptr.open(vpFileName.c_str());
        if (!flptr) throw ("Error: can not open the file [" + vpFileName + "] to read.");
        
        
        char tbuf[MAX_LINE_SIZE];
        int lineNum(0);
        vector<string> tmp_prbid;
        vector<double> tmp_var;
        
        while(!flptr.eof())
        {
            string tmpStr;
            flptr.getline(tbuf,MAX_LINE_SIZE);
            if(tbuf[0]!='\0'){
                istringstream iss(tbuf);
                iss>>tmpStr;
                tmp_prbid.push_back(tmpStr.c_str());
                iss>>tmpStr;
                tmp_var.push_back(atof(tmpStr.c_str()));
                lineNum++;
            }
        }
        flptr.close();
        
        
        eqtlinfo->_epi_var.resize(eqtlinfo->_epi_prbID.size());
        vector<int> idx;
        match_only(eqtlinfo->_epi_prbID, tmp_prbid, idx);
        if(idx.size()!=eqtlinfo->_epi_prbID.size())
        {
            cout<<"Some Probes in summary data are not in variance data!"<<endl;
            exit(1);
        }
        for(int i=0;i<idx.size();i++)
        {
            eqtlinfo->_epi_var[i]=tmp_var[idx[i]];
        }
        
        cout << idx.size() << " probes variance info to be included from [" + vpFileName + "]." << endl;

    }
    

    void iternal_test(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,int cis_itvl,char* smrFileName)
    {
        setNbThreads(thread_num);
        
        eqtlInfo etrait;
        eqtlInfo esdata;
        bInfo bdata;
        double threshold= chi_val(1,p_hetero);
        cis_itvl=cis_itvl*1000;
       
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(problstName != NULL) cout<<"WARNING: --extract-probe here presumes the probe list should contain both probes of exposure dataset and probes of outcome dataset.\n If you want to only extract probes from one dataset please include these probles in the file and all the probes of the other dataset as well.\n"<<endl;
        read_esifile(&etrait, string(eqtlFileName)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&etrait, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait, snplst2exclde);
        read_epifile(&etrait, string(eqtlFileName)+".epi");
        if(problstName != NULL) extract_prob(&etrait, problstName);
        if(problst2exclde != NULL) exclude_prob(&etrait, problst2exclde);
        if(oproblstName != NULL ) extract_prob(&etrait, oproblstName);
        if(oproblst2exclde != NULL) exclude_prob(&etrait, oproblst2exclde);
        
        
        read_besdfile(&etrait, string(eqtlFileName)+".besd");
        if(etrait._rowid.empty() && etrait._bxz.empty())
        {
            printf("No data included from %s under current condition.\n",eqtlFileName);
            exit(EXIT_FAILURE);
        }
        
        read_esifile(&esdata, string(eqtlFileName2)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&esdata, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
        
            read_famfile(&bdata, string(bFileName)+".fam");
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            read_bimfile(&bdata, string(bFileName)+".bim");
            if(snplstName != NULL) extract_snp(&bdata, snplstName);
            if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
            allele_check(&bdata, &etrait, &esdata);
            // if no snp left after check
            read_bedfile(&bdata, string(bFileName)+".bed");
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx(&bdata, &etrait, &esdata);
            }
            

        cout<<"Reading eQTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName2)+".epi");
        if(problstName != NULL) extract_prob(&esdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        if(eproblstName != NULL ) extract_prob(&esdata, eproblstName);
        if(eproblst2exclde != NULL) exclude_prob(&esdata, eproblst2exclde);
        
        read_besdfile(&esdata, string(eqtlFileName2)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included from %s under current condition.\n",eqtlFileName2);
            exit(EXIT_FAILURE);
        }
        
        
        
        
        for(int i=0;i<etrait._include.size();i++)
            etrait._probe_name_map.insert(pair<string, int>(etrait._epi_prbID[etrait._include[i]],etrait._include[i]));
        
        for(int i=0;i<esdata._include.size();i++)
            esdata._probe_name_map.insert(pair<string, int>(esdata._epi_prbID[esdata._include[i]],esdata._include[i]));
        
        cout<<endl<<"Performing interanl test..... "<<endl;
        vector<long> out_probid;
        vector<string> bxy;
        vector<string> sexy;
        vector<string> pxy;
        vector<string> bgwas;
        vector<string> segwas;
        vector<string> beqtl;
        vector<string> seeqtl;
        vector<string> pgwas;
        vector<string> peqtl;
        vector<string> rsid;
        vector<string> rschr;
        vector<string> rsbp;
        vector<string> rsfreq;
        
        vector<string> rsa1;
        vector<string> rsa2;
        vector<string> prb1;
        vector<string> nsnp_test1;
        
        vector<string> etrait_id;
        vector<string> etrait_bp;
        vector<string> etrait_gene;
        vector<string> etrait_chr;
        
        eSMRrlt erlt;
        read_esmr(&erlt,smrFileName);
        
        float progr0=0.0 , progr1;
        progress_print(progr0);
        
        for( int lid=0;lid<erlt.lineNum;lid++)
        {
            progr1=1.0*lid/erlt.lineNum;
            if(progr1-progr0-0.05>1e-6 || lid+1==erlt.lineNum)
            {
                if(lid+1==erlt.lineNum) progr1=1.0;
                progress_print(progr1);
                progr0=progr1;
            }

            
            string traitname=erlt._Outco_id[lid];
            
            map<string, int>::iterator iter;
            iter =  etrait._probe_name_map.find(traitname);
            if (iter == etrait._probe_name_map.end()) throw("Could not find probe "+traitname+" in "+eqtlFileName+" dataset, please check!");
            uint32_t ii = iter->second;
            
            gwasData gdata;
            gdata.allele_1.resize(etrait._esi_include.size());
            gdata.allele_2.resize(etrait._esi_include.size());
            gdata.byz=(double*)malloc(etrait._esi_include.size()*sizeof(double));
            gdata.seyz=(double*)malloc(etrait._esi_include.size()*sizeof(double));
            gdata.freq=(double*)malloc(etrait._esi_include.size()*sizeof(double));
            gdata.pvalue=(double*)malloc(etrait._esi_include.size()*sizeof(double));
            gdata.splSize=(uint32_t*)malloc(etrait._esi_include.size()*sizeof(uint32_t));
            
            
          
           // cout<<"\nPerforming analysis of eTrait [ "+traitname+" ]..."<<endl;
            memset(gdata.byz,0,etrait._esi_include.size()*sizeof(double));
            memset(gdata.seyz,0,etrait._esi_include.size()*sizeof(double));
            gdata._include.clear();
            gdata.snpName.clear();
            int count=0;
            if(etrait._rowid.empty())
            {
                for (int j = 0; j<etrait._esi_include.size(); j++)
                {
                    if (abs(etrait._bxz[ii][etrait._esi_include[j]] + 9) > 1e-6)
                    {
                        gdata.byz[count]=etrait._bxz[ii][etrait._esi_include[j]];
                        gdata.seyz[count]=etrait._sexz[ii][etrait._esi_include[j]];
                        gdata.snpName.push_back(etrait._esi_rs[etrait._esi_include[j]]);
                        gdata.allele_1[count]=etrait._esi_allele1[etrait._esi_include[j]];
                        gdata.allele_2[count]=etrait._esi_allele2[etrait._esi_include[j]];
                        gdata._include.push_back(etrait._esi_include[j]); // row id selected
                        count++;
                    }
                }
            }
            else
            {
                uint64_t beta_start=etrait._cols[ii<<1];
                uint64_t se_start=etrait._cols[1+(ii<<1)];
                uint64_t numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=etrait._rowid[beta_start+j];
                    if(binary_search(etrait._esi_include.begin(), etrait._esi_include.end(), ge_rowid))
                    {
                        gdata.byz[count]=etrait._val[beta_start+j];
                        gdata.seyz[count]=etrait._val[se_start+j];
                        gdata.snpName.push_back(etrait._esi_rs[ge_rowid]);
                        gdata.allele_1[count]=etrait._esi_allele1[ge_rowid];
                        gdata.allele_2[count]=etrait._esi_allele2[ge_rowid];
                        gdata._include.push_back(ge_rowid);
                        count++;
                    }
                    
                }
            }
            if(gdata.snpName.size()< m_hetero)
            {
                cout<<gdata.snpNum<<" common SNPs (less than parameter m_hetero: "+atos(m_hetero)+" ) are included from eTrait [ "+traitname+" ] summary."<<endl;
                continue;
            }
            gdata.snpNum=gdata.snpName.size();
            //cout<<gdata.snpNum<<" common SNPs are included from eTrait [ "+traitname+" ] summary."<<endl;
            
            vector<double> bxz;
            vector<double> sexz;
            vector<uint32_t> curId;
            vector<string> eName;
            vector<int> snpchrom;
            
            vector<string> allele1;
            vector<string> allele2;
            vector<uint32_t> bpsnp;
            vector<double> freq;
            
            vector<double> byz;
            vector<double> seyz;
            VectorXd zsxz;
            
            vector<int> sn_ids;
            
            VectorXd _byz;
            VectorXd _seyz;
            VectorXd _bxz;
            VectorXd _sexz;
            VectorXd _zsxz;
            MatrixXd _X;
            MatrixXd _LD;
            VectorXd ld_v;
            MatrixXd _LD_heidi;
            MatrixXd _X_heidi;
            
           
            string expo_name=erlt._Expo_id[lid];
            string refSNP=erlt._snp_rs[lid];
            map<string, int>::iterator iter1;
            iter1 =  esdata._probe_name_map.find(expo_name);
            if (iter1 == esdata._probe_name_map.end())
            {
                cout<<"Warning: Could not find probe "+expo_name+" in "+eqtlFileName2+" dataset, please check!"<<endl;
                continue;
            }
            int i = iter1->second;           
                
                //extract info from eqtl summary and gwas summary
                bxz.clear();
                sexz.clear();
                curId.clear(); // is the idxes of bfile._include not the values of
                eName.clear();
                snpchrom.clear();
                byz.clear();
                seyz.clear();
                allele1.clear();
                allele2.clear();
                bpsnp.clear();
                freq.clear();
                long maxid =-9;
                int probebp=esdata._epi_bp[i];
                int probechr=esdata._epi_chr[i];
                if(esdata._rowid.empty())
                {
                    for (int j = 0; j<esdata._esi_include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == etrait._esi_include.size()
                    {
                        if (abs(esdata._bxz[i][j] + 9) > 1e-6)
                        {
                            int snpbp=esdata._esi_bp[j];
                            int snpchr=esdata._esi_chr[j];
                            
                            int etrait_rid=etrait._esi_include[j];
                            long pos=find(gdata._include.begin(), gdata._include.end(), etrait_rid)-gdata._include.begin();
                            
                            if(snpchr==probechr && ABS(probebp-snpbp)<=cis_itvl && pos!=gdata._include.size())
                            {
                                bxz.push_back(esdata._bxz[i][j]);
                                sexz.push_back(esdata._sexz[i][j]);
                                byz.push_back(gdata.byz[pos]);
                                seyz.push_back(gdata.seyz[pos]);
                                curId.push_back(j);
                                eName.push_back(esdata._esi_rs[j]);
                                snpchrom.push_back(esdata._esi_chr[j]);
                                if(esdata._esi_rs[j]==refSNP) maxid=(eName.size()-1);
                                freq.push_back(bdata._mu[bdata._include[j]]/2);
                                allele1.push_back(esdata._esi_allele1[j]);
                                allele2.push_back(esdata._esi_allele2[j]);
                                bpsnp.push_back(esdata._esi_bp[j]);
                            }
                        }
                    }
                    
                }
                else{
                    uint64_t beta_start=esdata._cols[i<<1];
                    uint64_t se_start=esdata._cols[1+(i<<1)];
                    uint64_t numsnps=se_start-beta_start;
                    for(int j=0;j<numsnps;j++)
                    {
                        int ge_rowid=esdata._rowid[beta_start+j];
                        int snpbp=esdata._esi_bp[ge_rowid];
                        int snpchr=esdata._esi_chr[ge_rowid];
                        
                        int etrait_rid=etrait._esi_include[ge_rowid];
                        long pos=find(gdata._include.begin(), gdata._include.end(), etrait_rid)-gdata._include.begin();
                        
                        if(snpchr==probechr && ABS(probebp-snpbp)<=cis_itvl && pos!=gdata._include.size())
                        {
                            bxz.push_back(esdata._val[beta_start+j]);
                            sexz.push_back(esdata._val[se_start+j]);
                            byz.push_back(gdata.byz[pos]);
                            seyz.push_back(gdata.seyz[pos]);
                            curId.push_back(ge_rowid);
                            eName.push_back(esdata._esi_rs[ge_rowid]);
                            snpchrom.push_back(esdata._esi_chr[ge_rowid]);
                            if(esdata._esi_rs[ge_rowid]==refSNP) maxid=(eName.size()-1);
                            allele1.push_back(esdata._esi_allele1[ge_rowid]);
                            allele2.push_back(esdata._esi_allele2[ge_rowid]);
                            bpsnp.push_back(esdata._esi_bp[ge_rowid]);
                            freq.push_back(bdata._mu[bdata._include[ge_rowid]]/2);
                        }
                    }
                }
                if(maxid==-9) continue; //heidi SNP is not in selected SNPs
                if (bxz.size() == 0) continue;
                
                Map<VectorXd> ei_bxz(&bxz[0],bxz.size());
                Map<VectorXd> ei_sexz(&sexz[0],sexz.size());
                
                zsxz=ei_bxz.array()/ei_sexz.array();
                double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
            
            
                
                double bxy_val = byz[maxid] / bxz[maxid];
                double sexy_val = sqrt((seyz[maxid] * seyz[maxid] * bxz[maxid] * bxz[maxid] + sexz[maxid] * sexz[maxid] * byz[maxid] * byz[maxid]) / (bxz[maxid] * bxz[maxid] * bxz[maxid] * bxz[maxid]));
                double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
                
                
                double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
                
                double chisqyz = byz[maxid] / seyz[maxid];
                double pyz_val = pchisq(chisqyz*chisqyz, 1);
                
                out_probid.push_back(i);
                bxy.push_back(dtosf(bxy_val));
                sexy.push_back(dtosf(sexy_val));
                pxy.push_back(dtos(pxy_val));
                bgwas.push_back(dtosf(byz[maxid]));
                segwas.push_back(dtosf(seyz[maxid]));
                beqtl.push_back(dtosf(bxz[maxid]));
                seeqtl.push_back(dtosf(sexz[maxid]));
                pgwas.push_back(dtos(pyz_val));
                peqtl.push_back(dtos(pxz_val));
                rsid.push_back(eName[maxid]);
                rschr.push_back(atos(snpchrom[maxid]));
                rsbp.push_back(itos(bpsnp[maxid]));
                rsfreq.push_back(dtosf(freq[maxid]));
                
                rsa1.push_back(allele1[maxid]);
                rsa2.push_back(allele2[maxid]);
            
            etrait_id.push_back(etrait._epi_prbID[ii]);
            etrait_chr.push_back(atos(etrait._epi_chr[ii]));
            etrait_gene.push_back(etrait._epi_gene[ii]);
            etrait_bp.push_back(atos(etrait._epi_bp[ii]));
            
                    //extract info from reference
                    make_XMat(&bdata,curId, _X); //_X: one row one individual, one column one SNP
                    ld_calc_o2m(ld_v,maxid,_X);
                    
            
                    
                    sn_ids.clear(); //increase order
                    if(abs(ld_top-1)<1e-6) get_square_idxes(sn_ids,zsxz,threshold);
                    else get_square_ldpruning_idxes(sn_ids,zsxz,threshold,ld_v, maxid,ld_top);
                    
                    if(sn_ids.size() < m_hetero)
                    {
                        prb1.push_back(string("NA"));
                        nsnp_test1.push_back( string("NA"));
                        continue;
                    }
            
                    _byz.resize(sn_ids.size());
                    _seyz.resize(sn_ids.size());
                    _bxz.resize(sn_ids.size());
                    _sexz.resize(sn_ids.size());
                    _zsxz.resize(sn_ids.size());
                    _X_heidi.resize(_X.rows(), sn_ids.size());
                    
                    #pragma omp parallel for
                    for(int j=0;j<sn_ids.size();j++)
                    {
                        _byz[j]=byz[sn_ids[j]];
                        _seyz[j]=seyz[sn_ids[j]];
                        _bxz[j]=bxz[sn_ids[j]];
                        _sexz[j]=sexz[sn_ids[j]];
                        _zsxz[j]=zsxz[sn_ids[j]];
                        _X_heidi.col(j)=_X.col(sn_ids[j]);
                    }
                    _X.resize(0,0);
                    cor_calc(_LD_heidi, _X_heidi);
                    
            
                    _X_heidi.resize(0,0);
                    
                    long nsnp = sn_ids.size();
                    double pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, _LD_heidi, &nsnp);
                    
                    prb1.push_back(dtos(pdev));
                    if(nsnp>0) nsnp_test1.push_back(itos(nsnp));
                    else nsnp_test1.push_back(string("NA"));
            
            free_gwas_data( &gdata);
            
        }
        
        if(bxy.size()>0)
        {
            string smrfile = string(outFileName)+".smr";
            ofstream smr(smrfile.c_str());
            if (!smr) throw ("Error: can not open the file " + smrfile + " to save!");
            
            smr << "Expo_ID" <<'\t'<< "Expo_Chr" <<'\t' << "Expo_Gene"  << '\t' << "Expo_bp" << '\t'<< "Outco_ID" <<'\t'<< "Outco_Chr" <<'\t' << "Outco_Gene"  << '\t' << "Outco_bp" << '\t'<< "SNP"<< '\t' << "SNP_Chr"<< '\t' << "SNP_bp"<< '\t' << "A1"<< '\t'<< "A2"<< '\t'<<"Freq"<<'\t'<<"b_Outco"<<'\t'<<"se_Outco"<<'\t'<< "p_Outco" << '\t'<<"b_Expo"<<'\t'<<"se_Expo"<<'\t'<< "p_Expo" << '\t'<< "b_SMR" << '\t'<< "se_SMR"<< '\t' << "p_SMR" << "\t"<< "p_HET"<< "\t" << "nsnp" << '\n';
            
            for (int i = 0;i <bxy.size(); i++) {
                smr<<esdata._epi_prbID[out_probid[i]]<<'\t'<<esdata._epi_chr[out_probid[i]]<<'\t'<<esdata._epi_gene[out_probid[i]]<<'\t'<<esdata._epi_bp[out_probid[i]]<<'\t'<<etrait_id[i]<<'\t'<<etrait_chr[i]<<'\t'<<etrait_gene[i]<<'\t'<<etrait_bp[i]<<'\t'<<rsid[i]<<'\t'<<rschr[i]<<'\t'<<rsbp[i]<<'\t'<<rsa1[i]<<'\t'<<rsa2[i]<<'\t'<<rsfreq[i]<<'\t'<<bgwas[i]<<'\t'<<segwas[i]<<'\t'<<pgwas[i]<<'\t'<<beqtl[i]<<'\t'<<seeqtl[i]<<'\t'<<peqtl[i]<<'\t'<<bxy[i]<<'\t'<<sexy[i]<<'\t'<<pxy[i]<<'\t'<<prb1[i]<<'\t'<<nsnp_test1[i]<<'\n';
            }
            cout<<"Internal Test finished.\nInternal Test results of "<<bxy.size()<<" probes have been saved in the file [" + smrfile + "]."<<endl;
            smr.close();
            
        }else cout<<"Internal Test finished.\nInternal Test results of "<<bxy.size()<<" probes have been saved."<<endl;
        


    }
    
    void make_cojo(char* outFileName,char* eqtlFileName, char* snplstName,char* snplst2exclde, char* problstName, char* problst2exclde, char* genelistName, bool bFlag)
    {
        eqtlInfo eqtlinfo;
        cout<<endl<<"Reading eQTL summary data..."<<endl;
        if(eqtlFileName != NULL)
        {
            read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&eqtlinfo, snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&eqtlinfo, snplst2exclde);
            read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
            if(problstName != NULL) extract_prob(&eqtlinfo, problstName);
            if(genelistName != NULL) extract_prob_by_gene(&eqtlinfo, genelistName);
            if(problst2exclde != NULL) exclude_prob(&eqtlinfo, problst2exclde);
            read_besdfile(&eqtlinfo, string(eqtlFileName)+".besd");
            if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
            {
                printf("No data included from %s under current condition.\n",eqtlFileName);
                exit(EXIT_FAILURE);
            }
            
            
        }
        else throw ("Error: please input the eQTL summary information for the eQTL data files by the option --beqtl-summary.");
        
        vector<int> out_esi_id;
        vector<int> out_epi_id;
        vector<float> out_beta;
        vector<float> out_se;
        vector<double> out_pval;
        if(eqtlinfo._valNum==0)
        {
            for(uint32_t i=0;i<eqtlinfo._probNum;i++)
            {
                for(uint32_t j=0;j<eqtlinfo._snpNum;j++)
                {
                    double beta=eqtlinfo._bxz[i][j];
                    double se=eqtlinfo._sexz[i][j];
                    if(ABS(se+9)<1e-6) continue;
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                        out_esi_id.push_back(j);
                        out_epi_id.push_back(i);
                        out_beta.push_back(beta);
                        out_se.push_back(se);
                        out_pval.push_back(pxz);
                   
                }
            }
        }
        else
        {
            if(eqtlinfo._val.size()==0)
            {
                throw ("Error: No data extracted from the input, please check.");
            }
            
            for(uint32_t i=0;i<eqtlinfo._probNum;i++)
            {
                uint64_t proid=eqtlinfo._include[i];
                uint64_t pos=eqtlinfo._cols[proid<<1];
                uint64_t pos1=eqtlinfo._cols[(proid<<1)+1];
                uint64_t num=pos1-pos;
                for(int j=0;j<num;j++)
                {
                    double beta=eqtlinfo._val[pos+j];
                    double se=eqtlinfo._val[pos+j+num];
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                    
                        out_esi_id.push_back(eqtlinfo._rowid[pos+j]);
                        out_epi_id.push_back(i);
                        out_beta.push_back(beta);
                        out_se.push_back(se);
                        out_pval.push_back(pxz);
                   
                }
            }
        }
        
        string smrfile = string(outFileName)+".raw";
        ofstream smr(smrfile.c_str());
        if (!smr) throw ("Error: can not open the fam file " + smrfile + " to save!");
        
        smr << "ProbeID"<<'\t' << "SNP" <<'\t'<< "Chr" <<'\t' << "A1" << '\t'<< "A2"<< '\t' << "b"<<'\t'<< "SE" << '\t'<<"p"<<'\n';
        
        for (int i = 0;i <out_esi_id.size(); i++) {
            smr<<eqtlinfo._epi_prbID[out_epi_id[i]]<<'\t'<<eqtlinfo._esi_rs[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_chr[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele1[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele2[out_esi_id[i]]<<'\t'<<out_beta[i]<<'\t'<<out_se[i]<<'\t'<<out_pval[i]<< '\n';
        }
        
        smr.close();
        cout<<"Extracted results of "<<out_esi_id.size()<<" items have been saved in the plaint text file [" + smrfile + "]."<<endl;
        
    }
    
    void standardization(char* outFileName, char* eqtlFileName,bool bFlag,char* freqName, char* vpFileName)
    {
        eqtlInfo esdata;
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(freqName==NULL && vpFileName==NULL) throw("Error: please input feq data or variance data for standardisation by the flag --freq or --probe-var.");
        read_esifile(&esdata, string(eqtlFileName)+".esi");
        read_epifile(&esdata, string(eqtlFileName)+".epi");
         read_besdfile(&esdata, string(eqtlFileName)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included from %s under current condition.\n",eqtlFileName);
            exit(EXIT_FAILURE);
        }
        
        if(vpFileName!=NULL)
        {
            
            read_probevarfile(&esdata, string(vpFileName));
            for(int i=0;i<esdata._probNum;i++)
            {
                double prbvar_sqrt=sqrt(esdata._epi_var[i]);
                if(esdata._rowid.empty())
                {
                    for (int j = 0; j<esdata._esi_include.size(); j++)
                    {
                        if (abs(esdata._bxz[i][j] + 9) > 1e-6)
                        {
                            
                            float bxz=esdata._bxz[i][j];
                            float sexz=esdata._sexz[i][j];
                            esdata._bxz[i][j]=bxz/prbvar_sqrt;
                            esdata._sexz[i][j]=sexz/prbvar_sqrt;
                        }
                    }
                    
                }
                else{
                    uint64_t beta_start=esdata._cols[i<<1];
                    uint64_t se_start=esdata._cols[1+(i<<1)];
                    uint64_t numsnps=se_start-beta_start;
                    for(uint64_t j=0;j<numsnps;j++)
                    {
                        float bxz=esdata._val[beta_start+j];
                        float sexz=esdata._val[se_start+j];
                        esdata._val[beta_start+j]=bxz/prbvar_sqrt;
                        esdata._val[se_start+j]=sexz/prbvar_sqrt;
                    }
                }
                
            }

        }
        else if(freqName!=NULL)
        {
            
            int n=read_frqfile(&esdata, string(freqName));
            for(int i=0;i<esdata._probNum;i++)
            {
                if(esdata._rowid.empty())
                {
                    for (int j = 0; j<esdata._esi_include.size(); j++)
                    {
                        if (abs(esdata._bxz[i][j] + 9) > 1e-6)
                        {
                            
                            float bxz=esdata._bxz[i][j];
                            float sexz=esdata._sexz[i][j];
                            float p=esdata._esi_maf[j];
                            float z=bxz/sexz;
                            float b=z/sqrt(2*p*(1-p)*(n+z*z));
                            float se=1/sqrt(2*p*(1-p)*(n+z*z));
                            esdata._bxz[i][j]=b;
                            esdata._sexz[i][j]=se;
                        }
                    }
                    
                }
                else{
                    uint64_t beta_start=esdata._cols[i<<1];
                    uint64_t se_start=esdata._cols[1+(i<<1)];
                    uint64_t numsnps=se_start-beta_start;
                    for(uint64_t j=0;j<numsnps;j++)
                    {
                        uint64_t ge_rowid=esdata._rowid[beta_start+j];
                        float bxz=esdata._val[beta_start+j];
                        float sexz=esdata._val[se_start+j];
                        float p=esdata._esi_maf[ge_rowid];
                        float z=bxz/sexz;
                        float b=z/sqrt(2*p*(1-p)*(n+z*z));
                        float se=1/sqrt(2*p*(1-p)*(n+z*z));
                        esdata._val[beta_start+j]=b;
                        esdata._val[se_start+j]=se;
                    }
                }
                
            }

        }
       write_besd(outFileName, &esdata);
    }
   

    void lookup(char* outFileName,char* eqtlFileName, char* snplstName, char* problstName,char* genelistName, double plookup,bool bFlag, int chr,  int prbchr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs, char* prbname, char* fromprbname, char* toprbname,int snpWind, int prbWind,char* genename,int fromsnpkb, int tosnpkb, int fromprbkb, int toprbkb, bool snpwindFlag, bool prbwindFlag,bool cis_flag, int cis_itvl)
    {
        string logstr;
        int flag4chr=0;
        if(chr!=0) flag4chr++;
        if(prbchr!=0 || snpchr!=0) flag4chr++;
        if(flag4chr==2)
        {
            logstr="WARNING: --chr is not surpposed to use together with --probe-chr or --snp-chr. --chr will be disabled.\n";
            chr=0;
            fputs(logstr.c_str(), stdout);
        }       
        
        eqtlInfo eqtlinfo;
        cout<<endl<<"Reading eQTL summary data..."<<endl;
        if(eqtlFileName != NULL)
        {
            read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
            epi_man(&eqtlinfo, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
            
            read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
            esi_man(&eqtlinfo, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,  cis_itvl, prbname);
           
            
           read_besdfile(&eqtlinfo, string(eqtlFileName)+".besd");
            if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
            {
                printf("No data included from %s under current condition.\n",eqtlFileName);
                exit(EXIT_FAILURE);
            }

            
        }
        else throw ("Error: please input the eQTL summary information for the eQTL data files by the option --beqtl-summary.\n");
        
        vector<int> out_esi_id;
        vector<int> out_epi_id;
        vector<float> out_beta;
        vector<float> out_se;
        vector<double> out_pval;
        if(eqtlinfo._valNum==0)
        {
            for(uint32_t i=0;i<eqtlinfo._probNum;i++)
            {
                for(uint32_t j=0;j<eqtlinfo._snpNum;j++)
                {
                    double beta=eqtlinfo._bxz[i][j];
                    double se=eqtlinfo._sexz[i][j];
                    if(ABS(se+9)<1e-6) continue;
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                    if(pxz<=plookup)
                    {
                        out_esi_id.push_back(j);
                        out_epi_id.push_back(i);
                        out_beta.push_back(beta);
                        out_se.push_back(se);
                        out_pval.push_back(pxz);
                    }
                }
            }
        }
        else
        {
            if(eqtlinfo._val.size()==0)
            {
                throw ("Error: No data extracted from the input, please check.\n");
            }
            
            for(uint32_t i=0;i<eqtlinfo._probNum;i++)
            {
                uint64_t proid=eqtlinfo._include[i];
                uint64_t pos=eqtlinfo._cols[proid<<1];
                uint64_t pos1=eqtlinfo._cols[(proid<<1)+1];
                uint64_t num=pos1-pos;
                for(int j=0;j<num;j++)
                {
                    double beta=eqtlinfo._val[pos+j];
                    double se=eqtlinfo._val[pos+j+num];
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                    if(pxz<=plookup)
                    {
                        out_esi_id.push_back(eqtlinfo._rowid[pos+j]);
                        out_epi_id.push_back(i);
                        out_beta.push_back(beta);
                        out_se.push_back(se);
                        out_pval.push_back(pxz);
                    }
                }
            }
        }
        
        string smrfile = string(outFileName)+".txt";
        ofstream smr(smrfile.c_str());
        if (!smr) throw ("Error: can not open the fam file " + smrfile + " to save!");
        
        smr << "SNP" <<'\t'<< "Chr" <<'\t' << "BP"  << '\t' << "A1" << '\t'<< "A2"<< '\t' << "Probe"<< '\t' << "Probe_Chr"<< '\t'<< "Probe_bp"<< '\t'<<"Gene"<<'\t'<<"Orientation"<<'\t'<<"b"<<'\t'<< "SE" << '\t'<<"p"<<'\n';
        
        for (int i = 0;i <out_esi_id.size(); i++) {
            smr<<eqtlinfo._esi_rs[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_chr[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_bp[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele1[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele2[out_esi_id[i]]<<'\t'<<eqtlinfo._epi_prbID[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_chr[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_bp[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_gene[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_orien[out_epi_id[i]]<<'\t'<<out_beta[i]<<'\t'<<out_se[i]<<'\t'<<out_pval[i]<< '\n';
        }
        
        smr.close();
        cout<<"Extracted results of "<<out_esi_id.size()<<" SNPs have been saved in the file [" + smrfile + "]."<<endl;
        
    }
    
    void rm_unmatched_snp(gwasData* gdata, eqtlInfo* esdata)
    {
        cout<<"Allele check..."<<endl;
        // get the common SNPs
        vector<string> slctSNPs;
        vector<int> id_unmatched_esd;
        vector<int> id_unmatched_gwas;
        vector<int> gdId;
        vector<int> edId;
        edId.clear();
        int pre_size=esdata->_esi_include.size();
        vector<string> essnp;
        if(esdata->_esi_include.size()<esdata->_snpNum )
        {
            for(int i=0;i<esdata->_esi_include.size();i++) essnp.push_back(esdata->_esi_rs[esdata->_esi_include[i]]);
            StrFunc::match_only(essnp, gdata->snpName, gdId);
            if(gdId.empty()) throw("Error: no common SNPs found.");
            for(int i=0;i<gdId.size();i++) slctSNPs.push_back(gdata->snpName[gdId[i]]);
        }else
        {
            StrFunc::match_only(esdata->_esi_rs, gdata->snpName, gdId);
            if(gdId.empty()) throw("Error: no common SNPs found.");
            for(int i=0;i<gdId.size();i++) slctSNPs.push_back(gdata->snpName[gdId[i]]);
        }
        
        
        //alleles check
        StrFunc::match(slctSNPs, esdata->_esi_rs, edId);
        id_unmatched_esd.clear();
        id_unmatched_gwas.clear();
        for (int i = 0; i<edId.size(); i++)
        {
            string ga1, ga2, ea1, ea2;
            
            ga1 = gdata->allele_1[gdId[i]];
            ga2 = gdata->allele_2[gdId[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
            // use the allele in eQTL summary data as the reference allele. so we won't get the whole besd into memroy
            if( ea1 == ga1 && ea2 == ga2)
            {
                
            }
            else if(ea1 == ga2 && ea2 == ga1)
            {
                gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                
            }else {
                id_unmatched_esd.push_back(edId[i]);
                 id_unmatched_gwas.push_back(gdId[i]);
            }
        }
        if(id_unmatched_esd.size()>0)
        {
            vector<int> tmpvec;
            esdata->_esi_include.swap(tmpvec);
            set_complement(id_unmatched_esd,tmpvec, esdata->_esi_include);
        }
        if(id_unmatched_gwas.size()>0)
        {
            vector<int> tmpvec;
            gdata->_include.swap(tmpvec);
            set_complement(id_unmatched_gwas,tmpvec, gdata->_include);
        }
        cout<<id_unmatched_esd.size()<<" SNPs failed in allele check and have been excluded. Total "<<gdata->_include.size()<<" SNPs left in GWAS summary dataset and "<<esdata->_esi_include.size()<<" SNPs left in eQTL summary dtaset."<<endl;
       
    }

    void read_gene_anno(char* geneAnnoName,vector<int> &chr, vector<string> &genename,vector<int> &start,vector<int> &end)
    {
       
        ifstream flptr;
        flptr.open(geneAnnoName);
        if (!flptr) throw ("Error: can not open the file [" + string(geneAnnoName) + "] to read.");
      
        cout << "Reading gene annotation information from [" + string(geneAnnoName) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        
        while(!flptr.eof() )
        {
            string tmpStr;
            flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                istringstream iss(buf);
                iss>>tmpStr; //chr
                if(tmpStr=="X" || tmpStr=="x") chr.push_back(23);
                else if(tmpStr=="Y" || tmpStr=="y") chr.push_back(24);
                else chr.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr; // start
                start.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr; // end
                end.push_back(atoi(tmpStr.c_str()));
                iss>>tmpStr; // gene
                genename.push_back(tmpStr.c_str());
                lineNum++;
                
            }
        }
        cout << lineNum << " gene annotation infomation to be included from [" + string(geneAnnoName) + "]." << endl;
        flptr.close();

    }

    void plot(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero , char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, char* refSNP, bool heidioffFlag, int cis_itvl, char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag, char* geneAnnoName)
    {
        
        setNbThreads(thread_num);
        if(!prbwindFlag) throw("Error: please input probe window by the flag --probe-wind.");
        if(!heidioffFlag && bFileName == NULL ) throw("Error: please input Plink file by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data  by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data  by the flag --eqtl-summary.");
        if(geneAnnoName==NULL) throw("Error: please input gene annotation file by the flag --gene-list.");
        
        vector<int> gene_anno_chr;
        vector<string> gene_anno_genename;
        vector<int> gene_anno_start;
        vector<int> gene_anno_end;
        read_gene_anno(geneAnnoName,gene_anno_chr, gene_anno_genename,gene_anno_start,gene_anno_end);
        map<string, int> gene_anno_map;
        map<string, int>::iterator iter;
        for(int i=0;i<gene_anno_genename.size();i++) gene_anno_map.insert(pair<string,int>(gene_anno_genename[i], i));
        
        eqtlInfo prbhead;
        read_epifile(&prbhead, string(eqtlFileName)+".epi");
        if(problstName != NULL || genelistName != NULL)
        {
            if(problstName != NULL) extract_prob(&prbhead, problstName);
            if(genelistName != NULL) extract_prob_by_gene(&prbhead, genelistName);
        }
        else if(prbname!=NULL)
        {
            extract_eqtl_single_probe(&prbhead, prbname);
        }
  
        cout<<"\nThere would create "<<prbhead._include.size()<<" plot files..."<<endl;
        if(prbhead._include.size()>1024) cout<<"WARNING: Too many files!!! We strongly recommend using --probe or --extract-probe."<<endl;
        
        for(int plotid=0;plotid<prbhead._include.size();plotid++)
        {
            string plotprbname=prbhead._epi_prbID[prbhead._include[plotid]];
          
            gwasData gdata_;
            eqtlInfo esdata_;
            bInfo bdata;
            gwasData gdata;
            eqtlInfo esdata;
            double threshold= chi_val(1,p_hetero);
            bool heidiFlag=false;
            if(refSNP!=NULL) heidiFlag=true;
            
            
            
            cout<<"\nRetrieving SMR results for plot..."<<endl;
            read_gwas_data( &gdata, gwasFileName);
            read_esifile(&esdata, string(eqtlFileName)+".esi");
            esi_man(&esdata, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, false,  cis_itvl, plotprbname.c_str());
            if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
            if(!heidioffFlag)
            {
                read_famfile(&bdata, string(bFileName)+".fam");
                if(indilstName != NULL) keep_indi(&bdata,indilstName);
                if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
                read_bimfile(&bdata, string(bFileName)+".bim");
                if(snplstName != NULL) extract_snp(&bdata, snplstName);
                if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
                allele_check(&bdata, &gdata, &esdata);
                read_bedfile(&bdata, string(bFileName)+".bed");
                if (bdata._mu.empty()) calcu_mu(&bdata);
                if (maf > 0)
                {
                    filter_snp_maf(&bdata, maf);
                    update_geIndx(&bdata, &gdata, &esdata);
                }
                
            }else
            {
                allele_check(&gdata, &esdata);
            }
            
            update_gwas(&gdata);
            cout<<"Reading eQTL summary data..."<<endl;
            read_epifile(&esdata, string(eqtlFileName)+".epi");
            extract_prob(&esdata, plotprbname, prbWind);
           
            read_besdfile(&esdata, string(eqtlFileName)+".besd");
            if(esdata._rowid.empty() && esdata._bxz.empty())
            {
                printf("No data included from %s under current condition.\n",eqtlFileName);
                exit(EXIT_FAILURE);
            }
            
            
            unsigned int probNum = esdata._probNum;
            
            vector<double> bxz;
            vector<double> sexz;
            vector<uint32_t> curId;
            vector<string> eName;
            vector<int> snpchrom;
            
            vector<string> allele1;
            vector<string> allele2;
            vector<uint32_t> bpsnp;
            vector<double> freq;
            
            vector<double> byz;
            vector<double> seyz;
            vector<double> pyz;
            VectorXd zsxz;
            
            vector<int> sn_ids;
            
            VectorXd _byz;
            VectorXd _seyz;
            VectorXd _bxz;
            VectorXd _sexz;
            VectorXd _zsxz;
            MatrixXd _X;
            MatrixXd _LD;
            VectorXd ld_v;
            MatrixXd _LD_heidi;
            MatrixXd _X_heidi;
            int cis_itvl_bp=cis_itvl*1000;
            
            //for plot
            string plotdir="";
            string plotnm=outFileName;
            for(long j=strlen(outFileName)-1;j>=0;j--)
                if(outFileName[j]=='/')
                {
                    plotdir=string(outFileName).substr(0,j+1);
                    plotnm=string(outFileName).substr(j+1,strlen(outFileName));
                    break;
                }
            if(plotdir=="") plotdir="./";
            
            plotdir=string(plotdir)+"plot";
            struct stat st = {0};
            if (stat(plotdir.c_str(), &st) == -1) {
#if defined _WIN64 || defined _WIN32
                _mkdir(plotdir.c_str());
#else
                mkdir(plotdir.c_str(), 0755);
#endif
            }
            string plot_path= string(plotdir)+"/"+plotnm+"."+plotprbname+".txt";
            
            //probe info
            vector<string> epi_out_id;
            vector<int> epi_out_chr;
            vector<int> epi_out_bp;
            vector<string> epi_out_gene;
            vector<char> epi_out_orien;
            vector<double> epi_out_pheidi;
            vector<double> epi_out_psmr;
            vector<int> epi_out_start;
            vector<int> epi_out_end;
            
            //eqtl ld info
            vector<int> ldprbid;
            vector<string> ldprb;
            vector<int> ldnperprb;
            vector<string> ldrs;
            vector<double> outld;
            ldnperprb.push_back(0);
            
            
            long idx=find(esdata._epi_prbID.begin(), esdata._epi_prbID.end(), plotprbname)-esdata._epi_prbID.begin();
            if(idx==esdata._epi_prbID.size())
            {
                string logstr="ERROR: Can't find probe "+string(plotprbname)+".\n";
                fputs(logstr.c_str(),stdout);
                exit(1);
            }
            int minBP=esdata._epi_bp[idx]-prbWind*1000>0?(esdata._epi_bp[idx]-prbWind*1000):0;
            int maxBP=esdata._epi_bp[idx]+prbWind*1000;
            int refprbchr=esdata._epi_chr[idx];
            
            for(int i=0;i<probNum;i++)
            {
                
                //extract info from eqtl summary and gwas summary
                bxz.clear();
                sexz.clear();
                curId.clear(); // is the idxes of bfile._include not the values of
                eName.clear();
                snpchrom.clear();
                byz.clear();
                seyz.clear();
                pyz.clear();
                allele1.clear();
                allele2.clear();
                bpsnp.clear();
                freq.clear();
                long maxid =-9;
                int probebp=esdata._epi_bp[i];
                int probechr=esdata._epi_chr[i];
                string probenm=esdata._epi_prbID[i];
                int tmpminBP=probebp;
                int tmpmaxBP=probebp;
                if(esdata._rowid.empty())
                {
                    for (int j = 0; j<bdata._include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
                    {
                        if (abs(esdata._bxz[i][j] + 9) > 1e-6)
                        {
                            int snpbp=esdata._esi_bp[j];
                            int snpchr=esdata._esi_chr[j];
                            if(snpchr==probechr && ABS(probebp-snpbp)<=cis_itvl_bp)
                            {
                                tmpminBP=snpbp<tmpminBP?snpbp:tmpminBP;
                                tmpmaxBP=snpbp>tmpmaxBP?snpbp:tmpmaxBP;
                                bxz.push_back(esdata._bxz[i][j]);
                                sexz.push_back(esdata._sexz[i][j]);
                                byz.push_back(gdata.byz[j]);
                                seyz.push_back(gdata.seyz[j]);
                                pyz.push_back(gdata.pvalue[j]);
                                curId.push_back(j);
                                eName.push_back(esdata._esi_rs[j]);
                                snpchrom.push_back(esdata._esi_chr[j]);
                                if(heidiFlag && esdata._esi_rs[j]==string(refSNP)) maxid=(eName.size()-1);
                                if(!heidioffFlag) //if heidi off , bfile is not necessary to read.
                                {
                                    freq.push_back(bdata._mu[bdata._include[j]]/2);
                                }
                                allele1.push_back(esdata._esi_allele1[j]);
                                allele2.push_back(esdata._esi_allele2[j]);
                                bpsnp.push_back(esdata._esi_bp[j]);
                            }
                        }
                    }
                    
                }
                else{
                    uint64_t beta_start=esdata._cols[i<<1];
                    uint64_t se_start=esdata._cols[1+(i<<1)];
                    uint64_t numsnps=se_start-beta_start;
                    for(int j=0;j<numsnps;j++)
                    {
                        int ge_rowid=esdata._rowid[beta_start+j];
                        int snpbp=esdata._esi_bp[ge_rowid];
                        int snpchr=esdata._esi_chr[ge_rowid];
                        
                        if(snpchr==probechr && ABS(probebp-snpbp)<=cis_itvl_bp)
                        {
                            tmpminBP=snpbp<tmpminBP?snpbp:tmpminBP;
                            tmpmaxBP=snpbp>tmpmaxBP?snpbp:tmpmaxBP;
                            bxz.push_back(esdata._val[beta_start+j]);
                            sexz.push_back(esdata._val[se_start+j]);
                            byz.push_back(gdata.byz[ge_rowid]);
                            seyz.push_back(gdata.seyz[ge_rowid]);
                            pyz.push_back(gdata.pvalue[ge_rowid]);
                            curId.push_back(ge_rowid);
                            eName.push_back(esdata._esi_rs[ge_rowid]);
                            snpchrom.push_back(esdata._esi_chr[ge_rowid]);
                            if(heidiFlag && esdata._esi_rs[ge_rowid]==string(refSNP)) maxid=(eName.size()-1);
                            allele1.push_back(esdata._esi_allele1[ge_rowid]);
                            allele2.push_back(esdata._esi_allele2[ge_rowid]);
                            bpsnp.push_back(esdata._esi_bp[ge_rowid]);
                            if(!heidioffFlag){
                                freq.push_back(bdata._mu[bdata._include[ge_rowid]]/2);
                            }
                        }
                    }
                }
                if(heidiFlag && maxid==-9) continue; //heidi SNP is not in selected SNPs
                if (bxz.size() == 0) continue;
                
                epi_out_id.push_back(esdata._epi_prbID[i]);
                epi_out_chr.push_back(probechr);
                epi_out_bp.push_back(probebp);
                epi_out_gene.push_back(esdata._epi_gene[i]);
                epi_out_orien.push_back(esdata._epi_orien[i]);
                epi_out_pheidi.push_back(-9);
                epi_out_psmr.push_back(-9);
                iter = gene_anno_map.find((esdata._epi_gene[i]));
                if (iter != gene_anno_map.end())
                                          {
                                              epi_out_start.push_back(gene_anno_start[iter->second]);
                                              epi_out_end.push_back(gene_anno_end[iter->second]);
                                          } else {
                                              epi_out_start.push_back(-9);
                                              epi_out_end.push_back(-9);
                                          }
               
                
                Map<VectorXd> ei_bxz(&bxz[0],bxz.size());
                Map<VectorXd> ei_sexz(&sexz[0],sexz.size());
                
                zsxz=ei_bxz.array()/ei_sexz.array();
                if(!heidiFlag) maxid=max_abs_id(zsxz);
                double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
                
                if(!heidiFlag && pxz_val>p_smr) continue;
                
                
                if(tmpminBP<minBP) minBP=tmpminBP;
                if(tmpmaxBP>maxBP) maxBP=tmpmaxBP;
                
                double bxy_val = byz[maxid] / bxz[maxid];
                double sexy_val = sqrt((seyz[maxid] * seyz[maxid] * bxz[maxid] * bxz[maxid] + sexz[maxid] * sexz[maxid] * byz[maxid] * byz[maxid]) / (bxz[maxid] * bxz[maxid] * bxz[maxid] * bxz[maxid]));
                double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
                
                
                double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
                
                epi_out_psmr[i]=pxy_val;
                
                
                make_XMat(&bdata,curId, _X);
                ld_calc_o2m(ld_v,maxid,_X);
                
                // out ld
                ldprbid.push_back(i);
                ldprb.push_back(probenm);
                ldnperprb.push_back((int)curId.size()+ldnperprb[ldnperprb.size()-1]);
                for(int jj=0;jj<curId.size();jj++){
                    ldrs.push_back(esdata._esi_rs[curId[jj]]);
                    outld.push_back(ld_v(jj));
                }
                
                
                sn_ids.clear(); //increase order
                if(abs(ld_top-1)<1e-6) get_square_idxes(sn_ids,zsxz,threshold);
                else get_square_ldpruning_idxes(sn_ids,zsxz,threshold,ld_v, maxid,ld_top);
                
                if(sn_ids.size() < m_hetero)   continue;
                
                
                _byz.resize(sn_ids.size());
                _seyz.resize(sn_ids.size());
                _bxz.resize(sn_ids.size());
                _sexz.resize(sn_ids.size());
                _zsxz.resize(sn_ids.size());
                _X_heidi.resize(_X.rows(), sn_ids.size());
                
#pragma omp parallel for
                for(int j=0;j<sn_ids.size();j++)
                {
                    _byz[j]=byz[sn_ids[j]];
                    _seyz[j]=seyz[sn_ids[j]];
                    _bxz[j]=bxz[sn_ids[j]];
                    _sexz[j]=sexz[sn_ids[j]];
                    _zsxz[j]=zsxz[sn_ids[j]];
                    _X_heidi.col(j)=_X.col(sn_ids[j]);
                }
                _X.resize(0,0);
                cor_calc(_LD_heidi, _X_heidi);
                
                
                
                _X_heidi.resize(0,0);
                
                long nsnp = sn_ids.size();
                double pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, _LD_heidi, &nsnp);
                epi_out_pheidi[i]=pdev;
            }
            
            if(ldprbid.size()==0)
            {
                string logstr="No SMR results fetched.\n";
                fputs(logstr.c_str(),stdout);
                continue;
            }
            
            cout<<"\nRetrieving GWAS summary information and eQTL summary information for plot..."<<endl;
            read_gwas_data( &gdata_, gwasFileName);
            
            read_epifile(&esdata_, string(eqtlFileName)+".epi");
            extract_prob(&esdata_, plotprbname, prbWind);
            
            
            //get eQTL info in the region [minBP,maxBP]
            read_esifile(&esdata_, string(eqtlFileName)+".esi");
            vector<int> newIcld;
            for(int i=0;i<esdata_._esi_include.size();i++)
            {
                int tmpint=esdata_._esi_include[i];
                if( esdata_._esi_chr[tmpint]==refprbchr &&esdata_._esi_bp[tmpint]>=minBP && esdata_._esi_bp[tmpint]<=maxBP) newIcld.push_back(tmpint);
            }
            esdata_._esi_include.clear();
            esdata_._esi_include=newIcld;
            cout << esdata_._esi_include.size() << " SNPs are extracted from SNP BP: " +atos(minBP)+" bp to SNP BP: " + atos(maxBP) + " bp of eQTL summary dataset." << endl;
            //get RSs in the region [minBP,maxBP] to select gwas summary
            read_bimfile(&bdata, string(bFileName)+".bim");
            vector<string> bsnprs;
            vector<int> bsnpbp;
            for(int i=0;i<bdata._snp_num;i++)
                if(bdata._chr[i]==refprbchr && bdata._bp[i]<=maxBP && bdata._bp[i]>=minBP)
                {
                    bsnprs.push_back(bdata._snp_name[i]);
                    bsnpbp.push_back(bdata._bp[i]);
                }
            
            vector<int> idx1;
            gdata_.snpBp.clear();
            gdata_.snpBp.resize(gdata_.snpName.size());
            gdata_._include.clear();
            match(gdata_.snpName,bsnprs,idx1); //get gwas info
            for(int i=0;i<gdata_.snpName.size();i++)
            {
                if(idx1[i]!=-9)
                {
                    gdata_._include.push_back(i);
                    gdata_.snpBp[i]=bsnpbp[idx1[i]];
                }
                
            }
            cout << gdata_._include.size() << " SNPs are extracted from SNP BP: " +atos(minBP)+" bp to SNP BP: " + atos(maxBP) + " bp of GWAS summary dataset." << endl;
            
            update_gwas(&gdata_);
            rm_unmatched_snp(&gdata_, &esdata_); //allele check at the same time
            update_gwas(&gdata_);
            
            read_besdfile(&esdata_, string(eqtlFileName)+".besd");
            if(esdata_._rowid.empty() && esdata_._bxz.empty())
            {
                printf("No data included from %s under current condition.\n",eqtlFileName);
                exit(EXIT_FAILURE);
            }
            //SNP info, the union of gwas and eqtl
            vector<string> out_rs;
            vector<int> out_chr;
            vector<int> out_bp;
            vector<string> out_a1;
            vector<string> out_a2;
            map<string,int> snp_name_map;
            long mapsize=0;
            for(int i=0;i<esdata_._esi_include.size();i++)
            {
                snp_name_map.insert(pair<string, int>(esdata_._esi_rs[esdata_._esi_include[i]], mapsize));
                if (mapsize < snp_name_map.size()) {
                    out_rs.push_back(esdata_._esi_rs[esdata_._esi_include[i]]);
                    out_chr.push_back(esdata_._esi_chr[esdata_._esi_include[i]]);
                    out_bp.push_back(esdata_._esi_bp[esdata_._esi_include[i]]);
                    out_a1.push_back(esdata_._esi_allele1[esdata_._esi_include[i]]);
                    out_a2.push_back(esdata_._esi_allele2[esdata_._esi_include[i]]);
                    mapsize=snp_name_map.size();
                }
                
            }
            for(int i=0;i<gdata_._include.size();i++)
            {
                snp_name_map.insert(pair<string, int>(gdata_.snpName[gdata_._include[i]], mapsize));
                if (mapsize < snp_name_map.size()) {
                    out_rs.push_back(gdata_.snpName[gdata_._include[i]]);
                    out_chr.push_back(refprbchr);
                    out_bp.push_back(gdata_.snpBp[gdata_._include[i]]);
                    out_a1.push_back(gdata_.allele_1[gdata_._include[i]]);
                    out_a2.push_back(gdata_.allele_2[gdata_._include[i]]);
                    mapsize=snp_name_map.size();
                }
            }
            vector<int> bprank;
            getRank_norep(out_bp, bprank);
            vector<string> tmptmpstr;
            tmptmpstr.resize(out_rs.size());
            for(int i=0;i<out_rs.size();i++) tmptmpstr[bprank[i]]=out_rs[i];
            out_rs.swap(tmptmpstr);
            vector<int> tmptmpint;
            tmptmpint.resize(out_rs.size());
            for(int i=0;i<out_rs.size();i++) tmptmpint[bprank[i]]=out_bp[i];
            out_bp.swap(tmptmpint);
            vector<string> tmptmpchar;
            tmptmpchar.resize(out_rs.size());
            for(int i=0;i<out_rs.size();i++) tmptmpchar[bprank[i]]=out_a1[i];
            out_a1.swap(tmptmpchar);
            tmptmpchar.resize(out_rs.size());
            for(int i=0;i<out_rs.size();i++) tmptmpchar[bprank[i]]=out_a2[i];
            out_a2.swap(tmptmpchar);
            
            //gwas info
            vector<string> gwas_rs;
            vector<float> gwas_be;
            vector<float> gwas_se;
            for(int i=0;i<gdata_.snpNum;i++)
            {
                gwas_rs.push_back(gdata_.snpName[i]);
                gwas_be.push_back(gdata_.byz[i]);
                gwas_se.push_back(gdata_.seyz[i]);
            }
            //eqtl info
            
            
            vector<int> out_esi_id;
            vector<int> out_epi_id;
            vector<string> out_epi_name;
            vector<float> out_beta;
            vector<float> out_se;
            vector<double> out_pval;
            if(esdata_._valNum==0)
            {
                for(uint32_t ii=0;ii<ldprbid.size();ii++)
                {
                    int i=ldprbid[ii];
                    for(uint32_t j=0;j<esdata_._snpNum;j++)
                    {
                        double beta=esdata_._bxz[i][j];
                        double se=esdata_._sexz[i][j];
                        if(ABS(se+9)<1e-6) continue;
                        double zsxz=beta/se;
                        double pxz=pchisq(zsxz*zsxz, 1);
                        if(pxz<=1)
                        {
                            out_esi_id.push_back(j);
                            out_epi_id.push_back(i);
                            out_epi_name.push_back(esdata_._epi_prbID[i]);
                            out_beta.push_back(beta);
                            out_se.push_back(se);
                            out_pval.push_back(pxz);
                        }
                    }
                }
            }
            else
            {
                if(esdata_._val.size()==0)
                {
                    throw ("Error: No data extracted from the input, please check.\n");
                }
                
                for(uint32_t ii=0;ii<ldprbid.size();ii++)
                {
                    int i=ldprbid[ii];
                    uint64_t proid=esdata_._include[i];
                    uint64_t pos=esdata_._cols[proid<<1];
                    uint64_t pos1=esdata_._cols[(proid<<1)+1];
                    uint64_t num=pos1-pos;
                    for(int j=0;j<num;j++)
                    {
                        double beta=esdata_._val[pos+j];
                        double se=esdata_._val[pos+j+num];
                        double zsxz=beta/se;
                        double pxz=pchisq(zsxz*zsxz, 1);
                        if(pxz<=1)
                        {
                            out_esi_id.push_back(esdata_._rowid[pos+j]);
                            out_epi_id.push_back(i);
                            out_beta.push_back(beta);
                            out_se.push_back(se);
                            out_pval.push_back(pxz);
                        }
                    }
                }
            }
            
            vector<double> out_esi_ld;
            vector<string> out_esi_rs;
            vector<int> stend;
            stend.push_back(0);
            int curprid=out_epi_id[0];
            for(int i=0;i<out_esi_id.size();i++)
            {
                if(out_epi_id[i]!=curprid)
                {
                    stend.push_back(i);
                    curprid=out_epi_id[i];
                }
                out_esi_rs.push_back(esdata_._esi_rs[out_esi_id[i]]);
                out_esi_ld.push_back(-9);
            }
            stend.push_back((int)out_esi_id.size());
            
            
            for(uint32_t ii=0;ii<ldprbid.size();ii++)
            {
                vector<string> outrs;
                for(int j=stend[ii];j<stend[ii+1];j++) outrs.push_back(out_esi_rs[j]);
                vector<string> ldsnp;
                for(int j=ldnperprb[ii];j<ldnperprb[ii+1];j++) ldsnp.push_back(ldrs[j]);
                vector<int> idx;
                match(ldsnp,outrs,idx);
                for(int j=0;j<idx.size();j++) out_esi_ld[stend[ii]+idx[j]]=outld[ldnperprb[ii]+j];
            }
            cout<<"Total "<<out_esi_id.size()<<" eQTLs for "<<ldprbid.size()<<" probes are extracted."<<endl;
            free_gwas_data( &gdata_);
            
            
            
            
            FILE* plotfile=NULL;
            plotfile = fopen(plot_path.c_str(), "w");
            if (!(plotfile)) {
                printf("Open error %s\n", plot_path.c_str());
                exit(1);
            }
            
            string outstr="$probe "+atos(epi_out_id.size())+" "+plotprbname+'\n';
            if(fputs_checked(outstr.c_str(),plotfile))
            {
                printf("ERROR: in writing file %s .\n", plot_path.c_str());
                exit(EXIT_FAILURE);
            }
            for(int i=0;i<epi_out_id.size();i++)
            {
                outstr=epi_out_id[i]+' '+atos(epi_out_chr[i])+' '+atos(epi_out_bp[i])+' '+epi_out_gene[i]+' '+(epi_out_start[i]==-9?"NA":(epi_out_start[i]==23?"X":(epi_out_start[i]==24?"Y":(atos(epi_out_start[i])))))+' '+(epi_out_end[i]==-9?"NA":(epi_out_end[i]==23?"X":(epi_out_end[i]==24?"Y":(atos(epi_out_end[i])))))+' '+epi_out_orien[i]+' '+(epi_out_psmr[i]+9<1e-6?"NA":dtos(epi_out_psmr[i]))+' '+(epi_out_pheidi[i]+9<1e-6?"NA":dtos(epi_out_pheidi[i]))+'\n';
                if(fputs_checked(outstr.c_str(),plotfile))
                {
                    printf("ERROR: in writing file %s .\n", plot_path.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            outstr="$SNP "+atos(out_rs.size())+'\n';
            if(fputs_checked(outstr.c_str(),plotfile))
            {
                printf("ERROR: in writing file %s .\n", plot_path.c_str());
                exit(EXIT_FAILURE);
            }
            for(int i=0;i<out_rs.size();i++)
            {
                outstr=out_rs[i]+' '+atos(out_chr[i])+' '+atos(out_bp[i])+' '+out_a1[i]+' '+out_a2[i]+'\n';
                if(fputs_checked(outstr.c_str(),plotfile))
                {
                    printf("ERROR: in writing file %s .\n", plot_path.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            outstr="$GWAS "+atos(gwas_rs.size())+'\n';
            if(fputs_checked(outstr.c_str(),plotfile))
            {
                printf("ERROR: in writing file %s .\n", plot_path.c_str());
                exit(EXIT_FAILURE);
            }
            for(int i=0;i<gwas_rs.size();i++)
            {
                outstr=gwas_rs[i]+' '+atos(gwas_be[i])+' '+atos(gwas_se[i])+'\n';
                if(fputs_checked(outstr.c_str(),plotfile))
                {
                    printf("ERROR: in writing file %s .\n", plot_path.c_str());
                    exit(EXIT_FAILURE);
                }
            }
            outstr="$eQTL "+atos(ldprbid.size())+'\n';
            if(fputs_checked(outstr.c_str(),plotfile))
            {
                printf("ERROR: in writing file %s .\n", plot_path.c_str());
                exit(EXIT_FAILURE);
            }
            for(int i=0;i<stend.size()-1;i++)
            {
                outstr=esdata_._epi_prbID[out_epi_id[stend[i]]] +" "+atos(stend[i+1]-stend[i])+'\n';
                if(fputs_checked(outstr.c_str(),plotfile))
                {
                    printf("ERROR: in writing file %s .\n", plot_path.c_str());
                    exit(EXIT_FAILURE);
                }
                for(int j=stend[i];j<stend[i+1];j++)
                {
                    outstr=out_esi_rs[j]+' '+atos(out_beta[j])+' '+atos(out_se[j])+' '+(out_esi_ld[j]+9<1e-6?"NA":atos(out_esi_ld[j]))+'\n';
                    if(fputs_checked(outstr.c_str(),plotfile))
                    {
                        printf("ERROR: in writing file %s .\n", plot_path.c_str());
                        exit(EXIT_FAILURE);
                    }
                }
            }
            fclose(plotfile);
            free_gwas_data( &gdata);
            cout<<"Information for plot has been saved in "<<plot_path<<"."<<endl;
            
        }
    }
    
    void read_geneAnno(string gAnno_file, vector<string> &gene_name, vector<int> &gene_chr, vector<int> &gene_bp1, vector<int> &gene_bp2) {
        ifstream in_gAnno(gAnno_file.c_str());
        if (!in_gAnno) throw ("Error: can not open the file [" + gAnno_file + "] to read.");
        cout << "Reading physical positions of the genes from [" + gAnno_file + "]." << endl;
        string str_buf;
        vector<string> vs_buf;
        while (getline(in_gAnno, str_buf)) {
            if (StrFunc::split_string(str_buf, vs_buf) != 4) throw ("Error: in line \"" + str_buf + "\".");
            gene_chr.push_back(atoi(vs_buf[0].c_str()));
            gene_bp1.push_back(atoi(vs_buf[1].c_str()));
            gene_bp2.push_back(atoi(vs_buf[2].c_str()));
            gene_name.push_back(vs_buf[3]);
        }
        in_gAnno.close();
        cout << "Physical positions of " << gene_name.size() << " genes have been include." << endl;
    }
    
    void sbat_read_snpset(bInfo* bdata, char* snpset_file, vector<string> &set_name,  vector<int> &gene_chr, vector<int> &gene_bp1, vector<int> &gene_bp2, vector< vector<string> > &snpset)
    {
        ifstream in_snpset(snpset_file);
        if (!in_snpset) throw ("Error: can not open the file [" + string(snpset_file) + "] to read.");
        cout << "\nReading SNP sets from [" + string(snpset_file) + "]." << endl;
        string str_buf;
        vector<string> vs_buf, snpset_buf, snp_name;
        map<string, int>::iterator iter;
        int i = 0;
        map<int, int> chr_map;
        
        while (in_snpset>>str_buf) {
            if(str_buf!="END" && str_buf!="end") vs_buf.push_back(str_buf);
            else{
                if(vs_buf.empty()) continue;
                int bp1=INT_MAX;
                int bp2=0;
                int chr=0;
                chr_map.clear();
                for(i = 1; i < vs_buf.size(); i++){
                    iter = bdata->_snp_name_map.find(vs_buf[i]);
                    if (iter != bdata->_snp_name_map.end()){
                        int snpbp=bdata->_bp[iter->second];
                        bp1=bp1<snpbp?bp1:snpbp;
                        bp2=bp2>snpbp?bp2:snpbp;
                        chr=snpbp=bdata->_chr[iter->second];
                        chr_map.insert(pair<int,int>(chr,i));
                        if(chr_map.size()>1)
                        {
                            printf("ERROR: SNPs from different chromosomes found in SNP set: %s.\n", vs_buf[0].c_str());
                            exit(EXIT_FAILURE);
                        }
                        snpset_buf.push_back(vs_buf[i]);
                        snp_name.push_back(vs_buf[i]);
                    }
                }
                
                if(snpset_buf.size()>0)
                {
                    if(bp2-bp1>1000000)
                        printf("WARNING: The size of SNP set: %s is over 1Mb.\n", vs_buf[0].c_str());
                    set_name.push_back(vs_buf[0]);
                    gene_chr.push_back(chr);
                    gene_bp1.push_back(bp1);
                    gene_bp2.push_back(bp2);
                    snpset.push_back(snpset_buf);
                }
                vs_buf.clear();
                snpset_buf.clear();
            }
        }
        in_snpset.close();
        snp_name.erase(unique(snp_name.begin(), snp_name.end()), snp_name.end());        
        cout << snp_name.size() << " SNPs in " << snpset.size() << " sets have been matched and included." << endl;
    }

    void init_smr_wk(SMRWK* smrwk)
    {
        smrwk->bxz.clear(),smrwk->sexz.clear(),smrwk->curId.clear(),smrwk->eName.clear(),smrwk->snpchrom.clear(),smrwk->byz.clear();
        smrwk->seyz.clear(),smrwk->pyz.clear(),smrwk->bpsnp.clear();
    }
    long fill_smr_wk(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, bool refFlg, const char* refSNP, int lowerbp,int upperbp)
    {
        int i=smrwk->cur_prbid;
        long maxid=-9;
        if(esdata->_rowid.empty())
        {
            for (int j = 0; j<bdata->_include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (abs(esdata->_bxz[i][j] + 9) > 1e-6)
                {
                    int snpbp=esdata->_esi_bp[j];
                    int snpchr=esdata->_esi_chr[j];
                    if(snpchr==esdata->_epi_chr[i] && snpchr==smrwk->cur_chr && snpbp>=lowerbp && snpbp<=upperbp)
                    {
                        smrwk->bxz.push_back(esdata->_bxz[i][j]);
                        smrwk->sexz.push_back(esdata->_sexz[i][j]);
                        smrwk->byz.push_back(gdata->byz[j]);
                        smrwk->seyz.push_back(gdata->seyz[j]);
                        smrwk->pyz.push_back(gdata->pvalue[j]);
                        smrwk->curId.push_back(j);
                        smrwk->eName.push_back(esdata->_esi_rs[j]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                        if(refFlg && esdata->_esi_rs[j]==string(refSNP)) maxid=(smrwk->eName.size()-1);
                        smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                    }
                }
            }
            
        }
        else{
            uint64_t beta_start=esdata->_cols[i<<1];
            uint64_t se_start=esdata->_cols[1+(i<<1)];
            uint64_t numsnps=se_start-beta_start;
            for(int j=0;j<numsnps;j++)
            {
                int ge_rowid=esdata->_rowid[beta_start+j];
                int snpbp=esdata->_esi_bp[ge_rowid];
                int snpchr=esdata->_esi_chr[ge_rowid];
                
                if(snpchr==esdata->_epi_chr[i] && snpchr==smrwk->cur_chr && snpbp>=lowerbp && snpbp<=upperbp)
                {
                    smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                    smrwk->sexz.push_back(esdata->_val[se_start+j]);
                    smrwk->byz.push_back(gdata->byz[ge_rowid]);
                    smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                    smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                    smrwk->curId.push_back(ge_rowid);
                    smrwk->eName.push_back(esdata->_esi_rs[ge_rowid]);
                    smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                    if(refFlg && esdata->_esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->eName.size()-1);
                    smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                    
                }
            }
        }
        return maxid;
    }
    long fill_smr_wk(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, bool refFlg, const char* refSNP,int cis_itvl)
    {
        int i=smrwk->cur_prbid;
        long maxid =-9;
        if(esdata->_rowid.empty())
        {
            for (int j = 0; j<bdata->_include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
            {
                if (abs(esdata->_bxz[i][j] + 9) > 1e-6)
                {
                    int snpbp=esdata->_esi_bp[j];
                    int snpchr=esdata->_esi_chr[j];
                    if(snpchr==esdata->_epi_chr[i] && ABS(esdata->_epi_bp[i]-snpbp)<=cis_itvl)
                    {
                        smrwk->bxz.push_back(esdata->_bxz[i][j]);
                        smrwk->sexz.push_back(esdata->_sexz[i][j]);
                        smrwk->byz.push_back(gdata->byz[j]);
                        smrwk->seyz.push_back(gdata->seyz[j]);
                        smrwk->pyz.push_back(gdata->pvalue[j]);
                        smrwk->curId.push_back(j);
                        smrwk->eName.push_back(esdata->_esi_rs[j]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                        if(refFlg && esdata->_esi_rs[j]==string(refSNP)) maxid=(smrwk->eName.size()-1);
                        smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                    }
                }
            }
            
        }
        else{
            uint64_t beta_start=esdata->_cols[i<<1];
            uint64_t se_start=esdata->_cols[1+(i<<1)];
            uint64_t numsnps=se_start-beta_start;
            for(int j=0;j<numsnps;j++)
            {
                int ge_rowid=esdata->_rowid[beta_start+j];
                int snpbp=esdata->_esi_bp[ge_rowid];
                int snpchr=esdata->_esi_chr[ge_rowid];
                if(snpchr==esdata->_epi_chr[i] && ABS(esdata->_epi_bp[i]-snpbp)<=cis_itvl)
                {
                    smrwk->bxz.push_back(esdata->_val[beta_start+j]);
                    smrwk->sexz.push_back(esdata->_val[se_start+j]);
                    smrwk->byz.push_back(gdata->byz[ge_rowid]);
                    smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                    smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                    smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                    smrwk->eName.push_back(esdata->_esi_rs[ge_rowid]);
                    smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                    if(refFlg && esdata->_esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->eName.size()-1);
                    smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
					//smrwk->freq.push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                }
            }
        }
 
        return maxid;
    }
    void rm_cor_sbat(MatrixXd &R, double R_cutoff, int m, vector<int> &rm_ID1) {
        //Modified version of rm_cor_indi from grm.cpp
        
        int i = 0, j = 0, i_buf = 0;
        vector<int> rm_ID2;
        
        //float tmpr = 0;
        for (i = 0; i < m; i++) {
            for (j = 0; j < i; j++) {
                if (fabs(R(i,j)) > R_cutoff ) {
                    rm_ID1.push_back(i);
                    rm_ID2.push_back(j);
                }
            }
        }
        
        // count the number of appearance of each "position" in the vector, which involves a few steps
        vector<int> rm_uni_ID(rm_ID1);
        rm_uni_ID.insert(rm_uni_ID.end(), rm_ID2.begin(), rm_ID2.end());
        stable_sort(rm_uni_ID.begin(), rm_uni_ID.end());
        rm_uni_ID.erase(unique(rm_uni_ID.begin(), rm_uni_ID.end()), rm_uni_ID.end());
        map<int, int> rm_uni_ID_count;
        for (i = 0; i < rm_uni_ID.size(); i++) {
            i_buf = count(rm_ID1.begin(), rm_ID1.end(), rm_uni_ID[i]) + count(rm_ID2.begin(), rm_ID2.end(), rm_uni_ID[i]);
            rm_uni_ID_count.insert(pair<int, int>(rm_uni_ID[i], i_buf));
        }
        
        // swapping
        map<int, int>::iterator iter1, iter2;
        for (i = 0; i < rm_ID1.size(); i++) {
            iter1 = rm_uni_ID_count.find(rm_ID1[i]);
            iter2 = rm_uni_ID_count.find(rm_ID2[i]);
            if (iter1->second < iter2->second) {
                i_buf = rm_ID1[i];
                rm_ID1[i] = rm_ID2[i];
                rm_ID2[i] = i_buf;
            }
        }
        stable_sort(rm_ID1.begin(), rm_ID1.end());
        rm_ID1.erase(unique(rm_ID1.begin(), rm_ID1.end()), rm_ID1.end());
    }
    
    int maxabsid(vector<double> &zsxz, vector<int> &ids)
    {
        
        if(ids.size()==0)
        {
            printf("ERROR: no id values found!\n");
            exit(EXIT_FAILURE);
        }
        int id=ids[0];
        double tmpVal, cmpVal=abs(zsxz[ids[0]]);
        for( int i=1;i<ids.size();i++)
        {
            tmpVal=abs(zsxz[ids[i]]);
            if( cmpVal-tmpVal < 1e-6)
            {
                cmpVal=tmpVal;
                id=ids[i];
            }
        }
        return(id);
    }
    

    void rm_cor_sbat(MatrixXd &R, double R_cutoff, int m, vector<int> &rm_ID1,vector<double> &zxz4smr) {
        
        vector<int> remain,select;
        for(int i=0;i<zxz4smr.size();i++) remain.push_back(i);
        while(!remain.empty())
        {
            int maxid=maxabsid(zxz4smr,remain);
            select.push_back(maxid);
            vector<int> new_remain;
            for(int i=0;i<remain.size();i++)
            {
                if(remain[i]!= maxid )
                {
                    if( fabs(R(maxid,remain[i])) > R_cutoff){
                        rm_ID1.push_back(remain[i]);
                    }else {
                        new_remain.push_back(remain[i]);
                    }
                }
            }
            remain.swap(new_remain);
        }
         stable_sort(rm_ID1.begin(), rm_ID1.end());
    }
    
    void sbat_calcu_lambda(MatrixXd &X, VectorXd &eigenval,VectorXd &eigenvalxy, int &snp_count, double sbat_ld_cutoff, vector<int> &sub_indx,vector<double> &zxz4smr, vector<double> &zyz4smr)
    {
        int m = snp_count;
        
   
        vector<int> rm_ID1;
        double R_cutoff = sbat_ld_cutoff;
        int qi = 0; //alternate index
        
        VectorXd sumsq_x(m);
        for (int j = 0; j < m; j++) sumsq_x[j] = X.col(j).dot(X.col(j));
        
        MatrixXd C = X.transpose() * X;
        X.resize(0,0);
#pragma omp parallel for
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                double d_buf = sqrt(sumsq_x[i] * sumsq_x[j]);
                if(d_buf>0.0) C(i,j) /= d_buf;
                else C(i,j) = 0.0;
            }
        }
        
        if (sbat_ld_cutoff < 1) rm_cor_sbat(C, R_cutoff, m, rm_ID1,zxz4smr);
        //Create new index
        for (int i=0 ; i<m ; i++) {
            if (rm_ID1.size() == 0) sub_indx.push_back(i);
            else {
                if (rm_ID1[qi] == i) qi++; //Skip removed snp
                else sub_indx.push_back(i);
            }
        }
        snp_count = (int)sub_indx.size();
        if (sub_indx.size() < C.size()) { //Build new matrix
            MatrixXd D(sub_indx.size(),sub_indx.size());
            for (int i = 0 ; i < sub_indx.size() ; i++) {
                for (int j = 0 ; j < sub_indx.size() ; j++) {
                    D(i,j) = C(sub_indx[i],sub_indx[j]);
                }
            }
            C = D;
        }
        
        SelfAdjointEigenSolver<MatrixXd> saes(C);
        eigenval = saes.eigenvalues().cast<double>();
         /* save out ld*/
        /*
        cout<<C.rows()<<"*"<<C.cols()<<endl;
        FILE* bsefile=NULL;
        string ldfilename=atos(xh)+".ld.txt";
        bsefile = fopen(ldfilename.c_str(), "w");
        if (!(bsefile)) {
            printf("Open error %s\n", ldfilename.c_str());
            exit(1);
        }
        string outstr="";
        for(int j=0;j<C.rows();j++)
        {
            for(int k=0;k<C.cols()-1;k++) outstr+=atos(C(j,k))+'\t';
            outstr+=atos(C(j,C.rows()-1))+'\n';
        }
        
        if(fputs_checked(outstr.c_str(),bsefile))
        {
            printf("ERROR: in writing file %s .\n", ldfilename.c_str());
            exit(EXIT_FAILURE);
        }
        fclose(bsefile);
       
        cout<<eigenval<<endl;
        */
         /* end of saving*/
        
        VectorXd zs_xz(sub_indx.size()),zs_yz(sub_indx.size());
        for (int i = 0 ; i < sub_indx.size() ; i++){
            zs_xz(i)=zxz4smr[sub_indx[i]];
            zs_yz(i)=zyz4smr[sub_indx[i]];
        }
        MatrixXd xtx=zs_xz*zs_xz.transpose();
        MatrixXd yty=zs_yz*zs_yz.transpose();
        MatrixXd zz=yty.array()/xtx.array();
        MatrixXd numerator=C.array()*(xtx+yty).array()-zz.array();
        VectorXd zsq=zs_xz.array()*zs_xz.array()+zs_yz.array()*zs_yz.array();
        MatrixXd denominator=(zsq*zsq.transpose()).array().sqrt();
        C=numerator.array()/denominator.array();
        SelfAdjointEigenSolver<MatrixXd> saesxy(C);
        eigenvalxy = saesxy.eigenvalues().cast<double>();
    }
    
    double heidi_test(bInfo* bdata,vector<double> &slct_zsxz,vector<uint32_t> &slctId, long slct_maxid,double ld_top, double threshold, int m_hetero, vector<double> &slct_byz, vector<double> &slct_seyz,vector<double> &slct_bxz,vector<double> &slct_sexz, long &nsnp )
    {
        VectorXd ld_v;
        MatrixXd _X;
        vector<int> sn_ids;
        VectorXd tmp_zsxz(slct_zsxz.size());
        for(int j=0;j<slct_zsxz.size();j++) tmp_zsxz(j)=slct_zsxz[j];
        
        make_XMat(bdata,slctId, _X);
        ld_calc_o2m(ld_v,slct_maxid,_X);
        if(abs(ld_top-1)<1e-6) get_square_idxes(sn_ids,tmp_zsxz,threshold);
        else get_square_ldpruning_idxes(sn_ids,tmp_zsxz,threshold,ld_v, slct_maxid,ld_top);
        if(sn_ids.size() < m_hetero) return -9;
        
        VectorXd _byz,_seyz, _bxz,_sexz,_zsxz;
        MatrixXd _X_heidi, _LD_heidi;
        _byz.resize(sn_ids.size());
        _seyz.resize(sn_ids.size());
        _bxz.resize(sn_ids.size());
        _sexz.resize(sn_ids.size());
        _zsxz.resize(sn_ids.size());
        _X_heidi.resize(_X.rows(), sn_ids.size());
        
#pragma omp parallel for
        for(int j=0;j<sn_ids.size();j++)
        {
            _byz[j]=slct_byz[sn_ids[j]];
            _seyz[j]=slct_seyz[sn_ids[j]];
            _bxz[j]=slct_bxz[sn_ids[j]];
            _sexz[j]=slct_sexz[sn_ids[j]];
            _zsxz[j]=slct_zsxz[sn_ids[j]];
            _X_heidi.col(j)=_X.col(sn_ids[j]);
        }
        _X.resize(0,0);
        cor_calc(_LD_heidi, _X_heidi);
        
        _X_heidi.resize(0,0);
        
        nsnp = sn_ids.size();
        double pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, _LD_heidi, &nsnp);
        
        return pdev;
    }

    int smr_setbased_test(bInfo* bdata, vector<uint32_t> &slctId, vector<double> &slct_bxz,vector<double> &slct_sexz,vector<double> &slct_byz,vector<double> &slct_seyz, double p_smr, double ld_top,double &set_pval_smr, double &set_pval_gwas,double &set_pval_eqtl, vector<string> &snp4msmr)
    {
        /* step4: Filter out the SNPs with p-smr threshold and ld-pruning */
        vector<uint32_t> Id4smr;
        vector<double> bxz4smr;
        vector<double> sexz4smr;
        vector<double> byz4smr;
        vector<double> seyz4smr;
        vector<double> zxz4smr;
        vector<double> zyz4smr;
        double z_smr=sqrt(qchisq(p_smr,1));
        for(int j=0;j<slctId.size();j++)
        {
            double ztmp=abs(slct_bxz[j]/slct_sexz[j]);
            if(ztmp>=z_smr)
            {
                Id4smr.push_back(slctId[j]);
                bxz4smr.push_back(slct_bxz[j]);
                sexz4smr.push_back(slct_sexz[j]);
                zxz4smr.push_back(slct_bxz[j]/slct_sexz[j]);
                byz4smr.push_back(slct_byz[j]);
                seyz4smr.push_back(slct_seyz[j]);
                zyz4smr.push_back(slct_byz[j]/slct_seyz[j]);
            }
        }
        int snp_count=(int)Id4smr.size();
        if(snp_count==0) return -9;
        ////    printf("%ld SNPs are passed the threshold %6.2e.\n",Id4smr.size(), p_smr);
        /* step5: multiple-SNP SMR test */
        VectorXd eigenval;
        VectorXd eigenvalxy;
        vector<int> sub_indx;
        MatrixXd _X;
        make_XMat(bdata,Id4smr, _X); //_X: one row one individual, one column one SNP
        double sbat_ld_cutoff=sqrt(ld_top);
        sbat_calcu_lambda(_X, eigenval, eigenvalxy, snp_count,  sbat_ld_cutoff, sub_indx, zxz4smr, zyz4smr); //the index of slectId, snp_count can chage here
        ////     printf("%ld SNPs are passed LD-square threshold %6.2f.\n",sub_indx.size(), ld_top);
        vector<double> zsxysq_slct(sub_indx.size());
        double chisq_zy=0;
        double chisq_zx=0;
        for(int j=0;j<sub_indx.size();j++)
        {
            double z1=byz4smr[sub_indx[j]]/seyz4smr[sub_indx[j]];
            double z2=bxz4smr[sub_indx[j]]/sexz4smr[sub_indx[j]];
            chisq_zy += z1*z1;
            chisq_zx += z2*z2;
            zsxysq_slct[j]=z1*z1*z2*z2/(z1*z1+z2*z2);
        }
        double chisq_o = 0;
        for (int j = 0; j < sub_indx.size(); j++)  chisq_o += zsxysq_slct[j];
        for (int j = 0; j < sub_indx.size(); j++) snp4msmr.push_back(bdata->_snp_name[bdata->_include[Id4smr[sub_indx[j]]]]);
        
        /* save out byz,seyz,bxz,sexz*/
        /*
        FILE* bsefile=NULL;
        string bsefilename=atos(xh)+".bse.txt";
        bsefile = fopen(bsefilename.c_str(), "w");
        if (!(bsefile)) {
            printf("Open error %s\n", bsefilename.c_str());
            exit(1);
        }
         string outstr="";
        for(int j=0;j<sub_indx.size();j++)
        {
             outstr+=atos(byz4smr[sub_indx[j]])+'\t'+atos(seyz4smr[sub_indx[j]])+'\t'+atos(bxz4smr[sub_indx[j]])+'\t'+atos(sexz4smr[sub_indx[j]])+'\n';
        }
        
        if(fputs_checked(outstr.c_str(),bsefile))
        {
            printf("ERROR: in writing file %s .\n", bsefilename.c_str());
            exit(EXIT_FAILURE);
        }
        fclose(bsefile);
        xh++;
         */
        /* end of saving*/
   
        if(sub_indx.size() == 1)
        {
            set_pval_smr = pchisq(chisq_o, 1.0);
            set_pval_gwas= pchisq(chisq_zy, 1.0);
            set_pval_eqtl= pchisq(chisq_zx, 1.0);
            
        }
        else {
            set_pval_smr = chisq_o<1e-8?1:pchisqsum(chisq_o, eigenval);
            set_pval_gwas= chisq_zy<1e-8?1:pchisqsum(chisq_zy, eigenval);
            set_pval_eqtl= chisq_zx<1e-8?1:pchisqsum(chisq_zx, eigenval);
        }
        return snp_count;
    }
    void smr_multipleSNP(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero , char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, char* refSNP, bool heidioffFlag, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,char* setlstName, char* geneAnnoFileName, int expanWind)
    {
        
        setNbThreads(thread_num);
        
        bInfo bdata;
        gwasData gdata;
        eqtlInfo esdata;
        double threshold= chi_val(1,p_hetero);
        bool refFlg=false; // specify a SNP as the heidi SNP
        if(bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(refSNP!=NULL) refFlg=true;
        read_gwas_data( &gdata, gwasFileName);
        read_esifile(&esdata, string(eqtlFileName)+".esi");
        esi_man(&esdata, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, false,  cis_itvl, prbname);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        read_bimfile(&bdata, string(bFileName)+".bim");
        if(snplstName != NULL) extract_snp(&bdata, snplstName);
        if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        allele_check(&bdata, &gdata, &esdata);
        
        read_bedfile(&bdata, string(bFileName)+".bed");
        if (bdata._mu.empty()) calcu_mu(&bdata);
        if (maf > 0)
        {
            filter_snp_maf(&bdata, maf);
            update_geIndx(&bdata, &gdata, &esdata);
        }
        update_gwas(&gdata);
        cout<<"Reading eQTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName)+".epi");
        epi_man(&esdata, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        read_besdfile(&esdata, string(eqtlFileName)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included from %s under current condition.\n",eqtlFileName);
            exit(EXIT_FAILURE);
        }
     
        vector<string> set_name;
        vector< vector<string> > snpset;
        vector<int> gene_chr,gene_bp1,gene_bp2;
        if(setlstName!=NULL) sbat_read_snpset(&bdata,setlstName,set_name,gene_chr, gene_bp1,gene_bp2, snpset );
        else if(geneAnnoFileName!=NULL) read_geneAnno(geneAnnoFileName, set_name, gene_chr, gene_bp1, gene_bp2);
        
      
        unsigned int probNum = esdata._probNum;
        
        
        cout<<endl<<"Performing Multiple-SNP SMR and heterogeneity analysis..... "<<endl;
        float progr0=0.0 , progr1;
        progress_print(progr0);

        cis_itvl=cis_itvl*1000;
        if(expanWind!=-9) expanWind=expanWind*1000;
        string setlstfile = string(outFileName)+".snps4msmr.list";
        FILE* setlst=NULL;
        setlst = fopen(setlstfile.c_str(), "w");
        if (!(setlst)) {
            printf("Open error %s\n", setlstfile.c_str());
            exit(1);
        }
        
        string smrfile = string(outFileName)+".msmr";
        FILE* smr=NULL;
        smr = fopen(smrfile.c_str(), "w");
        if (!(smr)) {
            printf("Open error %s\n", smrfile.c_str());
            exit(1);
        }
        
        //string outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tp_GWAS_multi\tb_eQTL\tse_eQTL\tp_eQTL\tp_eQTL_multi\tb_SMR\tse_SMR\tp_SMR\tnsnp_msmr\tp_SMR_multi\tp_HET\tnsnp_heidi\n";
        
         string outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_SMR_multi\tp_HET\tnsnp_HET\n";
        if(fputs_checked(outstr.c_str(),smr))
        {
            printf("ERROR: in writing file %s .\n", smrfile.c_str());
            exit(EXIT_FAILURE);
        }
        long write_count=0;
        
        SMRWK smrwk;
        for(int i=0;i<probNum;i++)
        {
            progr1=1.0*i/probNum;
            if(progr1-progr0-0.05>1e-6 || i+1==probNum)
            {
                if(i+1==probNum) progr1=1.0;
                progress_print(progr1);
                progr0=progr1;
            }
            
            int probebp=esdata._epi_bp[i];
            int probechr=esdata._epi_chr[i];
            string probename=esdata._epi_prbID[i];
            string probegene=esdata._epi_gene[i];
            if(set_name.size()>0)
            {
                //gene list based or set list based
                for(int ii=0;ii<set_name.size();ii++)
                {
                    init_smr_wk(&smrwk);
                    smrwk.cur_prbid=i;
                    smrwk.cur_chr=gene_chr[ii];
                    
                    int lowerbp=gene_bp1[ii];
                    int upperbp=gene_bp2[ii];
                    long maxid =fill_smr_wk(&bdata, &gdata, &esdata, &smrwk, refFlg, refSNP,lowerbp, upperbp);
                    if(refFlg && maxid==-9) continue; //heidi SNP is not in selected SNPs
                    if (smrwk.bxz.size() == 0) continue;
                    
                    vector<uint32_t> slctId;
                    vector<double> slct_bxz, slct_sexz, slct_byz, slct_seyz, slct_zsxz;
                    vector<string> slct_snpName;
                    long slct_maxid;
                    if(snpset.size()==0)
                    {
                        // gene list
                        slctId.swap(smrwk.curId);
                        slct_bxz.swap(smrwk.bxz);
                        slct_sexz.swap(smrwk.sexz);
                        slct_byz.swap(smrwk.byz);
                        slct_seyz.swap(smrwk.seyz);
                        slct_snpName.swap(smrwk.eName);
                    }else {
                        // set list
                        vector<int> matchidx;
                        match_only(snpset[ii], smrwk.eName, matchidx);
                        if(refFlg && find(matchidx.begin(),matchidx.end(),maxid)==matchidx.end()) continue;
                        
                        for(int j=0;j<matchidx.size();j++)
                        {
                            slctId.push_back(smrwk.curId[matchidx[j]]);
                            slct_bxz.push_back(smrwk.bxz[matchidx[j]]);
                            slct_sexz.push_back(smrwk.sexz[matchidx[j]]);
                            slct_byz.push_back(smrwk.byz[matchidx[j]]);
                            slct_seyz.push_back(smrwk.seyz[matchidx[j]]);
                            slct_snpName.push_back(smrwk.eName[matchidx[j]]);
                        }
                    }
                    
                    Map<VectorXd> ei_bxz(&slct_bxz[0],slct_bxz.size());
                    Map<VectorXd> ei_sexz(&slct_sexz[0],slct_sexz.size());
                    VectorXd zsxz;
                    zsxz=ei_bxz.array()/ei_sexz.array();
                    if(!refFlg) maxid=max_abs_id(zsxz);
                    string topsnpname=slct_snpName[maxid];
                    slct_maxid=maxid;
                    for(int j=0;j<slct_bxz.size();j++) slct_zsxz.push_back(zsxz(j));
                    
					int out_raw_id = slctId[slct_maxid];
                    double bxz_max = slct_bxz[slct_maxid];
                    double sexz_max = slct_sexz[slct_maxid];
                    double byz_max = slct_byz[slct_maxid];
                    double seyz_max = slct_seyz[slct_maxid];
                    double bxy_max = byz_max / bxz_max;
                    double sexy_max = sqrt((seyz_max * seyz_max * bxz_max * bxz_max + sexz_max * sexz_max * byz_max * byz_max) / (bxz_max * bxz_max * bxz_max * bxz_max));
                    double chisqxy = bxy_max*bxy_max / (sexy_max*sexy_max);
                    double zsxz_max = bxz_max / sexz_max;
                    double zsyz_max = byz_max / seyz_max;
                    double pxz_max = pchisq(zsxz_max * zsxz_max, 1);
                    double pyz_max = pchisq(zsyz_max * zsyz_max, 1);
                    double pxy_max = pchisq(chisqxy, 1);

                    double set_pval_smr=-9;
                    double set_pval_gwas=-9;
                    double set_pval_eqtl=-9;
                     vector<string> snp4msmr;
                    int snp_count=smr_setbased_test(&bdata, slctId, slct_bxz,slct_sexz,slct_byz,slct_seyz, p_smr, ld_top,set_pval_smr, set_pval_gwas,set_pval_eqtl, snp4msmr);
                    if(snp_count==-9) continue;
                    
                    /* output snp set list of MSMR test*/
                    string setstr=probename+'\n';
                    if(fputs_checked(setstr.c_str(),setlst))
                    {
                        printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    for(int j=0;j<snp4msmr.size();j++)
                    {
                        setstr=snp4msmr[j]+'\n';
                        if(fputs_checked(setstr.c_str(),setlst))
                        {
                            printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                            exit(EXIT_FAILURE);
                        }
                    }
                    setstr="end\n";
                    if(fputs_checked(setstr.c_str(),setlst))
                    {
                        printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    /* end of output */
                    
                    /* step6: HEIDI test */
                    long nsnp=-9;
                    double pdev=-9;
                    if(!heidioffFlag)  pdev= heidi_test(&bdata,slct_zsxz,slctId, slct_maxid, ld_top,  threshold,  m_hetero, slct_byz, slct_seyz,slct_bxz,slct_sexz, nsnp );
                    
                    outstr = probename + '\t' + atos(probechr) + '\t' + probegene + '\t' + atos(probebp) + '\t' + topsnpname + '\t' + atos(esdata._esi_chr[out_raw_id]) + '\t' + atos(esdata._esi_bp[out_raw_id]) + '\t' + esdata._esi_allele1[out_raw_id] + '\t' + esdata._esi_allele2[out_raw_id] + '\t' + atos(bdata._mu[bdata._include[out_raw_id]] / 2) + '\t';
                    //outstr += atos(byz_max) + '\t' + atos(seyz_max) + '\t' + dtos(pyz_max) + '\t' + dtos(set_pval_gwas) + '\t';
                    //outstr += atos(bxz_max) + '\t' + atos(sexz_max) + '\t' + dtos(pxz_max) + '\t' + dtos(set_pval_eqtl) + '\t';
                    //outstr += atos(bxy_max) + '\t' + atos(sexy_max) + '\t' + dtos(pxy_max) + '\t' + atos(snp_count) + '\t' + dtos(set_pval_smr) + '\t' + (pdev > 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp) : "NA") + '\n';
                    outstr += atos(byz_max) + '\t' + atos(seyz_max) + '\t' + dtos(pyz_max) + '\t';
                    outstr += atos(bxz_max) + '\t' + atos(sexz_max) + '\t' + dtos(pxz_max) + '\t';
                    outstr += atos(bxy_max) + '\t' + atos(sexy_max) + '\t' + dtos(pxy_max) + '\t' + dtos(set_pval_smr) + '\t' + (pdev > 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp) : "NA") + '\n';
                    if(fputs_checked(outstr.c_str(),smr))
                    {
                        printf("ERROR: in writing file %s .\n", smrfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    write_count++;

                }
            }else{
                // top-SNP (/ ref-SNP) based
                init_smr_wk(&smrwk);
                smrwk.cur_prbid=i;
                /* step1: get cis-eQTLs */
                long maxid =fill_smr_wk(&bdata, &gdata, &esdata, &smrwk, refFlg, refSNP, cis_itvl);
                if(refFlg && maxid==-9) continue; //ref heidi SNP is not in selected SNPs
                if (smrwk.bxz.size() == 0) continue;
              ////  printf("%ld SNPs are included from the cis-region of probe %s.\n",smrwk.bxz.size(),probename.c_str());
                //now if sepcify reference SNP, maxid point to this SNP, otherwise maxid is -9
                
                /* step2: get top-SNP */
                Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
                Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
                VectorXd zsxz;
                zsxz=ei_bxz.array()/ei_sexz.array();
                if(!refFlg) maxid=max_abs_id(zsxz); // now maxid point to the sig eQTL SNP or ref SNP in the new datastruct(not the raw).
                
                string topsnpname=smrwk.eName[maxid];
                //cout<<maxid<<":"<<topsnpname<<":"<<esdata._esi_rs[smrwk.curId[maxid]]<<":"<<bdata._snp_name[bdata._include[smrwk.curId[maxid]]]<<":"<<gdata.snpName[smrwk.curId[maxid]]<<endl;
                /* step3: extract SNPs around the --set-wind around sig (or ref) SNP */
                vector<uint32_t> slctId;
                vector<double> slct_bxz;
                vector<double> slct_sexz;
                vector<double> slct_byz;
                vector<double> slct_seyz;
                vector<double> slct_zsxz;
                vector<string> slct_snpName;
                long slct_maxid=-9;
                if(expanWind!=-9){
                    for(int j=0;j<zsxz.size();j++)
                    {
                        int maxid_bp=smrwk.bpsnp[maxid];
                        int tmplower=maxid_bp-expanWind>0?maxid_bp-expanWind:0;
                        int tmpupper=maxid_bp+expanWind;
                        if(smrwk.bpsnp[j]>=tmplower && smrwk.bpsnp[j]<=tmpupper)
                        {
                            slctId.push_back(smrwk.curId[j]); // for get X
                            slct_bxz.push_back(smrwk.bxz[j]);
                            slct_sexz.push_back(smrwk.sexz[j]);
                            slct_byz.push_back(smrwk.byz[j]);
                            slct_seyz.push_back(smrwk.seyz[j]);
                            slct_zsxz.push_back(zsxz(j));
                            slct_snpName.push_back(smrwk.eName[j]);
                            if(j==maxid) slct_maxid=slctId.size()-1;
                        }
                    }
                }else {
                    slctId.swap(smrwk.curId);
                    slct_bxz.swap(smrwk.bxz);
                    slct_sexz.swap(smrwk.sexz);
                    slct_byz.swap(smrwk.byz);
                    slct_seyz.swap(smrwk.seyz);
                    slct_snpName.swap(smrwk.eName);
                    slct_maxid=maxid;
                    for(int j=0;j<slct_bxz.size();j++) slct_zsxz.push_back(zsxz(j));
                }
                
				int out_raw_id = slctId[slct_maxid];
				double bxz_max = slct_bxz[slct_maxid];
				double sexz_max = slct_sexz[slct_maxid];
				double byz_max = slct_byz[slct_maxid];
				double seyz_max = slct_seyz[slct_maxid];
				double bxy_max = byz_max / bxz_max;
				double sexy_max = sqrt((seyz_max * seyz_max * bxz_max * bxz_max + sexz_max * sexz_max * byz_max * byz_max) / (bxz_max * bxz_max * bxz_max * bxz_max));
				double chisqxy = bxy_max*bxy_max / (sexy_max*sexy_max);
				double zsxz_max = bxz_max / sexz_max;
				double zsyz_max = byz_max / seyz_max;
				double pxz_max = pchisq(zsxz_max * zsxz_max, 1);
				double pyz_max = pchisq(zsyz_max * zsyz_max, 1);
				double pxy_max = pchisq(chisqxy, 1);
               // cout<<slct_maxid<<":"<<slct_snpName[slct_maxid]<<endl;
             //// printf("%ld SNPs are included from the window %d kb around the top-SNP /ref-SNP %s.\n",slctId.size(),expanWind/1000,smrwk.eName[maxid].c_str());
                
                double set_pval_smr=-9;
                double set_pval_gwas=-9;
                double set_pval_eqtl=-9;
                vector<string> snp4msmr;
                int snp_count=smr_setbased_test(&bdata, slctId, slct_bxz,slct_sexz,slct_byz,slct_seyz, p_smr, ld_top,set_pval_smr, set_pval_gwas,set_pval_eqtl,snp4msmr);
                if(snp_count==-9) continue;
                /* output snp set list*/
                string setstr=probename+'\n';
                if(fputs_checked(setstr.c_str(),setlst))
                {
                    printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                for(int j=0;j<snp4msmr.size();j++)
                {
                    setstr=snp4msmr[j]+'\n';
                    if(fputs_checked(setstr.c_str(),setlst))
                    {
                        printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                }
                setstr="end\n";
                if(fputs_checked(setstr.c_str(),setlst))
                {
                    printf("ERROR: in writing file %s .\n", setlstfile.c_str());
                    exit(EXIT_FAILURE);
                }
                /* end of output */
                
                /* step6: HEIDI test */
                long nsnp=-9;
                double pdev=-9;
                
                if(!heidioffFlag)  pdev= heidi_test(&bdata,slct_zsxz,slctId, slct_maxid, ld_top,  threshold,  m_hetero, slct_byz, slct_seyz,slct_bxz,slct_sexz, nsnp );
                				
				outstr = probename + '\t' + atos(probechr) + '\t' + probegene + '\t' + atos(probebp) + '\t' + topsnpname + '\t' + atos(esdata._esi_chr[out_raw_id]) + '\t' + atos(esdata._esi_bp[out_raw_id]) + '\t' + esdata._esi_allele1[out_raw_id] + '\t' + esdata._esi_allele2[out_raw_id] + '\t' + atos(bdata._mu[bdata._include[out_raw_id]] / 2) + '\t';
				//outstr += atos(byz_max) + '\t' + atos(seyz_max) + '\t' + dtos(pyz_max) + '\t' + dtos(set_pval_gwas) + '\t';
				//outstr += atos(bxz_max) + '\t' + atos(sexz_max) + '\t' + dtos(pxz_max) + '\t' + dtos(set_pval_eqtl) + '\t';
				//outstr += atos(bxy_max) + '\t' + atos(sexy_max) + '\t' + dtos(pxy_max) + '\t' + atos(snp_count) + '\t' + dtos(set_pval_smr) + '\t' + (pdev > 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp) : "NA") + '\n';
                outstr += atos(byz_max) + '\t' + atos(seyz_max) + '\t' + dtos(pyz_max) + '\t';
                outstr += atos(bxz_max) + '\t' + atos(sexz_max) + '\t' + dtos(pxz_max) + '\t';
                outstr += atos(bxy_max) + '\t' + atos(sexy_max) + '\t' + dtos(pxy_max) + '\t' + dtos(set_pval_smr) + '\t' + (pdev > 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp) : "NA") + '\n';
                if(fputs_checked(outstr.c_str(),smr))
                {
                    printf("ERROR: in writing file %s .\n", smrfile.c_str());
                    exit(EXIT_FAILURE);
                }
                write_count++;

            }
        }
        cout<<"Multiple-SNP SMR and heterogeneity analysis finished.\nSMR and heterogeneity analysis results of "<<write_count<<" sets have been saved in the file [" + smrfile + "]."<<endl;
        cout<<"SNP sets included in Multiple SNPs based SMR test have been saved in the file [" + setlstfile + "]."<<endl;
        fclose(smr);
        fclose(setlst);
        free_gwas_data( &gdata);
        
    }
    
  
    
}