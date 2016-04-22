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
                erlt->_snp_a1.push_back(tmpStr.c_str()[0]);
                iss>>tmpStr;
                to_upper(tmpStr);
                erlt->_snp_a2.push_back(tmpStr.c_str()[0]);
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
        
        
        if(bFlag) read_besdfile(&etrait, string(eqtlFileName)+".besd");
        else      read_esdfile(&etrait, string(eqtlFileName)+".esd");
        
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
        
        if(bFlag) read_besdfile(&esdata, string(eqtlFileName2)+".besd");
        else      read_esdfile(&esdata, string(eqtlFileName2)+".esd");
        
        
        
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
        
        vector<char> rsa1;
        vector<char> rsa2;
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
            gdata.allele_1=(char*)malloc(etrait._esi_include.size()*sizeof(char));
            gdata.allele_2=(char*)malloc(etrait._esi_include.size()*sizeof(char));
            gdata.byz=(double*)malloc(etrait._esi_include.size()*sizeof(double));
            gdata.seyz=(double*)malloc(etrait._esi_include.size()*sizeof(double));
            gdata.freq=(double*)malloc(etrait._esi_include.size()*sizeof(double));
            gdata.pvalue=(double*)malloc(etrait._esi_include.size()*sizeof(double));
            gdata.splSize=(uint32_t*)malloc(etrait._esi_include.size()*sizeof(uint32_t));
            
            
          
           // cout<<"\nPerforming analysis of eTrait [ "+traitname+" ]..."<<endl;
            memset(gdata.byz,0,etrait._esi_include.size()*sizeof(double));
            memset(gdata.seyz,0,etrait._esi_include.size()*sizeof(double));
            memset(gdata.allele_1,0,etrait._esi_include.size()*sizeof(char));
            memset(gdata.allele_2,0,etrait._esi_include.size()*sizeof(char));
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
            
            vector<char> allele1;
            vector<char> allele2;
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
            if(bFlag) read_besdfile(&eqtlinfo, string(eqtlFileName)+".besd");
            else      read_esdfile(&eqtlinfo, string(eqtlFileName)+".esd");
            
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
        if(bFlag) read_besdfile(&esdata, string(eqtlFileName)+".besd");
        else      read_esdfile(&esdata, string(eqtlFileName)+".esd");
        
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
           
            
            if(bFlag) read_besdfile(&eqtlinfo, string(eqtlFileName)+".besd");
            else      read_esdfile(&eqtlinfo, string(eqtlFileName)+".esd");
            
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
        
        smr << "SNP" <<'\t'<< "Chr" <<'\t' << "BP"  << '\t' << "A1" << '\t'<< "A2"<< '\t' << "Probe"<< '\t' << "Probe_Chr"<< '\t'<< "Probe_bp"<< '\t'<<"Gene"<<'\t'<<"b"<<'\t'<< "SE" << '\t'<<"p"<<'\n';
        
        for (int i = 0;i <out_esi_id.size(); i++) {
            smr<<eqtlinfo._esi_rs[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_chr[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_bp[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele1[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele2[out_esi_id[i]]<<'\t'<<eqtlinfo._epi_prbID[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_chr[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_bp[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_gene[out_epi_id[i]]<<'\t'<<out_beta[i]<<'\t'<<out_se[i]<<'\t'<<out_pval[i]<< '\n';
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
            char ga1, ga2, ea1, ea2;
            
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

    void plot(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero , char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, char* refSNP, bool heidioffFlag, int cis_itvl, char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag)
    {
        
        setNbThreads(thread_num);
        if(!prbwindFlag) throw("Error: please input probe window by the flag --probe-wind.");
        if(!heidioffFlag && bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
       
        
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
            esi_man(&esdata, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,  cis_itvl, plotprbname.c_str());
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
           
            if(bFlag) read_besdfile(&esdata, string(eqtlFileName)+".besd");
            else      read_esdfile(&esdata, string(eqtlFileName)+".esd");
            
            
            unsigned int probNum = esdata._probNum;
            
            vector<double> bxz;
            vector<double> sexz;
            vector<uint32_t> curId;
            vector<string> eName;
            vector<int> snpchrom;
            
            vector<char> allele1;
            vector<char> allele2;
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
            
            for(long j=strlen(outFileName)-1;j>=0;j--)
                if(outFileName[j]=='/')
                {
                    plotdir=string(outFileName).substr(0,j+1);
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
            string plot_path= string(plotdir)+"/"+plotprbname+".txt";
            
            //probe info
            vector<string> epi_out_id;
            vector<int> epi_out_chr;
            vector<int> epi_out_bp;
            vector<string> epi_out_gene;
            vector<char> epi_out_orien;
            vector<double> epi_out_pheidi;
            vector<double> epi_out_psmr;
            
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
            //SNP info, the union of gwas and eqtl
            vector<string> out_rs;
            vector<int> out_chr;
            vector<int> out_bp;
            vector<char> out_a1;
            vector<char> out_a2;
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
            vector<char> tmptmpchar;
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
                outstr=epi_out_id[i]+' '+atos(epi_out_chr[i])+' '+atos(epi_out_bp[i])+' '+epi_out_gene[i]+' '+epi_out_orien[i]+' '+(epi_out_psmr[i]+9<1e-6?"NA":dtos(epi_out_psmr[i]))+' '+(epi_out_pheidi[i]+9<1e-6?"NA":dtos(epi_out_pheidi[i]))+'\n';
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
        
       
}