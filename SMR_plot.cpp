//
//  SMR_plot.cpp
//  SMR_CPP
//
//  Created by Futao Zhang on 29/06/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "SMR_plot.h"

namespace SMRDATA
{
    void smre2e(char* bFileName,char* eqtlFileName, char* eqtlFileName2, double maf,char* indilstName, char* snplstName,double p_hetero,double ld_top,int m_hetero, char* indilst2remove, char* snplst2exclde, double p_smr, int cis_itvl,int op_wind, const char* oprobe, vector<string> &mprobe, vector<double> &bsmr, vector<double> &sesmr, vector<double> &psmr, vector<double> &pheidi, vector<int> &nsnp)
    {
        setNbThreads(thread_num);
        string logstr;
         eqtlInfo etrait;
        eqtlInfo esdata;
        bInfo bdata;
        double threshold= chi_val(1,p_hetero);
        cis_itvl=cis_itvl*1000;
        
        if(bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
       
        read_esifile(&etrait, string(eqtlFileName)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&etrait, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait, snplst2exclde);
        read_epifile(&etrait, string(eqtlFileName)+".epi");
        if(oprobe != NULL) extract_eqtl_single_probe(&etrait, oprobe);
        
        
        
        read_besdfile(&etrait, string(eqtlFileName)+".besd");
        if(etrait._rowid.empty() && etrait._bxz.empty())
        {
            printf("No data included under current condition.\n");
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
            
      
        
        //the etrait is not updated, so from now on _esi_include should be used always.
        cout<<"Reading eQTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName2)+".epi");
      
        
        read_besdfile(&esdata, string(eqtlFileName2)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included under current condition.\n");
            exit(EXIT_FAILURE);
        }
        
        
        for( int ii=0;ii<etrait._probNum;ii++)
        {
            
            gwasData gdata;
            gdata.allele_1.resize(etrait._esi_include.size());
            gdata.allele_2.resize(etrait._esi_include.size());
            gdata.byz.resize(etrait._esi_include.size());
            gdata.seyz.resize(etrait._esi_include.size());
            gdata.freq.resize(etrait._esi_include.size());
            gdata.pvalue.resize(etrait._esi_include.size());
            gdata.splSize.resize(etrait._esi_include.size());
            
            
            string traitname=etrait._epi_prbID[ii];
            cout<<"\nPerforming analysis of eTrait [ "+traitname+" ]..."<<endl;
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
            cout<<gdata.snpNum<<" common SNPs are included from eTrait [ "+traitname+" ] summary."<<endl;
            
            int outcome_probe_wind=op_wind*1000;
            int traitbp=etrait._epi_bp[ii];
            int lowerbounder=(traitbp-outcome_probe_wind)>0?(traitbp-outcome_probe_wind):0;
            int upperbounder=traitbp+outcome_probe_wind;
            int traitchr=etrait._epi_chr[ii];
            vector<int> icld_tmp;
            for(int j=0;j<esdata._include.size();j++)
            {
                int bptmp=esdata._epi_bp[esdata._include[j]];
                if(esdata._epi_chr[esdata._include[j]]==traitchr && bptmp>=lowerbounder && bptmp<=upperbounder) icld_tmp.push_back(esdata._include[j]);
            }
            
            int outCount = -1;
            unsigned int probNum = icld_tmp.size();
            // unsigned int probNum = esdata._probNum;
            // vectors for output
            vector<long> out_probid(probNum);
            vector<string> bxy(probNum); // origin is double
            vector<string> sexy(probNum);// origin is double
            vector<string> pxy(probNum); // origin is double
            vector<string> bgwas(probNum); // origin is double
            vector<string> beqtl(probNum); // origin is double
            vector<string> segwas(probNum); // origin is double
            vector<string> seeqtl(probNum); // origin is double
            vector<string> pgwas(probNum); // origin is double
            vector<string> peqtl(probNum); // origin is double
            vector<string> rsid(probNum);
            vector<string> rschr(probNum);
            vector<string> rsa1(probNum);
            vector<string> rsa2(probNum);
            vector<string> rsbp(probNum); //origin is unsigned int
            vector<string> rsfreq(probNum); //origin is unsigned double
            vector<string> prb1(probNum); // origin is double
            vector<string> nsnp_test1(probNum);
            vector<string> top_match1(probNum);  // origin is int
            vector<string> ldrsq(probNum); // origin is double
            
            
            cout<<endl<<"Performing SMR and heterogeneity analysis..... "<<endl;
            float progr0=0.0 , progr1;
            progress_print(progr0);
            
            
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
            
            //for plot
            vector<string> plot_paths;
            for(int ti=0;ti<icld_tmp.size();ti++)
            {
                int i=icld_tmp[ti];
                progr1=1.0*ti/probNum;
                if(progr1-progr0-0.05>1e-6 || ti+1==probNum)
                {
                    if(ti+1==probNum) progr1=1.0;
                    progress_print(progr1);
                    progr0=progr1;
                }
                
                // for plot
                vector<uint32_t> plot_snpidx;
                vector<uint32_t> plot_probeidx;
                vector<double> plot_bxz;
                vector<double> plot_sexz;
                int plotdir_id=ti>>10;
                string plotdir="";
                
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
                           
                            allele1.push_back(esdata._esi_allele1[ge_rowid]);
                            allele2.push_back(esdata._esi_allele2[ge_rowid]);
                            bpsnp.push_back(esdata._esi_bp[ge_rowid]);
                            
                                freq.push_back(bdata._mu[bdata._include[ge_rowid]]/2);
                            
                        }
                    }
                }
                if( maxid==-9) continue; //heidi SNP is not in selected SNPs
                if (bxz.size() == 0) continue;
                
                Map<VectorXd> ei_bxz(&bxz[0],bxz.size());
                Map<VectorXd> ei_sexz(&sexz[0],sexz.size());
                
                zsxz=ei_bxz.array()/ei_sexz.array();
                maxid=max_abs_id(zsxz);
                double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
                
                
                if(pxz_val>p_smr) continue;
                else outCount++;
                
                double bxy_val = byz[maxid] / bxz[maxid];
                double sexy_val = sqrt((seyz[maxid] * seyz[maxid] * bxz[maxid] * bxz[maxid] + sexz[maxid] * sexz[maxid] * byz[maxid] * byz[maxid]) / (bxz[maxid] * bxz[maxid] * bxz[maxid] * bxz[maxid]));
                double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
                
                
                double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
                
                double chisqyz = byz[maxid] / seyz[maxid];
                double pyz_val = pchisq(chisqyz*chisqyz, 1);
                
                out_probid[outCount]= i;
                bxy[outCount]=dtosf(bxy_val);
                sexy[outCount]=dtosf(sexy_val);
                pxy[outCount]=dtos(pxy_val);
                bgwas[outCount]=dtosf(byz[maxid]);
                segwas[outCount]=dtosf(seyz[maxid]);
                beqtl[outCount]=dtosf(bxz[maxid]);
                seeqtl[outCount]=dtosf(sexz[maxid]);
                pgwas[outCount]=dtos(pyz_val);
                peqtl[outCount]=dtos(pxz_val);
                rsid[outCount]=eName[maxid];
                rschr[outCount]=atos(snpchrom[maxid]);
                rsbp[outCount]=itos(bpsnp[maxid]);
                 rsfreq[outCount]=dtosf(freq[maxid]);
                
                rsa1[outCount]=allele1[maxid];
                rsa2[outCount]=allele2[maxid];
                
                
                    //extract info from reference
                    make_XMat(&bdata,curId, _X); //_X: one row one individual, one column one SNP
                    //last vesion ref_snpData was used. row of ref_snpData is SNP, column of ref_snpData is individual
                    //cor_calc(_LD, _X);
                    
                    ld_calc_o2m(ld_v,maxid,_X);
                    
                    //test here
                    // for( int kk=0;kk<curId.size();kk++) if(abs(ld_v[kk]-_LD.col(maxid)[kk])>1e-6) cout<<"wrong"<<endl;
                    // test end
                    
                    sn_ids.clear(); //increase order
                    if(abs(ld_top-1)<1e-6) get_square_idxes(sn_ids,zsxz,threshold);
                    //else get_square_ldpruning_idxes(sn_ids,zsxz,threshold,_LD, maxid,ld_top);
                    else get_square_ldpruning_idxes(sn_ids,zsxz,threshold,ld_v, maxid,ld_top);
                    
                    if(sn_ids.size() < m_hetero)
                    {
                        prb1[outCount]= string("NA");
                        nsnp_test1[outCount]= string("NA");
                        top_match1[outCount]= string("NA");
                        ldrsq[outCount]= string("NA");
                        continue;
                    }
                    _byz.resize(sn_ids.size());
                    _seyz.resize(sn_ids.size());
                    _bxz.resize(sn_ids.size());
                    _sexz.resize(sn_ids.size());
                    _zsxz.resize(sn_ids.size());
                    // _LD_heidi.resize(sn_ids.size(),sn_ids.size());
                    _X_heidi.resize(_X.rows(), sn_ids.size());
                    
#pragma omp parallel for
                    for(int j=0;j<sn_ids.size();j++)
                    {
                        _byz[j]=byz[sn_ids[j]];
                        _seyz[j]=seyz[sn_ids[j]];
                        _bxz[j]=bxz[sn_ids[j]];
                        _sexz[j]=sexz[sn_ids[j]];
                        _zsxz[j]=zsxz[sn_ids[j]];
                        //   for(int k=0;k<=j;k++)_LD_heidi(j,k)=_LD_heidi(k,j)=_LD(sn_ids[j],sn_ids[k]);
                        _X_heidi.col(j)=_X.col(sn_ids[j]);
                    }
                    _X.resize(0,0);
                    cor_calc(_LD_heidi, _X_heidi);
                    
              
                    
                
                        _X_heidi.resize(0,0);
                        
                        long nsnp = sn_ids.size();
                        double pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, _LD_heidi, &nsnp);
                        
                        prb1[outCount] = dtos(pdev);
                        if(nsnp>0) nsnp_test1[outCount] = itos(nsnp);
                        else nsnp_test1[outCount] = string("NA");
                        
                        
                        // top GWAS SNP ?= top eQTL
                        int indx1 = max_abs_id(_zsxz);
                        _byz = _byz.array() / _seyz.array();
                        int indx2 = max_abs_id(_byz);
                        
                        if(indx1 == indx2) top_match1[outCount] =itos(1);
                        else top_match1[outCount] = itos(0);
                        double ldrsqVal = _LD_heidi(indx1, indx2) * _LD_heidi(indx1, indx2);
                        ldrsq[outCount]=dtosf(ldrsqVal);
            }
           
                if(outCount>=0)
                {
                    
                    for (int i = 0;i <=outCount; i++) {
                        mprobe.push_back(esdata._epi_prbID[out_probid[i]]);
                        bsmr.push_back(atof(bxy[i].c_str()));
                        sesmr.push_back(atof(sexy[i].c_str()));
                        psmr.push_back(atof(pxy[i].c_str()));
                        pheidi.push_back(atof(prb1[i].c_str()));
                        nsnp.push_back(atoi(nsnp_test1[i].c_str()));
                        
                    }
                }
            
            free_gwas_data( &gdata);
            
        }
    }
  
    void psudoclone(eqtlInfo* eqtlinfo, eqtlInfo* eqtlinfo2)
    {
        eqtlinfo2->_epi_bp=eqtlinfo->_epi_bp;
        eqtlinfo2->_epi_chr=eqtlinfo->_epi_chr;
        eqtlinfo2->_epi_prbID=eqtlinfo->_epi_prbID;
        eqtlinfo2->_epi_gd=eqtlinfo->_epi_gd;
        eqtlinfo2->_epi_gene=eqtlinfo->_epi_gene;
        eqtlinfo2->_epi_orien=eqtlinfo->_epi_orien;
        eqtlinfo2->_include=eqtlinfo->_include;
        eqtlinfo2->_probe_name_map=eqtlinfo->_probe_name_map;
        eqtlinfo2->_esi_chr=eqtlinfo->_esi_chr;
        eqtlinfo2->_esi_rs=eqtlinfo->_esi_rs;
        eqtlinfo2->_esi_gd=eqtlinfo->_esi_gd;
        eqtlinfo2->_esi_bp=eqtlinfo->_esi_bp;
        eqtlinfo2->_esi_allele1=eqtlinfo->_esi_allele1;
        eqtlinfo2->_esi_allele2=eqtlinfo->_esi_allele2;
        eqtlinfo2->_esi_include=eqtlinfo->_esi_include;
        eqtlinfo2->_snp_name_map=eqtlinfo->_snp_name_map;
        eqtlinfo2->_esi_freq=eqtlinfo->_esi_freq;
        eqtlinfo2->_snpNum=eqtlinfo->_snpNum;
        eqtlinfo2->_probNum=eqtlinfo->_probNum;
        eqtlinfo2->_valNum=eqtlinfo->_valNum;
    }
    void psudoclone(gwasData* gdata, gwasData* gdata2)
    {
        gdata2->snpName=gdata->snpName;
        gdata2->snpBp=gdata->snpBp;
        gdata2->allele_1=gdata->allele_1;
        gdata2->allele_2=gdata->allele_2;
        gdata2->freq=gdata->freq;
        gdata2->byz=gdata->byz;
        gdata2->seyz=gdata->seyz;
        gdata2->pvalue=gdata->pvalue;
        gdata2->splSize=gdata->splSize;
        gdata2->_include=gdata->_include;
        gdata2->snpNum=gdata->snpNum;
    }
    void psudoclone(bInfo* bdata, bInfo* bdata2)
    {
        bdata2->_chr=bdata->_chr;
        bdata2->_snp_name=bdata->_snp_name;
        bdata2->_snp_name_map=bdata->_snp_name_map;
        bdata2->_genet_dst=bdata->_genet_dst;
        bdata2->_bp=bdata->_bp;
        bdata2->_allele1=bdata->_allele1;
        bdata2->_allele2=bdata->_allele2;
        bdata2->_ref_A=bdata->_ref_A;
        bdata2->_other_A=bdata->_other_A;
        bdata2->_include=bdata->_include;
        bdata2->_snp_num=bdata->_snp_num;
        bdata2->_maf=bdata->_maf;
        bdata2->_fid=bdata->_fid;
        bdata2->_pid=bdata->_pid;
        bdata2->_id_map=bdata->_id_map;
        bdata2->_fa_id=bdata->_fa_id;
        bdata2->_mo_id=bdata->_mo_id;
        bdata2->_sex=bdata->_sex;
        bdata2->_pheno=bdata->_pheno;
        bdata2->_indi_num=bdata->_indi_num;
        bdata2->_keep=bdata->_keep;
        bdata2->_mu=bdata->_mu;
        bdata2->_autosome_num=bdata->_autosome_num;
        
    }
   
   
    void get_eTrait(gwasData* gdata, eqtlInfo* edata,int eTraitIdx)
    {
        gdata->allele_1.clear();
        gdata->allele_2.clear();
        gdata->byz.clear();
        gdata->seyz.clear();
        gdata->freq.clear();
        gdata->pvalue.clear();
        gdata->splSize.clear();
        gdata->_include.clear();
        gdata->snpName.clear();
        gdata->snpBp.clear();
        if(edata->_rowid.empty())
        {
            for (int j = 0; j<edata->_esi_include.size(); j++)
            {
                if (abs(edata->_bxz[eTraitIdx][edata->_esi_include[j]] + 9) > 1e-6)
                {
                    gdata->byz.push_back(edata->_bxz[eTraitIdx][edata->_esi_include[j]]);
                    gdata->seyz.push_back(edata->_sexz[eTraitIdx][edata->_esi_include[j]]);
                    double z=edata->_bxz[0][edata->_esi_include[j]]/edata->_sexz[eTraitIdx][edata->_esi_include[j]];
                    gdata->pvalue.push_back(pchisq(z*z,1));
                    gdata->snpName.push_back(edata->_esi_rs[edata->_esi_include[j]]);
                    gdata->allele_1.push_back(edata->_esi_allele1[edata->_esi_include[j]]);
                    gdata->allele_2.push_back(edata->_esi_allele2[edata->_esi_include[j]]);
                    gdata->freq.push_back(edata->_esi_freq[edata->_esi_include[j]]);
                    gdata->snpBp.push_back(edata->_esi_bp[edata->_esi_include[j]]);
                    gdata->_include.push_back(edata->_esi_include[j]); // row id selected
                    
                }
            }
        }
        else
        {
            int ii=0;
            uint64_t beta_start=edata->_cols[ii<<1];
            uint64_t se_start=edata->_cols[1+(ii<<1)];
            uint64_t numsnps=se_start-beta_start;
            for(int j=0;j<numsnps;j++)
            {
                int ge_rowid=edata->_rowid[beta_start+j];
                if(binary_search(edata->_esi_include.begin(), edata->_esi_include.end(), ge_rowid))
                {
                    gdata->byz.push_back(edata->_val[beta_start+j]);
                    gdata->seyz.push_back(edata->_val[se_start+j]);
                    double z=edata->_val[beta_start+j]/edata->_val[se_start+j];
                    gdata->pvalue.push_back(pchisq(z*z,1));
                    gdata->snpName.push_back(edata->_esi_rs[ge_rowid]);
                    gdata->allele_1.push_back(edata->_esi_allele1[ge_rowid]);
                    gdata->allele_2.push_back(edata->_esi_allele2[ge_rowid]);
                    gdata->freq.push_back(edata->_esi_freq[ge_rowid]);
                    gdata->snpBp.push_back(edata->_esi_bp[ge_rowid]);
                    gdata->_include.push_back(ge_rowid);
                    
                }
                
            }
        }
        gdata->snpNum=gdata->snpName.size();
    }
    void plot_triple(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName,char* meqtlFileName, double maf,char* indilstName, char* snplstName,double p_hetero,double ld_top,int m_hetero , char* indilst2remove, char* snplst2exclde, double p_smr, char* refSNP, int cis_itvl, char* prbname, int prbWind,bool prbwindFlag, int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag, char* geneAnnoName, double pthres_me2esmr,double threshpsmrest)
    {
        
        double threshpheidiest=0.05;
        setNbThreads(thread_num);
        if(prbname==NULL) throw("Error: please input probe to plot by the flag --probe.");
        if(!prbwindFlag) throw("Error: please input probe window by the flag --probe-wind.");
        if(bFileName == NULL ) throw("Error: please input Plink file by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data  by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data  by the flag --eqtl-summary.");
        if(meqtlFileName==NULL) throw("Error: please input eQTL summary data  by the flag --eqtl-summary.");
        if(geneAnnoName==NULL) throw("Error: please input gene annotation file by the flag --gene-list.");
        
        vector<int> gene_anno_chr;
        vector<string> gene_anno_genename;
        vector<int> gene_anno_start;
        vector<int> gene_anno_end;
        read_gene_anno(geneAnnoName,gene_anno_chr, gene_anno_genename,gene_anno_start,gene_anno_end);
        map<string, int> gene_anno_map;
        map<string, int>::iterator iter;
        for(int i=0;i<gene_anno_genename.size();i++) gene_anno_map.insert(pair<string,int>(gene_anno_genename[i], i));
        
        eqtlInfo edata;
        eqtlInfo mdata;
        bInfo bdata;
        gwasData gdata;
        read_epifile(&edata, string(eqtlFileName)+".epi");
        extract_prob(&edata, prbname, prbWind);
        
        
        long idx=-9;
        for(int i=0;i<edata._include.size();i++)
             if(edata._epi_prbID[edata._include[i]]== prbname) idx=edata._include[i];
        int prbbp=edata._epi_bp[idx];
        int curchr=edata._epi_chr[idx];
        int plotfrombp=(prbbp-prbWind*1000>0)?prbbp-prbWind*1000:0;
        int plottobp=prbbp+prbWind*1000;
        printf("The position of probe %s is %d. the plot region [%d , %d] is set.\n",prbname,prbbp,plotfrombp,plottobp);
        
        int plotfromkb=ceil(plotfrombp/1000.0);
        int plottokb=ceil(plottobp/1000.0);
        read_epifile(&mdata, string(meqtlFileName)+".epi");
        extract_eqtl_prob(&mdata, curchr, plotfromkb, plottokb);
        
        int from_esnpbp=(plotfrombp-cis_itvl*1000>0)?plotfrombp-cis_itvl*1000:0;
        int end_esnpbp=plottobp+cis_itvl*1000;
       
        printf("\nTo conduct SMR test and HEIDI test, the analysis region [%d , %d] is set.\n",from_esnpbp,end_esnpbp);
        int from_esnpkb=ceil(from_esnpbp/1000.0);
        int end_esnpkb=ceil(end_esnpbp/1000.0);
        
        read_esifile(&edata, string(eqtlFileName)+".esi");
        extract_eqtl_snp(&edata, curchr, from_esnpkb, end_esnpkb);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&edata, snplst2exclde);
        
        read_esifile(&mdata, string(meqtlFileName)+".esi");
        extract_eqtl_snp(&mdata, curchr, from_esnpkb, end_esnpkb);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&mdata, snplst2exclde);
        
        read_famfile(&bdata, string(bFileName)+".fam");
        if(indilstName != NULL) keep_indi(&bdata,indilstName);
        if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
        read_bimfile(&bdata, string(bFileName)+".bim");
        extract_region_bp(&bdata, curchr, from_esnpkb, end_esnpkb);
        if(snplstName != NULL) extract_snp(&bdata, snplstName);
        if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
        
        read_gwas_data( &gdata, gwasFileName);
        vector<int> idx1;
        gdata.snpBp.clear();
        gdata.snpBp.resize(gdata.snpName.size());
        gdata._include.clear();
        
        idx1.resize(gdata.snpName.size());
        #pragma omp parallel for private(iter)
        for (int l = 0; l<gdata.snpName.size(); l++){
            iter = bdata._snp_name_map.find(gdata.snpName[l]);
            if (iter != bdata._snp_name_map.end()) idx1[l]=iter->second;
            else {
                idx1[l]=-9;
            }
        }
        for(int i=0;i<gdata.snpName.size();i++)
        {
            if(idx1[i]!=-9)
            {
                gdata._include.push_back(i);
                gdata.snpBp[i]=bdata._bp[idx1[i]];
            }
        }
        printf("%ld SNPs are extracted from SNP BP: %d Kb to SNP BP: %d Kb of GWAS summary dataset.\n\n",gdata._include.size(),from_esnpkb,end_esnpkb);
        update_gwas(&gdata);
        
        // eSMR begins
        printf("\nPerforming eSMR test the and HEIDI test ...\n");
        eqtlInfo edata_clone;
        eqtlInfo mdata_clone;
        bInfo bdata_clone;
        gwasData gdata_clone;
        psudoclone(&edata, &edata_clone);
        psudoclone(&gdata, &gdata_clone);
        psudoclone(&bdata, &bdata_clone);
        
        allele_check(&bdata_clone, &gdata_clone, &edata_clone);
        read_bedfile(&bdata_clone, string(bFileName)+".bed");
        if (bdata_clone._mu.empty()) calcu_mu(&bdata_clone);
        if (maf > 0)
        {
            filter_snp_maf(&bdata_clone, maf);
            update_geIndx(&bdata_clone, &gdata_clone, &edata_clone);
        }
        
        update_gwas(&gdata_clone);
        read_besdfile(&edata_clone, string(eqtlFileName)+".besd");
        if(edata_clone._rowid.empty() && edata_clone._bxz.empty())
        {
            printf("No data included under current condition.\n");
            exit(EXIT_FAILURE);
        }
        
        vector<SMRRLT> smrrlts;
        smr_heidi_func(smrrlts,  NULL, &bdata_clone,&gdata_clone,&edata_clone,  cis_itvl,  false, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest);
        vector<int> egstart;
        vector<int> egend;
        for(int i=0;i<smrrlts.size();i++)
        {
            iter = gene_anno_map.find(smrrlts[i].Gene);
            if (iter != gene_anno_map.end())
            {
                egstart.push_back(gene_anno_start[iter->second]);
                egend.push_back(gene_anno_end[iter->second]);
            } else {
                egstart.push_back(-9);
                egend.push_back(-9);
            }
        }
        
        // mSMR begins
        printf("\nPerforming mSMR test the and HEIDI test ...\n");
        psudoclone(&gdata, &gdata_clone);
        psudoclone(&mdata, &mdata_clone);
        psudoclone(&bdata, &bdata_clone);
        
        allele_check(&bdata_clone, &gdata_clone, &mdata_clone);
        read_bedfile(&bdata_clone, string(bFileName)+".bed");
        if (bdata_clone._mu.empty()) calcu_mu(&bdata_clone);
        if (maf > 0)
        {
            filter_snp_maf(&bdata_clone, maf);
            update_geIndx(&bdata_clone, &gdata_clone, &mdata_clone);
        }
        
        update_gwas(&gdata_clone);
        read_besdfile(&mdata_clone, string(meqtlFileName)+".besd");
        if(mdata_clone._rowid.empty() && mdata_clone._bxz.empty())
        {
            printf("No data included under current condition.\n");
            exit(EXIT_FAILURE);
        }
        
        vector<SMRRLT> msmrrlts;
        smr_heidi_func(msmrrlts,  NULL, &bdata_clone,&gdata_clone,&mdata_clone,  cis_itvl,  false, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest);
        vector<int> mgstart;
        vector<int> mgend;
        for(int i=0;i<msmrrlts.size();i++)
        {
            iter = gene_anno_map.find(msmrrlts[i].Gene.substr(0,msmrrlts[i].Gene.length()-1));
            if (iter != gene_anno_map.end())
            {
                mgstart.push_back(gene_anno_start[iter->second]);
                mgend.push_back(gene_anno_end[iter->second]);
            } else {
                mgstart.push_back(-9);
                mgend.push_back(-9);
            }
        }

        
        // m2eSMR begins
        printf("\nPerforming m2eSMR test the and HEIDI test ...\n");
        psudoclone(&edata, &edata_clone);
        extract_eqtl_single_probe(&edata_clone, prbname);
        read_besdfile(&edata_clone, string(eqtlFileName)+".besd");
        if(edata_clone._rowid.empty() && edata_clone._bxz.empty())
        {
            printf("No data included under current condition.\n");
            exit(EXIT_FAILURE);
        }
        int eTraitIdx=0;
        get_eTrait(&gdata_clone, &edata_clone,eTraitIdx);
        if(gdata_clone.snpName.size()< m_hetero)
        {
            printf("ERROR: %ld  common SNPs (less than parameter m_hetero: %d are included from eTrait %s.\n",gdata_clone.snpNum,m_hetero,prbname);
            exit(EXIT_FAILURE);
        }
        psudoclone(&mdata, &mdata_clone);
        psudoclone(&bdata, &bdata_clone);
        allele_check(&bdata_clone, &gdata_clone, &mdata_clone);
        read_bedfile(&bdata_clone, string(bFileName)+".bed");
        if (bdata_clone._mu.empty()) calcu_mu(&bdata_clone);
        if (maf > 0)
        {
            filter_snp_maf(&bdata_clone, maf);
            update_geIndx(&bdata_clone, &gdata_clone, &mdata_clone);
        }
        
        update_gwas(&gdata_clone);
        read_besdfile(&mdata_clone, string(meqtlFileName)+".besd");
        if(mdata_clone._rowid.empty() && mdata_clone._bxz.empty())
        {
            printf("No data included under current condition.\n");
            exit(EXIT_FAILURE);
        }
        
        vector<SMRRLT> e2msmrrlts;
        smr_heidi_func(e2msmrrlts,  NULL, &bdata_clone,&gdata_clone,&mdata_clone,  cis_itvl,  false, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest);
        vector<int> e2mgstart;
        vector<int> e2mgend;
        for(int i=0;i<e2msmrrlts.size();i++)
        {
            iter = gene_anno_map.find(e2msmrrlts[i].Gene.substr(0,e2msmrrlts[i].Gene.length()-1));
            if (iter != gene_anno_map.end())
            {
                e2mgstart.push_back(gene_anno_start[iter->second]);
                e2mgend.push_back(gene_anno_end[iter->second]);
            } else {
                e2mgstart.push_back(-9);
                e2mgend.push_back(-9);
            }
        }
        
        //SNP info (combined gwas, eQTL and mQTL)
        vector<string> out_rs;
        vector<int> out_chr;
        vector<int> out_bp;
        vector<string> out_a1;
        vector<string> out_a2;
        map<string,int> snp_name_map;
        long mapsize=0;
        printf("\nRetrieving GWAS summary information from plot region...\n");
        //gwas info
        vector<string> gwas_rs;
        vector<float> gwas_be;
        vector<float> gwas_se;
        for(int i=0;i<gdata.snpNum;i++)
        {
            if(gdata.snpBp[i]>=plotfrombp && gdata.snpBp[i]<=plottobp){
                gwas_rs.push_back(gdata.snpName[i]);
                gwas_be.push_back(gdata.byz[i]);
                gwas_se.push_back(gdata.seyz[i]);
                snp_name_map.insert(pair<string,int>(gdata.snpName[i],mapsize));
                if (mapsize < snp_name_map.size()) {
                    out_rs.push_back(gdata.snpName[i]);
                    out_chr.push_back(curchr);
                    out_bp.push_back(gdata.snpBp[i]);
                    out_a1.push_back(gdata.allele_1[i]);
                    out_a2.push_back(gdata.allele_2[i]);
                    mapsize=snp_name_map.size();
                }
            }
        }
        printf("%ld SNPs are extracted from SNP BP: %d bp to SNP BP %d bp of GWAS summary dataset.\n", gwas_rs.size(),plotfrombp,plottobp);
        //eQTL info
         printf("\nRetrieving eQTL summary information from plot region...\n");
        vector<string> out_epi_name;
        for(int i=0;i<smrrlts.size();i++)
        {
            if(smrrlts[i].p_HET>=threshpheidiest){
                out_epi_name.push_back(smrrlts[i].ProbeID);
            }
        }
        printf("%ld eQTL probes are extracted by the estimated SMR test threshold %e and the estimated HEIDI test threshold %f.\n", out_epi_name.size(),threshpsmrest,threshpheidiest);
        vector<int> include;
        for(int i=0;i<out_epi_name.size();i++)
        {
            iter=edata._probe_name_map.find(out_epi_name[i]);
            if(iter!=edata._probe_name_map.end()) {
                include.push_back(iter->second);
            } else {
                printf("ERROR: Please report this bug even I don't this could happen.\n");
                exit(EXIT_FAILURE);
            }
        }
        stable_sort(include.begin(), include.end());
        edata._include.swap(include);
        extract_eqtl_snp(&edata, curchr, plotfromkb, plottokb);
        read_besdfile(&edata, string(eqtlFileName)+".besd");
        if(edata._rowid.empty() && edata._bxz.empty())
        {
            printf("No data included under current condition.\n");
            exit(EXIT_FAILURE);
        }
        vector<int> out_esi_id;
        vector<int> out_epi_id;
        out_epi_name.clear();
        vector<float> out_beta;
        vector<float> out_se;
        vector<double> out_pval;
        if(edata._valNum==0)
        {
            for(uint32_t i=0;i<edata._probNum;i++)
            {
                for(uint32_t j=0;j<edata._snpNum;j++)
                {
                    double beta=edata._bxz[i][j];
                    double se=edata._sexz[i][j];
                    if(ABS(se+9)<1e-6) continue;
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                    if(pxz<=1)
                    {
                        out_esi_id.push_back(j);
                        out_epi_id.push_back(i);
                        out_epi_name.push_back(edata._epi_prbID[i]);
                        out_beta.push_back(beta);
                        out_se.push_back(se);
                        out_pval.push_back(pxz);
                    }
                }
            }
        }
        else
        {
            if(edata._val.size()==0)
            {
                throw ("Error: No data extracted from the input, please check.\n");
            }
            
            for(uint32_t i=0;i<edata._probNum;i++)
            {
                uint64_t proid=edata._include[i];
                uint64_t pos=edata._cols[proid<<1];
                uint64_t pos1=edata._cols[(proid<<1)+1];
                uint64_t num=pos1-pos;
                for(int j=0;j<num;j++)
                {
                    double beta=edata._val[pos+j];
                    double se=edata._val[pos+j+num];
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                    if(pxz<=1)
                    {
                        out_esi_id.push_back(edata._rowid[pos+j]);
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
            out_esi_rs.push_back(edata._esi_rs[out_esi_id[i]]);
            out_esi_ld.push_back(-9);
            snp_name_map.insert(pair<string,int>(edata._esi_rs[out_esi_id[i]],mapsize));
            if (mapsize < snp_name_map.size()) {
                out_rs.push_back(edata._esi_rs[out_esi_id[i]]);
                out_chr.push_back(edata._esi_chr[out_esi_id[i]]);
                out_bp.push_back(edata._esi_bp[out_esi_id[i]]);
                out_a1.push_back(edata._esi_allele1[out_esi_id[i]]);
                out_a2.push_back(edata._esi_allele2[out_esi_id[i]]);
                mapsize=snp_name_map.size();
            }
        }
        stend.push_back((int)out_esi_id.size());
        
        //mQTL info
        
        printf("\nRetrieving mQTL summary information from plot region...\n");
        vector<string> out_epi_name_m;
        for(int i=0;i<msmrrlts.size();i++)
        {
            if(msmrrlts[i].p_HET>=threshpheidiest){
                out_epi_name_m.push_back(msmrrlts[i].ProbeID);
            }
        }
        printf("%ld mQTL probes are extracted by the estimated SMR test threshold %e and the estimated HEIDI test threshold %f.\n", out_epi_name_m.size(),threshpsmrest,threshpheidiest);
        vector<int> minclude;
        for(int i=0;i<out_epi_name_m.size();i++)
        {
            iter=mdata._probe_name_map.find(out_epi_name_m[i]);
            if(iter!=mdata._probe_name_map.end()) {
                minclude.push_back(iter->second);
            } else {
                printf("ERROR: Please report this bug even I don't this could happen.\n");
                exit(EXIT_FAILURE);
            }
        }
        stable_sort(minclude.begin(), minclude.end());
        mdata._include.swap(minclude);
        extract_eqtl_snp(&mdata, curchr, plotfromkb, plottokb);
        read_besdfile(&mdata, string(meqtlFileName)+".besd");
        if(mdata._rowid.empty() && mdata._bxz.empty())
        {
            printf("No data included under current condition.\n");
            exit(EXIT_FAILURE);
        }
        vector<int> out_esi_id_m;
        vector<int> out_epi_id_m;
        out_epi_name_m.clear();
        vector<float> out_beta_m;
        vector<float> out_se_m;
        vector<double> out_pval_m;
        if(mdata._valNum==0)
        {
            for(uint32_t i=0;i<mdata._probNum;i++)
            {
                for(uint32_t j=0;j<mdata._snpNum;j++)
                {
                    double beta=mdata._bxz[i][j];
                    double se=mdata._sexz[i][j];
                    if(ABS(se+9)<1e-6) continue;
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                    if(pxz<=1)
                    {
                        out_esi_id_m.push_back(j);
                        out_epi_id_m.push_back(i);
                        out_epi_name_m.push_back(mdata._epi_prbID[i]);
                        out_beta_m.push_back(beta);
                        out_se_m.push_back(se);
                        out_pval_m.push_back(pxz);
                    }
                }
            }
        }
        else
        {
            if(mdata._val.size()==0)
            {
                throw ("Error: No data extracted from the input, please check.\n");
            }
            
            for(uint32_t i=0;i<mdata._probNum;i++)
            {
                uint64_t proid=mdata._include[i];
                uint64_t pos=mdata._cols[proid<<1];
                uint64_t pos1=mdata._cols[(proid<<1)+1];
                uint64_t num=pos1-pos;
                for(int j=0;j<num;j++)
                {
                    double beta=mdata._val[pos+j];
                    double se=mdata._val[pos+j+num];
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                    if(pxz<=1)
                    {
                        out_esi_id_m.push_back(mdata._rowid[pos+j]);
                        out_epi_id_m.push_back(i);
                        out_beta_m.push_back(beta);
                        out_se_m.push_back(se);
                        out_pval_m.push_back(pxz);
                    }
                }
            }
        }
        
        vector<double> out_esi_ld_m;
        vector<string> out_esi_rs_m;
        vector<int> stend_m;
        stend_m.push_back(0);
        int curprid_m=out_epi_id_m[0];
        for(int i=0;i<out_esi_id_m.size();i++)
        {
            if(out_epi_id_m[i]!=curprid_m)
            {
                stend_m.push_back(i);
                curprid_m=out_epi_id_m[i];
            }
            out_esi_rs_m.push_back(mdata._esi_rs[out_esi_id_m[i]]);
            out_esi_ld_m.push_back(-9);
            snp_name_map.insert(pair<string,int>(mdata._esi_rs[out_esi_id_m[i]],mapsize));
            if (mapsize < snp_name_map.size()) {
                out_rs.push_back(mdata._esi_rs[out_esi_id_m[i]]);
                out_chr.push_back(mdata._esi_chr[out_esi_id_m[i]]);
                out_bp.push_back(mdata._esi_bp[out_esi_id_m[i]]);
                out_a1.push_back(mdata._esi_allele1[out_esi_id_m[i]]);
                out_a2.push_back(mdata._esi_allele2[out_esi_id_m[i]]);
                mapsize=snp_name_map.size();
            }

        }
        stend_m.push_back((int)out_esi_id_m.size());

        
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
        
       
        
        //for plot
        string plotdir="";
        string plotnm=outFileName;
        for(long j=strlen(outFileName)-1;j>=0;j--)
        {
            if(outFileName[j]=='/')
            {
                plotdir=string(outFileName).substr(0,j+1);
                plotnm=string(outFileName).substr(j+1,strlen(outFileName));
                break;
            }
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
        string plot_path= string(plotdir)+"/"+plotnm+"."+prbname+".txt";
        FILE* plotfile=NULL;
        plotfile = fopen(plot_path.c_str(), "w");
        if (!(plotfile)) {
            printf("Open error %s\n", plot_path.c_str());
            exit(1);
        }
        string outstr="$eprobe "+atos(smrrlts.size())+" "+prbname+'\n';
        if(fputs_checked(outstr.c_str(),plotfile))
        {
            printf("ERROR: in writing file %s .\n", plot_path.c_str());
            exit(EXIT_FAILURE);
        }
      
        for(int i=0;i<smrrlts.size();i++)
        {
            outstr=smrrlts[i].ProbeID+' '+atos(smrrlts[i].ProbeChr)+' '+atos(smrrlts[i].Probe_bp)+' '+smrrlts[i].Gene+' '+(egstart[i]==-9?"NA":atos(egstart[i]))+' '+(egend[i]==-9?"NA":atos(egend[i]))+' '+smrrlts[i].Orien+' '+(smrrlts[i].p_SMR+9<1e-6?"NA":dtos(smrrlts[i].p_SMR))+' '+(smrrlts[i].p_HET+9<1e-6?"NA":dtos(smrrlts[i].p_HET))+'\n';
            if(fputs_checked(outstr.c_str(),plotfile))
            {
                printf("ERROR: in writing file %s .\n", plot_path.c_str());
                exit(EXIT_FAILURE);
            }
        }
        outstr="$mprobe "+atos(msmrrlts.size())+'\n';
        if(fputs_checked(outstr.c_str(),plotfile))
        {
            printf("ERROR: in writing file %s .\n", plot_path.c_str());
            exit(EXIT_FAILURE);
        }
        
        for(int i=0;i<msmrrlts.size();i++)
        {
            outstr=msmrrlts[i].ProbeID+' '+atos(msmrrlts[i].ProbeChr)+' '+atos(msmrrlts[i].Probe_bp)+' '+msmrrlts[i].Gene+' '+(mgstart[i]==-9?"NA":atos(mgstart[i]))+' '+(mgend[i]==-9?"NA":atos(mgend[i]))+' '+msmrrlts[i].Orien+' '+(msmrrlts[i].p_SMR+9<1e-6?"NA":dtos(msmrrlts[i].p_SMR))+' '+(msmrrlts[i].p_HET+9<1e-6?"NA":dtos(msmrrlts[i].p_HET))+'\n';
            if(fputs_checked(outstr.c_str(),plotfile))
            {
                printf("ERROR: in writing file %s .\n", plot_path.c_str());
                exit(EXIT_FAILURE);
            }
        }
        outstr="$m2eprobe "+atos(e2msmrrlts.size())+'\n';
        if(fputs_checked(outstr.c_str(),plotfile))
        {
            printf("ERROR: in writing file %s .\n", plot_path.c_str());
            exit(EXIT_FAILURE);
        }
        
        for(int i=0;i<e2msmrrlts.size();i++)
        {
            outstr=e2msmrrlts[i].ProbeID+' '+atos(e2msmrrlts[i].ProbeChr)+' '+atos(e2msmrrlts[i].Probe_bp)+' '+e2msmrrlts[i].Gene+' '+(e2mgstart[i]==-9?"NA":atos(e2mgstart[i]))+' '+(e2mgend[i]==-9?"NA":atos(e2mgend[i]))+' '+e2msmrrlts[i].Orien+' '+(e2msmrrlts[i].p_SMR+9<1e-6?"NA":dtos(e2msmrrlts[i].p_SMR))+' '+(e2msmrrlts[i].p_HET+9<1e-6?"NA":dtos(e2msmrrlts[i].p_HET))+'\n';
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
        outstr="$eQTL "+atos(edata._probNum)+'\n';
        if(fputs_checked(outstr.c_str(),plotfile))
        {
            printf("ERROR: in writing file %s .\n", plot_path.c_str());
            exit(EXIT_FAILURE);
        }
        for(int i=0;i<stend.size()-1;i++)
        {
            outstr=edata._epi_prbID[out_epi_id[stend[i]]] +" "+atos(stend[i+1]-stend[i])+'\n';
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
        outstr="$mQTL "+atos(mdata._probNum)+'\n';
        if(fputs_checked(outstr.c_str(),plotfile))
        {
            printf("ERROR: in writing file %s .\n", plot_path.c_str());
            exit(EXIT_FAILURE);
        }
        for(int i=0;i<stend_m.size()-1;i++)
        {
            outstr=mdata._epi_prbID[out_epi_id_m[stend_m[i]]] +" "+atos(stend_m[i+1]-stend_m[i])+'\n';
            if(fputs_checked(outstr.c_str(),plotfile))
            {
                printf("ERROR: in writing file %s .\n", plot_path.c_str());
                exit(EXIT_FAILURE);
            }
            for(int j=stend_m[i];j<stend_m[i+1];j++)
            {
                outstr=out_esi_rs_m[j]+' '+atos(out_beta_m[j])+' '+atos(out_se_m[j])+' '+(out_esi_ld_m[j]+9<1e-6?"NA":atos(out_esi_ld_m[j]))+'\n';
                if(fputs_checked(outstr.c_str(),plotfile))
                {
                    printf("ERROR: in writing file %s .\n", plot_path.c_str());
                    exit(EXIT_FAILURE);
                }
            }
        }




        
           fclose(plotfile);
           
            cout<<"Information for plot has been saved in "<<plot_path<<"."<<endl;
        
    }
}