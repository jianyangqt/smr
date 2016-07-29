//
//  SMR_data.cpp
//  SRM_CPP
//
//  Created by Futao Zhang on 29/06/15.
//  Copyright (c) 2015 Futao Zhang. All rights reserved.
//

#include "SMR_data.h"
namespace SMRDATA
{
    int file_read_check(ifstream* in_file, const char* filename)
    {
        in_file->open(filename);
        if(!*in_file)  return 0;
        return 1;
    }
    void progress_print(float progress)
    {
        int barWidth = 70;
        cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i)
        {
            if (i < pos) cout << "=";
            else if (i == pos) cout << ">";
            else cout << " ";
        }
        cout << "] " << int(progress * 100.0) << " %\r";
        cout.flush();    
    }

    void read_famfile(bInfo* bdata, string famfile) {
		bdata->_autosome_num = 22;
        ifstream Fam(famfile.c_str());
        if (!Fam) throw ("Error: can not open the file [" + famfile + "] to read.");
        cout << "Reading PLINK FAM file from [" + famfile + "]." << endl;
        
        int i = 0;
        string str_buf;
        bdata->_fid.clear();
        bdata->_pid.clear();
        bdata->_fa_id.clear();
        bdata->_mo_id.clear();
        bdata->_sex.clear();
        bdata->_pheno.clear();
        while (Fam) {
            Fam >> str_buf;
            if (Fam.eof()) break;
            bdata->_fid.push_back(str_buf);
            Fam >> str_buf;
            bdata->_pid.push_back(str_buf);
            Fam >> str_buf;
            bdata->_fa_id.push_back(str_buf);
            Fam >> str_buf;
            bdata->_mo_id.push_back(str_buf);
            Fam >> str_buf;
            bdata->_sex.push_back(atoi(str_buf.c_str()));
            Fam >> str_buf;
            bdata->_pheno.push_back(atoi(str_buf.c_str()));
        }
        Fam.clear();
        Fam.close();
        bdata->_indi_num = bdata->_fid.size();
        cout << bdata->_indi_num << " individuals to be included from [" + famfile + "]." << endl;
        
        // Initialize _keep
        bdata->_keep.clear();
        bdata->_keep.resize(bdata->_indi_num);
        bdata->_id_map.clear();
        int size = 0;
        for (int i = 0; i < bdata->_indi_num; i++) {
            bdata->_keep[i] = i;
            bdata->_id_map.insert(pair<string, int>(bdata->_fid[i] + ":" + bdata->_pid[i], i));
            if (size == bdata->_id_map.size()) throw ("Error: Duplicate individual ID found: \"" + bdata->_fid[i] + "\t" + bdata->_pid[i] + "\".");
            size = bdata->_id_map.size();
        }
    }
    
    void update_bim(bInfo* bdata,vector<int> &rsnp) {
        int i = 0;
        
        //update bim information
        vector<int> chr_buf, bp_buf;
        vector<string> a1_buf, a2_buf, ref_A_buf, other_A_buf;
        vector<string> snp_name_buf;
        vector<double> genet_dst_buf, impRsq_buf;
        for (i = 0; i < bdata->_snp_num; i++) {
            if (!rsnp[i]) continue;
            chr_buf.push_back(bdata->_chr[i]);
            snp_name_buf.push_back(bdata->_snp_name[i]);
            genet_dst_buf.push_back(bdata->_genet_dst[i]);
            bp_buf.push_back(bdata->_bp[i]);
            a1_buf.push_back(bdata->_allele1[i]);
            a2_buf.push_back(bdata->_allele2[i]);
            ref_A_buf.push_back(bdata->_ref_A[i]);
            other_A_buf.push_back(bdata->_other_A[i]);
            if(bdata->_impRsq.size()>0) impRsq_buf.push_back(bdata->_impRsq[i]);
        }
        bdata->_chr.clear();
        bdata->_snp_name.clear();
        bdata->_genet_dst.clear();
        bdata->_bp.clear();
        bdata->_allele1.clear();
        bdata->_allele2.clear();
        bdata->_ref_A.clear();
        bdata->_other_A.clear();
        bdata->_impRsq.clear();
        bdata->_chr = chr_buf;
        bdata->_snp_name = snp_name_buf;
        bdata->_genet_dst = genet_dst_buf;
        bdata->_bp = bp_buf;
        bdata->_allele1 = a1_buf;
        bdata->_allele2 = a2_buf;
        bdata->_ref_A = ref_A_buf;
        bdata->_other_A = other_A_buf;
        bdata->_impRsq=impRsq_buf;
        bdata->_snp_num = bdata->_chr.size();
        bdata->_include.clear();
        bdata-> _include.resize(bdata->_snp_num);
        bdata->_snp_name_map.clear();
        
        for (i = 0; i < bdata->_snp_num; i++) {
            bdata->_include[i] = i;
            bdata->_snp_name_map.insert(pair<string, int>(bdata->_snp_name[i], i));
        }
    }
    
    void update_fam(bInfo* bdata,vector<int> &rindi) {
        //update fam information
        int i = 0;
        vector<string> fid_buf, pid_buf, fa_id_buf, mo_id_buf;
        vector<int> sex_buf;
        vector<double> pheno_buf;
        for (i = 0; i < bdata->_indi_num; i++) {
            if (!rindi[i]) continue;
            fid_buf.push_back(bdata->_fid[i]);
            pid_buf.push_back(bdata->_pid[i]);
            fa_id_buf.push_back(bdata->_fa_id[i]);
            mo_id_buf.push_back(bdata->_mo_id[i]);
            sex_buf.push_back(bdata->_sex[i]);
            pheno_buf.push_back(bdata->_pheno[i]);
        }
        bdata->_fid.clear();
        bdata->_pid.clear();
        bdata->_fa_id.clear();
        bdata->_mo_id.clear();
        bdata->_sex.clear();
        bdata->_pheno.clear();
        bdata->_fid = fid_buf;
        bdata->_pid = pid_buf;
        bdata->_fa_id = fa_id_buf;
        bdata->_mo_id = mo_id_buf;
        bdata->_sex = sex_buf;
        bdata->_pheno = pheno_buf;
        
        bdata->_indi_num = bdata->_fid.size();
        bdata->_keep.clear();
        bdata->_keep.resize(bdata->_indi_num);
        bdata->_id_map.clear();
        for (i = 0; i < bdata->_indi_num; i++) {
            bdata->_keep[i] = i;
            bdata->_id_map.insert(pair<string, int>(bdata->_fid[i] + ":" + bdata->_pid[i], i));
        }
    }
    
	void update_epi(eqtlInfo* eqtlinfo)
	{
		eqtlinfo->_probNum = eqtlinfo->_include.size();

		vector<int> chr_buf, gd_buf,bp_buf;
		vector<string> prbID_buf, gene_buf;
		vector<char> orien_buf;
		for (int i = 0; i < eqtlinfo->_probNum; i++)
		{
			chr_buf.push_back(eqtlinfo->_epi_chr[eqtlinfo->_include[i]]);
			gd_buf.push_back(eqtlinfo->_epi_gd[eqtlinfo->_include[i]]);
			bp_buf.push_back(eqtlinfo->_epi_bp[eqtlinfo->_include[i]]);
			prbID_buf.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]);
			gene_buf.push_back(eqtlinfo->_epi_gene[eqtlinfo->_include[i]]);
			orien_buf.push_back(eqtlinfo->_epi_orien[eqtlinfo->_include[i]]);
		}
		eqtlinfo->_epi_chr.clear();
		eqtlinfo->_epi_gd.clear();
		eqtlinfo->_epi_bp.clear();
		eqtlinfo->_epi_prbID.clear();
		eqtlinfo->_epi_gene.clear();
		eqtlinfo->_epi_orien.clear();
		eqtlinfo->_epi_chr.swap(chr_buf);
		eqtlinfo->_epi_gd.swap(gd_buf);
		eqtlinfo->_epi_bp.swap(bp_buf);
		eqtlinfo->_epi_prbID.swap(prbID_buf);
		eqtlinfo->_epi_gene.swap(gene_buf);
		eqtlinfo->_epi_orien.swap(orien_buf);
        
        eqtlinfo->_include.clear();
        eqtlinfo->_probe_name_map.clear();
        for (int i = 0; i < eqtlinfo->_probNum; i++)
        {
            eqtlinfo->_include.push_back(i);
            eqtlinfo->_probe_name_map.insert(pair<string,int>(eqtlinfo->_epi_prbID[i],i));
        }
	}

	void update_esi(eqtlInfo* eqtlinfo)
	{
		eqtlinfo->_snpNum = eqtlinfo->_esi_include.size();	

		vector<int> chr_buf, gd_buf, bp_buf;
		vector<string> rs_buf;
		vector<string> allele1_buf, allele2_buf;
        vector<float> freq_buf;
		for (int i = 0; i < eqtlinfo->_snpNum; i++)
		{
			chr_buf.push_back(eqtlinfo->_esi_chr[eqtlinfo->_esi_include[i]]);
			gd_buf.push_back(eqtlinfo->_esi_gd[eqtlinfo->_esi_include[i]]);
			bp_buf.push_back(eqtlinfo->_esi_bp[eqtlinfo->_esi_include[i]]);
			rs_buf.push_back(eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]]);
			allele1_buf.push_back(eqtlinfo->_esi_allele1[eqtlinfo->_esi_include[i]]);
			allele2_buf.push_back(eqtlinfo->_esi_allele2[eqtlinfo->_esi_include[i]]);
            freq_buf.push_back(eqtlinfo->_esi_freq[eqtlinfo->_esi_include[i]]);
		}
		eqtlinfo->_esi_chr.clear();
		eqtlinfo->_esi_gd.clear();
		eqtlinfo->_esi_bp.clear();
		eqtlinfo->_esi_rs.clear();
		eqtlinfo->_esi_allele1.clear();
		eqtlinfo->_esi_allele2.clear();
        eqtlinfo->_esi_freq.clear();
		eqtlinfo->_esi_chr.swap(chr_buf);
		eqtlinfo->_esi_gd.swap(gd_buf);
		eqtlinfo->_esi_bp.swap(bp_buf);
		eqtlinfo->_esi_rs.swap(rs_buf);
		eqtlinfo->_esi_allele1.swap(allele1_buf);
		eqtlinfo->_esi_allele2.swap(allele2_buf);
        eqtlinfo->_esi_freq.swap(freq_buf);
        
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_snp_name_map.clear();
        for (int i = 0; i < eqtlinfo->_snpNum; i++)
        {
            eqtlinfo->_esi_include.push_back(i);
            eqtlinfo->_snp_name_map.insert(pair<string, int>(eqtlinfo->_esi_rs[i], i));
        }
	}

    void read_bimfile(bInfo* bdata,string bimfile) {
        // Read bim file: recombination rate is defined between SNP i and SNP i-1
        int ibuf = 0;
        string cbuf = "0";
        double dbuf = 0.0;
        string str_buf;
        ifstream Bim(bimfile.c_str());
        if (!Bim) throw ("Error: can not open the file [" + bimfile + "] to read.");
        cout << "Reading PLINK BIM file from [" + bimfile + "]." << endl;
        bdata->_chr.clear();
        bdata->_snp_name.clear();
        bdata->_genet_dst.clear();
        bdata->_bp.clear();
        bdata->_allele1.clear();
        bdata->_allele2.clear();
        while (Bim) {
            Bim >> ibuf;
            if (Bim.eof()) break;
            bdata->_chr.push_back(ibuf);
            Bim >> str_buf;
            bdata->_snp_name.push_back(str_buf);
            Bim >> dbuf;
            bdata->_genet_dst.push_back(dbuf);
            Bim >> ibuf;
            bdata->_bp.push_back(ibuf);
            Bim >> cbuf;
            StrFunc::to_upper(cbuf);
            bdata->_allele1.push_back(cbuf.c_str());
            Bim >> cbuf;
            StrFunc::to_upper(cbuf);
            bdata->_allele2.push_back(cbuf.c_str());
        }
        Bim.close();
        bdata->_snp_num = bdata->_chr.size();
        bdata->_ref_A = bdata->_allele1;
        bdata->_other_A = bdata->_allele2;
        cout << bdata->_snp_num << " SNPs to be included from [" + bimfile + "]." << endl;
        
        // Initialize _include
         bdata->_include.clear();
         bdata->_include.resize( bdata->_snp_num);
         bdata->_snp_name_map.clear();
        
        for (int i = 0; i <  bdata->_snp_num; i++) {
            
            bdata->_include[i] = i;
            if( bdata->_snp_name_map.find(bdata->_snp_name[i]) != bdata->_snp_name_map.end()){
                cout << "Warning: Duplicated SNP ID \"" + bdata->_snp_name[i] + "\" ";
                stringstream ss;
                ss << bdata->_snp_name[i] << "_" << i + 1;
                bdata->_snp_name[i] = ss.str();
                cout<<"has been changed to \"" + bdata->_snp_name[i] + "\".\n";
            }
            bdata->_snp_name_map.insert(pair<string, int>(bdata->_snp_name[i], i));
        }
    }
    // some code are adopted from PLINK with modifications
    void read_bedfile(bInfo* bdata, string bedfile)
    {
        int i = 0, j = 0, k = 0;
        
        // Flag for reading individuals and SNPs
        vector<int> rindi, rsnp;
        //get_rindi
        rindi.clear();
        rindi.resize(bdata->_indi_num);
        for (int i = 0; i < bdata->_indi_num; i++) {
            if (bdata->_id_map.find(bdata->_fid[i] + ":" + bdata->_pid[i]) != bdata->_id_map.end()) rindi[i] = 1;
            else rindi[i] = 0;
        }
        //get_rsnp
        rsnp.clear();
        rsnp.resize(bdata->_snp_num);
        for (int i = 0; i < bdata->_snp_num; i++) {
            if (bdata->_snp_name_map.find(bdata->_snp_name[i]) != bdata->_snp_name_map.end()) rsnp[i] = 1;
            else rsnp[i] = 0;
        }
        
        if (bdata->_include.size() == 0) throw ("Error: No SNP is retained for analysis.");
        if (bdata->_keep.size() == 0) throw ("Error: No individual is retained for analysis.");
        
        // Read bed file
        char ch[1];
        bitset<8> b;
        bdata->_snp_1.resize(bdata->_include.size());
        bdata->_snp_2.resize(bdata->_include.size());
        for (i = 0; i < bdata->_include.size(); i++) {
			bdata->_snp_1[i].reserve(bdata->_keep.size());
			bdata->_snp_2[i].reserve(bdata->_keep.size());
        }
        fstream BIT(bedfile.c_str(), ios::in | ios::binary);
        if (!BIT) throw ("Error: can not open the file [" + bedfile + "] to read.");
        cout << "Reading PLINK BED file from [" + bedfile + "] in SNP-major format ..." << endl;
        for (i = 0; i < 3; i++) BIT.read(ch, 1); // skip the first three bytes
        int snp_indx = 0, indi_indx = 0;
        for (j = 0, snp_indx = 0; j < bdata->_snp_num; j++) { // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
            if (!rsnp[j]) {
                for (i = 0; i < bdata->_indi_num; i += 4) BIT.read(ch, 1);
                continue;
            }
            for (i = 0, indi_indx = 0; i < bdata->_indi_num;) {
                BIT.read(ch, 1);
                if (!BIT) throw ("Error: problem with the BED file ... has the FAM/BIM file been changed?");
                b = ch[0];                
                k = 0;
                while (k < 7 && i < bdata->_indi_num) { // change code: 11 for AA; 00 for BB;
                    if (!rindi[i]) k += 2;
                    else {						
                        bdata->_snp_2[snp_indx][indi_indx] = (!b[k++]);
                        bdata->_snp_1[snp_indx][indi_indx] = (!b[k++]);
                        indi_indx++;
                    }
                    i++;
                }
            }
            if (snp_indx == bdata->_include.size()) break;
            snp_indx++;
        }
        BIT.clear();
        BIT.close();
        cout << "Genotype data for " << bdata->_keep.size() << " individuals and " << bdata->_include.size() << " SNPs to be included from [" + bedfile + "]." << endl;
        
        update_fam(bdata, rindi);
        update_bim(bdata, rsnp);
    }
    


    void read_gwas_data(gwasData* gdata, char* gwasFileName)
    {
        ifstream gwasFile;
        if(!file_read_check(&gwasFile, gwasFileName))
        {
            fprintf (stderr, "%s: Couldn't open file %s\n",
                     gwasFileName, strerror (errno));
            exit (EXIT_FAILURE);
        }
        cout << "Reading GWAS summary-level statistics from [" + string(gwasFileName) + "]." << endl;
        gdata->_include.clear();
        gdata->snpName.clear();
        gdata->snpBp.clear();
        gdata->allele_1.clear();
        gdata->allele_2.clear();
        gdata->freq.clear();
        gdata->byz.clear();
        gdata->seyz.clear();
        gdata->pvalue.clear();
        gdata->splSize.clear();
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        gwasFile.getline(buf,MAX_LINE_SIZE);// the header
        vector<string> vs_buf;
        split_string(buf, vs_buf, " \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="SNP") {
            printf("ERROR: %s should have headers that start with \"snp\".\n", gwasFileName);
            exit(EXIT_FAILURE);
        }
        while(!gwasFile.eof())
        {
            gwasFile.getline(buf,MAX_LINE_SIZE);
            
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, " \t\n");
                if(col_num!=8) {
                    printf("ERROR: column number is not right in row %d!\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: \'NA\' of SNP name found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                gdata->snpName.push_back(vs_buf[0]);
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: \'NA\' of allele1 found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                to_upper(vs_buf[1]);
                gdata->allele_1.push_back(vs_buf[1]);
                
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: \'NA\' of allele2 found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                to_upper(vs_buf[2]);
                gdata->allele_2.push_back(vs_buf[2]);
                
                if(vs_buf[3]=="NA" || vs_buf[3]=="na"){
                    //printf("WARNING: \'NA\' of frequency found in row %d.\n", lineNum+2);
                }
                gdata->freq.push_back(atof(vs_buf[3].c_str()));
                
                if(vs_buf[4]=="NA" || vs_buf[4]=="na"){
                    printf("WARNING: \'NA\' of beta found in row %d.\n", lineNum+2);
                    gdata->byz.push_back(0);
                } else {
                    gdata->byz.push_back(atof(vs_buf[4].c_str()));
                }
                if(vs_buf[5]=="NA" || vs_buf[5]=="na"){
                    printf("WARNING: \'NA\' of se found in row %d.\n", lineNum+2);
                    gdata->seyz.push_back(-9);
                } else {
                    gdata->seyz.push_back(atof(vs_buf[5].c_str()));
                }

                gdata->pvalue.push_back(atof(vs_buf[6].c_str()));
                gdata->splSize.push_back(atoi(vs_buf[7].c_str()));
                gdata->_include.push_back(lineNum);
                lineNum++;
            }
        }
        gdata->snpNum=gdata->_include.size();
       cout <<"GWAS summary statistics of "<<gdata->snpNum << " SNPs to be included from [" + string(gwasFileName) + "]." << endl;
        gwasFile.close();
    }
	/*
    void update_include_map(eqtlInfo* eqtlinfo)
    {
        eqtlinfo->_incld_id_map.clear();
        long size=0;
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            eqtlinfo->_incld_id_map.insert(pair<int,int>(eqtlinfo->_esi_include[i],i));
            if (size == eqtlinfo->_incld_id_map.size()) throw ("Error: Duplicated SNP IDs found: \"" + eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]] + "\".");
            size = eqtlinfo->_incld_id_map.size();
        }
        
    }
	*/
    void read_esifile(eqtlInfo* eqtlinfo, string esifile, bool prtscr)
    {
        ifstream esi(esifile.c_str());
        if (!esi) throw ("Error: can not open the file [" + esifile + "] to read.");
        if(prtscr) cout << "Reading eQTL SNP information from [" + esifile + "]." << endl;
        eqtlinfo->_esi_chr.clear();
        eqtlinfo->_esi_rs.clear();
        eqtlinfo->_esi_gd.clear();
        eqtlinfo->_esi_bp.clear();
        eqtlinfo->_esi_allele1.clear();
        eqtlinfo->_esi_allele2.clear();
		eqtlinfo->_esi_include.clear();
        eqtlinfo->_snp_name_map.clear();
        eqtlinfo->_esi_freq.clear();

        char buf[MAX_LINE_SIZE];
        vector<string> vs_buf;
        int lineNum(0);
        while(!esi.eof())
        {
            esi.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, " \t\n");
                if(col_num!=6 && col_num!=7) {
                    printf("ERROR: column number is not right in row %d!\n", lineNum+1);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na" || vs_buf[0]=="0"){
                    printf("ERROR: \"NA\" chromosome found in row %d.\n",lineNum+1);
                    exit(EXIT_FAILURE);
                }
                int tmpchr;
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[0].c_str());
                eqtlinfo->_esi_chr.push_back(tmpchr);
                if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                    printf("ERROR: \'NA\' of SNP name found in row %d.\n", lineNum+1);
                    exit(EXIT_FAILURE);
                }
                if(eqtlinfo->_snp_name_map.find(vs_buf[1]) != eqtlinfo->_snp_name_map.end()){
                    cout << "Warning: Duplicated SNP ID \"" + vs_buf[1] + "\" ";
                    stringstream ss;
                    ss << vs_buf[1] << "_" << lineNum + 1;
                    vs_buf[1] = ss.str();
                    cout<<"has been changed to \"" + vs_buf[1] + "\".\n";
                }
                eqtlinfo->_snp_name_map.insert(pair<string, int>(vs_buf[1], lineNum));
                eqtlinfo->_esi_rs.push_back(vs_buf[1]);
                eqtlinfo->_esi_gd.push_back(atoi(vs_buf[2].c_str()));
                if(vs_buf[3]=="NA" || vs_buf[3]=="na"){
                    printf("ERROR: \'NA\' of SNP BP found in row %d.\n", lineNum+1);
                    exit(EXIT_FAILURE);
                }
                eqtlinfo->_esi_bp.push_back(atoi(vs_buf[3].c_str()));
                if(vs_buf[4]=="NA" || vs_buf[4]=="na") printf("WARNING: \"NA\" allele1 found in row %d.\n", lineNum+1);
                to_upper(vs_buf[4]);
                eqtlinfo->_esi_allele1.push_back(vs_buf[4]);
                if(vs_buf[5]=="NA" || vs_buf[5]=="na") printf("WARNING: \"NA\" allele2 found in row %d.\n", lineNum+1);
                to_upper(vs_buf[5]);
                eqtlinfo->_esi_allele2.push_back(vs_buf[5]);
                if(col_num==7)
                {
                    if(vs_buf[6]=="NA" || vs_buf[6]=="na"){
                        printf("WARNING: \"NA\" frequency found in row %d.\n", lineNum+1);
                        eqtlinfo->_esi_freq.push_back(-9);
                    } else {
                         eqtlinfo->_esi_freq.push_back(atof(vs_buf[6].c_str()));
                    }
                } else {
                    eqtlinfo->_esi_freq.push_back(-9);
                }
                eqtlinfo->_esi_include.push_back(lineNum);
                lineNum++;
            }
        }
        eqtlinfo->_snpNum=lineNum;
        if(prtscr) cout << eqtlinfo->_snpNum << " SNPs to be included from [" + esifile + "]." << endl;
        esi.close();
    }
    void read_epifile(eqtlInfo* eqtlinfo, string epifile, bool prtscr)
    {
        ifstream epi(epifile.c_str());
        if (!epi) throw ("Error: can not open the file [" + epifile + "] to read.");
        if(prtscr) cout << "Reading eQTL probe information from [" + epifile + "]." << endl;
        eqtlinfo->_epi_chr.clear();
        eqtlinfo->_epi_prbID.clear();
        eqtlinfo->_epi_gd.clear();
        eqtlinfo->_epi_bp.clear();
        eqtlinfo->_epi_gene.clear();
        eqtlinfo->_epi_orien.clear();
		eqtlinfo->_include.clear();
        eqtlinfo->_probe_name_map.clear();

        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        while(!epi.eof())
        {
            epi.getline(buf,MAX_LINE_SIZE);
            lineNum++;
        }
        if(buf[0]=='\0') lineNum--;
        eqtlinfo->_probNum=lineNum;
        if(prtscr) cout << eqtlinfo->_probNum << " Probes to be included from [" + epifile + "]." << endl;
        
        eqtlinfo->_epi_chr.resize(lineNum);
        eqtlinfo->_epi_prbID.resize(lineNum);
        eqtlinfo->_epi_gd.resize(lineNum);
        eqtlinfo->_epi_bp.resize(lineNum);
        eqtlinfo->_epi_gene.resize(lineNum);
        eqtlinfo->_epi_orien.resize(lineNum);
        eqtlinfo->_include.resize(lineNum);
        epi.clear(ios::goodbit);
        epi.seekg (0, ios::beg);
        for(int i=0;i<lineNum;i++)
        {
            string tmpStr;
            epi.getline(buf,MAX_LINE_SIZE);
            istringstream iss(buf);
            iss>>tmpStr;
            if(tmpStr=="NA" || tmpStr=="na" || tmpStr=="-9") {
                printf("ERROR: \"NA\" chromosome found in row %d.\n", i+1);
                exit(EXIT_FAILURE);
            }
            int tmpchr;
            if(tmpStr=="X" || tmpStr=="x") tmpchr=23;
            else if(tmpStr=="Y" || tmpStr=="y") tmpchr=24;
            else tmpchr=atoi(tmpStr.c_str());
            eqtlinfo->_epi_chr[i]=tmpchr;
            iss>>tmpStr;
            eqtlinfo->_include[i]=i;
            
            if(eqtlinfo->_probe_name_map.find(tmpStr) != eqtlinfo->_probe_name_map.end()){
                cout << "Warning: Duplicated probe ID \"" + tmpStr + "\" ";
                stringstream ss;
                ss << tmpStr << "_" << i + 1;
                tmpStr = ss.str();
                cout<<"has been changed to \"" + tmpStr + "\".\n";
            }
            eqtlinfo->_probe_name_map.insert(pair<string, int>(tmpStr, i));
            
            
            eqtlinfo->_epi_prbID[i]=tmpStr;
            iss>>tmpStr;
            eqtlinfo->_epi_gd[i]=atoi(tmpStr.c_str());
            iss>>tmpStr;
            if(tmpStr=="NA" || tmpStr=="na" || tmpStr=="0") {
                printf("ERROR: \"NA\" base position found in row %d.\n", i+1);
                exit(EXIT_FAILURE);
            }
            eqtlinfo->_epi_bp[i]=atoi(tmpStr.c_str());
            iss>>tmpStr;
            eqtlinfo->_epi_gene[i]=tmpStr.c_str();
            iss>>tmpStr;
            eqtlinfo->_epi_orien[i]=tmpStr.c_str()[0];
            
        }
        
        epi.close();
    }    
    
    void read_besdfile(eqtlInfo* eqtlinfo, string besdfile, bool prtscr)
    {
        if (eqtlinfo->_include.size() == 0) throw ("Error: No probe is retained for analysis.");
        if (eqtlinfo->_esi_include.size() == 0) throw ("Error: No SNP is retained for analysis.");
        
        eqtlinfo->_cols.clear();
        eqtlinfo->_rowid.clear();
        eqtlinfo->_val.clear();
        eqtlinfo->_valNum = 0;
        eqtlinfo->_bxz.clear();
        eqtlinfo->_sexz.clear();
        
        // the fastest way is using malloc and memcpy
        char SIGN[sizeof(uint64_t)+8];
        ifstream besd(besdfile.c_str(), ios::in|ios::binary);
        if(!besd)
        {
            fprintf (stderr, "%s: Couldn't open file %s\n",
                     besdfile.c_str(), strerror (errno));
            exit (EXIT_FAILURE);
        }
        if(prtscr)  cout << "Reading eQTL summary-level statistics from [" + besdfile + "]." << endl;
        
        besd.read(SIGN, 4);
        float* flag;
        flag=(float *)SIGN;
        if((int)*flag == SPARSE_FILE_TYPE_2){
            // clear datastruct for dense befor read sparse
            eqtlinfo->_bxz.clear();
            eqtlinfo->_sexz.clear();
            
            uint64_t colNum=(eqtlinfo->_probNum<<1)+1;
            uint64_t valNum;
            uint64_t lSize;
            char* buffer;
            besd.seekg(0,besd.end);
            lSize = besd.tellg();
            
            besd.seekg(4); // same as besd.seekg(4, besd.beg);
            besd.read(SIGN, sizeof(uint64_t));
            valNum=*(uint64_t *)SIGN;
            if( lSize - (sizeof(float) + sizeof(uint64_t) + (colNum+valNum)*sizeof(uint32_t) + valNum*sizeof(float)) != 0) {fputs ("wrong element number",stderr); exit (3);}
            
            
            buffer = (char*) malloc (sizeof(char)*(lSize));
            if (buffer == NULL) {fputs ("Memory error",stderr); exit (1);}
            besd.read(buffer,lSize);
            if (besd.gcount()+sizeof(float) + sizeof(uint64_t) != lSize) {fputs ("Reading error",stderr); exit (2);}
            
            
            uint32_t* ptr;
            ptr=(uint32_t *)buffer;
            
            if(eqtlinfo->_include.size()<eqtlinfo->_probNum || eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)
            {
                eqtlinfo->_cols.resize((eqtlinfo->_include.size()<<1)+1);
                eqtlinfo->_cols[0]=*ptr;
                uint32_t* row_ptr;
                row_ptr=ptr+colNum;
                float* val_ptr;
                val_ptr=(float*)(row_ptr+valNum);
                
                map<int, int > _incld_id_map;
                long size = 0;
                for (int i = 0; i<eqtlinfo->_esi_include.size(); i++)
                {
                    _incld_id_map.insert(pair<int, int>(eqtlinfo->_esi_include[i], i));
                    if (size == _incld_id_map.size()) throw ("Error: Duplicated SNP IDs found: \"" + eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]] + "\".");
                    size = _incld_id_map.size();
                }

                for(int i=0;i<eqtlinfo->_include.size();i++)
                {
                    uint32_t pid=eqtlinfo->_include[i];
                    uint32_t pos=*(ptr+(pid<<1));
                    uint32_t pos1=*(ptr+(pid<<1)+1);
                    uint32_t num=pos1-pos;
                    uint32_t real_num=0;
                    for(int j=0;j<num<<1;j++)
                    {
                        uint32_t rid=*(row_ptr+pos+j);
                        
                        map<int, int>::iterator iter;
                        iter=_incld_id_map.find(rid);
                        if(iter!=_incld_id_map.end())
                        {
                            int sid=iter->second;
                        
                       // long sid=find(eqtlinfo->_esi_include.begin(),eqtlinfo->_esi_include.end(),rid)-eqtlinfo->_esi_include.begin(); //slow
                      //  if(sid<eqtlinfo->_esi_include.size())
                      //  {
                            eqtlinfo->_rowid.push_back(sid);
                            eqtlinfo->_val.push_back(*(val_ptr+pos+j));
                            real_num++;
                        }
                       
                    }
                    eqtlinfo->_cols[(i<<1)+1]=(real_num>>1)+eqtlinfo->_cols[i<<1];
                    eqtlinfo->_cols[i+1<<1]=real_num+eqtlinfo->_cols[i<<1];
                }
                eqtlinfo->_valNum = eqtlinfo->_val.size();
                if(prtscr)  cout<<"eQTL summary-level statistics of "<<eqtlinfo->_include.size()<<" Probes and "<<eqtlinfo->_esi_include.size()<<" SNPs to be included from [" + besdfile + "]." <<endl;
                if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
                if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
            }
            else
            {
                eqtlinfo->_cols.resize(colNum);
                eqtlinfo->_rowid.resize(valNum);
                eqtlinfo->_val.resize(valNum);
                
                for(int i=0;i<colNum;i++) eqtlinfo->_cols[i]=*ptr++;
                for(int i=0;i<valNum;i++) eqtlinfo->_rowid[i]=*ptr++;
                float* val_ptr=(float*)ptr;
                for(int i=0;i<valNum;i++) eqtlinfo->_val[i]=*val_ptr++;
                eqtlinfo->_valNum = valNum;
               if(prtscr)  cout<<"eQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
            }
            // terminate
            free (buffer);
        }
        else if((int)*flag == DENSE_FILE_TYPE_1 )
        {
            // clear datastruct for sparse befor read dense
            eqtlinfo->_cols.clear();
            eqtlinfo->_rowid.clear();
            eqtlinfo->_val.clear();
            eqtlinfo->_valNum = 0;
            uint64_t memsize2use=eqtlinfo->_include.size()*eqtlinfo->_esi_include.size()*2*sizeof(float);
            if(memsize2use>0x200000000) printf("WARNING: %llu GB should be allocated for your besd file.\n",memsize2use>>30);
            eqtlinfo->_bxz.resize(eqtlinfo->_include.size());
            eqtlinfo->_sexz.resize(eqtlinfo->_include.size());
            for(unsigned int i=0;i<eqtlinfo->_include.size();i++)
            {
                eqtlinfo->_bxz[i].resize(eqtlinfo->_esi_include.size());
                eqtlinfo->_sexz[i].resize(eqtlinfo->_esi_include.size());
            }
            char* buffer;
            buffer = (char*) malloc (sizeof(char)*eqtlinfo->_snpNum<<3);
            if (buffer == NULL) {fputs ("Memory error",stderr); exit (1);}
            float* ft;
            float* se_ptr;
            if (eqtlinfo->_include.size()<eqtlinfo->_probNum || eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)  //means with the parameter --extract-probe. This also can read all the probes, but currently I don't think it is good for too many I/Os.
            {
                for(int i=0;i<eqtlinfo->_include.size();i++)
                {
                    unsigned long pid=eqtlinfo->_include[i];
                    besd.seekg(((pid<<1)*eqtlinfo->_snpNum+1)<<2);
                    memset(buffer,0,sizeof(char)*eqtlinfo->_snpNum<<3);
                    besd.read(buffer,eqtlinfo->_snpNum<<3);
                    ft=(float *)buffer;
                    for (int j = 0; j<eqtlinfo->_esi_include.size(); j++) eqtlinfo->_bxz[i][j] = *(ft + eqtlinfo->_esi_include[j]);
                    se_ptr = ft + eqtlinfo->_snpNum;
                    for (int j = 0; j<eqtlinfo->_esi_include.size(); j++) eqtlinfo->_sexz[i][j] = *(se_ptr + eqtlinfo->_esi_include[j]);
                }
                if(prtscr)  std::cout << "eQTL summary-level statistics of " << eqtlinfo->_include.size() << " Probes and " << eqtlinfo->_esi_include.size() << " SNPs to be included from [" + besdfile + "]." << endl;
                if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
                if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
            }
            else
            {
                //without --extract-probe, read with less I/O. and need not to update epi.
                //read with static buffer. If dynamic buffer, 2GB per I/O can be more efficient.
                /*
                 unsigned long long count=0;
                 while(!besd.eof())
                 {
                 besd.read(buf,MAX_LINE_NUM);
                 unsigned long Bread=besd.gcount();
                 buf[Bread]='\0';
                 char* ptr=buf;
                 //while(*ptr != '\0') //can not use this, too many 0x00 in buf
                 while(Bread)
                 {
                 unsigned long pid=count/eqtlinfo->_snpNum;
                 unsigned long sid=count++%eqtlinfo->_snpNum;
                 ft=(float *)ptr;
                 if(pid&1) eqtlinfo->_sexz[pid>>1][sid]=*ft;
                 else eqtlinfo->_bxz[pid>>1][sid]=*ft;
                 ptr+=4;
                 Bread-=4;
                 }
                 }
                 cout<<"eQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
                 */
                
                /*
                 
                 // besd.seekg(0,besd.end);
                 // uint64_t lSize = besd.tellg();
                 //  lSize-=4;
                 //  besd.seekg(4); // same as besd.seekg(4, besd.beg);
                 uint64_t alread=0;
                 
                 char* buff;
                 uint64_t buffszie=0x40000000;
                 buff = (char*) malloc (sizeof(char)*buffszie);
                 if (buff == NULL) {fputs ("Memory error",stderr); exit (1);}
                 memset(buff,0,sizeof(char)*buffszie);
                 uint64_t count=0;
                 while(!besd.eof())
                 {
                 besd.read(buff,buffszie);
                 unsigned long Bread=besd.gcount();
                 alread+=Bread;
                 char* ptr=buff;
                 while(Bread)
                 {
                 unsigned long pid=count/eqtlinfo->_snpNum;
                 unsigned long sid=count++%eqtlinfo->_snpNum;
                 ft=(float *)ptr;
                 if(pid&1) eqtlinfo->_sexz[pid>>1][sid]=*ft;
                 else eqtlinfo->_bxz[pid>>1][sid]=*ft;
                 ptr+=4;
                 Bread-=4;
                 }
                 cout<<alread<<":"<<(alread>>30)<<"GB "<<endl;
                 }
                 cout<<"eQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besd
                 file + "]." <<endl;
                 free(buff);
                 */
                
                /*
                char* buff;
                uint64_t buffszie=0x40000000;
                buff = (char*) malloc (sizeof(char)*buffszie);
                if (buff == NULL) {fputs ("Memory error",stderr); exit (1);}
                memset(buff,0,sizeof(char)*buffszie);
                
                uint64_t perbeta=(eqtlinfo->_snpNum<<2);
                uint64_t probonce=sizeof(char)*buffszie/perbeta;  //should be even number
                probonce>>=1;
                probonce<<=1;
                uint64_t readsize=perbeta*probonce;
                uint64_t probcount=0;
                while(!besd.eof())
                {
                    besd.read(buff,readsize);
                    unsigned long Bread=besd.gcount();
                    float* rptr=(float *)buff;
                    while(Bread)
                    {
                        memcpy(&eqtlinfo->_bxz[probcount][0],rptr,perbeta);
                        rptr+=eqtlinfo->_snpNum;
                        memcpy(&eqtlinfo->_sexz[probcount++][0],rptr,perbeta);
                        rptr+=eqtlinfo->_snpNum;
                        Bread-=(perbeta<<1);
                    }
                     cout<<probcount<<" done! "<<endl;
                }
                cout<<"eQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
                free(buff);
                 */
                
                char* buff;
                uint64_t buffszie=0x40000000;
                buff = (char*) malloc (sizeof(char)*buffszie);
                if (buff == NULL) {fputs ("Memory error when reading dense BESD file.",stderr); exit (1);}
                memset(buff,0,sizeof(char)*buffszie);
                
                uint64_t perbeta=(eqtlinfo->_snpNum<<2);
                uint64_t probonce=sizeof(char)*buffszie/perbeta;  //should be even number
                probonce>>=1;
                probonce<<=1;
                uint64_t readsize=perbeta*probonce;
                uint64_t probcount=0;
                while(!besd.eof())
                {
                    besd.read(buff,readsize);
                    uint64_t Bread=besd.gcount();
                    char* rptr=buff;
                    while(Bread)
                    {
                        memcpy(&eqtlinfo->_bxz[probcount][0],rptr,perbeta);
                        rptr+=perbeta;
                        memcpy(&eqtlinfo->_sexz[probcount++][0],rptr,perbeta);
                        rptr+=perbeta;
                        Bread-=(perbeta<<1);
                    }
                    printf("Redinging... %3.0f%%\r", 100.0*probcount/eqtlinfo->_probNum);
                    fflush(stdout);
                }
                if(prtscr) cout<<"\neQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
                free(buff);

            }
            free(buffer);
        }
        else if ((int)*flag == SPARSE_FILE_TYPE_1 )
        {
            // clear datastruct for dense befor read sparse
            eqtlinfo->_bxz.clear();
            eqtlinfo->_sexz.clear();
            
            uint64_t colNum=eqtlinfo->_probNum<<1;
            uint64_t valNum;
            uint64_t lSize;
            char* buffer;
            besd.seekg(0,besd.end);
            lSize = besd.tellg();
            
            besd.seekg(4); // same as besd.seekg(4, besd.beg);
            besd.read(SIGN, 4);
            valNum=(uint32_t)*(float *)SIGN; // int to float then float to int back can lose pricision. hence this clause and bellow are unbelievable
             if(lSize-((3+colNum+(valNum<<1))<<2) != 0) {fputs ("wrong element number",stderr); exit (3);}
            
            valNum=((lSize>>2)-3-colNum)>>1;
            
            buffer = (char*) malloc (sizeof(char)*(lSize-8));
            if (buffer == NULL) {fputs ("Memory error",stderr); exit (1);}
            besd.read(buffer,(lSize-8));
            if (besd.gcount()+8 != lSize) {fputs ("Reading error",stderr); exit (2);}
            float* ptr;
            ptr=(float *)buffer;
            
            if(eqtlinfo->_include.size()<eqtlinfo->_probNum || eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)
            {                
                eqtlinfo->_cols.resize((eqtlinfo->_include.size()<<1)+1);
                eqtlinfo->_cols[0]=(uint32_t)*ptr;
                float* row_ptr;
                row_ptr=ptr+colNum+1;
                float* val_ptr;
                val_ptr=row_ptr+valNum;

				map<int, int > _incld_id_map;
				long size = 0;
				for (int i = 0; i<eqtlinfo->_esi_include.size(); i++)
				{
					_incld_id_map.insert(pair<int, int>(eqtlinfo->_esi_include[i], i));
					if (size == _incld_id_map.size()) throw ("Error: Duplicated SNP IDs found: \"" + eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]] + "\".");
					size = _incld_id_map.size();
				}

                for(int i=0;i<eqtlinfo->_include.size();i++)
                {
                    unsigned long pid=eqtlinfo->_include[i];
                    uint32_t pos=(uint32_t)*(ptr+(pid<<1));
                    uint32_t pos1=(uint32_t)*(ptr+(pid<<1)+1);
                    uint32_t num=pos1-pos;
                    uint32_t real_num=0;
                    for(int j=0;j<num<<1;j++)
                    {
                        uint32_t rid=(uint32_t)*(row_ptr+pos+j);
						 
                        map<int, int>::iterator iter;
                        iter=_incld_id_map.find(rid);
                        if(iter!=_incld_id_map.end())
                        {
                            int sid=iter->second;							
                            eqtlinfo->_rowid.push_back(sid);
                            eqtlinfo->_val.push_back(*(val_ptr+pos+j));
                            real_num++;
                        }
                    }
                    eqtlinfo->_cols[(i<<1)+1]=(real_num>>1)+eqtlinfo->_cols[i<<1];
                    eqtlinfo->_cols[i+1<<1]=real_num+eqtlinfo->_cols[i<<1];
                }
                eqtlinfo->_valNum = eqtlinfo->_val.size();
                if(prtscr)  cout<<"eQTL summary-level statistics of "<<eqtlinfo->_include.size()<<" Probes and "<<eqtlinfo->_esi_include.size()<<" SNPs to be included from [" + besdfile + "]." <<endl;
                if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
                if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
            }
            else
            {
                eqtlinfo->_cols.resize(colNum+1);
                eqtlinfo->_rowid.resize(valNum);
                eqtlinfo->_val.resize(valNum);
                for(int i=0;i<=colNum;i++) eqtlinfo->_cols[i]=(uint32_t)*ptr++;
                for(int i=0;i<valNum;i++) eqtlinfo->_rowid[i]=(uint32_t)*ptr++;
                for(int i=0;i<valNum;i++) eqtlinfo->_val[i]=*ptr++;
                eqtlinfo->_valNum = valNum;
                if(prtscr)  cout<<"eQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
            }
            // terminate
            free (buffer);
        }
        else if ((int)*flag == SPARSE_FILE_TYPE_3)
        {
            // clear datastruct for dense befor read sparse
            eqtlinfo->_bxz.clear();
            eqtlinfo->_sexz.clear();
            char* buffer;
            uint64_t colNum=(eqtlinfo->_probNum<<1)+1;
            uint64_t valNum;
            uint64_t lSize;
           
            besd.seekg(0,besd.end);
            lSize = besd.tellg();
            
            besd.seekg(4); // same as besd.seekg(4, besd.beg);
            besd.read(SIGN, sizeof(uint64_t));
            valNum=*(uint64_t *)SIGN;
            if( lSize - (sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float)) != 0) {fputs ("wrong element number",stderr); exit (3);}
            
            
            if(eqtlinfo->_include.size()<eqtlinfo->_probNum || eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)
            {
                
                uint64_t colsize=colNum*sizeof(uint64_t);
                buffer = (char*) malloc (sizeof(char)*(colsize));
                if (buffer == NULL) {fputs ("Memory error when reading sparse BESD file.",stderr); exit (1);}
                besd.read(buffer,colsize);
                
                uint64_t* ptr;
                ptr=(uint64_t *)buffer;

                eqtlinfo->_cols.resize((eqtlinfo->_include.size()<<1)+1);
                eqtlinfo->_cols[0]=*ptr;
             
                map<int, int > _incld_id_map;
                long size = 0;
                for (int i = 0; i<eqtlinfo->_esi_include.size(); i++)
                {
                    _incld_id_map.insert(pair<int, int>(eqtlinfo->_esi_include[i], i));
                    if (size == _incld_id_map.size()) throw ("Error: Duplicated SNP IDs found: \"" + eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]] + "\".");
                    size = _incld_id_map.size();
                }

                uint64_t rowSTART=sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t);
                uint64_t valSTART=sizeof(float) + sizeof(uint64_t) + colNum*sizeof(uint64_t)+valNum*sizeof(uint32_t);
                for(int i=0;i<eqtlinfo->_include.size();i++)
                {
                    uint32_t pid=eqtlinfo->_include[i];
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
                    besd.seekg(rowSTART+pos*sizeof(uint32_t));
                    besd.read(row_char_ptr, 2*num*sizeof(uint32_t));
                    uint32_t* row_ptr=(uint32_t *)row_char_ptr;                    
                    besd.seekg(valSTART+pos*sizeof(float));
                    besd.read(val_char_ptr, 2*num*sizeof(float));
                    float* val_ptr=(float*)val_char_ptr;
                    for(int j=0;j<num<<1;j++)
                    {
                        uint32_t rid=*(row_ptr+j);
                        
                        map<int, int>::iterator iter;
                        iter=_incld_id_map.find(rid);
                        if(iter!=_incld_id_map.end())
                        {
                            int sid=iter->second;
                            
                            eqtlinfo->_rowid.push_back(sid);
                            eqtlinfo->_val.push_back(*(val_ptr+j));
                            real_num++;
                        }
                        
                    }
                    eqtlinfo->_cols[(i<<1)+1]=(real_num>>1)+eqtlinfo->_cols[i<<1];
                    eqtlinfo->_cols[i+1<<1]=real_num+eqtlinfo->_cols[i<<1];
                    free(row_char_ptr);
                    free(val_char_ptr);
                }
                eqtlinfo->_valNum = eqtlinfo->_val.size();
               
                if(prtscr)  cout<<"eQTL summary-level statistics of "<<eqtlinfo->_include.size()<<" Probes to be included from [" + besdfile + "]." <<endl;
                if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
                if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
            }
            else
            {
                buffer = (char*) malloc (sizeof(char)*(lSize));
                if (buffer == NULL) {fputs ("Memory error",stderr); exit (1);}
                besd.read(buffer,lSize);
                if (besd.gcount()+sizeof(float) + sizeof(uint64_t) != lSize) {fputs ("Reading error",stderr); exit (2);}
                uint64_t* ptr;
                ptr=(uint64_t *)buffer;

                eqtlinfo->_cols.resize(colNum);
                eqtlinfo->_rowid.resize(valNum);
                eqtlinfo->_val.resize(valNum);
                
                for(int i=0;i<colNum;i++) eqtlinfo->_cols[i]=*ptr++;
                uint32_t* ptr4B=(uint32_t *)ptr;
                for(int i=0;i<valNum;i++) eqtlinfo->_rowid[i]=*ptr4B++;
                float* val_ptr=(float*)ptr4B;
                for(int i=0;i<valNum;i++) eqtlinfo->_val[i]=*val_ptr++;
                eqtlinfo->_valNum = valNum;
               if(prtscr)  cout<<"eQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes to be included from [" + besdfile + "]." <<endl;
            }
            // terminate
            free (buffer);
        }
        besd.close();
        /*
        if(eqtlinfo->_rowid.empty() && eqtlinfo->_bxz.empty())
        {
            printf("NO data included from eQTL summary data %s, please check.\n",besdfile.c_str()); exit (EXIT_FAILURE);
        }
         */
    }
    
    void filter_probe_null(eqtlInfo* eqtlinfo)
    {
        vector<string> nullprobes;
        cout<<"\nfiltering out the probes with no value..."<<endl;
        if(eqtlinfo->_valNum==0)
        {
            eqtlinfo->_include.clear();
            for (int i = 0; i < eqtlinfo->_probNum; i++)
            {
                bool NA_flag = true;
                for (int j = 0; j<eqtlinfo->_snpNum; j++)
                {
                    if (abs(eqtlinfo->_sexz[i][j] + 9) > 1e-6)
                    {
                        NA_flag = false;
                        break;
                    }
                }
                if (!NA_flag) eqtlinfo->_include.push_back(i);
                else nullprobes.push_back(eqtlinfo->_epi_prbID[i]);
            }
        }
        else{
            eqtlinfo->_include.clear();
            for (int i = 0; i < eqtlinfo->_probNum; i++)
            {
                if (eqtlinfo->_cols[(i<<1)+1] > eqtlinfo->_cols[i<<1]) eqtlinfo->_include.push_back(i);
                else nullprobes.push_back(eqtlinfo->_epi_prbID[i]);
            }
        }
        eqtlinfo->_probe_name_map.clear();
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            eqtlinfo->_probe_name_map.insert(pair<string,int>(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]],eqtlinfo->_include[i]));
        }
        if(nullprobes.size()>0)
        {
            string fname="chr"+atos(eqtlinfo->_esi_chr[0])+".nullprobes.log";
            FILE* nullprobefile=fopen(fname.c_str(), "w");
            if (!(nullprobefile)) {
                printf("Error: Failed to open null probe log file.\n");
            }
            for(int i=0;i<nullprobes.size();i++)
            {
                string tmpstr=nullprobes[i]+'\n';
                fputs(tmpstr.c_str(),nullprobefile);
            }
            fclose(nullprobefile);
        }
        
        cout<<eqtlinfo->_include.size()<<" probes to be included."<<endl;
    }
    void filter_snp_null(eqtlInfo* eqtlinfo)
    {
        vector<string> nullSNPs;
        cout<<"\nfiltering out the SNPs with no value..."<<endl;
        vector<int> esi_include;
        if(eqtlinfo->_valNum==0)
        {
            for (int i = 0; i < eqtlinfo->_esi_include.size(); i++)
            {
                bool NA_flag = true;
                for (int j = 0; j<eqtlinfo->_include.size(); j++)
                {
                    if (abs(eqtlinfo->_sexz[eqtlinfo->_include[j]][eqtlinfo->_esi_include[i]] + 9) > 1e-6)
                    {
                        NA_flag = false;
                        break;
                    }
                }
                if (!NA_flag)esi_include.push_back(i);
                else nullSNPs.push_back(eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]]);
            }
        }
        else{
            vector<uint32_t> rowid;
            rowid=eqtlinfo->_rowid;
            getUnique(rowid);
            vector<int> idx;
           
            map<uint32_t, int> id_map;
            map<uint32_t, int>::iterator iter;
            for(int i=0; i<rowid.size(); i++) id_map.insert(pair<uint32_t,int>(rowid[i], i));
            for(int i=0; i< eqtlinfo->_esi_include.size(); i++){
                iter=id_map.find(eqtlinfo->_esi_include[i]);
                if(iter==id_map.end()) idx.push_back(-9);
                else idx.push_back(iter->second);
            }
        
            for (int i = 0; i < idx.size(); i++)
            {
                if (idx[i]!=-9) esi_include.push_back(eqtlinfo->_esi_include[i]);
                else nullSNPs.push_back(eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]]);
            }
        }
        eqtlinfo->_esi_include.swap(esi_include);
        eqtlinfo->_snp_name_map.clear();
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            eqtlinfo->_snp_name_map.insert(pair<string,int>(eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]],eqtlinfo->_esi_include[i]));
        }
        if(nullSNPs.size()>0)
        {
            string fname="chr"+atos(eqtlinfo->_esi_chr[0])+".nullSNPs.log";
            FILE* nullsnpfile=fopen(fname.c_str(), "w");
            if (!(nullsnpfile)) {
                printf("Error: Failed to open null probe log file.\n");
                exit(EXIT_FAILURE);
            }
            for(int i=0;i<nullSNPs.size();i++)
            {
                string tmpstr=nullSNPs[i]+'\n';
                fputs(tmpstr.c_str(),nullsnpfile);
            }
            fclose(nullsnpfile);
        }
        
        cout<<eqtlinfo->_esi_include.size()<<" SNPs to be included."<<endl;
    }
    


    bool has_suffix(const std::string &str, const std::string &suffix)
    {
        return str.size() >= suffix.size() &&
        str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }
    bool has_prefix(const std::string &str, const std::string &prefix)
    {
        return str.size() >= prefix.size() &&
        str.compare(0, prefix.size(), prefix) == 0;
    }

    void get_square_idxes(vector<int> &sn_ids,VectorXd &zsxz,double threshold)
    {
        
        for(int i=0;i<zsxz.size();i++)
        {
           if(zsxz[i]*zsxz[i]-threshold>1e-6) sn_ids.push_back(i);
        }
    }
    void get_square_idxes_sort(vector<int> &sn_ids,VectorXd &zsxz,double threshold)
    {
        vector<double> z;
        vector<int> rk;
        for(int i=0;i<zsxz.size();i++)
        {
            if(zsxz[i]*zsxz[i]-threshold>1e-6)
            {
                sn_ids.push_back(i);
                z.push_back(zsxz(i));
            }
        }
        getRank(z,rk);
        
    }
	void get_square_ldpruning_idxes(vector<int> &sn_ids, VectorXd &zsxz, double threshold, MatrixXd &LD, long maxid, double ld_top)
    {
        for(int i=0;i<zsxz.size();i++)
        {
            if(i!= maxid)
            {
                if((zsxz[i]*zsxz[i]-threshold)>1e-6 && (LD(maxid,i)*LD(maxid,i)-ld_top)<1e-6) sn_ids.push_back(i);
            }
            else{
                 if((zsxz[i]*zsxz[i]-threshold)>1e-6) sn_ids.push_back(i);
            }
            
        }
        
    }
    void get_square_ldpruning_idxes(vector<int> &sn_ids, VectorXd &zsxz, double threshold, VectorXd &ld_v, long maxid, double ld_top)
    {
        for(int i=0;i<zsxz.size();i++)
        {
            if(i!= maxid)
            {
                if((zsxz[i]*zsxz[i]-threshold)>1e-6 && (ld_v(i)*ld_v(i)-ld_top)<1e-6) sn_ids.push_back(i);
            }
            else{
                if((zsxz[i]*zsxz[i]-threshold)>1e-6) sn_ids.push_back(i);
            }
        }
    }
   void est_cov_bxy(MatrixXd &covbxy, VectorXd &_zsxz,VectorXd &_bxy,VectorXd &_seyz,VectorXd &_bxz, MatrixXd &_LD_heidi)
    {
        long nsnp =_zsxz.size();
        if(nsnp>1)
        {          
           MatrixXd bxytbxy= _bxy*_bxy.transpose();
           MatrixXd zsxztzsxz= _zsxz*_zsxz.transpose();
           covbxy=_LD_heidi.array()*((_seyz*_seyz.transpose()).array()/(_bxz*_bxz.transpose()).array() + bxytbxy.array()/zsxztzsxz.array()) - bxytbxy.array()/(zsxztzsxz.array()*zsxztzsxz.array());
             
        }
    }

  
   float bxy_hetero3(VectorXd &_byz, VectorXd &_bxz, VectorXd &_seyz, VectorXd &_sexz, VectorXd &_zsxz, MatrixXd &_LD_heidi, long* snp_num)
    {
        VectorXd _bxy;
        VectorXd _sexy;
        VectorXd dev;
        VectorXd tmp3;
        long nsnp=*snp_num;
        int maxid;        
        float pdev=-1.0;
        MatrixXd covbxy(nsnp,nsnp);
        MatrixXd vdev(nsnp-1,nsnp-1);
       
        VectorXd bxz2;
        dev.resize(nsnp-1);
        if(nsnp>1)
        {
            _bxy=_byz.array()/_bxz.array();
            bxz2=_bxz.array()*_bxz.array();
            _sexy=(_seyz.array()*_seyz.array()*bxz2.array()+_sexz.array()*_sexz.array()*_byz.array()*_byz.array())/(bxz2.array()*bxz2.array()).sqrt();
            
            maxid=max_abs_id(_zsxz);
            
            for(int j=0;j<maxid;j++) dev[j]=_bxy[maxid]-_bxy[j];
            for(int j=maxid+1;j<nsnp;j++) dev[j-1]=_bxy[maxid]-_bxy[j];			

            est_cov_bxy(covbxy, _zsxz,_bxy,_seyz,_bxz,_LD_heidi);            
                       
            double tmp1=covbxy(maxid,maxid);           
            tmp3.resize(nsnp-1);
            for(int i=0; i<maxid; i++) tmp3[i]=covbxy(maxid,i);
            for(int i=maxid+1; i<nsnp; i++) tmp3[i-1]=covbxy(maxid,i);
            // vdev as tmp2
			vdev.block(0, 0, maxid, maxid) = covbxy.block(0, 0, maxid, maxid);
			vdev.block(0, maxid, maxid, nsnp - maxid - 1) = covbxy.block(0, maxid + 1, maxid, nsnp - maxid - 1);
			vdev.block(maxid, 0, nsnp - maxid - 1, maxid) = covbxy.block(maxid + 1, 0, nsnp - maxid - 1, maxid);
			vdev.block(maxid, maxid, nsnp - maxid - 1, nsnp - maxid - 1) = covbxy.block(maxid + 1, maxid + 1, nsnp - maxid - 1, nsnp - maxid - 1);            
			
			// get vdev
			VectorXd v1 = VectorXd::Zero(nsnp - 1);
			v1 = v1.array() + 1.0;
		
			vdev = tmp1 + vdev.array() - (v1*tmp3.transpose()).array() - (tmp3*v1.transpose()).array();
			for (int i = 0; i<nsnp - 1; i++)  vdev(i,i) += 1e-8; // in R code  		

            //tmp3 as vardev            
            for(int i=0; i<maxid; i++) tmp3[i]=tmp1+covbxy(i,i)-2*tmp3[i]+ 1e-8;
            for(int i=maxid+1; i<nsnp; i++) tmp3[i-1]=tmp1+covbxy(i,i)-2*tmp3[i-1]+ 1e-8;
			
            //dev as chisq_dev
            
            for(int i=0;i<nsnp-1;i++) dev[i]=dev[i]*dev[i]/tmp3[i];
            
            double sumChisq_dev=0.0;
            for(int i=0;i<nsnp-1;i++)sumChisq_dev+=dev[i];
			
            //using covbxy to store corr_dev
			covbxy.resize(nsnp - 1, nsnp - 1);
			covbxy = vdev.array() / sqrt((vdev.diagonal()*vdev.diagonal().transpose()).array());
           
			
            // using Eigen Library
            
			SelfAdjointEigenSolver<MatrixXd> es(covbxy);			
            VectorXd lambda;
            lambda=es.eigenvalues();
            
            /*
             EigenSolver<MatrixXd> es(A);
             cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
             cout<< es.eigenvalues().transpose()<<endl;
             */
            
            /*
             MatrixXd D = es.pseudoEigenvalueMatrix();
             cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
            */
            
            pdev= pchisqsum(sumChisq_dev,lambda);
            *snp_num=lambda.size();
         
        }else *snp_num=-9;      
        
        return(pdev);
    }
    

    string dtos(double value)
    {
        stringstream ss;
        ss << scientific<< value;
        // ss << fixed << setprecision(400) << __value;
        return(ss.str());
    }
    string dtosf(double value)
    {
        stringstream ss;        
        ss << fixed << value;
        return(ss.str());
    }
    string itos(int value)
    {
        stringstream ss;
        ss << value;      
        return(ss.str());
    }
    string ltos(long value)
    {
        stringstream ss;
        ss << value;
        return(ss.str());
    }
    
    void free_gwas_data(gwasData* gdata)
    {
        gdata->snpName.clear();
        gdata->allele_1.clear();
        gdata->allele_2.clear();
        gdata->freq.clear();
        gdata->byz.clear();
        gdata->seyz.clear();
        gdata->pvalue.clear();
        gdata->splSize.clear();
    }
            
    void mu_func(bInfo* bdata, int j, vector<double> &fac) {
        int i = 0;
		bdata->_dosage_flag = 0;
        double fcount = 0.0, f_buf = 0.0;
        if (bdata->_dosage_flag) {
            for (i = 0; i < bdata->_keep.size(); i++) {
                if (bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]] < 1e5) {
                    bdata->_mu[bdata->_include[j]] += fac[i] * bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]];
                    fcount += fac[i];
                }
            }
        } else {
            for (i = 0; i < bdata->_keep.size(); i++) {
                if (!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]) {
                    f_buf = (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                    if (bdata->_allele2[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) f_buf = 2.0 - f_buf;
                    bdata->_mu[bdata->_include[j]] += fac[i] * f_buf;
                    fcount += fac[i];
                }
            }
        }
        
        if (fcount > 0.0)bdata->_mu[bdata->_include[j]] /= fcount;
    }
    

    void calcu_mu(bInfo* bdata, bool ssq_flag)
    {
        int i = 0, j = 0;
        
        vector<double> auto_fac(bdata->_keep.size()), xfac(bdata->_keep.size()), fac(bdata->_keep.size());
        for (i = 0; i < bdata->_keep.size(); i++)
        {
            auto_fac[i] = 1.0;
            if (bdata->_sex[bdata->_keep[i]] == 1) xfac[i] = 0.5;
            else if (bdata->_sex[bdata->_keep[i]] == 2) xfac[i] = 1.0;
            fac[i] = 0.5;
        }
        
        cout << "Calculating allele frequencies ..." << endl;
        bdata->_mu.clear();
        bdata->_mu.resize(bdata->_snp_num);
        
#pragma omp parallel for
        for (j = 0; j < bdata->_include.size(); j++)
        {
            if (bdata->_chr[bdata->_include[j]]<(bdata->_autosome_num + 1)) mu_func(bdata, j, auto_fac);
            else if (bdata->_chr[bdata->_include[j]] == (bdata->_autosome_num + 1)) mu_func(bdata,j, xfac);
            else mu_func(bdata, j, fac);
        }
    }

    eigenMatrix reg(vector<double> &y, vector<double> &x, vector<double> &rst, bool table=false)
    {
        int N = x.size();
        if (N != y.size() || N < 1) throw ("Error: The lengths of x and y do not match.");
        
        int i = 0;
        double d_buf = 0.0, y_mu = 0.0, x_mu = 0.0, x_var = 0.0, y_var = 0.0, cov = 0.0;
        for (i = 0; i < N; i++) {
            x_mu += x[i];
            y_mu += y[i];
        }
        x_mu /= (double) N;
        y_mu /= (double) N;
        for (i = 0; i < N; i++) {
            d_buf = (x[i] - x_mu);
            x_var += d_buf*d_buf;
            d_buf = (y[i] - y_mu);
            y_var += d_buf*d_buf;
        }
        x_var /= (double) (N - 1.0);
        y_var /= (double) (N - 1.0);
        for (i = 0; i < N; i++) cov += (x[i] - x_mu)*(y[i] - y_mu);
        cov /= (double) (N - 1);
        double a = 0.0, b = 0.0, sse = 0.0, a_se = 0.0, b_se = 0.0, p = 0.0, rsq = 0.0, r = 0.0;
        if (x_var > 0.0) b = cov / x_var;
        a = y_mu - b*x_mu;
        for (i = 0; i < N; i++) {
            d_buf = y[i] - a - b * x[i];
            sse += d_buf*d_buf;
        }
        if (x_var > 0.0) {
            a_se = sqrt((sse / (N - 2.0))*(1.0 / N + x_mu * x_mu / (x_var * (N - 1.0))));
            b_se = sqrt(sse / x_var / (N - 1.0) / (N - 2.0));
        }
        if (x_var > 0.0 && y_var > 0.0) {
            r = cov / sqrt(y_var * x_var);
            rsq = r*r;
        }
        double t = 0.0;
        if (b_se > 0.0) t = fabs(b / b_se);
        p = StatFunc::t_prob(N - 2.0, t, true);
        rst.clear();
        rst.push_back(b);
        rst.push_back(b_se);
        rst.push_back(p);
        rst.push_back(rsq);
        rst.push_back(r);
        
        eigenMatrix reg_sum(3, 3);
        if (table) {
            reg_sum(2, 0) = rsq;
            reg_sum(1, 0) = b;
            reg_sum(1, 1) = b_se;
            reg_sum(1, 2) = p;
            if (a_se > 0.0) t = fabs(a / a_se);
            p = StatFunc::t_prob(N - 2.0, t, true);
            reg_sum(0, 0) = a;
            reg_sum(0, 1) = a_se;
            reg_sum(0, 2) = p;
            return (reg_sum);
        }
        return (reg_sum);
    }    

  
    bool make_XMat(bInfo* bdata, MatrixXf &X)
    {
        if (bdata->_mu.empty()) calcu_mu(bdata);
        
        cout << "Recoding genotypes (individual major mode) ..." << endl;
        bool have_mis = false;
        unsigned long i = 0, j = 0, n = bdata->_keep.size(), m = bdata->_include.size();
        
        X.resize(0,0);
        X.resize(n, m);
//#pragma omp parallel for private(j)
        for (i = 0; i < n; i++) {
            if (bdata->_dosage_flag) {
                for (j = 0; j < m; j++) {
                    if (bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]] < 1e5) {
                        if (bdata->_allele1[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) X(i,j) = bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]];
                        else X(i,j) = 2.0 - bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]];
                    }
                    else {
                        X(i,j) = 1e6;
                        have_mis = true;
                    }
                }
                bdata->_geno_dose[i].clear();
            }
            else {
                for (j = 0; j < bdata->_include.size(); j++) {
                    if (!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]) {
                        if (bdata->_allele1[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) X(i,j) = bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]];
                        else X(i,j) = 2.0 - (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                    }
                    else {
                        X(i,j) = 1e6;
                        have_mis = true;
                    }
                }
            }
        }
        return have_mis;
    }
   
    void makex_eigenVector(bInfo* bdata,int j, eigenVector &x, bool resize, bool minus_2p)
    {
        int i = 0;
        if (resize) x.resize(bdata->_keep.size());
        for (i = 0; i < bdata->_keep.size(); i++) {
            if (!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]) {
                if (bdata->_allele1[bdata->_include[j]] == bdata->_ref_A[bdata->_include[j]]) x[i] = (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                else x[i] = 2.0 - (bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
            } else x[i] = bdata->_mu[bdata->_include[j]];
            if (minus_2p) x[i] -= bdata->_mu[bdata->_include[j]];
        }
    }
    
       
    void makeptrx(bInfo* bdata,int bsnpid,int cursnpid, float* X, bool minus_2p)
    {
        int i = 0;
        for (i = 0; i < bdata->_keep.size(); i++) {
            if (!bdata->_snp_1[bdata->_include[bsnpid]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[bsnpid]][bdata->_keep[i]]) {
                if (bdata->_allele1[bdata->_include[bsnpid]] == bdata->_ref_A[bdata->_include[bsnpid]]) X[cursnpid*bdata->_indi_num+i] = (bdata->_snp_1[bdata->_include[bsnpid]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[bsnpid]][bdata->_keep[i]]);
                else X[cursnpid*bdata->_indi_num+i] = 2.0 - (bdata->_snp_1[bdata->_include[bsnpid]][bdata->_keep[i]] + bdata->_snp_2[bdata->_include[bsnpid]][bdata->_keep[i]]);
            } else X[cursnpid*bdata->_indi_num+i] = bdata->_mu[bdata->_include[bsnpid]];
            if (minus_2p) X[cursnpid*bdata->_indi_num+i] -= bdata->_mu[bdata->_include[bsnpid]];
        }
    }
    
    void read_indi_list(string indi_list_file, vector<string> &indi_list)
    {
        ifstream i_indi_list(indi_list_file.c_str());
        if(!i_indi_list) throw("Error: can not open the file ["+indi_list_file+"] to read.");
        string str_buf, id_buf;
        indi_list.clear();
        while(i_indi_list){
            i_indi_list>>str_buf;
            if(i_indi_list.eof()) break;
            id_buf=str_buf+":";
            i_indi_list>>str_buf;
            id_buf+=str_buf;
            indi_list.push_back(id_buf);
            getline(i_indi_list, str_buf);
        }
        i_indi_list.close();
    }
    
    void read_msglist(string msglistfile, vector<string> &msglist, string msg)
    {
        // Read msglist file
        msglist.clear();
        string StrBuf;
        ifstream i_msglist(msglistfile.c_str());
        if(!i_msglist) throw("Error: can not open the file ["+msglistfile+"] to read.");
        cout<<"Reading a list of "<<msg<<" from ["+msglistfile+"]."<<endl;
        while(i_msglist>>StrBuf){
            msglist.push_back(StrBuf);
            getline(i_msglist, StrBuf);
        }
        i_msglist.close();
    }
    
    
    void update_id_map_kp(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep)
    {
        int i=0;
        map<string, int> id_map_buf(id_map);
        for(i=0; i<id_list.size(); i++) id_map_buf.erase(id_list[i]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) id_map.erase(iter->first);
        
        keep.clear();
        for(iter=id_map.begin(); iter!=id_map.end(); iter++) keep.push_back(iter->second);
        stable_sort(keep.begin(), keep.end());
    }
    void update_id_map_rm(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep)
    {
        int i = 0;
        for (i = 0; i < id_list.size(); i++) id_map.erase(id_list[i]);
        
        keep.clear();
        map<string, int>::iterator iter;
        for (iter = id_map.begin(); iter != id_map.end(); iter++) keep.push_back(iter->second);
        stable_sort(keep.begin(), keep.end());
    }
    

    void keep_indi(bInfo* bdata,string indi_list_file)
    {
        vector<string> indi_list;
        read_indi_list(indi_list_file, indi_list);
        update_id_map_kp(indi_list, bdata->_id_map, bdata->_keep);
        cout<<bdata->_keep.size()<<" individuals are kept from ["+indi_list_file+"]."<<endl;
    }
    
    void remove_indi(bInfo* bdata, string indi_list_file) {
        vector<string> indi_list;
        read_indi_list(indi_list_file, indi_list);
        int prev_size = bdata->_keep.size();
        update_id_map_rm(indi_list, bdata->_id_map, bdata->_keep);
        cout << prev_size - bdata->_keep.size() << " individuals are removed from [" + indi_list_file + "] and there are " << bdata->_keep.size() << " individuals remaining." << endl;
    }
    void extract_region_bp(bInfo* bdata, int chr, int fromkb, int tokb)
    {
        int frombp=fromkb*1000;
        int tobp=tokb*1000;
        vector<string> snplist;
        for(int i = 0; i < bdata->_include.size(); i++){
            int j = bdata->_include[i];
            if(bdata->_chr[j] == chr && bdata->_bp[j]<=tobp && bdata->_bp[j]>=frombp) snplist.push_back(bdata->_snp_name[j]);
        }
        if(snplist.empty()) throw ("Error: on SNP found in this region.");
        update_id_map_kp(snplist, bdata->_snp_name_map, bdata->_include);
        cout << bdata->_include.size() << " SNPs are extracted from SNP BP: "<<fromkb<<"Kb"<<"to SNP BP: "<<tokb<<"Kb."<< endl;
    }
    void extract_snp(bInfo* bdata,string snplistfile)
    {
        vector<string> snplist;
        string msg="SNPs";
        read_msglist(snplistfile, snplist,msg);
        update_id_map_kp(snplist, bdata->_snp_name_map, bdata->_include);
        cout<<bdata->_include.size()<<" SNPs are extracted from ["+snplistfile+"]."<<endl;
    }
    void exclude_snp(bInfo* bdata,string snplistfile)
    {
        vector<string> snplist;
         string msg="SNPs";
        read_msglist(snplistfile, snplist,msg);
        int prev_size = bdata->_include.size();
        update_id_map_rm(snplist, bdata->_snp_name_map, bdata->_include);
        cout << prev_size - bdata->_include.size() << " SNPs are excluded from [" + snplistfile + "] and there are " << bdata->_include.size() << " SNPs remaining." << endl;
    }
   
    void extract_eqtl_by_chr(eqtlInfo* eqtlinfo, int snpchr)
    {
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            int tmpint=eqtlinfo->_esi_include[i];
            if(eqtlinfo->_esi_chr[tmpint]==snpchr ) newIcld.push_back(tmpint);
        }
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_esi_include=newIcld;
        cout << eqtlinfo->_esi_include.size() << " SNPs are extracted from chromosome [" + atos(snpchr) + "]." << endl;
    }
    void extract_epi_by_chr(eqtlInfo* eqtlinfo, int prbchr)
    {
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            int tmpint=eqtlinfo->_include[i];
            if(eqtlinfo->_epi_chr[tmpint]==prbchr ) newIcld.push_back(tmpint);
        }
        eqtlinfo->_include.clear();
        eqtlinfo->_include=newIcld;
        cout << eqtlinfo->_include.size() << " probes are extracted from chromosome [" + atos(prbchr) + "]." << endl;
    }

	void extract_eqtl_snp(eqtlInfo* eqtlinfo, string snplstName)
	{
		vector<string> snplist;
		string msg = "SNPs";
		read_msglist(snplstName, snplist, msg);
        if(eqtlinfo->_esi_include.size()==eqtlinfo->_snpNum)
        {
            eqtlinfo->_esi_include.clear();
            StrFunc::match_only(snplist, eqtlinfo->_esi_rs, eqtlinfo->_esi_include);
            stable_sort(eqtlinfo->_esi_include.begin(), eqtlinfo->_esi_include.end());
        }
        else
        {
            vector<int> icld_step1;
            StrFunc::match_only(snplist, eqtlinfo->_esi_rs, icld_step1);
            
            vector<int> common_probes;
            set_intersect(icld_step1, eqtlinfo->_esi_include, common_probes);
            eqtlinfo->_esi_include.clear();
            eqtlinfo->_esi_include=common_probes;
            stable_sort(eqtlinfo->_esi_include.begin(), eqtlinfo->_esi_include.end());
        }
        
		cout << eqtlinfo->_esi_include.size() << " SNPs are extracted from [" + snplstName + "]." << endl;
	}
   
    void extract_eqtl_snp(eqtlInfo* eqtlinfo, string snporprb, int Wind, string msg)
    {
        string logstr;
        
        int bp=-9;
        int chr=-9;
        if(msg=="SNP")
        {
            long idx=find(eqtlinfo->_esi_rs.begin(), eqtlinfo->_esi_rs.end(), snporprb)-eqtlinfo->_esi_rs.begin();
            if(idx==eqtlinfo->_esi_rs.size())
            {
                logstr="ERROR: Can't find SNP "+snporprb+" in the dataset. Please check.\n";
                fputs(logstr.c_str(),stdout);
                exit(1);
            }
            bp=eqtlinfo->_esi_bp[idx];
            chr=eqtlinfo->_esi_chr[idx];
        }
        else if(msg=="probe")
        {
            long idx=-9;
            for(int i=0;i<eqtlinfo->_include.size();i++)
            {
                if(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]==snporprb) idx=eqtlinfo->_include[i];
            }
            if(idx==-9)
            {
                logstr="ERROR: Can't find probe "+snporprb+" in the dataset. Please check.\n";
                fputs(logstr.c_str(),stdout);
                exit(1);
            }
            bp=eqtlinfo->_epi_bp[idx];
            chr=eqtlinfo->_epi_chr[idx];
        }
        
        int upbound=bp+Wind*1000;
        int tmpint=bp-Wind*1000;
        int lowbound=tmpint>0?tmpint:0;
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            tmpint=eqtlinfo->_esi_include[i];
            if(eqtlinfo->_esi_chr[tmpint]==chr && eqtlinfo->_esi_bp[tmpint]>=lowbound && eqtlinfo->_esi_bp[tmpint]<=upbound) newIcld.push_back(tmpint);
        }
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_esi_include=newIcld;
        cout << eqtlinfo->_esi_include.size() << " SNPs are extracted from the region: " +atos(Wind)+" Kb around [" + snporprb + "]." << endl;
    }
    void extract_eqtl_single_snp(eqtlInfo* eqtlinfo, string snprs)
    {
        string logstr;
        long idx=find(eqtlinfo->_esi_rs.begin(), eqtlinfo->_esi_rs.end(), snprs)-eqtlinfo->_esi_rs.begin();
        if(idx==eqtlinfo->_esi_rs.size())
        {
            logstr="ERROR: Can't find SNP "+snprs+" in the dataset. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_esi_include.push_back((int)idx);
        cout << snprs << " is extracted. " << endl;
    }
    void extract_eqtl_snp(eqtlInfo* eqtlinfo, string fromsnprs, string tosnprs)
    {
        string logstr;
        long fromidx=find(eqtlinfo->_esi_rs.begin(), eqtlinfo->_esi_rs.end(), fromsnprs)-eqtlinfo->_esi_rs.begin();
        if(fromidx==eqtlinfo->_esi_rs.size())
        {
            logstr="ERROR: Can't find SNP "+fromsnprs+" in the dataset. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        int fromsnpbp=eqtlinfo->_esi_bp[fromidx];
        int snpchr=eqtlinfo->_esi_chr[fromidx];
        
        long toidx=find(eqtlinfo->_esi_rs.begin(), eqtlinfo->_esi_rs.end(), tosnprs)-eqtlinfo->_esi_rs.begin();
        if(toidx==eqtlinfo->_esi_rs.size())
        {
            logstr="ERROR: Can't find SNP "+tosnprs+" in the dataset. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        int tosnpbp=eqtlinfo->_esi_bp[toidx];
        int tosnpchr=eqtlinfo->_esi_chr[toidx];
        if(tosnpchr != snpchr)
        {
            logstr="ERROR: SNP "+fromsnprs+" and SNP "+tosnprs +" are not from the same chromosome. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        
        if(fromsnpbp>tosnpbp)
        {
            int tmp=fromsnpbp;
            fromsnpbp=tosnpbp;
            tosnpbp=tmp;
        }
     
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            int tmpint=eqtlinfo->_esi_include[i];
            if(eqtlinfo->_esi_chr[tmpint]==snpchr && eqtlinfo->_esi_bp[tmpint]>=fromsnpbp && eqtlinfo->_esi_bp[tmpint]<=tosnpbp) newIcld.push_back(tmpint);
        }
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_esi_include=newIcld;
        cout << eqtlinfo->_esi_include.size() << " SNPs are extracted from SNP " +fromsnprs+" to SNP " + tosnprs + "." << endl;
    }
    void extract_eqtl_snp(eqtlInfo* eqtlinfo, int chr, int fromsnpkb, int tosnpkb)
    {
        int fromsnpbp=fromsnpkb*1000;
        int tosnpbp=tosnpkb*1000;
        
        if(fromsnpbp>tosnpbp)
        {
            int tmp=fromsnpbp;
            fromsnpbp=tosnpbp;
            tosnpbp=tmp;
        }
        
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_esi_include.size();i++)
        {
            int tmpint=eqtlinfo->_esi_include[i];
            if( eqtlinfo->_esi_chr[tmpint]==chr &&eqtlinfo->_esi_bp[tmpint]>=fromsnpbp && eqtlinfo->_esi_bp[tmpint]<=tosnpbp) newIcld.push_back(tmpint);
        }
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_esi_include=newIcld;
        cout << eqtlinfo->_esi_include.size() << " SNPs are extracted from SNP BP: " +atos(fromsnpkb)+" Kb to SNP BP: " + atos(tosnpkb) + " Kb." << endl;
    }



    void exclude_eqtl_snp(eqtlInfo* eqtlinfo, string snplstName)
    {
        vector<string> snplist;
        vector<string> mapstr;
         vector<int> tmp;
        string msg = "SNPs";
        read_msglist(snplstName, snplist, msg);
        int pre_num=eqtlinfo->_esi_include.size();
        mapstr.resize(pre_num);
        tmp.resize(pre_num);
        for(int i=0;i<pre_num;i++)
        {
            mapstr[i]=eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]];
            tmp[i]=eqtlinfo->_esi_include[i];
        }
        
        eqtlinfo->_esi_include.clear();
        StrFunc::set_complement(snplist, mapstr, tmp, eqtlinfo->_esi_include); //sorted
        stable_sort(eqtlinfo->_esi_include.begin(), eqtlinfo->_esi_include.end());       
        cout << pre_num-eqtlinfo->_esi_include.size() << " SNPs are excluded from [" + snplstName + "] and there are " << eqtlinfo->_esi_include.size() << " SNPs remaining." << endl;
    }
    
    void extract_gwas_snp(gwasData* gdata, string snplstName)
    {
        vector<string> snplist;
        string msg = "SNPs";
        read_msglist(snplstName, snplist, msg);
        gdata->_include.clear();
        StrFunc::match_only(snplist, gdata->snpName, gdata->_include);
        stable_sort(gdata->_include.begin(), gdata->_include.end());
        cout << gdata->_include.size() << " SNPs are extracted from [" + snplstName + "]." << endl;
    }
    
    void exclude_gwas_snp(gwasData* gdata, string snplstName)
    {
        vector<string> snplist;
        vector<string> mapstr;
         vector<int> tmp;
        string msg = "SNPs";
        read_msglist(snplstName, snplist, msg);
        int pre_num=gdata->_include.size();
        mapstr.resize(pre_num);
         tmp.resize(pre_num);
        for(int i=0;i<pre_num;i++)
        {
            mapstr[i]=gdata->snpName[gdata->_include[i]];
            tmp[i]=gdata->_include[i];
        }
        
        gdata->_include.clear();
        StrFunc::set_complement(snplist, mapstr, tmp,gdata->_include); //sorted
        stable_sort(gdata->_include.begin(), gdata->_include.end());
        cout << pre_num-gdata->_include.size() << " SNPs are excluded from [" + snplstName + "] and there are " << gdata->_include.size() << " SNPs remaining." << endl;
    }

    
    void extract_prob(eqtlInfo* eqtlinfo,string problstName)
    {
        vector<string> problist;
        string msg="probes";
        read_msglist(problstName, problist,msg);
        vector<string> raw_problist;
        for(int i=0;i<eqtlinfo->_include.size();i++) raw_problist.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]);
        vector<string> common_probes;
        set_intersect(problist, raw_problist, common_probes);
        eqtlinfo->_include.clear();
		StrFunc::match_only(common_probes, eqtlinfo->_epi_prbID, eqtlinfo->_include);
        stable_sort(eqtlinfo->_include.begin(), eqtlinfo->_include.end());
        cout<<eqtlinfo->_include.size()<<" probes are extracted from ["+problstName+"]."<<endl;
    }
    
    void extract_prob_by_gene(eqtlInfo* eqtlinfo, string genelistName)
    {
        vector<string> genelist;
        string msg="genes";
        read_msglist(genelistName, genelist,msg);
        vector<string> raw_problist;
        vector<string> raw_genelist;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            raw_problist.push_back(eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]);
            raw_genelist.push_back(eqtlinfo->_epi_gene[eqtlinfo->_include[i]]);
        }
        vector<int> idx;
        for(int i=0;i<genelist.size();i++)
        {
            string tmpname1=genelist[i];
            for(int j=0;j<raw_genelist.size();j++)
            {
                string tmpname2=raw_genelist[j];
                if(tmpname1==tmpname2)  idx.push_back(j);
                else
                {
                    vector<string> substrs;
                    int tmpnum=split_string(tmpname2,substrs);
                    if(tmpnum>1)
                        for(int k=0;k<tmpnum;k++)
                            if(tmpname1==substrs[k])
                            {
                                idx.push_back(j);
                                break;
                            }
                }
            }
            out2:;
        }
        vector<string> common_probes;
        for(int i=0;i<idx.size();i++) common_probes.push_back(raw_problist[idx[i]]);
        eqtlinfo->_include.clear();
        StrFunc::match_only(common_probes, eqtlinfo->_epi_prbID, eqtlinfo->_include);
        stable_sort(eqtlinfo->_include.begin(), eqtlinfo->_include.end());
        cout<<eqtlinfo->_include.size()<<" probes are extracted from ["+genelistName+"]."<<endl;
    }
    
    void extract_prob(eqtlInfo* eqtlinfo, string prbname, int prbWind)
    {
        string logstr;
        long idx=find(eqtlinfo->_epi_prbID.begin(), eqtlinfo->_epi_prbID.end(), prbname)-eqtlinfo->_epi_prbID.begin();
        if(idx==eqtlinfo->_epi_prbID.size())
        {
            logstr="ERROR: Can't find probe "+prbname+" in the dataset. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        int prbbp=eqtlinfo->_epi_bp[idx];
        int prbchr=eqtlinfo->_epi_chr[idx];
        int upbound=prbbp+prbWind*1000;
        int tmpint=prbbp-prbWind*1000;
        int lowbound=tmpint>0?tmpint:0;
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            tmpint=eqtlinfo->_include[i];
            if(eqtlinfo->_epi_chr[tmpint]==prbchr && eqtlinfo->_epi_bp[tmpint]>=lowbound && eqtlinfo->_epi_bp[tmpint]<=upbound) newIcld.push_back(tmpint);
        }
        eqtlinfo->_include.clear();
        eqtlinfo->_include=newIcld;
        cout << eqtlinfo->_include.size() << " probes are extracted from the region: " +atos(prbWind)+" Kb around [" + prbname + "]." << endl;
    }
    void extract_eqtl_single_probe(eqtlInfo* eqtlinfo, string prbname, bool prtscr)
    {
        string logstr;
        int idx=find(eqtlinfo->_epi_prbID.begin(), eqtlinfo->_epi_prbID.end(), prbname)-eqtlinfo->_epi_prbID.begin();
        if(idx==eqtlinfo->_epi_prbID.size())
        {
            logstr="ERROR: Can't find probe "+prbname+" in the dataset. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }

        eqtlinfo->_include.clear();
        eqtlinfo->_include.push_back(idx);
        if(prtscr) cout << prbname << " is extracted. " << endl;
    }
    void extract_eqtl_prob(eqtlInfo* eqtlinfo, string fromprbname, string toprbname)
    {
        string logstr;
        long fromidx=find(eqtlinfo->_epi_prbID.begin(), eqtlinfo->_epi_prbID.end(), fromprbname)-eqtlinfo->_epi_prbID.begin();
        if(fromidx==eqtlinfo->_epi_prbID.size())
        {
            logstr="ERROR: Can't find SNP "+fromprbname+" in the dataset. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        int fromprbbp=eqtlinfo->_epi_bp[fromidx];
        int prbchr=eqtlinfo->_epi_chr[fromidx];
        
        long toidx=find(eqtlinfo->_epi_prbID.begin(), eqtlinfo->_epi_prbID.end(), toprbname)-eqtlinfo->_epi_prbID.begin();
        if(toidx==eqtlinfo->_epi_prbID.size())
        {
            logstr="ERROR: Can't find SNP "+toprbname+" in the dataset. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        int toprbbp=eqtlinfo->_epi_bp[toidx];
        int toprbchr=eqtlinfo->_epi_chr[toidx];
        if(toprbchr != prbchr)
        {
            logstr="ERROR: probe "+fromprbname+" and probe "+toprbname +" are not from the same chromosome. Please check.\n";
            fputs(logstr.c_str(),stdout);
            exit(1);
        }
        
        if(fromprbbp>toprbbp)
        {
            int tmp=fromprbbp;
            fromprbbp=toprbbp;
            toprbbp=tmp;
        }
        
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            int tmpint=eqtlinfo->_include[i];
            if(eqtlinfo->_epi_chr[tmpint]==prbchr && eqtlinfo->_epi_bp[tmpint]>=fromprbbp && eqtlinfo->_epi_bp[tmpint]<=toprbbp) newIcld.push_back(tmpint);
        }
        eqtlinfo->_include.clear();
        eqtlinfo->_include=newIcld;
        cout << eqtlinfo->_include.size() << " probes are extracted from probe : " +fromprbname+" to probe " + toprbname + "." << endl;
    }
    void extract_eqtl_prob(eqtlInfo* eqtlinfo, int chr, int fromprbkb, int toprbkb)
    {
        int fromprbbp=fromprbkb*1000;
        int toprbbp=toprbkb*1000;
        
        if(fromprbbp>toprbbp)
        {
            int tmp=fromprbbp;
            fromprbbp=toprbbp;
            toprbbp=tmp;
        }
        
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            int tmpint=eqtlinfo->_include[i];
            if( eqtlinfo->_epi_chr[tmpint]==chr && eqtlinfo->_epi_bp[tmpint]>=fromprbbp && eqtlinfo->_epi_bp[tmpint]<=toprbbp) newIcld.push_back(tmpint);
        }
        eqtlinfo->_include.clear();
        eqtlinfo->_include=newIcld;
        cout << eqtlinfo->_include.size() << " probes are extracted from probe BP: " +atos(fromprbkb)+"Kb to probe BP: " + atos(toprbkb) + "Kb." << endl;
    }
    void extract_prob_by_single_gene(eqtlInfo* eqtlinfo, string genename)
    {
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            int tmpint=eqtlinfo->_include[i];
            if( eqtlinfo->_epi_gene[tmpint]==genename) newIcld.push_back(tmpint);
        }
        eqtlinfo->_include.clear();
        eqtlinfo->_include=newIcld;
        cout << eqtlinfo->_include.size() << " probes are extracted that are mapping to gene: " +genename+ "." << endl;
    }

    void exclude_prob(eqtlInfo* eqtlinfo,string problstName)
    {
        vector<string> problist;
        vector<string> mappro;
        vector<int> tmp;
        string msg="probes";
        read_msglist(problstName, problist,msg);
        int pre_num=eqtlinfo->_include.size();
        mappro.resize(pre_num);
        tmp.resize(pre_num);
        for(int i=0;i<pre_num;i++)
        {
            mappro[i]=eqtlinfo->_epi_prbID[eqtlinfo->_include[i]];
            tmp[i]=eqtlinfo->_include[i];
        }
       
        eqtlinfo->_include.clear();
        StrFunc::set_complement(problist, mappro, tmp, eqtlinfo->_include);
        
        cout<<pre_num-eqtlinfo->_include.size()<<" probes are excluded from ["+problstName+"]and there are "<<eqtlinfo->_include.size()<<" probes remaining."<<endl;
    }
    
    void exclude_eqtl_single_probe(eqtlInfo* eqtlinfo, string prbname)
    {
        string logstr;
        vector<int> newIcld;
        for(int i=0;i<eqtlinfo->_include.size();i++)
        {
            int idx=eqtlinfo->_include[i];
            if(eqtlinfo->_epi_prbID[idx] != prbname) newIcld.push_back(idx);
        }
       
        if(newIcld.size()==eqtlinfo->_include.size())
        {
            logstr="WARNING: Can't find probe "+prbname+" in the dataset. Nothing to exclude.\n";
            fputs(logstr.c_str(),stdout);
            
        }else
        {
            eqtlinfo->_include.clear();
            eqtlinfo->_include=newIcld;
            cout << prbname << " is excluded from the dataset. " << endl;
        }
        
    }

    void filter_snp_maf(bInfo* bdata,double maf)
    {
        if(bdata->_mu.empty()) calcu_mu(bdata);
        
        cout<<"Pruning SNPs with MAF > "<<maf<<" ..."<<endl;
        map<string, int> id_map_buf(bdata->_snp_name_map);
        map<string, int>::iterator iter, end=id_map_buf.end();
        int prev_size=bdata->_include.size();
        double fbuf=0.0;
        bdata->_include.clear();
        bdata->_snp_name_map.clear();
        for(iter=id_map_buf.begin(); iter!=end; iter++){
            fbuf=bdata->_mu[iter->second]*0.5;
            if(fbuf<=maf || (1.0-fbuf)<=maf) continue;
            bdata->_snp_name_map.insert(*iter);
            bdata->_include.push_back(iter->second);
        }
        if(bdata->_include.size()==0) throw("Error: No SNP is retained for analysis.");
        else{
            stable_sort(bdata->_include.begin(), bdata->_include.end());
            cout<<"After pruning SNPs with MAF > "<<maf<<", there are "<<bdata->_include.size()<<" SNPs ("<<prev_size-bdata->_include.size()<<" SNPs with MAF < "<<maf<<")."<<endl;
        }
        
    }
    
    void allele_check(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Performing allele check between GWAS summary dataset, eQTL summary dataset and Plink dataset. ";
        cout<<logstr<<endl;
        
        vector<string> bsnp;
        vector<string> essnp;
        if(bdata->_include.size()< bdata->_snp_num || esdata->_esi_include.size()<esdata->_snpNum )
        {
            for(int i=0;i<bdata->_include.size();i++) bsnp.push_back(bdata->_snp_name[bdata->_include[i]]);
            for(int i=0;i<esdata->_esi_include.size();i++) essnp.push_back(esdata->_esi_rs[esdata->_esi_include[i]]);
            StrFunc::match_only(bsnp, gdata->snpName, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
            for(int i=0;i<edId.size();i++) cmmnSNPs.push_back(gdata->snpName[edId[i]]);
            StrFunc::match_only(cmmnSNPs, essnp, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            for(int i=0;i<edId.size();i++) edId[i]=esdata->_esi_include[edId[i]];
            for(int i=0;i<edId.size();i++) slctSNPs.push_back(esdata->_esi_rs[edId[i]]);
        }else
        {
            StrFunc::match_only(bdata->_snp_name, gdata->snpName, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
            for(int i=0;i<edId.size();i++) cmmnSNPs.push_back(gdata->snpName[edId[i]]);
            edId.clear();
            StrFunc::match_only(cmmnSNPs, esdata->_esi_rs, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            for(int i=0;i<edId.size();i++) slctSNPs.push_back(esdata->_esi_rs[edId[i]]);
        }
        
        //slctSNPs is in the order as bdata. so bdId is in increase order. edId and gdId may not.
        
        //alleles check
        StrFunc::match(slctSNPs, bdata->_snp_name, bdId);
        StrFunc::match(slctSNPs, gdata->snpName, gdId);
        cmmnSNPs.clear();
        bdata->_include.clear();
        gdata->_include.clear();
        esdata->_esi_include.clear();
        
        for (int i = 0; i<edId.size(); i++)
        {
            string a1, a2, ga1, ga2, ea1, ea2;
            a1 = bdata->_allele1[bdId[i]];
            a2 = bdata->_allele2[bdId[i]];
            ga1 = gdata->allele_1[gdId[i]];
            ga2 = gdata->allele_2[gdId[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
            // use the allele in eQTL summary data as the reference allele. so we won't get the whole besd into memroy
            if(ea1 == a1 &&  ea2 == a2)
            {
                if( ea1 == ga1 && ea2 == ga2)
                {
                    
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                }
                else if(ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    
                }
            }
            else if(ea1 == a2 &&  ea2 == a1)
            {
                
                if( ea1 == ga1 && ea2 == ga2)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                   
                }
                else if(ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
            }

            
        }
        
        logstr=itos(bdata->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
        
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }
    
    bool allele_check_(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Performing allele check between GWAS summary dataset, eQTL summary dataset and Plink dataset. ";
        cout<<logstr<<endl;
        
        vector<string> bsnp;
        vector<string> essnp;
        if(bdata->_include.size()< bdata->_snp_num || esdata->_esi_include.size()<esdata->_snpNum )
        {
            for(int i=0;i<bdata->_include.size();i++) bsnp.push_back(bdata->_snp_name[bdata->_include[i]]);
            for(int i=0;i<esdata->_esi_include.size();i++) essnp.push_back(esdata->_esi_rs[esdata->_esi_include[i]]);
            StrFunc::match_only(bsnp, gdata->snpName, edId);
            if(edId.empty()) return false;
            for(int i=0;i<edId.size();i++) cmmnSNPs.push_back(gdata->snpName[edId[i]]);
            StrFunc::match_only(cmmnSNPs, essnp, edId);
            if(edId.empty()) return false;
            for(int i=0;i<edId.size();i++) edId[i]=esdata->_esi_include[edId[i]];
            for(int i=0;i<edId.size();i++) slctSNPs.push_back(esdata->_esi_rs[edId[i]]);
        }else
        {
            StrFunc::match_only(bdata->_snp_name, gdata->snpName, edId);
            if(edId.empty()) return false;
            for(int i=0;i<edId.size();i++) cmmnSNPs.push_back(gdata->snpName[edId[i]]);
            edId.clear();
            StrFunc::match_only(cmmnSNPs, esdata->_esi_rs, edId);
            if(edId.empty()) return false;
            for(int i=0;i<edId.size();i++) slctSNPs.push_back(esdata->_esi_rs[edId[i]]);
        }
        
        //slctSNPs is in the order as bdata. so bdId is in increase order. edId and gdId may not.
        
        //alleles check
        StrFunc::match(slctSNPs, bdata->_snp_name, bdId);
        StrFunc::match(slctSNPs, gdata->snpName, gdId);
        cmmnSNPs.clear();
        bdata->_include.clear();
        gdata->_include.clear();
        esdata->_esi_include.clear();
        
        for (int i = 0; i<edId.size(); i++)
        {
            string a1, a2, ga1, ga2, ea1, ea2;
            a1 = bdata->_allele1[bdId[i]];
            a2 = bdata->_allele2[bdId[i]];
            ga1 = gdata->allele_1[gdId[i]];
            ga2 = gdata->allele_2[gdId[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
            // use the allele in eQTL summary data as the reference allele. so we won't get the whole besd into memroy
            if(ea1 == a1 &&  ea2 == a2)
            {
                if( ea1 == ga1 && ea2 == ga2)
                {
                    
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                }
                else if(ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    
                }
            }
            else if(ea1 == a2 &&  ea2 == a1)
            {
                
                if( ea1 == ga1 && ea2 == ga2)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
                else if(ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
            }
            
            
        }
        
        logstr=itos(bdata->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
        
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
        
        return true;
        
    }
    
    void allele_check(gwasData* gdata, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> cmmnSNPs;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Performing allele check between GWAS summary dataset and eQTL summary dataset. ";
        cout<<logstr<<endl;
       
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
        cmmnSNPs.clear();
        gdata->_include.clear();
        esdata->_esi_include.clear();
        
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
               
                gdata->_include.push_back(gdId[i]);
                esdata->_esi_include.push_back(edId[i]);
            }
            else if(ea1 == ga2 && ea2 == ga1)
            {
                gdata->_include.push_back(gdId[i]);
                esdata->_esi_include.push_back(edId[i]);
                
                gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                
            }
        }
        logstr=itos(gdata->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
        
    }
    
    void allele_check(eqtlInfo* etrait, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> cmmnSNPs;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Performing allele check between GWAS summary dataset and eQTL summary dataset. ";
        cout<<logstr<<endl;
         vector<string> etsnp;
        vector<string> essnp;
        if(esdata->_esi_include.size()<esdata->_snpNum || etrait->_esi_include.size()<etrait->_snpNum)
        {
            for(int i=0;i<esdata->_esi_include.size();i++) essnp.push_back(esdata->_esi_rs[esdata->_esi_include[i]]);
            for(int i=0;i<etrait->_esi_include.size();i++) etsnp.push_back(etrait->_esi_rs[etrait->_esi_include[i]]);
            StrFunc::match_only(essnp, etsnp, gdId);
            if(gdId.empty()) throw("Error: no common SNPs found.");
            for(int i=0;i<gdId.size();i++) slctSNPs.push_back(etsnp[gdId[i]]);
        }else
        {
            StrFunc::match_only(esdata->_esi_rs, etrait->_esi_rs, gdId);
            if(gdId.empty()) throw("Error: no common SNPs found.");
            for(int i=0;i<gdId.size();i++) slctSNPs.push_back(etrait->_esi_rs[gdId[i]]);
        }
        
        
        //alleles check
        StrFunc::match(slctSNPs, esdata->_esi_rs, edId);
        gdId.clear();
        StrFunc::match(slctSNPs, etrait->_esi_rs, gdId);
        cmmnSNPs.clear();
        etrait->_esi_include.clear();
        esdata->_esi_include.clear();
        
        for (int i = 0; i<edId.size(); i++)
        {
            string ta1, ta2, ea1, ea2;
            
            ta1 = etrait->_esi_allele1[gdId[i]];
            ta2 = etrait->_esi_allele2[gdId[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
            // use the allele in eQTL summary data as the reference allele. so we won't get the whole besd into memroy
            if( ea1 == ta1 && ea2 == ta2)
            {
                
                etrait->_esi_include.push_back(gdId[i]);
                esdata->_esi_include.push_back(edId[i]);
            }
            else if(ea1 == ta2 && ea2 == ta1)
            {
                etrait->_esi_include.push_back(gdId[i]);
                esdata->_esi_include.push_back(edId[i]);
                
                if(etrait->_val.size()>0)
                {
                    int count=0;
                    for(int j=0;j<etrait->_rowid.size();j++)
                    {
                        if(etrait->_rowid[j]==gdId[i])
                        {
                            count++;
                            if(count & 1)
                                etrait->_val[j]=-etrait->_val[j];
                        }
                    }
                }
                else
                {
                    for(int j=0;j<etrait->_include.size();j++) if( etrait->_bxz[j][gdId[i]]+9 > 1e-6 ) etrait->_bxz[j][gdId[i]]=-etrait->_bxz[j][gdId[i]];
                }
                
            }
        }
        logstr=itos(etrait->_esi_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
        
    }
    
    void allele_check(gwasData* gdata1, gwasData* gdata2)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> cmmnSNPs;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Performing allele check between GWAS summary dataset and eQTL summary dataset. ";
        cout<<logstr<<endl;
        
        StrFunc::match_only(gdata2->snpName, gdata1->snpName, gdId);
        if(gdId.empty()) throw("Error: no common SNPs found.");
        
        for(int i=0;i<gdId.size();i++) slctSNPs.push_back(gdata1->snpName[gdId[i]]);
        
        //alleles check
        StrFunc::match(slctSNPs, gdata2->snpName, edId);
        cmmnSNPs.clear();
        gdata1->_include.clear();
        gdata2->_include.clear();
        
        for (int i = 0; i<edId.size(); i++)
        {
            string ga1, ga2, ea1, ea2;
            
            ga1 = gdata1->allele_1[gdId[i]];
            ga2 = gdata1->allele_2[gdId[i]];
            ea1 = gdata2->allele_1[edId[i]];
            ea2 = gdata2->allele_2[edId[i]];
            // use the allele in eQTL summary data as the reference allele. so we won't get the whole besd into memroy
            if( ea1 == ga1 && ea2 == ga2)
            {
                
                gdata1->_include.push_back(gdId[i]);
                gdata2->_include.push_back(edId[i]);
            }
            else if(ea1 == ga2 && ea2 == ga1)
            {
                gdata1->_include.push_back(gdId[i]);
                gdata2->_include.push_back(edId[i]);
                
                gdata1->byz[gdId[i]]=-gdata1->byz[gdId[i]];
            }
        }
        logstr=atos(gdata1->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
    }
    void allele_check(bInfo* bdata, eqtlInfo* etrait, eqtlInfo* esdata)
    {
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        string logstr="Performing allele check among GWAS summary dataset, eQTL summary dataset and Plink dataset. ";
        cout<<logstr<<endl;
        
        vector<string> bsnp(bdata->_include.size());
        vector<string> etsnp(etrait->_esi_include.size());
        vector<string> essnp(esdata->_esi_include.size());
        if(bdata->_include.size()< bdata->_snp_num || esdata->_esi_include.size()<esdata->_snpNum || etrait->_esi_include.size()<etrait->_snpNum )
        {
            #pragma omp parallel for
            for(int i=0;i<bdata->_include.size();i++)
                bsnp[i]=bdata->_snp_name[bdata->_include[i]];
            #pragma omp parallel for
            for(int i=0;i<etrait->_esi_include.size();i++)
                etsnp[i]=etrait->_esi_rs[etrait->_esi_include[i]];
            #pragma omp parallel for
            for(int i=0;i<esdata->_esi_include.size();i++)
                essnp[i]=esdata->_esi_rs[esdata->_esi_include[i]];
            match_only(bsnp, etsnp, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and eTrait data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etsnp[edId[i]];
            match_only(cmmnSNPs, essnp, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                edId[i]=esdata->_esi_include[edId[i]];
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=esdata->_esi_rs[edId[i]];
        }else
        {
            match_only(bdata->_snp_name, etrait->_esi_rs, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
            cmmnSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                cmmnSNPs[i]=etrait->_esi_rs[edId[i]];
            edId.clear();
            match_only(cmmnSNPs, esdata->_esi_rs, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
            slctSNPs.resize(edId.size());
            #pragma omp parallel for
            for(int i=0;i<edId.size();i++)
                slctSNPs[i]=esdata->_esi_rs[edId[i]];
        }
        
        //slctSNPs is in the order as bdata. so bdId is in increase order. edId and gdId may not.
        
        //alleles check
        match(slctSNPs, bdata->_snp_name, bdId);
        match(slctSNPs, etrait->_esi_rs, gdId);
        cmmnSNPs.clear();
        bdata->_include.clear();
        etrait->_esi_include.clear();
        esdata->_esi_include.clear();
        
        for (int i = 0; i<edId.size(); i++)
        {
            printf("Checking...  %3.0f%%\r", 100.0*i/edId.size());
            fflush(stdout);
            string a1, a2, ta1, ta2, ea1, ea2;
            a1 = bdata->_allele1[bdId[i]];
            a2 = bdata->_allele2[bdId[i]];
            ta1 = etrait->_esi_allele1[gdId[i]];
            ta2 = etrait->_esi_allele2[gdId[i]];
            ea1 = esdata->_esi_allele1[edId[i]];
            ea2 = esdata->_esi_allele2[edId[i]];
         
            // use the allele in eQTL summary data "esdata" as the reference allele. so we won't get the whole besd into memroy
            if(ea1 == a1 &&  ea2 == a2)
            {
                if( ea1 == ta1 && ea2 == ta2)
                {
                    
                    bdata->_include.push_back(bdId[i]);
                    etrait->_esi_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                }
                else if(ea1 == ta2 && ea2 == ta1)
                {
                    bdata->_include.push_back(bdId[i]);
                    etrait->_esi_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    if(etrait->_val.size()>0)
                    {
                        /*
                        int count=0;
                        for(int j=0;j<etrait->_rowid.size();j++)
                        {
                            if(etrait->_rowid[j]==gdId[i])
                            {
                                count++;
                                if(count & 1)
                                    etrait->_val[j]=-1*etrait->_val[j];
                            }
                        }
                         */
                        #pragma omp parallel for private(i)
                        for(int j=0;j<etrait->_probNum;j++)
                        {
                            uint64_t beta_start=etrait->_cols[j<<1];
                            uint64_t se_start=etrait->_cols[1+(j<<1)];
                            uint64_t numsnps=se_start-beta_start;
                            for(int k=0;k<numsnps;k++)
                            {
                                uint32_t ge_rowid=etrait->_rowid[beta_start+k];
                                if(ge_rowid==gdId[i])
                                {
                                    etrait->_val[beta_start+k]=-1*etrait->_val[beta_start+k];
                                    break; //the rowid in a probe is unique.
                                }
                            }
                        }
                    }
                    else
                    {
                        #pragma omp parallel for private(i)
                        for(int j=0;j<etrait->_include.size();j++)
                            if( etrait->_bxz[j][gdId[i]]+9 > 1e-6 )
                                etrait->_bxz[j][gdId[i]]=-1*etrait->_bxz[j][gdId[i]];
                    }
                }
            }
            else if(ea1 == a2 &&  ea2 == a1)
            {
                
                if( ea1 == ta1 && ea2 == ta2)
                {
                    bdata->_include.push_back(bdId[i]);
                    etrait->_esi_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
                else if(ea1 == ta2 && ea2 == ta1)
                {
                    bdata->_include.push_back(bdId[i]);
                    etrait->_esi_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    if(etrait->_val.size()>0)
                    {
                        /*
                        int count=0;
                        for(int j=0;j<etrait->_rowid.size();j++)
                        {
                            if(etrait->_rowid[j]==gdId[i])
                            {
                                count++;
                                if(count & 1) etrait->_val[j]=-etrait->_val[j];
                            }
                        }
                         */
                        #pragma omp parallel for private(i)
                        for(int j=0;j<etrait->_probNum;j++)
                        {
                            uint64_t beta_start=etrait->_cols[j<<1];
                            uint64_t se_start=etrait->_cols[1+(j<<1)];
                            uint64_t numsnps=se_start-beta_start;
                            for(int k=0;k<numsnps;k++)
                            {
                                uint32_t ge_rowid=etrait->_rowid[beta_start+k];
                                if(ge_rowid==gdId[i])
                                {
                                    etrait->_val[beta_start+k]=-1*etrait->_val[beta_start+k];
                                    break; //the rowid in a probe is unique.
                                }
                            }
                        }

                    }
                    else
                    {
                        #pragma omp parallel for private(i)
                        for(int j=0;j<etrait->_include.size();j++)
                            if( etrait->_bxz[j][gdId[i]]+9 > 1e-6 )
                                etrait->_bxz[j][gdId[i]]=-1*etrait->_bxz[j][gdId[i]];
                    }
                    string tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                    
                }
            }
            
            
        }
        
        logstr=atos(bdata->_include.size())+" SNPs are included after allele check. ";
        cout<<logstr<<endl;
        
        // only update _snp_name_map which would be used in read bed file.
        map<string, int> id_map_buf(bdata->_snp_name_map);
        for(int i=0; i<bdata->_include.size(); i++) id_map_buf.erase(bdata->_snp_name[bdata->_include[i]]);
        map<string, int>::iterator iter;
        for(iter=id_map_buf.begin(); iter!=id_map_buf.end(); iter++) bdata->_snp_name_map.erase(iter->first);
    }

    void update_gwas(gwasData* gdata){
        
        bool hasBP=false;
        if(gdata->snpBp.size()>0) hasBP=true;
        vector<int> snpBp;
        if(hasBP) snpBp.resize(gdata->_include.size());
        vector<string> snpName(gdata->_include.size());
        vector<string> allele_1(gdata->_include.size());
        vector<string> allele_2(gdata->_include.size());
        vector<double> freq(gdata->_include.size());
        vector<double> byz(gdata->_include.size());
        vector<double> seyz(gdata->_include.size());
        vector<double> pvalue(gdata->_include.size());
        vector<uint32_t> splSize(gdata->_include.size());
        
        gdata->snpNum=gdata->_include.size();
        for(int i=0;i<gdata->_include.size();i++ )
        {
            snpName[i]=gdata->snpName[gdata->_include[i]];
            allele_1[i]=gdata->allele_1[gdata->_include[i]];
            allele_2[i]=gdata->allele_2[gdata->_include[i]];
            freq[i]=gdata->freq[gdata->_include[i]];
            byz[i]=gdata->byz[gdata->_include[i]];
            seyz[i]=gdata->seyz[gdata->_include[i]];
            pvalue[i]=gdata->pvalue[gdata->_include[i]];
            splSize[i]=gdata->splSize[gdata->_include[i]];
            if(hasBP) snpBp[i]=gdata->snpBp[gdata->_include[i]];
        }
        
        gdata->allele_1.clear();
        gdata->allele_2.clear();
        gdata->freq.clear();
        gdata->byz.clear();
        gdata->seyz.clear();
        gdata->pvalue.clear();
        gdata->splSize.clear();
        
        gdata->snpName.swap(snpName);
        gdata->allele_1.swap(allele_1);
        gdata->allele_2.swap(allele_2);
        gdata->freq.swap(freq);
        gdata->byz.swap(byz);
        gdata->seyz.swap(seyz);
        gdata->pvalue.swap(pvalue);
        gdata->splSize.swap(splSize);
        
        if(hasBP) gdata->snpBp=snpBp;
        for(int i=0;i<gdata->snpNum;i++) gdata->_include[i]=i;
        
    }
    
    void make_XMat(bInfo* bdata,vector<uint32_t> &snpids, MatrixXd &X, bool minus_2p) {
        // Eigen is column-major by default. here row of X is individual, column of X is SNP.
        uint64_t snpNum=snpids.size();
        X.resize(bdata->_keep.size(),snpNum);
        #pragma omp parallel for
        for (int i = 0; i < snpNum ; i++)
        {
            uint32_t snpid=snpids[i];
            for (int j = 0; j < bdata->_keep.size() ; j++)
            {
                
                if (!bdata->_snp_1[bdata->_include[snpid]][bdata->_keep[j]] || bdata->_snp_2[bdata->_include[snpid]][bdata->_keep[j]])
                {
                    if (bdata->_allele1[bdata->_include[snpid]] == bdata->_ref_A[bdata->_include[snpid]]) X(j,i)= bdata->_snp_1[bdata->_include[snpid]][bdata->_keep[j]] + bdata->_snp_2[bdata->_include[snpid]][bdata->_keep[j]];
                    else X(j,i)= 2.0 - (bdata->_snp_1[bdata->_include[snpid]][bdata->_keep[j]] + bdata->_snp_2[bdata->_include[snpid]][bdata->_keep[j]]);
                } else X(j,i) = bdata->_mu[bdata->_include[snpid]];
                if (minus_2p) X(j,i) -= bdata->_mu[bdata->_include[snpid]];
                
            }
        }
    }
    void cor_calc(MatrixXd &LD, MatrixXd &X)
    {
        long size=X.cols();
        long n=X.rows();
        VectorXd tmpX(size);
        VectorXd tmpX2(size);
        #pragma omp parallel for
        for(int i=0;i<size;i++){
            tmpX[i]=X.col(i).sum();
            tmpX2[i]=X.col(i).dot(X.col(i));
        }
        LD.noalias()=X.transpose()*X;
        VectorXd tmpXX=(sqrt(tmpX2.array()*n-tmpX.array()*tmpX.array())).matrix();
        LD = (LD*n-tmpX*tmpX.transpose()).array()/ (tmpXX*tmpXX.transpose()).array();
        
    }
    void update_geIndx(bInfo* bdata, gwasData* gdata, eqtlInfo* esdata)
    {
        vector<int> tmpIdx1;
        vector<int> tmpIdx2;
        for (int i = 0; i < bdata->_include.size(); i++)
        {
            tmpIdx1.push_back(gdata->_include[bdata->_include[i]]);
            tmpIdx2.push_back(esdata->_esi_include[bdata->_include[i]]);
        }
        gdata->_include.clear();
        esdata->_esi_include.clear();
        gdata->_include = tmpIdx1;
        esdata->_esi_include = tmpIdx2;
    }
    void update_geIndx(bInfo* bdata, eqtlInfo* etrait, eqtlInfo* esdata)
    {
        vector<int> tmpIdx1;
        vector<int> tmpIdx2;
        for (int i = 0; i < bdata->_include.size(); i++)
        {
            tmpIdx1.push_back(etrait->_esi_include[bdata->_include[i]]);
            tmpIdx2.push_back(esdata->_esi_include[bdata->_include[i]]);
        }
        etrait->_esi_include.clear();
        esdata->_esi_include.clear();
        etrait->_esi_include = tmpIdx1;
        esdata->_esi_include = tmpIdx2;
    }


    void read_smaslist(vector<string> &smasNames, string eqtlsmaslstName)
    {
        ifstream smas(eqtlsmaslstName.c_str());
        if (!smas) throw ("Error: can not open the file [" + eqtlsmaslstName + "] to read.");
        cout << "Reading eQTL summary file names from [" + eqtlsmaslstName + "]." << endl;
        char buf[MAX_LINE_SIZE];
        map<string, int> probe_map;
        long mapsize=0;
        while (smas.getline(buf, MAX_LINE_SIZE))
        {
            string tmpStr;
            istringstream iss(buf);
            iss >> tmpStr;
            probe_map.insert(pair<string,int>(tmpStr,mapsize));
            if(mapsize<probe_map.size())
            {
                smasNames.push_back(tmpStr);
                mapsize=probe_map.size();
            }
            else
            {
                printf("WARNING: duplicate summary file name %s found and skipped.\n",tmpStr.c_str());
            }
            
        }
        cout << smasNames.size()<<" eQTL summary file names are included from [" + eqtlsmaslstName + "]." << endl;
        smas.close();
    }
  
    void ld_calc_o2m(VectorXd &ld_v,long targetid, MatrixXd &X)
    {
        long size=X.cols();
        long n=X.rows();
        
        VectorXd tmpX(size);
        VectorXd tmpX2(size);
        VectorXd tmpXY(size);
        #pragma omp parallel for
        for(int i=0;i<size;i++){
            tmpX[i]=X.col(i).sum();
            tmpX2[i]=X.col(i).dot(X.col(i));
            tmpXY[i]=X.col(targetid).dot(X.col(i));
        }
        float tmpY=X.col(targetid).sum();;
        float tmpY2=X.col(targetid).dot(X.col(targetid));
        ld_v=(tmpXY*n-tmpX*tmpY).array()/sqrt((tmpX2.array()*n-tmpX.array()*tmpX.array())*(tmpY2*n-tmpY*tmpY));
    }
    void init_smr_wk(SMRWK* smrwk)
    {
        smrwk->bxz.clear(),smrwk->sexz.clear(),smrwk->curId.clear(),smrwk->rs.clear(),smrwk->snpchrom.clear(),smrwk->byz.clear();
        smrwk->seyz.clear(),smrwk->pyz.clear(),smrwk->bpsnp.clear(),smrwk->allele1.clear(),smrwk->allele2.clear(),smrwk->freq.clear(),smrwk->zxz.clear();
    }
    long fill_smr_wk(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, const char* refSNP, int lowerbp,int upperbp)
    {
        int i=smrwk->cur_prbidx;
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
                        smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                        smrwk->byz.push_back(gdata->byz[j]);
                        smrwk->seyz.push_back(gdata->seyz[j]);
                        smrwk->pyz.push_back(gdata->pvalue[j]);
                        smrwk->curId.push_back(j);
                        smrwk->rs.push_back(esdata->_esi_rs[j]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                        if(refSNP!=NULL && esdata->_esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
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
                    smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                    smrwk->byz.push_back(gdata->byz[ge_rowid]);
                    smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                    smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                    smrwk->curId.push_back(ge_rowid);
                    smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                    smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                    if(refSNP!=NULL && esdata->_esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                    smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                    
                }
            }
        }
        return maxid;
    }
    long fill_smr_wk(bInfo* bdata,gwasData* gdata,eqtlInfo* esdata,SMRWK* smrwk, const char* refSNP,int cis_itvl,bool heidioffFlag)
    {
        int i=smrwk->cur_prbidx;
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
                        smrwk->zxz.push_back(esdata->_bxz[i][j]/esdata->_sexz[i][j]);
                        smrwk->byz.push_back(gdata->byz[j]);
                        smrwk->seyz.push_back(gdata->seyz[j]);
                        smrwk->pyz.push_back(gdata->pvalue[j]);
                        smrwk->curId.push_back(j);
                        smrwk->rs.push_back(esdata->_esi_rs[j]);
                        smrwk->snpchrom.push_back(esdata->_esi_chr[j]);
                        smrwk->allele1.push_back(esdata->_esi_allele1[j]);
                        smrwk->allele2.push_back(esdata->_esi_allele2[j]);
                        if(refSNP!=NULL && esdata->_esi_rs[j]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                        smrwk->bpsnp.push_back(esdata->_esi_bp[j]);
                        if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[j]] / 2);
                        else smrwk->freq.push_back(esdata->_esi_freq[j]);
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
                    smrwk->zxz.push_back(esdata->_val[beta_start+j]/esdata->_val[se_start+j]);
                    smrwk->byz.push_back(gdata->byz[ge_rowid]);
                    smrwk->seyz.push_back(gdata->seyz[ge_rowid]);
                    smrwk->pyz.push_back(gdata->pvalue[ge_rowid]);
                    smrwk->curId.push_back(ge_rowid); //save snp id of the raw datastruct
                    smrwk->rs.push_back(esdata->_esi_rs[ge_rowid]);
                    smrwk->snpchrom.push_back(esdata->_esi_chr[ge_rowid]);
                    smrwk->allele1.push_back(esdata->_esi_allele1[ge_rowid]);
                    smrwk->allele2.push_back(esdata->_esi_allele2[ge_rowid]);
                    if(refSNP!=NULL && esdata->_esi_rs[ge_rowid]==string(refSNP)) maxid=(smrwk->rs.size()-1);
                    smrwk->bpsnp.push_back(esdata->_esi_bp[ge_rowid]);
                    if(!heidioffFlag) smrwk->freq.push_back(bdata->_mu[bdata->_include[ge_rowid]] / 2);
                    else smrwk->freq.push_back(esdata->_esi_freq[ge_rowid]);
                }
            }
        }
        
        return maxid;
    }
    double heidi_test(bInfo* bdata,SMRWK* smrwk, long maxid,double ld_top, double threshold, int m_hetero,long &nsnp )
    {
        VectorXd ld_v;
        MatrixXd _X;
        vector<int> sn_ids;
        Map<VectorXd> ei_bxz(&smrwk->bxz[0],smrwk->bxz.size());
        Map<VectorXd> ei_sexz(&smrwk->sexz[0],smrwk->sexz.size());
        
        VectorXd zsxz=ei_bxz.array()/ei_sexz.array();
        
        
        make_XMat(bdata,smrwk->curId, _X);
        ld_calc_o2m(ld_v,maxid,_X);
        if(abs(ld_top-1)<1e-6) get_square_idxes(sn_ids,zsxz,threshold);
        else get_square_ldpruning_idxes(sn_ids,zsxz,threshold,ld_v, maxid,ld_top);
        if(sn_ids.size() < m_hetero) {
            
            printf("INFO: HEIDI test is skipped becasuse the number of SNPs (%ld) is less than a threshold (%d).\n", sn_ids.size(), m_hetero);
            return -9;
        }
        
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
            _byz[j]=smrwk->byz[sn_ids[j]];
            _seyz[j]=smrwk->seyz[sn_ids[j]];
            _bxz[j]=smrwk->bxz[sn_ids[j]];
            _sexz[j]=smrwk->sexz[sn_ids[j]];
            _zsxz[j]=zsxz(sn_ids[j]);
            _X_heidi.col(j)=_X.col(sn_ids[j]);
        }
        _X.resize(0,0);
        cor_calc(_LD_heidi, _X_heidi);
        
        _X_heidi.resize(0,0);
        
        nsnp = sn_ids.size();
        double pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, _LD_heidi, &nsnp);
        
        return pdev;
    }
    void update_snidx(SMRWK* smrwk,vector<int> &sn_ids,int max_snp_slct, string forwhat)
    {
        if(sn_ids.size()>max_snp_slct)
        {
            printf("Top %d SNPs are selected for %s...\n",max_snp_slct,forwhat.c_str());
            priority_queue< pair<double, int> > q;
            vector<int> sn_slct_ids;
            for(int i=0;i<sn_ids.size();i++) q.push((pair<double, int>(abs(smrwk->zxz[sn_ids[i]]), sn_ids[i])));
            for(int i=0;i<max_snp_slct;i++){
                sn_slct_ids.push_back(q.top().second);
                q.pop();
            }
            sn_ids.swap(sn_slct_ids);
        }
        
    }
    void extract_smrwk(SMRWK* smrwk,vector<int> &sn_ids,SMRWK* smrwk2)
    {
        smrwk2->cur_chr=smrwk->cur_chr;
        smrwk2->cur_prbidx=smrwk->cur_prbidx;
        for(int i=0;i<sn_ids.size();i++)
        {
            smrwk2->bxz.push_back(smrwk->bxz[sn_ids[i]]);
            smrwk2->sexz.push_back(smrwk->sexz[sn_ids[i]]);
            smrwk2->freq.push_back(smrwk->freq[sn_ids[i]]);
            smrwk2->byz.push_back(smrwk->byz[sn_ids[i]]);
            smrwk2->seyz.push_back(smrwk->seyz[sn_ids[i]]);
            smrwk2->pyz.push_back(smrwk->pyz[sn_ids[i]]);
            smrwk2->zxz.push_back(smrwk->zxz[sn_ids[i]]);
            smrwk2->curId.push_back(smrwk->curId[sn_ids[i]]);
            smrwk2->bpsnp.push_back(smrwk->bpsnp[sn_ids[i]]);
            smrwk2->snpchrom.push_back(smrwk->snpchrom[sn_ids[i]]);
            smrwk2->rs.push_back(smrwk->rs[sn_ids[i]]);
            smrwk2->allele1.push_back(smrwk->allele1[sn_ids[i]]);
            smrwk2->allele2.push_back(smrwk->allele2[sn_ids[i]]);
        }

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
    void update_smrwk_x(SMRWK* smrwk,vector<int> &sn_ids,MatrixXd &X)
    {
        vector<double> byz,seyz, bxz,sexz,zxz,freq,pyz;
        vector<uint32_t> curId;
        vector<int>  bpsnp, snpchrom;
        vector<string> rs,allele1, allele2;
        MatrixXd _X;
        byz.resize(sn_ids.size());
        seyz.resize(sn_ids.size());
        bxz.resize(sn_ids.size());
        sexz.resize(sn_ids.size());
        zxz.resize(sn_ids.size());
        freq.resize(sn_ids.size());
        pyz.resize(sn_ids.size());
        curId.resize(sn_ids.size());
        bpsnp.resize(sn_ids.size());
        snpchrom.resize(sn_ids.size());
        rs.resize(sn_ids.size());
        allele1.resize(sn_ids.size());
        allele2.resize(sn_ids.size());
        _X.resize(X.rows(), sn_ids.size());
        
        #pragma omp parallel for
        for(int j=0;j<sn_ids.size();j++)
        {
            byz[j]=smrwk->byz[sn_ids[j]];
            seyz[j]=smrwk->seyz[sn_ids[j]];
            bxz[j]=smrwk->bxz[sn_ids[j]];
            sexz[j]=smrwk->sexz[sn_ids[j]];
            zxz[j]=smrwk->zxz[sn_ids[j]];
            freq[j]=smrwk->freq[sn_ids[j]];
            pyz[j]=smrwk->pyz[sn_ids[j]];
            curId[j]=smrwk->curId[sn_ids[j]];
            bpsnp[j]=smrwk->bpsnp[sn_ids[j]];
            snpchrom[j]=smrwk->snpchrom[sn_ids[j]];
            rs[j]=smrwk->rs[sn_ids[j]];
            allele1[j]=smrwk->allele1[sn_ids[j]];
            allele2[j]=smrwk->allele2[sn_ids[j]];
            _X.col(j)=X.col(sn_ids[j]);
        }
        smrwk->bxz.swap(bxz);
        smrwk->sexz.swap(sexz);
        smrwk->freq.swap(freq);
        smrwk->byz.swap(byz);
        smrwk->seyz.swap(seyz);
        smrwk->pyz.swap(pyz);
        smrwk->zxz.swap(zxz);
        smrwk->curId.swap(curId);
        smrwk->bpsnp.swap(bpsnp);
        smrwk->snpchrom.swap(snpchrom);
        smrwk->rs.swap(rs);
        smrwk->allele1.swap(allele1);
        smrwk->allele2.swap(allele2);
        X=_X;
    }


    double heidi_test_new(bInfo* bdata,SMRWK* smrwk, long maxid,double ldr2_top, double threshold, int m_hetero,long &nsnp )
    {
        VectorXd ld_v;
        MatrixXd _X;
        vector<int> sn_ids;
        double pthres=pchisq(threshold,1);
        printf("Filtering SNPs (total %ld) with a p-value threshold for HEIDI test %e.\n",smrwk->zxz.size(), pthres);
        for(int i=0;i<smrwk->zxz.size();i++)
        {
            if(smrwk->zxz[i]*smrwk->zxz[i]-threshold>1e-6) sn_ids.push_back(i);
        }
        if(sn_ids.size() < m_hetero) {
            
            printf("INFO: HEIDI test is skipped becasuse the number of SNPs (%ld) is less than a threshold (%d).\n", sn_ids.size(), m_hetero);
            return -9;
        }
        printf("%ld SNPs left after filtration.\n",sn_ids.size());
        update_snidx(smrwk,sn_ids,500,"LD pruning");
        
        SMRWK smrwk_heidi;
        extract_smrwk(smrwk,sn_ids,&smrwk_heidi);
        long maxid_heidi=max_abs_id(smrwk_heidi.zxz);
        
        make_XMat(bdata,smrwk_heidi.curId, _X);
        printf("Removing SNPs with LD r-square > the threshold %f between top-SNP %s ...\n",ldr2_top,smrwk_heidi.rs[maxid_heidi].c_str());
        ld_calc_o2m(ld_v,maxid_heidi,_X);
        sn_ids.clear();
        if(abs(ldr2_top-1)>1e-6) {
            
            for(int i=0;i<smrwk_heidi.zxz.size();i++)
            {
                if(i!= maxid_heidi)
                {
                    if((ld_v(i)*ld_v(i)-ldr2_top)<1e-6) sn_ids.push_back(i);
                }
                else{
                    sn_ids.push_back(i);
                }
            }
        }
        printf("%ld SNPs are removed and %ld SNPs are left.\n",smrwk_heidi.zxz.size()-sn_ids.size(),sn_ids.size());
        if(sn_ids.size() < m_hetero) {
            printf("INFO: HEIDI test is skipped becasuse the number of SNPs (%ld) is less than a threshold (%d).\n", sn_ids.size(), m_hetero);
            return -9;
        }
        update_smrwk_x(&smrwk_heidi,sn_ids,_X);
        maxid_heidi=max_abs_id(smrwk_heidi.zxz);
        printf("Removing one of each pair of remaining SNPs with LD r-square > the threshold %f...\n",ldr2_top);
        int m = (int)smrwk_heidi.bxz.size();
        vector<int> rm_ID1;
        MatrixXd C;
        cor_calc(C, _X);
        double ld_top=sqrt(ldr2_top);
        if (ld_top < 1) rm_cor_sbat(C, ld_top, m, rm_ID1);
        printf("%ld SNPs are removed and %ld SNPs (including the top SNP %s) are left.\n",rm_ID1.size(),m-rm_ID1.size(),smrwk_heidi.rs[maxid_heidi].c_str());
        if(m-rm_ID1.size() < m_hetero) {
            
            printf("INFO: HEIDI test is skipped becasuse the number of SNPs (%ld) is less than a threshold (%d).\n", sn_ids.size(), m_hetero);
            return -9;
        }
        //Create new index
        sn_ids.clear();
        int qi=0;
        for (int i=0 ; i<m ; i++) {
            if (rm_ID1.size() == 0) sn_ids.push_back(i);
            else {
                if (rm_ID1[qi] == i) qi++; //Skip removed snp
                else sn_ids.push_back(i);
            }
        }
        
        update_snidx(&smrwk_heidi,sn_ids,20,"HEIDI test");
        if (sn_ids.size() < C.size()) { //Build new matrix
            MatrixXd D(sn_ids.size(),sn_ids.size());
            for (int i = 0 ; i < sn_ids.size() ; i++) {
                for (int j = 0 ; j < sn_ids.size() ; j++) {
                    D(i,j) = C(sn_ids[i],sn_ids[j]);
                }
            }
            C = D;
        }        
        
        VectorXd _byz,_seyz, _bxz,_sexz,_zsxz;
        _byz.resize(sn_ids.size());
        _seyz.resize(sn_ids.size());
        _bxz.resize(sn_ids.size());
        _sexz.resize(sn_ids.size());
        _zsxz.resize(sn_ids.size());
        #pragma omp parallel for
        for(int j=0;j<sn_ids.size();j++)
        {
            _byz[j]=smrwk_heidi.byz[sn_ids[j]];
            _seyz[j]=smrwk_heidi.seyz[sn_ids[j]];
            _bxz[j]=smrwk_heidi.bxz[sn_ids[j]];
            _sexz[j]=smrwk_heidi.sexz[sn_ids[j]];
            _zsxz[j]=smrwk_heidi.zxz[sn_ids[j]];
        }
       
        nsnp = sn_ids.size();
        double pdev=bxy_hetero3(_byz,  _bxz, _seyz, _sexz, _zsxz, C, &nsnp);
        
        return pdev;
    }
    void smr_heidi_func(vector<SMRRLT> &smrrlts, char* outFileName, bInfo* bdata,gwasData* gdata,eqtlInfo* esdata, int cis_itvl, bool heidioffFlag, const char* refSNP,double p_hetero,double ld_top,int m_hetero , double p_smr,double threshpsmrest)
    {
        
        uint64_t probNum = esdata->_probNum;
        double thresh_heidi= chi_val(1,p_hetero);
        VectorXd _byz,_seyz,_bxz,_sexz,_zsxz,ld_v,zsxz;
        MatrixXd _X,_LD,_LD_heidi,_X_heidi;
        
        FILE* smr=NULL;
        long write_count=0;
        string outstr="";
        string smrfile="";
        if(outFileName!=NULL){
            smrfile = string(outFileName)+".smr";
            smr = fopen(smrfile.c_str(), "w");
            if (!(smr)) {
                printf("Open error %s\n", smrfile.c_str());
                exit(1);
            }
            outstr="probeID\tProbeChr\tGene\tProbe_bp\ttopSNP\ttopSNP_chr\ttopSNP_bp\tA1\tA2\tFreq\tb_GWAS\tse_GWAS\tp_GWAS\tb_eQTL\tse_eQTL\tp_eQTL\tb_SMR\tse_SMR\tp_SMR\tp_HET\tnsnp_HET\n";
            if(fputs_checked(outstr.c_str(),smr))
            {
                printf("ERROR: in writing file %s .\n", smrfile.c_str());
                exit(EXIT_FAILURE);
            }
        } else {
            smrrlts.clear();
        }
        
        
        SMRWK smrwk;
        cout<<endl<<"Performing SMR and heterogeneity analysis..... "<<endl;
        float progr0=0.0 , progr1;
        progress_print(progr0);
        
        cis_itvl=cis_itvl*1000;
        for(int i=0;i<probNum;i++)
        {
            
            progr1=1.0*i/probNum;
            if(progr1-progr0-0.05>1e-6 || i+1==probNum)
            {
                if(i+1==probNum) progr1=1.0;
                progress_print(progr1);
                progr0=progr1;
            }
            
            int probebp=esdata->_epi_bp[i];
            int probechr=esdata->_epi_chr[i];
            string probename=esdata->_epi_prbID[i];
            string probegene=esdata->_epi_gene[i];
            char probeorien=esdata->_epi_orien[i];
            
            init_smr_wk(&smrwk);
            smrwk.cur_prbidx=i;
            smrwk.cur_chr=probechr;
            long maxid =fill_smr_wk(bdata, gdata, esdata, &smrwk, refSNP, cis_itvl,heidioffFlag);
            
            if(refSNP!=NULL && maxid==-9)
            {
                printf("WARNING: can't find target SNP %s in probe %s.\n",refSNP, probename.c_str());
                continue;
            }
            if (smrwk.bxz.size() == 0) {
                
                printf("WARNING: No SNP fetched in probe %s.\n", probename.c_str());
                continue;
            }
            
            Map<VectorXd> ei_bxz(&smrwk.bxz[0],smrwk.bxz.size());
            Map<VectorXd> ei_sexz(&smrwk.sexz[0],smrwk.sexz.size());
            
            zsxz=ei_bxz.array()/ei_sexz.array();
            if(refSNP==NULL) maxid=max_abs_id(zsxz);
            double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
            
            if(refSNP==NULL && pxz_val>p_smr){
                printf("WARNING: no SNP passed a p-value threshold %e for SMR analysis in probe %s.\n", p_smr, probename.c_str());
                continue;
            } else {
                printf("Dealing with probe %s...\n", probename.c_str());
            }
            double bxy_val = smrwk.byz[maxid] / smrwk.bxz[maxid];
            double sexy_val = sqrt((smrwk.seyz[maxid] * smrwk.seyz[maxid] * smrwk.bxz[maxid] * smrwk.bxz[maxid] + smrwk.sexz[maxid] * smrwk.sexz[maxid] * smrwk.byz[maxid] * smrwk.byz[maxid]) / (smrwk.bxz[maxid] * smrwk.bxz[maxid] * smrwk.bxz[maxid] * smrwk.bxz[maxid]));
            double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
            double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
            
            SMRRLT currlt;
            if(smr)
            {
                outstr = probename + '\t' + atos(probechr) + '\t' + probegene + '\t' + atos(probebp) + '\t' + smrwk.rs[maxid] + '\t' + atos(smrwk.snpchrom[maxid]) + '\t' + atos(smrwk.bpsnp[maxid]) + '\t' + smrwk.allele1[maxid] + '\t' + smrwk.allele2[maxid] + '\t' + atos(smrwk.freq[maxid]) + '\t';
                outstr += atos(smrwk.byz[maxid]) + '\t' + atos(smrwk.seyz[maxid]) + '\t' + dtos(smrwk.pyz[maxid]) + '\t';
                outstr += atos(smrwk.bxz[maxid]) + '\t' + atos(smrwk.sexz[maxid]) + '\t' + dtos(pxz_val) + '\t';
                outstr += atos(bxy_val) + '\t' + atos(sexy_val) + '\t' + dtos(pxy_val) + '\t';
                
            } else {
                currlt.ProbeID=probename;
                currlt.ProbeChr=probechr;
                currlt.Gene=probegene;
                currlt.Probe_bp=probebp;
                currlt.Orien=probeorien;
                currlt.SNP=smrwk.rs[maxid];
                currlt.SNP_Chr=smrwk.snpchrom[maxid];
                currlt.SNP_bp=smrwk.bpsnp[maxid];
                currlt.A1=smrwk.allele1[maxid];
                currlt.A2=smrwk.allele2[maxid];
                currlt.Freq=smrwk.freq[maxid];
                currlt.b_GWAS=smrwk.byz[maxid];
                currlt.se_GWAS=smrwk.seyz[maxid];
                currlt.p_GWAS=smrwk.pyz[maxid];
                currlt.b_eQTL=smrwk.bxz[maxid];
                currlt.se_eQTL=smrwk.sexz[maxid];
                currlt.p_eQTL=pxz_val;
                currlt.b_SMR=bxy_val;
                currlt.se_SMR=sexy_val;
                currlt.p_SMR=pxy_val;
            }
            
            if(heidioffFlag || pxy_val>threshpsmrest)
            {
                printf("INFO: HEIDI test for probe %s is skipped becasuse HEIDI test turns off or p_SMR does not pass a threshold %e.\n", probename.c_str(),threshpsmrest);
                if(smr)
                {
                    outstr+= string("NA") + '\t' + string("NA") + '\n';
                    if(fputs_checked(outstr.c_str(),smr))
                    {
                        printf("ERROR: in writing file %s .\n", smrfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    write_count++;
                } else {
                    currlt.p_HET=-9;
                    currlt.nsnp=-9;
                    smrrlts.push_back(currlt);
                }
                
            }
            else
            {
                
                /****Generate the test data***/
                /*
                 MatrixXd _X;
                 make_XMat(bdata,smrwk.curId, _X);
                 string xfilename =string(outFileName)+".mat";
                 FILE* xfile = fopen(xfilename.c_str(), "w");
                 if (!(xfile)) {
                 printf("Open error %s\n", xfilename.c_str());
                 exit(1);
                 }
                 string str=smrwk.rs[0];
                 for( int j=1;j<smrwk.rs.size();j++) str+='\t'+smrwk.rs[j];
                 str+='\n';
                 fputs_checked(str.c_str(), xfile);
                 for(int i=0;i<_X.rows();i++)
                 {
                 string str=atos(_X(i,0));
                 for( int j=1;j<_X.cols();j++) str+='\t'+atos(_X(i,j));
                 str+='\n';
                 fputs_checked(str.c_str(), xfile);
                 }
                 fclose(xfile);
                
                string gwasfname =string(outFileName)+".gwas";
                FILE* gfile = fopen(gwasfname.c_str(), "w");
                if (!(gfile)) {
                    printf("Open error %s\n", gwasfname.c_str());
                    exit(1);
                }
                for( int j=0;j<smrwk.rs.size();j++)
                {
                    string str=smrwk.rs[j]+'\t'+atos(smrwk.byz[j])+'\t'+atos(smrwk.seyz[j])+'\n';
                    fputs_checked(str.c_str(), gfile);
                }
                fclose(gfile);
                
                string eqtlfname =string(outFileName)+".eqtl";
                FILE* efile = fopen(eqtlfname.c_str(), "w");
                if (!(efile)) {
                    printf("Open error %s\n", eqtlfname.c_str());
                    exit(1);
                }
                for( int j=0;j<smrwk.rs.size();j++)
                {
                    string str=smrwk.rs[j]+'\t'+atos(smrwk.bxz[j])+'\t'+atos(smrwk.sexz[j])+'\n';
                    fputs_checked(str.c_str(), efile);
                }
                fclose(efile);
                 */
                /*****/

                long nsnp=-9;
                double pdev=-9;
                if(!heidioffFlag)  pdev= heidi_test_new(bdata,&smrwk, maxid, ld_top,  thresh_heidi,  m_hetero, nsnp );
                if(smr)
                {
                    outstr+= (pdev > 0 ? dtos(pdev) : "NA") + '\t' + (nsnp > 0 ? atos(nsnp) : "NA") + '\n';
                    if(fputs_checked(outstr.c_str(),smr))
                    {
                        printf("ERROR: in writing file %s .\n", smrfile.c_str());
                        exit(EXIT_FAILURE);
                    }
                    write_count++;
                } else {
                    currlt.p_HET=pdev;
                    currlt.nsnp=(int)nsnp;
                    smrrlts.push_back(currlt);
                }
            }
        }
        
        if(smr){
            printf("Results of %ld probes have been saved in file %s.\n",write_count,smrfile.c_str());
            fclose(smr);
        } else {
            printf("Results of %ld probes have been returned.\n",smrrlts.size());
        }
        
    }
    void smr(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero , char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, char* refSNP, bool heidioffFlag, int cis_itvl,bool plotflg, char* genelistName, int chr,int prbchr, const char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag,double threshpsmrest)
    {
        
        setNbThreads(thread_num);
        
        bInfo bdata;
        gwasData gdata;
        eqtlInfo esdata;
        bool heidiFlag=false;
        
        if(!heidioffFlag && bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(refSNP!=NULL) heidiFlag=true;
        read_gwas_data( &gdata, gwasFileName);
        read_esifile(&esdata, string(eqtlFileName)+".esi");
        esi_man(&esdata, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, false,  cis_itvl, prbname);
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
        epi_man(&esdata, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde); //no switch the place ahead of extract_prob()
        read_besdfile(&esdata, string(eqtlFileName)+".besd");
       if(esdata._rowid.empty() && esdata._bxz.empty())
       {
           printf("No data included under current condition.\n");
           exit(EXIT_FAILURE);
       }
       vector<SMRRLT> smrrlts;
       smr_heidi_func(smrrlts,  outFileName, &bdata,&gdata,&esdata,  cis_itvl,  heidioffFlag, refSNP,p_hetero,ld_top, m_hetero , p_smr, threshpsmrest);
    }
    
    // test all the meQTL or eQTL <=threshold.  toooo many tests and toooo slow
    void smr_trans(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf, char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero , char* indilst2remove, char* snplst2exclde, char* problst2exclde, double p_trans,char* refSNP, bool heidioffFlag,int trans_itvl,bool plotflg)
    {
        setNbThreads(thread_num);
        
        bInfo bdata;
        gwasData gdata;
        eqtlInfo esdata;
        double threshold= chi_val(1,p_hetero);
        bool heidiFlag=false;
        
        if(!heidioffFlag && bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(refSNP!=NULL) heidiFlag=true;
        read_gwas_data( &gdata, gwasFileName);
        read_esifile(&esdata, string(eqtlFileName)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&esdata, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde); //no switch the place ahead of extract_eqtl_snp()
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
        if(problstName != NULL) extract_prob(&esdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde); //no switch the place ahead of extract_prob()
        read_besdfile(&esdata, string(eqtlFileName)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included under current condition.\n");
            exit(EXIT_FAILURE);
        }
       
        int outCount = -1;
        unsigned int probNum = esdata._probNum;
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
        
        // for plot
        vector<uint32_t> plot_snpidx;
        vector<uint32_t> plot_probeidx;
        vector<double> plot_bxz;
        vector<double> plot_sexz;
        
        cout<<endl<<"Performing SMR and heterogeneity analysis..... "<<endl;
        float progr0=0.0 , progr1;
        progress_print(progr0);
        
        
        vector<double> bxz;
        vector<double> sexz;
        vector<uint32_t> curId;
        vector<string> eName;
        vector<int> echr;
        
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
        trans_itvl=trans_itvl*1000;
        for(int i=0;i<probNum;i++)
        {
            
            progr1=1.0*i/probNum;
            if(progr1-progr0-0.05>1e-6 || i+1==probNum)
            {
                if(i+1==probNum) progr1=1.0;
                progress_print(progr1);
                progr0=progr1;
            }
            
            //extract info from eqtl summary and gwas summary
            bxz.clear();
            sexz.clear();
            curId.clear(); // is the idxes of bfile._include not the values of
            eName.clear();
            echr.clear();
            byz.clear();
            seyz.clear();
            allele1.clear();
            allele2.clear();
            bpsnp.clear();
            freq.clear();
            long maxid =-9;
            int probebp=esdata._epi_bp[i];
            int probechr=esdata._epi_chr[i];
            //get top trans SNPs
            vector<int> topTransBP;
            vector<int> topTransRowid;
            vector<float> topbeta;
            if(esdata._rowid.empty())
            {
                for (int j = 0; j<bdata._include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
                {
                    if (abs(esdata._bxz[i][j] + 9) > 1e-6)
                    {
                        int snpbp=esdata._esi_bp[j];
                        int snpchr=esdata._esi_chr[j];
                        float beta=esdata._bxz[i][j];
                        float se=esdata._sexz[i][j];
                        if(snpchr!=probechr || (snpchr==probechr && ABS(probebp-snpbp)>5000000))
                        {
                            float zeqtl=beta/se;
                            float pval=pchisq(zeqtl*zeqtl, 1);
                            if(pval<=p_trans)
                            {
                                topTransBP.push_back(snpbp);
                                topTransRowid.push_back(j);
                            }
                            
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
                    float beta=esdata._val[beta_start+j];
                    float se=esdata._val[se_start+j];

                    if(snpchr!=probechr || (snpchr==probechr && ABS(probebp-snpbp)>5000000))
                    {
                        float zeqtl=beta/se;
                        float pval=pchisq(zeqtl*zeqtl, 1);
                        if(pval<=p_trans)
                        {
                            topTransBP.push_back(snpbp);
                            topTransRowid.push_back(j); // here j is not row id in dense format. it is the value id of current probe
                            topbeta.push_back(beta);
                            
                        }
                        
                    }
                    
                }
            }
           
            for(int j=0;j<topTransBP.size();j++)
            {
                int curbp=topTransBP[j]; // =esdata._esi_bp[esdata._rowid[esdata._cols[i<<1]+topTransRowid[j]]];
                if(esdata._rowid.empty())
                 {
                     
                 }
                 else{
                     uint64_t beta_start=esdata._cols[i<<1];
                     uint64_t se_start=esdata._cols[1+(i<<1)];
                     uint64_t numsnps=se_start-beta_start;
                     for(int k=0;k<numsnps;k++)
                     {
                         int ge_rowid=esdata._rowid[beta_start+k];
                         int snpbp=esdata._esi_bp[ge_rowid];
                         int snpchr=esdata._esi_chr[ge_rowid];
                         if(ABS(curbp-snpbp)<=trans_itvl)
                         {
                              bxz.push_back(esdata._val[beta_start+k]);
                              sexz.push_back(esdata._val[se_start+k]);
                              byz.push_back(gdata.byz[ge_rowid]);
                              seyz.push_back(gdata.seyz[ge_rowid]);
                              curId.push_back(ge_rowid);
                              eName.push_back(esdata._esi_rs[ge_rowid]);
                             echr.push_back(snpchr);
                             if(heidiFlag && esdata._esi_rs[ge_rowid]==string(refSNP)) maxid=(eName.size()-1);
                             else if(curbp==snpbp) maxid=(eName.size()-1);
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
                
                outCount++;
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
                rschr[outCount]=atos(echr[maxid]);
                rsbp[outCount]=itos(bpsnp[maxid]);
                if(!heidioffFlag){
                    rsfreq[outCount]=dtosf(freq[maxid]);
                }
                rsa1[outCount]=allele1[maxid];
                rsa2[outCount]=allele2[maxid];
                
                if(heidioffFlag)
                {
                    prb1[outCount]= string("NA");
                    nsnp_test1[outCount]= string("NA");

                }
                else
                {
                    make_XMat(&bdata,curId, _X);
                    ld_calc_o2m(ld_v,maxid,_X);
                    
                    sn_ids.clear(); //increase order
                    if(abs(ld_top-1)<1e-6) get_square_idxes(sn_ids,zsxz,threshold);
                    else get_square_ldpruning_idxes(sn_ids,zsxz,threshold,ld_v, maxid,ld_top);
                    
                    if(sn_ids.size() < m_hetero)
                    {
                        prb1[outCount]= string("NA");
                        nsnp_test1[outCount]= string("NA");
                        top_match1[outCount]= string("NA");
                        ldrsq[outCount]= string("NA");
                        continue;
                    }
                    
                    if(plotflg){
                        for(int j=0;j<sn_ids.size();j++) {
                            plot_probeidx.push_back(i);
                            plot_snpidx.push_back(curId[sn_ids[j]]);
                            plot_bxz.push_back(bxz[sn_ids[j]]);
                            plot_sexz.push_back(sexz[sn_ids[j]]);
                        }
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
                
            }
            
        }
        
        //genesn in R is just probidx here. can refer to probinfo to get probeid, probenm, chr, gene, bp
        if(plotflg && !heidioffFlag)
        {
            string ldfile = string(outFileName)+".engs";
            ofstream ldio(ldfile.c_str());
            if (!ldio) throw ("Error: can not open the file " + ldfile + " to save!");
            ldio << "Chr" <<'\t' << "ProbeID"  << '\t' << "Prob_bp" << '\t'<< "SNP"<< '\t'  << "SNP_bp"<< '\t'<< "A1"<< '\t'<< "A2"<< '\t'<<"b_GWAS"<<'\t'<<"se_GWAS"<<'\t'<<"b_eQTL"<<'\t'<<"se_eQTL"<<'\n';
            for (int i = 0;i <plot_snpidx.size(); i++) {
                ldio<<esdata._epi_chr[plot_probeidx[i]]<<'\t'<<esdata._epi_prbID[plot_probeidx[i]]<<'\t'<<esdata._epi_bp[plot_probeidx[i]]<<'\t'<<esdata._esi_rs[plot_snpidx[i]]<<'\t'<<esdata._esi_bp[plot_snpidx[i]]<<'\t'<<esdata._esi_allele1[plot_snpidx[i]]<<'\t'<<esdata._esi_allele2[plot_snpidx[i]]<<'\t'<<gdata.byz[plot_snpidx[i]]<<'\t'<<gdata.seyz[plot_snpidx[i]]<<'\t'<<plot_bxz[i]<<'\t'<<plot_sexz[i]<<'\n';
            }
            cout<<"informantion of "<<plot_snpidx.size()<<" probes for plot have been saved in the file [" + ldfile + "]."<<endl;
            ldio.close();
            
            getUnique(plot_snpidx);
            make_XMat(&bdata,plot_snpidx, _X);
            cor_calc(_LD, _X);
            
            string ldidfile = string(outFileName)+".ld.id";
            ldio.open(ldidfile.c_str());
            if (!ldio) throw ("Error: can not open the file " + ldidfile + " to save!");
            for (int i = 0;i <plot_snpidx.size(); i++)
                ldio<<esdata._esi_rs[plot_snpidx[i]]<<'\n';
            cout<<"rs id of "<<plot_snpidx.size()<<" unique SNPs for plot have been saved in the file [" + ldidfile + "]."<<endl;
            ldio.close();
            
            string ldinfofile = string(outFileName)+".ld.gz";
            gzFile gz_outfile = gzopen(ldinfofile.c_str(), "wb");
            string ldstrs="";
            for(int i=0;i<_LD.cols();i++)
                for(int j=0;j<i;j++) ldstrs+=atos(_LD(i,j))+"\n";
            
            if(gzputs(gz_outfile, ldstrs.c_str()) == -1) cout<<"error"<<endl;
            gzclose(gz_outfile);
            cout<< "Lower triangle of LD score matrix of "<<plot_snpidx.size()<<" SNPs have been saved by row-major in the file ["+ ldinfofile + "]."<<endl;
        }
        string smrfile = string(outFileName)+".smr";
        ofstream smr(smrfile.c_str());
        if (!smr) throw ("Error: can not open the file " + smrfile + " to save!");
        
            smr << "ProbeID" <<'\t'<< "ProbeChr" <<'\t' << "Gene"  << '\t' << "Prob_bp" << '\t'<< "SNP"<< '\t'<<"SNP_Chr"<<'\t' << "SNP_bp"<< '\t' << "A1"<< '\t'<< "A2"<< '\t'<<"Freq"<<'\t'<<"b_GWAS"<<'\t'<<"se_GWAS"<<'\t'<< "p_GWAS" << '\t'<<"b_eQTL"<<'\t'<<"se_eQTL"<<'\t'<< "p_eQTL" << '\t'<< "b_SMR" << '\t'<< "se_SMR"<< '\t' << "p_SMR" << "\t"<< "p_HET"<< "\t" << "nsnp" << '\n';
            
            for (int i = 0;i <=outCount; i++) {
                smr<<esdata._epi_prbID[out_probid[i]]<<'\t'<<esdata._epi_chr[out_probid[i]]<<'\t'<<esdata._epi_gene[out_probid[i]]<<'\t'<<esdata._epi_bp[out_probid[i]]<<'\t'<<rsid[i]<<'\t'<<rschr[i]<<'\t'<<rsbp[i]<<'\t'<<rsa1[i]<<'\t'<<rsa2[i]<<'\t'<<rsfreq[i]<<'\t'<<bgwas[i]<<'\t'<<segwas[i]<<'\t'<<pgwas[i]<<'\t'<<beqtl[i]<<'\t'<<seeqtl[i]<<'\t'<<peqtl[i]<<'\t'<<bxy[i]<<'\t'<<sexy[i]<<'\t'<<pxy[i]<<'\t'<<prb1[i]<<'\t'<<nsnp_test1[i]<<'\n';
            }
            cout<<"SMR and heterogeneity analysis for trans regions finished.\nSMR and heterogeneity analysis results of "<<outCount+1<<" probes have been saved in the file [" + smrfile + "]."<<endl;
       
        smr.close();
        
        free_gwas_data( &gdata);
        

    }
      
    void smr_trans_wholeInOne(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf, char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero , char* indilst2remove, char* snplst2exclde, char* problst2exclde, double p_trans,char* refSNP, bool heidioffFlag,int trans_itvl,bool plotflg)
    {
        setNbThreads(thread_num);
        
        bInfo bdata;
        gwasData gdata;
        eqtlInfo esdata;
        double threshold= chi_val(1,p_hetero);
        bool heidiFlag=false;
        
        if(!heidioffFlag && bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(refSNP!=NULL) heidiFlag=true;
        read_gwas_data( &gdata, gwasFileName);
        read_esifile(&esdata, string(eqtlFileName)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&esdata, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde); //no switch the place ahead of extract_eqtl_snp()
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
        if(problstName != NULL) extract_prob(&esdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde); //no switch the place ahead of extract_prob()
       
        read_besdfile(&esdata, string(eqtlFileName)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included under current condition.\n");
            exit(EXIT_FAILURE);
        }
        
        int outCount = -1;
        unsigned int probNum = esdata._probNum;
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
        
        // for plot
        vector<uint32_t> plot_snpidx;
        vector<uint32_t> plot_probeidx;
        vector<double> plot_bxz;
        vector<double> plot_sexz;
        
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
        vector<int> bpsnp;
        vector<double> freq;
        
        vector<double> byz;
        vector<double> seyz;
        VectorXd zsxz;
        
        vector<int> sn_ids;
        
        VectorXi _bpsnp;
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
        trans_itvl=trans_itvl*1000;
        for(int i=0;i<probNum;i++)
        {
            
            progr1=1.0*i/probNum;
            if(progr1-progr0-0.05>1e-6 || i+1==probNum)
            {
                if(i+1==probNum) progr1=1.0;
                progress_print(progr1);
                progr0=progr1;
            }
            
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
            //get top trans SNPs
            vector<float> topbeta;
            if(esdata._rowid.empty())
            {
                for (int j = 0; j<bdata._include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
                {
                    if (abs(esdata._bxz[i][j] + 9) > 1e-6)
                    {
                        int snpbp=esdata._esi_bp[j];
                        int snpchr=esdata._esi_chr[j];
                        float beta=esdata._bxz[i][j];
                        float se=esdata._sexz[i][j];
                        if(snpchr!=probechr || (snpchr==probechr && ABS(probebp-snpbp)>5000000))
                        {
                            float zeqtl=beta/se;
                            float pval=pchisq(zeqtl*zeqtl, 1);
                            if(pval<=p_trans)
                            {
                                bxz.push_back(esdata._bxz[i][j]);
                                sexz.push_back(esdata._sexz[i][j]);
                                byz.push_back(gdata.byz[j]);
                                seyz.push_back(gdata.seyz[j]);
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
                
            }
            else
            {
                uint64_t beta_start=esdata._cols[i<<1];
                uint64_t se_start=esdata._cols[1+(i<<1)];
                uint64_t numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=esdata._rowid[beta_start+j];
                    int snpbp=esdata._esi_bp[ge_rowid];
                    int snpchr=esdata._esi_chr[ge_rowid];
                    float beta=esdata._val[beta_start+j];
                    float se=esdata._val[se_start+j];
                    
                    if(snpchr!=probechr || (snpchr==probechr && ABS(probebp-snpbp)>5000000))
                    {
                        float zeqtl=beta/se;
                        float pval=pchisq(zeqtl*zeqtl, 1);
                        if(pval<=p_trans)
                        {
                            bxz.push_back(esdata._val[beta_start+j]);
                            sexz.push_back(esdata._val[se_start+j]);
                            byz.push_back(gdata.byz[ge_rowid]);
                            seyz.push_back(gdata.seyz[ge_rowid]);
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
            }
            
       
                
            if(heidiFlag && maxid==-9) continue; //heidi SNP is not in selected SNPs
            if (bxz.size() == 0) continue;
                
            Map<VectorXd> ei_bxz(&bxz[0],bxz.size());
            Map<VectorXd> ei_sexz(&sexz[0],sexz.size());
            zsxz=ei_bxz.array()/ei_sexz.array();
            if(!heidiFlag) maxid=max_abs_id(zsxz);
            double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
                
                double bxy_val = byz[maxid] / bxz[maxid];
                double sexy_val = sqrt((seyz[maxid] * seyz[maxid] * bxz[maxid] * bxz[maxid] + sexz[maxid] * sexz[maxid] * byz[maxid] * byz[maxid]) / (bxz[maxid] * bxz[maxid] * bxz[maxid] * bxz[maxid]));
                double chisqxy = bxy_val*bxy_val / (sexy_val*sexy_val);
                
                double pxy_val = pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
                
                double chisqyz = byz[maxid] / seyz[maxid];
                double pyz_val = pchisq(chisqyz*chisqyz, 1);
                
                outCount++;
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
                if(!heidioffFlag){
                    rsfreq[outCount]=dtosf(freq[maxid]);
                }
                rsa1[outCount]=allele1[maxid];
                rsa2[outCount]=allele2[maxid];
                
                if(heidioffFlag)
                {
                    prb1[outCount]= string("NA");
                    nsnp_test1[outCount]= string("NA");
                    
                }
                else
                {
                    //distance vector
                    VectorXd vdistance(curId.size());
                    for(int j=0;j<curId.size();j++)
                        if(ABS(bpsnp[maxid]-bpsnp[j])>1e7) vdistance(j)=0;
                        else vdistance(j)=1;
                    //                   
                    
                    make_XMat(&bdata,curId, _X);
                    ld_calc_o2m(ld_v,maxid,_X);
                   
                    ld_v=ld_v.array()*vdistance.array();
                    
                    sn_ids.clear(); //increase order
                    if(abs(ld_top-1)<1e-6) get_square_idxes(sn_ids,zsxz,threshold);
                    else get_square_ldpruning_idxes(sn_ids,zsxz,threshold,ld_v, maxid,ld_top);
                    
                    if(sn_ids.size() < m_hetero)
                    {
                        prb1[outCount]= string("NA");
                        nsnp_test1[outCount]= string("NA");
                        top_match1[outCount]= string("NA");
                        ldrsq[outCount]= string("NA");
                        continue;
                    }
                    
                    if(plotflg){
                        for(int j=0;j<sn_ids.size();j++) {
                            plot_probeidx.push_back(i);
                            plot_snpidx.push_back(curId[sn_ids[j]]);
                            plot_bxz.push_back(bxz[sn_ids[j]]);
                            plot_sexz.push_back(sexz[sn_ids[j]]);
                        }
                    }
                    _bpsnp.resize(sn_ids.size());
                    _byz.resize(sn_ids.size());
                    _seyz.resize(sn_ids.size());
                    _bxz.resize(sn_ids.size());
                    _sexz.resize(sn_ids.size());
                    _zsxz.resize(sn_ids.size());
                    _X_heidi.resize(_X.rows(), sn_ids.size());
                    
                    #pragma omp parallel for
                    for(int j=0;j<sn_ids.size();j++)
                    {
                        _bpsnp[j]=bpsnp[sn_ids[j]];
                        _byz[j]=byz[sn_ids[j]];
                        _seyz[j]=seyz[sn_ids[j]];
                        _bxz[j]=bxz[sn_ids[j]];
                        _sexz[j]=sexz[sn_ids[j]];
                        _zsxz[j]=zsxz[sn_ids[j]];
                        _X_heidi.col(j)=_X.col(sn_ids[j]);
                    }
                    
                    _X.resize(sn_ids.size(),sn_ids.size()); //as MatrixXd mdist;
                    for(int j=0;j<sn_ids.size();j++)
                        for(int k=0;k<=j;k++)
                            if(ABS(_bpsnp[j]-_bpsnp[k])>1e7) _X(j,k)=_X(k,j)=0;
                            else _X(j,k)=_X(k,j)=1;
                  
                    
                    cor_calc(_LD_heidi, _X_heidi);
               
                    _LD_heidi=_LD_heidi.array()*_X.array();
                    
                    _X_heidi.resize(0,0);
                    _X.resize(0,0); // free space
                    
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
                
            }

        
        //genesn in R is just probidx here. can refer to probinfo to get probeid, probenm, chr, gene, bp
        if(plotflg && !heidioffFlag)
        {
            string ldfile = string(outFileName)+".engs";
            ofstream ldio(ldfile.c_str());
            if (!ldio) throw ("Error: can not open the file " + ldfile + " to save!");
            ldio << "Chr" <<'\t' << "ProbeID"  << '\t' << "Prob_bp" << '\t'<< "SNP"<< '\t'  << "SNP_bp"<< '\t'<< "A1"<< '\t'<< "A2"<< '\t'<<"b_GWAS"<<'\t'<<"se_GWAS"<<'\t'<<"b_eQTL"<<'\t'<<"se_eQTL"<<'\n';
            for (int i = 0;i <plot_snpidx.size(); i++) {
                ldio<<esdata._epi_chr[plot_probeidx[i]]<<'\t'<<esdata._epi_prbID[plot_probeidx[i]]<<'\t'<<esdata._epi_bp[plot_probeidx[i]]<<'\t'<<esdata._esi_rs[plot_snpidx[i]]<<'\t'<<esdata._esi_bp[plot_snpidx[i]]<<'\t'<<esdata._esi_allele1[plot_snpidx[i]]<<'\t'<<esdata._esi_allele2[plot_snpidx[i]]<<'\t'<<gdata.byz[plot_snpidx[i]]<<'\t'<<gdata.seyz[plot_snpidx[i]]<<'\t'<<plot_bxz[i]<<'\t'<<plot_sexz[i]<<'\n';
            }
            cout<<"informantion of "<<plot_snpidx.size()<<" probes for plot have been saved in the file [" + ldfile + "]."<<endl;
            ldio.close();
            
            getUnique(plot_snpidx);
            make_XMat(&bdata,plot_snpidx, _X);
            cor_calc(_LD, _X);
            
            string ldidfile = string(outFileName)+".ld.id";
            ldio.open(ldidfile.c_str());
            if (!ldio) throw ("Error: can not open the file " + ldidfile + " to save!");
            for (int i = 0;i <plot_snpidx.size(); i++)
                ldio<<esdata._esi_rs[plot_snpidx[i]]<<'\n';
            cout<<"rs id of "<<plot_snpidx.size()<<" unique SNPs for plot have been saved in the file [" + ldidfile + "]."<<endl;
            ldio.close();
            
            string ldinfofile = string(outFileName)+".ld.gz";
            gzFile gz_outfile = gzopen(ldinfofile.c_str(), "wb");
            string ldstrs="";
            for(int i=0;i<_LD.cols();i++)
                for(int j=0;j<i;j++) ldstrs+=atos(_LD(i,j))+"\n";
            
            if(gzputs(gz_outfile, ldstrs.c_str()) == -1) cout<<"error"<<endl;
            gzclose(gz_outfile);
            cout<< "Lower triangle of LD score matrix of "<<plot_snpidx.size()<<" SNPs have been saved by row-major in the file ["+ ldinfofile + "]."<<endl;
        }
        string smrfile = string(outFileName)+".smr";
        ofstream smr(smrfile.c_str());
        if (!smr) throw ("Error: can not open the file " + smrfile + " to save!");
        
        smr << "ProbeID" <<'\t'<< "ProbeChr" <<'\t' << "Gene"  << '\t' << "Prob_bp" << '\t'<< "SNP"<< '\t'<<"SNP_Chr"<<'\t' << "SNP_bp"<< '\t' << "A1"<< '\t'<< "A2"<< '\t'<<"Freq"<<'\t'<<"b_GWAS"<<'\t'<<"se_GWAS"<<'\t'<< "p_GWAS" << '\t'<<"b_eQTL"<<'\t'<<"se_eQTL"<<'\t'<< "p_eQTL" << '\t'<< "b_SMR" << '\t'<< "se_SMR"<< '\t' << "p_SMR" << "\t"<< "p_HET"<< "\t" << "nsnp" << '\n';
        
        for (int i = 0;i <=outCount; i++) {
            smr<<esdata._epi_prbID[out_probid[i]]<<'\t'<<esdata._epi_chr[out_probid[i]]<<'\t'<<esdata._epi_gene[out_probid[i]]<<'\t'<<esdata._epi_bp[out_probid[i]]<<'\t'<<rsid[i]<<'\t'<<rschr[i]<<'\t'<<rsbp[i]<<'\t'<<rsa1[i]<<'\t'<<rsa2[i]<<'\t'<<rsfreq[i]<<'\t'<<bgwas[i]<<'\t'<<segwas[i]<<'\t'<<pgwas[i]<<'\t'<<beqtl[i]<<'\t'<<seeqtl[i]<<'\t'<<peqtl[i]<<'\t'<<bxy[i]<<'\t'<<sexy[i]<<'\t'<<pxy[i]<<'\t'<<prb1[i]<<'\t'<<nsnp_test1[i]<<'\n';
        }
        cout<<"SMR and heterogeneity analysis for trans regions finished.\nSMR and heterogeneity analysis results of "<<outCount+1<<" probes have been saved in the file [" + smrfile + "]."<<endl;
        
        smr.close();
        
        free_gwas_data( &gdata);
        
        
    }
    
   void make_esd_file(char* outFileName, char* eqtlFileName, char* snplstName,char* problstName,bool bFlag,bool make_besd_flag, char* snplst2exclde, char* problst2exclde, int cis_itvl,char* genelistName, int chr,int prbchr, char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag)
    {
        
        eqtlInfo eqtlinfo;
        
        cout<<endl<<"Reading eQTL summary data..."<<endl;
        if(eqtlFileName != NULL)
        {
            
            read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
            epi_man(&eqtlinfo, problstName, genelistName,  chr, prbchr,  prbname,  fromprbname,  toprbname, prbWind, fromprbkb,  toprbkb, prbwindFlag,  genename);
            if(problst2exclde != NULL) exclude_prob(&eqtlinfo, problst2exclde);
            
            read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
            esi_man(&eqtlinfo, snplstName, chr, snpchr,  snprs,  fromsnprs,  tosnprs, snpWind, fromsnpkb,  tosnpkb, snpwindFlag, cis_flag,  cis_itvl, prbname);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&eqtlinfo, snplst2exclde);
            
           read_besdfile(&eqtlinfo, string(eqtlFileName)+".besd");
            if(eqtlinfo._rowid.empty() && eqtlinfo._bxz.empty())
            {
                printf("No data included under current condition.\n");
                exit(EXIT_FAILURE);
            }
            
        }
        else throw ("Error: please input the eQTL summary information for the eQTL data files by the option --beqtl-summary.");
        
      
                filter_probe_null(&eqtlinfo); // at the same time, reset the vector _include
                filter_snp_null(&eqtlinfo);
                cout<<"\nsaving eQTL data..."<<endl;
                string esdfile = string(outFileName)+".esi";
                ofstream smr(esdfile.c_str());
                if (!smr) throw ("Error: can not open the esi file " + esdfile + " to save!");
                for (int i = 0;i <eqtlinfo._snpNum; i++) {
                    smr<<eqtlinfo._esi_chr[i]<<'\t'<<eqtlinfo._esi_rs[i]<<'\t'<<eqtlinfo._esi_gd[i]<<'\t'<<eqtlinfo._esi_bp[i]<<'\t'<<eqtlinfo._esi_allele1[i]<<'\t'<<eqtlinfo._esi_allele2[i]<<'\t'<<eqtlinfo._esi_freq[i]<<'\n';
                }
                smr.close();
                cout<<eqtlinfo._snpNum<<" SNPs have been saved in the file [" + esdfile + "]."<<endl;
                
                esdfile = string(outFileName)+".epi";
                smr.open(esdfile.c_str());
                if (!smr) throw ("Error: can not open the epi file " + esdfile + " to save!");
                for (int i = 0;i <eqtlinfo._include.size(); i++) {
                    smr<<eqtlinfo._epi_chr[eqtlinfo._include[i]]<<'\t'<<eqtlinfo._epi_prbID[eqtlinfo._include[i]]<<'\t'<<eqtlinfo._epi_gd[eqtlinfo._include[i]]<<'\t'<<eqtlinfo._epi_bp[eqtlinfo._include[i]]<<'\t'<<eqtlinfo._epi_gene[eqtlinfo._include[i]]<<'\t'<<eqtlinfo._epi_orien[eqtlinfo._include[i]]<<'\n';
                }
                smr.close();
                cout<<eqtlinfo._include.size()<<" probes have been saved in the file [" + esdfile + "]."<<endl;
            
                if(make_besd_flag)
                {
                    esdfile = string(outFileName)+".besd";
                    FILE * smr;
                    smr = fopen (esdfile.c_str(), "wb");
                    if(eqtlinfo._valNum==0)
                    {
                        uint64_t bsize=(eqtlinfo._include.size()*eqtlinfo._snpNum<<1)+1;
                        float* buffer=(float*)malloc (sizeof(float)*bsize);
                        memset(buffer,0,sizeof(float)*bsize);
                        float* ptr=buffer;
                        *ptr++=0.0;
                        uint64_t pro_num=eqtlinfo._include.size();
                        uint64_t snp_num=eqtlinfo._snpNum;
                        for(int i=0;i<pro_num;i++)
                        {
                            memcpy(ptr+(i<<1)*snp_num,&eqtlinfo._bxz[eqtlinfo._include[i]][0],sizeof(float)*snp_num);
                            memcpy(ptr+((i<<1)+1)*snp_num,&eqtlinfo._sexz[eqtlinfo._include[i]][0],sizeof(float)*snp_num);
                        }
                        fwrite (buffer,sizeof(float), bsize, smr);
                        free(buffer);
                    }
                    else
                    {
                        uint64_t colSize=sizeof(uint64_t)*((eqtlinfo._include.size()<<1)+1);
                        uint64_t rowSize=sizeof(uint32_t)*eqtlinfo._valNum;
                        uint64_t valSize=sizeof(float)*eqtlinfo._valNum;
                        uint64_t valNum=eqtlinfo._valNum;
                        uint64_t bufsize=sizeof(float)+sizeof(uint64_t)+colSize+rowSize+valSize;
                        
                        char* buffer=(char*)malloc (sizeof(char)*bufsize);
                        memset(buffer,0,sizeof(char)*bufsize);
                        float ftype=SPARSE_FILE_TYPE_3;
                        memcpy(buffer,&ftype,sizeof(float));
                        char* wptr=buffer+sizeof(float);
                        memcpy(wptr,&valNum,sizeof(uint64_t));
                        wptr+=sizeof(uint64_t);
                        uint64_t* uptr=(uint64_t*)wptr; *uptr++=0;
                        for(int i=0;i<eqtlinfo._include.size();i++)
                        {
                            *uptr++=eqtlinfo._cols[(eqtlinfo._include[i]<<1)+1];
                            *uptr++=eqtlinfo._cols[eqtlinfo._include[i]+1<<1];
                        }
                        wptr+=colSize;
                        memcpy(wptr,&eqtlinfo._rowid[0],rowSize);
                        wptr+=rowSize;
                        memcpy(wptr,&eqtlinfo._val[0],valSize);
                        fwrite (buffer,sizeof(char), bufsize, smr);
                        free(buffer);
                    }
                    fclose (smr);
                    
                    cout<<"Beta values and SE values for "<<eqtlinfo._include.size()<<" Probes and "<<eqtlinfo._snpNum<<" SNPs have been saved in the binary file [" + esdfile + "]." <<endl;
                }
    
    }
    
 
    int get_besd_format(string besdfName)
    {
        string besdfile = besdfName + ".besd";
        char buf[8];
        ifstream besd(besdfile.c_str(), ios::in|ios::binary);
        if(!besd)
        {
            fprintf (stderr, "%s: Couldn't open file %s\n",
                     besdfile.c_str(), strerror (errno));
            exit (EXIT_FAILURE);
        }
        
        besd.read(buf, 4);
        besd.close();
        
        float* flag=(float *)buf;
        return  (int)*flag;

    }
  

    //sort in ascend order
    int comp(const void *a,const void *b){ return (((*(probeinfolst *)a).probechr>(*(probeinfolst *)b).probechr) || ( ((*(probeinfolst *)a).probechr ==(*(probeinfolst *)b).probechr) && ((*(probeinfolst *)a).bp > (*(probeinfolst *)b).bp) ))?1:-1; }
    int comp2(const void *a,const void *b){ return (((*(probeinfolst2 *)a).probechr>(*(probeinfolst2 *)b).probechr) || ( ((*(probeinfolst2 *)a).probechr ==(*(probeinfolst2 *)b).probechr) && ((*(probeinfolst2 *)a).bp > (*(probeinfolst2 *)b).bp) ))?1:-1; }

   int comp_esi(const void *a,const void *b){ return (((*(snpinfolst *)a).snpchr >(*(snpinfolst *)b).snpchr) || ( ((*(snpinfolst *)a).snpchr ==(*(snpinfolst *)b).snpchr) && ((*(snpinfolst *)a).bp > (*(snpinfolst *)b).bp) ))?1:-1; }
   int comp_estn(const void *a,const void *b){ return ((*(snpinfolst *)a).estn >(*(snpinfolst *)b).estn)?1:-1; }
    
     void smr_e2e(char* outFileName, char* bFileName,char* eqtlFileName, char* eqtlFileName2, double maf,char* indilstName, char* snplstName,char* problstName, char* oproblstName,char* eproblstName,bool bFlag,double p_hetero,double ld_top,int m_hetero, char* indilst2remove, char* snplst2exclde, char* problst2exclde, char* oproblst2exclde,char* eproblst2exclde,double p_smr,char* refSNP, bool heidioffFlag,int cis_itvl,char* traitlstName,int op_wind, char* oprobe, char* eprobe, char* oprobe2rm, char* eprobe2rm, double prepsmrthres)
    {
        setNbThreads(thread_num);
        string logstr;
        if(oproblstName!=NULL && oprobe!=NULL)
        {
            logstr="WARNING: --extract-single-outcome-probe is not surpposed to use together with --extract-outcome-probe. --extract-single-outcome-probe will be disabled.\n";
            oprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblstName!=NULL && eprobe!=NULL)
        {
            logstr="WARNING: --extract-single-exposure-probe is not surpposed to use together with --extract-exposure-probe. --extract-single-exposure-probe will be disabled.\n";
            eprobe=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(oproblst2exclde!=NULL && oprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-outcome-probe is not surpposed to use together with --exclude-outcome-probe. --exclude-single-outcome-probe will be disabled.\n";
            oprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        if(eproblst2exclde!=NULL && eprobe2rm!=NULL)
        {
            logstr="WARNING: --exclude-single-exposure-probe is not surpposed to use together with --exclude-exposure-probe. --exclude-single-exposure-probe will be disabled.\n";
            eprobe2rm=NULL;
            fputs(logstr.c_str(), stdout);
        }
        
        eqtlInfo etrait;
        eqtlInfo esdata;
        bInfo bdata;
        double threshold= chi_val(1,p_hetero);
        double pre_psmr_thresh = chi_val(1,prepsmrthres);
        bool heidiFlag=false;
        cis_itvl=cis_itvl*1000;
        
        if(!heidioffFlag && bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(refSNP!=NULL) heidiFlag=true;
        if(problstName != NULL) cout<<"WARNING: --extract-probe here presumes the probe list should contain both probes of exposure dataset and probes of outcome dataset.\n If you want to only extract probes from one dataset please include these probles in the file and all the probes of the other dataset as well.\n"<<endl;
        read_esifile(&etrait, string(eqtlFileName)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&etrait, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&etrait, snplst2exclde);
        read_epifile(&etrait, string(eqtlFileName)+".epi");
        if(problstName != NULL) extract_prob(&etrait, problstName);
        if(problst2exclde != NULL) exclude_prob(&etrait, problst2exclde);
        if(oproblstName != NULL ) extract_prob(&etrait, oproblstName);
        else if(oprobe != NULL) extract_eqtl_single_probe(&etrait, oprobe);
        if(oproblst2exclde != NULL) exclude_prob(&etrait, oproblst2exclde);
        else if(oprobe2rm != NULL) exclude_eqtl_single_probe(&etrait, oprobe2rm);
        
        read_besdfile(&etrait, string(eqtlFileName)+".besd");
        if(etrait._rowid.empty() && etrait._bxz.empty())
        {
            printf("No data included under current condition.\n");
            exit(EXIT_FAILURE);
        }

        
        read_esifile(&esdata, string(eqtlFileName2)+".esi");
        //if (snplstName != NULL) extract_eqtl_snp(&esdata, snplstName); //no need here, already extracted in etrait
        //if(snplst2exclde != NULL) exclude_eqtl_snp(&esdata, snplst2exclde);
        if(!heidioffFlag)
        {
            read_famfile(&bdata, string(bFileName)+".fam");
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            read_bimfile(&bdata, string(bFileName)+".bim");
            //if(snplstName != NULL) extract_snp(&bdata, snplstName);
            //if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
            allele_check(&bdata, &etrait, &esdata);
            // if no snp left after check
            read_bedfile(&bdata, string(bFileName)+".bed");
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if (maf > 0)
            {
                filter_snp_maf(&bdata, maf);
                update_geIndx(&bdata, &etrait, &esdata);
            }
            
        }else
        {
            allele_check(&etrait, &esdata);
        }
       
        //the etrait is not updated, so from now on _esi_include should be used always.
        cout<<"Reading eQTL summary data..."<<endl;
        read_epifile(&esdata, string(eqtlFileName2)+".epi");
        if(problstName != NULL) extract_prob(&esdata, problstName);
        if(problst2exclde != NULL) exclude_prob(&esdata, problst2exclde);
        if(eproblstName != NULL ) extract_prob(&esdata, eproblstName);
        else if(eprobe != NULL) extract_eqtl_single_probe(&esdata, eprobe);
        if(eproblst2exclde != NULL) exclude_prob(&esdata, eproblst2exclde);       
        else if(eprobe2rm != NULL) exclude_eqtl_single_probe(&esdata, eprobe2rm);
        
        read_besdfile(&esdata, string(eqtlFileName2)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included under current condition.\n");
            exit(EXIT_FAILURE);
        }
        

        
        /********************/
        /*
        for(int z=0;z<esdata._snpNum ;z++)
        {
            string esrs=esdata._esi_rs[z];
            string etrs=etrait._esi_rs[etrait._esi_include[z]];
            string brs=bdata._snp_name[bdata._include[z]];
            bool flg=(esrs==etrs) && (esrs==brs) && (etrs==brs);
            if(!flg) cout<<"err"<<endl;
        }
         */
        /********************/
         /********************/
        /*
        string trs="rs393270"; //"rs413592";
        string prnm="cg05202042";
        long pos=find(esdata._esi_rs.begin(), esdata._esi_rs.end(), trs)-esdata._esi_rs.begin();
        long pos2=find(esdata._epi_prbID.begin(), esdata._epi_prbID.end(), prnm)-esdata._epi_prbID.begin();
        uint64_t beta_start=esdata._cols[pos2<<1];
        uint64_t se_start=esdata._cols[1+(pos2<<1)];
        uint64_t numsnps=se_start-beta_start;
        for(int j=0;j<numsnps;j++)
        {
            int ge_rowid=esdata._rowid[beta_start+j];
            if(ge_rowid==pos) cout<<esdata._val[beta_start+j]<<":"<<esdata._val[se_start+j]<<endl;
        }
         */
         /********************/
        
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
            
            vector<string> allele1, allele2;
            vector<uint32_t> bpsnp;
            vector<double> freq, byz, seyz;
            vector<int> sn_ids;
            VectorXd zsxz,_byz,_seyz,_bxz,_sexz,_zsxz,ld_v;
            MatrixXd _X, _LD,_LD_heidi,_X_heidi;
  
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
                        
                        int etrait_rid=etrait._esi_include[ge_rowid];
                        long pos=find(gdata._include.begin(), gdata._include.end(), etrait_rid)-gdata._include.begin();
                        /****the same **/
                        //string ers=esdata._esi_rs[ge_rowid];
                        //string trs=etrait._esi_rs[etrait._esi_include[ge_rowid]];
                        //long pos=find(gdata.snpName.begin(), gdata.snpName.end(), trs)-gdata.snpName.begin();
                        /******/
                        if(snpchr==probechr && ABS(probebp-snpbp)<=cis_itvl && pos!=gdata._include.size())
                        {
                            bxz.push_back(esdata._val[beta_start+j]);
                            sexz.push_back(esdata._val[se_start+j]);
                            byz.push_back(gdata.byz[pos]);
                            seyz.push_back(gdata.seyz[pos]);
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
                
                Map<VectorXd> ei_bxz(&bxz[0],bxz.size());
                Map<VectorXd> ei_sexz(&sexz[0],sexz.size());
                
                zsxz=ei_bxz.array()/ei_sexz.array();
                if(!heidiFlag) maxid=max_abs_id(zsxz);
                double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
                
                
                if(!heidiFlag && pxz_val>p_smr) continue;
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
                if(heidioffFlag)   rsfreq[outCount]=string("NA");
                else rsfreq[outCount]=dtosf(freq[maxid]);
                
                rsa1[outCount]=allele1[maxid];
                rsa2[outCount]=allele2[maxid];
                
                if(heidioffFlag || chisqxy<pre_psmr_thresh)
                {
                    prb1[outCount]= string("NA");
                    nsnp_test1[outCount]= string("NA");
                }
                else
                {
                    make_XMat(&bdata,curId, _X);
                    ld_calc_o2m(ld_v,maxid,_X);
                    
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
                
            }
            
            //genesn in R is just probidx here. can refer to probinfo to get probeid, probenm, chr, gene, bp
         
                if(outCount>=0)
                {
                    string smrfile = string(outFileName)+"."+traitname+".smr";
                    ofstream smr(smrfile.c_str());
                    if (!smr) throw ("Error: can not open the file " + smrfile + " to save!");
                    
                    smr << "Expo_ID" <<'\t'<< "Expo_Chr" <<'\t' << "Expo_Gene"  << '\t' << "Expo_bp" << '\t'<< "Outco_ID" <<'\t'<< "Outco_Chr" <<'\t' << "Outco_Gene"  << '\t' << "Outco_bp" << '\t'<< "SNP"<< '\t' << "SNP_Chr"<< '\t' << "SNP_bp"<< '\t' << "A1"<< '\t'<< "A2"<< '\t'<<"Freq"<<'\t'<<"b_Outco"<<'\t'<<"se_Outco"<<'\t'<< "p_Outco" << '\t'<<"b_Expo"<<'\t'<<"se_Expo"<<'\t'<< "p_Expo" << '\t'<< "b_SMR" << '\t'<< "se_SMR"<< '\t' << "p_SMR" << "\t"<< "p_HET"<< "\t" << "nsnp" << '\n';
                    
                    for (int i = 0;i <=outCount; i++) {                        
                        smr<<esdata._epi_prbID[out_probid[i]]<<'\t'<<esdata._epi_chr[out_probid[i]]<<'\t'<<esdata._epi_gene[out_probid[i]]<<'\t'<<esdata._epi_bp[out_probid[i]]<<'\t'<<traitname<<'\t'<<etrait._epi_chr[ii]<<'\t'<<etrait._epi_gene[ii]<<'\t'<<etrait._epi_bp[ii]<<'\t'<<rsid[i]<<'\t'<<rschr[i]<<'\t'<<rsbp[i]<<'\t'<<rsa1[i]<<'\t'<<rsa2[i]<<'\t'<<rsfreq[i]<<'\t'<<bgwas[i]<<'\t'<<segwas[i]<<'\t'<<pgwas[i]<<'\t'<<beqtl[i]<<'\t'<<seeqtl[i]<<'\t'<<peqtl[i]<<'\t'<<bxy[i]<<'\t'<<sexy[i]<<'\t'<<pxy[i]<<'\t'<<prb1[i]<<'\t'<<nsnp_test1[i]<<'\n';
                    }
                    cout<<"SMR and heterogeneity analysis finished.\nSMR and heterogeneity analysis results of "<<outCount+1<<" probes have been saved in the file [" + smrfile + "]."<<endl;
                    smr.close();
                    
                }else cout<<"SMR and heterogeneity analysis finished.\nSMR and heterogeneity analysis results of "<<outCount+1<<" probes have been saved."<<endl;
            
            
            free_gwas_data( &gdata);
            
        }
    }
    
    void update_freq(char* eqtlFileName, char* frqfile)
    {
        eqtlInfo eqtlinfo;
         if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        
        read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
        ifstream frq(frqfile);
        if (!frq) throw ("Error: can not open the file [" + string(frqfile) + "] to read.");
        cout << "Reading eQTL probe information from [" + string(frqfile) + "]." << endl;
        FILE* fsnpfptr=fopen((string(frqfile)+".unmatched.snp.list").c_str(), "w");
        if (!(fsnpfptr)) {
            printf("Error: Failed to open unmatched SNP log file.\n");
        }

        vector<string> vs_buf;
        map<string, int>::iterator iter;
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        int fsnpNum(0);
        int falleleNum(0);
        int hitNum(0);
        while(!frq.eof())
        {
            frq.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, " \t\n");
                if(col_num!=4) {
                    printf("ERROR: column number is not right in row %d!\n", lineNum+1);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na" || vs_buf[0]=="0"){
                    printf("ERROR: \"NA\" SNP found in row %d.\n",lineNum+1);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[3]=="NA" || vs_buf[3]=="na" ){
                    printf("ERROR: \"NA\" SNP freq found in row %d.\n",lineNum+1);
                    exit(EXIT_FAILURE);
                }
                if(atof(vs_buf[3].c_str())>=1 || atof(vs_buf[3].c_str())<=0)
                {
                    printf("ERROR: freqency should be between (0 , 1) in row %d.\n",lineNum+1);
                    exit(EXIT_FAILURE);
                }
                iter=eqtlinfo._snp_name_map.find(vs_buf[0]);
                if(iter!=eqtlinfo._snp_name_map.end())
                {
                    int snpidx=iter->second;
                    if(eqtlinfo._esi_allele1[snpidx]==vs_buf[1] && eqtlinfo._esi_allele2[snpidx]==vs_buf[2]) {
                        eqtlinfo._esi_freq[snpidx]=atof(vs_buf[3].c_str());
                        hitNum++;
                    }
                    else if(eqtlinfo._esi_allele1[snpidx]==vs_buf[2] && eqtlinfo._esi_allele2[snpidx]==vs_buf[1])
                    {
                        eqtlinfo._esi_freq[snpidx]=1-atof(vs_buf[3].c_str());
                        hitNum++;
                    }
                    else {
                        falleleNum++;
                        string logstr=vs_buf[0]+'\t'+"failedAlleleCheck"+'\n';
                        fputs(logstr.c_str(),fsnpfptr);
                    }
                    
                } else {
                    fsnpNum++;
                    string logstr=vs_buf[0]+'\t'+"failedSNPMatch"+'\n';
                    fputs(logstr.c_str(),fsnpfptr);
                }
                lineNum++;
            }
        }
        fclose(fsnpfptr);
        frq.close();
        cout << lineNum << " SNPs to be included from [" + string(frqfile) + "]." << endl;
        printf("%d SNPs hitted, %d SNPs failed in SNP match and %d SNPs failed in allele check.\n",hitNum,fsnpNum,falleleNum);
        
        string esifile = string(eqtlFileName)+".esi";
        ofstream smr(esifile.c_str());
        if (!smr) throw ("Error: can not open the esi file " + esifile + " to save!");
        for (int i = 0;i <eqtlinfo._snpNum; i++) {
            smr<<eqtlinfo._esi_chr[i]<<'\t'<<eqtlinfo._esi_rs[i]<<'\t'<<eqtlinfo._esi_gd[i]<<'\t'<<eqtlinfo._esi_bp[i]<<'\t'<<eqtlinfo._esi_allele1[i]<<'\t'<<eqtlinfo._esi_allele2[i]<<'\t'<<eqtlinfo._esi_freq[i]<<'\n';
        }
        smr.close();
        cout<<eqtlinfo._snpNum<<" SNPs have been saved in the file [" + esifile + "]."<<endl;
        
       
    }
    
    void write_besd(string outFileName, eqtlInfo* eqtlinfo)
    {
        filter_probe_null(eqtlinfo); // at the same time, reset the vector _include
        filter_snp_null(eqtlinfo);
        cout<<"\nsaving eQTL data..."<<endl;
        string esdfile = string(outFileName)+".esi";
        ofstream smr(esdfile.c_str());
        if (!smr) throw ("Error: can not open the esi file " + esdfile + " to save!");
        for (int i = 0;i <eqtlinfo->_snpNum; i++) {
            smr<<eqtlinfo->_esi_chr[i]<<'\t'<<eqtlinfo->_esi_rs[i]<<'\t'<<eqtlinfo->_esi_gd[i]<<'\t'<<eqtlinfo->_esi_bp[i]<<'\t'<<eqtlinfo->_esi_allele1[i]<<'\t'<<eqtlinfo->_esi_allele2[i]<<'\t'<<eqtlinfo->_esi_freq[i]<<'\n';
        }
        smr.close();
        cout<<eqtlinfo->_snpNum<<" SNPs have been saved in the file [" + esdfile + "]."<<endl;
        
        esdfile = string(outFileName)+".epi";
        smr.open(esdfile.c_str());
        if (!smr) throw ("Error: can not open the epi file " + esdfile + " to save!");
        for (int i = 0;i <eqtlinfo->_include.size(); i++) {
            smr<<eqtlinfo->_epi_chr[eqtlinfo->_include[i]]<<'\t'<<eqtlinfo->_epi_prbID[eqtlinfo->_include[i]]<<'\t'<<eqtlinfo->_epi_gd[eqtlinfo->_include[i]]<<'\t'<<eqtlinfo->_epi_bp[eqtlinfo->_include[i]]<<'\t'<<eqtlinfo->_epi_gene[eqtlinfo->_include[i]]<<'\t'<<eqtlinfo->_epi_orien[eqtlinfo->_include[i]]<<'\n';
        }
        smr.close();
        cout<<eqtlinfo->_include.size()<<" probes have been saved in the file [" + esdfile + "]."<<endl;
        
        esdfile = string(outFileName)+".besd";
        FILE * smrbesd;
        smrbesd = fopen (esdfile.c_str(), "wb");
        if(eqtlinfo->_valNum==0)
        {
            uint64_t bsize=(eqtlinfo->_include.size()*eqtlinfo->_snpNum<<1)+1;
            float* buffer=(float*)malloc (sizeof(float)*bsize);
            memset(buffer,0,sizeof(float)*bsize);
            float* ptr=buffer;
            *ptr++=0.0;
            uint64_t pro_num=eqtlinfo->_include.size();
            uint64_t snp_num=eqtlinfo->_snpNum;
            for(int i=0;i<pro_num;i++)
            {
                memcpy(ptr+(i<<1)*snp_num,&eqtlinfo->_bxz[eqtlinfo->_include[i]][0],sizeof(float)*snp_num);
                memcpy(ptr+((i<<1)+1)*snp_num,&eqtlinfo->_sexz[eqtlinfo->_include[i]][0],sizeof(float)*snp_num);
            }
            fwrite(buffer,sizeof(float), bsize, smrbesd);
            free(buffer);
        }
        else
        {
            uint64_t colSize=sizeof(uint64_t)*((eqtlinfo->_include.size()<<1)+1);
            uint64_t rowSize=sizeof(uint32_t)*eqtlinfo->_valNum;
            uint64_t valSize=sizeof(float)*eqtlinfo->_valNum;
            uint64_t valNum=eqtlinfo->_valNum;
            uint64_t bufsize=sizeof(float)+sizeof(uint64_t)+colSize+rowSize+valSize;
            
            char* buffer=(char*)malloc (sizeof(char)*bufsize);
            memset(buffer,0,sizeof(char)*bufsize);
            float ftype=SPARSE_FILE_TYPE_3;
            memcpy(buffer,&ftype,sizeof(float));
            char* wptr=buffer+sizeof(float);
            memcpy(wptr,&valNum,sizeof(uint64_t));
            wptr+=sizeof(uint64_t);
            uint64_t* uptr=(uint64_t*)wptr; *uptr++=0;
            for(int i=0;i<eqtlinfo->_include.size();i++)
            {
                *uptr++=eqtlinfo->_cols[(eqtlinfo->_include[i]<<1)+1];
                *uptr++=eqtlinfo->_cols[eqtlinfo->_include[i]+1<<1];
            }
            wptr+=colSize;
            memcpy(wptr,&eqtlinfo->_rowid[0],rowSize);
            wptr+=rowSize;
            memcpy(wptr,&eqtlinfo->_val[0],valSize);
            fwrite (buffer,sizeof(char), bufsize, smrbesd);
            free(buffer);

        }

        
        fclose (smrbesd);
        
        cout<<"Beta values and SE values for "<<eqtlinfo->_include.size()<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs have been saved in the binary file [" + esdfile + "]." <<endl;
        

    }
    
    
    void meta(char* outFileName,char* eqtlFileName, char* eqtlFileName2)
    {
        setNbThreads(thread_num);
        
        eqtlInfo etrait;
        eqtlInfo esdata;
        eqtlInfo metadata;
        
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        if(eqtlFileName2==NULL) throw("Error: please input another eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        
        read_esifile(&etrait, string(eqtlFileName)+".esi");
        read_esifile(&esdata, string(eqtlFileName2)+".esi");
        vector<int> idx;
        vector<string> cmmnSymbs;
        match_only(esdata._esi_rs, etrait._esi_rs, idx);
        etrait._esi_include=idx;
        for(int i=0;i<idx.size();i++) cmmnSymbs.push_back(etrait._esi_rs[idx[i]]);
        idx.clear();
        match_only(cmmnSymbs,esdata._esi_rs,idx);
        esdata._esi_include=idx;
        cout<<"There are "<<idx.size()<<" SNPs in common."<<endl;
        read_epifile(&etrait, string(eqtlFileName)+".epi");
        read_epifile(&esdata, string(eqtlFileName2)+".epi");
        idx.clear();
        match_only(esdata._epi_prbID, etrait._epi_prbID, idx);
        etrait._include=idx;
        cmmnSymbs.clear();
        for(int i=0;i<idx.size();i++) cmmnSymbs.push_back(etrait._epi_prbID[idx[i]]);
        idx.clear();
        match_only(cmmnSymbs,esdata._epi_prbID,idx);
        esdata._include=idx;
        cout<<"There are "<<idx.size()<<" Probes in common."<<endl;
        read_besdfile(&etrait, string(eqtlFileName)+".besd");
        if(etrait._rowid.empty() && etrait._bxz.empty())
        {
            printf("No data included from %s under current condition.\n",eqtlFileName);
            exit(EXIT_FAILURE);
        }
        read_besdfile(&esdata, string(eqtlFileName2)+".besd");
        if(esdata._rowid.empty() && esdata._bxz.empty())
        {
            printf("No data included from %s under current condition.\n",eqtlFileName2);
            exit(EXIT_FAILURE);
        }
        metadata._cols.push_back(0);
        cout<<"Performing Meta analysis..."<<endl;
        map<string, int> unmatched_rs_map;
        int unmatched_map_size=0;
        for(int i=0;i<etrait._probNum;i++)
        {
            string ref_pid=etrait._epi_prbID[i];
            string alt_pid=esdata._epi_prbID[i];
            if(ref_pid!=alt_pid)
            {
                cout<<"Some bugs here, please help to report!"<<endl;
                exit(1);
            }
            
            vector<float> ref_byz;
            vector<float> ref_seyz;
            vector<string> ref_a1;
            vector<string> ref_a2;
            vector<uint32_t> ref_rowid;
            vector<float> meta_beta;
            vector<float> meta_se;
            vector<uint32_t> meta_rowid;
            
            uint64_t beta_start=etrait._cols[i<<1];
            uint64_t se_start=etrait._cols[1+(i<<1)];
            uint64_t numsnps=se_start-beta_start;
            for(uint64_t j=0;j<numsnps;j++)
            {
                uint32_t ge_rowid=etrait._rowid[beta_start+j];
                ref_rowid.push_back(ge_rowid);
                ref_a1.push_back(etrait._esi_allele1[ge_rowid]);
                ref_a2.push_back(etrait._esi_allele2[ge_rowid]);
                ref_byz.push_back(etrait._val[beta_start+j]);
                ref_seyz.push_back(etrait._val[se_start+j]);
            }
            
            vector<float> byz;
            vector<float> seyz;
            vector<string> a1;
            vector<string> a2;
            vector<uint32_t> rowid;
          
            
            uint64_t beta_stt=esdata._cols[i<<1];
            uint64_t se_stt=esdata._cols[1+(i<<1)];
            uint64_t nums=se_stt-beta_stt;
            for(uint64_t j=0;j<nums;j++)
            {
                uint32_t ge_rid=esdata._rowid[beta_stt+j];
                rowid.push_back(ge_rid);
                a1.push_back(esdata._esi_allele1[ge_rid]);
                a2.push_back(esdata._esi_allele2[ge_rid]);
                byz.push_back(esdata._val[beta_stt+j]);
                seyz.push_back(esdata._val[se_stt+j]);
               
            }
            
            vector<int> idx;
            match(ref_rowid,rowid,idx);
            for(int j=0;j<idx.size();j++)
            {
                if(idx[j]==-9) continue;
                uint32_t tmp_rowid=ref_rowid[j];
                string tmp_rs=etrait._esi_rs[tmp_rowid];
                if(tmp_rowid!=rowid[idx[j]])
                {
                    cout<<"Some bugs here, please help to report!"<<endl;
                    exit(1);
                }
                
                string tmp_ref_a1=ref_a1[j];
                string tmp_ref_a2=ref_a2[j];
                string tmp_alt_a1=a1[idx[j]];
                string tmp_alt_a2=a2[idx[j]];
                float se1=ref_seyz[j];
                float se2=seyz[idx[j]];
                float beta1=ref_byz[j];
                float beta2=byz[idx[j]];
                
                
                if (tmp_ref_a1 == tmp_alt_a1 && tmp_ref_a2 == tmp_alt_a2)
                {
                    float tmpSE=se1*se2/sqrt(se1*se1+se2*se2);
                    meta_rowid.push_back(tmp_rowid);
                    meta_se.push_back(tmpSE);
                    meta_beta.push_back((beta1/(se1*se1) + beta2/(se2*se2))*tmpSE*tmpSE);
                }
                else if(tmp_ref_a1 == tmp_alt_a2 && tmp_ref_a2 == tmp_alt_a1)
                {
                    beta2=-beta2;
                    float tmpSE=se1*se2/sqrt(se1*se1+se2*se2);
                    meta_rowid.push_back(tmp_rowid);
                    meta_se.push_back(tmpSE);
                    meta_beta.push_back((beta1/(se1*se1) + beta2/(se2*se2))*tmpSE*tmpSE);
                }
                else {
                    unmatched_rs_map.insert(pair<string, int>(tmp_rs, unmatched_map_size));
                   
                }
            }
        
            for(int j=0;j<meta_beta.size();j++) metadata._val.push_back(meta_beta[j]);
            for(int j=0;j<meta_beta.size();j++) metadata._val.push_back(meta_se[j]);
            for(int j=0;j<meta_beta.size();j++) metadata._rowid.push_back(meta_rowid[j]);
            for(int j=0;j<meta_beta.size();j++) metadata._rowid.push_back(meta_rowid[j]);
            metadata._cols.push_back(meta_beta.size()+metadata._cols[metadata._cols.size()-1]);
            metadata._cols.push_back(meta_beta.size()+metadata._cols[metadata._cols.size()-1]);
           
        }
        metadata._esi_allele1=etrait._esi_allele1;
        metadata._esi_allele2=etrait._esi_allele2;
        metadata._esi_bp=etrait._esi_bp;
        metadata._esi_chr=etrait._esi_chr;
        metadata._esi_gd=etrait._esi_gd;
        metadata._esi_include=etrait._esi_include;
        metadata._esi_rs=etrait._esi_rs;
        metadata._epi_bp=etrait._epi_bp;
        metadata._epi_chr=etrait._epi_chr;
        metadata._epi_gd=etrait._epi_gd;
        metadata._epi_gene=etrait._epi_gene;
        metadata._epi_orien=etrait._epi_orien;
        metadata._epi_prbID=etrait._epi_prbID;
        metadata._include=etrait._include;
        metadata._probNum=etrait._probNum;
        metadata._snpNum=etrait._snpNum;
        metadata._valNum=metadata._val.size();
        
        string unmatchedsnpfname = string(outFileName)+".unmatched.snp.list";
        FILE* unmatchedsnpfile=fopen(unmatchedsnpfname.c_str(), "w");
        if (!(unmatchedsnpfile)) {
            printf("Error: Failed to open unmatchedsnpfile file.\n");
            exit(1);
        }
        for(std::map<string,int>::iterator it=unmatched_rs_map.begin(); it!=unmatched_rs_map.end(); ++it){
            string tmpstr=it->first+'\n';
            fputs(tmpstr.c_str(),unmatchedsnpfile);
        }
        fclose(unmatchedsnpfile);
        write_besd(outFileName, &metadata);
        
    }
    
    void esi_man(eqtlInfo* eqtlinfo,char* snplstName,int chr,int snpchr, char* snprs, char* fromsnprs, char* tosnprs,int snpWind,int fromsnpkb, int tosnpkb,bool snpwindFlag,bool cis_flag, int cis_itvl,const char* prbname)
    {
        string logstr;
        int flags4snp=0;
        if(snplstName != NULL) flags4snp++;
        if(snprs != NULL) flags4snp++;
        if(fromsnprs!=NULL) flags4snp++;
        if(fromsnpkb>=0) flags4snp++;
        if(flags4snp>1)
        {
            logstr="WARNING: Flags for SNPs in this section are mutual exclusive. The priority order (from high to low) is: --extract-snp, --snp-wind, --snp, --from(to)--snp, --from(to)-snp-kb.\n";
            fputs(logstr.c_str(), stdout);
        }
        if(snpchr!=0)
        {
            extract_eqtl_by_chr(eqtlinfo, snpchr);
        }
        else if(chr!=0)
        {
            extract_eqtl_by_chr(eqtlinfo, chr);
        }
        
        if(prbname!=NULL && cis_flag)
        {
            extract_eqtl_snp(eqtlinfo, prbname, cis_itvl, "probe"); // extract cis eQTLs
        }
        else if (snplstName != NULL) extract_eqtl_snp(eqtlinfo, snplstName);
        else if (snpwindFlag)
        {
            if(snprs==NULL)
            {
                logstr="ERROR: Please identify the SNP name by --snp when using --snp-wind.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            extract_eqtl_snp(eqtlinfo, snprs, snpWind, "SNP");
        }
        else if(snprs!=NULL)
        {
            extract_eqtl_single_snp(eqtlinfo, snprs);
        }
        else if(fromsnprs!=NULL)
        {
            if(tosnprs==NULL)
            {
                logstr="ERROR: Please identify the SNP name by --to-snp.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            extract_eqtl_snp(eqtlinfo, fromsnprs, tosnprs);
        }
        else if(fromsnpkb>=0)
        {
            
            if(fromsnpkb>=0 && chr==0 && snpchr==0) {
                logstr="ERROR: Please identify the chromosome by --snp-chr or --chr.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            
            if(tosnpkb<0)
            {
                logstr="ERROR: snp BP can't be negative.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            extract_eqtl_snp(eqtlinfo, snpchr, fromsnpkb, tosnpkb);
        }
        
    }
    void epi_man(eqtlInfo* eqtlinfo,char* problstName,char* genelistName, int chr,int prbchr, const char* prbname, char* fromprbname, char* toprbname,int prbWind,int fromprbkb, int toprbkb,bool prbwindFlag, char* genename)
    {
        string logstr;
        int flags4prb=0;
        if(problstName != NULL) flags4prb++;
        if(prbname != NULL) flags4prb++;
        if(fromprbname!=NULL) flags4prb++;
        if(fromprbkb>=0) flags4prb++;
        if(genename != NULL) flags4prb++;
        if(flags4prb>1)
        {
            logstr="WARNING: Flags for probes in this section are mutual exclusive. The priority order (from high to low) is: --extract-probe, --gene-list, --probe-wind, --probe, --from(to)--probe, --from(to)-probe-kb, --gene.\n";
            fputs(logstr.c_str(), stdout);
        }
        
        if(prbchr!=0)
        {
            extract_epi_by_chr(eqtlinfo, prbchr);
        }
        else if(chr!=0)
        {
            extract_epi_by_chr(eqtlinfo, chr);
        }
        
        if(problstName != NULL || genelistName != NULL)
        {
            if(problstName != NULL) extract_prob(eqtlinfo, problstName);
            if(genelistName != NULL) extract_prob_by_gene(eqtlinfo, genelistName);
        }
        else if(prbwindFlag)
        {
            if(prbname==NULL)
            {
                logstr="ERROR: Please identify the probe name by --probe when using --probe-wind.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            extract_prob(eqtlinfo, prbname, prbWind);
        }
        else if(prbname!=NULL)
        {
            extract_eqtl_single_probe(eqtlinfo, prbname);
        }
        else if(fromprbname!=NULL)
        {
            if(toprbname==NULL)
            {
                logstr="ERROR: Please identify the probe name by --to-probe.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            extract_eqtl_prob(eqtlinfo, fromprbname, toprbname);
        }
        else if(fromprbkb>=0)
        {
            if(fromprbkb>=0 && chr==0 && prbchr==0) {
                logstr="ERROR: Please identify the chromosome by --probe-chr or --chr.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            if(toprbkb<0)
            {
                logstr="ERROR: probe BP can't be negative.\n";
                fputs(logstr.c_str(), stdout);
                exit(1);
            }
            extract_eqtl_prob(eqtlinfo, prbchr, fromprbkb, toprbkb);
        }
        else if(genename!=NULL)
        {
            
            extract_prob_by_single_gene(eqtlinfo, genename);
        }
        
    }
    
}
