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
        vector<char> a1_buf, a2_buf, ref_A_buf, other_A_buf;
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

		vector<uint32_t> chr_buf, gd_buf,bp_buf;
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
        for (int i = 0; i < eqtlinfo->_probNum; i++) eqtlinfo->_include.push_back(i);

		/*		
		eqtlinfo->_epi_chr.assign(chr_buf.begin(),chr_buf.end());
		eqtlinfo->_epi_gd.assign(gd_buf.begin(), gd_buf.end());
		eqtlinfo->_epi_bp.assign(bp_buf.begin(), bp_buf.end());
		eqtlinfo->_epi_prbID.assign(prbID_buf.begin(), prbID_buf.end());
		eqtlinfo->_epi_gene.assign(gene_buf.begin(), gene_buf.end());
		eqtlinfo->_epi_orien.assign(orien_buf.begin(), orien_buf.end());
		*/
	}

	void update_esi(eqtlInfo* eqtlinfo)
	{
		eqtlinfo->_snpNum = eqtlinfo->_esi_include.size();	

		vector<uint32_t> chr_buf, gd_buf, bp_buf;
		vector<string> rs_buf;
		vector<char> allele1_buf, allele2_buf;
		for (int i = 0; i < eqtlinfo->_snpNum; i++)
		{
			chr_buf.push_back(eqtlinfo->_esi_chr[eqtlinfo->_esi_include[i]]);
			gd_buf.push_back(eqtlinfo->_esi_gd[eqtlinfo->_esi_include[i]]);
			bp_buf.push_back(eqtlinfo->_esi_bp[eqtlinfo->_esi_include[i]]);
			rs_buf.push_back(eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]]);
			allele1_buf.push_back(eqtlinfo->_esi_allele1[eqtlinfo->_esi_include[i]]);
			allele2_buf.push_back(eqtlinfo->_esi_allele2[eqtlinfo->_esi_include[i]]);
		}
		eqtlinfo->_esi_chr.clear();
		eqtlinfo->_esi_gd.clear();
		eqtlinfo->_esi_bp.clear();
		eqtlinfo->_esi_rs.clear();
		eqtlinfo->_esi_allele1.clear();
		eqtlinfo->_esi_allele2.clear();
		eqtlinfo->_esi_chr.swap(chr_buf);
		eqtlinfo->_esi_gd.swap(gd_buf);
		eqtlinfo->_esi_bp.swap(bp_buf);
		eqtlinfo->_esi_rs.swap(rs_buf);
		eqtlinfo->_esi_allele1.swap(allele1_buf);
		eqtlinfo->_esi_allele2.swap(allele2_buf);
        
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_snp_name_map.clear();
        long size=0;
        for (int i = 0; i < eqtlinfo->_snpNum; i++)
        {
            eqtlinfo->_esi_include.push_back(i);
            /*
             // disabled currently
            eqtlinfo->_snp_name_map.insert(pair<string, int>(eqtlinfo->_esi_rs[i], i));
            if (size == eqtlinfo->_snp_name_map.size()) throw ("Error: Duplicated SNP IDs found: \"" + eqtlinfo->_esi_rs[i] + "\".");
            size = eqtlinfo->_snp_name_map.size();
             */
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
            bdata->_allele1.push_back(cbuf.c_str()[0]);
            Bim >> cbuf;
            StrFunc::to_upper(cbuf);
            bdata->_allele2.push_back(cbuf.c_str()[0]);
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
        int size = 0;
        for (int i = 0; i <  bdata->_snp_num; i++) {
             bdata->_include[i] = i;
             bdata->_snp_name_map.insert(pair<string, int>( bdata->_snp_name[i], i));
            if (size ==  bdata->_snp_name_map.size()) throw ("Error: Duplicated SNP IDs found: \"" +  bdata->_snp_name[i] + "\".");
            size =  bdata->_snp_name_map.size();
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
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        gwasFile.getline(buf,MAX_LINE_SIZE);// the header
        while(!gwasFile.eof())
        {
            gwasFile.getline(buf,MAX_LINE_SIZE);
            lineNum++;
        }
        if(buf[0]=='\0') lineNum--;
        
        gdata->snpNum=lineNum;
        gdata->snpName.resize(lineNum);
        gdata->allele_1=(char*)malloc(lineNum*sizeof(char));
        gdata->allele_2=(char*)malloc(lineNum*sizeof(char));
        gdata->freq=(double*)malloc(lineNum*sizeof(double));
        gdata->byz=(double*)malloc(lineNum*sizeof(double));
        gdata->seyz=(double*)malloc(lineNum*sizeof(double));
        gdata->pvalue=(double*)malloc(lineNum*sizeof(double));
        gdata->splSize=(uint32_t*)malloc(lineNum*sizeof(uint32_t));
        gdata->_include.resize(lineNum);
        
        gwasFile.clear(ios::goodbit);
        gwasFile.seekg (0, ios::beg);
        gwasFile.getline(buf,MAX_LINE_SIZE);
        for(int i=0;i<lineNum;i++)
        {
            string tmpStr;
            gwasFile.getline(buf,MAX_LINE_SIZE);
            istringstream iss(buf);
            iss>>tmpStr;
            gdata->snpName[i]=tmpStr;
            iss>>tmpStr;
            memcpy(&gdata->allele_1[i],tmpStr.c_str(),1);
            iss>>tmpStr;
            memcpy(&gdata->allele_2[i],tmpStr.c_str(),1);
            iss>>tmpStr;
            gdata->freq[i]=atof(tmpStr.c_str());
            iss>>tmpStr;
            gdata->byz[i]=atof(tmpStr.c_str());
            iss>>tmpStr;
            gdata->seyz[i]=atof(tmpStr.c_str());
            iss>>tmpStr;
            gdata->pvalue[i]=atof(tmpStr.c_str());
            iss>>tmpStr;
            gdata->splSize[i]=atoi(tmpStr.c_str());
            gdata->_include[i]=i;
        }
        to_upper(gdata->allele_1,lineNum);
        to_upper(gdata->allele_2,lineNum);
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
    void read_esifile(eqtlInfo* eqtlinfo, string esifile)
    {
        ifstream esi(esifile.c_str());
        if (!esi) throw ("Error: can not open the file [" + esifile + "] to read.");
        cout << "Reading eQTL SNP information from [" + esifile + "]." << endl;
        eqtlinfo->_esi_chr.clear();
        eqtlinfo->_esi_rs.clear();
        eqtlinfo->_esi_gd.clear();
        eqtlinfo->_esi_bp.clear();
        eqtlinfo->_esi_allele1.clear();
        eqtlinfo->_esi_allele2.clear();
		eqtlinfo->_esi_include.clear();

        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        while(!esi.eof())
        {
            esi.getline(buf,MAX_LINE_SIZE);
            lineNum++;
        }
        if(buf[0]=='\0') lineNum--;
        eqtlinfo->_snpNum=lineNum;        
        cout << eqtlinfo->_snpNum << " SNPs to be included from [" + esifile + "]." << endl;
        
        eqtlinfo->_esi_chr.resize(lineNum);
        eqtlinfo->_esi_rs.resize(lineNum);
        eqtlinfo->_esi_gd.resize(lineNum);
        eqtlinfo->_esi_bp.resize(lineNum);
        eqtlinfo->_esi_allele1.resize(lineNum);
        eqtlinfo->_esi_allele2.resize(lineNum);
		eqtlinfo->_esi_include.resize(lineNum);
        esi.clear(ios::goodbit);
        esi.seekg (0, ios::beg);
        long size=0;
        for(int i=0;i<lineNum;i++)
        {
            string tmpStr;
            esi.getline(buf,MAX_LINE_SIZE);
            istringstream iss(buf);
            iss>>tmpStr;
            eqtlinfo->_esi_chr[i]=atoi(tmpStr.c_str());
            iss>>tmpStr;
            eqtlinfo->_esi_rs[i]=tmpStr;
            /*
             //disabled currently
            eqtlinfo->_snp_name_map.insert(pair<string, int>(tmpStr, i));
            if (size == eqtlinfo->_snp_name_map.size()) throw ("Error: Duplicated SNP IDs found: \"" + tmpStr + "\".");
            size = eqtlinfo->_snp_name_map.size();
            */
            iss>>tmpStr;
            eqtlinfo->_esi_gd[i]=atoi(tmpStr.c_str());
            iss>>tmpStr;
            eqtlinfo->_esi_bp[i]=atoi(tmpStr.c_str());
            iss>>tmpStr;
            StrFunc::to_upper(tmpStr);
            eqtlinfo->_esi_allele1[i]=tmpStr.c_str()[0];
            iss>>tmpStr;
            StrFunc::to_upper(tmpStr);
            eqtlinfo->_esi_allele2[i]=tmpStr.c_str()[0];
			eqtlinfo->_esi_include[i] = i;
            
        }       
        esi.close();
       
    }
    void read_epifile(eqtlInfo* eqtlinfo, string epifile)
    {
        ifstream epi(epifile.c_str());
        if (!epi) throw ("Error: can not open the file [" + epifile + "] to read.");
        cout << "Reading eQTL probe information from [" + epifile + "]." << endl;
        eqtlinfo->_epi_chr.clear();
        eqtlinfo->_epi_prbID.clear();
        eqtlinfo->_epi_gd.clear();
        eqtlinfo->_epi_bp.clear();
        eqtlinfo->_epi_gene.clear();
        eqtlinfo->_epi_orien.clear();
		eqtlinfo->_include.clear();

        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        while(!epi.eof())
        {
            epi.getline(buf,MAX_LINE_SIZE);
            lineNum++;
        }
        if(buf[0]=='\0') lineNum--;
        eqtlinfo->_probNum=lineNum;
        cout << eqtlinfo->_probNum << " Probes to be included from [" + epifile + "]." << endl;
        
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
            eqtlinfo->_epi_chr[i]=atoi(tmpStr.c_str());
            iss>>tmpStr;
            eqtlinfo->_epi_prbID[i]=tmpStr;
            iss>>tmpStr;
            eqtlinfo->_epi_gd[i]=atoi(tmpStr.c_str());
            iss>>tmpStr;
            eqtlinfo->_epi_bp[i]=atoi(tmpStr.c_str());
            iss>>tmpStr;
            eqtlinfo->_epi_gene[i]=tmpStr.c_str();
            iss>>tmpStr;
            eqtlinfo->_epi_orien[i]=tmpStr.c_str()[0];
            eqtlinfo->_include[i]=i;
        }
        
        epi.close();
    }    
    
    void read_besdfile(eqtlInfo* eqtlinfo, string besdfile)
    {
        if (eqtlinfo->_include.size() == 0) throw ("Error: No probe is retained for analysis.");
        if (eqtlinfo->_esi_include.size() == 0) throw ("Error: No SNP is retained for analysis.");
        
        // the fastest way is using malloc and memcpy
        char buf[MAX_LINE_SIZE+4];
        ifstream besd(besdfile.c_str(), ios::in|ios::binary);
        if(!besd)
        {
            fprintf (stderr, "%s: Couldn't open file %s\n",
                     besdfile.c_str(), strerror (errno));
            exit (EXIT_FAILURE);
        }
        cout << "Reading eQTL summary-level statistics from [" + besdfile + "]." << endl;
        
        besd.read(buf, 4);
        float* flag;
        flag=(float *)buf;
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
            besd.read(buf, sizeof(uint64_t));
            valNum=*(uint64_t *)buf;
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
                cout<<"eQTL summary-level statistics of "<<eqtlinfo->_include.size()<<" Probes and "<<eqtlinfo->_esi_include.size()<<" SNPs to be included from [" + besdfile + "]." <<endl;
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
                cout<<"eQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
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
                std::cout << "eQTL summary-level statistics of " << eqtlinfo->_include.size() << " Probes and " << eqtlinfo->_esi_include.size() << " SNPs to be included from [" + besdfile + "]." << endl;
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
                cout<<"\neQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
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
            besd.read(buf, 4);
            valNum=(uint32_t)*(float *)buf; // int to float then float to int back can lose pricision. hence this clause and bellow are unbelievable
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
                cout<<"eQTL summary-level statistics of "<<eqtlinfo->_include.size()<<" Probes and "<<eqtlinfo->_esi_include.size()<<" SNPs to be included from [" + besdfile + "]." <<endl;
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
                cout<<"eQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
            }
            // terminate
            free (buffer);
        }
        besd.close();
    }
    
    void read_esdfile(eqtlInfo* eqtlinfo, string esdfile)
    {
		// clear datastruct for sparse befor read dense
		eqtlinfo->_cols.clear();
		eqtlinfo->_rowid.clear();
		eqtlinfo->_val.clear();
		eqtlinfo->_valNum = 0;

        ifstream esd;
        if(!file_read_check(&esd, esdfile.c_str()))
        {
            fprintf (stderr, "%s: Couldn't open file %s\n",
                     esdfile.c_str(), strerror (errno));
            exit (EXIT_FAILURE);
        }
        cout << "Reading eQTL summary-level statistics from [" + esdfile + "]." << endl;
        char buf[MAX_LINE_SIZE+4];
        int lineNum(0);      
        while(!esd.eof())
        {
            esd.getline(buf,MAX_LINE_SIZE);
            lineNum++;
        }
        if(buf[0]=='\0') lineNum--;
        if(lineNum != eqtlinfo->_snpNum) throw ("Error: no match of SNP numbers.");
               
		if (eqtlinfo->_include.size()<eqtlinfo->_probNum || eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum)
		{
			eqtlinfo->_bxz.resize(eqtlinfo->_include.size());
			eqtlinfo->_sexz.resize(eqtlinfo->_include.size());
			for (unsigned int i = 0; i<eqtlinfo->_include.size(); i++)
			{
				eqtlinfo->_bxz[i].reserve(eqtlinfo->_esi_include.size());
				eqtlinfo->_sexz[i].reserve(eqtlinfo->_esi_include.size());
			}

			esd.clear(ios::goodbit);
			esd.seekg(0, ios::beg);
			int cur_line = 0;
			for (int i = 0; i<eqtlinfo->_esi_include.size(); i++)
			{
				int sid=eqtlinfo->_esi_include[i];
				for (int j = cur_line; j < sid; j++) esd.getline(buf, MAX_LINE_SIZE);
				string tmpStr;				
				esd.getline(buf, MAX_LINE_SIZE);
				istringstream iss(buf);
				int pos = 0;
				for (int j = 0; j < eqtlinfo->_include.size(); j++)
				{
					for (int k = pos; k < eqtlinfo->_include[j]<<1; k++) iss >> tmpStr;
					iss >> tmpStr;
                    if(!tmpStr.compare("NA")) eqtlinfo->_bxz[j][i]=-9;
                    else eqtlinfo->_bxz[j][i] = atof(tmpStr.c_str());
					iss >> tmpStr;
                    if(!tmpStr.compare("NA")) eqtlinfo->_sexz[j][i] = -9;
                    else eqtlinfo->_sexz[j][i] = atof(tmpStr.c_str());
					pos = (eqtlinfo->_include[j]+1)<<1;
				}
				cur_line = sid + 1;
			}
            if(eqtlinfo->_include.size()<eqtlinfo->_probNum ) update_epi(eqtlinfo);
            if(eqtlinfo->_esi_include.size()<eqtlinfo->_snpNum) update_esi(eqtlinfo);
			/*
			for (int i = 0; i < eqtlinfo->_include.size(); i++)
			{
				for (int j = 0; j < eqtlinfo->_snpNum; j++)
					if (abs(eqtlinfo->_sexz[i][j] + 9) > 1e-6) cout << eqtlinfo->_bxz[i][j] << "," << eqtlinfo->_sexz[i][j] << endl;				
				cout << endl; cout << endl;
			}
			*/			
		}
		else
		{
			eqtlinfo->_bxz.resize(eqtlinfo->_probNum);
			eqtlinfo->_sexz.resize(eqtlinfo->_probNum);
			for (unsigned int i = 0; i<eqtlinfo->_probNum; i++)
			{
				eqtlinfo->_bxz[i].reserve(eqtlinfo->_snpNum);
				eqtlinfo->_sexz[i].reserve(eqtlinfo->_snpNum);
			}

			esd.clear(ios::goodbit);
			esd.seekg(0, ios::beg);
			for (int i = 0; i<lineNum; i++)
			{
				string tmpStr;
				esd.getline(buf, MAX_LINE_SIZE);
				istringstream iss(buf);
				for (int j = 0; j<eqtlinfo->_probNum; j++)
				{
					iss >> tmpStr;
                    if(!tmpStr.compare("NA")) eqtlinfo->_bxz[j][i] = -9;
					else eqtlinfo->_bxz[j][i] = atof(tmpStr.c_str());
					iss >> tmpStr;
                    if(!tmpStr.compare("NA")) eqtlinfo->_sexz[j][i] = -9;
					else eqtlinfo->_sexz[j][i] = atof(tmpStr.c_str());
				}
			}
		}
		     
		cout << "eQTL summary statistics of " << eqtlinfo->_include.size() << " Probes and " << eqtlinfo->_snpNum << " SNPs to be included from [" + esdfile + "]." << endl;
        esd.close();
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

    
  //this function derived from the source code of FastGCN
    void ld_calcualte(double* ref_ld,float* ref_snpData,int rsize,int csize)
    {
        vector<float> tmpX(rsize);
        vector<float> tmpX2(rsize);
        
        //
        //	int Num_CPU_Cores;
        //#ifdef WIN32
        //		SYSTEM_INFO si;
        //		GetSystemInfo(&si);
        //		Num_CPU_Cores=si.dwNumberOfProcessors;
        //#elif
        //		Num_CPU_Cores=sysconf(_SC_NPROCESSORS_CONF);
        //#endif
        //	omp_set_num_threads(Num_CPU_Cores);
        //
        //#pragma omp parallel for
        for(int i=0;i<rsize;i++){
            for(int j=0;j<csize;j++){
                tmpX[i] += ref_snpData[i*csize+j];
                tmpX2[i] += ref_snpData[i*csize+j]*ref_snpData[i*csize+j];
            }
        }
        float tmp0;
        // use lower triangular, because of i>j
        //#pragma omp parallel for private(tmp0)
        for(int i=0;i<rsize;i++){
            for(int j=0;j<i;j++){
                tmp0=0;
                for(int k=0;k<csize;k++) tmp0 += ref_snpData[i*csize+k]*ref_snpData[j*csize+k];
                ref_ld[i*(i-1)/2+j] = 
				(tmp0-tmpX[i]*tmpX[j]/csize)/sqrtf((tmpX2[i]-tmpX[i]*tmpX[i]/csize)*(tmpX2[j]-tmpX[j]*tmpX[j]/csize));
				
            }
        }
        
    }
    
    void get_square_idxes(vector<int> &sn_ids,VectorXd &zsxz,double threshold)
    {
        
        for(int i=0;i<zsxz.size();i++)
        {
           if(zsxz[i]*zsxz[i]-threshold>1e-6) sn_ids.push_back(i);
        }
    }
	void get_square_ldpruning_idxes(vector<int> &sn_ids, VectorXd &zsxz, double threshold, MatrixXd &LD, long maxid, double ld_top)
    {
       
        for(int i=0;i<zsxz.size();i++)
        {
            if(i!= maxid)
            {
                if((zsxz[i]*zsxz[i]-threshold)>1e-6 && (LD(maxid,i)-ld_top)<1e-6) sn_ids.push_back(i);
            }
            else{
                 if((zsxz[i]*zsxz[i]-threshold)>1e-6) sn_ids.push_back(i);
            }
            
        }
        
        
        
    }
   void est_cov_bxy(MatrixXd &covbxy, VectorXd &_zsxz,VectorXd &_bxy,VectorXd &_seyz,VectorXd &_bxz, MatrixXd &_LD_heidi)
    {
        int nsnp =_zsxz.size();
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
       
        dev.resize(nsnp-1);
        if(nsnp>1)
        {
            _bxy=_byz.array()/_bxz.array();
            _sexy=sqrt((_seyz.array()*_seyz.array()*_bxz.array()*_bxz.array()+_sexz.array()*_sexz.array()*_byz.array()*_byz.array())/(_bxz.array()*_bxz.array()*_bxz.array()*_bxz.array()));
            
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
            
            //unsigned int m=lambda.size();
            
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
        free(gdata->allele_1);
        free(gdata->allele_2);
        free(gdata->freq);
        free(gdata->byz);
        free(gdata->seyz);
        free(gdata->pvalue);
        free(gdata->splSize);
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
        
//#pragma omp parallel for
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

    void EstLD(bInfo* bdata, vector<int> &smpl, double wind_size, vector< vector<string> > &snp, vector< vector<double> > &r, vector<double> &r2, vector<double> &md_r2, vector<double> &max_r2, vector<string> &max_r2_snp, vector<double> &dL, vector<double> &dR, vector<int> &K, vector<string> &L_SNP, vector<string> &R_SNP, double alpha, bool IncldQ)
    {
        int i = 0, j = 0, L = 0, R = 0, maxL = 0, maxR = 0, i_buf = 0;
        map<int, int> smpl_snp_map;
        for (i = 0; i < smpl.size(); i++) smpl_snp_map.insert(pair<int, int>(bdata->_include[smpl[i]], i));
        
        cout << "Parameters used to search SNPs in LD with the given SNPs: window size=" << (int) (wind_size * 0.001) << "Kb, significant level=" << alpha << endl;
        vector<double> rst, y, x;
        for (i = 0; i < smpl.size(); i++) {
            vector<int> buf;
            vector<double> r_buf, rsq_buf;
            maxL = maxR = L = R = smpl[i];
            cout<<bdata->_snp_name[bdata->_include[R]]<<endl;
            makex(bdata,L, y);
            if (IncldQ) {
                buf.push_back(L);
                rsq_buf.push_back(1.0);
                r_buf.push_back(1.0);
            }
            while (1) {
                if (R == bdata->_include.size() - 1) break;
                if (bdata->_chr[bdata->_include[R]] != bdata->_chr[bdata->_include[R + 1]]) break;
                if (bdata->_bp[bdata->_include[R + 1]] - bdata->_bp[bdata->_include[smpl[i]]] > wind_size) break;
                R++;
                cout<<bdata->_snp_name[bdata->_include[R]]<<endl;
                if (smpl_snp_map.find(bdata->_include[R]) != smpl_snp_map.end()) continue;
                makex(bdata,R, x);
                reg(y, x, rst);
                if (rst[2] < alpha) {
                    maxR = R;
                    buf.push_back(R);
                    rsq_buf.push_back(rst[3]);
                    r_buf.push_back(rst[4]);
                }
            }
            while (1) {
                if (L == 0) break;
                if (bdata->_chr[bdata->_include[L]] != bdata->_chr[bdata->_include[L - 1]]) break;
                if (bdata->_bp[bdata->_include[smpl[i]]] - bdata->_bp[bdata->_include[L - 1]] > wind_size) break;
                L--;
                cout<<bdata->_snp_name[bdata->_include[L]]<<endl;
                if (smpl_snp_map.find(bdata->_include[L]) != smpl_snp_map.end()) continue;
                makex(bdata, L, x);
                reg(y, x, rst);
                if (rst[2] < alpha) {
                    maxL = L;
                    buf.insert(buf.begin(), L);
                    rsq_buf.insert(rsq_buf.begin(), rst[3]);
                    r_buf.insert(r_buf.begin(), rst[4]);
                }
            }
            if (buf.size() == 0) {
                K[i] = 0;
                dL[i] = 0;
                dR[i] = 0;
                L_SNP[i] = "NA";
                R_SNP[i] = "NA";
                r[i].push_back(0.0);
                snp[i].push_back("NA");
                r2[i] = 0.0;
                md_r2[i] = 0.0;
                max_r2[i] = 0.0;
                max_r2_snp[i] = "NA";
            } else {
                K[i] = buf.size();
                dL[i] = bdata->_bp[bdata->_include[smpl[i]]] - bdata->_bp[bdata->_include[maxL]];
                dR[i] = bdata->_bp[bdata->_include[maxR]] - bdata->_bp[bdata->_include[smpl[i]]];
                L_SNP[i] = bdata->_snp_name[bdata->_include[maxL]];
                R_SNP[i] = bdata->_snp_name[bdata->_include[maxR]];
                for (j = 0; j < K[i]; j++) {
                    r[i].push_back(r_buf[j]);
                    snp[i].push_back(bdata->_snp_name[bdata->_include[buf[j]]]);
                }
                r2[i] = CommFunc::mean(rsq_buf);
                md_r2[i] = CommFunc::median(rsq_buf);
                i_buf = max_element(rsq_buf.begin(), rsq_buf.end()) - rsq_buf.begin();
                max_r2[i] = rsq_buf[i_buf];
                max_r2_snp[i] = snp[i][i_buf];
            }
            cout << i + 1 << " of " << smpl.size() << " target SNPs.\r";
        }
    }
    // est LD values for one SNP versus the SNPs around it.  
    void LD_Blocks(bInfo* bdata, double wind_size, double alpha, bool IncldQ, vector<string> &_ld_target_snp)
    {
        int i = 0, j = 0;
        
        // Read snplist file
        vector<int> smpl_buf, smpl;
        vector<string> uni_snp;
        for (i = 0; i < bdata->_include.size(); i++) uni_snp.push_back(bdata->_snp_name[bdata->_include[i]]);
        StrFunc::match(_ld_target_snp, uni_snp, smpl_buf);
        for (i = 0; i < smpl_buf.size(); i++) {
            if (smpl_buf[i]>-1) smpl.push_back(smpl_buf[i]);
        }
        int SNP_SmplNum = smpl.size(); // smpl is the position of _include        
        // Calculate LD structure
        cout << "Estimating LD structure..." << endl;
        vector<int> K(SNP_SmplNum);
        vector<double> r2(SNP_SmplNum), md_r2(SNP_SmplNum), max_r2(SNP_SmplNum), dL(SNP_SmplNum), dR(SNP_SmplNum);
        vector<string> max_r2_snp(SNP_SmplNum);
        vector< vector<double> > r(SNP_SmplNum);
        vector<string> L_SNP(SNP_SmplNum), R_SNP(SNP_SmplNum);
        vector< vector<string> > snp_ls(SNP_SmplNum);
        
        
        EstLD(bdata, smpl, wind_size, snp_ls, r, r2, md_r2, max_r2, max_r2_snp, dL, dR, K, L_SNP, R_SNP, alpha, IncldQ);        
        
        
       
        cout << "target_SNP\tfreq\tL_region\tR_region\tL_snp\tR_snp\tnSNPs\tmean_rsq\tmedian_rsq\tmax_rsq\tmax_rsq_snp" << endl;
        for (i = 0; i < SNP_SmplNum; i++) cout << bdata->_snp_name[bdata->_include[smpl[i]]] << "\t" << "\t" << dL[i] << "\t" << dR[i] << "\t" << L_SNP[i] << "\t" << R_SNP[i] << "\t" << K[i] << "\t" << r2[i] << "\t" << md_r2[i] << "\t" << max_r2[i] << "\t" << max_r2_snp[i] << endl;
       
        
        for (i = 0; i < SNP_SmplNum; i++) {
            for (j = 0; j < r[i].size(); j++) cout << r[i][j] << " ";
            cout << endl;
        }
        
        for (i = 0; i < SNP_SmplNum; i++) {
            for (j = 0; j < snp_ls[i].size(); j++) cout << snp_ls[i][j] << " ";
            cout << endl;
        }       
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
   
    // SNP mode: each vector is the genotype of a SNP with all of individuals
    void make_XMat_SNPs(bInfo* bdata, vector< vector<float> > &X, bool miss_with_mu)
    {
        if(bdata->_mu.empty() && miss_with_mu) calcu_mu(bdata);
        
        cout<<"Recoding genotypes (SNP major mode) ..."<<endl;
        int i=0, j=0;
        X.clear();
        X.resize(bdata->_include.size());
        for(i=0; i<bdata->_include.size(); i++) X[i].resize(bdata->_keep.size());
        for(i=0; i<bdata->_keep.size(); i++){
            bool need2fill=false;
            if(bdata->_dosage_flag){
                for(j=0; j<bdata->_include.size(); j++){
                    if(bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]]<1e5){
                        if(bdata->_allele1[bdata->_include[j]]==bdata->_ref_A[bdata->_include[j]]) X[j][i]=bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]];
                        else X[j][i]=2.0-bdata->_geno_dose[bdata->_keep[i]][bdata->_include[j]];
                    }
                    else{ X[j][i]=1e6; need2fill=true; }
                }
                bdata->_geno_dose[i].clear();
            }
            else{
                for(j=0; j<bdata->_include.size(); j++){
                    if(!bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]] || bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]){
                        if(bdata->_allele1[bdata->_include[j]]==bdata->_ref_A[bdata->_include[j]]) X[j][i]=(bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]]+bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                        else X[j][i]=2.0-(bdata->_snp_1[bdata->_include[j]][bdata->_keep[i]]+bdata->_snp_2[bdata->_include[j]][bdata->_keep[i]]);
                    }
                    else{ X[j][i]=1e6; need2fill=true; }
                }
            }
            // Fill the missing genotype with the mean of x (2p)
            for(j=0; j<bdata->_include.size() && miss_with_mu && need2fill; j++){
                if(X[j][i]>1e5) X[j][i]=bdata->_mu[bdata->_include[j]];
            }
        }
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

	void extract_eqtl_snp(eqtlInfo* eqtlinfo, string snplstName)
	{
		vector<string> snplist;
		string msg = "SNPs";
		read_msglist(snplstName, snplist, msg);
		eqtlinfo->_esi_include.clear();
		StrFunc::match_only(snplist, eqtlinfo->_esi_rs, eqtlinfo->_esi_include);
		stable_sort(eqtlinfo->_esi_include.begin(), eqtlinfo->_esi_include.end());        
		cout << eqtlinfo->_esi_include.size() << " SNPs are extracted from [" + snplstName + "]." << endl;
	}
    
    void exclude_eqtl_snp(eqtlInfo* eqtlinfo, string snplstName)
    {
        vector<string> snplist;
        vector<string> mapstr;
        string msg = "SNPs";
        read_msglist(snplstName, snplist, msg);
        int pre_num=eqtlinfo->_esi_include.size();
        mapstr.resize(pre_num);
        for(int i=0;i<pre_num;i++)
            mapstr[i]=eqtlinfo->_esi_rs[eqtlinfo->_esi_include[i]];
        eqtlinfo->_esi_include.clear();
        StrFunc::set_complement(snplist, mapstr, eqtlinfo->_esi_include); //sorted
        stable_sort(eqtlinfo->_esi_include.begin(), eqtlinfo->_esi_include.end());       
        cout << pre_num-eqtlinfo->_esi_include.size() << " SNPs are excluded from [" + snplstName + "] and there are " << eqtlinfo->_esi_include.size() << " SNPs remaining." << endl;
    }
    
    void extract_prob(eqtlInfo* eqtlinfo,string problstName)
    {
        vector<string> problist;
        string msg="probes";
        read_msglist(problstName, problist,msg);        
        eqtlinfo->_include.clear();
		StrFunc::match_only(problist, eqtlinfo->_epi_prbID, eqtlinfo->_include);
        stable_sort(eqtlinfo->_include.begin(), eqtlinfo->_include.end());
        cout<<eqtlinfo->_include.size()<<" probes are extracted from ["+problstName+"]."<<endl;
    }
    void exclude_prob(eqtlInfo* eqtlinfo,string problstName)
    {
        vector<string> problist;
        vector<string> mappro;
        string msg="probes";
        read_msglist(problstName, problist,msg);
        int pre_num=eqtlinfo->_include.size();
        mappro.resize(pre_num);
        for(int i=0;i<pre_num;i++) mappro[i]=eqtlinfo->_epi_prbID[eqtlinfo->_include[i]];
        eqtlinfo->_include.clear();
        StrFunc::set_complement(problist, mappro, eqtlinfo->_include);//sorted
        //stable_sort(eqtlinfo->_include.begin(), eqtlinfo->_include.end());
        cout<<pre_num-eqtlinfo->_include.size()<<" probes are excluded from ["+problstName+"]and there are "<<eqtlinfo->_include.size()<<" probes remaining."<<endl;
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
    void filter_probe_null(eqtlInfo* eqtlinfo)
    {
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
            }
        }
        else{
            eqtlinfo->_include.clear();
            for (int i = 0; i < eqtlinfo->_probNum; i++)
            {
                if (eqtlinfo->_cols[(i<<1)+1] > eqtlinfo->_cols[i<<1]) eqtlinfo->_include.push_back(i);
            }
        }
        cout<<eqtlinfo->_include.size()<<" probes to be included."<<endl;
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
        }else
        {
            StrFunc::match_only(bdata->_snp_name, gdata->snpName, edId);
            if(edId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
            for(int i=0;i<edId.size();i++) cmmnSNPs.push_back(gdata->snpName[edId[i]]);
            edId.clear();
            StrFunc::match_only(cmmnSNPs, esdata->_esi_rs, edId);
            if(edId.empty()) throw("Error: no common SNPs found.");
        }
        for(int i=0;i<edId.size();i++) slctSNPs.push_back(esdata->_esi_rs[edId[i]]);
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
            char a1, a2, ga1, ga2, ea1, ea2;
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
                    
                    char tmpch=bdata->_ref_A[bdId[i]];
                    bdata->_ref_A[bdId[i]]=bdata->_other_A[bdId[i]];
                    bdata->_other_A[bdId[i]]=tmpch;
                   
                }
                else if(ea1 == ga2 && ea2 == ga1)
                {
                    bdata->_include.push_back(bdId[i]);
                    gdata->_include.push_back(gdId[i]);
                    esdata->_esi_include.push_back(edId[i]);
                    
                    gdata->byz[gdId[i]]=-gdata->byz[gdId[i]];
                    char tmpch=bdata->_ref_A[bdId[i]];
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
        }else
        {
            StrFunc::match_only(esdata->_esi_rs, gdata->snpName, gdId);
            if(gdId.empty()) throw("Error: no common SNPs found.");
        }
        for(int i=0;i<gdId.size();i++) slctSNPs.push_back(gdata->snpName[gdId[i]]);
        
        //alleles check
        StrFunc::match(slctSNPs, esdata->_esi_rs, edId);
        cmmnSNPs.clear();
        gdata->_include.clear();
        esdata->_esi_include.clear();
        
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
    void update_gwas(gwasData* gdata){
        
        
        vector<string> snpName(gdata->_include.size());
        char* allele_1=(char*)malloc((gdata->_include.size()+1)*sizeof(char));
        char* allele_2=(char*)malloc((gdata->_include.size()+1)*sizeof(char));
        double* freq=(double*)malloc(gdata->_include.size()*sizeof(double));
        double* byz=(double*)malloc(gdata->_include.size()*sizeof(double));
        double* seyz=(double*)malloc(gdata->_include.size()*sizeof(double));
        double* pvalue=(double*)malloc(gdata->_include.size()*sizeof(double));
        uint32_t* splSize=(uint32_t*)malloc(gdata->_include.size()*sizeof(uint32_t));
        
        gdata->snpNum=gdata->_include.size();
        for(int i=0;i<gdata->_include.size();i++ )
        {
            snpName[i]=gdata->snpName[gdata->_include[i]];
            allele_1[i]=gdata->allele_1[gdata->_include[i]];
            allele_2[i]=gdata->allele_2[gdata->_include[i]];
            freq[i]=gdata->freq[gdata->_include[i]];
            byz[i]=gdata->byz[gdata->_include[i]];
            seyz[i]=gdata->seyz[gdata->_include[i]];
        }
        
        free(gdata->allele_1);
        free(gdata->allele_2);
        free(gdata->freq);
        free(gdata->byz);
        free(gdata->seyz);
        free(gdata->pvalue);
        free(gdata->splSize);
        
        gdata->snpName=snpName;
        gdata->allele_1=allele_1;
        gdata->allele_2=allele_2;
        gdata->freq=freq;
        gdata->byz=byz;
        gdata->seyz=seyz;
        gdata->pvalue = pvalue;
        gdata->splSize=splSize;
             
        gdata->allele_1[gdata->_include.size()]='\0';
        gdata->allele_2[gdata->_include.size()]='\0';
        
        for(int i=0;i<gdata->snpNum;i++) gdata->_include[i]=i;
        
    }

    void make_XMat(bInfo* bdata,vector<uint32_t> &snpids, MatrixXd &X, bool minus_2p = false) {
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
        vector<uint32_t> tmpIdx1;
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

    void read_smaslist(vector<string> &smasNames, string eqtlsmaslstName)
    {
        ifstream smas(eqtlsmaslstName.c_str());
        if (!smas) throw ("Error: can not open the file [" + eqtlsmaslstName + "] to read.");
        cout << "Reading eQTL summary file names from [" + eqtlsmaslstName + "]." << endl;
        char buf[MAX_LINE_SIZE];
        while (smas.getline(buf, MAX_LINE_SIZE))
        {
            string tmpStr;
            istringstream iss(buf);
            iss >> tmpStr;
            smasNames.push_back(tmpStr);
        }
        cout << smasNames.size()<<" eQTL summary file names are included from [" + eqtlsmaslstName + "]." << endl;
        smas.close();
    }
    void combine_esi(eqtlInfo* eqtlinfo, vector<string> &smasNames)
    {
        eqtlinfo->_esi_chr.clear();
        eqtlinfo->_esi_rs.clear();
        eqtlinfo->_esi_gd.clear();
        eqtlinfo->_esi_bp.clear();
        eqtlinfo->_esi_allele1.clear();
        eqtlinfo->_esi_allele2.clear();
        eqtlinfo->_esi_include.clear();
        eqtlinfo->_snp_name_map.clear();
        long size = 0, counter = 0;
        
        for (int i = 0; i < smasNames.size(); i++)
        {
            string esifile = smasNames[i]+".esi";
            ifstream esi(esifile.c_str());
            if (!esi) throw ("Error: can not open the file [" + esifile + "] to read.");
            cout << "Reading eQTL SNP information from [" + esifile + "]." << endl;
            
            char buf[MAX_LINE_SIZE];
            int lineNum(0);
            while (!esi.eof())
            {
                esi.getline(buf, MAX_LINE_SIZE);
                lineNum++;
            }
            if (buf[0] == '\0') lineNum--;
            cout << lineNum << " SNPs to be included from [" + esifile + "]." << endl;
            
            esi.clear(ios::goodbit);
            esi.seekg(0, ios::beg);
            for (int j = 0; j<lineNum; j++)
            {
                string tmpStr;
                esi.getline(buf, MAX_LINE_SIZE);
                istringstream iss(buf);
                iss >> tmpStr;
                int chrtmp = atoi(tmpStr.c_str());
                iss >> tmpStr;
                
                eqtlinfo->_snp_name_map.insert(pair<string, int>(tmpStr,counter));
                if (size < eqtlinfo->_snp_name_map.size())
                {
                    size = eqtlinfo->_snp_name_map.size();
                    eqtlinfo->_esi_chr.push_back(chrtmp);
                    eqtlinfo->_esi_rs.push_back(tmpStr);
                    iss >> tmpStr;
                    eqtlinfo->_esi_gd.push_back(atoi(tmpStr.c_str()));
                    iss >> tmpStr;
                    eqtlinfo->_esi_bp.push_back(atoi(tmpStr.c_str()));
                    iss >> tmpStr;
                    StrFunc::to_upper(tmpStr);
                    eqtlinfo->_esi_allele1.push_back(tmpStr.c_str()[0]);
                    iss >> tmpStr;
                    StrFunc::to_upper(tmpStr);
                    eqtlinfo->_esi_allele2.push_back(tmpStr.c_str()[0]);
                    eqtlinfo->_esi_include.push_back(counter++);
                }
            }
            esi.close();
        }
        eqtlinfo->_snpNum =eqtlinfo->_esi_include.size();
        cout <<"Total " <<eqtlinfo->_snpNum << " SNPs to be included! " << endl;
    }
    
    void combine_epi_esd(eqtlInfo* eqtlinfo, vector<string> &smasNames, string outFileName)
    {
        eqtlinfo->_cols.push_back(0);
        vector<uint32_t> duplicatedid;
        for (int i = 0; i < smasNames.size(); i++)
        {
            string esifile = smasNames[i] + ".esi";
            string epifile = smasNames[i] + ".epi";
            string esdfile = smasNames[i] + ".besd";
            eqtlInfo tmp_info;
            read_esifile(&tmp_info, esifile);
            read_epifile(&tmp_info, epifile);
            read_besdfile(&tmp_info, esdfile);
            for (int j = 0; j < tmp_info._probNum; j++)
            {
                duplicatedid.clear();
                uint32_t beta_start = tmp_info._cols[j << 1];
                uint32_t se_start = tmp_info._cols[1 + (j << 1)];
                int numsnps = se_start - beta_start;
                for (int k = 0; k < numsnps; k++)
                {
                    uint32_t ge_rowid = tmp_info._rowid[beta_start + k];
                    string snpname = tmp_info._esi_rs[ge_rowid];
                    map<string, int>::iterator iter;
                    iter = eqtlinfo->_snp_name_map.find(snpname);
                    if (iter == eqtlinfo->_snp_name_map.end()) throw("Unknown error!");
                    uint32_t sid = iter->second;
                    eqtlinfo->_rowid.push_back(sid);
                    duplicatedid.push_back(sid);
                    eqtlinfo->_val.push_back(tmp_info._val[beta_start + k]);
                }
                for (int k = 0; k < numsnps; k++)
                {
                    uint32_t ge_rowid = tmp_info._rowid[se_start + k];
                    uint32_t sid = duplicatedid[k];
                    eqtlinfo->_rowid.push_back(sid);
                    eqtlinfo->_val.push_back(tmp_info._val[se_start + k]);
                }
                uint32_t lastcolnum = eqtlinfo->_cols[eqtlinfo->_cols.size() - 1] ;
                eqtlinfo->_cols.push_back(lastcolnum + numsnps);
                eqtlinfo->_cols.push_back(lastcolnum + (numsnps<<1));
                
                eqtlinfo->_epi_bp.push_back(tmp_info._epi_bp[j]);
                eqtlinfo->_epi_chr.push_back(tmp_info._epi_chr[j]);
                eqtlinfo->_epi_gd.push_back(tmp_info._epi_gd[j]);
                eqtlinfo->_epi_gene.push_back(tmp_info._epi_gene[j]);
                eqtlinfo->_epi_orien.push_back(tmp_info._epi_orien[j]);
                eqtlinfo->_epi_prbID.push_back(tmp_info._epi_prbID[j]);
            }
        }
        eqtlinfo->_valNum = eqtlinfo->_val.size();
        eqtlinfo->_probNum = eqtlinfo->_epi_prbID.size();
        for (int i = 0; i < eqtlinfo->_probNum; i++) eqtlinfo->_include.push_back(i);
        cout << "Total "<<eqtlinfo->_probNum << " probes and "<< eqtlinfo->_snpNum << " SNPs have been combined." << endl;
        
        filter_probe_null(eqtlinfo); // at the same time, reset the vector _include
        cout << "\nsaving eQTL data..." << endl;
        string esdfile = outFileName + ".esi";
        ofstream smr(esdfile.c_str());
        if (!smr) throw ("Error: can not open the esi file " + esdfile + " to save!");
        for (int i = 0; i <eqtlinfo->_snpNum; i++) {
            smr << eqtlinfo->_esi_chr[i] << '\t' << eqtlinfo->_esi_rs[i] << '\t' << eqtlinfo->_esi_gd[i] << '\t' << eqtlinfo->_esi_bp[i] << '\t' << eqtlinfo->_esi_allele1[i] << '\t' << eqtlinfo->_esi_allele2[i] << '\n';
        }
        smr.close();
        cout << eqtlinfo->_snpNum << " SNPs have been saved in the file [" + esdfile + "]." << endl;
        
        esdfile = string(outFileName) + ".epi";
        smr.open(esdfile.c_str());
        if (!smr) throw ("Error: can not open the epi file " + esdfile + " to save!");
        for (int i = 0; i <eqtlinfo->_include.size(); i++) {
            smr << eqtlinfo->_epi_chr[eqtlinfo->_include[i]] << '\t' << eqtlinfo->_epi_prbID[eqtlinfo->_include[i]] << '\t' << eqtlinfo->_epi_gd[eqtlinfo->_include[i]] << '\t' << eqtlinfo->_epi_bp[eqtlinfo->_include[i]] << '\t' << eqtlinfo->_epi_gene[eqtlinfo->_include[i]] << '\t' << eqtlinfo->_epi_orien[eqtlinfo->_include[i]] << '\n';
        }
        smr.close();
        cout << eqtlinfo->_probNum << " probes have been saved in the file [" + esdfile + "]." << endl;
        
        esdfile = string(outFileName) + ".besd";
        FILE * smrbesd;
        smrbesd = fopen(esdfile.c_str(), "wb");
        uint64_t colSize = sizeof(uint32_t)*((eqtlinfo->_include.size() << 1) + 1);
        uint64_t rowSize = sizeof(uint32_t)*eqtlinfo->_valNum;
        uint64_t valSize = sizeof(float)*eqtlinfo->_valNum;
        uint64_t valNum = eqtlinfo->_valNum;
        uint64_t bufsize = sizeof(float) + sizeof(uint64_t) + colSize + rowSize + valSize;
        
        char* buffer = (char*)malloc(sizeof(char)*bufsize);
        memset(buffer, 0, sizeof(char)*bufsize);
        float ftype = SPARSE_FILE_TYPE_2;
        memcpy(buffer, &ftype, sizeof(float));
        char* wptr = buffer + sizeof(float);
        memcpy(wptr, &valNum, sizeof(uint64_t));
        wptr += sizeof(uint64_t);
        uint32_t* uptr = (uint32_t*)wptr; *uptr++ = 0;
        for (int i = 0; i<eqtlinfo->_include.size(); i++)
        {
            *uptr++ = eqtlinfo->_cols[(eqtlinfo->_include[i] << 1) + 1];
            *uptr++ = eqtlinfo->_cols[eqtlinfo->_include[i] + 1 << 1];
        }
        wptr += colSize;
        memcpy(wptr, &eqtlinfo->_rowid[0], rowSize);
        wptr += rowSize;
        memcpy(wptr, &eqtlinfo->_val[0], valSize);
        fwrite(buffer, sizeof(char), bufsize, smrbesd);
        free(buffer);
        fclose(smrbesd);
        
        cout << "Beta values and SE values for " << eqtlinfo->_include.size() << " Probes and " << eqtlinfo->_snpNum << " SNPs have been saved in the binary file [" + esdfile + "]." << endl;
        
    }

    
    
    void smr(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero , char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr, char* refSNP, bool heidioffFlag)
    {
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
        if(bFlag) read_besdfile(&esdata, string(eqtlFileName)+".besd");
        else      read_esdfile(&esdata, string(eqtlFileName)+".esd");
        //test the SNP alignement
        //for(int i=0;i<bdata._include.size();i++)
        //    if(gdata.snpName[i]!=esdata._esi_rs[i] || bdata._snp_name[bdata._include[i]] != gdata.snpName[i] || bdata._snp_name[bdata._include[i]] !=esdata._esi_rs[i] ) cout<<i<<endl;
        
        // till now. SNPs in the three dataset are aligned in bfile order. Alleles are checked by eQTL data as the reference
        // in sparse besd, rowid may be out-of-order
        // bdata._include is the SNPs id to test in case of MAF pruning and HWE pruning being conducted.
        // use bdata._include to extract genotype data for LD
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
        vector<char> rsa1(probNum);
        vector<char> rsa2(probNum);
        vector<string> rsbp(probNum); //origin is unsigned int
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
        
        vector<char> allele1;
        vector<char> allele2;
        vector<uint32_t> bpsnp;
        
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
        MatrixXd _LD_heidi;
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
            byz.clear();
            seyz.clear();
            allele1.clear();
            allele2.clear();
            bpsnp.clear();
            long maxid =-9;
            if(esdata._rowid.empty())
            {
                for (int j = 0; j<bdata._include.size(); j++)// bdata._include.size() == esdata._esi_include.size() == gdata._include.size()
                {
                    if (abs(esdata._bxz[i][j] + 9) > 1e-6)
                    {
                        bxz.push_back(esdata._bxz[i][j]);
                        sexz.push_back(esdata._sexz[i][j]);
                        byz.push_back(gdata.byz[j]);
                        seyz.push_back(gdata.seyz[j]);
                        curId.push_back(j);
                        eName.push_back(esdata._esi_rs[j]);
                        if(heidiFlag && esdata._esi_rs[j]==string(refSNP)) maxid=(eName.size()-1);
                        if(!heidioffFlag) bpsnp.push_back(bdata._bp[bdata._include[j]]);
                        allele1.push_back(esdata._esi_allele1[j]);
                        allele2.push_back(esdata._esi_allele2[j]);
                    }
                }
                
            }
            else{
                int beta_start=esdata._cols[i<<1];
                int se_start=esdata._cols[1+(i<<1)];
                int numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    int ge_rowid=esdata._rowid[beta_start+j];
                    bxz.push_back(esdata._val[beta_start+j]);
                    sexz.push_back(esdata._val[se_start+j]);
                    byz.push_back(gdata.byz[ge_rowid]);
                    seyz.push_back(gdata.seyz[ge_rowid]);
                    curId.push_back(ge_rowid);
                    eName.push_back(esdata._esi_rs[ge_rowid]);
                    if(heidiFlag && esdata._esi_rs[ge_rowid]==string(refSNP)) maxid=(eName.size()-1);
                    allele1.push_back(esdata._esi_allele1[ge_rowid]);
                    allele2.push_back(esdata._esi_allele2[ge_rowid]);
                    if(!heidioffFlag) bpsnp.push_back(bdata._bp[bdata._include[ge_rowid]]);
                   
                }
            }
            if(heidiFlag && maxid==-9) continue;
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
            if(!heidioffFlag) rsbp[outCount]=itos(bpsnp[maxid]);
            rsa1[outCount]=allele1[maxid];
            rsa2[outCount]=allele2[maxid];
            
            if(!heidioffFlag)
            {
                //extract info from reference
                make_XMat(&bdata,curId, _X); //_X: one row one individual, one column one SNP
                
                //last vesion ref_snpData was used. row of ref_snpData is SNP, column of ref_snpData is individual
                cor_calc(_LD, _X);
                
                sn_ids.clear(); //increase order
                if(abs(ld_top-1)<1e-6) get_square_idxes(sn_ids,zsxz,threshold);
                else get_square_ldpruning_idxes(sn_ids,zsxz,threshold,_LD, maxid,ld_top);
                
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
                _LD_heidi.resize(sn_ids.size(),sn_ids.size());
                
                #pragma omp parallel for
                for(int j=0;j<sn_ids.size();j++)
                {
                    _byz[j]=byz[sn_ids[j]];
                    _seyz[j]=seyz[sn_ids[j]];
                    _bxz[j]=bxz[sn_ids[j]];
                    _sexz[j]=sexz[sn_ids[j]];
                    _zsxz[j]=zsxz[sn_ids[j]];
                    for(int k=0;k<=j;k++)_LD_heidi(j,k)=_LD_heidi(k,j)=_LD(sn_ids[j],sn_ids[k]);
                }
                
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
        
        string smrfile = string(outFileName)+".smr";
        ofstream smr(smrfile.c_str());
        if (!smr) throw ("Error: can not open the fam file " + smrfile + " to save!");
        if(heidioffFlag)
        {
            smr << "ProbeID" <<'\t'<< "Chr" <<'\t' << "Gene"  << '\t' << "Prob_bp" << '\t'<< "SNP"<< '\t' << "A1"<< '\t'<< "A2"<< '\t'<<"b_GWAS"<<'\t'<<"se_GWAS"<<'\t'<< "p_GWAS" << '\t'<<"b_eQTL"<<'\t'<<"se_eQTL"<<'\t'<< "p_eQTL" << '\t'<< "b_SMR" << '\t'<< "se_SMR"<< '\t' << "p_SMR" << '\n';
            
            for (int i = 0;i <=outCount; i++) {
                smr<<esdata._epi_prbID[out_probid[i]]<<'\t'<<esdata._epi_chr[out_probid[i]]<<'\t'<<esdata._epi_gene[out_probid[i]]<<'\t'<<esdata._epi_bp[out_probid[i]]<<'\t'<<rsid[i]<<'\t'<<rsa1[i]<<'\t'<<rsa2[i]<<'\t'<<bgwas[i]<<'\t'<<segwas[i]<<'\t'<<pgwas[i]<<'\t'<<beqtl[i]<<'\t'<<seeqtl[i]<<'\t'<<peqtl[i]<<'\t'<<bxy[i]<<'\t'<<sexy[i]<<'\t'<<pxy[i]<<'\n';
            }
            cout<<"SMR analysis finished.\nSMR analysis results of "<<outCount+1<<" probes have been saved in the file [" + smrfile + "]."<<endl;
        }
        else
        {
            smr << "ProbeID" <<'\t'<< "Chr" <<'\t' << "Gene"  << '\t' << "Prob_bp" << '\t'<< "SNP"<< '\t' << "SNP_bp"<< '\t' << "A1"<< '\t'<< "A2"<< '\t'<<"b_GWAS"<<'\t'<<"se_GWAS"<<'\t'<< "p_GWAS" << '\t'<<"b_eQTL"<<'\t'<<"se_eQTL"<<'\t'<< "p_eQTL" << '\t'<< "b_SMR" << '\t'<< "se_SMR"<< '\t' << "p_SMR" << "\t"<< "p_HET"<< "\t" << "nsnp" << '\n';
            
            for (int i = 0;i <=outCount; i++) {
                smr<<esdata._epi_prbID[out_probid[i]]<<'\t'<<esdata._epi_chr[out_probid[i]]<<'\t'<<esdata._epi_gene[out_probid[i]]<<'\t'<<esdata._epi_bp[out_probid[i]]<<'\t'<<rsid[i]<<'\t'<<rsbp[i]<<'\t'<<rsa1[i]<<'\t'<<rsa2[i]<<'\t'<<bgwas[i]<<'\t'<<segwas[i]<<'\t'<<pgwas[i]<<'\t'<<beqtl[i]<<'\t'<<seeqtl[i]<<'\t'<<peqtl[i]<<'\t'<<bxy[i]<<'\t'<<sexy[i]<<'\t'<<pxy[i]<<'\t'<<prb1[i]<<'\t'<<nsnp_test1[i]<<'\n';
            }
             cout<<"SMR and heterogeneity analysis finished.\nSMR and heterogeneity analysis results of "<<outCount+1<<" probes have been saved in the file [" + smrfile + "]."<<endl;
        }
        smr.close();
       
        free_gwas_data( &gdata);
        
    }
    
    void make_esd_file(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf, char* indilstName, char* snplstName,char* problstName,bool bFlag,bool make_besd_flag,bool make_esd_flag, char* indilst2remove, char* snplst2exclde, char* problst2exclde, bool cis_flag, int cis, float transThres, float restThres)
    {
        bInfo bdata;
        gwasData gdata;
        eqtlInfo eqtlinfo;
        
        if(bFileName != NULL)
        {
            read_famfile(&bdata, string(bFileName)+".fam");
            if(indilstName != NULL) keep_indi(&bdata,indilstName);
            if(indilst2remove != NULL) remove_indi(&bdata, indilst2remove);
            read_bimfile(&bdata, string(bFileName)+".bim");
            if(snplstName != NULL) extract_snp(&bdata, snplstName);
            if(snplst2exclde != NULL) exclude_snp(&bdata, snplst2exclde);
            read_bedfile(&bdata, string(bFileName)+".bed");
            if (bdata._mu.empty()) calcu_mu(&bdata);
            if(maf>0) filter_snp_maf(&bdata,maf);
            cout<<endl;
        }
        
        if(gwasFileName != NULL) read_gwas_data( &gdata, gwasFileName);
        
        
        cout<<endl<<"Reading eQTL summary data..."<<endl;
        if(eqtlFileName != NULL)
        {
            read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&eqtlinfo, snplstName);
            if(snplst2exclde != NULL) exclude_eqtl_snp(&eqtlinfo, snplst2exclde); //no switch the place ahead of extract_eqtl_snp()
            read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
            if(problstName != NULL) extract_prob(&eqtlinfo, problstName);
            if(problst2exclde != NULL) exclude_prob(&eqtlinfo, problst2exclde); //no switch the place ahead of extract_prob()
            if(bFlag) read_besdfile(&eqtlinfo, string(eqtlFileName)+".besd");
            else      read_esdfile(&eqtlinfo, string(eqtlFileName)+".esd");
            
        }
        else throw ("Error: please input the eQTL summary information for the eQTL data files by the option --beqtl-summary.");
        
        // get common SNPs
        vector<string> slctSNPs;
        vector<char> A1,A2;
        slctSNPs.clear();
        if(bFileName != NULL && gwasFileName != NULL)
        {
            cout<<endl<<"extracting common SNPs among PLink data,GWAS data and eQTL data ..."<<endl;
            vector<int> cmmnId;
            vector<string> cmmnSNPs;
            vector<string> b_snp_name;
            cmmnSNPs.clear();
            cmmnId.clear();
            if(maf>0)
            {
                b_snp_name.resize(bdata._include.size());
                for(int i=0;i<bdata._include.size();i++)
                    b_snp_name[i]=bdata._snp_name[bdata._include[i]];
                StrFunc::match_only(b_snp_name, gdata.snpName, cmmnId);
            }
            else
                StrFunc::match_only(bdata._snp_name, gdata.snpName, cmmnId);
            if(cmmnId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
            for(int i=0;i<cmmnId.size();i++) cmmnSNPs.push_back(gdata.snpName[cmmnId[i]]);
            //allele check
            vector<int> bdId;
            StrFunc::match(cmmnSNPs, bdata._snp_name, bdId);
            for (int i = 0; i<cmmnSNPs.size(); i++)
            {
                char a1, a2, ga1, ga2;
                a1 = bdata._allele1[bdId[i]];
                a2 = bdata._allele2[bdId[i]];
                ga1 = gdata.allele_1[cmmnId[i]];
                ga2 = gdata.allele_2[cmmnId[i]];
                if (a1 != '0' && a2 != '0' && ga1 != '0' && ga2 != '0' )
                    if ((a1 == ga1 && a2 == ga2) || (a1 == ga2 && a2 == ga1) )
                    {
                        slctSNPs.push_back(cmmnSNPs[i]);
                        A1.push_back(ga1);
                        A2.push_back(ga2);
                    }
            }
            if(slctSNPs.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
        }
        else if(bFileName != NULL && gwasFileName == NULL)
        {
            cout<<endl<<"extracting common SNPs between PLink data and eQTL data ..."<<endl;
            
            if(maf>0)
            {
                slctSNPs.resize(bdata._include.size());
                A1.resize(bdata._include.size());
                A2.resize(bdata._include.size());
                for(int i=0;i<bdata._include.size();i++)
                {
                    slctSNPs[i]=bdata._snp_name[bdata._include[i]];
                    A1[i]=bdata._allele1[bdata._include[i]];
                    A2[i]=bdata._allele2[bdata._include[i]];
                }
                
            }
            else
            {
                slctSNPs.assign(bdata._snp_name.begin(),bdata._snp_name.end());
                A1.assign(bdata._allele1.begin(),bdata._allele1.end());
                A2.assign(bdata._allele2.begin(),bdata._allele2.end());
            }
        }
        else if(bFileName == NULL && gwasFileName != NULL)
        {
            cout<<endl<<"extracting common SNPs between GWAS data and eQTL data ..."<<endl;
            slctSNPs.assign(gdata.snpName.begin(),gdata.snpName.end());
            uint32_t length=strlen(gdata.allele_1);
            A1.assign(gdata.allele_1,gdata.allele_1+length);
            A2.assign(gdata.allele_2,gdata.allele_2+length);
        }
        
        if(slctSNPs.empty())
        {
            //output (no Plink data or GWAS data input)
            if(cis_flag)
            {
                cout<<"Transforming dense file to sparse file ..."<<endl;
                if(eqtlinfo._valNum>0) throw("Error: please input dense format eQTL summary data file.");
                //if got flag (--extract-snp), filter_probe_null() should be invoked, then _include.
                vector<uint32_t> cols((eqtlinfo._probNum<<1)+1);
                vector<uint32_t> rowids;
                vector<float> val;
                FILE* logfile=NULL;
                char cCurrentPath[FILENAME_MAX];
                
                if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
                {
                   throw ("Error: can not get current directory!");
                }
                
                cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */

                
                string fname=string(cCurrentPath)+"/smr.log";
                logfile=fopen(fname.c_str(), "w");
                if (!(logfile)) {
                    sprintf("Error: Failed to open %s.\n", fname.c_str());
                    exit(1);
                }
                fputs(string("--cis-itvl:\t"+itos(cis)+"Mb\n").c_str(),logfile);
                 fputs(string("--trans-thres:\t"+dtos(transThres)+"\n").c_str(),logfile);
                 fputs(string("--rest-thres:\t"+dtos(restThres)+"\n").c_str(),logfile);
                 fputs(string("-----------------------------------------\n").c_str(),logfile);
                
                fputs("ProbeID\tProbeChr\tType\tStartBP\tEndBP\tnsnp\n",logfile);
                cis=cis*1e6;
                cols[0]=0;
                //in case of too many values over transThres in the trans region
                string pretransPrbId="";
                for(uint32_t i=0;i<eqtlinfo._probNum;i++)
                {
                    vector<int> esi_include;
                    int probchr=eqtlinfo._epi_chr[i];
                    int probbp= eqtlinfo._epi_bp[i];
                    uint32_t uperBounder=probbp+cis;
                    uint32_t lowerBounder=((probbp-cis>0)?(probbp-cis):0);
                    uint32_t cisNum=0;
                    uint32_t otherSlctNum=0;
                    for(int j=0;j<eqtlinfo._snpNum;j++)
                    {
                        double zsxz=eqtlinfo._bxz[i][j]/eqtlinfo._sexz[i][j];
                        double pxz=pchisq(zsxz*zsxz, 1);
                        if(eqtlinfo._esi_chr[j] == probchr && eqtlinfo._esi_bp[j]<=uperBounder && eqtlinfo._esi_bp[j]>=lowerBounder && abs(eqtlinfo._sexz[i][j]+9)>1e-6)
                        {
                            esi_include.push_back(j);
                            cisNum++;
                        }
                       
                        else if(pxz<transThres)
                        {
                            uint32_t transNum=0;
                            string logstr="";
                            esi_include.push_back(j);
                            long transbp=eqtlinfo._esi_bp[j];
                            long translowerBounder=((transbp-cis>0)?(transbp-cis):0);
                            long transuperBounder=transbp+cis;
                            //if esi is not sorted.                            
                           // for(int k=0;k<eqtlinfo._snpNum;k++)
                           //     if(k!=j && eqtlinfo._esi_chr[j] == eqtlinfo._esi_chr[k]  && abs(transbp-eqtlinfo._esi_bp[k])<=cis)   esi_include.push_back(k);
                            
                            // if esi is sorted
                            int startptr=j-1;
                            while(startptr>=0 && eqtlinfo._esi_chr[j] == eqtlinfo._esi_chr[startptr] && transbp-eqtlinfo._esi_bp[startptr]<=cis)
                            {
                                translowerBounder=eqtlinfo._esi_bp[startptr];
                                esi_include.push_back(startptr);
                                startptr--;
                                transNum++;
                                
                            }
                            startptr=j+1;
							while (startptr<eqtlinfo._snpNum && eqtlinfo._esi_chr[j] == eqtlinfo._esi_chr[startptr] && eqtlinfo._esi_bp[startptr] - transbp <= cis)
                            {
                                transuperBounder=eqtlinfo._esi_bp[startptr];
                                esi_include.push_back(startptr);
                                startptr++;
                                transNum++;
                                
                                //j=startptr-1; wrong!!! is there are trans upper, then the trans region should extend.
                            }
                            // for log
                            if( eqtlinfo._epi_prbID[i] != pretransPrbId )
                            {
                                logstr=eqtlinfo._epi_prbID[i]+"\t"+itos(eqtlinfo._epi_chr[i])+"\ttrans"+"\t"+ltos(translowerBounder)+"\t"+ltos(transuperBounder)+"\t"+itos(transNum)+"\n";
                                fputs(logstr.c_str(),logfile);
                                pretransPrbId=eqtlinfo._epi_prbID[i];
                            }                           
                        }
                        else if(pxz<restThres)
                        {
                            esi_include.push_back(j);
                            otherSlctNum++;
                        }
                        
                    }
                    if(cisNum)
                    {
                        string logstr=eqtlinfo._epi_prbID[i]+"\t"+itos(eqtlinfo._epi_chr[i])+"\tcis"+"\t"+ltos(lowerBounder)+"\t"+ltos(uperBounder)+"\t"+itos(cisNum)+"\n";
                        fputs(logstr.c_str(),logfile);
                    }
                    sort( esi_include.begin(), esi_include.end() );
                    vector<int> ::iterator it=unique(esi_include.begin(),esi_include.end());
                    esi_include.erase( it, esi_include.end() );
                    
                    for(int j=0;j<esi_include.size();j++)
                    {
                        val.push_back(eqtlinfo._bxz[i][esi_include[j]]);
                        rowids.push_back(esi_include[j]);
                    }
                    for(int j=0;j<esi_include.size();j++)
                    {
                        val.push_back(eqtlinfo._sexz[i][esi_include[j]]);
                        rowids.push_back(esi_include[j]);
                    }
                    uint64_t real_num=esi_include.size();
                    cols[(i<<1)+1]=real_num+cols[i<<1];
                    cols[i+1<<1]=(real_num<<1)+cols[i<<1];
                }
                
                uint64_t colSize=sizeof(uint32_t)*cols.size();
                uint64_t rowSize=sizeof(uint32_t)*rowids.size();
                uint64_t valSize=sizeof(float)*val.size();
                uint64_t valNum=val.size();
                uint64_t bufsize=sizeof(float)+sizeof(uint64_t)+colSize+rowSize+valSize;
                char* buffer=(char*)malloc (bufsize*sizeof(char));
                char* wptr=buffer;
                float ftype=SPARSE_FILE_TYPE_2;
                memcpy(wptr,&ftype,sizeof(float));
                wptr+=sizeof(float);
                memcpy(wptr,&valNum,sizeof(uint64_t));
                wptr+=sizeof(uint64_t);
                memcpy(wptr,&cols[0],colSize);
                wptr+=colSize;
                memcpy(wptr,&rowids[0],rowSize);
                wptr+=rowSize;
                memcpy(wptr,&val[0],valSize);
                
                cout<<"\nsaving eQTL data..."<<endl;
                string esdfile = string(outFileName)+".esi";
                ofstream smr(esdfile.c_str());
                if (!smr) throw ("Error: can not open the esi file " + esdfile + " to save!");
                for (int i = 0;i <eqtlinfo._snpNum; i++) {
                    smr<<eqtlinfo._esi_chr[i]<<'\t'<<eqtlinfo._esi_rs[i]<<'\t'<<eqtlinfo._esi_gd[i]<<'\t'<<eqtlinfo._esi_bp[i]<<'\t'<<eqtlinfo._esi_allele1[i]<<'\t'<<eqtlinfo._esi_allele2[i]<<'\n';
                }
                smr.close();
                cout<<eqtlinfo._snpNum<<" SNPs have been saved in the file [" + esdfile + "]."<<endl;
                
                esdfile = string(outFileName)+".epi";
                smr.open(esdfile.c_str());
                if (!smr) throw ("Error: can not open the epi file " + esdfile + " to save!");
                for (int i = 0;i <eqtlinfo._probNum; i++) {
                    smr<<eqtlinfo._epi_chr[i]<<'\t'<<eqtlinfo._epi_prbID[i]<<'\t'<<eqtlinfo._epi_gd[i]<<'\t'<<eqtlinfo._epi_bp[i]<<'\t'<<eqtlinfo._epi_gene[i]<<'\t'<<eqtlinfo._epi_orien[i]<<'\n';
                }
                smr.close();
                cout<<eqtlinfo._probNum<<" probes have been saved in the file [" + esdfile + "]."<<endl;
                
                esdfile = string(outFileName)+".besd";
                FILE * besd;
                besd = fopen (esdfile.c_str(), "wb");
                fwrite (buffer,sizeof(char), bufsize, besd);
                fclose(besd);
                if(logfile) fclose(logfile);
                free(buffer);
                cout<<"Beta values and SE values for "<<eqtlinfo._include.size()<<" Probes and "<<eqtlinfo._snpNum<<" SNPs have been saved in the binary file [" + esdfile + "]." <<endl;
                cout<<"Log information has been saved in the file [" + fname + "]." <<endl;
            }
            else
            {
                filter_probe_null(&eqtlinfo); // at the same time, reset the vector _include
                cout<<"\nsaving eQTL data..."<<endl;
                string esdfile = string(outFileName)+".esi";
                ofstream smr(esdfile.c_str());
                if (!smr) throw ("Error: can not open the esi file " + esdfile + " to save!");
                for (int i = 0;i <eqtlinfo._snpNum; i++) {
                    smr<<eqtlinfo._esi_chr[i]<<'\t'<<eqtlinfo._esi_rs[i]<<'\t'<<eqtlinfo._esi_gd[i]<<'\t'<<eqtlinfo._esi_bp[i]<<'\t'<<eqtlinfo._esi_allele1[i]<<'\t'<<eqtlinfo._esi_allele2[i]<<'\n';
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
                if(make_esd_flag)
                {
                    esdfile = string(outFileName)+".esd";
                    smr.open(esdfile.c_str());
                    if (!smr) throw ("Error: can not open the esd file " + esdfile + " to save!");
                    if(eqtlinfo._valNum==0)
                    {
                        for (int i = 0;i <eqtlinfo._snpNum; i++) {
                            for(int j=0;j<eqtlinfo._include.size();j++)
                            {
                                double v_tmp=eqtlinfo._sexz[eqtlinfo._include[j]][i];
                                if(abs(v_tmp+9)<1e-6) smr<<"NA"<<'\t'<<"NA"<<'\t';
                                else smr<<eqtlinfo._bxz[eqtlinfo._include[j]][i]<<'\t'<<eqtlinfo._sexz[eqtlinfo._include[j]][i]<<'\t';
                            }
                            smr<<'\n';
                        }
                    }
                    else
                    {   // if output transposed file, should be more easy.
                        eqtlinfo._bxz.resize(eqtlinfo._include.size());
                        eqtlinfo._sexz.resize(eqtlinfo._include.size());
                        for( uint32_t i=0;i<eqtlinfo._include.size();i++)
                        {
                            eqtlinfo._bxz[i].reserve(eqtlinfo._esi_include.size());
                            eqtlinfo._sexz[i].reserve(eqtlinfo._esi_include.size());
                        }
                        for(uint32_t i=0;i<eqtlinfo._include.size();i++)
                            for(uint32_t j=0;j<eqtlinfo._esi_include.size();j++)
                            {
                                eqtlinfo._bxz[i][j]=-9;
                                eqtlinfo._sexz[i][j]=-9;
                            }
                        for(uint32_t i=0;i<eqtlinfo._include.size();i++)
                        {
                            int proid=eqtlinfo._include[i];
                            int pos=eqtlinfo._cols[proid<<1];
                            int pos1=eqtlinfo._cols[(proid<<1)+1];
                            int num=pos1-pos;
                            for(int j=0;j<num;j++)
                            {
                                eqtlinfo._bxz[i][eqtlinfo._rowid[pos+j]]=eqtlinfo._val[pos+j];
                                eqtlinfo._sexz[i][eqtlinfo._rowid[pos+j]]=eqtlinfo._val[pos+j+num];
                            }
                        }
                        for (uint32_t i = 0;i <eqtlinfo._snpNum; i++) {
                            for(uint32_t j=0;j<eqtlinfo._include.size();j++)
                            {
                                double v_tmp=eqtlinfo._sexz[j][i];
                                if(abs(v_tmp+9)<1e-6) smr<<"NA"<<'\t'<<"NA"<<'\t';
                                else smr<<eqtlinfo._bxz[j][i]<<'\t'<<eqtlinfo._sexz[j][i]<<'\t';
                            }
                            smr<<'\n';
                        }
                    }
                    smr.close();
                    cout<<"Beta values and SE values for "<<eqtlinfo._include.size()<<" Probes and "<<eqtlinfo._snpNum<<" SNPs have been saved in the file [" + esdfile + "]." <<endl;
                    /*
                     esdfile = string(outFileName)+".tesd";
                     smr.open(esdfile.c_str());
                     if (!smr) throw ("Error: can not open the fam file " + esdfile + " to save!");
                     for (int i = 0;i <eqtlinfo._include.size(); i++) {
                     for(int j=0;j<eqtlinfo._snpNum;j++)
                     {
                     double v_tmp=eqtlinfo._sexz[eqtlinfo._include[i]][j];
                     if(abs(v_tmp+9)<1e-6) smr<<"NA"<<'\t';
                     else smr<<eqtlinfo._bxz[eqtlinfo._include[i]][j]<<'\t';
                     }
                     smr<<'\n';
                     for(int j=0;j<eqtlinfo._snpNum;j++)
                     {
                     double v_tmp=eqtlinfo._sexz[eqtlinfo._include[i]][j];
                     if(abs(v_tmp+9)<1e-6) smr<<"NA"<<'\t';
                     else smr<<eqtlinfo._sexz[eqtlinfo._include[i]][j]<<'\t';
                     }
                     smr<<'\n';
                     }
                     smr.close();
                     cout<<"Beta values and SE values for "<<eqtlinfo._include.size()<<" Probes and "<<eqtlinfo._snpNum<<" SNPs have been saved in the transoped file [" + esdfile + "]." <<endl;
                     */
                }
                
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
                        uint64_t colSize=sizeof(uint32_t)*((eqtlinfo._include.size()<<1)+1);
                        uint64_t rowSize=sizeof(uint32_t)*eqtlinfo._valNum;
                        uint64_t valSize=sizeof(float)*eqtlinfo._valNum;
                        uint64_t valNum=eqtlinfo._valNum;
                        uint64_t bufsize=sizeof(float)+sizeof(uint64_t)+colSize+rowSize+valSize;
                        
                        char* buffer=(char*)malloc (sizeof(char)*bufsize);
                        memset(buffer,0,sizeof(char)*bufsize);
                        float ftype=SPARSE_FILE_TYPE_2;
                        memcpy(buffer,&ftype,sizeof(float));
                        char* wptr=buffer+sizeof(float);
                        memcpy(wptr,&valNum,sizeof(uint64_t));
                        wptr+=sizeof(uint64_t);
                        uint32_t* uptr=(uint32_t*)wptr; *uptr++=0;
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
        }
        else
        {
            vector<int> cmmnId, slctId;
            vector<string> cmmnSNPs;
            cmmnSNPs.clear();
            cmmnId.clear();
            StrFunc::match_only(slctSNPs, eqtlinfo._esi_rs, cmmnId);
            if(cmmnId.empty()) throw("Error: no common SNPs.");
            for(int i=0;i<cmmnId.size();i++) cmmnSNPs.push_back(eqtlinfo._esi_rs[cmmnId[i]]);
            
            vector<int> slId;
            StrFunc::match(cmmnSNPs, slctSNPs, slId);
            slctSNPs.clear();
            for (int i = 0; i<cmmnSNPs.size(); i++)
            {
                char a1, a2, ea1, ea2;
                a1 = A1[slId[i]];
                a2 = A2[slId[i]];
                ea1 = eqtlinfo._esi_allele1[cmmnId[i]];
                ea2 = eqtlinfo._esi_allele2[cmmnId[i]];
                if (a1 != '0' && a2 != '0' && ea1 != '0' && ea2 != '0')
                    if ((a1 == ea1 && a2 == ea2) || (a1 == ea2 && a2 == ea1) || (a1 != ea1 && a1 != ea2 && a2 != ea2 && a2 != ea1))
                    {
                        slctSNPs.push_back(cmmnSNPs[i]);
                        slctId.push_back(cmmnId[i]);  // something like _esi_include
                    }
                
            }
            if(slctSNPs.empty()) throw("Error: no common SNPs after allele check.");
            cout<<slctSNPs.size()<<" SNPs in common to be included after allele check."<<endl;
            
            cout<<"\nfiltering out the probes with no value..."<<endl; // should be invoked always here
            vector<int> include; // something like _include
            if(eqtlinfo._valNum==0)
            {
                for (int i = 0; i < eqtlinfo._probNum; i++)
                {
                    bool NA_flag = true;
                    for (int j = 0; j<slctId.size(); j++)
                    {
                        if (abs(eqtlinfo._sexz[i][slctId[j]] + 9) > 1e-6)
                        {
                            NA_flag = false;
                            break;
                        }
                    }
                    if (!NA_flag) include.push_back(i);
                }
            }
            else
            {
                map<int, int > _incld_id_map;
                for (int i = 0; i<slctId.size(); i++)
                {
                    _incld_id_map.insert(pair<int, int>(slctId[i], i));
                    
                }

                for (int i = 0; i < eqtlinfo._probNum; i++)
                {
                    uint32_t pos1=eqtlinfo._cols[(i<<1)+1];
                    uint32_t pos=eqtlinfo._cols[i<<1];
                    uint32_t num=pos1-pos;
                    if (num>0)
                    {
                        for(int j=0;j<num;j++)
                        {
                            uint32_t rid=eqtlinfo._rowid[pos+j];
                            map<int, int>::iterator iter;
                            iter=_incld_id_map.find(rid);
                            if(iter!=_incld_id_map.end())
                            {
                                include.push_back(i);
                                break;
                            }
                        }
                    }
                }
            }
            
            cout<<include.size()<<" probes to be included."<<endl;
            
            //output
            cout<<"\nsaving eQTL data..."<<endl;
            string esdfile = string(outFileName) + ".esi";
            ofstream smr(esdfile.c_str());
            if (!smr) throw ("Error: can not open the fam file " + esdfile + " to save!");
            for (int i = 0; i <slctId.size(); i++) {
                smr << eqtlinfo._esi_chr[slctId[i]] << '\t' << eqtlinfo._esi_rs[slctId[i]] << '\t' << eqtlinfo._esi_gd[slctId[i]] << '\t' << eqtlinfo._esi_bp[slctId[i]] << '\t' << eqtlinfo._esi_allele1[slctId[i]] << '\t' << eqtlinfo._esi_allele2[slctId[i]] << '\n';
            }
            smr.close();
            cout << slctId.size() << " SNPs have been saved in the file [" + esdfile + "]." << endl;
            
            esdfile = string(outFileName) + ".epi";
            smr.open(esdfile.c_str());
            if (!smr) throw ("Error: can not open the fam file " + esdfile + " to save!");
            for (int i = 0; i <include.size(); i++) {
                smr << eqtlinfo._epi_chr[include[i]] << '\t' << eqtlinfo._epi_prbID[include[i]] << '\t' << eqtlinfo._epi_gd[include[i]] << '\t' << eqtlinfo._epi_bp[include[i]] << '\t' << eqtlinfo._epi_gene[include[i]] << '\t' << eqtlinfo._epi_orien[include[i]] << '\n';
            }
            smr.close();
            cout << include.size() << " probes have been saved in the file [" + esdfile + "]." << endl;
            if(make_esd_flag)
            {
                esdfile = string(outFileName) + ".esd";
                smr.open(esdfile.c_str());
                if (!smr) throw ("Error: can not open the fam file " + esdfile + " to save!");
                if(eqtlinfo._valNum==0)
                {
                    for (int i = 0; i <slctId.size(); i++) {
                        for (int j = 0; j<include.size(); j++)
                        {
                            double v_tmp = eqtlinfo._sexz[include[j]][slctId[i]];
                            if (abs(v_tmp + 9)<1e-6) smr << "NA" << '\t' << "NA" << '\t';
                            else smr << eqtlinfo._bxz[include[j]][slctId[i]] << '\t' << eqtlinfo._sexz[include[j]][slctId[i]] << '\t';
                        }
                        smr << '\n';
                    }
                }
                else
                {
                    vector<int> rk;
                    getRank(slctId, rk);
                    map<int, int > _incld_id_map;
                    for (int i = 0; i<slctId.size(); i++)
                    {
                        _incld_id_map.insert(pair<int, int>(slctId[i], i));
                        
                    }
                    eqtlinfo._bxz.resize(include.size());
                    eqtlinfo._sexz.resize(include.size());
                    for( uint32_t i=0;i<include.size();i++)
                    {
                        eqtlinfo._bxz[i].reserve(slctId.size());
                        eqtlinfo._sexz[i].reserve(slctId.size());
                    }
                    for(uint32_t i=0;i<include.size();i++)
                        for(uint32_t j=0;j<slctId.size();j++)
                        {
                            eqtlinfo._bxz[i][j]=-9;
                            eqtlinfo._sexz[i][j]=-9;
                        }
                    for(uint32_t i=0;i<include.size();i++)
                    {
                        int proid=include[i];
                        int pos=eqtlinfo._cols[proid<<1];
                        int pos1=eqtlinfo._cols[(proid<<1)+1];
                        int num=pos1-pos;
                        for(int j=0;j<num;j++)
                        {
                            uint32_t rid=eqtlinfo._rowid[pos+j];
                            map<int, int>::iterator iter;
                            iter=_incld_id_map.find(rid);
                            if(iter!=_incld_id_map.end())
                            {
                                int sid=iter->second;
                                eqtlinfo._bxz[i][rk[sid]]=eqtlinfo._val[pos+j];
                                eqtlinfo._sexz[i][rk[sid]]=eqtlinfo._val[pos+j+num];
                            }
                        }
                    }
                    for (uint32_t i = 0;i <slctId.size(); i++) {
                        for(uint32_t j=0;j<include.size();j++)
                        {
                            double v_tmp=eqtlinfo._sexz[j][i];
                            if(abs(v_tmp+9)<1e-6) smr<<"NA"<<'\t'<<"NA"<<'\t';
                            else smr<<eqtlinfo._bxz[j][i]<<'\t'<<eqtlinfo._sexz[j][i]<<'\t';
                        }
                        smr<<'\n';
                    }
                }
                smr.close();
                cout << "Beta values and SE values for " << include.size() << " Probes and " << slctId.size() << " SNPs have been saved in the file [" + esdfile + "]." << endl;
            }
            
            if(make_besd_flag)
            {
                esdfile = string(outFileName)+".besd";
                FILE * smr;
                smr = fopen (esdfile.c_str(), "wb");
                if(eqtlinfo._valNum==0)
                {
                    unsigned long bsize=(include.size()*slctId.size()<<1)+1;
                    float* buffer=(float*)malloc (sizeof(float)*bsize);
                    memset(buffer,0,sizeof(float)*bsize);
                    float* ptr=buffer;
                    *ptr++=0.0; //dense flag
                    unsigned long pro_num=include.size();
                    unsigned long snp_num=slctId.size();
                    for(int i=0;i<pro_num;i++)
                    {
                        for(int j=0;j<snp_num;j++)
                            ptr[(i<<1)*snp_num+j]=eqtlinfo._bxz[include[i]][slctId[j]];
                        
                        for(int j=0;j<snp_num;j++)
                            ptr[((i<<1)+1)*snp_num+j]=eqtlinfo._sexz[include[i]][slctId[j]];
                    }
                    fwrite (buffer,sizeof(float), bsize, smr);
                    free(buffer);
                }
                else{
                    
                    vector<uint32_t> cols((include.size()<<1)+1);
                    vector<uint32_t> rowid;
                    vector<float> val;
                    vector<int> rk;
                    vector<float> tmp;
                    getRank(slctId, rk);
                    map<int, int > _incld_id_map;
                    for (int i = 0; i<slctId.size(); i++)
                    {
                        _incld_id_map.insert(pair<int, int>(slctId[i], i));
                        
                    }
                    cols[0]=0.0;
                    for(int i=0;i<include.size();i++)
                    {
                        tmp.clear();
                        uint32_t pos=eqtlinfo._cols[include[i]<<1];
                        uint32_t pos1=eqtlinfo._cols[(include[i]<<1)+1];
                        uint32_t num=pos1-pos;
                        uint32_t real_num=0;
                        for(uint32_t j=0;j<num;j++)
                        {
                            uint32_t rid=eqtlinfo._rowid[pos+j];
                            map<int, int>::iterator iter;
                            iter=_incld_id_map.find(rid);
                            if(iter!=_incld_id_map.end())
                            {
                                int it=iter->second;
                                val.push_back(eqtlinfo._val[pos+j]);
                                tmp.push_back(eqtlinfo._val[pos+j+num]);
                                rowid.push_back(rk[it]);
                                real_num++;
                            }
                        }
                        int ris=rowid.size()-real_num;;
                        for(uint32_t j=0;j<real_num;j++)
                        {
                            val.push_back(tmp[j]);
                            rowid.push_back(rowid[ris+j]);
                        }
                        cols[(i<<1)+1]=real_num+cols[i<<1];
                        cols[i+1<<1]=(real_num<<1)+cols[i<<1];
                    }
                    
                    uint64_t colSize=sizeof(uint32_t)*cols.size();
                    uint64_t rowSize=sizeof(uint32_t)*rowid.size();
                    uint64_t valNum=val.size();
                    uint64_t valSize=sizeof(float)*valNum;
                    uint64_t bufsize=sizeof(float)+sizeof(uint64_t)+colSize+rowSize+valSize;
                    
                    char* buffer=(char*)malloc (sizeof(char)*bufsize);
                    memset(buffer,0,sizeof(char)*bufsize);
                    float ftype=SPARSE_FILE_TYPE_2;
                    memcpy(buffer,&ftype,sizeof(float));
                    char* wptr=buffer+sizeof(float);
                    memcpy(wptr,&valNum,sizeof(uint64_t));
                    wptr+=sizeof(uint64_t);
                    memcpy(wptr,&cols[0],colSize);
                    wptr+=colSize;
                    memcpy(wptr,&rowid[0],rowSize);
                    wptr+=rowSize;
                    memcpy(wptr,&val[0],valSize);
                    fwrite (buffer,sizeof(char), bufsize, smr);
                    free(buffer);
                }
                fclose (smr);
               
                cout<<"Beta values and SE values for "<<include.size()<<" Probes and "<<slctId.size()<<" SNPs have been saved in the binary file [" + esdfile + "]." <<endl;
            }
        }
        
    }
    
    void lookup(char* outFileName,char* eqtlFileName, char* snplstName, char* problstName, float plookup, bool bFlag)
    {
        eqtlInfo eqtlinfo;
        cout<<endl<<"Reading eQTL summary data..."<<endl;
        if(eqtlFileName != NULL)
        {
            read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
            if (snplstName != NULL) extract_eqtl_snp(&eqtlinfo, snplstName);
            read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
            if(problstName != NULL) extract_prob(&eqtlinfo, problstName);           
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
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                    if(pxz<plookup)
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
            for(uint32_t i=0;i<eqtlinfo._probNum;i++)
            {
                int proid=eqtlinfo._include[i];
                int pos=eqtlinfo._cols[proid<<1];
                int pos1=eqtlinfo._cols[(proid<<1)+1];
                int num=pos1-pos;
                for(int j=0;j<num;j++)
                {
                    double beta=eqtlinfo._val[pos+j];
                    double se=eqtlinfo._val[pos+j+num];
                    double zsxz=beta/se;
                    double pxz=pchisq(zsxz*zsxz, 1);
                    if(pxz<plookup)
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
        
        string smrfile = string(outFileName)+".lkp";
        ofstream smr(smrfile.c_str());
        if (!smr) throw ("Error: can not open the fam file " + smrfile + " to save!");
        
        smr << "SNP" <<'\t'<< "Chr" <<'\t' << "BP"  << '\t' << "A1" << '\t'<< "A2"<< '\t' << "Probe"<< '\t' << "Probe_Chr"<< '\t'<< "Probe_bp"<< '\t'<<"Gene"<<'\t'<<"b"<<'\t'<< "SE" << '\t'<<"p"<<'\n';
        
        for (int i = 0;i <out_esi_id.size(); i++) {
            smr<<eqtlinfo._esi_rs[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_chr[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_bp[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele1[out_esi_id[i]]<<'\t'<<eqtlinfo._esi_allele2[out_esi_id[i]]<<'\t'<<eqtlinfo._epi_prbID[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_chr[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_bp[out_epi_id[i]]<<'\t'<<eqtlinfo._epi_gene[out_epi_id[i]]<<'\t'<<out_beta[i]<<'\t'<<out_se[i]<<'\t'<<out_pval[i]<< '\n';
        }
       
        smr.close();
        cout<<"Extracted results of "<<out_esi_id.size()<<" SNPs have been saved in the file [" + smrfile + "]."<<endl;
        
    }
    
    void combineCis(char* eqtlsmaslstName, char* outFileName)
    {
        
        eqtlInfo eqtlinfo;
        vector<string> smasNames;
        
        read_smaslist(smasNames, string(eqtlsmaslstName));
        combine_esi(&eqtlinfo, smasNames);
        combine_epi_esd(&eqtlinfo, smasNames, string(outFileName));
        
    }
}
