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
        for (int i = 0; i < eqtlinfo->_snpNum; i++) eqtlinfo->_esi_include.push_back(i);

		
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
            //StrFunc::to_upper(cbuf);
            bdata->_allele1.push_back(cbuf.c_str()[0]);
            Bim >> cbuf;
            //StrFunc::to_upper(cbuf);
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
        }
       cout <<"GWAS summary statistics of "<<gdata->snpNum << " SNPs to be included from [" + string(gwasFileName) + "]." << endl;
        gwasFile.close();
    }
    
    
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
        for(int i=0;i<lineNum;i++)
        {
            string tmpStr;
            esi.getline(buf,MAX_LINE_SIZE);
            istringstream iss(buf);
            iss>>tmpStr;
            eqtlinfo->_esi_chr[i]=atoi(tmpStr.c_str());
            iss>>tmpStr;
            eqtlinfo->_esi_rs[i]=tmpStr;
            iss>>tmpStr;
            eqtlinfo->_esi_gd[i]=atoi(tmpStr.c_str());
            iss>>tmpStr;
            eqtlinfo->_esi_bp[i]=atoi(tmpStr.c_str());
            iss>>tmpStr;
            eqtlinfo->_esi_allele1[i]=tmpStr.c_str()[0];
            iss>>tmpStr;
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
        if((int)*flag){
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
                vector<int> rk;
                getRank(eqtlinfo->_esi_include, rk); // if eqtlinfo->_esi_include is sorted, rk is the indices.
                eqtlinfo->_cols.resize((eqtlinfo->_include.size()<<1)+1);
                eqtlinfo->_cols[0]=*ptr;
                uint32_t* row_ptr;
                row_ptr=ptr+colNum;
                float* val_ptr;
                val_ptr=(float*)(row_ptr+valNum);
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
                        
                        long sid=find(eqtlinfo->_esi_include.begin(),eqtlinfo->_esi_include.end(),rid)-eqtlinfo->_esi_include.begin();
                        if(sid<eqtlinfo->_esi_include.size())
                        {
                            eqtlinfo->_rowid.push_back(rk[sid]);
                            eqtlinfo->_val.push_back(*(val_ptr+pos+j));
                            real_num++;
                        }
                        
                        /*
                         if(find(eqtlinfo->_esi_include.begin(),eqtlinfo->_esi_include.end(),rid)!=eqtlinfo->_esi_include.end()) //can be optimized
                         {
                         eqtlinfo->_rowid.push_back(rid);
                         eqtlinfo->_val.push_back(*(val_ptr+pos+j));
                         real_num++;
                         }
                         */
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
        else
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
                
                char* buff;
                uint64_t buffszie=0x40000000;
                buff = (char*) malloc (sizeof(char)*buffszie);
                if (buff == NULL) {fputs ("Memory error",stderr); exit (1);}
                memset(buff,0,sizeof(char)*buffszie);
                
                uint64_t perbeta=(eqtlinfo->_snpNum<<2);
                uint64_t probonce=sizeof(char)*buffszie/perbeta;
                uint64_t readsize=perbeta*probonce;
                uint64_t probcount=0;
                while(!besd.eof())
                {
                    besd.read(buff,readsize);
                    unsigned long Bread=besd.gcount();
                    char* ptr=buff;
                    while(Bread)
                    {
                        memcpy(&eqtlinfo->_bxz[probcount][0],ptr,perbeta);
                        ptr+=perbeta;
                        memcpy(&eqtlinfo->_sexz[probcount++][0],ptr,perbeta);
                        ptr+=perbeta;
                        Bread-=(perbeta<<1);
                    }
                     cout<<probcount<<" done! "<<endl;
                }
                cout<<"eQTL summary-level statistics of "<<eqtlinfo->_probNum<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs to be included from [" + besdfile + "]." <<endl;
                free(buff);
            }
            free(buffer);
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
        
        float* tmpX=new float[rsize];
        float* tmpX2=new float[rsize];
        ::memset(tmpX,0,rsize*sizeof(float));
        ::memset(tmpX2,0,rsize*sizeof(float));
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
        
        delete[] tmpX;
        delete[] tmpX2;
    }
    
    void get_square_idxes(vector<int> &sn_ids,vector<double> &zsxz,double threshold)
    {
        for(int i=0;i<zsxz.size();i++)
        {
           if(zsxz[i]*zsxz[i]-threshold>1e-6) sn_ids.push_back(i);
        }
    }
    void get_square_ldpruning_idxes(vector<int> &sn_ids,vector<double> &zsxz,double threshold,double* ref_ld, int maxid,double ld_top)
    {
        
        for(int i=0;i<zsxz.size();i++)
        {
            if(i<maxid)
            {
                if((zsxz[i]*zsxz[i]-threshold)>1e-6 && (ref_ld[maxid*(maxid-1)/2+i]-ld_top)<1e-6) sn_ids.push_back(i);
            }
            else if(i>maxid)
            {
                if((zsxz[i]*zsxz[i]-threshold)>1e-6 && (ref_ld[i*(i-1)/2+maxid]-ld_top)<1e-6) sn_ids.push_back(i);
            }
            else{
                 if((zsxz[i]*zsxz[i]-threshold)>1e-6) sn_ids.push_back(i);
            }
            
        }
        
        
        
    }
    void est_cov_bxy(double* covbxy, vector<double> zsxz1,vector<double> bxy1,vector<double> seyz1,vector<double> bxz1,double* ref_ld1)
    {
        int nsnp =zsxz1.size();
        if(nsnp>1)
        {
            double tmpval1, tmpval2;
            
            for(int i=0;i<nsnp;i++)
            {
                for(int j=0;j<i;j++)
                {
                    tmpval1=bxy1[i]*bxy1[j];
                    tmpval2=zsxz1[i]*zsxz1[j];
                    covbxy[i*nsnp+j]=covbxy[j*nsnp+i]=ref_ld1[i*(i-1)/2+j]*seyz1[i]*seyz1[j]/(bxz1[i]*bxz1[j]) + ref_ld1[i*(i-1)/2+j]*tmpval1/tmpval2 - tmpval1/(tmpval2*tmpval2);
                }
            }
            for(int i=0;i<nsnp;i++)
            {
                tmpval1=bxy1[i]*bxy1[i];
                tmpval2=zsxz1[i]*zsxz1[i];
                covbxy[i*nsnp+i]=seyz1[i]*seyz1[i]/(bxz1[i]*bxz1[i]) + tmpval1/tmpval2 - tmpval1/(tmpval2*tmpval2);
            }            
           
        }
    }

    double bxy_hetero3(vector<double> byz1, vector<double> bxz1,vector<double> seyz1,vector<double> sexz1,vector<double> zsxz1,double* ref_ld1, long* snp_num)
    {
        vector<double> bxy1;
        vector<double> sexy1;
        vector<double> dev;
        long nsnp=*snp_num;
        int maxid;        
        double pdev=-1.0;
        double* covbxy=new double[nsnp*nsnp];
        ::memset(covbxy,0,nsnp*nsnp*sizeof(double));
        double* vdev=new double[(nsnp-1)*(nsnp-1)];
        ::memset(vdev,0,(nsnp-1)*(nsnp-1)*sizeof(double));
        dev.resize(nsnp-1);
        if(nsnp>1)
        {
            for(int j=0;j<nsnp;j++)
            {   bxy1.push_back(byz1[j]/bxz1[j]);
                sexy1.push_back(sqrt((seyz1[j]*seyz1[j]*bxz1[j]*bxz1[j]+ sexz1[j]*sexz1[j]*byz1[j]*byz1[j]) / (bxz1[j]*bxz1[j]*bxz1[j]*bxz1[j])));
            }
            maxid=max_abs_id(zsxz1);
            
            for(int j=0;j<maxid;j++) dev[j]=bxy1[maxid]-bxy1[j];
            for(int j=maxid+1;j<nsnp;j++) dev[j-1]=bxy1[maxid]-bxy1[j];
            
            est_cov_bxy(covbxy, zsxz1,bxy1,seyz1,bxz1,ref_ld1);
            
            double tmp1=covbxy[maxid*nsnp+maxid];
            double* tmp3=new double[nsnp-1];
            ::memset(tmp3,0,(nsnp-1)*sizeof(double));
            for(int i=0; i<maxid; i++) tmp3[i]=covbxy[maxid*nsnp+i];
            for(int i=maxid+1; i<nsnp; i++) tmp3[i-1]=covbxy[maxid*nsnp+i];
            // vdev as tmp2
            for(int i=0; i<maxid; i++)
            {
                for(int j=0;j<=i;j++) vdev[i*(nsnp-1)+j]=vdev[j*(nsnp-1)+i]= covbxy[i*nsnp+j];
            }
            
            for(int i=maxid+1; i<nsnp; i++)
            {
                for(int j=0;j<maxid;j++) vdev[(i-1)*(nsnp-1)+j]=vdev[j*(nsnp-1)+i-1]= covbxy[i*nsnp+j];
                for(int j=maxid+1;j<=i;j++) vdev[(i-1)*(nsnp-1)+j-1]=vdev[(j-1)*(nsnp-1)+i-1]= covbxy[i*nsnp+j];
            }
            // get vdev
            for(int i=0;i<nsnp-1;i++)
            {
                for(int j=0;j<=i;j++) vdev[i*(nsnp-1)+j]=vdev[j*(nsnp-1)+i]=tmp1+vdev[i*(nsnp-1)+j]-tmp3[i]-tmp3[j];
            }
            
            for(int i=0;i<nsnp-1;i++)  vdev[i*(nsnp-1)+i]+=1e-8; // in R code
            
            //tmp3 as vardev
            
            for(int i=0; i<maxid; i++) tmp3[i]=tmp1+covbxy[i*nsnp+i]-2*tmp3[i]+ 1e-8;
            for(int i=maxid+1; i<nsnp; i++) tmp3[i-1]=tmp1+covbxy[i*nsnp+i]-2*tmp3[i-1]+ 1e-8;
            
            //dev as chisq_dev
            
            for(int i=0;i<nsnp-1;i++) dev[i]=dev[i]*dev[i]/tmp3[i];
            
            double sumChisq_dev=0.0;
            for(int i=0;i<nsnp-1;i++)sumChisq_dev+=dev[i];
            
            //using covbxy to store corr_dev
            memset(covbxy,0,nsnp*nsnp*sizeof(double));
            for( int i=0;i<nsnp-1;i++)
            {
                for( int j=0;j<=i;j++) {
                    covbxy[i*(nsnp-1)+j] = covbxy[j*(nsnp-1)+i] = vdev[i*(nsnp-1)+j]/sqrt(vdev[i*(nsnp-1)+i]*vdev[j*(nsnp-1)+j]);
                }
            }
            
            // using Eigen Library
            
            Map<MatrixXd> A(covbxy,nsnp-1,nsnp-1);
            SelfAdjointEigenSolver<MatrixXd> es(A);
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
            
            delete[] tmp3;          
         
        }else *snp_num=-9;
        
        delete[] vdev;
        delete[] covbxy;
        
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
        //stable_sort(eqtlinfo->_esi_include.begin(), eqtlinfo->_esi_include.end());
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

    void smr(char* outFileName, char* bFileName,char* gwasFileName, char* eqtlFileName, double maf,char* indilstName, char* snplstName,char* problstName,bool bFlag,double p_hetero,double ld_top,int m_hetero , char* indilst2remove, char* snplst2exclde, char* problst2exclde,double p_smr)
    {
        bInfo bdata;
        gwasData gdata;
        eqtlInfo eqtlinfo;
        double threshold= chi_val(1,p_hetero);
        
        if(bFileName == NULL ) throw("Error: please input Plink file for SMR analysis by the flag --bfile.");
        if(gwasFileName==NULL) throw("Error: please input GWAS summary data for SMR analysis by the flag --gwas-summary.");
        if(eqtlFileName==NULL) throw("Error: please input eQTL summary data for SMR analysis by the flag --eqtl-summary.");
        
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
        read_gwas_data( &gdata, gwasFileName);
        cout<<endl<<"Reading eQTL summary data..."<<endl;
        read_esifile(&eqtlinfo, string(eqtlFileName)+".esi");
        if (snplstName != NULL) extract_eqtl_snp(&eqtlinfo, snplstName);
        if(snplst2exclde != NULL) exclude_eqtl_snp(&eqtlinfo, snplst2exclde); //no switch the place ahead of extract_eqtl_snp()
        read_epifile(&eqtlinfo, string(eqtlFileName)+".epi");
        if(problstName != NULL) extract_prob(&eqtlinfo, problstName);
        if(problst2exclde != NULL) exclude_prob(&eqtlinfo, problst2exclde); //no switch the place ahead of extract_prob()
        if(bFlag) read_besdfile(&eqtlinfo, string(eqtlFileName)+".besd");
        else      read_esdfile(&eqtlinfo, string(eqtlFileName)+".esd");
        
        // get the common SNPs
        vector<string> slctSNPs;
        vector<string> b_snp_name;
        vector<string> cmmnSNPs;
        vector<int> bdId;
        vector<int> gdId;
        vector<int> edId;
        cmmnSNPs.clear();
        edId.clear();
        if(maf>0)
        {
            b_snp_name.resize(bdata._include.size());
            for(int i=0;i<bdata._include.size();i++)
                b_snp_name[i]=bdata._snp_name[bdata._include[i]];
            StrFunc::match_only(b_snp_name, gdata.snpName, edId);
        }
        else
            StrFunc::match_only(bdata._snp_name, gdata.snpName, edId);
        if(edId.empty()) throw("Error: no common SNPs between PLink data and GWAS data.");
        for(int i=0;i<edId.size();i++) cmmnSNPs.push_back(gdata.snpName[edId[i]]);
        edId.clear();
        StrFunc::match_only(cmmnSNPs, eqtlinfo._esi_rs, edId);
        if(edId.empty()) throw("Error: no common SNPs found.");
        for(int i=0;i<edId.size();i++) slctSNPs.push_back(eqtlinfo._esi_rs[edId[i]]);
        
        
        //alleles check
        StrFunc::match(slctSNPs, bdata._snp_name, bdId);
        StrFunc::match(slctSNPs, gdata.snpName, gdId);
        cmmnSNPs.clear();
        for (int i = 0; i<slctSNPs.size(); i++)
        {
            char a1, a2, ga1, ga2, ea1, ea2;
            a1 = bdata._allele1[bdId[i]];
            a2 = bdata._allele2[bdId[i]];
            ga1 = gdata.allele_1[gdId[i]];
            ga2 = gdata.allele_2[gdId[i]];
            ea1 = eqtlinfo._esi_allele1[edId[i]];
            ea2 = eqtlinfo._esi_allele2[edId[i]];
            if (a1 != '0' && a2 != '0' && ga1 != '0' && ga2 != '0' && ea1 != '0' && ea2 != '0')
                if ((a1 == ga1 && a2 == ga2) || (a1 == ga2 && a2 == ga1) || (a1 != ga1 && a1 != ga2 && a2 != ga2 && a2 != ga1))
                    if ((a1 == ea1 && a2 == ea2) || (a1 == ea2 && a2 == ea1) || (a1 != ea1 && a1 != ea2 && a2 != ea2 && a2 != ea1))
                        cmmnSNPs.push_back(slctSNPs[i]);
        }
        
        // get indices
        bdId.clear();
        gdId.clear();
        edId.clear();
        StrFunc::match(cmmnSNPs, bdata._snp_name, bdId);
        StrFunc::match(cmmnSNPs, gdata.snpName, gdId);
        StrFunc::match(cmmnSNPs, eqtlinfo._esi_rs, edId);
        
        
        int outCount = -1;
        unsigned int probNum = eqtlinfo._probNum;
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
        for(int i=0;i<probNum;i++)
        {
            
            progr1=1.0*i/probNum;
            if(progr1-progr0-0.05>1e-6 || i+1==probNum)
            {
                if(i+1==probNum) progr1=1.0;
                progress_print(progr1);
                progr0=progr1;
            }
            
            //extract info from eqtl
            vector<float> bxz;
            vector<float> sexz;
            vector<uint32_t> curId;
            vector<string> eName;
            if(eqtlinfo._rowid.empty())
            {
                for (int j = 0; j<edId.size(); j++)
                    if (abs(eqtlinfo._bxz[i][edId[j]] + 9) > 1e-6)
                    {
                        bxz.push_back(eqtlinfo._bxz[i][edId[j]]);
                        sexz.push_back(eqtlinfo._sexz[i][edId[j]]);
                        curId.push_back(j);
                        eName.push_back(eqtlinfo._esi_rs[edId[j]]);
                    }
            }
            else{
                int beta_start=eqtlinfo._cols[i<<1];
                int se_start=eqtlinfo._cols[1+(i<<1)];
                int numsnps=se_start-beta_start;
                for(int j=0;j<numsnps;j++)
                {
                    int pos=find(edId.begin(),edId.end(),eqtlinfo._rowid[beta_start+j])-edId.begin();
                    if(pos<edId.size())
                    {
                        bxz.push_back(eqtlinfo._val[beta_start+j]);
                        sexz.push_back(eqtlinfo._val[se_start+j]);
                        curId.push_back(pos);
                        eName.push_back(eqtlinfo._esi_rs[edId[pos]]);
                    }
                }
            }
            
            if (curId.size() == 0) continue;
            
            //extract info from reference
            vector<char> allele1(curId.size());
            vector<char> allele2(curId.size());
            vector<string> reName(curId.size());
            vector<string> tstName(curId.size());
            vector<unsigned int> bpsnp(curId.size());
            float* ref_snpData; //row major: snps belonging to the the same individual store adjacently
            ref_snpData=new float[bdata._indi_num*curId.size()];
            ::memset(ref_snpData,0,bdata._indi_num*curId.size()*sizeof(float));
            for (int j = 0; j<curId.size(); j++)
            {
                allele1[j] = bdata._allele1[bdId[curId[j]]]; // string to char
                allele2[j] = bdata._allele2[bdId[curId[j]]];
                bpsnp[j] = bdata._bp[bdId[curId[j]]];
                reName[j]=bdata._snp_name[bdId[curId[j]]];
                tstName[j]=cmmnSNPs[curId[j]];
                makeptrx(&bdata,bdId[curId[j]], j, ref_snpData);
            }
            
            
            // extract info from gwas data. NOTE: the gwas data is not same as the one that R code uses.
            vector<char> allele_1(curId.size());
            vector<char> allele_2(curId.size());
            vector<double> byz(curId.size());
            vector<double> seyz(curId.size());
            vector<string> gName(curId.size());
            for (int j = 0; j<curId.size(); j++)
            {
                allele_1[j] = gdata.allele_1[gdId[curId[j]]];
                allele_2[j] = gdata.allele_2[gdId[curId[j]]];
                byz[j] = gdata.byz[gdId[curId[j]]];
                seyz[j] = gdata.seyz[gdId[curId[j]]];
                gName[j]=gdata.snpName[gdId[curId[j]]];
            }
            //get LD info from reference
            double* ref_ld;
            ref_ld=new double[curId.size()*(curId.size()-1)/2];
            ::memset(ref_ld,0,curId.size()*(curId.size()-1)*sizeof(double)/2);
            ld_calcualte(ref_ld,ref_snpData,curId.size(),bdata._indi_num);
            
            
            // top_snp_test();
            vector<double> zsxz(curId.size());
            for(int j=0;j<curId.size();j++) zsxz[j]=bxz[j]/sexz[j];
            int maxid=max_abs_id(zsxz);
            
            double bxy_val = byz[maxid]/bxz[maxid];
            double sexy_val = sqrt((seyz[maxid]*seyz[maxid]*bxz[maxid]*bxz[maxid]+ sexz[maxid]*sexz[maxid]*byz[maxid]*byz[maxid]) / (bxz[maxid]*bxz[maxid]*bxz[maxid]*bxz[maxid]) );
            double chisqxy = bxy_val*bxy_val/(sexy_val*sexy_val);
            
            
            double pxy_val=pchisq(chisqxy, 1);  //   double pxy=chi_prob(1,chisqxy); //works
            double chisqyz=byz[maxid]/seyz[maxid];
            double pyz_val = pchisq(chisqyz*chisqyz, 1);
            double pxz_val = pchisq(zsxz[maxid]*zsxz[maxid], 1);
            
            if(pxz_val>p_smr) continue;
            else outCount++;
            
            
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
            rsid[outCount]=cmmnSNPs[curId[maxid]];
            rsbp[outCount]=itos(bpsnp[maxid]);
            rsa1[outCount]=allele1[maxid];
            rsa2[outCount]=allele2[maxid];
            
            
            vector<int> sn_ids; //increase order
            if(abs(ld_top-1)<1e-6) get_square_idxes(sn_ids,zsxz,threshold);
            else get_square_ldpruning_idxes(sn_ids,zsxz,threshold,ref_ld, maxid,ld_top);
            
            vector<double> byz1;
            vector<double> seyz1;
            vector<double> bxz1;
            vector<double> sexz1;
            vector<double> zsxz1;
            double* ref_ld1;
            ref_ld1=new double[sn_ids.size()*(sn_ids.size()-1)/2];
            ::memset(ref_ld1,0,sn_ids.size()*(sn_ids.size()-1)*sizeof(double)/2);
            
            for(int j=0;j<sn_ids.size();j++)
            {
                byz1.push_back(byz[sn_ids[j]]);
                seyz1.push_back(seyz[sn_ids[j]]);
                bxz1.push_back(bxz[sn_ids[j]]);
                sexz1.push_back(sexz[sn_ids[j]]);
                zsxz1.push_back(zsxz[sn_ids[j]]);
                for(int k=0;k<j;k++)ref_ld1[j*(j-1)/2+k]=ref_ld[sn_ids[j]*(sn_ids[j]-1)/2+sn_ids[k]];
            }
            
            if(sn_ids.size() < m_hetero)
            {
                prb1[outCount]= string("NA");
                nsnp_test1[outCount]= string("NA");
                top_match1[outCount]= string("NA");
                ldrsq[outCount]= string("NA");
                continue;
            }
            
            long nsnp = sn_ids.size();
            double pdev=bxy_hetero3(byz1,  bxz1, seyz1, sexz1, zsxz1,ref_ld1, &nsnp);
            
            prb1[outCount] = dtos(pdev);
            if(nsnp>0) nsnp_test1[outCount] = itos(nsnp);
            else nsnp_test1[outCount] = string("NA");
            
            
            // top GWAS SNP ?= top eQTL
            int indx1=max_abs_id(zsxz1);
            for(int j=0;j<byz1.size();j++) byz1[j]=byz1[j]/seyz1[j];
            int indx2=max_abs_id(byz1);
            
            
            if(indx1 == indx2) top_match1[outCount] =itos(1);
            else top_match1[outCount] = itos(0);
            double ldrsqVal;
            if(indx1 > indx2) ldrsqVal = ref_ld1[indx1*(indx1-1)/2+indx2]*ref_ld1[indx1*(indx1-1)/2+indx2];
            else if(indx1 < indx2) ldrsqVal = ref_ld1[indx2*(indx2-1)/2+indx1]*ref_ld1[indx2*(indx2-1)/2+indx1];
            else ldrsqVal=1;
            ldrsq[outCount]=dtosf(ldrsqVal);
            
            delete[] ref_ld;
            delete[] ref_ld1;
            delete[] ref_snpData;
            
        }
        
        //genesn in R is just probidx here. can refer to probinfo to get probeid, probenm, chr, gene, bp
        
        string smrfile = string(outFileName)+".smr";
        ofstream smr(smrfile.c_str());
        if (!smr) throw ("Error: can not open the fam file " + smrfile + " to save!");
        
        smr << "ProbeID" <<'\t'<< "Chr" <<'\t' << "Gene"  << '\t' << "Prob_bp" << '\t'<< "SNP"<< '\t' << "SNP_bp"<< '\t' << "A1"<< '\t'<< "A2"<< '\t'<<"b_GWAS"<<'\t'<<"se_GWAS"<<'\t'<< "p_GWAS" << '\t'<<"b_eQTL"<<'\t'<<"se_eQTL"<<'\t'<< "p_eQTL" << '\t'<< "b_SMR" << '\t'<< "se_SMR"<< '\t' << "p_SMR" << "\t"<< "p_HET"<< "\t" << "nsnp" << '\n';
        
        for (int i = 0;i <=outCount; i++) {
            smr<<eqtlinfo._epi_prbID[out_probid[i]]<<'\t'<<eqtlinfo._epi_chr[out_probid[i]]<<'\t'<<eqtlinfo._epi_gene[out_probid[i]]<<'\t'<<eqtlinfo._epi_bp[out_probid[i]]<<'\t'<<rsid[i]<<'\t'<<rsbp[i]<<'\t'<<rsa1[i]<<'\t'<<rsa2[i]<<'\t'<<bgwas[i]<<'\t'<<segwas[i]<<'\t'<<pgwas[i]<<'\t'<<beqtl[i]<<'\t'<<seeqtl[i]<<'\t'<<peqtl[i]<<'\t'<<bxy[i]<<'\t'<<sexy[i]<<'\t'<<pxy[i]<<'\t'<<prb1[i]<<'\t'<<nsnp_test1[i]<<'\n';
        }
        
        /*
         smr << "probeid" << "\t" << "chr" << "\t" << "gene" << "\t" << "prob_bp" << "\t" << "snp"<< "\t" << "snp_bp"<< "\t"<< "a1"<< "\t"<< "a2"<< "\t"<< "pyz" << "\t"<< "pxz" << "\t"<< "bxy" << "\t"<< "sexy"<< "\t" << "pxy" << "\t"<< "prb1"<< "\t" << "nsnp_test1" << "\t"<< "top_match1"<< "\t" << "ldrsq" << endl;
         
         for (int i = 0;i <=outCount; i++) {
         smr<<eqtlinfo._epi_prbID[out_probid[i]]<<"\t"<<eqtlinfo._epi_chr[out_probid[i]]<<"\t"<<eqtlinfo._epi_gene[out_probid[i]]<<"\t"<<eqtlinfo._epi_bp[out_probid[i]]<<"\t"<<rsid[i]<<"\t"<<rsbp[i]<<"\t"<<rsa1[i]<<"\t"<<rsa2[i]<<"\t"<<pgwas[i]<<"\t"<<peqtl[i]<<"\t"<<bxy[i]<<"\t"<<sexy[i]<<"\t"<<pxy[i]<<"\t"<<prb1[i]<<"\t"<<nsnp_test1[i]<<"\t"<<top_match1[i]<<"\t"<<ldrsq[i]<<'\n';
         }
         */
        smr.close();
        cout<<"SMR and heterogeneity analysis finished.\nSMR and heterogeneity analysis results of "<<outCount+1<<" probes have been saved in the file [" + smrfile + "]."<<endl;
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
                    if ((a1 == ga1 && a2 == ga2) || (a1 == ga2 && a2 == ga1) || (a1 != ga1 && a1 != ga2 && a2 != ga2 && a2 != ga1))
                    {
                        slctSNPs.push_back(cmmnSNPs[i]);
                        A1.push_back(a1);
                        A2.push_back(a2);
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
                if(eqtlinfo._valNum>0) throw("Error: please input dense format eQTL summary data file.");
                //if got flag (--extract-snp), filter_probe_null() should be invoked, then _include.
                vector<uint32_t> cols((eqtlinfo._probNum<<1)+1);
                vector<uint32_t> rowids;
                vector<float> val;
                
                cis=cis*1e6;
                cols[0]=0;
                for(uint32_t i=0;i<eqtlinfo._probNum;i++)
                {
                    vector<int> esi_include;
                    int probchr=eqtlinfo._epi_chr[i];
                    int probbp= eqtlinfo._epi_bp[i];
                    uint32_t uperBounder=probbp+cis;
                    uint32_t lowerBounder=((probbp-cis>0)?(probbp-cis):0);
                    for(int j=0;j<eqtlinfo._snpNum;j++)
                    {
                        double zsxz=eqtlinfo._bxz[i][j]/eqtlinfo._sexz[i][j];
                        double pxz=pchisq(zsxz*zsxz, 1);
                        if(eqtlinfo._esi_chr[j] == probchr && eqtlinfo._esi_bp[j]<=uperBounder && eqtlinfo._esi_bp[j]>=lowerBounder && abs(eqtlinfo._sexz[i][j]+9)>1e-6)
                        {
                            esi_include.push_back(j);
                        }
                        else if(pxz<transThres)
                        {
                            esi_include.push_back(j);
                            long transbp=eqtlinfo._esi_bp[j];
                           
                            //if esi is not sorted.                            
                           // for(int k=0;k<eqtlinfo._snpNum;k++)
                           //     if(k!=j && eqtlinfo._esi_chr[j] == eqtlinfo._esi_chr[k]  && abs(transbp-eqtlinfo._esi_bp[k])<=cis)   esi_include.push_back(k);
                            
                            // if esi is sorted
                           
                            int startptr=j-1;
                            while(startptr>=0 && eqtlinfo._esi_chr[j] == eqtlinfo._esi_chr[startptr] && abs(transbp-eqtlinfo._esi_bp[startptr])<=cis)
                            {
                                esi_include.push_back(startptr);
                                startptr--;
                                
                            }
                            startptr=j+1;
                            while(startptr<=eqtlinfo._snpNum && eqtlinfo._esi_chr[j] == eqtlinfo._esi_chr[startptr] && abs(transbp-eqtlinfo._esi_bp[startptr])<=cis)
                            {
                                esi_include.push_back(startptr);
                                startptr++;
                                //j=startptr-1; wrong!!! is there are trans upper, then the trans region should extend.
                            }
                           
                        }
                        else if(pxz<restThres) esi_include.push_back(j);
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
                float ftype=SPARSE_FILE_TYPE_1;
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
                
                free(buffer);
                cout<<"Beta values and SE values for "<<eqtlinfo._include.size()<<" Probes and "<<eqtlinfo._snpNum<<" SNPs have been saved in the binary file [" + esdfile + "]." <<endl;
                
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
                        uint64_t colSize=sizeof(uint32_t)*(eqtlinfo._cols.size());
                        uint64_t rowSize=sizeof(uint32_t)*eqtlinfo._valNum;
                        uint64_t valSize=sizeof(float)*eqtlinfo._valNum;
                        uint64_t valNum=eqtlinfo._valNum;
                        uint64_t bufsize=sizeof(float)+sizeof(uint64_t)+colSize+rowSize+valSize;
                        
                        char* buffer=(char*)malloc (sizeof(char)*bufsize);
                        memset(buffer,0,sizeof(char)*bufsize);
                        float ftype=SPARSE_FILE_TYPE_1;
                        memcpy(buffer,&ftype,sizeof(float));
                        char* wptr=buffer+sizeof(float);
                        memcpy(wptr,&valNum,sizeof(uint64_t));
                        wptr+=sizeof(uint64_t);
                        memcpy(wptr,&eqtlinfo._cols[0],colSize);
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
                            if(find(slctId.begin(),slctId.end(),rid)!=slctId.end())
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
                            long sid=find(slctId.begin(),slctId.end(),rid)-slctId.begin();
                            if(sid<slctId.size())
                            {
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
                            int it=find(slctId.begin(),slctId.end(),rid)-slctId.begin();
                            if(it<slctId.size())
                            {
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
                    uint64_t bufsize=1+sizeof(uint64_t)+colSize+rowSize+valSize;
                    
                    char* buffer=(char*)malloc (sizeof(char)*bufsize);
                    memset(buffer,0,sizeof(char)*bufsize);
                    float ftype=SPARSE_FILE_TYPE_1;
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
    
}
