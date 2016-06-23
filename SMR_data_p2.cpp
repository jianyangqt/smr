//
//  SMR_data_p2.cpp
//  SMR_CPP
//
//  Created by Futao Zhang on 10/06/2016.
//  Copyright (c) 2016 Futao Zhang. All rights reserved.
//

#include "SMR_data_p2.h"

namespace SMRDATA
{
    // below for txt 2 besd
    int read_probeinfolst(vector<probeinfolst> &prbiflst,char* syllabusName)
    {
        ifstream flptr(syllabusName);
        if (!flptr) {
            printf("Error: can not open the file %s to read.\n",syllabusName);
            exit(EXIT_FAILURE);
        }
        printf( "Reading eQTL probe information from \"%s\".\n", syllabusName);
        map<string, int> probe_map;
        map<string, int> probe_cp_map;
        map<string, int> probe_cpb_map;
        long mapsize=0;
        char buf[MAX_LINE_SIZE];
        flptr.getline(buf, MAX_LINE_SIZE); //header
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, " \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: %s should have headers that start with \"Chr\".\n", syllabusName);
            exit(EXIT_FAILURE);
        }
        int lineNum(0);
        while(!flptr.eof())
        {
            flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0') lineNum++;
        }
        flptr.clear(ios::goodbit);
        flptr.seekg (0, ios::beg);
        flptr.getline(buf, MAX_LINE_SIZE); //header
        int fcount=0;
        while(!flptr.eof())
        {
            flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                col_num = split_string(buf, vs_buf, " \t\n");
                if(col_num!=7 && col_num!=8) {
                    printf("ERROR: column number is not right in row %d with probe ID %s!\n", fcount+2, vs_buf[1].c_str());
                    exit(EXIT_FAILURE);
                }
                to_upper(vs_buf[0]);
                to_upper(vs_buf[3]);
                if(vs_buf[0]=="NA" || vs_buf[3]=="NA")
                {
                    printf("ERROR: \"NA\" found in Choromosome or probe position of probe %s!\n", vs_buf[1].c_str());
                    exit(EXIT_FAILURE);
                }
                string cpstr=vs_buf[0]+":"+vs_buf[1];
                string cpbstr=vs_buf[0]+":"+vs_buf[1]+":"+vs_buf[3];
                probe_cp_map.insert(pair<string,int>(cpstr,fcount));
                probe_cpb_map.insert(pair<string,int>(cpbstr,fcount));
                probe_map.insert(pair<string,int>(vs_buf[6],fcount));
                if(probe_cp_map.size()!=probe_cpb_map.size())
                {
                    printf("ERROR: Different BPs of probe %s found.\n",vs_buf[1].c_str() );
                    exit(EXIT_FAILURE);
                }
                if(mapsize==probe_map.size())
                {
                    printf("ERROR: Duplicate esd file name %s found.\n",vs_buf[6].c_str() );
                    exit(EXIT_FAILURE);
                } else {
                    mapsize=probe_map.size();
                }
                probeinfolst tmpinfo;
                int tmpchr;
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[0].c_str());
                tmpinfo.probechr=tmpchr;
                tmpinfo.probeId=vs_buf[1].c_str();
                tmpinfo.gd=atoi(vs_buf[2].c_str());
                tmpinfo.bp=atoi(vs_buf[3].c_str());
                tmpinfo.genename=vs_buf[4].c_str();
                tmpinfo.orien=vs_buf[5].c_str()[0];
                tmpinfo.esdpath=vs_buf[6].c_str();
                if(col_num==8) tmpinfo.bfilepath=vs_buf[7].c_str();
                prbiflst.push_back(tmpinfo);
                fcount++;
            }
        }
        flptr.close();
        printf("%d esd files of total %ld probes infomation to be included from \"%s\".\n",lineNum, probe_cp_map.size(), syllabusName);
        return lineNum;
    }
    void read_smr_sa(vector<string> &rs,vector<int> &chr, vector<int> &bp, vector<string> &a1, vector<string> &a2, vector<float> &beta,vector<float> &se,string esdpath)
    {
        gzFile gzfile=NULL;
        ifstream flptr;
        bool gzflag=has_suffix(esdpath, "gz");
        if(gzflag)
        {
            gzfile = gzopen(esdpath.c_str(), "rb");
            if (!(gzfile)) {
                fprintf (stderr, "%s: Couldn't open file %s\n",
                         esdpath.c_str(), strerror (errno));
                exit (EXIT_FAILURE);
            }
        }
        else
        {
            flptr.open(esdpath.c_str());
            if (!flptr) throw ("Error: can not open the file [" + esdpath + "] to read.");
        }
        cout << "Reading eQTL information from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
        else flptr.getline(buf,MAX_LINE_SIZE);// the header
        /* check headers */
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, " \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"Chr\"");
            exit(EXIT_FAILURE);
        }
        map<string, int> currs_map;
        long currssize=0;
        while(!flptr.eof() && !gzeof(gzfile))
        {
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
            else flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                col_num = split_string(buf, vs_buf, " \t\n");
                if(col_num!=9) {
                    printf("ERROR: column number is not right in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                string tmpstr=vs_buf[1]+":"+vs_buf[3]+":"+vs_buf[4];
                currs_map.insert(pair<string,int>(tmpstr,lineNum));
                if(currssize<currs_map.size())
                {
                    if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                        printf("ERROR: \'NA\' of chromosome found in row %d.\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }
                    int tmpchr;
                    if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                    else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                    else tmpchr=atoi(vs_buf[0].c_str());
                    chr.push_back(tmpchr);
                    rs.push_back(vs_buf[1]);
                    if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                        printf("ERROR: \'NA\' of SNP position found in row %d.\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }

                    bp.push_back(atoi(vs_buf[2].c_str()));
                    to_upper(vs_buf[3]);
                    a1.push_back(vs_buf[3]);
                    to_upper(vs_buf[4]);
                    a2.push_back(vs_buf[4]);
                    to_upper(vs_buf[6]);
                    to_upper(vs_buf[7]);
                    if(vs_buf[6].compare("NA") && vs_buf[7].compare("NA")){
                        beta.push_back(atof(vs_buf[6].c_str()));
                        se.push_back(atof(vs_buf[7].c_str()));
                    } else {
                        beta.push_back(-9);
                        se.push_back(-9);
                    }
                    lineNum++;
                    currssize=currs_map.size();
                } else {
                    printf("WARNING: duplicate SNP %s with the same alleles found in file [%s]. This row was skiped:\n", vs_buf[1].c_str(),esdpath.c_str());
                    printf("%s.\n",buf);
                }
            }
        }
        
        if(gzflag) gzclose(gzfile);
        else flptr.close();
        cout << lineNum << " SNPs infomation to be included from [" + string(esdpath) + "]." << endl;
    }
    void read_smr_sa(vector<snpinfolst> &snpinfo,string esdpath)
    {
        gzFile gzfile=NULL;
        ifstream flptr;
        bool gzflag=has_suffix(esdpath, "gz");
        if(gzflag)
        {
            gzfile = gzopen(esdpath.c_str(), "rb");
            if (!(gzfile)) {
                fprintf (stderr, "%s: Couldn't open file %s\n",
                         esdpath.c_str(), strerror (errno));
                exit (EXIT_FAILURE);
            }
        }
        else
        {
            flptr.open(esdpath.c_str());
            if (!flptr) throw ("Error: can not open the file [" + esdpath + "] to read.");
        }
        cout << "Reading eQTL information from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
        else flptr.getline(buf,MAX_LINE_SIZE);// the header
        /* check headers */
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, " \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"Chr\"");
            exit(EXIT_FAILURE);
        }
        map<string, int> currs_map;
        long currssize=0;
        while(!flptr.eof() && !gzeof(gzfile))
        {
            snpinfolst tmpesd;
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
            else flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                col_num = split_string(buf, vs_buf, " \t\n");
                if(col_num!=9) {
                    printf("ERROR: column number is not right in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                string tmpstr=vs_buf[1]+":"+vs_buf[3]+":"+vs_buf[4];
                currs_map.insert(pair<string,int>(tmpstr,lineNum));
                if(currssize<currs_map.size())
                {
                    if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                        printf("ERROR: \'NA\' of chromosome found in row %d.\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }
                    if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpesd.snpchr=23;
                    else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpesd.snpchr=24;
                    else tmpesd.snpchr=atoi(vs_buf[0].c_str());
                    tmpesd.snprs=vs_buf[1].c_str();
                    if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                        printf("ERROR: \'NA\' of SNP position found in row %d.\n", lineNum+2);
                        exit(EXIT_FAILURE);
                    }
                    tmpesd.bp=atoi(vs_buf[2].c_str());
                     to_upper(vs_buf[3]);
                    tmpesd.a1=vs_buf[3].c_str();
                     to_upper(vs_buf[4]);
                    tmpesd.a2=vs_buf[4].c_str();
                     to_upper(vs_buf[6]);
                     to_upper(vs_buf[7]);
                    if(vs_buf[6].compare("NA") && vs_buf[7].compare("NA")){
                        tmpesd.beta=atof(vs_buf[6].c_str());
                        tmpesd.se=atof(vs_buf[7].c_str());
                    } else {
                        tmpesd.beta=-9;
                        tmpesd.se=-9;
                    }
                    snpinfo.push_back(tmpesd);
                    lineNum++;
                    currssize=currs_map.size();
                } else {
                    printf("WARNING: duplicate SNP %s with the same alleles found in file [%s]. This row was skiped:\n", vs_buf[1].c_str(),esdpath.c_str());
                    printf("%s.\n",buf);
                }
            }
        }
        
        if(gzflag) gzclose(gzfile);
        else flptr.close();
        cout << lineNum << " SNPs infomation to be included from [" + string(esdpath) + "]." << endl;
    }
    void read_plink_qassoc_gz(vector<string> &rs,vector<int> &chr,vector<int> &bp, vector<float> &beta,vector<float> &se,string esdpath){
        
        cout << "Reading summary information from [" + esdpath + "]." << endl;
        char tbuf[MAX_LINE_SIZE];
        gzFile gzfile = gzopen(esdpath.c_str(), "rb");
        if (!(gzfile)) {
            fprintf (stderr, "%s: Couldn't open file %s\n",
                     esdpath.c_str(), strerror (errno));
            exit (EXIT_FAILURE);
        }
        int lineNum(0);
        gzgets(gzfile, tbuf, MAX_LINE_SIZE); //head
        vector<string> vs_buf;
        int col_num = split_string(tbuf, vs_buf, " \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"CHR\"");
            exit(EXIT_FAILURE);
        }
        
        while(!gzeof(gzfile))
        {
            gzgets(gzfile, tbuf, MAX_LINE_SIZE);
            if(tbuf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(tbuf, vs_buf, " \t\n");
                if(col_num!=9) {
                    printf("ERROR: column number is not right in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: \'NA\' of chromosome found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                int tmpchr;
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[0].c_str());
                chr.push_back(tmpchr);
                rs.push_back(vs_buf[1]);
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: \'NA\' of SNP position found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }

                bp.push_back(atoi(vs_buf[2].c_str()));
                 to_upper(vs_buf[4]);
                 to_upper(vs_buf[5]);
                if(vs_buf[4].compare("NA") && vs_buf[5].compare("NA"))
                {
                    beta.push_back(atof(vs_buf[4].c_str()));
                    se.push_back(atof(vs_buf[5].c_str()));
                } else {
                    beta.push_back(-9);
                    se.push_back(-9);
                }
                lineNum++;
            }
        }
        gzclose(gzfile);
        cout << lineNum << " SNPs summary info to be included from [" + esdpath + "]." << endl;
    }
    
    void read_plink_qassoc(vector<string> &rs,vector<int> &chr, vector<int> &bp, vector<float> &beta,vector<float> &se,string esdpath){
        
        ifstream flptr(esdpath.c_str());
        if (!flptr) throw ("Error: can not open the file [" + string(esdpath) + "] to read.");
        cout << "Reading eQTL information from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        flptr.getline(buf,MAX_LINE_SIZE);
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, " \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"CHR\"");
            exit(EXIT_FAILURE);
        }
        while(!flptr.eof())
        {
            flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, " \t\n");
                if(col_num!=9) {
                    printf("ERROR: column number is not right in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: \'NA\' of chromosome found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                int tmpchr;
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[0].c_str());
                chr.push_back(tmpchr);
                rs.push_back(vs_buf[1]);
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: \'NA\' of SNP position found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                bp.push_back(atoi(vs_buf[2].c_str()));
                 to_upper(vs_buf[4]);
                 to_upper(vs_buf[5]);
                if(vs_buf[4].compare("NA") && vs_buf[5].compare("NA"))
                {
                    beta.push_back(atof(vs_buf[4].c_str()));
                    se.push_back(atof(vs_buf[5].c_str()));
                } else {
                    beta.push_back(-9);
                    se.push_back(-9);
                }
                lineNum++;
            }
        }
        flptr.close();
        cout << lineNum << " SNPs summary info to be included from [" + esdpath + "]." << endl;
    }
    void read_plink_qassoc(vector<snpinfolst> &snpinfo,string esdpath)
    {
        gzFile gzfile=NULL;
        ifstream flptr;
        bool gzflag=has_suffix(esdpath, "gz");
        if(gzflag)
        {
            gzfile = gzopen(esdpath.c_str(), "rb");
            if (!(gzfile)) {
                fprintf (stderr, "%s: Couldn't open file %s\n",
                         esdpath.c_str(), strerror (errno));
                exit (EXIT_FAILURE);
            }
        }
        else
        {
            flptr.open(esdpath.c_str());
            if (!flptr) throw ("Error: can not open the file [" + esdpath + "] to read.");
        }
        cout << "Reading eQTL information from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
        else flptr.getline(buf,MAX_LINE_SIZE);// the header
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, " \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"CHR\"");
            exit(EXIT_FAILURE);
        }
        while(!flptr.eof())
        {
            snpinfolst tmpesd;
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
            else flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, " \t\n");
                if(col_num!=9) {
                    printf("ERROR: column number is not right in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: \'NA\' of chromosome found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpesd.snpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpesd.snpchr=24;
                else tmpesd.snpchr=atoi(vs_buf[0].c_str());
                
                tmpesd.snprs=vs_buf[1].c_str();
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: \'NA\' of SNP position found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                tmpesd.bp=atoi(vs_buf[2].c_str());
                 to_upper(vs_buf[4]);
                 to_upper(vs_buf[5]);
                if(vs_buf[4].compare("NA") && vs_buf[5].compare("NA"))
                {
                   tmpesd.beta=atof(vs_buf[4].c_str());
                   tmpesd.se=atof(vs_buf[5].c_str());
                } else {
                   tmpesd.beta=-9;
                   tmpesd.se=-9;
                }
                 snpinfo.push_back(tmpesd);
                lineNum++;
            }
        }
        if(gzflag) gzclose(gzfile);
        else flptr.close();
        cout << lineNum << " SNPs summary info to be included from [" + esdpath + "]." << endl;
    }
    void read_gemma_qassoc(vector<string> &rs,vector<int> &chr, vector<int> &bp, vector<string> &a1, vector<string> &a2, vector<float> &beta,vector<float> &se,string esdpath)
    {
        gzFile gzfile=NULL;
        ifstream flptr;
        bool gzflag=has_suffix(esdpath, "gz");
        if(gzflag)
        {
            gzfile = gzopen(esdpath.c_str(), "rb");
            if (!(gzfile)) {
                fprintf (stderr, "%s: Couldn't open file %s\n",
                         esdpath.c_str(), strerror (errno));
                exit (EXIT_FAILURE);
            }
        }
        else
        {
            flptr.open(esdpath.c_str());
            if (!flptr) throw ("Error: can not open the file [" + esdpath + "] to read.");
        }
        cout << "Reading eQTL information from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
        else flptr.getline(buf,MAX_LINE_SIZE);// the header
        /* check headers */
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, " \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"chr\"");
            exit(EXIT_FAILURE);
        }
        
        while(!flptr.eof() && !gzeof(gzfile))
        {
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
            else flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, " \t\n");
                if(col_num!=11) {
                    printf("ERROR: column number is not right in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: \'NA\' of chromosome found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                int tmpchr;
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                else tmpchr=atoi(vs_buf[0].c_str());
                chr.push_back(tmpchr);
                
                rs.push_back(vs_buf[1]);
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: \'NA\' of SNP position found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                bp.push_back(atoi(vs_buf[2].c_str()));
                 to_upper(vs_buf[4]);
                a1.push_back(vs_buf[4]);
                 to_upper(vs_buf[5]);
                a2.push_back(vs_buf[5]);
                 to_upper(vs_buf[7]);
                 to_upper(vs_buf[8]);
                if(vs_buf[7].compare("NA") && vs_buf[8].compare("NA")) {
                    beta.push_back(atof(vs_buf[7].c_str()));
                    se.push_back(atof(vs_buf[8].c_str()));
                } else {
                    beta.push_back(-9);
                    se.push_back(-9);
                }
                lineNum++;
                
            }
        }
        
        cout << lineNum << " SNPs infomation to be included from [" + string(esdpath) + "]." << endl;
        if(gzflag) gzclose(gzfile);
        else flptr.close();
    }
    void read_gemma_qassoc(vector<snpinfolst> &snpinfo,string esdpath)
    {
        gzFile gzfile=NULL;
        ifstream flptr;
        bool gzflag=has_suffix(esdpath, "gz");
        if(gzflag)
        {
            gzfile = gzopen(esdpath.c_str(), "rb");
            if (!(gzfile)) {
                fprintf (stderr, "%s: Couldn't open file %s\n",
                         esdpath.c_str(), strerror (errno));
                exit (EXIT_FAILURE);
            }
        }
        else
        {
            flptr.open(esdpath.c_str());
            if (!flptr) throw ("Error: can not open the file [" + esdpath + "] to read.");
        }
        cout << "Reading eQTL information from [" + string(esdpath) + "]." << endl;
        
        char buf[MAX_LINE_SIZE];
        int lineNum(0);
        if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
        else flptr.getline(buf,MAX_LINE_SIZE);// the header
        /* check headers */
        vector<string> vs_buf;
        int col_num = split_string(buf, vs_buf, " \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="CHR") {
            printf("ERROR: the headers should start with \"chr\"");
            exit(EXIT_FAILURE);
        }
        
        while(!flptr.eof() && !gzeof(gzfile))
        {
            snpinfolst tmpesd;
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
            else flptr.getline(buf,MAX_LINE_SIZE);
            if(buf[0]!='\0'){
                vs_buf.clear();
                int col_num = split_string(buf, vs_buf, " \t\n");
                if(col_num!=11) {
                    printf("ERROR: column number is not right in row %d.!", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                    printf("ERROR: \'NA\' of chromosome found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpesd.snpchr=23;
                else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpesd.snpchr=24;
                else tmpesd.snpchr=atoi(vs_buf[0].c_str());
                
                tmpesd.snprs=vs_buf[1].c_str();
                if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                    printf("ERROR: \'NA\' of SNP position found in row %d.\n", lineNum+2);
                    exit(EXIT_FAILURE);
                }
                tmpesd.bp=atoi(vs_buf[2].c_str());
                 to_upper(vs_buf[4]);
                tmpesd.a1=vs_buf[4].c_str();
                 to_upper(vs_buf[5]);
                tmpesd.a2=vs_buf[5].c_str();
                 to_upper(vs_buf[7]);
                 to_upper(vs_buf[8]);
                if(vs_buf[7].compare("NA") && vs_buf[8].compare("NA")) {
                    tmpesd.beta=atof(vs_buf[7].c_str());
                    tmpesd.se=atof(vs_buf[8].c_str());
                } else {
                    tmpesd.beta=-9;
                    tmpesd.se=-9;
                }
                snpinfo.push_back(tmpesd);
                lineNum++;
            }
        }
        
        cout << lineNum << " SNPs infomation to be included from [" + string(esdpath) + "]." << endl;
        if(gzflag) gzclose(gzfile);
        else flptr.close();
    }
    uint64_t get_esi_info(vector<snpinfolst> &snpinfo, probeinfolst** prbiflst,int prb_num,int fformat)
    {
        probeinfolst* locinfolst=*prbiflst;
        map<string, int> rs_map;
        //map<string, int> rsbp_map;
        //map<string, int> rsaa_map;
        long size = 0;
        uint64_t ttl_value_num =0;
        char buf[MAX_LINE_SIZE];
        for(int j=0;j<prb_num;j++)
        {
            printf("Scanning... %3.0f%%\r", 100.0*j/prb_num);
            fflush(stdout);
            //map<string, int> currs_map;
            //long currssize=0;
            string esdfilename=(locinfolst+j)->esdpath;
            gzFile gzfile = NULL;
            ifstream flptr;
            bool gzflag=has_suffix(esdfilename, "gz");
            if(gzflag)
            {
                gzfile = gzopen(esdfilename.c_str(), "rb");
                if (!(gzfile)) {
                    fprintf (stderr, "%s: Couldn't open file %s\n",
                             esdfilename.c_str(), strerror (errno));
                    exit (EXIT_FAILURE);
                }
            }
            else
            {
                flptr.open((locinfolst+j)->esdpath.c_str());
                if (!flptr) throw ("Error: can not open the file [" + string((locinfolst+j)->esdpath) + "] to read.");
            }
            
            if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE); //head
            else flptr.getline(buf,MAX_LINE_SIZE);// the header
            vector<string> vs_buf;
            int col_num = split_string(buf, vs_buf, " \t\n");
            to_upper(vs_buf[0]);
            if(vs_buf[0]!="CHR") {
                printf("ERROR: the headers of file %s should be start with \"Chr\".\n",esdfilename.c_str());
                exit(EXIT_FAILURE);
            }
            
            if(fformat==0)
            {
                int lineNum=0;
                while(!flptr.eof() && !gzeof(gzfile))
                {
                    if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
                    else flptr.getline(buf,MAX_LINE_SIZE);
                    if(buf[0]!='\0')
                    {
                        ttl_value_num++;
                        vs_buf.clear();
                        col_num = split_string(buf, vs_buf, " \t\n");
                        if(col_num!=9) {
                            printf("ERROR: column number is not right of row %d in esd file \"%s\"!\n", lineNum+2, esdfilename.c_str());
                            exit(EXIT_FAILURE);
                        }
                        if(vs_buf[0]=="NA" ||vs_buf[0]=="na" || vs_buf[2]=="NA" ||vs_buf[2]=="na")
                        {
                            printf("ERROR: chromosome or SNP position of row %d can't be \"NA\" in esd file \"%s\"!\n", lineNum+2, esdfilename.c_str());
                            exit(EXIT_FAILURE);
                        }
                        /*
                        currs_map.insert(pair<string, int>(vs_buf[1].c_str(), lineNum));
                        if(currs_map.size()==currssize)
                        {
                            printf("WARING: duplicate SNP %s found in file %s.\n", vs_buf[1].c_str(),esdfilename.c_str());
                        } else {
                            currssize=currs_map.size();
                        }
                         */
                         to_upper(vs_buf[3]);
                         to_upper(vs_buf[4]);
                        if(vs_buf[3]=="NA" ||vs_buf[4]=="NA")
                        {
                            printf("WARNING: Allele infomation of row %d is \"NA\" in esd file \"%s\"! This may cause incorrect positive/negative sign of effect size when conducting allele check.\n", lineNum+2, esdfilename.c_str());
                           
                        }
                        string crsstr=vs_buf[0]+":"+vs_buf[1];
                        //string crsbpstr=vs_buf[0]+":"+vs_buf[1]+":"+vs_buf[2];
                        //string crsbpaastr=vs_buf[0]+":"+vs_buf[1]+":"+vs_buf[2]+":"+vs_buf[3]+":"+vs_buf[4];
                        //string crsbpaastr_iv=vs_buf[0]+":"+vs_buf[1]+":"+vs_buf[2]+":"+vs_buf[4]+":"+vs_buf[3];
                        rs_map.insert(pair<string, int>(crsstr.c_str(), lineNum));
                        //rsbp_map.insert(pair<string, int>(crsbpstr.c_str(), lineNum));
                        //rsaa_map.insert(pair<string, int>(crsbpaastr.c_str(), lineNum));
                        //if(rs_map.size() != rsbp_map.size())
                        //{
                        //    printf("ERROR: SNP %s on Chromosome %s has multiple BPs, please check.\n", vs_buf[1].c_str(), vs_buf[0].c_str());
                        //    exit(EXIT_FAILURE);
                        //}
                        bool newsnp=false;
                        if (size < rs_map.size())
                        {
                            newsnp=true;
                            size = rs_map.size();
                        }/*else {
                            //check multi-allelic SNPs
                            if(rs_map.size() != rsaa_map.size())
                            {
                                long tmpsize=rsaa_map.size();
                                rsaa_map.insert(pair<string, int>(crsbpaastr_iv.c_str(), lineNum));
                                if(tmpsize==rsaa_map.size())
                                {
                                    printf("WARING: switched ref allele with alt allele of SNP %s found.\n", vs_buf[1].c_str());
                                } else {
                                    newsnp=true;
                                    printf("WARING: multi-allelic SNPs with duplicate SNP ID %s found.\n", vs_buf[1].c_str());
                                    rs_map.insert(pair<string, int>((crsbpaastr+"_m").c_str(), lineNum));
                                    rsbp_map.insert(pair<string, int>((crsbpaastr+"_m").c_str(), lineNum));
                                }
                                rs_map.insert(pair<string, int>((crsbpaastr+"_").c_str(), lineNum));
                                rsbp_map.insert(pair<string, int>((crsbpaastr+"_").c_str(), lineNum));
                                size = rs_map.size();
                            }
                        }*/
                        if(newsnp)
                        {
                            snpinfolst snptmp;
                            int tmpchr;
                            if(vs_buf[0]=="X" || vs_buf[0]=="x") tmpchr=23;
                            else if(vs_buf[0]=="Y" || vs_buf[0]=="y") tmpchr=24;
                            else tmpchr=atoi(vs_buf[0].c_str());
                            snptmp.snpchr=tmpchr;
                            snptmp.snprs=vs_buf[1].c_str();
                            snptmp.bp=atoi(vs_buf[2].c_str());
                            snptmp.a1=vs_buf[3].c_str();
                            snptmp.a2=vs_buf[4].c_str();
                            snptmp.gd=0;
                            snpinfo.push_back(snptmp);
                        }
                        lineNum++;
                    }
                }
            }
            else if(fformat==1)
            {
                bInfo binfo;
                string bimfname=string((locinfolst+j)->bfilepath)+".bim";
                read_bimfile(&binfo, bimfname);
                
                while(!flptr.eof() && !gzeof(gzfile))
                {
                    if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
                    else flptr.getline(buf,MAX_LINE_SIZE);
                    if(buf[0]!='\0')
                    {
                        string tmpStr;
                        int tmpchr;
                        snpinfolst snptmp;
                        ttl_value_num++;
                        istringstream iss(buf);
                        iss>>tmpStr; //chr
                        if(tmpStr=="NA" ||tmpStr=="na" )
                        {
                            printf("ERROR: chromosome  can't be \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                            exit(EXIT_FAILURE);
                        }
                        if(tmpStr=="X" || tmpStr=="x") tmpchr=23;
                        else if(tmpStr=="Y" || tmpStr=="y") tmpchr=24;
                        else tmpchr=atoi(tmpStr.c_str());
                        iss>>tmpStr; //rs
                        if(tmpStr[0]!='\0') rs_map.insert(pair<string, int>(tmpStr.c_str(), tmpchr));
                        if (size < rs_map.size())
                        {
                            snptmp.snpchr=tmpchr;
                            snptmp.snprs=tmpStr.c_str();
                            map<string, int>::iterator iter;
                            iter=binfo._snp_name_map.find(tmpStr);
                            if( iter != binfo._snp_name_map.end())
                            {
                                snptmp.a1=binfo._allele1[iter->second];
                                snptmp.a2=binfo._allele2[iter->second];
                            }else {
                                printf("ERROR: can't find SNP %s in the bim file %s.", tmpStr.c_str(), bimfname.c_str());
                                exit(EXIT_FAILURE);
                            }
                            iss>>tmpStr; //BP
                            snptmp.bp=atoi(tmpStr.c_str());
                            snptmp.gd=0;
                            size = rs_map.size();
                            snpinfo.push_back(snptmp);
                        }
                    }
                }
            }
            else if(fformat==2)
            {
                while(!flptr.eof() && !gzeof(gzfile) )
                {
                    if(gzflag) gzgets(gzfile, buf, MAX_LINE_SIZE);
                    else flptr.getline(buf,MAX_LINE_SIZE);
                    if(buf[0]!='\0')
                    {
                        string tmpStr;
                        int tmpchr;
                        snpinfolst snptmp;
                        ttl_value_num++;
                        istringstream iss(buf);
                        iss>>tmpStr; //chr
                        if(tmpStr=="NA" ||tmpStr=="na" )
                        {
                            printf("ERROR: chromosome can't be \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                            exit(EXIT_FAILURE);
                        }
                        if(tmpStr=="X" || tmpStr=="x") tmpchr=23;
                        else if(tmpStr=="Y" || tmpStr=="y") tmpchr=24;
                        else tmpchr=atoi(tmpStr.c_str());
                        iss>>tmpStr; //rs
                        if(tmpStr[0]!='\0') rs_map.insert(pair<string, int>(tmpStr.c_str(), tmpchr));
                        if (size < rs_map.size())
                        {
                            snptmp.snpchr=tmpchr;
                            snptmp.snprs=tmpStr.c_str();
                            iss>>tmpStr; //BP
                            if(tmpStr=="NA" ||tmpStr=="na" )
                            {
                                printf("ERROR: SNP position can't be \"NA\" in  file \"%s\"!\n", esdfilename.c_str());
                                exit(EXIT_FAILURE);
                            }
                            snptmp.bp=atoi(tmpStr.c_str());
                            iss>>tmpStr;//n_miss
                            iss>>tmpStr; //A1 minor allele
                             to_upper(tmpStr);
                            snptmp.a1=tmpStr.c_str();
                            iss>>tmpStr; //A0 major allele
                             to_upper(tmpStr);
                            snptmp.a2=tmpStr.c_str();
                            snptmp.gd=0;
                            size = rs_map.size();
                            snpinfo.push_back(snptmp);
                        }
                    }
                }
            }
            
        }
        return ttl_value_num;
    }
    void save_txts_dbesd(char* outFileName, long esiNum, long epiNum,vector<int> &epi2esd,probeinfolst* prbiflst, int fformat,vector<string> &esi_rs,vector<string> &esi_a1,vector<string> &esi_a2)
    {
        // get esd info
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
        map<string, int> esi_map;
        for(int j=0;j<esi_rs.size();j++)
        {
            esi_map.insert(pair<string,int>(esi_rs[j],j));
        }
        
        map<string, int>::iterator iter;
        for(int j=0;j<epiNum;j++)
        {
            printf("Saving... %3.0f%%\r", 100.0*j/epiNum);
            fflush(stdout);
            int esdnumcurprb=epi2esd[j+1]-epi2esd[j];
            printf("Reading %d text file(s) of probe %s.\n", esdnumcurprb, prbiflst[epi2esd[j]].probeId.c_str());
            for(int k=0;k<bsize;k++) buffer[k]=-9; //init
            for(int k=epi2esd[j];k<epi2esd[j+1];k++)
            {
                vector<string> _rs;
                vector<float> _beta;
                vector<float> _se;
                vector<int> _chr;
                vector<string> _a1;
                vector<string> _a2;
                vector<int> _bp;
                string zname=prbiflst[k].esdpath;
                bool gzflag=has_suffix(zname, "gz");
                switch (fformat)
                {
                    case 0:
                        read_smr_sa(_rs,_chr,_bp, _a1, _a2, _beta, _se, prbiflst[k].esdpath);
                        break;
                    case 1:
                        if(gzflag) read_plink_qassoc_gz(_rs, _chr,_bp, _beta, _se, prbiflst[k].esdpath);
                        else read_plink_qassoc(_rs,_chr, _bp,_beta, _se, prbiflst[k].esdpath);
                        break;
                    case 2:
                        read_gemma_qassoc(_rs,_chr, _bp,_a1, _a2, _beta, _se, prbiflst[k].esdpath);
                        break;
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
                    if(fformat==0)
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
                    } else {
                        buffer[rsid[l]]=_beta[l];
                        buffer[esiNum+rsid[l]]=_se[l];
                    }
                }
            }
            fwrite (buffer,sizeof(float), bsize, smr1);
        }
        fclose (smr1);
        free(buffer);
        cout<<"Beta values and SE values for "<<epiNum<<" Probes and "<<esiNum<<" SNPs have been saved in the dense binary file [" + esdfile + "]." <<endl;

    }
    
     void save_full_txts_sbesd(char* outFileName, long esiNum, long epiNum,vector<int> &epi2esd,probeinfolst* prbiflst, int fformat,vector<string> &esi_rs,vector<string> &esi_a1,vector<string> &esi_a2)
    {
        map<string, int> esi_map;
        for(int j=0;j<esi_rs.size();j++)
        {
            esi_map.insert(pair<string,int>(esi_rs[j],j));
        }
        
        // get esd info
        
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
        
        map<string, int>::iterator iter;
        for(int j=0;j<epiNum;j++)
        {
            printf("Saving... %3.0f%%\r", 100.0*j/epiNum);
            fflush(stdout);
            int esdnumcurprb=epi2esd[j+1]-epi2esd[j];
            printf("\nReading %d text file(s) of probe %s.\n", esdnumcurprb, prbiflst[epi2esd[j]].probeId.c_str());
            
            map<string, int> rsa_map;
            long rsNum=0;
            
            vector<uint32_t> tmprid;
            vector<float> tmpse;
            for(int k=epi2esd[j];k<epi2esd[j+1];k++)
            {
                vector<string> _rs;
                vector<float> _beta;
                vector<float> _se;
                vector<int> _chr;
                vector<string> _a1;
                vector<string> _a2;
                vector<int> _bp;
                string zname=prbiflst[k].esdpath;
                bool gzflag=has_suffix(zname, "gz");
                switch (fformat)
                {
                    case 0:
                        read_smr_sa(_rs,_chr,_bp, _a1, _a2, _beta, _se, prbiflst[k].esdpath);
                        break;
                    case 1:
                        if(gzflag) read_plink_qassoc_gz(_rs, _chr,_bp, _beta, _se, prbiflst[k].esdpath);
                        else read_plink_qassoc(_rs,_chr, _bp,_beta, _se, prbiflst[k].esdpath);
                        break;
                    case 2:
                        read_gemma_qassoc(_rs, _chr, _bp, _a1, _a2, _beta, _se, prbiflst[k].esdpath);
                        break;
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
                        string chckstr=(fformat==1?_rs[l]:(_rs[l]+":"+_a1[l]+":"+_a2[l]));
                        rsa_map.insert(pair<string,int>(chckstr,l));
                        if(rsNum<rsa_map.size())
                        {
                            if(fformat==0)
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
                            } else {
                                    val.push_back(_beta[l]);
                                    rowids.push_back(rsid[l]);
                                    tmpse.push_back(_se[l]);
                                    tmprid.push_back(rsid[l]);
                            }
                            rsNum=rsa_map.size();
                        } else {
                            printf("WARNING: duplicate SNP %s with the same alleles belonging to the same probe %s found in current esd file\"%s\" and the esd file(s) read before. \n",_rs[l].c_str(),prbiflst[k].probeId.c_str(),zname.c_str());
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

    void slct_sparse_per_prb(vector<int> &slct_idx, probeinfolst* prbifo, vector<snpinfolst> &snpinfo, long cis_itvl, long trans_itvl,double transThres,double restThres,FILE* logfile)
    {
        
        // for log file
        map<string, int> rsa_map;
        map<string, int>::iterator iter;
        long rsNum=0;

        vector<int> tran_chr;
        vector<uint64_t> lowerBp;
        vector<uint64_t> upperBp;
        vector<uint64_t> nsnp;
        vector<bool> extend;
        vector<bool> merge;
        vector<int> cis_idx;
        map<string, int> other_idx;
        
        //extract info
        string probid=prbifo->probeId;
        long probbp=prbifo->bp;
        long probchr=prbifo->probechr;
        long cisuperBounder= probbp+cis_itvl;
        long cislowerBounder=((probbp-cis_itvl>0)?(probbp-cis_itvl):0);
       
        printf("Extracting information of probe %s...\n", probid.c_str());
        for(int l=0;l<snpinfo.size();l++)
        {
            if(abs(snpinfo[l].se+9)>1e-6)
            {
                double zsxz=snpinfo[l].beta/snpinfo[l].se;
                double pxz=pchisq(zsxz*zsxz, 1);
                
                if(snpinfo[l].snpchr == probchr && snpinfo[l].bp<=cisuperBounder && snpinfo[l].bp>=cislowerBounder)
                {
                    string chckstr=snpinfo[l].snprs+":"+snpinfo[l].a1+":"+snpinfo[l].a2;
                    rsa_map.insert(pair<string,int>(chckstr,rsNum));
                    if(rsNum<rsa_map.size()){
                        slct_idx.push_back(l);
                        rsNum=rsa_map.size();
                        cis_idx.push_back(l);
                    }
                    
                }
                else if(pxz<=transThres)
                {
                    uint64_t transNum=0;
                    int curChr=snpinfo[l].snpchr;
                    string chckstr=snpinfo[l].snprs+":"+snpinfo[l].a1+":"+snpinfo[l].a2;
                    rsa_map.insert(pair<string, int>(chckstr,rsNum));
                    if(rsNum<rsa_map.size()){
                        slct_idx.push_back(l);
                        rsNum=rsa_map.size();
                        transNum++;
                    }
                    
                    long transbp=snpinfo[l].bp;
                    long translowerBounder=((transbp-trans_itvl>0)?(transbp-trans_itvl):0);
                    long transuperBounder=transbp+trans_itvl;
                    bool extended=false;
                    bool merged=false;
                    
                    int startptr=l-1;
                    while(startptr>=0 && curChr == snpinfo[startptr].snpchr && transbp-snpinfo[startptr].bp<=trans_itvl && abs(snpinfo[startptr].se+9)>1e-6)
                    {
                        translowerBounder=snpinfo[startptr].bp;
                        if(upperBp.size()>0 && translowerBounder<=upperBp[upperBp.size()-1] && tran_chr[upperBp.size()-1]==snpinfo[startptr].snpchr) //trans region merges
                        {
                            merged=true;
                            transNum=transNum+nsnp[nsnp.size()-1];
                            break;
                        }
                        if(snpinfo[startptr].snpchr == probchr && translowerBounder<=cisuperBounder && translowerBounder>=cislowerBounder) // trans touches cis region
                        {
                            break;
                        }
                        string chckstr=snpinfo[startptr].snprs+":"+snpinfo[startptr].a1+":"+snpinfo[startptr].a2;
                        iter=other_idx.find(chckstr);
                        if(iter!=other_idx.end()) // trans merges other eqtl
                        {
                            other_idx.erase(iter->first);
                            transNum++;
                        } else {
                            rsa_map.insert(pair<string, int>(chckstr,rsNum));
                            if(rsNum<rsa_map.size()){
                                slct_idx.push_back(startptr);
                                rsNum=rsa_map.size();
                                transNum++;
                            }
                        }
                        
                        startptr--;
                    }
                    startptr=l+1;
                    while (startptr<snpinfo.size() && curChr == snpinfo[startptr].snpchr && snpinfo[startptr].bp - transbp <= trans_itvl && abs(snpinfo[startptr].se+9)>1e-6)
                    {
                        transuperBounder=snpinfo[startptr].bp;
                        if(snpinfo[startptr].snpchr == probchr && transuperBounder>=cislowerBounder && transuperBounder<=cisuperBounder) // trans touches cis region
                        {
                            break;
                            
                        }  else {
                            
                            string chckstr=snpinfo[startptr].snprs+":"+snpinfo[startptr].a1+":"+snpinfo[startptr].a2;
                            rsa_map.insert(pair<string, int>(chckstr,rsNum));
                            if(rsNum<rsa_map.size()){
                                slct_idx.push_back(startptr);
                                rsNum=rsa_map.size();
                                transNum++;
                                double zsxz_tmp=snpinfo[startptr].beta/snpinfo[startptr].se;
                                double pxz_tmp=pchisq(zsxz_tmp*zsxz_tmp, 1);
                                if(pxz_tmp<transThres) // trans region extends
                                {
                                    extended=true;
                                    transbp=snpinfo[startptr].bp;
                                }
                                l=startptr;
                            }
                        }
                        startptr++;
                    }
                    if(merged)
                    {
                        upperBp[upperBp.size()-1]=transuperBounder;
                        nsnp[nsnp.size()-1]=transNum;
                        extend[nsnp.size()-1]=extended;
                        merge[nsnp.size()-1]=merged;
                    }
                    else
                    {
                        tran_chr.push_back(curChr);
                        upperBp.push_back(transuperBounder);
                        lowerBp.push_back(translowerBounder);
                        nsnp.push_back(transNum);
                        extend.push_back(extended);
                        merge.push_back(merged);
                    }
                    
                }
                else if(pxz<restThres)
                {
                    string chckstr=snpinfo[l].snprs+":"+snpinfo[l].a1+":"+snpinfo[l].a2;
                    rsa_map.insert(pair<string,int>(chckstr,rsNum));
                    if(rsNum<rsa_map.size()){
                        slct_idx.push_back(l);
                        rsNum=rsa_map.size();
                        other_idx.insert(pair<string,int>(chckstr,l));
                    }
                }
            } else {
                printf("WARNING: SNP %s  with \"NA\" value found and skipped.\n",snpinfo[l].snprs.c_str());
            }
        }
        long cissnpnum=cis_idx.size();
        long transnum=nsnp.size();
        long transsnpnum=0;
        long othersnpnum=other_idx.size();
        for(int l=0;l<nsnp.size();l++) transsnpnum+=nsnp[l];
        if(slct_idx.size()-cissnpnum-transsnpnum-othersnpnum != 0)
        {
            printf("Something is wrong with this selection methold. Please report this bug.\n");
            exit(EXIT_FAILURE);
        }
        printf("%ld SNPs in the cis-region, %ld SNPs in total %ld trans-region(s), %ld SNPs beyond cis-region and trans-region(s) have been extracted.\n",cissnpnum,transsnpnum,transnum,othersnpnum);
        //log
       
        if(cis_idx.size() || nsnp.size() || other_idx.size())
        {
            string logstr="{"+probid+","+atos(probchr)+","+atos(probbp)+"}\t";
            if(cis_idx.size()>0)
            {                
                logstr+= "["+atos(probchr)+","+ atos(snpinfo[cis_idx[0]].bp)+","+atos(snpinfo[cis_idx[cis_idx.size()-1]].bp)+","+atos(cis_idx.size())+"]\t";
            }else logstr+="[]\t";
            if(nsnp.size()>0)
            {
                for(int h=0;h<nsnp.size();h++)
                {
                    logstr+= "<"+atos(tran_chr[h])+","+ atos(lowerBp[h])+","+atos(upperBp[h])+","+atos(nsnp[h])+">\t";
                }
            }else logstr+="<>\t";
            
            if(other_idx.size()>0)
            {
                //for(iter=other_idx.begin();iter!=other_idx.end();iter++) logstr+="("+atos(snpinfo[iter->second].bp)+"), ";
                //logstr+="\n";
                logstr+="("+atos(other_idx.size())+")\n";
            }
            else logstr+="()\n";
            
            fputs(logstr.c_str(),logfile);
            fflush(logfile);
        }
    }
    void save_slct_txts_sbesd(char* outFileName, long esiNum, long epiNum,vector<int> &epi2esd,probeinfolst* prbiflst, int fformat,vector<string> &esi_rs,vector<string> &esi_a1,vector<string> &esi_a2,int cis_itvl, int trans_itvl, float transThres, float restThres)
    {
        map<string, int> esi_map;
        for(int j=0;j<esi_rs.size();j++)
        {
            esi_map.insert(pair<string,int>(esi_rs[j],j));
        }
        
        // get esd info
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
            vector<uint32_t> tmprid;
            vector<float> tmpse;
            int esdnumcurprb=epi2esd[j+1]-epi2esd[j];
            printf("\nReading %d text file(s) of probe %s.\n", esdnumcurprb, prbiflst[epi2esd[j]].probeId.c_str());
            snpinfo.clear();
            for(int k=epi2esd[j];k<epi2esd[j+1];k++)
            {
                string zname=prbiflst[k].esdpath;
                switch (fformat)
                {
                    case 0:
                        read_smr_sa(snpinfo, prbiflst[k].esdpath);
                        break;
                    case 1:
                        read_plink_qassoc(snpinfo, prbiflst[k].esdpath);
                        break;
                    case 2:
                        read_gemma_qassoc(snpinfo, prbiflst[k].esdpath);
                        break;
                }
            }
            snpinfolst* sortptr=&snpinfo[0];
            qsort(sortptr,snpinfo.size(),sizeof(snpinfolst),comp_esi);
            
            probeinfolst prbifo=prbiflst[epi2esd[j]];
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
                        if(fformat==0)
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
                        } else {
                            val.push_back(_beta[l]);
                            rowids.push_back(rsid[l]);
                            tmpse.push_back(_se[l]);
                            tmprid.push_back(rsid[l]);
                        }
                        rsNum=rsa_map.size();
                    } else {
                        printf("WARNING: duplicate SNP %s with the same alleles belonging to the same probe %s. \n",_rs[l].c_str(),prbiflst[epi2esd[j]].probeId.c_str());
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
    void make_besd(char*outFileName, char* syllabusName, bool gctaflag,bool plinkflag,bool gemmaflag,bool merlinflag,bool save_dense_flag,int cis_itvl, int trans_itvl, float transThres, float restThres)
    {
        if(syllabusName==NULL) throw("Error: please input eQTL syllabus file by the flag --eqtl-flist.");
        int flagcount=0, fformat=0;
        if(gctaflag){ flagcount++; fformat=0;}
        if(plinkflag){ flagcount++;fformat=1;}
        if(gemmaflag){ flagcount++; fformat=2;}
        if(merlinflag){ flagcount++; fformat=3;}
        if(flagcount==0 || flagcount>1) throw("Error: please verify the file format flags. Only one flag can be specified.");
        vector<probeinfolst> prbiflstinfo;
        int esdfileNum=read_probeinfolst(prbiflstinfo,syllabusName);
        probeinfolst* prbiflst=&prbiflstinfo[0];
        qsort(prbiflst,esdfileNum,sizeof(probeinfolst),comp);
        
        //get epi
        printf("Generating epi file...\n");
        map<string, int> epi_map;
        long epiNum=0;
        string epifile = string(outFileName)+string(".epi");
        ofstream epi(epifile.c_str());
        if (!epi)
        {
            printf("Error: can not open the epi file %s to save!", epifile.c_str());
            exit(EXIT_FAILURE);
        }
        vector<int> epi2esd;
        for (int j = 0;j <esdfileNum; j++) {
            epi_map.insert(pair<string,int>(prbiflst[j].probeId,epiNum));
            if(epiNum<epi_map.size())
            {
                epi2esd.push_back(j);
                epi<<prbiflst[j].probechr<<'\t'<<prbiflst[j].probeId<<'\t'<<0<<'\t'<<prbiflst[+j].bp<<'\t'<<prbiflst[j].genename<<'\t'<<prbiflst[j].orien<<'\n';
                epiNum++;
            }
        }
        epi2esd.push_back((int)esdfileNum);
        epi.close();
        printf("%ld probes have been saved in the file \"%s\".\n",epiNum, epifile.c_str());
        
        //get esi: from bim file
        vector<snpinfolst> snpinfo;
        printf("\nScaning %d esd files to generate esi file...\n", esdfileNum);
        uint64_t ttlv = get_esi_info(snpinfo,&prbiflst,esdfileNum,fformat);
        long esiNum=snpinfo.size();
        snpinfolst* esiptr=&snpinfo[0];
        qsort(esiptr,esiNum,sizeof(snpinfolst),comp_esi);
        vector<string> esi_rs(esiNum);
        vector<string> esi_a1(esiNum);
        vector<string> esi_a2(esiNum);
        string esifile =  string(outFileName)+string(".esi");
        ofstream esi(esifile.c_str());
        if (!esi) throw ("Error: can not open the esi file to save!");
        for (int j = 0;j <esiNum; j++) {
            esi<<snpinfo[j].snpchr<<'\t'<<snpinfo[j].snprs<<'\t'<<snpinfo[j].gd<<'\t'<<snpinfo[j].bp<<'\t'<<snpinfo[j].a1<<'\t'<<snpinfo[j].a2<<'\n';
            esi_rs[j]=snpinfo[j].snprs;
            esi_a1[j]=snpinfo[j].a1;
            esi_a2[j]=snpinfo[j].a2;
        }
        esi.close();
        cout<<esiNum<<" SNPs have been saved in the file [" + esifile + "]."<<endl;
        printf("\nScaning %d esd files to generate besd file...\n", esdfileNum);
        if(save_dense_flag)
        {
            double sparsity=1.0*ttlv/(esiNum*epiNum);
            if(sparsity>=0.4)
            {
                printf("The sparsity of your data matrix is %f. We are going to save it in dense format!\n", sparsity);
                save_txts_dbesd(outFileName, esiNum, epiNum,epi2esd, prbiflst,fformat, esi_rs, esi_a1,esi_a2);
            } else {
                printf("The sparsity of your data matrix is %f. We are going to save it in sparse format!\n", sparsity);
                save_full_txts_sbesd( outFileName,  esiNum,  epiNum,epi2esd,prbiflst,fformat,esi_rs,esi_a1,esi_a2);
            }
            
        } else {
            save_slct_txts_sbesd(outFileName, esiNum, epiNum,epi2esd,prbiflst,fformat,esi_rs,esi_a1,esi_a2, cis_itvl,  trans_itvl,  transThres,  restThres);
        }        
    }
    
    void make_besd_byQfile(char* qfileName, char* outFileName,bool save_dense_flag,int cis_itvl, int trans_itvl, float transThres, float restThres)
    {
       
        FILE* qfile=fopen(qfileName, "r");
        if (!qfile) {
            printf("ERROR: Can't open file %s.\n ", qfileName);
            exit(EXIT_FAILURE);
        }
        
        printf("Reading data from %s ...\n", qfileName);
        vector<probeinfolst> prbiflst;
        vector<snpinfolst> snpinfo;
        vector<string> _ttl_rs;
        vector<string> _ttl_probe;
        vector<float> _ttl_beta;
        vector<float> _ttl_se;
        vector<string> _ttl_a1;
        vector<string> _ttl_a2;
        vector<int> _ttl_chr;
        vector<int> _ttl_bp;
        bool warningdone=false;
        char buf[MAX_LINE_SIZE];
        fgets(buf, MAX_LINE_SIZE, qfile); //header
        vector<string> vs_buf;
        split_string(buf, vs_buf, " \t\n");
        to_upper(vs_buf[0]);
        if(vs_buf[0]!="SNP") {
            printf("ERROR: the headers of query file should start with \"SNP\".\n");
            exit(EXIT_FAILURE);
        }
        long lineNum(0);
        long esiNum(0);
        long epiNum(0);
        snpinfolst tmpsnpinfo;
        probeinfolst tmpprbinfo;
        map<string, int> esi_map;
         long size=0;
        map<string, int> epi_map;
        map<string, int> rsbp_map;
        map<string, int> rsaa_map;
        map<string, int> probe_cp_map;
        map<string, int> probe_cpb_map;
        while(fgets(buf, MAX_LINE_SIZE, qfile))
        {
            vs_buf.clear();
            int col_num = split_string(buf, vs_buf, " \t\n");
            if(col_num!=13) {
                printf("ERROR: column number is not right in row %ld.\n", lineNum+2);
                exit(EXIT_FAILURE);
            }
            if(!save_dense_flag && (vs_buf[1]=="NA" || vs_buf[1]=="na"))
            {
                printf("ERROR: chromosome infomation of SNPs can't be \"NA\".\n");
                exit(EXIT_FAILURE);
            }
            if(!save_dense_flag && (vs_buf[6]=="NA" || vs_buf[6]=="na"))
            {
                printf("ERROR: chromosome infomation of probes can't be \"NA\". \n");
                exit(EXIT_FAILURE);
            }
            if(!save_dense_flag && (vs_buf[2]=="NA" || vs_buf[2]=="na"))
            {
                printf("ERROR: position of SNPs can't be \"NA\". \n");
                exit(EXIT_FAILURE);
            }
            if(!save_dense_flag && (vs_buf[7]=="NA" || vs_buf[7]=="na"))
            {
                printf("ERROR: position of probes can't be \"NA\". \n");
                exit(EXIT_FAILURE);
            }
            if((vs_buf[3]=="NA" || vs_buf[4]=="NA" || vs_buf[3]=="na" || vs_buf[4]=="na") && !warningdone)
            {
               printf("WARING: Allele information of one or more SNPs is \"NA\". This may cause incorrect positive/negative sign of effect size when conducting allele check.\n");
                warningdone=true;
               
            }
            int tmpchr;
            if(vs_buf[1]=="X" || vs_buf[1]=="x") tmpchr=23;
            else if(vs_buf[1]=="Y" || vs_buf[1]=="y") tmpchr=24;
            else tmpchr=atoi(vs_buf[1].c_str());
            to_upper(vs_buf[3]);
            to_upper(vs_buf[4]);
            to_upper(vs_buf[10]);
            to_upper(vs_buf[11]);
            if(!vs_buf[10].compare("NA") || !vs_buf[11].compare("NA")) {
                
                printf("WARNING: \"NA\" value of beta or se found, this row is skipped.\n");
                printf("%s",buf);
                continue;
            }
            string crsstr=vs_buf[0]+":"+vs_buf[1];
            string crsbpstr=vs_buf[0]+":"+vs_buf[1]+":"+vs_buf[2];
            string crsbpaastr=vs_buf[0]+":"+vs_buf[1]+":"+vs_buf[2]+":"+vs_buf[3]+":"+vs_buf[4];
            string crsbpaastr_iv=vs_buf[0]+":"+vs_buf[1]+":"+vs_buf[2]+":"+vs_buf[4]+":"+vs_buf[3];
            esi_map.insert(pair<string, int>(crsstr.c_str(), lineNum));
            rsbp_map.insert(pair<string, int>(crsbpstr.c_str(), lineNum));
            rsaa_map.insert(pair<string, int>(crsbpaastr.c_str(), lineNum));
            if(esi_map.size() != rsbp_map.size())
            {
                printf("ERROR: SNP %s on Chromosome %s has multiple BPs, please check.\n", vs_buf[0].c_str(), vs_buf[1].c_str());
                exit(EXIT_FAILURE);
            }
            bool newsnp=false;
            if (size < esi_map.size())
            {
                newsnp=true;
                size = esi_map.size();
            }else {
                //check multi-allelic SNPs
                if(esi_map.size() != rsaa_map.size())
                {
                    long tmpsize=rsaa_map.size();
                    rsaa_map.insert(pair<string, int>(crsbpaastr_iv.c_str(), lineNum));
                    if(tmpsize==rsaa_map.size())
                    {
                        printf("WARING: switched ref allele with alt allele of SNP %s found.\n", vs_buf[0].c_str());
                        vs_buf[10]=atos(-1.0*atof(vs_buf[10].c_str()));
                        string a0=vs_buf[3];
                        vs_buf[3]=vs_buf[4];
                        vs_buf[4]=a0;
                    } else {
                        newsnp=true;
                        printf("WARING: multi-allelic SNP %s found.\n", vs_buf[0].c_str());
                        esi_map.insert(pair<string, int>((crsbpaastr+"_m").c_str(), lineNum));
                        rsbp_map.insert(pair<string, int>((crsbpaastr+"_m").c_str(), lineNum));
                    }
                    esi_map.insert(pair<string, int>((crsbpaastr+"_").c_str(), lineNum));
                    rsbp_map.insert(pair<string, int>((crsbpaastr+"_").c_str(), lineNum));
                    size = esi_map.size();
                }
            }
            if(newsnp)
            {
                snpinfolst snptmp;
                snptmp.snpchr=tmpchr;
                snptmp.snprs=vs_buf[0].c_str();
                snptmp.bp=atoi(vs_buf[2].c_str());
                snptmp.a1=vs_buf[3].c_str();
                snptmp.a2=vs_buf[4].c_str();
                snptmp.gd=0;
                snpinfo.push_back(snptmp);
                esiNum++;
            }
            string cpstr=vs_buf[5]+":"+vs_buf[6];
            string cpbstr=vs_buf[5]+":"+vs_buf[6]+":"+vs_buf[7];
            probe_cp_map.insert(pair<string,int>(cpstr,lineNum));
            probe_cpb_map.insert(pair<string,int>(cpbstr,lineNum));
            if(probe_cp_map.size()!=probe_cpb_map.size())
            {
                printf("ERROR: Different BPs of probe %s found.\n",vs_buf[5].c_str() );
                exit(EXIT_FAILURE);
            }
            epi_map.insert(pair<string,int>(vs_buf[5],epiNum));
            if(epiNum<epi_map.size())
            {
                tmpprbinfo.probeId=vs_buf[5].c_str();
                tmpprbinfo.probechr=atoi(vs_buf[6].c_str());
                tmpprbinfo.bp=atoi(vs_buf[7].c_str());
                tmpprbinfo.genename=vs_buf[8].c_str();
                tmpprbinfo.orien=vs_buf[9][0];
                tmpprbinfo.gd=0;
                prbiflst.push_back(tmpprbinfo);
                epiNum++;
            }
            
            _ttl_rs.push_back(vs_buf[0].c_str());
            _ttl_a1.push_back(vs_buf[3].c_str());
            _ttl_a2.push_back(vs_buf[4].c_str());
            _ttl_probe.push_back(vs_buf[5]);
            _ttl_beta.push_back(atof(vs_buf[10].c_str()));
            _ttl_se.push_back(atof(vs_buf[11].c_str()));
            _ttl_chr.push_back(tmpchr);
            _ttl_bp.push_back(atoi(vs_buf[2].c_str()));
            lineNum++;
        }
        printf("%ld SNPs summary info to be included from %s.\n", lineNum, qfileName);
        fclose(qfile);
        
        printf("\nGenerating epi file...\n");
        probeinfolst* epiptr=&prbiflst[0];
        qsort(epiptr,epiNum,sizeof(probeinfolst),comp);
        string epifile = string(outFileName)+string(".epi");
        ofstream epi(epifile.c_str());
        if (!epi) throw ("Error: can not open the epi file " + epifile + " to save!");
        for (int j = 0;j <epiNum; j++) {
            epi<<(prbiflst[j].probechr==0?"NA":atos(prbiflst[j].probechr))<<'\t'<<prbiflst[j].probeId<<'\t'<<0<<'\t'<<(prbiflst[j].bp==0?"NA":atos(prbiflst[j].bp))<<'\t'<<prbiflst[j].genename<<'\t'<<prbiflst[j].orien<<'\n';
        }
        epi.close();
        printf("%ld probes have been saved in the file %s.\n",epiNum,epifile.c_str());
        
        
        printf("\nGenerating esi file...\n");
        snpinfolst* esiptr=&snpinfo[0];
        qsort(esiptr,esiNum,sizeof(snpinfolst),comp_esi);
        epi_map.clear();
        esi_map.clear();
        vector<string> esi_rs(esiNum);
        vector<string> esi_a1(esiNum);
        vector<string> esi_a2(esiNum);
        string esifile =  string(outFileName)+string(".esi");
        ofstream esi(esifile.c_str());
        if (!esi) throw ("Error: can not open the esi file to save!");
        for (int j = 0;j <esiNum; j++) {
            esi<<(snpinfo[j].snpchr==0?"NA":atos(snpinfo[j].snpchr))<<'\t'<<snpinfo[j].snprs<<'\t'<<snpinfo[j].gd<<'\t'<<(snpinfo[j].bp==0?"NA":atos(snpinfo[j].bp))<<'\t'<<snpinfo[j].a1<<'\t'<<snpinfo[j].a2<<'\n';
            esi_rs[j]=snpinfo[j].snprs;
            esi_map.insert(pair<string,int>(snpinfo[j].snprs,j));
            esi_a1[j]=snpinfo[j].a1;
            esi_a2[j]=snpinfo[j].a2;
        }
        esi.close();
        printf("%ld SNPs have been saved in the file %s.\n",esiNum,esifile.c_str());
        
        double sparsity=1.0*lineNum/(esiNum*epiNum);
        rsbp_map.clear();
        rsaa_map.clear();
        probe_cp_map.clear();
        probe_cpb_map.clear();
        printf("\nGenerating besd file...\n");
       
        
        if(save_dense_flag)
        {
            if(sparsity>0.4)
            {
                printf("The sparsity of your data matrix is %f. We are going to save it in dense format!\n", sparsity);
                
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
                map<string, int>::iterator iter;
                for(int j=0;j<epiNum;j++)
                {
                    printf("Redinging... %3.0f%%\r", 100.0*j/epiNum);
                    fflush(stdout);
                    for(int k=0;k<bsize;k++) buffer[k]=-9; //init
                    vector<string> _rs;
                    vector<float> _beta;
                    vector<float> _se;
                    vector<string> _a1;
                    vector<string> _a2;
                    string curprb=prbiflst[j].probeId;
                    for( int k=0;k<_ttl_probe.size();k++)
                        if(curprb==_ttl_probe[k])
                        {
                            _rs.push_back(_ttl_rs[k]);
                            _a1.push_back(_ttl_a1[k]);
                            _a2.push_back(_ttl_a2[k]);
                            _beta.push_back(_ttl_beta[k]);
                            _se.push_back(_ttl_se[k]);
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
                        string chckstr=_rs[l]+":"+_a1[l]+":"+_a2[l];
                        rsa_map.insert(pair<string,int>(chckstr,l));
                        if(rsNum<rsa_map.size())
                        {
                            if(esi_a1[rsid[l]]==_a1[l] && esi_a2[rsid[l]]==_a2[l]){
                                buffer[rsid[l]]=_beta[l];
                                buffer[esiNum+rsid[l]]=_se[l];
                            } else if(esi_a1[rsid[l]]==_a2[l] && esi_a2[rsid[l]]==_a1[l]){
                                // can't enter here, cos it has done in reading.
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
                            rsNum=rsa_map.size();
                        } else {
                            printf("WARNING: duplicate SNP %s with the same alleles belonging to the same probe %s found. \n",_rs[l].c_str(),curprb.c_str());
                        }
                    }
                    fwrite (buffer,sizeof(float), bsize, smr1);
                }
                fclose (smr1);
                free(buffer);
                cout<<"Beta values and SE values for "<<epiNum<<" Probes and "<<esiNum<<" SNPs have been saved in the dense binary file [" + esdfile + "]." <<endl;
                
            } else {
                
                printf("The sparsity of your data matrix is %f. We are going to save it in sparse format!\n", sparsity);
                
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
                
                map<string, int>::iterator iter;
                for(int j=0;j<epiNum;j++)
                {
                    printf("Saving... %3.0f%%\r", 100.0*j/epiNum);
                    fflush(stdout);
                    vector<string> _rs;
                    vector<float> _beta;
                    vector<float> _se;
                    vector<string> _a1;
                    vector<string> _a2;
                    string curprb=prbiflst[j].probeId;
                    for( int k=0;k<_ttl_probe.size();k++)
                        if(curprb==_ttl_probe[k])
                        {
                            _rs.push_back(_ttl_rs[k]);
                            _a1.push_back(_ttl_a1[k]);
                            _a2.push_back(_ttl_a2[k]);
                            _beta.push_back(_ttl_beta[k]);
                            _se.push_back(_ttl_se[k]);
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
                  
                    if(rsid.size()!=_rs.size())
                    {
                        printf("ERROR: SNP names don't match in %s. Please report this bug.", curprb.c_str());
                        exit(EXIT_FAILURE);
                    }

                    map<string, int> rsa_map;
                    long rsNum=0;
                    vector<uint32_t> tmprid;
                    vector<float> tmpse;
                    for(int l=0;l<rsid.size();l++)
                    {
                        if(abs(_se[l]+9)>1e-6)
                        {
                            string chckstr=_rs[l]+":"+_a1[l]+":"+_a2[l];
                            rsa_map.insert(pair<string,int>(chckstr,l));
                            if(rsNum<rsa_map.size())
                            {
                                if(esi_a1[rsid[l]]==_a1[l] && esi_a2[rsid[l]]==_a2[l] )
                                {
                                    val.push_back(_beta[l]);
                                    rowids.push_back(rsid[l]);
                                    tmpse.push_back(_se[l]);
                                    tmprid.push_back(rsid[l]);
                                } else if(esi_a1[rsid[l]]==_a2[l] && esi_a2[rsid[l]]==_a1[l] )
                                {
                                     // can't enter here, cos it has done in reading.
                                    val.push_back(-1.0*_beta[l]);
                                    rowids.push_back(rsid[l]);
                                    tmpse.push_back(_se[l]);
                                    tmprid.push_back(rsid[l]);
                                } else
                                {
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
                            } else
                            {
                                printf("WARNING: duplicate SNP %s with the same alleles belonging to the same probe %s found. \n",_rs[l].c_str(),curprb.c_str());
                            }
                        } else
                        {
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
                }
                uint64_t valNum=val.size();
                fwrite (&valNum,sizeof(uint64_t), 1, smr1);
                fwrite (&cols[0],sizeof(uint64_t), cols.size(), smr1);
                fwrite (&rowids[0],sizeof(uint32_t), rowids.size(), smr1);
                fwrite (&val[0],sizeof(float), val.size(), smr1);
                fclose (smr1);
                
                cout<<"Beta values and SE values for "<<epiNum<<" Probes and "<<esiNum<<" SNPs have been saved in the sparse binary file [" + esdfile + "]." <<endl;
            }
        } else {
            // get esd info
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
                printf("Redinging... %3.0f%%\r", 100.0*j/epiNum);
                fflush(stdout);
                vector<uint32_t> tmprid;
                vector<float> tmpse;
                string curprb=prbiflst[j].probeId;
                snpinfo.clear();
                
                for( int k=0;k<_ttl_probe.size();k++)
                    if(curprb==_ttl_probe[k])
                    {
                        snpinfolst tmpinfo;
                        tmpinfo.snprs=_ttl_rs[k];
                        tmpinfo.a1=_ttl_a1[k];
                        tmpinfo.a2=_ttl_a2[k];
                        tmpinfo.beta=_ttl_beta[k];
                        tmpinfo.se=_ttl_se[k];
                        tmpinfo.snpchr=_ttl_chr[k];
                        tmpinfo.bp=_ttl_bp[k];
                        snpinfo.push_back(tmpinfo);
                    }
                
                snpinfolst* sortptr=&snpinfo[0];
                qsort(sortptr,snpinfo.size(),sizeof(snpinfolst),comp_esi);
                
                probeinfolst prbifo=prbiflst[j];
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
                if(rsid.size()!=_rs.size())
                {
                    printf("ERROR: SNP names don't match in probe %s. Please report this bug.", prbiflst[j].probeId.c_str());
                    exit(EXIT_FAILURE);
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
                            printf("WARNING: duplicate SNP %s with the same alleles belonging to the same probe %s. \n",_rs[l].c_str(),prbiflst[j].probeId.c_str());
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
        
    }
}