//
//  SMR.h
//  SMR_CPP
//
//  Created by Futao Zhang on 1/09/2015.
//  Copyright (c) 2015 Futao Zhang. All rights reserved.
//

#ifndef SMR_CPP_SMR_h
#define SMR_CPP_SMR_h
void option(int option_num, char* option_str[]);
static inline bool not_in_flags(vector<string> &flags, string str)
{
    return find(flags.begin(),flags.end(),str) == flags.end();
}
static inline void FLAGS_VALID_CK(int option_num, char* option_str[])
{
    const char *flgs[] = { "--bfile","--gwas-summary","--beqtl-summary","--maf","--keep","--remove","--extract-snp","--exclude-snp","--extract-probe","--extract-outcome-probe","--extract-exposure-probe",
        "--exclude-probe","--exclude-outcome-probe","--exclude-exposure-probe","--eqtl-summary","--ld-pruning","--peqtl-heidi","--heidi-m","--make-besd","--out", "--peqtl-smr","--smr",
        "--cis-wind","--peqtl-trans","--peqtl-other","--efile","--query","--heidi-off","--target-snp","--extract-trait","--thread-num","--combine-besd","--beqtl-flist",
        "--trans-wind","--plot","--trans","--eqtl-flist","--smr-format","--merlin-fastassoc-format","--plink-qassoc-format","--gemma-format","--freq","--esd-std",
        "--meta","--genes","--smr-wind","--est-n","--internal-test", "--smr-file","--recode","--probe-var","--make-sparse","--chr","--probe-chr","--snp-chr", "--snp", "--from-snp", "--to-snp", "--probe", "--from-probe", "--to-probe","--snp-wind","--probe-wind", "--gene","--from-snp-kb","--to-snp-kb","--from-probe-kb","--to-probe-kb","--extract-single-exposure-probe","--extract-single-outcome-probe","--exclude-single-exposure-probe","--exclude-single-outcome-probe","--gene-list", "--set-list", "--set-wind","--smr-multi","--qfile","--save-dense"
    };
    
    vector<string> flags(flgs, flgs + sizeof(flgs)/sizeof(flgs[0]));
    
    if(option_num<3)
    {
        cout<<"flags include:"<<endl;
        int cur_mark=0;
        for(int i=0;i<flags.size();i++)
        {
            int tmp=i>>2;
            if(tmp>cur_mark)
            {
                cout<<endl;
                cur_mark=tmp;
            }
            cout<<flags[i]<<",";
        }
        cout<<endl;
        exit (EXIT_FAILURE);
    }
    for(int i=0;i<option_num;i++)
    {
        if(SMRDATA::has_prefix(option_str[i],"--"))
            if(not_in_flags(flags, option_str[i]))
            {
                
                fprintf (stderr, "%s: Invalid option\n",
                         option_str[i]);
                exit (EXIT_FAILURE);
            }
    }
    
}
static inline void FLAG_VALID_CK(string str, char* flag)
{
    if(flag==NULL || SMRDATA::has_prefix(flag, "--"))
    {
        fprintf (stderr, "Please verify the flag %s!: \n",
                 str.c_str());
        exit (EXIT_FAILURE);
    }
}


#endif
