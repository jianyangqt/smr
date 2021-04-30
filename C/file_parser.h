#ifndef FILE_FARSER_CALMT
#define FILE_FARSER_CALMT

struct MSMR_DATA_eqtl {
    struct MSMR_NODE_eqtl {
        char expo_id[20];
        unsigned char expo_chr;
        char expo_gene[20];
        unsigned int expo_bp;
        char outco_id[20];
        unsigned char outco_chr;
        char outco_gene[20];
        unsigned int outco_bp;
        char top_snp[20];
        unsigned char top_snp_chr;
        unsigned int top_snp_bp;
        char a1[10];
        char a2[10];
        unsigned char frequence_na;
        double frequence;
        unsigned char b_outco_na;
        double b_outco;
        unsigned char se_outco_na;
        double se_outco;
        unsigned char p_outco_na;
        double p_outco;
        unsigned char b_expo_na;
        double b_expo;
        unsigned char se_expo_na;
        double se_expo;
        unsigned char p_expo_na;
        double p_expo;
        unsigned char b_smr_na;
        double b_smr;
        unsigned char se_smr_na;
        double se_smr;
        unsigned char p_smr_na;
        double p_smr;
        unsigned char p_smr_multi_na;
        double p_smr_multi;
        unsigned char p_heidi_na;
        double p_heidi;
        unsigned char nsnp_heidi_na;
        unsigned short  nsnp_heidi;
        
        struct MSMR_NODE_eqtl * prev;
        struct MSMR_NODE_eqtl * next;

    } * msmr_node;

    unsigned int node_len;
};


extern struct MSMR_DATA_eqtl parse_msmr_file_eqtl(const char *, bool);


struct MSMR_DATA_gwas {
    struct MSMR_NODE_gwas {
        char probe_id[50];
        unsigned char probe_chr;
        char gene[50];
        unsigned int probe_bp;
        char top_snp[50];
        unsigned char top_snp_chr;
        unsigned int top_snp_bp;
        char a1[10];
        char a2[10];

        unsigned char frequence_na;
        double frequence;
        unsigned char b_gwas_na;
        double b_gwas;
        unsigned char se_gwas_na;
        double se_gwas;
        unsigned char p_gwas_na;
        double p_gwas;
        unsigned char b_eqtl_na;
        double b_eqtl;
        unsigned char se_eqtl_na;
        double se_eqtl;
        unsigned char p_eqtl_na;
        double p_eqtl;
        unsigned char b_smr_na;
        double b_smr;
        unsigned char se_smr_na;
        double se_smr;
        unsigned char p_smr_na;
        double p_smr;
        unsigned char p_smr_multi_na;
        double p_smr_multi;
        unsigned char p_heidi_na;
        double p_heidi;
        unsigned char nsnp_heidi_na;
        unsigned short nsnp_heidi;
        
        struct MSMR_NODE_gwas * prev;
        struct MSMR_NODE_gwas * next;

    } * msmr_node;

    
    unsigned int node_len;

};


extern struct MSMR_DATA_gwas parse_msmr_file_gwas(const char *, bool);


struct ANNO_DATA {
    struct ANNO_NODE {
        unsigned char chr;
        char gene_id[50];
        char gene_name[50];

        unsigned char gene_type_na;
        char gene_type[20];
        unsigned char gene_start_na;
        unsigned int gene_start;
        unsigned char gene_end_na;
        unsigned int gene_end;
        unsigned char orientation_na;
        char orientation;

        struct ANNO_NODE * prev;
        struct ANNO_NODE * next;
    } * data_node;

    unsigned int node_len;
};


extern struct ANNO_DATA parse_anno_file_msmr(const char *, bool);



#endif
