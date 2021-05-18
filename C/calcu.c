#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


struct Anno_data{
    char gene_id[256];
    char gene_name[256];
    struct Anno_data * left;
    struct Anno_data * right;
};


//data structure used to contain eqtl and gwas data.
struct eqtl_gwas_data{
    //to eqtl feature.
    char Expo_ID[256];
    unsigned char Outco_Chr;
    char Outco_Gene[256];
    unsigned int Outco_bp;
    float p_SMR_multi;
    float p_HEIDI;
    char gene_name[256];
    // to gwas freature.
    char probeID[256];
    float b_SMR;
    float p_SMR;
    float p_SMR_multi_gwas;
    float p_ACAT;
    float b_SMR_weight;
    struct eqtl_gwas_data * next;
};


//data structure used to contain eqtl and pqtl data.
//eqtl need search anno to get gene name but pqtl needed.
struct gwas_data{
    char probeID[256];
    char Gene[256];
    float b_SMR;
    float p_SMR;
    float p_SMR_multi;
    float p_ACAT;
    float b_SMR_weight;
    char gene_name[256];
    struct gwas_data * next;
};


struct linker{
    struct linker * next;
    void * ptr;
};


struct gene_mapper{
    char gene_name[256];
    struct linker * caqtl;
    struct linker * mqtl;
    struct linker * h27ac;
    struct linker * k4me1;
    struct linker * eqtl;
    struct linker * pqtl;
    struct gene_mapper * next;
};


struct value_list{
    float value;
    struct value_list * next;

};


static struct Anno_data * read_anno_file(const char *);
static int add_tree_node(struct Anno_data *, struct Anno_data *);
static void split_gene_id(char *);
static struct eqtl_gwas_data * read_eqtl_file(const char *);
static struct gwas_data * read_gwas_file(const char *);
static struct eqtl_gwas_data * filter_eqtl_data(struct eqtl_gwas_data *, float);
static struct Anno_data * search_anno_tree(struct Anno_data *, const char *);
static void add_gene_name_to_eqtl(struct eqtl_gwas_data *, struct Anno_data *);
static void add_gene_name_to_gwas(struct gwas_data *, struct Anno_data *);
static void calculate_p_ACAT(struct gwas_data *);
static void union_eqtl_and_gwas(struct eqtl_gwas_data *, struct gwas_data **);
static struct Anno_data * find_all_gene_name(struct eqtl_gwas_data *, struct eqtl_gwas_data *, \
    struct eqtl_gwas_data *, struct eqtl_gwas_data *, struct gwas_data *, struct gwas_data *);
static int add_tree_node_2(struct Anno_data *, struct Anno_data *);
static struct gene_mapper * make_gene_list(struct Anno_data *);
static struct gene_mapper * add_gene(struct gene_mapper *, struct Anno_data *);
static void collect_gene_blongs(struct gene_mapper *, void *, int);
static struct eqtl_gwas_data * wash_eqtl(struct eqtl_gwas_data *);
static struct gwas_data * wash_gwas(struct gwas_data *);
static void print_res(struct gene_mapper *, const char *);
static float trav_linker(struct linker *, int);
static void trav_linker_2(struct linker *, int, float);
static void value_link(struct linker *, int *, unsigned int *, float *, float *, \
    struct value_list **, struct value_list **, struct value_list **);
static void value_link_2(struct linker *, float *, float *, struct value_list **, \
    struct value_list **, struct value_list **);
static float acat_func(struct value_list *);
static float acat_weight_func(struct value_list *, struct value_list *);




static void
trav_tree(struct Anno_data * tree)
{
    if (tree){
        printf("%s\n", tree -> gene_name);
        trav_tree(tree -> left);
        trav_tree(tree -> right);
    }
    return;
}

int
main(int argc, char * argv[])
{
    char * file_list;
    if (argc != 2){
        fprintf(stderr, "error, argument error.\n");
        exit(1);
    }
    file_list = argv[1];

    FILE * f_in = fopen(file_list, "r");
    if (!f_in){
        fprintf(stderr, "error, open file list error.\n");
        exit(1);
    }

    char line_buffer[1024];
    char anno_f[512];
    char caqtl_eqtl_f[512];
    char mqtl_eqtl_f[512];
    char h27ac_eqtl_f[512];
    char k4me1_eqtl_f[512];
    char caqtl_gwas_f[512];
    char mqtl_gwas_f[512];
    char h27ac_gwas_f[512];
    char k4me1_gwas_f[512];
    char eqtl_gwas_f[512];
    char pqtl_gwas_f[512];
    char key[512];
    char value[512];
    char chr = '\0';
    int i = 0;
    int j = 0;
    char sw = 1;
    while(fgets(line_buffer, 1023, f_in)){
        //printf("%s", line_buffer);
        j = 0;
        sw = 1;
        for (i = 0; i < strlen(line_buffer); i++){
            chr = line_buffer[i];
            if (chr != ' ' && chr != '\t' && chr != '\n'){
                if (sw){
                    if (chr == '='){
                        key[j] = '\0';
                        sw = 0;
                        j = 0;
                    } else{
                        key[j] = chr;
                        j++;
                    }
                } else{
                    value[j] = chr;
                    j++;
                }
            }
        }
        value[j] = '\0';
        //printf(">%s\n", key);
        //printf(">%s\n", value);
        if (!strcmp(key, "anno_file")){
            strcpy(anno_f, value);
        }else if (!strcmp(key, "caqtl_eqtl")){
            strcpy(caqtl_eqtl_f, value);
        } else if (!strcmp(key, "mqtl_eqtl")){
            strcpy(mqtl_eqtl_f, value);
        } else if (!strcmp(key, "h27ac_eqtl")){
            strcpy(h27ac_eqtl_f, value);
        } else if (!strcmp(key, "k4me1_eqtl")){
            strcpy(k4me1_eqtl_f, value);
        } else if (!strcmp(key, "caqtl_gwas")){
            strcpy(caqtl_gwas_f, value);
        } else if (!strcmp(key, "mqtl_gwas")){
            strcpy(mqtl_gwas_f, value);
        } else if (!strcmp(key, "h27ac_gwas")){
            strcpy(h27ac_gwas_f, value);
        } else if(!strcmp(key, "k4me1_gwas")){
            strcpy(k4me1_gwas_f, value);
        } else if (!strcmp(key, "eqtl_gwas")){
            strcpy(eqtl_gwas_f, value);
        } else if (!strcmp(key, "pqtl_gwas")){
            strcpy(pqtl_gwas_f, value);
        } else{
            fprintf(stderr, "not recognized file list key.\n");
            exit(1);
        }

    }
    fclose(f_in);

    struct Anno_data * anno_dt;
    anno_dt = read_anno_file(anno_f);
    //trav_tree(anno_dt);

    struct eqtl_gwas_data * caqtl_eqtl = NULL;
    struct gwas_data * caqtl_gwas = NULL;
    caqtl_eqtl = read_eqtl_file(caqtl_eqtl_f);
    caqtl_gwas = read_gwas_file(caqtl_gwas_f);
    caqtl_eqtl = filter_eqtl_data(caqtl_eqtl, 0.01);
    add_gene_name_to_eqtl(caqtl_eqtl, anno_dt);
    calculate_p_ACAT(caqtl_gwas);
    union_eqtl_and_gwas(caqtl_eqtl, &caqtl_gwas);


    struct eqtl_gwas_data * mqtl_eqtl = NULL;
    struct gwas_data * mqtl_gwas = NULL;
    mqtl_eqtl = read_eqtl_file(mqtl_eqtl_f);
    mqtl_gwas = read_gwas_file(mqtl_gwas_f);
    mqtl_eqtl = filter_eqtl_data(mqtl_eqtl, 0.01);
    add_gene_name_to_eqtl(mqtl_eqtl, anno_dt);
    calculate_p_ACAT(mqtl_gwas);
    union_eqtl_and_gwas(mqtl_eqtl, &mqtl_gwas);


    struct eqtl_gwas_data * h27ac_eqtl = NULL;
    struct gwas_data * h27ac_gwas = NULL;
    h27ac_eqtl = read_eqtl_file(h27ac_eqtl_f);
    h27ac_gwas = read_gwas_file(h27ac_gwas_f);
    h27ac_eqtl = filter_eqtl_data(h27ac_eqtl, 0.01);
    add_gene_name_to_eqtl(h27ac_eqtl, anno_dt);
    calculate_p_ACAT(h27ac_gwas);
    union_eqtl_and_gwas(h27ac_eqtl, &h27ac_gwas);


    struct eqtl_gwas_data * k4me1_eqtl;
    struct gwas_data * k4me1_gwas;
    k4me1_eqtl = read_eqtl_file(k4me1_eqtl_f);
    k4me1_gwas = read_gwas_file(k4me1_gwas_f);
    k4me1_eqtl = filter_eqtl_data(k4me1_eqtl, 0.01);
    add_gene_name_to_eqtl(k4me1_eqtl, anno_dt);
    calculate_p_ACAT(k4me1_gwas);
    union_eqtl_and_gwas(k4me1_eqtl, &k4me1_gwas);



    struct gwas_data * eqtl_gwas = NULL;
    eqtl_gwas = read_gwas_file(eqtl_gwas_f);
    add_gene_name_to_gwas(eqtl_gwas, anno_dt);
    calculate_p_ACAT(eqtl_gwas);




    struct gwas_data * pqtl_gwas = NULL, * pqtl_gwas_tmp;
    pqtl_gwas = read_gwas_file(pqtl_gwas_f);
    calculate_p_ACAT(pqtl_gwas);
    pqtl_gwas_tmp = pqtl_gwas;
    while(pqtl_gwas_tmp){
        strcpy(pqtl_gwas_tmp -> gene_name, pqtl_gwas_tmp -> Gene);
        pqtl_gwas_tmp = pqtl_gwas_tmp -> next;
    }



    caqtl_eqtl = wash_eqtl(caqtl_eqtl);
    mqtl_eqtl = wash_eqtl(mqtl_eqtl);
    h27ac_eqtl = wash_eqtl(h27ac_eqtl);
    k4me1_eqtl = wash_eqtl(k4me1_eqtl);
    eqtl_gwas = wash_gwas(eqtl_gwas);
    pqtl_gwas = wash_gwas(pqtl_gwas);

/*
    struct eqtl_gwas_data * ttp;
    ttp = k4me1_eqtl;
    while(ttp){
        printf("%s %s %s %e %e\n", ttp -> Expo_ID, ttp -> gene_name, ttp -> probeID, \
            ttp -> p_SMR_multi, ttp -> p_ACAT);
        ttp = ttp -> next;
    }
    exit(0);
*/

/*
    struct gwas_data * ppt;
    ppt = pqtl_gwas;
    while(ppt){
        printf("%s %s %s %e\n", ppt -> probeID, ppt -> Gene, ppt -> gene_name, ppt -> p_ACAT);
        ppt = ppt -> next;
    }
*/


    struct Anno_data * all_gene;
    all_gene = find_all_gene_name(caqtl_eqtl, mqtl_eqtl, h27ac_eqtl, k4me1_eqtl, \
        eqtl_gwas, pqtl_gwas);

    //trav_tree(all_gene);
    struct gene_mapper * gene_mapper_chain = NULL;
    gene_mapper_chain = make_gene_list(all_gene);

/*
    while(gene_mapper_chain){
        printf("%s\n", gene_mapper_chain -> gene_name);
        gene_mapper_chain = gene_mapper_chain -> next;
    }
*/



    collect_gene_blongs(gene_mapper_chain, caqtl_eqtl, 1);
    collect_gene_blongs(gene_mapper_chain, mqtl_eqtl, 2);
    collect_gene_blongs(gene_mapper_chain, h27ac_eqtl, 3);
    collect_gene_blongs(gene_mapper_chain, k4me1_eqtl, 4);
    collect_gene_blongs(gene_mapper_chain, eqtl_gwas, 5);
    collect_gene_blongs(gene_mapper_chain, pqtl_gwas, 6);

/*
    struct linker * linker_tmp;
    struct eqtl_gwas_data * eqtl_tmp;
    struct gwas_data * gwas_tmp;
    while(gene_mapper_chain){
        printf(">%s\n", gene_mapper_chain -> gene_name);
        printf("--caqtl dt\n");
        linker_tmp = gene_mapper_chain -> caqtl;
        while (linker_tmp){
            eqtl_tmp = (struct eqtl_gwas_data *)linker_tmp -> ptr;
            printf("%s\n",  eqtl_tmp -> gene_name);
            linker_tmp = linker_tmp -> next;
        }

        printf("--mqtl dt\n");
        linker_tmp = gene_mapper_chain -> mqtl;
        while (linker_tmp){
            eqtl_tmp = (struct eqtl_gwas_data *)linker_tmp -> ptr;
            printf("%s\n",  eqtl_tmp -> gene_name);
            linker_tmp = linker_tmp -> next;
        }

        printf("--h27ac dt\n");
        linker_tmp = gene_mapper_chain -> h27ac;
        while (linker_tmp){
            eqtl_tmp = (struct eqtl_gwas_data *)linker_tmp -> ptr;
            printf("%s\n",  eqtl_tmp -> gene_name);
            linker_tmp = linker_tmp -> next;
        }

        printf("--k4me1 dt\n");
        linker_tmp = gene_mapper_chain -> k4me1;
        while (linker_tmp){
            eqtl_tmp = (struct eqtl_gwas_data *)linker_tmp -> ptr;
            printf("%s\n",  eqtl_tmp -> gene_name);
            linker_tmp = linker_tmp -> next;
        }

        printf("--eqtl dt\n");
        linker_tmp = gene_mapper_chain -> eqtl;
        while (linker_tmp){
            gwas_tmp = (struct gwas_data *)linker_tmp -> ptr;
            printf("%s\n",  gwas_tmp -> gene_name);
            linker_tmp = linker_tmp -> next;
        }

        printf("--pqtl dt\n");
        linker_tmp = gene_mapper_chain -> pqtl;
        while (linker_tmp){
            gwas_tmp = (struct gwas_data *)linker_tmp -> ptr;
            printf("%s\n",  gwas_tmp -> gene_name);
            linker_tmp = linker_tmp -> next;
        }

        gene_mapper_chain = gene_mapper_chain -> next;
    }
*/


    print_res(gene_mapper_chain, "res.csv");
}


static struct Anno_data *
read_anno_file(const char * file_name)
{
    char chromo[5];
    char gene_id[256];
    char gene_name[256];
    char gene_type[256];
    char start[32];
    char end[32];
    char strand[32];
    FILE * f_in = fopen(file_name, "r");
    int chr;
    while((chr = fgetc(f_in)) != '\n'){
        ;
    }
    struct Anno_data * first_node = NULL, * new_node = NULL;
    while(fscanf(f_in, "%s %s %s %s %s %s %s", chromo, gene_id, gene_name, gene_type, start, end, strand) == 7){
        new_node = (struct Anno_data *)malloc(sizeof(struct Anno_data));
        new_node -> left = NULL;
        new_node -> right = NULL;
        split_gene_id(gene_id);
        strcpy(new_node -> gene_id, gene_id);
        strcpy(new_node -> gene_name, gene_name);
        if (first_node){
            if (add_tree_node(first_node, new_node)){
                free(new_node);
            }
        } else{
            first_node = new_node;
        }
    }

    return first_node;
}


static int
add_tree_node(struct Anno_data * des_node, struct Anno_data * new_node)
{
    int cmp_res = 0;
    cmp_res = strcmp(des_node -> gene_id, new_node -> gene_id);

    if (cmp_res < 0){
        if (des_node -> left){
            add_tree_node(des_node -> left, new_node);
        } else{
            des_node -> left = new_node;
        }
    } else if (cmp_res > 0){
        if (des_node -> right){
            add_tree_node(des_node -> right, new_node);
        } else{
            des_node -> right = new_node;
        }
    } else{
        fprintf(stderr, "waring, duplicated gene_id found\n");
        return 1;
    }
    return 0;
}


static void
split_gene_id(char * gene_id)
{
    char tmp[256];
    int i = 0;
    for (i = 0; i < strlen(gene_id); i++){
        if (gene_id[i] != '.'){
            tmp[i] = gene_id[i];
        } else{
            break;
        }
    }
    tmp[i] = '\0';
    strcpy(gene_id, tmp);
    return;
}


static struct eqtl_gwas_data *
read_eqtl_file(const char * file_name)
{
    FILE * f_in = fopen(file_name, "r");
    if (!f_in){
        fprintf(stderr, "open %s failed\n", file_name);
        exit(1);
    }
    int chr;
    //remove head line
    while((chr = fgetc(f_in)) != '\n'){
        ;
    }
    char Expo_ID[256];
    char Expo_chr[5];
    char Expo_Gene[256];
    char Expo_bp[32];
    char Outco_ID[256];
    char Outco_Chr[5];
    char Outco_Gene[256];
    char Outco_bp[32];
    char top_SNP[256];
    char top_SNP_chr[256];
    char top_SNP_bp[256];
    char A1[512];
    char A2[512];
    char Freq[32];
    char b_Outco[32];
    char se_Outco[32];
    char p_Outco[32];
    char b_Expo[32];
    char se_Expo[32];
    char p_Expo[32];
    char b_SMR[32];
    char se_SMR[32];
    char p_SMR[32];
    char p_SMR_multi[32];
    char p_HEIDI[32];
    char nsnp_HEIDI[32];
    struct eqtl_gwas_data * first_node = NULL, * new_node = NULL, * ptr = NULL;

    while(fscanf(f_in, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \
        %s %s %s %s %s %s %s %s %s", Expo_ID, Expo_chr, Expo_Gene, Expo_bp, Outco_ID, \
        Outco_Chr, Outco_Gene, Outco_bp, top_SNP, top_SNP_chr, top_SNP_bp, \
        A1, A2, Freq, b_Outco, se_Outco, p_Outco, b_Expo, se_Expo, p_Expo, \
        b_SMR, se_SMR, p_SMR, p_SMR_multi, p_HEIDI, nsnp_HEIDI) == 26){

        new_node = (struct eqtl_gwas_data *)malloc(sizeof(struct eqtl_gwas_data));
        new_node -> next = NULL;
        strcpy(new_node -> gene_name, "NA");
        strcpy(new_node -> probeID, "NA");
        new_node -> b_SMR = -1;
        new_node -> p_SMR = -1;
        new_node -> p_SMR_multi_gwas = -1;
        new_node -> p_ACAT = -1;
        new_node -> b_SMR_weight = -1;
        strcpy(new_node -> Expo_ID, Expo_ID);
        if (!strcmp(Outco_Chr, "NA")){
            if (strlen(Outco_Chr) == 1){
                if (Outco_Chr[0] == 'X' || Outco_Chr[0] == 'x')
                    new_node -> Outco_Chr = 23;
                else if (Outco_Chr[0] == 'Y' || Outco_Chr[0] == 'y')
                    new_node -> Outco_Chr = 24;
                else if (isdigit(Outco_Chr[0]))
                    new_node -> Outco_Chr = atoi(Outco_Chr);
                else{
                    fprintf(stderr, "Outco_Chr not recognized\n");
                    exit(1);
                }
            } else if (isdigit(Outco_Chr[0])){
                new_node -> Outco_Chr = atoi(Outco_Chr);
            } else{
                fprintf(stderr, "Outco_Chr not recognized\n");
                exit(1); //actually, befor exit, should close file and free mem.
            }
        } else{
            new_node -> Outco_Chr = 0;
        }
        strcpy(new_node -> Outco_Gene, Outco_Gene);
        if (!strcmp(Outco_bp, "NA")){
            new_node -> Outco_bp = 0;
        } else{
            new_node -> Outco_bp = (unsigned int)atol(Outco_bp);
        }

        if (!strcmp(p_SMR_multi, "NA")){
            new_node -> p_SMR_multi = -1;
        } else{
            new_node -> p_SMR_multi = atof(p_SMR_multi);
        }

        if (!strcmp(p_HEIDI, "NA")){
            new_node -> p_HEIDI = -1;
        } else{
            new_node -> p_HEIDI = atof(p_HEIDI);
        }

        if (first_node){
            ptr -> next = new_node;
            ptr = new_node;
        } else{
            first_node = ptr = new_node;
        }
    }
        fclose(f_in);
    return first_node;
}


static struct gwas_data *
read_gwas_file(const char * file_name)
{
    FILE * f_in = fopen(file_name, "r");
    if (!f_in){
        fprintf(stderr, "error, %s open failed.\n", file_name);
        exit(1);
    }
    int chr;
    while((chr = fgetc(f_in)) != '\n'){
        ;
    }
    char probeID[256];
    char ProbeChr[5];
    char Gene[256];
    char Probe_bp[32];
    char topSNP[256];
    char topSNP_chr[5];
    char topSNP_bp[32];
    char A1[512];
    char A2[512];
    char Freq[32];
    char b_GWAS[32];
    char se_GWAS[32];
    char p_GWAS[32];
    char b_eQTL[32];
    char se_eQTL[32];
    char p_eQTL[32];
    char b_SMR[32];
    char se_SMR[32];
    char p_SMR[32];
    char p_SMR_multi[32];
    char p_HEIDI[32];
    char nsnp_HEIDI[32];
    struct gwas_data * first_node = NULL, * new_node = NULL, * ptr = NULL;
    while(fscanf(f_in, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \
        %s %s %s %s", probeID, ProbeChr, Gene, Probe_bp, topSNP, topSNP_chr, \
        topSNP_bp, A1, A2, Freq, b_GWAS, se_GWAS, p_GWAS, b_eQTL, se_eQTL, p_eQTL, \
        b_SMR, se_SMR, p_SMR, p_SMR_multi, p_HEIDI, nsnp_HEIDI) == 22){

        new_node = (struct gwas_data *)malloc(sizeof(struct gwas_data));
        new_node -> next = NULL;
        new_node -> p_ACAT = -1;
        new_node -> b_SMR_weight = -1;
        strcpy(new_node -> gene_name, "NA");
        strcpy(new_node -> probeID, probeID);
        strcpy(new_node -> Gene, Gene);
        if (!strcmp(b_SMR, "NA")){
            new_node -> b_SMR = -1;
        } else{
            new_node -> b_SMR = atof(b_SMR);
        }

        if (!strcmp(p_SMR, "NA")){
            new_node -> p_SMR = -1;
        } else{
            new_node -> p_SMR = atof(p_SMR);
        }

        if (!strcmp(p_SMR_multi, "NA")){
            new_node -> p_SMR_multi = -1;
        } else{
            new_node -> p_SMR_multi = atof(p_SMR_multi);
        }

        if (first_node){
            ptr -> next = new_node;
            ptr = new_node;
        } else{
            first_node = ptr = new_node;
        }

    }

    return first_node;
}


static struct eqtl_gwas_data *
filter_eqtl_data(struct eqtl_gwas_data * eqtl_dt, float p_HEIDI_thr)
{
    struct eqtl_gwas_data * dt_out, * ptr, * tmp;
    unsigned int counter = 0;
    dt_out = eqtl_dt;
    while(dt_out){
        counter++;
        dt_out = dt_out -> next;
    }
    float p_SMR_multi_thr = 0.5 / counter;
    dt_out = ptr = eqtl_dt;
    eqtl_dt = eqtl_dt -> next;
    dt_out -> next = NULL;

    while(eqtl_dt){
        if (eqtl_dt -> p_SMR_multi != -1 && eqtl_dt -> p_SMR_multi < p_SMR_multi_thr \
                && eqtl_dt -> p_HEIDI != -1 && eqtl_dt -> p_HEIDI > p_HEIDI_thr){
            ptr -> next = eqtl_dt;
            ptr = eqtl_dt;
            tmp = eqtl_dt -> next;
            eqtl_dt -> next = NULL;
        } else{
            tmp = eqtl_dt -> next;
            free(eqtl_dt);
        }
        eqtl_dt = tmp;
    }

    if (dt_out -> p_SMR_multi != -1 && dt_out -> p_SMR_multi < p_SMR_multi_thr \
        && dt_out -> p_HEIDI != -1 && dt_out -> p_HEIDI > p_HEIDI_thr){
        ;
    } else{
        dt_out = dt_out -> next;
    }
    return dt_out;
}


static struct Anno_data *
search_anno_tree(struct Anno_data * anno_dt, const char * search_key)
{
    struct Anno_data * dt_out = NULL;
    if (!anno_dt){
        return anno_dt;
    }
    int cmp_res = strcmp(anno_dt -> gene_id, search_key);
    if (!cmp_res){
        dt_out = anno_dt;
    } else if(cmp_res < 0){
        dt_out = search_anno_tree(anno_dt -> left, search_key);
    } else{
        dt_out = search_anno_tree(anno_dt -> right, search_key);
    }
    return dt_out;
}


static void
add_gene_name_to_eqtl(struct eqtl_gwas_data * eqtl_dt, struct Anno_data * anno_dt)
{
    struct Anno_data * search_res;
    while(eqtl_dt){
        search_res = search_anno_tree(anno_dt, eqtl_dt -> Outco_Gene);
        if (search_res){
            strcpy(eqtl_dt -> gene_name, search_res -> gene_name);
        } else{
            fprintf(stderr, "Warning, %s not search out its match gene name\n", \
                eqtl_dt -> Outco_Gene);
        }
        eqtl_dt = eqtl_dt -> next;
    }
    return;
}


static void
add_gene_name_to_gwas(struct gwas_data * gwas_dt, struct Anno_data * anno_dt)
{
    struct Anno_data * search_res;
    while(gwas_dt){
        search_res = search_anno_tree(anno_dt, gwas_dt -> probeID);
        if (search_res){
            strcpy(gwas_dt -> gene_name, search_res -> gene_name);
        } else{
            fprintf(stderr, "Warning, %s not search out its match gene name\n", \
                gwas_dt -> probeID);
        }
        gwas_dt = gwas_dt -> next;
    }
    return;
}


static void
calculate_p_ACAT(struct gwas_data * gwas_dt)
{
    float p_SMR = 0;
    float p_SMR_multi = 0;
    while(gwas_dt){
        p_SMR = gwas_dt -> p_SMR;
        p_SMR_multi = gwas_dt -> p_SMR_multi;
        if (p_SMR != -1 && p_SMR_multi != -1){
            gwas_dt -> p_ACAT = 0.5 - atan((tan((0.5 - p_SMR) * M_PI) + tan((0.5 - p_SMR_multi) * M_PI)) / 2) / M_PI;
        }
        gwas_dt = gwas_dt -> next;
    }
}



//merge two kind of data.
static void
union_eqtl_and_gwas(struct eqtl_gwas_data * eqtl_dt, struct gwas_data ** gwas_dt)
{
    struct gwas_data * gwas_dt_pt, * tmp;
    struct eqtl_gwas_data * eqtl_pt;
    gwas_dt_pt = *gwas_dt;
    *gwas_dt = NULL;

    while(gwas_dt_pt){
        //printf("%s\n", gwas_dt_pt -> probeID);
        eqtl_pt = eqtl_dt;
        while(eqtl_pt){
            //printf("%s\n", eqtl_pt -> Expo_ID);
            if (!strcmp(eqtl_pt -> Expo_ID, gwas_dt_pt -> probeID)){
                //printf(">heer\n");
                strcpy(eqtl_pt -> probeID, gwas_dt_pt -> probeID);
                eqtl_pt -> b_SMR = gwas_dt_pt -> b_SMR;
                eqtl_pt -> p_SMR = gwas_dt_pt -> p_SMR;
                eqtl_pt -> p_SMR_multi_gwas = gwas_dt_pt -> p_SMR_multi;
                eqtl_pt -> p_ACAT = gwas_dt_pt -> p_ACAT;
            }
            eqtl_pt = eqtl_pt -> next;
        }
        tmp = gwas_dt_pt -> next;
        free(gwas_dt_pt);
        gwas_dt_pt = tmp;
    }
    return;
}


//here, using Anno data data structure to build tree to unique gene name.
static struct Anno_data *
find_all_gene_name(struct eqtl_gwas_data * caqtl_eqtl, struct eqtl_gwas_data * mqtl_eqtl, \
    struct eqtl_gwas_data * h27ac_eqtl, struct eqtl_gwas_data * k4me1_eqtl, \
    struct gwas_data * eqtl_gwas, struct gwas_data * pqtl_gwas)
{
    struct Anno_data * first_node = NULL, * new_node = NULL;
    while(caqtl_eqtl){
        new_node = (struct Anno_data *)malloc(sizeof(struct Anno_data));
        new_node -> left = NULL;
        new_node -> right = NULL;
        strcpy(new_node -> gene_id, "NA");
        strcpy(new_node -> gene_name, caqtl_eqtl -> gene_name);
        if (first_node){
            if (add_tree_node_2(first_node, new_node)){
                free(new_node);
            }

        } else{
            first_node = new_node;
        }

        caqtl_eqtl = caqtl_eqtl -> next;
    }

    while (mqtl_eqtl){
        new_node = (struct Anno_data *)malloc(sizeof(struct Anno_data));
        new_node -> left = NULL;
        new_node -> right = NULL;
        strcpy(new_node -> gene_id, "NA");
        strcpy(new_node -> gene_name, mqtl_eqtl -> gene_name);
        if (add_tree_node_2(first_node, new_node)){
            free(new_node);
        }
        mqtl_eqtl = mqtl_eqtl -> next;
    }

    while (h27ac_eqtl){
        new_node = (struct Anno_data *)malloc(sizeof(struct Anno_data));
        new_node -> left = NULL;
        new_node -> right = NULL;
        strcpy(new_node -> gene_id, "NA");
        strcpy(new_node -> gene_name, h27ac_eqtl -> gene_name);
        if (add_tree_node_2(first_node, new_node)){
            free(new_node);
        }
        h27ac_eqtl = h27ac_eqtl -> next;
    }

    while(k4me1_eqtl){
        new_node = (struct Anno_data *)malloc(sizeof(struct Anno_data));
        new_node -> left = NULL;
        new_node -> right = NULL;
        strcpy(new_node -> gene_id, "NA");
        strcpy(new_node -> gene_name, k4me1_eqtl -> gene_name);
        if (add_tree_node_2(first_node, new_node)){
            free(new_node);
        }
        k4me1_eqtl = k4me1_eqtl -> next;
    }

    while(eqtl_gwas){
        new_node = (struct Anno_data *)malloc(sizeof(struct Anno_data));
        new_node -> left = NULL;
        new_node -> right = NULL;
        strcpy(new_node -> gene_id, "NA");
        strcpy(new_node -> gene_name, eqtl_gwas -> gene_name);
        if (add_tree_node_2(first_node, new_node)){
            free(new_node);
        }

        eqtl_gwas = eqtl_gwas -> next;
    }

    while(pqtl_gwas){
        new_node = (struct Anno_data *)malloc(sizeof(struct Anno_data));
        new_node -> left = NULL;
        new_node -> right = NULL;
        strcpy(new_node -> gene_id, "NA");
        strcpy(new_node -> gene_name, pqtl_gwas -> gene_name);
        if (add_tree_node_2(first_node, new_node)){
            free(new_node);
        }
        pqtl_gwas = pqtl_gwas -> next;
    }

    return first_node;
}

static int
add_tree_node_2(struct Anno_data * des_node, struct Anno_data * new_node)
{
    int cmp_res = 0;
    cmp_res = strcmp(des_node -> gene_name, new_node -> gene_name);

    if (cmp_res < 0){
        if (des_node -> left){
            add_tree_node_2(des_node -> left, new_node);
        } else{
            des_node -> left = new_node;
        }
    } else if (cmp_res > 0){
        if (des_node -> right){
            add_tree_node_2(des_node -> right, new_node);
        } else{
            des_node -> right = new_node;
        }
    } else{
        return 1;
    }
    return 0;
}


static struct gene_mapper *
make_gene_list(struct Anno_data * gene_tree)
{
    struct gene_mapper * head = NULL, * tmp = NULL, * first_node = NULL, \
        * ptr = NULL;
    head = (struct gene_mapper *)malloc(sizeof(struct gene_mapper));
    strcpy(head -> gene_name, "NA");
    head -> caqtl = NULL;
    head -> mqtl = NULL;
    head -> h27ac = NULL;
    head -> k4me1 = NULL;
    head -> eqtl = NULL;
    head -> pqtl = NULL;
    add_gene(head, gene_tree);

    head = head -> next;
    char gene_name[256];
    while(head){
        strcpy(gene_name, head -> gene_name);
        //printf(">%s\n", gene_name);
        if (!strcmp(gene_name, "Y_RNA") || !strcmp(gene_name, "U6") || \
            !strcmp(gene_name, "U4") || !strcmp(gene_name, "7SK") || \
            !strcmp(gene_name, "Metazoa_SRP")){
            tmp = head -> next;
            free(head);
        } else{
            if (first_node){
                tmp = head -> next;
                head -> next = NULL;
                ptr -> next = head;
                ptr = head;

            } else{
                tmp = head -> next;
                head -> next = NULL;
                first_node = ptr = head;
            }
        }
        head = tmp;
    }

    return first_node;
}

static struct gene_mapper *
add_gene(struct gene_mapper * gene_list, struct Anno_data * gene_name_tree)
{
    struct gene_mapper * new_node;
    if (gene_name_tree){
        new_node = (struct gene_mapper *)malloc(sizeof(struct gene_mapper));
        strcpy(new_node -> gene_name, gene_name_tree -> gene_name);
        new_node -> caqtl = NULL;
        new_node -> mqtl = NULL;
        new_node -> h27ac = NULL;
        new_node -> k4me1 = NULL;
        new_node -> eqtl = NULL;
        new_node -> pqtl = NULL;
        new_node -> next = NULL;
        gene_list -> next = new_node;
        gene_list = new_node;
        gene_list = add_gene(gene_list, gene_name_tree -> left);
        gene_list = add_gene(gene_list, gene_name_tree -> right);
    }
    return gene_list;
}


static void
collect_gene_blongs(struct gene_mapper * gene_chain, void * data_in, int type)
{
    struct eqtl_gwas_data * eqtl_dt;
    struct gwas_data * gwas_dt;
    struct linker * first_node = NULL, * new_node = NULL, * ptr = NULL;
    if (type == 1){
        while(gene_chain){
            first_node = new_node = ptr = NULL;
            eqtl_dt = data_in;
            while(eqtl_dt){
                if (!strcmp(gene_chain -> gene_name, eqtl_dt -> gene_name)){
                    new_node = (struct linker *)malloc(sizeof(struct linker));
                    new_node -> next = NULL;
                    new_node -> ptr = eqtl_dt;
                    if (first_node){
                        ptr -> next = new_node;
                        ptr = new_node;
                    } else{
                        first_node = ptr = new_node;
                    }
                }

                eqtl_dt = eqtl_dt -> next;
            }
            gene_chain -> caqtl = first_node;
            gene_chain = gene_chain -> next;
        }
    } else if (type == 2){
        while(gene_chain){
            first_node = new_node = ptr = NULL;
            eqtl_dt = data_in;
            while(eqtl_dt){
                if (!strcmp(gene_chain -> gene_name, eqtl_dt -> gene_name)){
                    new_node = (struct linker *)malloc(sizeof(struct linker));
                    new_node -> next = NULL;
                    new_node -> ptr = eqtl_dt;
                    if (first_node){
                        ptr -> next = new_node;
                        ptr = new_node;
                    } else{
                        first_node = ptr = new_node;
                    }
                }

                eqtl_dt = eqtl_dt -> next;
            }
            gene_chain -> mqtl = first_node;
            gene_chain = gene_chain -> next;
        }
    } else if(type == 3){
        while(gene_chain){
            first_node = new_node = ptr = NULL;
            eqtl_dt = data_in;
            while(eqtl_dt){
                if (!strcmp(gene_chain -> gene_name, eqtl_dt -> gene_name)){
                    new_node = (struct linker *)malloc(sizeof(struct linker));
                    new_node -> next = NULL;
                    new_node -> ptr = eqtl_dt;
                    if (first_node){
                        ptr -> next = new_node;
                        ptr = new_node;
                    } else{
                        first_node = ptr = new_node;
                    }
                }

                eqtl_dt = eqtl_dt -> next;
            }
            gene_chain -> h27ac = first_node;
            gene_chain = gene_chain -> next;
        }
    } else if (type == 4){
        while(gene_chain){
            first_node = new_node = ptr = NULL;
            eqtl_dt = data_in;
            while(eqtl_dt){
                if (!strcmp(gene_chain -> gene_name, eqtl_dt -> gene_name)){
                    new_node = (struct linker *)malloc(sizeof(struct linker));
                    new_node -> next = NULL;
                    new_node -> ptr = eqtl_dt;
                    if (first_node){
                        ptr -> next = new_node;
                        ptr = new_node;
                    } else{
                        first_node = ptr = new_node;
                    }
                }

                eqtl_dt = eqtl_dt -> next;
            }
            gene_chain -> k4me1 = first_node;
            gene_chain = gene_chain -> next;
        }
    } else if (type == 5){
        while(gene_chain){
            first_node = new_node = ptr = NULL;
            gwas_dt = data_in;
            while(gwas_dt){
                if (!strcmp(gene_chain -> gene_name, gwas_dt -> gene_name)){
                    new_node = (struct linker *)malloc(sizeof(struct linker));
                    new_node -> next = NULL;
                    new_node -> ptr = gwas_dt;
                    if (first_node){
                        ptr -> next = new_node;
                        ptr = new_node;
                    } else{
                        first_node = ptr = new_node;
                    }
                }

                gwas_dt = gwas_dt -> next;
            }
            gene_chain -> eqtl = first_node;
            gene_chain = gene_chain -> next;
        }
    } else if (type == 6){
        while(gene_chain){
            first_node = new_node = ptr = NULL;
            gwas_dt = data_in;
            while(gwas_dt){
                if (!strcmp(gene_chain -> gene_name, gwas_dt -> gene_name)){
                    new_node = (struct linker *)malloc(sizeof(struct linker));
                    new_node -> next = NULL;
                    new_node -> ptr = gwas_dt;
                    if (first_node){
                        ptr -> next = new_node;
                        ptr = new_node;
                    } else{
                        first_node = ptr = new_node;
                    }
                }

                gwas_dt = gwas_dt -> next;
            }
            gene_chain -> pqtl = first_node;
            gene_chain = gene_chain -> next;
        }
    } else {
        fprintf(stderr, "type not recognized\n");
    }

    return;
}


static struct eqtl_gwas_data *
wash_eqtl(struct eqtl_gwas_data * dt_in)
{
    struct eqtl_gwas_data * dt_out = NULL, * ptr = NULL, * tmp = NULL;
    while(dt_in){
        if (strcmp(dt_in -> gene_name, "NA") && strcmp(dt_in -> probeID, "NA")){
            if (dt_out){
                tmp = dt_in -> next;
                dt_in -> next = NULL;
                ptr -> next = dt_in;
                ptr = dt_in;
            } else{
                tmp = dt_in -> next;
                dt_in -> next = NULL;
                dt_out = ptr = dt_in;
            }
        } else{
            tmp = dt_in -> next;
            free(dt_in);
        }
        dt_in = tmp;
    }

    return dt_out;
}


static struct gwas_data *
wash_gwas(struct gwas_data * dt_in)
{
    struct gwas_data * dt_out = NULL, * ptr = NULL, * tmp = NULL;

    while(dt_in){
        if (strcmp(dt_in -> gene_name, "NA")){
            if (dt_out){
                tmp = dt_in -> next;
                dt_in -> next = NULL;
                ptr -> next = dt_in;
                ptr = dt_in;
            } else{
                tmp = dt_in -> next;
                dt_in -> next = NULL;
                dt_out = ptr = dt_in;
            }
        } else{
            tmp = dt_in -> next;
            free(dt_in);
        }
        dt_in = tmp;
    }
    return dt_out;
}


static void
print_res(struct gene_mapper * gene_chain, const char * file_name)
{
    struct gene_mapper * keeper = NULL;
    keeper = gene_chain;
    FILE * f_out = fopen(file_name, "w");
    if (!f_out){
        fprintf(stderr, "error, %s open failed.\n", file_name);
        exit(1);
    }
    fprintf(f_out, "gene\tChr\tPos\tp_ACAT_gene\tp_ACAT_multi\tp_ACAT_gene_weight \
        \tp_ACAT_multi_weight\tP_caQTL_SMR_top\tP_mQTL_SMR_top\tP_eQTL_SMR_top\t \
        P_pQTL_SMR_top\tP_caQTL_SMR_multiP_mQTL_SMR_multi\tP_eQTL_SMR_multi\tP_pQTL_SMR_multi\n");


    float sum_bxy2 = 0;
    int gene_chr = 0;
    unsigned int gene_bp = 0;
    float caqtl_smr_top = 2;
    float caqtl_smr_multi = 2;
    float mqtl_smr_top = 2;
    float mqtl_smr_multi = 2;
    float h27ac_smr_top = 2;
    float h27ac_smr_multi = 2;
    float k4me1_smr_top = 2;
    float k4me1_smr_multi = 2;
    float eqtl_smr_top = 2;
    float eqtl_smr_multi = 2;
    float pqtl_smr_top = 2;
    float pqtl_smr_multi = 2;

    struct linker * linker_ptr = NULL;
    //add sum_bxy2 up.
    while(gene_chain){
        sum_bxy2 = 0;
        linker_ptr = gene_chain -> caqtl;
        sum_bxy2 += trav_linker(linker_ptr, 1);
        linker_ptr = gene_chain -> mqtl;
        sum_bxy2 += trav_linker(linker_ptr, 1);
        linker_ptr = gene_chain -> h27ac;
        sum_bxy2 += trav_linker(linker_ptr, 1);
        linker_ptr = gene_chain -> k4me1;
        sum_bxy2 += trav_linker(linker_ptr, 1);
        linker_ptr = gene_chain -> eqtl;
        sum_bxy2 += trav_linker(linker_ptr, 2);
        linker_ptr = gene_chain -> pqtl;
        sum_bxy2 += trav_linker(linker_ptr, 2);
        gene_chain = gene_chain -> next;

    }

    //calculate weight value.
    gene_chain = keeper;
    while(gene_chain){
        sum_bxy2 = 0;
        linker_ptr = gene_chain -> caqtl;
        trav_linker_2(linker_ptr, 1, sum_bxy2);
        linker_ptr = gene_chain -> mqtl;
        trav_linker_2(linker_ptr, 1, sum_bxy2);
        linker_ptr = gene_chain -> h27ac;
        trav_linker_2(linker_ptr, 1, sum_bxy2);
        linker_ptr = gene_chain -> k4me1;
        trav_linker_2(linker_ptr, 1, sum_bxy2);
        linker_ptr = gene_chain -> eqtl;
        trav_linker_2(linker_ptr, 2, sum_bxy2);
        linker_ptr = gene_chain -> pqtl;
        trav_linker_2(linker_ptr, 2, sum_bxy2);

        gene_chain = gene_chain -> next;

    }

    gene_chain = keeper;
    struct value_list weight_list;
    weight_list.value = -1;
    weight_list.next = NULL;
    struct value_list p_ACAT_list;
    p_ACAT_list.value = -1;
    p_ACAT_list.next = NULL;
    struct value_list multi_p_list;
    multi_p_list.value = -1;
    multi_p_list.next = NULL;
    struct value_list * ptr_1 = NULL, * ptr_2 = NULL, * ptr_3 = NULL;
    ptr_1 = &weight_list;
    ptr_2 = &p_ACAT_list;
    ptr_3 = &multi_p_list;
    while(gene_chain){
        gene_chr = 0;
        gene_bp = 0;

        caqtl_smr_top = 2;
        caqtl_smr_multi = 2;
        linker_ptr = gene_chain -> caqtl;
        value_link(linker_ptr, &gene_chr, &gene_bp, &caqtl_smr_top, &caqtl_smr_multi, \
            &ptr_1, &ptr_2, &ptr_3);

        mqtl_smr_top = 2;
        mqtl_smr_multi = 2;
        linker_ptr = gene_chain -> mqtl;
        value_link(linker_ptr, &gene_chr, &gene_bp, &mqtl_smr_top, &mqtl_smr_multi, \
            &ptr_1, &ptr_2, &ptr_3);

        h27ac_smr_top = 2;
        h27ac_smr_multi = 2;
        linker_ptr = gene_chain -> h27ac;
        value_link(linker_ptr, &gene_chr, &gene_bp, &h27ac_smr_top, &h27ac_smr_multi, \
            &ptr_1, &ptr_2, &ptr_3);

        k4me1_smr_top = 2;
        k4me1_smr_multi = 2;
        linker_ptr = gene_chain -> k4me1;
        value_link(linker_ptr, &gene_chr, &gene_bp, &k4me1_smr_top, &k4me1_smr_multi, \
            &ptr_1, &ptr_2, &ptr_3);

        eqtl_smr_top = 2;
        eqtl_smr_multi = 2;
        linker_ptr = gene_chain -> eqtl;
        value_link_2(linker_ptr, &eqtl_smr_top, &eqtl_smr_multi, &ptr_1, &ptr_2, &ptr_3);

        pqtl_smr_top = 2;
        pqtl_smr_multi = 2;
        linker_ptr = gene_chain -> pqtl;
        value_link_2(linker_ptr, &pqtl_smr_top, &pqtl_smr_multi, &ptr_1, &ptr_2, &ptr_3);


        float p_ACAT_p = acat_func(p_ACAT_list.next);
        float p_multi_p = acat_func(multi_p_list.next);
        float p_ACAT_p_w = acat_weight_func(weight_list.next, p_ACAT_list.next);
        float p_multi_p_w = acat_weight_func(weight_list.next, multi_p_list.next);
        char chr_str[5];
        if (gene_chr == 23){
            strcpy(chr_str, "X");
        } else if(gene_chr == 24){
            strcpy(chr_str, "Y");
        } else{
            sprintf(chr_str, "%u", gene_chr);
        }
        fprintf(f_out, "%s\t%s\t%u\t%f\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", \
            gene_chain -> gene_name, chr_str, gene_bp, p_ACAT_p, p_multi_p, \
            p_ACAT_p_w, p_multi_p_w, caqtl_smr_top, mqtl_smr_top, eqtl_smr_top, pqtl_smr_top, \
            caqtl_smr_multi, mqtl_smr_multi, eqtl_smr_multi, pqtl_smr_multi);

        gene_chain = gene_chain -> next;

    }



    printf("%f\n", sum_bxy2);

    fclose(f_out);
    return;
}


static float
acat_func(struct value_list * dt_in)
{
    float dt_out = 0;
    int i = 0;
    while(dt_in){
        i++;
        dt_out += tan((0.5 - dt_in -> value) * M_PI);
        dt_in = dt_in -> next;
    }
    dt_out = dt_out / i;
    dt_out = 0.5 - atan(dt_out) / M_PI;

    return dt_out;
}


static float
acat_weight_func(struct value_list * dt_in_w, struct value_list * dt_in_p)
{
    float dt_out = 0;
    while(dt_in_w && dt_in_p){
        dt_out += dt_in_w -> value * tan((0.5 - dt_in_p -> value) * M_PI);
        dt_in_w = dt_in_w -> next;
        dt_in_p = dt_in_p -> next;
    }

    dt_out  = 0.5 - atan(dt_out) / M_PI;

    return dt_out;
}



static void
value_link(struct linker * linker_ptr, int * gene_chr, unsigned int * gene_bp, \
    float * top, float * multi, struct value_list ** ptr_1, struct value_list ** ptr_2, \
    struct value_list ** ptr_3)
{
    struct value_list * new_node_1 = NULL, * new_node_2 = NULL, * new_node_3 = NULL;
    struct eqtl_gwas_data * eqtl_dt;
    while(linker_ptr){
        eqtl_dt = (struct eqtl_gwas_data *)linker_ptr -> ptr;
        *gene_chr = eqtl_dt -> Outco_Chr;
        *gene_bp = eqtl_dt -> Outco_bp;
        if (eqtl_dt -> p_SMR < *top)
            *top = eqtl_dt -> p_SMR;
        if (eqtl_dt -> p_SMR_multi_gwas < *multi)
            *multi = eqtl_dt -> p_SMR_multi_gwas;
        new_node_1 = (struct value_list *)malloc(sizeof(struct value_list));
        new_node_2 = (struct value_list *)malloc(sizeof(struct value_list));
        new_node_3 = (struct value_list *)malloc(sizeof(struct value_list));
        new_node_1 -> value = eqtl_dt -> b_SMR_weight;
        new_node_1 -> next = NULL;
        new_node_2 -> value = eqtl_dt -> p_ACAT;
        new_node_2 -> next = NULL;
        new_node_3 -> value = eqtl_dt -> p_SMR_multi_gwas;
        new_node_3 -> next = NULL;
        (*ptr_1) -> next = new_node_1;
        (*ptr_2) -> next = new_node_2;
        (*ptr_3) -> next = new_node_3;
        *ptr_1 = new_node_1;
        *ptr_2 = new_node_2;
        *ptr_3 = new_node_3;
        linker_ptr = linker_ptr -> next;
    }

    return;
}


static void
value_link_2(struct linker * linker_ptr, float * top, float * multi, \
    struct value_list ** ptr_1, struct value_list ** ptr_2, struct value_list ** ptr_3)
{

    struct value_list * new_node_1 = NULL, * new_node_2 = NULL, * new_node_3 = NULL;
    struct gwas_data * gwas_dt;

    while(linker_ptr){
        gwas_dt = (struct gwas_data *)linker_ptr -> ptr;
        if (gwas_dt -> p_SMR < *top)
            *top = gwas_dt -> p_SMR;
        if (gwas_dt -> p_SMR_multi < *multi)
            *multi = gwas_dt -> p_SMR_multi;
        new_node_1 = (struct value_list *)malloc(sizeof(struct value_list));
        new_node_2 = (struct value_list *)malloc(sizeof(struct value_list));
        new_node_3 = (struct value_list *)malloc(sizeof(struct value_list));
        new_node_1 -> value = gwas_dt -> b_SMR_weight;
        new_node_1 -> next = NULL;
        new_node_2 -> value = gwas_dt -> p_ACAT;
        new_node_2 -> next = NULL;
        new_node_3 -> value = gwas_dt -> p_SMR_multi;
        new_node_3 -> next = NULL;
        (*ptr_1) -> next = new_node_1;
        (*ptr_2) -> next = new_node_2;
        (*ptr_3) -> next = new_node_3;
        *ptr_1 = new_node_1;
        *ptr_2 = new_node_2;
        *ptr_3 = new_node_3;
        linker_ptr = linker_ptr -> next;
    }

    return;

}



static float
trav_linker(struct linker * linker_ptr, int type)
{
    float dt_out = 0;
    struct eqtl_gwas_data * eqtl_dt = NULL;
    struct gwas_data * gwas_dt = NULL;

    if (type == 1){
        while(linker_ptr){
            eqtl_dt = (struct eqtl_gwas_data *)linker_ptr -> ptr;
            dt_out += pow(eqtl_dt -> b_SMR, 2);
            linker_ptr = linker_ptr -> next;
        }

    } else if (type == 2){
        while(linker_ptr){
            gwas_dt = (struct gwas_data *)linker_ptr -> ptr;
            dt_out += pow(gwas_dt -> b_SMR, 2);
            linker_ptr = linker_ptr -> next;
        }

    } else{
        fprintf(stderr, "type not recognized.\n");
        exit(1);
    }

    return dt_out;
}



static void
trav_linker_2(struct linker * linker_ptr, int type, float sum_bxy2)
{
    struct eqtl_gwas_data * eqtl_dt = NULL;
    struct gwas_data * gwas_dt = NULL;

    if (type == 1){
        while(linker_ptr){
            eqtl_dt = (struct eqtl_gwas_data *)linker_ptr -> ptr;
            eqtl_dt -> b_SMR_weight = pow(eqtl_dt -> b_SMR, 2) / sum_bxy2;
            linker_ptr = linker_ptr -> next;
        }

    } else if (type == 2){
        while(linker_ptr){
            gwas_dt = (struct gwas_data *)linker_ptr -> ptr;
            gwas_dt -> b_SMR_weight = pow(gwas_dt -> b_SMR, 2) / sum_bxy2;
            linker_ptr = linker_ptr -> next;
        }

    } else{
        fprintf(stderr, "type not recognized.\n");
        exit(1);
    }

    return;
}
//
