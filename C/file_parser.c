#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>

//#define TEST_FUNC
#define ELE_PER_BYTE 4
#define MASK1 3
#define MASK2 12
#define MASK3 48
#define MASK4 192
#define SNP_NUM_INDEX 2
#define PROB_NUM_INDEX 3

//define for msmr parsing
#define LINE_BUF_LEN 500

// uint32 + floats
#define DENSE_FILE_TYPE_1 0

// 16*uint32s + floats (indicator+samplesize+snpnumber+probenumber+ 6*-9s + values) [default]]
#define DENSE_FILE_TYPE_3 5

// uint32 + uint64_t + uint64_ts + uint32_ts + floats
#define SPARSE_FILE_TYPE_3F 0x40400000

// 16*uint32s + uint64_t + uint64_ts + uint32_ts +
// floats (indicator+samplesize+snpnumber+probenumber+ 6*-9s +
// valnumber+cols+rowids+betases) [default]]
#define SPARSE_FILE_TYPE_3 3

#define  ALLEL_TYPE(al) ((al == 1)? 'N': ((al == 0)? '2': ((al == 2)? '1': '0')))


/*
    fam file contain first six column of ped file, which is:
        Family_ID, Individual_ID, Paternal_ID, Maternal_ID,
        Sex(1=male, 2=female, other=unknow), Phenotype
 */
// data structure used to contain fam data.
struct FAM_RES {
    struct FAM_NODE {
        char family_id[10];
        char individual_id[10];
        char paternal_id[10];
        char maternal_id[10];
        char sex;
        signed char phenotype;
        struct FAM_NODE * next;

    } * fam_data;

    unsigned long fam_len;

};


// fam file parsing function.
static struct FAM_RES
parse_fam(const char * fam_file_name)
{
    char family_id[10];
    char individual_id[10];
    char paternal_id[10];
    char maternal_id[10];
    int sex;
    int phenotype;
    struct FAM_NODE * fam_node_data = NULL, * new_node = NULL, * pt = NULL;
    struct FAM_RES data_out;
    unsigned long fam_len = 0;

    FILE * f_in = fopen(fam_file_name, "r");
    if (!f_in){
        fprintf(stderr, "error, fam file open failed.\n");
        exit(EXIT_FAILURE);
    }

    while(fscanf(f_in, "%s %s %s %s %d %d", family_id, individual_id, \
        paternal_id, maternal_id, &sex, &phenotype) == 6){
        new_node = (struct FAM_NODE *)malloc(sizeof(struct FAM_NODE));
        fam_len++;
        strcpy(new_node -> family_id, family_id);
        strcpy(new_node -> individual_id, individual_id);
        strcpy(new_node -> paternal_id, paternal_id);
        strcpy(new_node -> maternal_id, maternal_id);
        new_node ->  sex = (char)sex;
        new_node -> phenotype = (signed char)phenotype;
        new_node -> next = NULL;

        if (fam_node_data){
            pt -> next = new_node;
            pt = new_node;
        } else{
            fam_node_data = pt = new_node;
        }
    }
    data_out.fam_data = fam_node_data;
    data_out.fam_len = fam_len;
    return data_out;
}


/*
    bim file is extended MAP file, two extra columns is allels name,
    the columns is described as fllowing:
        chromosome, rs or snp, Genetic distance, Base pair position, ref allel, alt allel
    here I conver chromosome to a unsigned int number, 1 - 22 will be itself, and X is
    set to 99, Y set to 98, XY set to 97.
 */
struct BIM_RES {
    struct BIM_NODE {
        unsigned char chrome;
        char snp_name[25];
        double distance;
        unsigned long position;
        unsigned char al1;
        unsigned char al2;
        struct BIM_NODE * next;

    } * bim_data;

    unsigned long bim_len;

};


//parse bim file.
static struct BIM_RES
parse_bim(const char * bim_file_name)
{
    char chrome_tmp[5];
    char snp_name[25];
    double distance;
    unsigned long position;
    unsigned char al1;
    unsigned char al2;
    struct BIM_NODE * bim_node_data = NULL, * new_node = NULL, * pt = NULL;
    struct BIM_RES data_out;
    unsigned long bim_len = 0;

    FILE * bim_f_in = fopen(bim_file_name, "r");
    if (!bim_f_in){
        fprintf(stderr, "error, bim file not open successesfully.\n");
        exit(EXIT_FAILURE);
    }

    while(fscanf(bim_f_in, "%s %s %lf %ld %c %c", chrome_tmp, snp_name, \
        &distance, &position, &al1, &al2) == 6){
        new_node = (struct BIM_NODE *)malloc(sizeof(struct BIM_NODE));
        bim_len++;
        if (isdigit(chrome_tmp[0])){
            new_node -> chrome = atoi(chrome_tmp);
        } else if (strcmp(chrome_tmp, "X")){
            new_node -> chrome = 99;
        } else if (strcmp(chrome_tmp, "Y")) {
            new_node -> chrome = 98;
        } else if (strcmp(chrome_tmp, "XY")){
            new_node -> chrome = 97;
        } else {
            fprintf(stderr, "error, %s is not recognized\n",chrome_tmp);
        }

        strcpy(new_node -> snp_name, snp_name);
        new_node -> distance = distance;
        new_node -> position = position;
        new_node -> al1 = al1;
        new_node -> al2 = al2;
        new_node -> next = NULL;

        if (bim_node_data){
            pt -> next = new_node;
            pt = new_node;

        } else {
            bim_node_data = pt = new_node;
        }
    }
    data_out.bim_data = bim_node_data;
    data_out.bim_len = bim_len;

    return data_out;
}


/*
 * .bed file was generate by plink, it was is a binary file and used to record genotype data.
 * every allel was recorded using 2bit, the meaning is discriped as follow:
 *  00 homozygote
 *  10 heterozygous
 *  11 homozygote
 *  01 missing
 *
 * When readin the bed file, the very first 3 byte should trimed, it is not belong to item data.
 * And, the row of data is represent for all individual of a snp/other, if the number of individual
 * cann't divided by 4 integely, zero using to fill left bit of a byte.
 *
 * As reference to Plinks docs and Plink2R, I using follwing character to represent decode results:
 *
 *  00 homozygote     -> '2'
 *  10 heterozygous   -> '1'
 *  11 homozygote     -> '0'
 *  01 missing        -> 'N'
 *
 * The function reture a 2 dimension char array; row_length = individual_num,
 * column_length = snp_amount. In other words, row_number = snp_amount,
 * column_number = individual_amount.
 *
 * data_out[snp_amount][individual_amount]
 */
static void unpack_byte(unsigned char, unsigned char *);
static void *
parse_bed(const char * bed_file_name, unsigned long row_num, unsigned long clo_num)
{
    unsigned char (* data_out)[clo_num] = (unsigned char (*)[clo_num]) \
        malloc(sizeof(unsigned char) * row_num * clo_num);
    //byte length for each raw.
    unsigned long packed_len = (int)ceil((double)clo_num / ELE_PER_BYTE);
    unsigned char packed_buffer[packed_len];
    unsigned char unpacked_buffer[packed_len * ELE_PER_BYTE];
    unsigned long row_c = 0;
    unsigned long clo_c = 0;
    unsigned long i = 0, j = 0, k = 0;
    unsigned char unpacked_byte[ELE_PER_BYTE];

    FILE * f_in = fopen(bed_file_name, "r");
    if (!f_in){
        fprintf(stderr, "error, open bed file failed.\n");
        exit(1);
    }
    //remove first 3 byte.
    i = 0;
    while(fgetc(f_in)){
        i++;
        if (i == 3)
            break;
    }

    while (fread(packed_buffer, sizeof(unsigned char), packed_len, f_in) == packed_len){
        j = 0;
        for (i = 0; i < packed_len; i++){
            unpack_byte(packed_buffer[i], unpacked_byte);
            k = 0;
            while (k < ELE_PER_BYTE){
                unpacked_buffer[j] = unpacked_byte[k];
                j++;
                k++;
            }
        }

        for (clo_c = 0; clo_c < clo_num; clo_c++){
            data_out[row_c][clo_c] = unpacked_buffer[clo_c];
        }
        row_c++;
    }

    return data_out;
}


static void
unpack_byte(unsigned char packed_char, unsigned char * unpacked_byte)
{
    //the length is ELE_PER_BYTE, which is 4 here.
    unsigned char f1, f2, f3, f4;
    f1 = packed_char & MASK1;
    f1 = ALLEL_TYPE(f1);

    f2 = (packed_char & MASK2) >> 2;
    f2 = ALLEL_TYPE(f2);

    f3 = (packed_char & MASK3) >> 4;
    f3 = ALLEL_TYPE(f3);

    f4 = (packed_char & MASK4) >> 6;
    f4 = ALLEL_TYPE(f4);

    unpacked_byte[0] = f1;
    unpacked_byte[1] = f2;
    unpacked_byte[2] = f3;
    unpacked_byte[3] = f4;
    return;
}




/*
    esi file, which is same as bmi file.
    The column is following:
        chromosome, snp name, genetic distance, position,
        reference allel, alternative allel, frequency of reference allel(effect allel).

    if a field can be NA, then the bigest number of its type will be used to represent
    the NA condition. For example, ESI_NODE.frequency is unsigned long (minimux is uint32),
    then (2**32 - 1) will be used represent NA.

 */
struct ESI_NODE {
    unsigned char chromosome;
    char snp_name[25];
    double distance;
    unsigned long position;
    char al1;
    char al2;
    unsigned long frequency;

};


/* The besd file is a kind of packed file, which store data
 * as binary instead of ASCII to reduce size of file. Here we
 * describe the file as Unit_array.
 * The filestore Beta and SE data. For each prob, there are multi snp which
 * asociated with it, and every snp contain a Beta value and SE value
 * (Values is used to describ this two value).
 * In order to associated Values to Prob(epi file) and SNP(esi file) data.
 * Some mate date was stored in Unit_array in front of Values_array. There
 * are offset of each prob's Beta and SE data block, and Each Value related SNP
 * index of SNP file(esi file).
 *
 * Here is digram to show the data array layout:
 * [first 16 int, <First 16 int>],
 * [Number of Values, 1 Uint64, <Values_arrary_length>],
 * [offset of Beta(block) and SE(block) value of each prob, Prob_num * 2 + 1 Uint64, <Beta_SE_offset>],
 * [SNP index(line number of esi file) of Values, segma(Prob_num, snp) * 2 Uint32, <snp_index>],
 * [Values array, segma(Prob_num, snp) * 2 Float, <Values_array>]
 *
 * Here are details of Unit_array_layout(Unit_array count from zero):
 * int32(0)[file_type], int(1)[-9], int(2)[snp_amount], int(3)[prob_amount],
 * int(4)[-9], ..., int(15)[-9];
 * uint64(16)[number of Values, Values array length];
 * uint64(17)[0], uint64(18)[offset of first prob's Beta values],
 * uint64(20)[offset of first prob's SE value], ...;
 * uint32(m)[], uint32(m+1)[], ...;
 * float(n)[], float(n+1), ... float(n+p1), float(n+p1+1), float(n+p1+2), ... float(n+2*p1);
 * ...
 */
struct BESD_DATA {
    int32_t first_16_int[16];
    uint64_t value_num;
    uint64_t * value_beta_se_pos;
    uint32_t * snp_pos;
    float * value_beta_se;

};

struct TMP_S {
    unsigned long prob_index;
    unsigned int snp_contained;
    struct SNP_DT {
        unsigned int snp_index;
        float beta;
        float se;
    } * SNP_data_pt;
};


static struct BESD_DATA
decode_besd_file(const char * besd_file_name)
{
    struct BESD_DATA data_out;
    data_out.value_beta_se_pos = NULL;
    data_out.snp_pos = NULL;
    data_out.value_beta_se = NULL;
    int file_pack_type = 0;
    int snp_amount = 0;
    int prob_amount = 0;
    int32_t first_16[16];
    uint64_t Value_num = 0;
    unsigned long i = 0;
    unsigned long prob_index = 0;
    unsigned long snp_index = 0;
    unsigned long offset_first = 0, offset_second = 0, offset_third = 0;
    unsigned int beta_num_per_prob = 0, se_num_per_prob = 0;
    float beta_value = 0.0;
    float se_value = 0.0;
    unsigned int snp_index_beta;
    unsigned int snp_index_se;



    FILE * f_in = fopen(besd_file_name, "r");
    if (!f_in){
        fprintf(stderr, "error, open best file failed.\n");
        exit(EXIT_FAILURE);
    }

    if (fread(&first_16, sizeof(int32_t), 16, f_in) != 16){
        exit(1);
    }
    file_pack_type = first_16[0];
    snp_amount = first_16[SNP_NUM_INDEX];
    prob_amount = first_16[PROB_NUM_INDEX];
    if (file_pack_type != 3){
        fprintf(stderr, "error, From now, this prog only can unpack type 3 file.\n");
        exit(1);
    }

    if (fread(&Value_num, sizeof(uint64_t), 1, f_in) != 1){
        exit(1);
    }
    printf("snp:%ld prob:%ld Value:%ld\n",snp_amount, prob_amount, Value_num);

    unsigned long beta_se_offset_length = prob_amount * 2 + 1;
    uint64_t * beta_se_offset = (uint64_t *)malloc(sizeof(uint64_t) * beta_se_offset_length);
    if (fread(beta_se_offset, sizeof(uint64_t), beta_se_offset_length, f_in) != beta_se_offset_length){
        exit(1);
    }
    if (beta_se_offset[beta_se_offset_length - 1] != Value_num){
        fprintf(stderr, "error, the last offset number should equal to value number.\n");
        exit(1);
    }

    unsigned long snp_index_block_len = Value_num;
    uint32_t * snp_index_beta_se = (uint32_t *)malloc(sizeof(uint32_t) * snp_index_block_len);
    if(fread(snp_index_beta_se, sizeof(uint32_t), snp_index_block_len, f_in) != snp_index_block_len){
        exit(1);
    }

    float * Values = (float *)malloc(sizeof(float) * Value_num);
    if (fread(Values, sizeof(float), Value_num, f_in) != Value_num){
        exit(1);
    }

    for (i = 0; i < snp_index_block_len; i++){
        printf("%u ", snp_index_beta_se[i]);
    }
    printf("\n");

    for (i = 0; i < Value_num; i++){
        printf("%f ", Values[i]);
    }
    printf("\n");

    offset_first = beta_se_offset[0];
    for (prob_index = 1; prob_index <= prob_amount; prob_index += 2 ){
        offset_second = beta_se_offset[prob_index];
        offset_third = beta_se_offset[prob_index + 1];
        beta_num_per_prob = offset_second - offset_first;
        se_num_per_prob = offset_third - offset_second;
        if (beta_num_per_prob != se_num_per_prob){
            fprintf(stderr, "beta value number should equal se value number for same prob.\n");
            exit(1);
        }
        printf(">> %ld,%ld\n", beta_num_per_prob, se_num_per_prob);

        for (i = 0; i < beta_num_per_prob; i++){
            snp_index_beta = snp_index_beta_se[offset_first + i];
            snp_index_se = snp_index_beta_se[offset_second + i];
            beta_value = Values[offset_first + i];
            se_value = Values[offset_second + i];
            printf("%u %u %.8f %.8f \n", snp_index_beta, snp_index_se, beta_value, se_value);
        }
        offset_first = offset_third;
    }

    printf("go end\n");
    return data_out;
}




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


static int
convert_chromosome(const char * chromo_str)
{
    if (isdigit(chromo_str[0])){
        return atoi(chromo_str);
    } else if(strcmp(chromo_str, "X") == 0){
        return 23;
    } else if (strcmp(chromo_str, "Y") == 0){
        return 24;
    } else{
        perror("chromosome cann't be covert to intege");
    }
}


static void
assign_na_value(unsigned char * na_value, double * value, const char * field)
{
    if (strcmp(field, "NA") == 0){
        *na_value = 0;
        *value = 0;

    } else{
        *na_value = 1;
        *value = atof(field);

    }
    return;
}


struct MSMR_DATA_eqtl
parse_msmr_file_eqtl(const char * file_name, bool head)
{
    struct MSMR_DATA_eqtl data_out;
    data_out.msmr_node = NULL;
    unsigned int node_len = 0;
    const int field_len = 26;
    char line_buffer[LINE_BUF_LEN];

    char expo_id[50];
    char expo_chr_str[5];
    char expo_gene[50];
    unsigned int expo_bp;

    char outco_id[50];
    char outco_chr_str[5];
    char outco_gene[50];
    unsigned int outco_bp;

    char top_snp[50];
    char top_snp_chr_str[5];
    unsigned int top_snp_bp;

    char a1[10];
    char a2[10];

    char frequence[20];
    char b_outco[20];
    char se_outco[20];
    char p_outco[20];
    char b_expo[20];
    char se_expo[20];
    char p_expo[20];
    char b_smr[20];
    char se_smr[20];
    char p_smr[20];
    char p_smr_multi[20];
    char p_heidi[20];
    char nsnp_heidi[20];

    struct MSMR_NODE_eqtl * new_node = NULL, * first_node = NULL, * pt = NULL;


    FILE * f_in = fopen(file_name, "r");
    if (!f_in){
        fprintf(stderr, "error, open msmr file failed.\n");
        exit(1);
    }

    if (head){
        fgets(line_buffer, LINE_BUF_LEN, f_in);
    }

    node_len = 0;
    while(fscanf(f_in, "%s %s %s %u %s %s %s %u %s %s %u %s %s \
        %s %s %s %s %s %s %s %s %s %s %s %s %s", \
        expo_id, expo_chr_str, expo_gene, &expo_bp, \
        outco_id, outco_chr_str, outco_gene, &outco_bp, \
        top_snp, top_snp_chr_str, &top_snp_bp, \
        a1, a2,\
        frequence, b_outco, se_outco, p_outco, \
        b_expo, se_expo, p_expo, b_smr, se_smr, p_smr, \
        p_smr_multi, p_heidi, nsnp_heidi) == field_len){
        node_len++;

        new_node = (struct MSMR_NODE_eqtl *)malloc(sizeof(struct MSMR_NODE_eqtl));
        new_node -> prev = NULL;
        new_node -> next = NULL;
        strcpy(new_node -> expo_id, expo_id);
        new_node -> expo_chr = convert_chromosome(expo_chr_str);
        strcpy(new_node -> expo_gene, expo_gene);
        new_node -> expo_bp = expo_bp;
        strcpy(new_node -> outco_id, outco_id);
        new_node -> outco_chr = convert_chromosome(outco_chr_str);
        strcpy(new_node -> outco_gene, outco_gene);
        new_node -> outco_bp = outco_bp;
        strcpy(new_node -> top_snp, top_snp);
        new_node -> top_snp_chr = convert_chromosome(top_snp_chr_str);
        new_node -> top_snp_bp = top_snp_bp;
        strcpy(new_node -> a1, a1);
        strcpy(new_node -> a2, a2);

        assign_na_value(&(new_node -> frequence_na), &(new_node -> frequence), frequence);
        assign_na_value(&(new_node -> b_outco_na), &(new_node -> b_outco), b_outco);
        assign_na_value(&(new_node -> se_outco_na), &(new_node -> se_outco), se_outco);
        assign_na_value(&(new_node -> p_outco_na), &(new_node -> p_outco), p_outco);
        assign_na_value(&(new_node -> b_expo_na), &(new_node -> b_expo), b_expo);
        assign_na_value(&(new_node -> se_expo_na), &(new_node -> se_expo), se_expo);
        assign_na_value(&(new_node -> p_expo_na), &(new_node -> p_expo), p_expo);
        assign_na_value(&(new_node -> b_smr_na), &(new_node -> b_smr), b_smr);
        assign_na_value(&(new_node -> se_smr_na), &(new_node -> se_smr), se_smr);
        assign_na_value(&(new_node -> p_smr_na), &(new_node -> p_smr), p_smr);
        assign_na_value(&(new_node -> p_smr_multi_na), &(new_node -> p_smr_multi), p_smr_multi);
        assign_na_value(&(new_node -> p_heidi_na), &(new_node -> p_heidi), p_heidi);
        if (strcmp(nsnp_heidi, "NA") == 0){
            new_node -> nsnp_heidi_na = 0;
            new_node -> nsnp_heidi = 0;

        } else{
            new_node -> nsnp_heidi_na = 1;
            new_node -> nsnp_heidi = (unsigned short)atoi(nsnp_heidi);
        }

        if (first_node){
            new_node -> prev = pt;
            pt -> next = new_node;
            pt = new_node;

        } else{
            first_node = pt = new_node;
        }
    }
    if (!feof(f_in)){
        fclose(f_in);
        fprintf(stderr, "error, msmr read line exit before read all lines.\n");
        exit(1);
    }

    data_out.node_len = node_len;
    data_out.msmr_node = first_node;

    return data_out;
}


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


struct MSMR_DATA_gwas
parse_msmr_file_gwas(const char * file_name, bool head)
{
    struct MSMR_DATA_gwas data_out;
    data_out.msmr_node = NULL;
    char line_buffer[LINE_BUF_LEN];
    const int field_len = 22;
    unsigned int node_len = 0;

    char probe_id[50];
    char probe_chr_str[5];
    char gene[50];
    unsigned int probe_bp;
    char top_snp[50];
    char top_snp_chr_str[5];
    unsigned int top_snp_bp;
    char a1[10];
    char a2[10];

    char frequence[20];
    char b_gwas[20];
    char se_gwas[20];
    char p_gwas[20];
    char b_eqtl[20];
    char se_eqtl[20];
    char p_eqtl[20];
    char b_smr[20];
    char se_smr[20];
    char p_smr[20];
    char p_smr_multi[20];
    char p_heidi[20];
    char nsnp_heidi[20];

    struct MSMR_NODE_gwas * first_node = NULL, * new_node = NULL, * pt = NULL;

    FILE * f_in = fopen(file_name, "r");
    if(!f_in){
        fprintf(stderr, "error, open file error.\n");
        exit(1);
    }
    if (head){
        fgets(line_buffer, LINE_BUF_LEN, f_in);
    }

    node_len = 0;
    while(fscanf(f_in, "%s %s %s %u %s %s %u %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s", \
        probe_id, probe_chr_str, gene, &probe_bp, top_snp, top_snp_chr_str, &top_snp_bp, \
        a1, a2, frequence, b_gwas, se_gwas, p_gwas, b_eqtl, se_eqtl, p_eqtl, b_smr, se_smr, \
        p_smr, p_smr_multi, p_heidi, nsnp_heidi) == field_len){
        
        node_len ++;
        //printf("%u\n", node_len);

        new_node = (struct MSMR_NODE_gwas *)malloc(sizeof(struct MSMR_NODE_gwas));
        //this is the most saft way to store string.
        new_node -> prev = NULL;
        new_node -> next = NULL;
        strncpy(new_node -> probe_id, probe_id, 49);
        (new_node -> probe_id)[49] = '\0';
        new_node -> probe_chr = convert_chromosome(probe_chr_str);
        strncpy(new_node -> gene, gene, 49);
        (new_node -> gene)[49] = '\0';
        new_node -> probe_bp = probe_bp;
        strncpy(new_node -> top_snp, top_snp, 49);
        (new_node -> top_snp)[49] = '\0';
        new_node -> top_snp_chr = convert_chromosome(top_snp_chr_str);
        new_node -> top_snp_bp = top_snp_bp;
        strcpy(new_node -> a1, a1);
        strcpy(new_node -> a2, a2);

        assign_na_value(&(new_node -> frequence_na), &(new_node -> frequence), frequence);
        assign_na_value(&(new_node -> b_gwas_na), &(new_node -> b_gwas), b_gwas);
        assign_na_value(&(new_node -> se_gwas_na), &(new_node -> se_gwas), se_gwas);
        assign_na_value(&(new_node -> p_gwas_na), &(new_node -> p_gwas), p_gwas);
        assign_na_value(&(new_node -> b_eqtl_na), &(new_node -> b_eqtl), b_eqtl);
        assign_na_value(&(new_node -> se_eqtl_na), &(new_node -> se_eqtl), se_eqtl);
        assign_na_value(&(new_node -> p_eqtl_na), &(new_node -> p_eqtl), p_eqtl);
        assign_na_value(&(new_node -> b_smr_na), &(new_node -> b_smr), b_smr);
        assign_na_value(&(new_node -> se_smr_na), &(new_node -> se_smr), se_smr);
        assign_na_value(&(new_node -> p_smr_na), &(new_node -> p_smr), p_smr);
        assign_na_value(&(new_node -> p_smr_multi_na), &(new_node -> p_smr_multi), p_smr_multi);
        assign_na_value(&(new_node -> p_heidi_na), &(new_node -> p_heidi), p_heidi);
        if (strcmp(nsnp_heidi, "NA") == 0){
            new_node -> nsnp_heidi_na = 0;
            new_node -> nsnp_heidi = 0;
        } else{
            new_node -> nsnp_heidi_na = 1;
            new_node -> nsnp_heidi = (unsigned short)atoi(nsnp_heidi);
        }

        if (first_node){
            new_node -> prev = pt;
            pt -> next = new_node;
            pt = new_node;
        } else{
            first_node = pt = new_node;
        }
    }

    if (!feof(f_in)){
        fprintf(stderr, "error, not all line was readed.\n");
        exit(1);
    }

    data_out.node_len = node_len;
    data_out.msmr_node = first_node;

    return data_out;
}


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


static void
assign_na_value_int(unsigned char * na_ptr, unsigned int * value_ptr, const char * field)
{
    if (strcmp(field, "NA") == 0){
        *na_ptr = 0;
        *value_ptr = 0;
    } else{
        *na_ptr = 1;
        *value_ptr = atoi(field);
    }
    return;
}


struct ANNO_DATA
parse_anno_file_msmr(const char * file_name, bool head)
{
    struct ANNO_DATA data_out;
    data_out.data_node = NULL;
    const unsigned char field_len = 7;
    unsigned int node_len = 0;
    
    char chr_str[5];
    char gene_id[50];
    char gene_name[50];
    
    char gene_type[20];
    char gene_start[20];
    char gene_end[20];
    char orientation[5];

    FILE * f_in = fopen(file_name, "r");
    if (!f_in){
        fprintf(stderr, "file open error.\n");
        exit(1);
    }

    struct ANNO_NODE * first_node, * pt, * new_node;
    node_len = 0;
    while(fscanf(f_in, "%s %s %s %s %s %s %s", \
        chr_str, gene_id, gene_name, gene_type, gene_start, gene_end, orientation) == field_len){
        
        node_len ++;
        new_node = (struct ANNO_NODE *)malloc(sizeof(struct ANNO_NODE));
        new_node -> prev = NULL;
        new_node -> next = NULL;

        new_node -> chr = convert_chromosome(chr_str);
        //here has danger of line buffer overflow.
        strcpy(new_node -> gene_id, gene_id);
        strcpy(new_node -> gene_name, gene_name);
        if (strcmp(gene_type, "NA") == 0){
            new_node -> gene_type_na = 0;
            strcpy(new_node -> gene_type, "");

        } else{
            new_node -> gene_type_na = 1;
            strcpy(new_node -> gene_type, gene_type);
        }
        assign_na_value_int(&(new_node -> gene_start_na), &(new_node -> gene_start), gene_start);
        assign_na_value_int(&(new_node -> gene_end_na), &(new_node -> gene_end), gene_end);
        if (strcmp(gene_type, "NA") == 0){
            new_node -> orientation_na = 0;
            new_node -> orientation = 0;
        } else{
            new_node -> orientation_na = 1;
            new_node -> orientation = orientation[0];
        }

        if (first_node){
            new_node -> prev = pt;
            pt -> next = new_node;

        } else{
            first_node = pt = new_node;
        }

    }

    if (!feof(f_in)){
        fprintf(stderr, "read line exit early.\n");
        exit(1);
    }

    data_out.node_len = node_len;
    data_out.data_node = first_node;

    return data_out;

}



//Hash function using short.
#define HASH_LEN_STR 10
unsigned short
hash_func(const char * short_str)
{
    unsigned short dt_out = 0;
    unsigned char j = 0;
    for (int i = 0; i < strlen(short_str); i++){
        if (j == 0){
            j = 1;
            dt_out ^= short_str[i] << 8;

        } else{
            j = 0;
            dt_out ^= short_str[i];
        }
    }

    return dt_out;
}


#ifdef TEST_FUNC
int
main(int argc, char * argv[])
{
    FILE * f_in = fopen(argv[1], "r");
    unsigned short hv;
    char a[50], b[50], c[50], d[50], e[50], f[50], g[50];

    while(fscanf(f_in, "%s %s %s %s %s %s %s", a, b, c, d, e, f, g) == 7){

        hv = hash_func(b);
        printf("%hu\n", hv);


    }

    return 0;
}
#endif
