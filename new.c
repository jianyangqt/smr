#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>

#define UTILITY
#define ELE_PER_BYTE 4
#define MASK1 3
#define MASK2 12
#define MASK3 48
#define MASK4 192
#define SNP_NUM_INDEX 2
#define PROB_NUM_INDEX 3

#define  ALLEL_TYPE(al) ((al == 1)? 'N': ((al == 0)? '2': ((al == 2)? '1': '0')))

struct BIM_NODE{
    unsigned char chrome;
    char name[25];
    double distance;
    long position;
    unsigned char al1;
    unsigned char al2;
    struct BIM_NODE * next;

};

struct BIM_RES {
    struct BIM_NODE * bim_data;
    unsigned long bim_len;

};

struct FAM_NODE{
    char family_id[10];
    char individual_id[10];
    char paternal_id[10];
    char maternal_id[10];
    char sex;
    signed char phenotype;
    struct FAM_NODE * next;

};

struct FAM_RES {
    struct FAM_NODE * fam_data;
    unsigned long fam_len;

};


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



static struct BIM_RES  parse_bim(const char *);
static struct FAM_RES  parse_fam(const char *);
static void * parse_bed(const char *, unsigned long, unsigned long);
static void unpack_byte(unsigned char, unsigned char *);
static struct BESD_DATA decode_besd_file(const char *);


#ifdef UTILITY
int
main(int argc, char ** argv)
{
    struct BIM_RES  bim_data;
    struct FAM_RES  fam_data;
    const char * bim_file = argv[1];
    const char * fam_file = argv[2];
    const char * bed_file = argv[3];
    const char * besd_file = argv[4];
    unsigned char (*bed_data)[fam_data.fam_len];
    struct BESD_DATA besd_data;

    bim_data = parse_bim(bim_file);
    fam_data = parse_fam(fam_file);
    //printf("%ld %ld\n", bim_data.bim_len, fam_data.fam_len);
    bed_data = parse_bed(bed_file, bim_data.bim_len, fam_data.fam_len); 
    /*
    unsigned long i, j;
    for (i = 0; i < bim_data.bim_len; i++){
        for (j = 0; j < fam_data.fam_len; j++){
            printf("%c", bed_data[i][j]);
        }
        printf("\n");
    }
    */

    besd_data = decode_besd_file(besd_file); 


    return 0;

}
#endif

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






static struct BIM_RES
parse_bim(const char * bim_file_name)
{
    char chrome_tmp[5];
    char name[25];
    double distance;
    long position;
    unsigned char al1;
    unsigned char al2;
    struct BIM_NODE * data_out_ = NULL, * new_n = NULL, * pt = NULL;
    struct BIM_RES data_out;
    unsigned long bim_len = 0;

    FILE * bim_f_in = fopen(bim_file_name, "r");
    if (!bim_f_in){
        fprintf(stderr, "error, bim file not open successesfully.\n");
        exit(EXIT_FAILURE);
    }
    
    while(fscanf(bim_f_in, "%s %s %lf %ld %c %c", chrome_tmp, name, &distance, &position, &al1, &al2) == 6){
        new_n = (struct BIM_NODE *)malloc(sizeof(struct BIM_NODE));
        bim_len++;
        if (isdigit(chrome_tmp[0])){
            new_n -> chrome = atoi(chrome_tmp);
        } else if (strcmp(chrome_tmp, "X")){
            new_n -> chrome = 99;
        } else if (strcmp(chrome_tmp, "Y")) {
            new_n -> chrome = 98;

        } else if (strcmp(chrome_tmp, "XY")){

            new_n -> chrome = 97;
        } else {

            fprintf(stderr, "error, %s is not recognized\n",chrome_tmp);
        }

        strcpy(new_n -> name, name);
        new_n -> distance = distance;
        new_n -> position = position;
        new_n -> al1 = al1;
        new_n -> al2 = al2;
        new_n -> next =NULL;

        if (data_out_){
            pt -> next = new_n;
            pt = new_n;

        } else {
            data_out_ = pt = new_n;
        }
    }
    data_out.bim_data = data_out_;
    data_out.bim_len = bim_len;

    return data_out;
}



static struct FAM_RES
parse_fam(const char * fam_file_name)
{
    char family_id[10];
    char individual_id[10];
    char paternal_id[10];
    char maternal_id[10];
    int sex;
    int phenotype;
    struct FAM_NODE * data_out_ = NULL, * new_ = NULL, * pt = NULL;
    struct FAM_RES data_out;
    unsigned long fam_len = 0;

    FILE * f_in = fopen(fam_file_name, "r");
    if (!f_in){
        fprintf(stderr, "error, fam file open failed.\n");
        exit(EXIT_FAILURE);

    }

    while(fscanf(f_in, "%s %s %s %s %d %d", family_id, individual_id, \
        paternal_id, maternal_id, &sex, &phenotype) == 6){
        new_ = (struct FAM_NODE *)malloc(sizeof(struct FAM_NODE));
        fam_len++;
        strcpy(new_ -> family_id, family_id);
        strcpy(new_ -> individual_id, individual_id);
        strcpy(new_ -> paternal_id, paternal_id);
        strcpy(new_ -> maternal_id, maternal_id);
        new_ ->  sex = (char)sex;
        new_ -> phenotype = (signed char)phenotype;
        new_ -> next = NULL;

        if (data_out_){
            pt -> next = new_;
            pt = new_;
        } else{
            data_out_ = pt = new_;
        }
    }
    data_out.fam_data = data_out_;
    data_out.fam_len = fam_len;
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
 * The function reture a char 2d array; row_num = indivadual_num, clolum_num = snp_amount.
 * data_out[individual_amount][snp_amount]
 */
static void *
parse_bed(const char * bed_file_name, unsigned long row_num, unsigned long clo_num)
{
    unsigned char (* data_out)[clo_num] = (unsigned char (*)[clo_num]) \
        malloc(sizeof(unsigned char) * row_num * clo_num);
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
