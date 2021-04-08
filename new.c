#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>

#define UTILITY
#define ELE_PER_BYTE 4
#define MASK1 3
#define MASK2 12
#define MASK3 48
#define MASK4 192

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


static struct BIM_RES  parse_bim(char *);
static struct FAM_RES  parse_fam(char *);
static void * parse_bed(char *, unsigned long, unsigned long);
static void unpack_byte(unsigned char, unsigned char *);


#ifdef UTILITY
int
main(int argc, char ** argv)
{
    
    struct BIM_RES  bim_data;
    struct FAM_RES  fam_data;
    char * bim_file;
    char * fam_file;
    char * bed_file;

    bim_file = argv[1];
    fam_file = argv[2];
    bed_file = argv[3];
    bim_data = parse_bim(bim_file);
    fam_data = parse_fam(fam_file);
    printf("%ld %ld\n", bim_data.bim_len, fam_data.fam_len);

    unsigned char (*bed_data)[fam_data.fam_len];
    bed_data = parse_bed(bed_file, bim_data.bim_len, fam_data.fam_len); 

    unsigned long i, j;
    for (i = 0; i < bim_data.bim_len; i++){
        for (j = 0; j < fam_data.fam_len; j++){
            printf("%c", bed_data[i][j]);
        }
        printf("\n");
    }

    return 0;

}
#endif


static struct BIM_RES
parse_bim(char * bim_file_name)
{
    char chrome_tmp[5];
    char name[25];
    double distance;
    long position;
    unsigned char al1;
    unsigned char al2;
    struct BIM_NODE * data_out_ = NULL, * new_n = NULL, * pt = NULL;
    struct BIM_RES data_out;
    unsigned long bim_len;

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
parse_fam(char * fam_file_name)
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
 */
static void *
parse_bed(char * bed_file_name, unsigned long row_num, unsigned long clo_num)
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
