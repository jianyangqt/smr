#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


#define ELE_PER_BYTE 4
#define MASK1 3
#define MASK2 12
#define MASK3 48
#define MASK4 192
#define SNP_NUM_INDEX 2
#define PROB_NUM_INDEX 3


#define  ALLEL_TYPE(al) ((al == 1)? 'N': ((al == 0)? '2': ((al == 2)? '1': '0')))


/*
    fam file contain first six column of ped file, which is:
        Family_ID, Individual_ID, Paternal_ID, Maternal_ID,
        Sex(1=male, 2=female, other=unknow), Phenotype
 */
// data structure used to contain fam data.
struct FAM_RES {
    struct FAM_NODE {
        char family_id[32];
        char individual_id[32];
        char paternal_id[32];
        char maternal_id[32];
        char sex;
        int phenotype;
        struct FAM_NODE * next;

    } * fam_data;

    unsigned long fam_len;

};


// fam file parsing function.
static struct FAM_RES
parse_fam(const char * fam_file_name)
{
    char family_id[32];
    char individual_id[32];
    char paternal_id[32];
    char maternal_id[32];
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
        new_node -> next = NULL;
        fam_len++;
        strncpy(new_node -> family_id, family_id, 31);
        (new_node -> family_id)[31] = '\0';
        strncpy(new_node -> individual_id, individual_id, 31);
        (new_node -> family_id)[31] = '\0';
        strncpy(new_node -> paternal_id, paternal_id, 31);
        (new_node -> paternal_id)[31] = '\0';
        strcpy(new_node -> maternal_id, maternal_id);
        (new_node -> maternal_id)[31] = '\0';
        new_node ->  sex = (char)sex;
        new_node -> phenotype = phenotype;

        if (fam_node_data){
            pt -> next = new_node;
            pt = new_node;
        } else{
            fam_node_data = pt = new_node;
        }
    }

    if (!feof(f_in)){
        fprintf(stderr, "error, line of file not all readed.\n");
        exit(1);
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
    set to 23, Y set to 24.
 */
struct BIM_RES {
    struct BIM_NODE {
        unsigned char chrome;
        char snp_name[32];
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
    char snp_name[32];
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
            new_node -> chrome = 23;
        } else if (strcmp(chrome_tmp, "Y")) {
            new_node -> chrome = 24;
        } else {
            fprintf(stderr, "error, %s is not recognized\n",chrome_tmp);
        }

        strncpy(new_node -> snp_name, snp_name, 31);
        (new_node -> snp_name)[31] = '\0';
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

    if (!feof(bim_f_in)){
        fprintf(stderr, "error, line of file not all readed.\n");
        exit(1);

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
parse_bed(const char * bed_file_name, unsigned long row_num, unsigned long row_len)
{
    unsigned char (* data_out)[row_len] = (unsigned char (*)[row_len]) \
        malloc(sizeof(unsigned char) * row_num * row_len);
    //byte length for each raw.
    unsigned long packed_len = (int)ceil((double)row_len / ELE_PER_BYTE);
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

        for (clo_c = 0; clo_c < row_len; clo_c++){
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


//here the matrix element is char.
void *
transformation_matrix(void * matrix, unsigned int row_num, unsigned int row_len)
{
    char (* matrix_now)[row_len] = (char (*)[row_len])matrix;
    char * dt_out = (char *)malloc(sizeof(char) * row_num * row_len);
    unsigned int i, j;


    for (i = 0; i  < row_len; i++){
        for (j = 0; j < row_num; j++){
            dt_out[i * j] = matrix_now[j][i];
        }
    } 

    return dt_out;
}


static void
print_ped_and_map_file(struct FAM_RES fam_dt, struct BIM_RES bim_dt, \
    void * bed_dt, const char * out_file)
{
    unsigned long i, j;
    char ped_file[256];
    char map_file[256];
    strncpy(ped_file, out_file, 251);
    strncpy(map_file, out_file, 251);
    strcat(ped_file, ".ped");
    strcat(map_file, ".map");
    struct FAM_NODE * fam_node;
    struct BIM_NODE * bim_node;

    FILE * f_out_ped = fopen(ped_file, "w");
    FILE * f_out_map = fopen(map_file, "w");

    char (* bed_dt_transed)[bim_dt.bim_len] = \
        transformation_matrix(bed_dt, bim_dt.bim_len, fam_dt.fam_len);
    fam_node = fam_dt.fam_data;
    bim_node = bim_dt.bim_data;
    for (i = 0; i < fam_dt.fam_len; i++){
        fprintf(f_out_ped, "%s\t%s\t%s\t%s\t%d\t%d\t", fam_node -> family_id, fam_node -> individual_id, \
            fam_node -> paternal_id, fam_node -> maternal_id, fam_node -> sex, \
            fam_node -> phenotype);
        for (j = 0; j < bim_dt.bim_len; j++){
                                                        
            


        }
        fam_node = fam_node -> next;


    }

}



void
plink_decode(const char * plink_file_prefix, const char * out_file)
{
    char fam_file_name[256];
    strncpy(fam_file_name, plink_file_prefix, 251);
    strcat(fam_file_name, ".fam");
    char bim_file_name[256];
    strncpy(bim_file_name, plink_file_prefix, 251);
    strcat(bim_file_name, ".bim");
    char bed_file_name[256];
    strncpy(bed_file_name, plink_file_prefix, 251);
    strcat(bed_file_name, ".bed");

    printf("%s %s %s \n", fam_file_name, bim_file_name, bed_file_name);
    struct FAM_RES fam_dt = parse_fam(fam_file_name);
    struct BIM_RES bim_dt = parse_bim(bim_file_name);
    void * bed_dt = parse_bed(bed_file_name, bim_dt.bim_len, fam_dt.fam_len);
    print_ped_and_map_file(fam_dt, bim_dt, bed_dt, out_file);

}







void
plink_encode(const char * plink_file_prefix, const char * out_file)
{
    printf("The plink file encode is under development.\n");
    return;

}



int
main(int argc, char * argv[])
{

    char * pref = argv[1];

    plink_decode(pref, "hhh");

    return 0;
}




