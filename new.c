#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#define UTILITY
#define MASK1 3
#define MASK2 12
#define MASK3 48
#define MASK4 192

struct BIM_NODE{
    unsigned char chrome;
    char name[25];
    double distance;
    long position;
    unsigned char al1;
    unsigned char al2;
    struct BIM_NODE * next;

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

struct BED_NODE{
    unsigned char * genotype_raw;
    struct BED_NODE * next;

};


static struct BIM_NODE * parse_bim(char *);
static struct FAM_NODE * parse_fam(char *);
static struct BED_NODE * parse_bed(char *, unsigned int);

#ifdef UTILITY
int
main(int argc, char ** argv)
{
    
    struct BIM_NODE * bim_data;
    struct BIM_NODE * bim_data_tmp;
    struct FAM_NODE * fam_data;
    struct FAM_NODE * fam_data_tmp;
    struct BED_NODE * bed_data;
    struct BED_NODE * bed_data_tmp;
    char * bim_file;
    char * fam_file;
    char * bed_file;
    unsigned int bim_len = 0;
    unsigned int i = 0;
    bim_file = argv[1];
    fam_file = argv[2];
    bed_file = argv[3];
    bim_data = parse_bim(bim_file);
    fam_data = parse_fam(fam_file);

    
    bim_data_tmp = bim_data;
    while(bim_data_tmp){
        bim_len ++;
        bim_data_tmp = bim_data_tmp -> next;
    }
    //printf("%d\n", bim_len);

    fam_data_tmp = fam_data;
    while(fam_data_tmp){
        //printf("%s %s\n", fam_data_tmp -> family_id, fam_data_tmp -> individual_id);
        fam_data_tmp = fam_data_tmp -> next;
    }


    bed_data = parse_bed(bed_file, bim_len);
    bed_data_tmp = bed_data;

    while (bed_data_tmp){
        i = 0;
        while (i < bim_len){

            printf("%d ", (bed_data_tmp -> genotype_raw)[i]);
            i++;
        }

        printf("\n");

        bed_data_tmp = bed_data_tmp -> next;

    }    
    return 0;

}
#endif


static struct BIM_NODE *
parse_bim(char * bim_file_name)
{
    char chrome_tmp[5];
    char name[25];
    double distance;
    long position;
    unsigned char al1;
    unsigned char al2;
    struct BIM_NODE * data_out = NULL, * new_n = NULL, * pt = NULL;

    FILE * bim_f_in = fopen(bim_file_name, "r");
    if (!bim_f_in){
        fprintf(stderr, "error, bim file not open successesfully.\n");
        exit(EXIT_FAILURE);
    }
    
    while(fscanf(bim_f_in, "%s %s %lf %ld %c %c", chrome_tmp, name, &distance, &position, &al1, &al2) == 6){
        new_n = (struct BIM_NODE *)malloc(sizeof(struct BIM_NODE));
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

        if (data_out){
            pt -> next = new_n;
            pt = new_n;

        } else {
            data_out = new_n;
            pt = new_n;
        }
    }

    return data_out;
}



static struct FAM_NODE *
parse_fam(char * fam_file_name)
{
    char family_id[10];
    char individual_id[10];
    char paternal_id[10];
    char maternal_id[10];
    int sex;
    int phenotype;
    struct FAM_NODE * data_out = NULL, * new_ = NULL, * pt = NULL;

    FILE * f_in = fopen(fam_file_name, "r");
    if (!f_in){
        fprintf(stderr, "error, fam file open failed.\n");
        exit(EXIT_FAILURE);

    }

    while(fscanf(f_in, "%s %s %s %s %d %d", family_id, individual_id, paternal_id, maternal_id, &sex, &phenotype) == 6){
        new_ = (struct FAM_NODE *)malloc(sizeof(struct FAM_NODE));
        strcpy(new_ -> family_id, family_id);
        strcpy(new_ -> individual_id, individual_id);
        strcpy(new_ -> paternal_id, paternal_id);
        strcpy(new_ -> maternal_id, maternal_id);
        new_ ->  sex = (char)sex;
        new_ -> phenotype = (signed char)phenotype;
        new_ -> next = NULL;

        if (data_out){
            pt -> next = new_;
            pt = new_;
        } else{
            data_out = pt = new_;
        }
    }
    return data_out;
}

/*
 * Bed file was generate by plink, it was is a binary file and used to record genotype data.
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
static struct BED_NODE *
parse_bed(char * bed_file_name, unsigned int raw_len)
{

    /*
        00 -> 0 -> 0
        10 -> 2 -> 1
        11 -> 3 -> 2
        01 -> 1 -> 3 missing

     */

    struct BED_NODE * data_out = NULL, * new_ = NULL, * pt = NULL;
    unsigned char * raw_buffer;
    int c;
    unsigned int i;
    unsigned char j;
    unsigned char gtype[4];

    FILE * f_in = fopen(bed_file_name, "r");
    if (!f_in){
        fprintf(stderr, "error, open bed file failed.\n");
        exit(EXIT_FAILURE);
    }


    //remove first 3 char, they are magic number which not belong to data.
    j = 0;
    while((c = fgetc(f_in)) != EOF){
        j++;
        if (j == 3)
            break;
    }
    

    i = raw_len;
    //printf("raw len %d", raw_len);
    while((c = fgetc(f_in)) != EOF){
        gtype[0] = (unsigned char)c & MASK1;
        gtype[1] = ((unsigned char)c & MASK2) >> 2;
        gtype[2] = ((unsigned char)c & MASK3) >> 4;
        gtype[3] = ((unsigned char)c & MASK4) >> 6;
        //printf("%d %d %d %d\n", gtype[0], gtype[1], gtype[2], gtype[3]);
        j = 0; 
        while(j < 4){
            if (i == raw_len){
                new_ = (struct BED_NODE *) malloc(sizeof(struct BED_NODE));
                new_ -> next = NULL;
                raw_buffer = (unsigned char *)malloc(raw_len);
                new_ -> genotype_raw = raw_buffer;
                if (data_out){
                    pt -> next = new_;
                    pt = new_;
                } else {
                    data_out = pt = new_;
                }
                i = 0;
            }
            raw_buffer[i] = gtype[j];
            i ++;
            //printf("%d\n", i);
            j ++;
        }

    }

    return data_out;
}




