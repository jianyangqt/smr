#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>

#include "my_functions.h"
#include "avl_tree.h"


#define MAX_POSITION_VALUE 0x7fffffff
#define MAX_PROBE_ID_LEN 63
#define MAX_GENE_ID_LEN 62

#define HASH_TABLE_LEN 65536
#define MASK_INT 2147483648

#define DENSE_FILE_TYPE_3 5
#define SPARSE_FILE_TYPE_3 3


struct EPI_DATA_LIST_NODE {
    //this field is required.
    unsigned char chromo;
    //this field is required.
    char probe_id[64];
    //can any string.
    char gene_distance[32];
    //max position value is 2**31 - 1, first left bit for NA. 0 indicate NA.
    uint32_t position;
    //can contain 62 character, most left for NA, 0 indicate NA.
    char gene_id[64];
    //0 for NA, '-' for -, '+' for +
    char orientation;
    struct EPI_DATA_LIST_NODE * next;

};


struct ESI_DATA_LIST_NODE{
    //required
    unsigned char chromo;
    //rerquired
    char snp_id[64];
    //any string
    char gene_distance[32];
    //max position value is 2**31 - 1, first left bit for NA. 0 indicate NA.
    uint32_t position;
    //can contain 30 characters, most lesf for NA, 0 indicate NA.
    char al1[32];
    char al2[32];
    //-1 for NA
    double frequence;
    struct ESI_DATA_LIST_NODE * next;

};


struct BESD_DATA{
    int file_type;
    long addn;
    unsigned long esi_num;
    unsigned long epi_num;
    struct EPI_DATA_LIST_NODE * epi_dt;
    struct ESI_DATA_LIST_NODE * esi_dt;
    struct BESD_PROBE_NODE {
        uint32_t probe_index;
        void * probe_infor_node;
        struct PROBE_SNP_NODE {
            uint32_t snp_index;
            void * snp_infor_node;
            float beta;
            float se;
            struct PROBE_SNP_NODE * next;

        } * snp_contained;

        struct BESD_PROBE_NODE * next;

    } * probe_node;

};


struct snp_index_btree{
        unsigned long index;
        uint64_t pointer_num;
        struct snp_index_btree * left;
        struct snp_index_btree * right;
};


static struct EPI_DATA_LIST_NODE * read_epi_file(const char *);
static struct ESI_DATA_LIST_NODE * read_esi_file(const char *);
static struct BESD_DATA read_besd_file_dense(const char *);
static struct BESD_DATA read_besd_file_sparse(const char *);
static struct BESD_DATA read_besd_file(const char *);
static void union_besd_data(struct EPI_DATA_LIST_NODE *, struct ESI_DATA_LIST_NODE *, struct BESD_DATA);
static struct BESD_DATA structure_besd_data(const char *);
static void add_besd_data(void **, struct BESD_DATA *);
static void ** merge_multi_besd_data(const char *);
static void reindex_probe(void **);
static struct snp_index_btree * make_btree(struct snp_index_btree *, uint64_t, unsigned long *);
static void count_treenode(struct snp_index_btree *, void *);
static void make_btree_hash(struct snp_index_btree *, void *);
static void look_up_btree(struct snp_index_btree *, void (*)(struct snp_index_btree *, void *), void *);
static struct tree_node * reindex_snp(void **);
static void print_epi_res(void **, const char *);
static void print_esi_res(struct snp_index_btree *, const char *);
static void print_besd_res(void **, const char *);
static void add_probe_node_of_slot(void **, unsigned short, struct BESD_PROBE_NODE *);
static void merge_besd_probe_node(struct BESD_PROBE_NODE *, struct BESD_PROBE_NODE *);


int
main(int argc, char * argv[])
{
/*
    struct EPI_DATA_LIST_NODE * epi_dt;
    epi_dt = read_epi_file(argv[1]);

    while(epi_dt){
        printf("%d\t%s\t%s\t%u\t%s\t%c\n", epi_dt -> chromo, epi_dt -> probe_id, epi_dt -> gene_distance, epi_dt -> position ^= MASK_INT, epi_dt -> gene_id, epi_dt -> orientation);
        epi_dt = epi_dt -> next;
    }
*/

/*
    struct ESI_DATA_LIST_NODE * esi_dt;
    esi_dt = read_esi_file(argv[1]);

    while(esi_dt){
        printf("%d\t%s\t%s\t%u\t%s\t%s\t%lf\n", esi_dt -> chromo, esi_dt -> snp_id, \
            esi_dt -> gene_distance, esi_dt -> position ^= MASK_INT, esi_dt -> al1, esi_dt -> al2, esi_dt -> frequence);
        esi_dt = esi_dt -> next;
    }
*/
/*
    struct BESD_DATA besd_dt;
    besd_dt = read_besd_file(argv[1]);
    printf("%d %ld %lu %lu\n", besd_dt.file_type, besd_dt.addn, besd_dt.esi_num, besd_dt.epi_num);
    struct BESD_PROBE_NODE * probe_dt;
    struct PROBE_SNP_NODE * snp_dt;

    probe_dt = besd_dt.probe_node;
    while(probe_dt){
        printf(">>%u\n", probe_dt ->probe_index);
        snp_dt = probe_dt -> snp_contained;
        while(snp_dt){
            printf("%u %f %f\n", snp_dt -> snp_index, snp_dt -> beta, snp_dt -> se);

            snp_dt = snp_dt -> next;
        }

        probe_dt = probe_dt -> next;

    }
*/
/*
    unsigned int i = 0;
    struct BESD_DATA besd_dt = structure_besd_data(argv[1]);
    void * hash_table[HASH_TABLE_LEN];
    for (i = 0; i < HASH_TABLE_LEN; i++){
        hash_table[i] = NULL;
    }

    add_besd_data(hash_table, besd_dt);
*/
    void ** hash_table;
    struct tree_node * index_tree;
    struct BESD_PROBE_NODE * probe_ptr;
    struct PROBE_SNP_NODE * snp_ptr;
    struct ESI_DATA_LIST_NODE * esi_node;

    hash_table = merge_multi_besd_data(argv[1]);
    reindex_probe(hash_table);

    index_tree = reindex_snp(hash_table);
    //traverse_tree_h(index_tree);

    print_epi_res(hash_table, "test.epi");
    //print_esi_res(index_btree, "test.esi");

    return 0;
}


static void
print_epi_res(void ** hash_table, const char * file_name)
{
    unsigned int i = 0;
    struct BESD_PROBE_NODE * probe_ptr;
    struct EPI_DATA_LIST_NODE * epi_ptr;
    char chromo_str[5];
    unsigned char chromo;
    char distance[32];
    char position_str[32];
    char gene_id[64];
    char orientation[32];


    FILE * f_out = fopen(file_name, "w");
    if (!f_out){
        fprintf(stderr, "open epi out file error\n");
        exit(1);
    }

    for (i = 0; i < HASH_TABLE_LEN; i++){
        probe_ptr = (struct BESD_PROBE_NODE *)hash_table[i];
        while(probe_ptr){
            epi_ptr = (struct EPI_DATA_LIST_NODE *)probe_ptr -> probe_infor_node;
            chromo = epi_ptr -> chromo;
            if (chromo == 23){
                strcpy(chromo_str, "X");
            } else if (chromo == 24){
                strcpy(chromo_str, "Y");
            } else{
                sprintf(chromo_str, "%d", chromo);
            }

            if ((epi_ptr -> gene_distance)[0] == 0){
                strcpy(distance, "NA");
            } else{
                strcpy(distance, epi_ptr -> gene_distance);
            }
            if ((epi_ptr -> position & MASK_INT) == 0){
                strcpy(position_str, "NA");
            } else{
                sprintf(position_str, "%d", epi_ptr -> position);
            }

            if ((epi_ptr -> gene_id)[0] == 0){
                strcpy(gene_id, "NA");
            } else{
                strcpy(gene_id, epi_ptr -> gene_id);
            }

            if (epi_ptr -> orientation == 0){
                strcpy(orientation, "NA");
            } else{
                orientation[0] = epi_ptr -> orientation;
            }

            fprintf(f_out, "%s\t%s\t%s\t%s\t%s\t%s\n", chromo_str, epi_ptr -> probe_id, distance, \
                position_str, gene_id, orientation);
            probe_ptr = probe_ptr -> next;
        }
    }
    fclose(f_out);

}


static void
print_esi_res(struct snp_index_btree * btree, const char * file_name)
{
    unsigned long tree_node_amount = 0;
    unsigned long i = 0;
    struct ESI_DATA_LIST_NODE * esi_infor;

    FILE * f_out = fopen(file_name, "w");
    if (!f_out){
        fprintf(stderr, "error, open esi out file failed.\n");
        exit(1);
    }

    look_up_btree(btree, count_treenode, &tree_node_amount);

/*
    uint64_t * tree_node_hash_list = (uint64_t *)malloc(sizeof(uint64_t) * tree_node_amount);
    look_up_btree(btree, make_btree_hash, tree_node_hash_list);
    for (i = 1; i < tree_node_amount; i++){
        esi_infor = (struct ESI_DATA_LIST_NODE *)tree_node_hash_list[i];
        fprintf(f_out, "%s", esi_infor -> snp_id);
    }
    free(tree_node_hash_list);

*/
}

static void
print_besd_res(void ** hash_table, const char * file_name)
{

}


static struct tree_node *
reindex_snp(void ** hash_table)
{
    struct tree_node ** tree_pt_pt = NULL;
    struct tree_node * tree_pt = NULL;
    struct tree_node * res_node = NULL;
    struct BESD_PROBE_NODE * probe_ptr = NULL;
    struct PROBE_SNP_NODE * snp_ptr = NULL;
    unsigned int i = 0;
    unsigned int j = 0;
    uint64_t snp_infor_pointer_num = 0;
    unsigned long node_index = 0;

    tree_pt_pt = &tree_pt;

    for (i = 0; i < HASH_TABLE_LEN; i++){
        probe_ptr = (struct BESD_PROBE_NODE *)hash_table[i];
        while (probe_ptr){
            j++;
            snp_ptr = probe_ptr -> snp_contained;
            while(snp_ptr){
                snp_infor_pointer_num = (uint64_t)snp_ptr -> snp_infor_node;
                //printf(">>%lu\n", snp_infor_pointer_num);
                res_node = make_avl_tree(tree_pt_pt, snp_infor_pointer_num, &node_index);
                snp_ptr -> snp_index = res_node -> index;
                snp_ptr = snp_ptr -> next;
            }
            probe_ptr = probe_ptr -> next;
        }
    }
    //printf(">>%lu\n", node_index);
    return tree_pt;
}


static void
look_up_btree(struct snp_index_btree * btree, void (* func)(struct snp_index_btree *, void *), void * container)
{
    struct snp_index_btree * left_node = NULL, * right_node = NULL;
    if (btree){
        left_node = btree -> left;
        right_node = btree -> right;
        func(btree, container);
        look_up_btree(left_node, func, container);
        look_up_btree(right_node, func, container);
    }
    return;
}


static void
count_treenode(struct snp_index_btree * btree, void * container)
{
    unsigned long * tmp;
    tmp = (unsigned long *)container;
    (*tmp)++;
    return;
}


static void
make_btree_hash(struct snp_index_btree * btree, void * container)
{
    uint64_t * tmp;
    tmp = (uint64_t *)container;
    tmp[btree -> index] = btree -> pointer_num;
    return;
}


static struct snp_index_btree *
make_btree(struct snp_index_btree * btree, uint64_t snp_pointer_num, unsigned long * node_index)
{
    struct snp_index_btree * keep = NULL;
    //printf("+++++++++++++++++++++++\n");
    while (btree){
        keep = btree;
        //printf("%lu\n", btree -> pointer_num);
        if (snp_pointer_num > (btree -> pointer_num)){
            btree = btree -> right;
        } else if (snp_pointer_num < (btree -> pointer_num)){
            btree = btree -> left;
        } else{
            printf("here\n");
            return btree;
        }
    }

    struct snp_index_btree * new_tree_node = (struct snp_index_btree *)malloc(sizeof(struct snp_index_btree));
    new_tree_node -> left = NULL;
    new_tree_node -> right = NULL;
    new_tree_node -> index = *node_index;
    (*node_index)++;
    new_tree_node -> pointer_num = snp_pointer_num;
    if (keep){
        if (snp_pointer_num > (keep -> pointer_num)){
            keep -> right = new_tree_node;
        } else{
            keep -> left = new_tree_node;
        }
    }
    return new_tree_node;
}


static void
reindex_probe(void ** hash_table)
{
    struct BESD_PROBE_NODE * probe_node_ptr = NULL;
    unsigned int i = 0;
    unsigned long j = 0;

    j = 0;
    for (i = 0; i < HASH_TABLE_LEN; i++){
        //printf(">>%u\n", i);
        probe_node_ptr = (struct BESD_PROBE_NODE *)hash_table[i];
        while(probe_node_ptr){
                probe_node_ptr -> probe_index = j;
                j++;
                probe_node_ptr = probe_node_ptr -> next;
        }
    }
    return;
}


static void **
merge_multi_besd_data(const char * besd_file_name_list)
{
    FILE * f_in = fopen(besd_file_name_list, "r");
    char besd_file_name[1024];
    struct BESD_DATA besd_dt;
    unsigned int i = 0;
    void ** hash_table = (void **)malloc(sizeof(void *) * HASH_TABLE_LEN);
    for (i = 0; i < HASH_TABLE_LEN; i++){
        hash_table[i] = NULL;
    }
    while(fscanf(f_in, "%s", besd_file_name) == 1){
        printf("%s\n", besd_file_name);
        besd_dt = structure_besd_data(besd_file_name);
        add_besd_data(hash_table, &besd_dt);
    }

    return hash_table;
}


static struct BESD_DATA
structure_besd_data(const char * file_name_prefix)
{
    char epi_file_name[128];
    char esi_file_name[128];
    char besd_file_name[128];

    if(strlen(file_name_prefix) > 123){
        fprintf(stderr, "error, besd file name prefix is too long.\n");
        exit(1);
    }

    strcpy(epi_file_name, file_name_prefix);
    strcpy(esi_file_name, file_name_prefix);
    strcpy(besd_file_name, file_name_prefix);

    strcat(epi_file_name, ".epi");
    strcat(esi_file_name, ".esi");
    strcat(besd_file_name, ".besd");

    struct EPI_DATA_LIST_NODE * epi_dt;
    struct ESI_DATA_LIST_NODE * esi_dt;
    struct BESD_DATA besd_dt;

    epi_dt = read_epi_file(epi_file_name);
    esi_dt = read_esi_file(esi_file_name);
    besd_dt = read_besd_file(besd_file_name);

/*
    struct BESD_PROBE_NODE * tmp;
    tmp = besd_dt.probe_node;
    while(tmp){
        printf("%u\n", tmp -> probe_index);
        tmp = tmp -> next;
    }
*/

    union_besd_data(epi_dt, esi_dt, besd_dt);

    return besd_dt;
}


static struct EPI_DATA_LIST_NODE *
read_epi_file(const char * epi_file_name)
{
    char chromo_str[5];
    unsigned char chromo;
    char probe_id[64];
    char gene_distance[32];
    char position_str[32];
    uint32_t position;
    char gene_id[64];
    char orientation_str[5];
    char orientation;
    struct EPI_DATA_LIST_NODE * first_node = NULL, * new_node = NULL, * ptr = NULL;


    FILE * f_in = fopen(epi_file_name, "r");
    if (!f_in){
        fprintf(stderr, "epi file open failed.\n");
        exit(1);

    }

    while(fscanf(f_in, "%s %s %s %s %s %s", chromo_str, probe_id, \
        gene_distance, position_str, gene_id, orientation_str) == 6){
        if (strcmp(chromo_str, "X") == 0 || strcmp(chromo_str, "x") == 0){
            chromo = 23;
        } else if (strcmp(chromo_str, "Y") == 0 || strcmp(chromo_str, "y") == 0){
            chromo = 24;
        } else if (isdigit(chromo_str[0])){
            chromo = (unsigned char)atoi(chromo_str);
        }

        if (strlen(probe_id) > 63){
            probe_id[63] = '\0';
            fprintf(stderr, "Warning, probe name too long. cut to 63 char\n, ");
        }

        if (strcmp(position_str, "NA") == 0){
            position = 0;
        } else if (isdigit(position_str[0])){
            position = atoi(position_str);
            if (position > MAX_POSITION_VALUE){
                fprintf(stderr, "error, the probe position to big.\n");
                exit(1);
            }
            position |= MASK_INT;

        } else{
            fprintf(stderr, "probe position should be integer.\n");
            exit(1);
        }

        if (strlen(gene_id) > 63){
            gene_id[63] = '\0';
            fprintf(stderr, "Warning, gene id too long, and trim to 63 char.\n");
        }

        if (strcmp(gene_id, "NA") == 0){
            gene_id[0] = 0;

        }

        if (strcmp(orientation_str, "NA") == 0){
            orientation  = 0;

        } else if (strlen(orientation_str) == 1 && orientation_str[0] == '+'){
            orientation = '+';

        } else if (strlen(orientation_str) == 1 && orientation_str[0] == '-'){
            orientation  = '-';

        } else{
            fprintf(stderr, "orientation not recognized, set to NA");
            orientation  = 0;

        }

        new_node = (struct EPI_DATA_LIST_NODE *)malloc(sizeof(struct EPI_DATA_LIST_NODE));
        new_node -> next = NULL;
        new_node -> chromo = chromo;
        strcpy(new_node -> probe_id, probe_id);
        strcpy(new_node -> gene_distance, gene_distance);
        new_node -> position = position;
        strcpy(new_node -> gene_id, gene_id);
        new_node -> orientation = orientation;

        if(first_node){
            ptr -> next = new_node;
            ptr = new_node;

        } else{
            first_node = ptr = new_node;

        }
    }

    if (!feof(f_in)){
        fprintf(stderr, "file not completely readed.\n");
        fclose(f_in);
        exit(1);
    }
    fclose(f_in);

    return first_node;

}


static struct ESI_DATA_LIST_NODE *
read_esi_file(const char * esi_file_name)
{
    FILE * f_in = fopen(esi_file_name, "r");
    if (!f_in){
        fprintf(stderr, "esi open failed\n");
        exit(1);
    }

    char chromo_str[5];
    unsigned char chromo;
    char snp_id[64];
    char gene_distance[32];
    char position_str[32];
    uint32_t position;
    char al1[516];
    char al2[516];
    char frequence_str[32];
    double frequence;
    unsigned long i = 0;

    struct ESI_DATA_LIST_NODE * first_node = NULL, * new_node = NULL, * ptr = NULL;

    while(fscanf(f_in, "%s %s %s %s %s %s %s", chromo_str, snp_id, gene_distance, \
        position_str, al1, al2, frequence_str) == 7){
            //fprintf(stderr, ">>%lu  %s %s %s %s %s %s %s\n", i, chromo_str, snp_id, gene_distance, position_str, al1, al2, frequence_str);
            i++;
            if (strcmp(chromo_str, "X") == 0 || strcmp(chromo_str, "x") == 0){
                chromo = 23;
            } else if (strcmp(chromo_str, "Y") == 0 || strcmp(chromo_str, "y") == 0){
                chromo = 24;
            } else if(isdigit(chromo_str[0])){
                chromo = (unsigned char)atoi(chromo_str);
            } else{
                fprintf(stderr, "chromosome is not recognized in esi file\n");
                exit(1);
            }

            if (strlen(snp_id) > 63){
                snp_id[63] = '\0';
                fprintf(stderr, "Warning, snp id is too long, and is trimed into \
                    63 char.\n");
            }

            if (strlen(gene_distance) > 31){
                fprintf(stderr, "gene distance is trimed. %s\n", gene_distance);
                gene_distance[31] = '\0';
            }
            if (strcmp(gene_distance, "NA") == 0){
                gene_distance[0] = 0;

            }

            if (strcmp(position_str, "NA") == 0){
                position = 0;
            } else if(isdigit(position_str[0])){
                position = atoi(position_str);
                if (position > MAX_POSITION_VALUE){
                    fprintf(stderr, "snp position exect max value.\n");
                    exit(1);
                }
                position |= MASK_INT;
            }

            if (strlen(al1) > 31){
                fprintf(stderr, "al1 is trimed. %s\n", al1);
                al1[31] = '\0';
            }
            if (strlen(al2) > 31){
                fprintf(stderr, "al2 is trimed. %s\n", al2);
                al2[31] = '\0';
            }

            if(strcmp(frequence_str, "NA") == 0){
                frequence = -1;
            } else{
                frequence = atof(frequence_str);
            }

            new_node = (struct ESI_DATA_LIST_NODE *)malloc(sizeof(struct ESI_DATA_LIST_NODE));
            new_node -> next = NULL;
            new_node -> chromo = chromo;
            strcpy(new_node -> snp_id, snp_id);
            strcpy(new_node -> gene_distance, gene_distance);
            new_node -> position = position;
            strcpy(new_node -> al1, al1);
            strcpy(new_node -> al2, al2);
            new_node -> frequence = frequence;

            if(first_node){
                ptr -> next = new_node;
                ptr = new_node;
            } else{
                first_node = ptr = new_node;

            }
        }

        if(!feof(f_in)){
            fprintf(stderr, "esi file not completely readed.\n");
            fclose(f_in);
            exit(1);
        }
        fclose(f_in);

        return first_node;
}


static struct BESD_DATA
read_besd_file(const char * besd_file_name)
{
    FILE * f_in = fopen(besd_file_name, "r");
    if (!f_in){
        fprintf(stderr, "open besd file failed.\n");
        exit(1);
    }

    int file_type;
    if (fread(&file_type, sizeof(int), 1, f_in) != 1){
        fprintf(stderr, "read besd file type failed\n");
        exit(1);
    }
    fclose(f_in);

    struct BESD_DATA data_out = {0, 0, 0, 0, 0};
    if (file_type == DENSE_FILE_TYPE_3){
        data_out = read_besd_file_dense(besd_file_name);
    } else if (file_type == SPARSE_FILE_TYPE_3){
        data_out = read_besd_file_sparse(besd_file_name);
    } else{
        fprintf(stderr, "besd file type not recognized.\n");
        exit(1);
    }

    return data_out;
}


struct BESD_DATA
read_besd_file_dense(const char * besd_file_name)
{
    struct BESD_DATA data_out;
    data_out.probe_node = NULL;
    FILE * f_in = fopen(besd_file_name, "r");
    int first_16_int[16];
    int file_type, addn, esi_num, epi_num;
    uint32_t i, j;

    if (fread(first_16_int, sizeof(int), 16, f_in) == 16){
        file_type = first_16_int[0];
        addn = first_16_int[1];
        esi_num = first_16_int[2];
        epi_num = first_16_int[3];
    } else{
        fprintf(stderr, "read besd file first 16 int failed.\n");
        exit(1);
    }
    data_out.file_type = file_type;
    data_out.addn = addn;
    data_out.esi_num = esi_num;
    data_out.epi_num = epi_num;

    struct BESD_PROBE_NODE * first_node = NULL, * new_node = NULL, * ptr = NULL;
    struct PROBE_SNP_NODE * head_node = NULL, * birth_node = NULL, * ptr_esi = NULL;
    //float block_buffer[esi_num * 2];
    float * block_buffer = (float *)malloc(sizeof(float) * esi_num * 2);
    const uint32_t block_len = esi_num * 2;
    float beta, se;

    for (i = 0; i < epi_num; i++){
        head_node = NULL;
        birth_node = NULL;
        ptr_esi = NULL;
        new_node = (struct BESD_PROBE_NODE *)malloc(sizeof(struct BESD_PROBE_NODE));
        new_node -> next = NULL;
        new_node -> snp_contained = NULL;
        new_node -> probe_index = i;
        if (fread(block_buffer, sizeof(float), block_len, f_in) == block_len){
            for (j = 0; j < esi_num; j++){
                beta = block_buffer[j];
                se = block_buffer[esi_num + j];

                if (beta != -9 && se != -9){
                    birth_node = (struct PROBE_SNP_NODE *)malloc(sizeof(struct PROBE_SNP_NODE));
                    birth_node -> next = NULL;
                    birth_node -> snp_index = j;
                    birth_node -> beta = beta;
                    birth_node -> se = se;
                    if (head_node){
                        ptr_esi -> next = birth_node;
                        ptr_esi = birth_node;
                    } else{
                        head_node = ptr_esi = birth_node;
                    }
                } else if (beta == -9 && se == -9){
                    ;
                } else{
                    fprintf(stderr, "besd file error, beta and se value should \
                        both be -9 or a both a none -9 value.\n");
                    exit(1);
                }
            }
        } else{
            fprintf(stderr, "read besd block data failed.\n");
            exit(1);
        }
        new_node -> snp_contained = head_node;
        if (first_node){
            ptr -> next = new_node;
            ptr = new_node;
        } else{
            first_node = ptr = new_node;
        }
    }
    fgetc(f_in); //remove EOF char
    if(!feof(f_in)){
        fprintf(stderr, "file not completely readed.\n");
        exit(1);
    }
    fclose(f_in);
    data_out.probe_node = first_node;

    return data_out;
}


struct BESD_DATA
read_besd_file_sparse(const char * besd_file_name)
{
    //printf("spase\n");
    struct BESD_DATA data_out;
    data_out.probe_node = NULL;
    FILE * f_in = fopen(besd_file_name, "r");
    int first_16_int[16];
    int file_type, addn, esi_num, epi_num;
    unsigned long i, j;
    if (fread(first_16_int, sizeof(int), 16, f_in) == 16){
        file_type = first_16_int[0];
        addn = first_16_int[1];
        esi_num = first_16_int[2];
        epi_num = first_16_int[3];
    } else{
        fprintf(stderr, "error, read besd file first 16 int failed\n");
        exit(1);
    }
    data_out.file_type = file_type;
    data_out.addn = addn;
    data_out.esi_num = esi_num;
    data_out.epi_num = epi_num;

    struct BESD_PROBE_NODE * first_node = NULL, * new_node = NULL, * ptr = NULL;
    struct PROBE_SNP_NODE * head_node = NULL, * birth_node = NULL, * ptr_esi = NULL;
    uint64_t val_num;
    unsigned long offset_len = epi_num * 2 + 1;
    uint64_t probe_offset[offset_len];

    if (fread(&val_num, sizeof(uint64_t), 1, f_in) != 1){
        fprintf(stderr, "error, read value number bytes error\n");
        exit(1);
    }

    if (fread(probe_offset, sizeof(uint64_t), offset_len, f_in) != offset_len){
        fprintf(stderr, "error, read probe value offset error.\n");
        exit(1);
    }

    unsigned long max_value_block_len = 0;
    unsigned long block_len = 0;
    for (i = 0; i < offset_len - 1; i += 2){
        //printf("%lu %lu\n", probe_offset[i], probe_offset[i + 1]);
        block_len = probe_offset[i + 1] - probe_offset[i];
        if (block_len > max_value_block_len)
            max_value_block_len = block_len;
    }
    //printf("%lu\n", max_value_block_len);

    unsigned long max_value_block_len_x2 = max_value_block_len * 2;
    //uint32_t snp_index_block_buffer[max_value_block_len_x2];
    uint32_t * snp_index_block_buffer = (uint32_t *)malloc(sizeof(uint32_t) * max_value_block_len_x2);
    //float value_block_buffer[max_value_block_len_x2];
    float * value_block_buffer = (float *)malloc(sizeof(float) * max_value_block_len_x2);
    fpos_t cur_pos;

    for (i = 0; i < offset_len - 1; i += 2){
        new_node = (struct BESD_PROBE_NODE *)malloc(sizeof(struct BESD_PROBE_NODE));
        new_node -> next = NULL;
        new_node -> snp_contained = NULL;
        head_node = NULL;
        birth_node = NULL;
        ptr_esi = NULL;
        new_node -> probe_index = i / 2;

        block_len = probe_offset[i + 2] - probe_offset[i];
        if (fread(snp_index_block_buffer, sizeof(uint32_t), block_len, f_in) != block_len){
            fprintf(stderr, "read erroe\n");
            exit(1);
        }
        fgetpos(f_in, &cur_pos);
        fseek(f_in, val_num * sizeof(uint32_t) - block_len * sizeof(uint32_t), SEEK_CUR);
        if (fread(value_block_buffer, sizeof(float), block_len, f_in) != block_len){
            fprintf(stderr, "read error\n");
            exit(1);
        }
        fsetpos(f_in, &cur_pos);
        block_len = block_len / 2;
        for (j = 0; j < block_len; j++){
            birth_node = (struct PROBE_SNP_NODE *)malloc(sizeof(struct PROBE_SNP_NODE));
            birth_node -> next = NULL;
            birth_node -> snp_index = snp_index_block_buffer[j];
            birth_node -> beta = value_block_buffer[j];
            birth_node -> se = value_block_buffer[j + block_len];

            if (head_node){
                ptr_esi -> next = birth_node;
                ptr_esi = birth_node;
            } else{
                head_node = ptr_esi = birth_node;
            }
        }

        new_node -> snp_contained = head_node;

        if (first_node){
            ptr -> next = new_node;
            ptr = new_node;
        } else{
            first_node = ptr = new_node;
        }
    }

    fseek(f_in, val_num * sizeof(uint32_t), SEEK_CUR);
    fgetc(f_in);
    if (!feof(f_in)){
        fprintf(stderr, "error, file not read completely\n");
        exit(1);
    }
    fclose(f_in);
    data_out.probe_node = first_node;

    return data_out;
}


static void
union_besd_data(struct EPI_DATA_LIST_NODE * epi_dt, \
    struct ESI_DATA_LIST_NODE * esi_dt, struct BESD_DATA besd_dt)
{
    struct BESD_PROBE_NODE * besd_probe_ptr = NULL;
    struct PROBE_SNP_NODE * probe_snp_ptr = NULL;

    unsigned long epi_len = 0;
    unsigned long esi_len = 0;
    unsigned long i = 0;

    struct EPI_DATA_LIST_NODE * epi_node_ptr = NULL;
    struct ESI_DATA_LIST_NODE * esi_node_ptr = NULL;

    epi_node_ptr = epi_dt;
    epi_len = 0;
    while(epi_node_ptr){
        epi_len ++;
        epi_node_ptr = epi_node_ptr -> next;
    }

    esi_node_ptr = esi_dt;
    esi_len = 0;
    while(esi_node_ptr){
        esi_len ++;
        esi_node_ptr = esi_node_ptr -> next;
    }

    struct EPI_DATA_LIST_NODE ** epi_index_table = (struct EPI_DATA_LIST_NODE **)\
        malloc(sizeof(struct EPI_DATA_LIST_NODE *) * epi_len);
    char * epi_refed = (char *)malloc(sizeof(char) * epi_len);
    epi_node_ptr = epi_dt;
    for (i = 0; i < epi_len; i++){
        epi_refed[i] = 0;
        epi_index_table[i] = epi_node_ptr;
        epi_node_ptr = epi_node_ptr -> next;
    }

    struct ESI_DATA_LIST_NODE ** esi_index_table = (struct ESI_DATA_LIST_NODE **)\
        malloc(sizeof(struct ESI_DATA_LIST_NODE *) * esi_len);
    char * esi_refed = (char *)malloc(sizeof(char) * esi_len);
    esi_node_ptr = esi_dt;
    for (i = 0; i < esi_len; i++){
        esi_refed[i] = 0;
        esi_index_table[i] = esi_node_ptr;
        esi_node_ptr = esi_node_ptr -> next;
    }

    besd_probe_ptr = besd_dt.probe_node;
    while(besd_probe_ptr){
        besd_probe_ptr -> probe_infor_node = (void *)epi_index_table[besd_probe_ptr -> probe_index];
        epi_refed[besd_probe_ptr -> probe_index] = 1;
        probe_snp_ptr = besd_probe_ptr -> snp_contained;
        while(probe_snp_ptr){
            probe_snp_ptr -> snp_infor_node = (void *)esi_index_table[probe_snp_ptr -> snp_index];
            esi_refed[probe_snp_ptr -> snp_index] = 1;
            probe_snp_ptr = probe_snp_ptr -> next;
        }
        besd_probe_ptr = besd_probe_ptr -> next;
    }

    //besd_dt.epi_dt = epi_dt;
    for (int i = 0; i < epi_len; i++){
        if (epi_refed[i] == 1){
            epi_index_table[i] -> next = NULL;
        } else{
            free(epi_index_table[i]);
        }
    }

    //besd_dt.esi_dt = esi_dt;
    for (int i = 0; i < esi_len; i++){
        if (esi_refed[i] == 1){
            esi_index_table[i] -> next = NULL;
        } else{
            free(esi_index_table[i]);
        }
    }

    free(epi_index_table);
    free(epi_refed);
    free(esi_index_table);
    free(esi_refed);
    return;
}


static void
add_besd_data(void ** hash_table, struct BESD_DATA * besd_dt)
{
    struct BESD_PROBE_NODE * besd_probe_ptr, * probe_ptr;
    struct EPI_DATA_LIST_NODE * epi_node;
    unsigned short probe_hash_index = 0;
    char * probe_id;

    besd_probe_ptr = besd_dt -> probe_node;
    besd_dt -> probe_node = NULL;
    while (besd_probe_ptr){
        probe_ptr = besd_probe_ptr -> next;
        besd_probe_ptr -> next = NULL;
        epi_node = (struct EPI_DATA_LIST_NODE *)(besd_probe_ptr -> probe_infor_node);
        probe_id = epi_node -> probe_id;
        probe_hash_index = hash_func(probe_id);
        //printf("hash_index:%s %hu\n", probe_id, probe_hash_index);
        add_probe_node_of_slot(hash_table, probe_hash_index, besd_probe_ptr);
        besd_probe_ptr = probe_ptr;
        //exit(0);
    }
    return;
}


static void
add_probe_node_of_slot(void ** hash_table, unsigned short probe_hash_index, struct BESD_PROBE_NODE * besd_probe_ptr)
{
    char * probe_id_a, * probe_id_b;
    uint32_t probe_position_a, probe_position_b;
    unsigned char probe_chromo_a, probe_chromo_b;
    //unsigned char linked_marke = 0;
    struct BESD_PROBE_NODE * probe_ptr, * probe_ptr_keep;
    struct EPI_DATA_LIST_NODE * epi_node;
    epi_node = besd_probe_ptr -> probe_infor_node;
    probe_id_a = epi_node -> probe_id;
    probe_position_a = epi_node -> position;
    probe_chromo_a = epi_node -> chromo;

    if (hash_table[probe_hash_index] == NULL){
        hash_table[probe_hash_index] = (void *)besd_probe_ptr;
        //printf("new\n");
        return;
    } else{
        probe_ptr = (struct BESD_PROBE_NODE *)hash_table[probe_hash_index];
        //printf("slot merged\n");
        while(probe_ptr){
            probe_ptr_keep = probe_ptr;
            epi_node = probe_ptr -> probe_infor_node;
            probe_id_b = epi_node -> probe_id;
            probe_position_b = epi_node -> position;
            probe_chromo_b = epi_node -> chromo;

            if (strcmp(probe_id_a, probe_id_b) == 0){
                if ((probe_chromo_a == probe_chromo_b) && (probe_position_a == probe_position_b || probe_position_a == 0 || probe_position_b == 0)){
                    besd_probe_ptr -> probe_infor_node = NULL;
                    free(epi_node);
                    merge_besd_probe_node(besd_probe_ptr, probe_ptr);
                    free(besd_probe_ptr);
                    return;
                } else{
                    fprintf(stderr, "error, same probe id shoude located at same \
                        chromosome, and have same position. This probe will be ignored.\n");
                    return;
                }
            }
            probe_ptr = probe_ptr -> next;
        }
        probe_ptr_keep -> next = besd_probe_ptr;
    }

    return;
}


// merge information of node a into node b.
static void
merge_besd_probe_node(struct BESD_PROBE_NODE * node_a, struct BESD_PROBE_NODE * node_b)
{
    struct PROBE_SNP_NODE * snp_node_a, * snp_node_a_tmp, * snp_node_b, * snp_node_b_keep;
    struct ESI_DATA_LIST_NODE * esi_node_a, * esi_node_b;
    char * snp_id_a, * snp_id_b;
    uint32_t snp_positon_a, snp_positon_b;
    unsigned char snp_chromo_a, snp_chromo_b;
    unsigned char merged_marke = 0;

    snp_node_a = node_a -> snp_contained;
    node_a -> snp_contained = NULL;
    snp_node_b = node_b -> snp_contained;
    snp_node_b_keep = snp_node_b;
    //here is really a time consumming step.
    while(snp_node_a){
        snp_node_a_tmp = snp_node_a -> next;
        snp_node_a -> next = NULL;
        esi_node_a = (struct ESI_DATA_LIST_NODE *)snp_node_a -> snp_infor_node;
        snp_id_a = esi_node_a -> snp_id;
        snp_chromo_a = esi_node_a -> chromo;
        snp_positon_a = esi_node_a -> position;
        merged_marke = 0;
        while(snp_node_b){
            esi_node_b = (struct ESI_DATA_LIST_NODE *)snp_node_b -> snp_infor_node;
            snp_id_b = esi_node_b -> snp_id;
            snp_chromo_b = esi_node_b -> chromo;
            snp_positon_b = esi_node_b -> position;
            if (strcmp(snp_id_a, snp_id_b) == 0){
                if ((snp_chromo_a == snp_chromo_b) && (snp_positon_a == snp_positon_b || snp_positon_a == 0 || snp_positon_b == 0)){
                    merged_marke = 1;
                    snp_node_a -> snp_infor_node = NULL;
                    free(esi_node_a);
                    free(snp_node_a);
                    break;
                } else{
                    fprintf(stderr, "snp have same id shoude located in same \
                        chromosome, and same position. The snp %s will be obmit.\n", snp_id_a);
                    break;
                }
            }

            if (snp_node_b -> next == NULL && !merged_marke){
                snp_node_b -> next = snp_node_a;
                break;
            }
            snp_node_b = snp_node_b -> next;
        }

        snp_node_b = snp_node_b_keep;
        snp_node_a = snp_node_a_tmp;
    }

    return;
}
