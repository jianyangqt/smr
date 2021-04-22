#include <stdio.h>
#include <argp.h>
#include <stdlib.h>
#include <stdbool.h>


#include "file_parser.h"

#define TEST_FUNC

struct arguments {
    const char * bfile;
    const char * gwas_summary;
    const char * beqtl_summary;
    const char * out;
    const char * thread;

};


static struct argp_option options[] = {
    {"bfile", 'b', "FILE", 0, "bfile prefix", 0},
    {"gwas_summary", 'g', "FILE", 0, "gwas summary file name", 1},
    {"beqtl_summary", 'e', "FILE", 0, "beqtl summary file prefix", 2},
    {"out", 'o', "FILE", 0, "out file name", 3},
    {"thread", 't', "INT", 0, "thread number using", 4},
    {0}
};

static error_t
opt_parser(int key, char * arg, struct argp_state * state)
{
    struct arguments * arguments =  (struct arguments *)state -> input;
    switch (key){
    case 'b':
        arguments -> bfile = arg;
        break;
    case 'g':
        arguments -> gwas_summary = arg;
        break;
    case 'e':
        arguments -> beqtl_summary = arg;
        break;
    case 'o':
        arguments -> out = arg;
        break;
    case 't':
        arguments -> thread = arg;
        break;
    case ARGP_KEY_ARG:
        argp_usage(state);
        break;
    case ARGP_KEY_END:
        if (state -> arg_num > 0)
            argp_usage(state);
        break;
    default:
        return ARGP_ERR_UNKNOWN;
        break;
    }

    return 0;
}


static char arg_doc[] = "<no position arguments>";
static char doc[] = "calculate molecular and trait relation.";
static struct argp argp = {options, opt_parser, arg_doc, doc};
static void filter_msmr_eqtl(struct MSMR_DATA_eqtl *, double, double);


struct arguments
get_args(int argc, char * argv[])
{
    struct arguments args;
    argp_parse(&argp, argc, argv, 0, 0, &args);
    return args;
}


int
calculate_cauchy(int argc, char * argv[])
{
    if (argc < 2){
        printf("smr calmt --help for more information\n");
        exit(1);
    }
    //struct arguments args = get_args(argc, argv);
    
    const char * k4me1_eqtl = argv[1];
    const char * h27ac_eqtl = argv[2];
    const char * caqtl_eqtl = argv[3];
    const char * mqtl_eqtl = argv[4];

    struct MSMR_DATA_eqtl k4me1_eqtl_dt = parse_msmr_file_eqtl(k4me1_eqtl, true);
    struct MSMR_DATA_eqtl h27ac_eqtl_dt = parse_msmr_file_eqtl(h27ac_eqtl, true);
    struct MSMR_DATA_eqtl caqtl_eqtl_dt = parse_msmr_file_eqtl(caqtl_eqtl, true);
    struct MSMR_DATA_eqtl mqtl_eqtl_dt = parse_msmr_file_eqtl(mqtl_eqtl, true);

    double max_p_smr_multi = 0;
    double min_p_heidi = 0;
    max_p_smr_multi = 0.05 / k4me1_eqtl_dt.node_len;
    min_p_heidi = 0.01;
    filter_msmr_eqtl(&k4me1_eqtl_dt, max_p_smr_multi, min_p_heidi);
    max_p_smr_multi = 0.05 / h27ac_eqtl_dt.node_len;
    filter_msmr_eqtl(&h27ac_eqtl_dt, max_p_smr_multi, min_p_heidi);
    max_p_smr_multi = 0.05 / caqtl_eqtl_dt.node_len;
    filter_msmr_eqtl(&caqtl_eqtl_dt, max_p_smr_multi, min_p_heidi);
    max_p_smr_multi = 0.05 / mqtl_eqtl_dt.node_len;
    filter_msmr_eqtl(&mqtl_eqtl_dt, max_p_smr_multi, min_p_heidi);


    return 0;
}



#ifdef TEST_FUNC
int
main(int argc, char * argv[])
{

    calculate_cauchy(argc, argv);


    return 0;
}
#endif


static void
filter_msmr_eqtl(struct MSMR_DATA_eqtl * eqtl_dt, double max_p_smr_multi, double min_p_heidi)
{
    struct MSMR_NODE_eqtl * pt = NULL, * tmp = NULL;
    pt = eqtl_dt -> msmr_node;

    while(pt){
        if (pt -> p_smr_multi >= max_p_smr_multi || pt -> p_heidi <= min_p_heidi){
            eqtl_dt -> node_len -= 1;
            tmp = pt;
            if (pt -> prev){
                (pt -> prev) -> next = pt -> next;
            } 

            if (pt -> next){
                (pt -> next) -> prev = pt -> prev;
            }
        }
        pt = pt -> next;
        free(tmp);
    }
    //here is a danger of geneate dangling pointer
    if (eqtl_dt -> node_len == 0){
        eqtl_dt -> msmr_node = NULL;
    }
    return;
}





