#include <stdio.h>
#include <argp.h>
#include <stdlib.h>


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
    struct arguments args = get_args(argc, argv);

    return 0;
}


