#include <string>
#include <iostream>
#include <vector>

#include "SMR_main.hpp"
#include "calmt.hpp"


using namespace std;

string get_args(int, char *[]);


int
main(int argc, char * argv[])
{
    string usage = "Usage: smr [--version/-v] [--help/-h] <command> [<args>]";
    string version = "smr version 1.1.0";

    vector<string> help_mesg;
    string prog_name = "smr";
    string prog_description = 
        "Tool which implements the SMR & HEIDI methods to test for\n"
        "pleiotropic association between the expression level\n"
        "of a gene and a complex trait of interest using\n" 
        "summary-level data from GWAS and expression quantitative\n"
        "trait loci (eQTL) studies.";

    string mesg_help = "--help/-h:    print help message and exit.";
    string mesg_version = "--version/-v:     print software version and exit.";
    string mesg_smr_main = "smr_main:    smr main program.";
    string mesg_calmt = "calmt:    calculate cauchy between moulecular and trait.";

    help_mesg.push_back(version);
    help_mesg.push_back("\n");
    help_mesg.push_back(usage);
    help_mesg.push_back("\n");
    help_mesg.push_back("\n");
    help_mesg.push_back(prog_description);
    help_mesg.push_back("\n");
    help_mesg.push_back("\n");
    help_mesg.push_back(mesg_help);
    help_mesg.push_back("\n");
    help_mesg.push_back(mesg_version);
    help_mesg.push_back("\n");
    help_mesg.push_back(mesg_smr_main);
    help_mesg.push_back("\n");
    help_mesg.push_back(mesg_calmt);
    help_mesg.push_back("\n");


    string argument = "";
    argument = get_args(argc, argv);
    if (argc > 2){
        --argc;
        ++argv;
    }
    
    if (argument == ""){
        cout << usage << endl;
        exit(0);
    } else if (argument.substr(0, 1) == "-"){
        if (argument == "--help" || argument == "-h"){
            for (string mesg: help_mesg){
                cout << mesg;
            }
        } else if (argument == "--version" || argument == "-v"){
            cout << version << endl;
        } else {
            cout << "unknown option: " << argument << endl;
            cout << usage << endl;
        }
    } else {
        if (argument == "smr_main") {
            smr_main(argc, argv);
            exit(0);

        } else if (argument == "calmt"){
            calculate_cauchy(argc, argv);
            exit(0);

        } else {
            cout << argument << " is not a smr command" << endl;
            exit(0);
        }
    }

    
    return 0;
}


string
get_args(int argc, char * argv[])
{
    
    if (argc < 2){
        return "";
    } else {
        return argv[1];
    }
}


