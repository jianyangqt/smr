#include <string>
#include <iostream>
#include <vector>

#include "SMR_main.hpp"
#include "../C/calmt.h"


using namespace std;


int
main(int argc, char * argv[])
{
    string usage = "Usage: smr [--version/-v] [--help/-h] <command> [<args>]\n"
                   "    \"smr command --help\" for command help message.";
    string version = "smr version 1.1.0";

    vector<string> help_mesg;
    string prog_name = "smr";
    string prog_description = 
        "Tool which implements the SMR & HEIDI methods to test for\n"
        "pleiotropic association between the expression level\n"
        "of a gene and a complex trait of interest using\n" 
        "summary-level data from GWAS and expression quantitative\n"
        "trait loci (eQTL) studies.";

    string mesg_smr_main = "main            smr main program.";
    string mesg_calmt =    "calmt           calculate cauchy between moulecular and trait.";
    string mesg_help =     "--help/-h       print help message and exit.";
    string mesg_version =  "--version/-v    print software version and exit.";

    help_mesg.push_back(version);
    help_mesg.push_back("\n");
    help_mesg.push_back(usage);
    help_mesg.push_back("\n");
    help_mesg.push_back("\n");
    help_mesg.push_back(prog_description);
    help_mesg.push_back("\n");
    help_mesg.push_back("\n");
    help_mesg.push_back(mesg_smr_main);
    help_mesg.push_back("\n");
    help_mesg.push_back(mesg_calmt);
    help_mesg.push_back("\n");
    help_mesg.push_back(mesg_help);
    help_mesg.push_back("\n");
    help_mesg.push_back(mesg_version);
    help_mesg.push_back("\n");


    if (argc == 1){
        cout << usage << endl;
        exit(0);
    }

    string first_arg = string(argv[1]);

    if (first_arg.substr(0, 1) == "-"){
        if (first_arg == "--help" || first_arg == "-h"){
            for (string mesg: help_mesg){
                cout << mesg;
            }
            exit(0);

        } else if(first_arg == "--version" || first_arg == "-h"){
            cout << version << endl;
            exit(0);

        } else if (first_arg == "--usage" || first_arg == "-u"){
            cout << usage << endl;
            exit(0);

        } else {
            cerr << "unknow option: " << first_arg << endl;
            cerr << usage << endl;
            exit(0);
        }

    } else {
        
        if (first_arg == "main"){
            smr_main(--argc, ++argv);
            exit(0);

        } else if (first_arg == "calmt"){
            calculate_cauchy(--argc, ++argv);
            exit(0);

        } else {
            cerr << first_arg << " is not smr command." << endl;
            cerr << usage << endl;
        }
    }

    return 0;
}

