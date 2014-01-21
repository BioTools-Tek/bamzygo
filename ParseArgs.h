#ifndef PARSEARGS_H
#define PARSEARGS_H

#include <QString>
#include <QStringList>
#include <iostream>

using namespace std;

class ParseArgs
{
public:
    ParseArgs(int argc, char *argv[]);
    //PARSE ARGS
    const char * colfile;
    const char * bamfile;
    const char * unkfile;

    bool extra;

    int range; //Default ~ Range of 1. If this was 10, then while the next line's position is less than 10 from the origin, it will include it in the same pileup call.

    int hetlim; // Default ~ 60
    int qualim; // Default ~ 10
//    int bad_depth;
//    bool outliers, headers, verbose, keep_chr;
};

void usage(){
    cerr << "v2013_07_16" << endl;
    cerr << "Takes in a two (or more) column file of <chr> <pos> for single positions (e.g. VCF, VCF-A, BED, etc) and appends the percentage of the max base, and the zygosity as two seperate columns at the end from an input bam file matching those positions" << endl;
    cerr << "\nusage: bam_indel_zygo <in.bam> <in.colfile> <out.unknowns> [OPTIONS]" << endl;
    cerr << endl;
    cerr << "   OPTIONS:"<<endl;
    cerr << "    --qualim N   Ignores bases with quality scores less than N" << endl;
    cerr << "    --hetlim N   Heterozygous if max base percent is less than N" << endl;
    cerr << "    --range N    Read multiple lines from colfile that fall within" << endl;
    cerr << "                 <last_position>+<N> per pileup call." << endl;
    cerr << "    --extra      include <bases> <quality> columns too in output" << endl;
    cerr << endl;
    exit(0);
}

ParseArgs::ParseArgs(int argc, char *argv[]){
    //PARSE ARGS
    bamfile="";colfile=""; unkfile="";

    //Default opts
    hetlim=60;qualim=0; range=1;
    extra=false;

    if (argc<4){
        usage();
    }
    if (argc >=4 ){
        bamfile = argv[1]; colfile = argv[2]; unkfile= argv[3];

        //Opts
        bool errs=false;
        QString errorsT = "Cannot parse: ";
        QStringList okay_vars;

        for (int i=4; i < argc; i++){
            QString tmp = argv[i];
            if(tmp=="--extra") extra=true;
            else if(tmp=="--hetlim"){
                if (argc <= i+1){
                    cerr << "Give a het limit value. Default is 60" << endl;
                    exit(-1);
                }
                hetlim = QString(argv[i+1]).toInt();
                okay_vars.append(argv[i+1]);
            }
            else if(tmp=="--range"){
                if (argc <= i+1){
                    cerr << "Give range value. Default is 1. (i.e. process line-by-line)" << endl;
                    exit(-1);
                }
                range = QString(argv[i+1]).toInt();
                okay_vars.append(argv[i+1]);
            }
            else if(tmp=="--qualim"){
                if (argc <= i+1){
                    cerr << "Give a quality limit value. Default is 10" << endl;
                    exit(-1);
                }
                qualim = QString(argv[i+1]).toInt();
                okay_vars.append(argv[i+1]);
            }

            else {
                bool error_yes=true;
                for (int i=0; i< okay_vars.length(); i++){
                    if (tmp==okay_vars.at(i)){
                        error_yes=false;
                        break;
                    }
                }
                if (error_yes) { errorsT.append(tmp+" "); errs=true;}
            }
        }
        if (errs) { cerr << errorsT.toUtf8().data() << endl; exit(-1);}
    }
}

#endif // PARSEARGS_H
