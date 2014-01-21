#include "Mpileup.h"
#include "ParseArgs.h"

#include <QFile>
#include <QTextStream>

//#define DEBUG
//#define BUFFERLINE_APPROACH

#ifdef DEBUG
#define DEBUG_BASES
#define DEBUG_SCORES
#endif

/**
  Takes in a two (or more) column file of <chr> <pos> for single positions (e.g. VCF, VCF-A, BED, etc) and appends
  the percentage of the max base, and the zygosity as two seperate columns at the end from an input bam file matching
  those positions
**/


uint countlines(const char * file){
    uint lines=0;
    //For lines in infile
    QFile inputFile(file);
    if (inputFile.open(QIODevice::ReadOnly))
    {
        while (!inputFile.atEnd()){
            inputFile.readLine();
            lines ++;
        }
    }
    return lines;
}

ofstream unk;


int main(int argc, char *argv[])
{
    ParseArgs *pa = new ParseArgs(argc, argv);

    unk.open(pa->unkfile); //Open unknowns

    quint64 maxval = -1;
   // cerr <<  "maxval=~" << maxval << endl;
    //exit(-1);

    //PSEUDO:
    //For lines in infile:
    // 1. Read and store N lines
    // 2. Prepare N lines for single MPileUp Call --> Pileup Output
    // 3. Match pileup output with the input positions:
    // a. returns <chr> <pos> <bases> <quality> array if successful
    //    else cerr << "<chr> <pos> NOT IN BAM" << endl. If all empty, return []

    // 4. Read stored lines of infile and for each matching <chr> <pos>:
    //  a. parse bases and quality --> determine maxbase, zygosity.
    //  b. append: <bases> <quality> <base> <zygo> (if --extra flag) else <base> <zygo>
    //  c. Print line to cout

    //Find num lines
    uint num_lines = countlines(pa->colfile);
#ifdef DEBUG
    cerr << "Num lines:" << num_lines << endl;
#endif

    //For lines in infile
    QFile inputFile(pa->colfile);

    if (inputFile.open(QIODevice::ReadOnly))
    {
        QTextStream in(&inputFile);
        uint count_lines = 0;
        bool double_stop = false;

        QString line="##";

        handleHeaders(line, in, count_lines, unk, pa->extra);

        bool alreadychecked = false;

        while ( !double_stop )
        {
            QStringList lines;
            QList<PileUp*> line_data;
            QString first_chrome="";
            uint first_index=0;

            bool stop = false;

#ifdef BUFFERLINE_APPROACH
            QString bufferline = "";   //Holds last position before line becomes out of range
            QStringList available_lines;
#else
            qint64 bufferpos=0;
#endif

            while( !stop )
            {
                if ( in.atEnd() ) {double_stop=true; break;}
                QString line = in.readLine();

#ifdef BUFFERLINE_APPROACH
                available_lines.append(line);

                while (available_lines.length() > 0){
                    QString line = available_lines.first();
#endif

                    if (line.length()>0)
                    {
                        QStringList tokens = line.split('\t');

                        //Check if file has already been processed
                        if (!alreadychecked){
                            QStringList format_col = tokens.at(FORMAT_INDEX).split(':'); //"AL:GY:ZYG:ETC"
                            if (format_col.contains(ZYG_ID)){
                                cerr << pa->colfile << " has already been processed! Printing out everything." << endl;
                                cout << line.toUtf8().data() << endl;
                                while (!in.atEnd()){
                                    cout << in.readLine().toUtf8().data() << endl;
                                }
                                exit(0);
                            }
                            alreadychecked = true;
                        }

                        //Print Header(s)
                        if (tokens[0].length() > 5) {
                            cout << line.toUtf8().data() << endl;
#ifdef BUFFERLINE_APPROACH
                            available_lines.removeFirst();
#endif
                            continue;
                        }

                        QString chrom = tokens[0].split("chr")[1];
                        uint index = tokens[1].toUInt();

//                        cerr << "\r line pos:" << chrom.toUtf8().data() << '\t' << index << "   " << flush;


                        //Assign first chromosome and index if not set;
                        if (first_chrome == ""){
                            first_chrome = chrom;
                            first_index = index;
                        }

                        //Check Range
                        if ( (first_chrome == chrom ) && ( (index - first_index) <= pa->range ))
                        {
                            lines.append(line);
                            line_data.append(new PileUp(chrom, index));
#ifndef BUFFERLINE_APPROACH
                            bufferpos = in.pos();
#endif
                        }
                        else
                        {

#ifdef BUFFERLINE_APPROACH
                            available_lines.append(bufferline);
#else
                            in.seek(bufferpos);
#endif

                            stop = true;
                            break;
                        }
                        count_lines += 1;
                    }
#ifdef BUFFERLINE_APPROACH
                    bufferline = line;                  // Update bufferline while true
                    available_lines.removeFirst();      // Remove line
                }
#else
                    bufferpos = in.pos();
#endif
            }
            cerr << "\r[" << ((int)(100*count_lines)/num_lines) << "%],P=" << line_data.length() <<  flush;

            //2. Prepare N lines for single MPileUp Call
            // Do in-place assignments to line_data
            MPileup *mp = new MPileup(QString(pa->bamfile), line_data);
            delete mp;

            // 3. Match pileup output with the input positions:
            // < MPileup updates the line_data with bases and quality data,
            //   and declares whether the line has_data or not. >

            // 4. Read stored lines of infile and for each matching <chr> <pos>:           

            for (int xx=0; xx < line_data.length(); xx++){

//                cout << lines[xx].toUtf8().data();
                PileUp *p = line_data[xx];

                //  a. parse bases and quality --> determine maxbase, zygosity.
                if (p->has_data)
                {
                    Zygosity *z = new Zygosity(p, pa->hetlim, pa->qualim);

                    QStringList tokens = lines[xx].split('\t');

                    tokens[FORMAT_INDEX].append(':').append(ZYG_ID);
                    QString build = z->max_base+"("+QString::number(z->percent)+"%)"
                            +","+(z->individual_scores)+",#"+QString::number(z->num_goodstr)
                            +","+(z->homozygous);
                    tokens[INDIV_START_INDEX].append(':').append(&build);

                    if (pa->extra){
                        tokens[FORMAT_INDEX].append(':').append(BQ_ID);
                        QString build2 = p->bases+" "+p->qualities;
                        tokens[INDIV_START_INDEX].append(':').append(&build2);
                    }

                    for (int p=0; p < tokens.length(); p++){
                        cout << tokens[p].toUtf8().data() << '\t' << flush;
                    }
                    cout << endl;


                    delete z;
                }
                else {
                    unk << lines[xx].toUtf8().data() << "\tNO BAM DATA" << endl;
                }
            }
        }
        cerr << "\nEND" << endl;
    }
    return 0;
}
