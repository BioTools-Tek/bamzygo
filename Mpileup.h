#ifndef MPILEUP_H
#define MPILEUP_H

#include <samtools/sam.h>
#include <QString>
#include <QList>
#include <QProcess>
#include <iostream>
#include <fstream>
#include <QTextStream>


#include "Zygosity.h"

//Sam buffer
typedef struct {
    int beg, end;
    samfile_t *in;
} tmpstruct_t;

/** Returns QList of pileup data **/
class MPileup {
public:
    QList<PileUp*> pileups;
    QString filer;

    MPileup(QString file_arg, QList<PileUp*> &piles)
    {
        this->filer = file_arg;
        this->pileups = piles;

        //Find range of buffered positions, and split upon change of chromosome
        QString current_chrom="";
        uint min_index = 0xFFFFFFFF; //max possible value for unsigned
        uint max_index= 0;

        for (int i=0; i< piles.length(); i++){
            PileUp *p = piles[i];

            uint pos = p->position;
            max_index = (pos >= max_index)?pos:max_index;
            min_index = (pos <  min_index)?pos:min_index;

            if ( current_chrom == "" ) current_chrom = p->chr;
            if ( (current_chrom != p->chr) || (i==(piles.length()-1) )) //New chrom or last index
            {
                //New chromosome across lines!
                QString region = current_chrom.append(":")
                        .append(QString::number(min_index)).append('-')
                        .append(QString::number(max_index));

                const char * reg = region.toUtf8().data();

                cerr << '[' << reg << "]  " << flush;

                performPileup(region);

                //Clear all
                current_chrom = p->chr;
                min_index = 0xFFFFFFFF;
                max_index = 0;
            }
        }
    }
private:
    void performPileup(QString &locus);
};

QProcess *qp = 0;

void MPileup::performPileup(QString &locus)
{
    qp = new QProcess;
    QString ret("samtools mpileup -A -r ");
    ret.append(locus).append(' ').append(this->filer);

//    cout << ret.toUtf8().data() << endl;

    qp->start(ret);
    qp->waitForFinished();

    QString result(qp->readAllStandardOutput().trimmed());
    QStringList outputlist = result.split('\n');

    for (int j=0; j < this->pileups.length(); j++)
    {
        PileUp *inpile = this->pileups[j];
        inpile->has_data = false;

        while (outputlist.length()>0)
        {
            QStringList out = outputlist.first().split('\t');

            if (out.length()<2) {
                //cerr << "BLANK AT:" << out.at(0).toUtf8().data() << " " << out.at(1).toUtf8().data() <<  endl;
                outputlist.removeFirst();
                continue;
            }

            PileUp *outpile = new PileUp(QString(out.at(0)), QString(out.at(1)).toUInt() );

            if ( inpile->chr == outpile->chr ){

                if ( outpile->position < inpile->position ){

//                    cerr << "outpilepos=" << outpile->position << " < " << inpile->position << " " << outputlist.length() << endl;
                    //Remove if input pile has greater position than output pile for speedup
                    outputlist.removeFirst();
                    continue;
                }

                if (outpile->position == inpile->position ){

                    inpile->depth = QString(out.at(3)).toInt();
                    inpile->bases = QString(out.at(4));
                    inpile->qualities = QString(out.at(5));

//                    cerr << inpile->chr.toUtf8().data() << " "
//                         << inpile->position << " -- "
//                         << inpile->depth << ":"
//                         << inpile->bases.toUtf8().data() << ":"
//                         << inpile->qualities.toUtf8().data()
//                         << endl;

                    inpile->has_data = true;

                    //Speed up future searches
                    outputlist.removeFirst();
                    continue;
                }
                else { //greater than
                    break;
                }
            }
            else {
                //current output chrome != input chrome
                cerr << "\nAnythin here" << endl;
                break;
            }
        }
    }
    qp->terminate();
    qp->close();
    qp->kill();
}


#define HFORMAT "##FORMAT="
#define FILLER ",Number=.,Type=String,Description="
#define ZYG_ID "ZYG"
#define BQ_ID "BQ"

#define HEADER_ZYGO HFORMAT "<ID=" ZYG_ID
#define HEADER_BQ HFORMAT "<ID=" BQ_ID

#define HEADER_ZYGO_FULL HEADER_ZYGO FILLER "\"Zygosity (Obtained from upstream BAM data): Max base(percent), Frequencies, Number of total stretches, Zygosity \">"
#define HEADER_BQ_FULL HEADER_BQ FILLER "\"Bases and Qualities (Obtained from upstream BAM data)\">"

static int FORMAT_INDEX = 8; // Default
static int INDIV_START_INDEX = 9; // Default


static void handleHeaders(QString &line, QTextStream &in, uint &countline, ofstream &unfound, bool baseAndQ)
{
    line="##";

    qint64 bufferpos = 0;
    bool found_zygo_header = false, found_bq_header = false;
    short found_format = -1;

    //Print headers
    do{
        line = in.readLine();
        if (line[0]!='#') break;

        //Look for ##
        if (line.startsWith(HFORMAT)){
            found_format = 0; //Found
            if (line.startsWith(HEADER_ZYGO)) found_zygo_header = true;
            else if (baseAndQ && line.startsWith(HEADER_BQ)) found_bq_header = true;
        }
        else if (line.startsWith("#CHROM")){
            QStringList tokes = line.split('\t');
            bool form = false;
            for (int p=0; p < tokes.length(); p++){
                if (tokes[p].toUpper().contains("FORMAT")) {
                    FORMAT_INDEX = p;
                    INDIV_START_INDEX = p+1;
                    form = true;
                    break;
                }
            }
            if (!form) cerr << "Could not find FORMAT column! Assuming column:" << (FORMAT_INDEX + 1) << endl;

        }


        else {
            if(found_format==0){ // Found format but now it's ended
                if (!found_zygo_header) cout << HEADER_ZYGO_FULL << endl; //Never found header, print new one
                if (baseAndQ && !found_bq_header) cout << HEADER_BQ_FULL << endl; //Never found header, print new one
                found_format=-2; // so that we know that it at least exists
            }
        }
        cout << line.toUtf8().data() << endl; // Output to rejects and cout too
        unfound <<  line.toUtf8().data() << endl;

        countline++;
        bufferpos = in.pos();
    } while(!in.atEnd());

    //Never found ##FORMAT in header
    if(found_format==-1){
        cerr << "No Format line in header!\nPrinting new one anyway..." << endl;
        if (!found_zygo_header) cout << HEADER_ZYGO_FULL << endl; //Never found header, print new one
        if (baseAndQ && !found_bq_header) cout << HEADER_BQ_FULL << endl; //Never found header, print new one
    }

    //Go to end of headers
    in.seek(bufferpos);
}



////STATIC
//// callback for bam_fetch()
//static int fetch_func(const bam1_t *b, void *data)
//{
//    bam_plbuf_t *buf = (bam_plbuf_t*)data;
//    bam_plbuf_push(b, buf);
//    return 0;
//}


////callback for bam_plbuf_init()--returns: chr pos <pileup -A>
//static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
//{
//    tmpstruct_t *tmp = (tmpstruct_t*)data;
//    if (print_all) {
////        tmp->in->header->hash

//        if ( n < bad_depth ) printf("%s\t%d\t%d -- OUT\n", tmp->in->header->target_name[tid], pos + 1, n);
//        else printf("%s\t%d\t%d\n", tmp->in->header->target_name[tid], pos + 1, n);

//    } else{
//        int current_index = (pos+1);
//        cerr << "\r" << current_index << "   " << flush;

//        PileUp *p = new PileUp(tmp->in->header->target_name[tid], current_index);
//        results.append(p);
//    }
//    return 0;
//}


//void MPileup::performPileup(const char * locus)
//{
//    const char * file = this->filer.toUtf8().data();

//    tmpstruct_t tmp;
//    tmp.beg = 0; tmp.end = 0x7fffffff;
//    tmp.in = samopen(file, "rb", 0);
//    if (tmp.in == 0) {
//        fprintf(stderr, "Fail to open BAM file %s\n", file);
//        exit(-1);
//    }
//    if (locus=="") {
//        cerr << "Bad region" << endl;
//        exit(-1);
//        //        // if a region is not specified
//        //        sampileup(tmp.in, -1, pileup_func, &tmp);
//    } else {
//        int ref;
//        bam_index_t *idx;
//        bam_plbuf_t *buf;
//        idx = bam_index_load(file); // BAM index
//        if (idx == 0) {
//            fprintf(stderr, "BAM indexing file is not available.\n");
//            exit(-1);
//        }
//        bam_parse_region(tmp.in->header, locus, &ref,
//                         &tmp.beg, &tmp.end); // parse the region
//        if (ref < 0) {
//            fprintf(stderr, "Invalid region %s\n", locus);
//            exit(-1);
//        }
//        //        cerr << "initialising bam pileup buffer" << endl;
//        buf = bam_plbuf_init(pileup_func, &tmp); // initialize pileup

//        //        cout << "fetching region = " << region << " from buffer" << endl;
//        bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, buf, fetch_func);

//        bam_plbuf_push(0, buf);
//        bam_index_destroy(idx);
//        bam_plbuf_destroy(buf);
//    }
//    samclose(tmp.in);
//}



#endif // MPILEUP_H
