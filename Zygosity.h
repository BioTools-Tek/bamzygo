#ifndef ZYGOSITY_H
#define ZYGOSITY_H

#include <qmap.h>

using namespace std;

//When number Insertions ~ Deletions, this is the percentage assigned
#define FIFTYFIFTY 911

/** Simply a container for pileup data
  **/
class PileUp{
public:
    QString chr;
    uint position; uint depth;
    QString bases;
    QString qualities;

    bool has_data; //

    PileUp(QString chr, uint position)
    {
        this->chr = chr; this->position = position;
    }
};


class Zygosity
{
public:
    QString max_base;
    QString homozygous;
    int percent;
    bool indel;
    QString individual_scores;
    int num_goodstr;

    Zygosity(PileUp *p, int hetlim, int qualim)
    {
        QString bases = p->bases;
        QString quality = p->qualities;
        QMap<QString, int> map_scores;
        QMap<QString, int> map_freq;

        individual_scores="";

        int qlength = quality.length();
        int blength = bases.length();
        num_goodstr=0;

        int num_insertions = 0;
        int num_deletions = 0;

        //        cerr << bases.toUtf8().data() << "  " << quality.toUtf8().data() << endl;

        int qindex = 0;
        int bindex = 0;


        while( qindex < qlength )
        {
            int score = quality[qindex].toLatin1()-33;
            if ( score > qualim ) num_goodstr ++;

            // Skip Base-Specific Scores
            if (bases[bindex]=='^') bindex +=2;
            else if (bases[bindex]=='$') bindex += 1;

            //Look ahead
            QString change = QString(bases[bindex]);
            bool insertion = false;
            bool deletion = false;

            int old_bindex = bindex;
            if (bindex+1 < blength)
            {
                change = QString(bases[bindex+1]);
                // We have insertion or deletion!

// Things to note:
// Actual deletions are represented by '*' in the pileup, i.e. they display empty spaces.
// VCF, instead of bisecting deletions randomnly in the middle, tend to display just where they start and
// how long they are.
// So if chr1 56789  has TTttTT-4NNNNt-4NNNNttT-4NNNN, then these three deletions do not actually appear for position 56789
// but instead appear at position 56790 as AAaA***a*A***, etc.
// So technically position 56789 is homozygous T's, but refers to the beginning of a heterozygous indel

                if (change=="+" || change=="-")
                {
                    insertion = !(deletion = (change=="-"));
                    bindex ++; //We now move into the scope of the indel

                    //Find number
                    int num_length = 0;
                    bool ok = true;
                    while (ok){
                        QString(bases[(bindex+1)+num_length]).toInt(&ok);
                        num_length ++;
                    }
                    int num = QString(bases.mid(bindex+1, num_length-1)).toInt();

                    //Extract the previous and the symbol
                    change.prepend(QChar(bases[old_bindex])); // Add previous. Change = a+ now
                    //              cout << "change1=" << change.toUtf8().data() << endl;

                    ///// BINDEX is now pointing at '+'

                    //Place number
                    int yy = 1;
                    while (yy < num_length) {change.append(bases[bindex+yy]); yy ++;}  // Change = a+1 now
                    bindex += yy;
                    //              cout << "change2=" << change.toUtf8().data() << endl;

                    //Place insertion bases
                    int jj = 0;
                    while (jj < num){
                        change.append(bases[(bindex+jj)]);
                        jj ++;
                    }
                    //            cout << "change3=" << change.toUtf8().data() << endl;

                    bindex += jj-1;
                }
                else{
                    //No insertions or deletions
                    bindex = old_bindex;
                    change = bases[bindex];
                }
            }

            if ( score > qualim ) {
#ifdef DEBUG_BASES
                cout << change.toUtf8().data() << "  " << score << ", " << flush;
#endif
                if(insertion) num_insertions ++;
                if(deletion) num_deletions ++;

                //Store score
                map_scores[change] = score;

                //Store Frequencies
                change = change.toUpper();
                if (map_freq.contains(change)) map_freq[change] +=1;
                else map_freq[change] = 1;
            }

            bindex ++;
            qindex ++;
        }

        //Check they both end at the end of their strings;
        if ( (bindex < (blength-5)) || (qindex < qlength)){
            cerr << "\nbindex=" << bindex << "  " << blength << endl;
            cerr << "qindex=" << qindex << "  " << qlength << endl;
            cerr << "AT: chr" << p->chr.toUtf8().data() << "\t" << p->position << endl;
//            exit(-1);
        }

#ifdef DEBUG_BASES
        cout << endl;
#endif
        ////Find and Set Max Base

        QList<QString> keys = map_freq.keys();
        int max=0;

        for (int i=0; i< keys.length(); i++)
        {
            QString joke = keys[i];
            int freq = map_freq[joke];

#ifdef DEBUG_SCORES
            cout << joke.toUtf8().data() << "[]" << freq << endl;
#endif

            individual_scores.append(joke+"["+QString::number(freq)+"] ");

            if (joke.length()>1){ // e.g. C-4NNNN

                //Mark indels sperately when finding max base
                if (joke[1]=='-' || joke[1]=='+'){
                    //keys.removeAt(i);
                    continue;
                }
            }
            // Otherwise find max SNV
            else if ( freq > max ){
                max = freq;
                this->max_base = joke;
            }
        }
        //// Find Percentage and Homozygosity
        // Check for indels
        if ((num_insertions>0) || (num_deletions>0))
        {
            indel = true;

            int perc_insertions = (int)(100*((float)(num_insertions))/((float)(num_goodstr)));
            int perc_deletions = (int)(100*((float)(num_deletions))/((float)(num_goodstr)));

#ifdef DEBUG_BASES
            cout << "Num_GoodSTR  " << num_goodstr << endl;
            cerr << "\ninsertions:" << num_insertions << ", perc:" << perc_insertions << endl;
            cerr << "deletions :" << num_deletions << ", perc:" << perc_deletions << endl;
#endif
            //Pick majority
            if ( num_insertions > num_deletions ){
                this->percent = perc_insertions;
                this->max_base = "+";
            } else {
                this->percent = perc_deletions;
                this->max_base = "-";
            }

            //Check that insertions and deletions aren't of equivalent majority (~20% diff)
            float ratio = (
                   (perc_insertions>perc_deletions)?
                   ((float)(perc_deletions)/(float)(perc_insertions)):
                   ((float)(perc_insertions)/(float)(perc_deletions))
            );
            if ( (ratio >= 0.7) ){
                cerr << "\t[" << p->chr.toUtf8().data() << ':' << p->position << ']'
                    << "#Ins(" << perc_insertions << ")~#Del("<< perc_deletions << ')' << flush;
                this->percent = FIFTYFIFTY;
                this->max_base = "+/-";
            }
        }
        // Dealing with snvs
        else {
            indel = false;
            this->percent =  (int)(100*((float)(max))/((float)(num_goodstr)));
        }

        //Determine Zygosity String
        if (this->percent == FIFTYFIFTY ) this->homozygous = "CHECK";
        else if (this->percent < (100-hetlim)) this->homozygous = "HETEROZY?";
        else if (this->percent >= hetlim ) this->homozygous = "HOMOZY";
        else this->homozygous = "HETEROZY";
    }
};


#endif // ZYGOSITY_H
