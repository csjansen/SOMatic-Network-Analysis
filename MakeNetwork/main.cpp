#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<map>
#include<string>
#include<math.h>
#include<algorithm>


using namespace std;

class TSSsite {
public:
    string chrom;
    int pos;
    string strand;
    TSSsite() {
        chrom = "";
        pos = 0;
        strand = "";
    }
};

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()


class genomicRegion {
public:
    string gene;
    string chrom;
    int start;
    int stop;
    string strand;
};


template <typename T>
  string NumberToString ( T Number )
  {
     ostringstream ss;
     ss << Number;
     return ss.str();
  }

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
        return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
        split(s, delim, elems);
        return elems;
}
map<string, TSSsite>* parseGtfFile(string gtfFileName) {
    ifstream gtffile(gtfFileName.c_str());
    map<string, TSSsite>* TSSsites=new map<string, TSSsite>();
    string line;

    while(getline(gtffile,line)) {
        if(line[0]=='#') continue;
        vector<string> splitz = split(line, '\t');
        if(splitz[2].compare("exon")==0) {
            string geneName;
            vector<string> attributes = split(splitz[8], ';');
            for(int i = 0; i < attributes.size(); i++) {
                vector<string> pairItems = split(attributes[i],' ');
//                              cout<<pairItems[0]<<endl;
                if(pairItems[1].compare("gene_name")==0) {
//                                      cout<<pairItems[1]<<endl;
                    geneName = pairItems[2].substr(1,pairItems[2].size()-2);
                    break;
                }
            }
                        //cout<<geneName<<endl;
            TSSsite record = TSSsites->operator[](geneName);
            record.strand = splitz[6];
            record.chrom = splitz[0];
            int pos;
            if(record.chrom.compare("+")==0) {
                istringstream(splitz[3])>>pos;
                if(pos < record.pos || record.pos == 0) record.pos = pos;
            } else {
                istringstream(splitz[4])>>pos;
                if(pos > record.pos || record.pos == 0) record.pos = pos;
            }
            TSSsites->operator[](geneName) = record;
        }
    }

    return TSSsites;
}
vector<genomicRegion> GetRegRegions(map<string, TSSsite>* TSSsites, vector<string> genes, string CompareType) {
    vector<genomicRegion> regions;
    for(int i = 0; i < genes.size(); i++) {
                if(i%10000==0) cout<<i<<" Done"<<endl;
        map<string, TSSsite>::iterator it = TSSsites->find(genes[i]);
        TSSsite TSS;
        if(it!=TSSsites->end()) {
            TSS = it->second;
        } else {
            cout<<"Could not find: "<<genes[i]<<".  Skipping..."<<endl;
            continue;
        }
        int start;
        int stop;
        if(CompareType.compare("TwoClosest")==0) {
            int prevdist = 50000;
            if(TSS.pos < 50000) start = 0;
            else start = TSS.pos - 50000;
            int nextdist = 50000;
            stop = TSS.pos + 50000;

            for(map<string, TSSsite>::iterator looper = TSSsites->begin(); looper != TSSsites->end(); looper++) {
                TSSsite tester = looper->second;
                if(TSS.chrom.compare(tester.chrom)==0 && ((TSS.pos<tester.pos&&tester.pos-TSS.pos<nextdist)||(TSS.pos>tester.pos&&TSS.pos-tester.pos<prevdist))) {
                    if(TSS.pos<tester.pos) {
                        stop = tester.pos;
                        nextdist = tester.pos-TSS.pos;
                    } else {
                        start = tester.pos;
                        prevdist = TSS.pos - tester.pos;
                    }
                }
            }
            genomicRegion temp;
            temp.gene = genes[i];
            temp.chrom = TSS.chrom;
            temp.start = start;
            temp.stop = stop;
            regions.push_back(temp);
        }
                if(CompareType.compare("OneClosest")==0) {
                        float distance1=50000;
            float distance2=50000;
            for(map<string, TSSsite>::iterator looper = TSSsites->begin(); looper != TSSsites->end(); looper++) {
                TSSsite tester = looper->second;
                if(TSS.chrom.compare(tester.chrom)!=0) continue;
                if(TSS.pos < tester.pos && tester.pos-TSS.pos < distance1) distance1=tester.pos-TSS.pos;
                if(TSS.pos > tester.pos && TSS.pos-tester.pos < distance2) distance2=TSS.pos-tester.pos;
            }
            if(distance2!=50000)
                start = TSS.pos-distance2/2;
            else
                start = TSS.pos-distance2;
            if(distance1!=50000)
                stop = TSS.pos+distance1/2;
            else
                stop = TSS.pos+distance1;
            genomicRegion temp;
            temp.gene = genes[i];
            temp.chrom = TSS.chrom;
            temp.start = start;
            temp.stop = stop;
            regions.push_back(temp);

                }
    }
    return regions;
}
struct Combo {
        string chr;
        int start;
        int stop;
        string gene;
};

int main(int argc, char* argv[]) {
        if(argc < 2) {
        cout << "Usage: ./MakeNetwork [options] -Metacluster1 <number of metaclusters from SOM1> -Metacluster2 <number of metaclusters from SOM2> -FusionBreakup <Fusion Breakup File Location> -ZScoreFile <Input Z File Location> -MotifFile <tsv from HOCOMOCO with all motifs> -Output <Output File Location>" <<endl;
        cout << "Options: <default>" <<endl;
        cout << "-ChIP: Use ChIP peak beds for the Zscore file.  Make MotifFile into sample names file."<<endl;
        return 0;
    }
        string inputprefix;
        int row1;
    int col1;
    string FusionBreakupFileName;
        string outputprefix;
        string motifFileName;
        string XenoGeneFileName;
        bool ChIP=false;
        for(int i = 0; i < argc; i++) {
        string temp = argv[i];
                if(temp.compare("-Metacluster1")==0)
            istringstream(argv[i+1])>>row1;
                if(temp.compare("-Metacluster2")==0)
            istringstream(argv[i+1])>>col1;
                if(temp.compare("-FusionBreakup")==0)
                        FusionBreakupFileName = argv[i+1];
                if(temp.compare("-Output")==0)
                        outputprefix = argv[i+1];
                if(temp.compare("-ZScoreFile")==0)
                    inputprefix = argv[i+1];
                if(temp.compare("-MotifFile")==0)
                        motifFileName = argv[i+1];
                if(temp.compare("-XenoGene")==0)
                        XenoGeneFileName=argv[i+1];
                if(temp.compare("-ChIP")==0)
                        ChIP=true;
        }
        string ChIPName;
        if(ChIP) {
                ChIPName = motifFileName;
        }
        map<string,string> XenoGene;
        if(XenoGeneFileName.compare("")!=0) {
                string line;
                ifstream XenoGeneFile(XenoGeneFileName.c_str());
                while(getline(XenoGeneFile,line)) {
                        vector<string> splitz = split(line,'\t');
                        transform(splitz[2].begin(), splitz[2].end(), splitz[2].begin(), ::toupper);
                        XenoGene[splitz[0]]=splitz[2];
                }
                XenoGeneFile.close();
        }

        cout<<"Parsing Motif File"<<endl;
        ifstream motifFile(motifFileName.c_str());
        string line;
        map<string,string> IDtoGeneName;
	bool first = true;
        while(getline(motifFile,line)) {
		if(first) {
			first = false;
			continue;
		}
                if(!ChIP) {
                        vector<string> splitz = split(line,'\t');
//			cout<<splitz[0]<<'\t'<<splitz[2]<<endl;
//vector<string> splitz2 = split(splitz[0],'.');
                        IDtoGeneName[splitz[0]]=splitz[2];
                        cout<<splitz[0]<<'\t'<<splitz[2]<<endl;
                }
        }
        cout<<"Parsing Fusion Combo Files"<<endl;
//        cout<<"FOXI1_MOUSE.H10MO.B"<<'\t'<<IDtoGeneName["FOXI1_MOUSE.H10MO.B"]<<endl;
        vector<vector<vector<Combo> > > Combos;
        for(int i = 0; i < row1; i++) {
                vector<vector<Combo> > ComboRow;
                for(int j = 0; j < col1; j++) {
                        vector<Combo> ComboEle;
                        cout<<(FusionBreakupFileName+"/Combo_"+SSTR(i)+"_"+SSTR(j))<<endl;
                        ifstream infile((FusionBreakupFileName+"/Combo_"+SSTR(i)+"_"+SSTR(j)).c_str());
                        string line;
                        while(getline(infile,line)) {
//                              cout<<line<<endl;
                                vector<string> splitz = split(line,'\t');
				//vector<string> splitz2 = split(splitz[0],':');
				//vector<string> splitz3 = split(splitz2[1],'-');
                                Combo temp;
                                temp.chr = splitz[0];
                                istringstream(splitz[1])>>temp.start;
                                istringstream(splitz[2])>>temp.stop;
                                //vector<string> splitz2 = split(splitz[3],'_');
                                //temp.gene=splitz2[1];
                                if(XenoGene[splitz[3]].compare("")==0)
                                        temp.gene = splitz[3];
                                else
                                        temp.gene = XenoGene[splitz[3]];
                                ComboEle.push_back(temp);
                        }
                        ComboRow.push_back(ComboEle);
                }
                Combos.push_back(ComboRow);
        }
        ifstream BedFile(inputprefix.c_str());
cout<<inputprefix<<endl;
    ofstream OutBedFile(outputprefix.c_str());
        while(getline(BedFile,line)) {
                vector<string> splitz = split(line,'\t');
                if(splitz.size() < 6) continue;
                if(!ChIP) {
//			cout<<line<<endl;
                        string TFName = splitz[0];
                        int row;
                        int col;
                        istringstream(splitz[1])>>row;
                        istringstream(splitz[2])>>col;
                        for(int i = 6; i < splitz.size(); i++) {
//                              cout<<splitz[i]<<endl;
                                genomicRegion temp;
                                vector<string> splitz2 = split(splitz[i],':');
//				cout<<"1"<<endl;
                                vector<string> splitz3 = split(splitz2[1],'-');
//				cout<<"1"<<endl;
                                temp.chrom = splitz2[0];
//				cout<<"1"<<endl;
                                istringstream(splitz3[0])>>temp.start;
//				cout<<"1"<<endl;
                                istringstream(splitz3[1])>>temp.stop;
//				cout<<"1"<<endl;
                                temp.strand=(splitz2[1])[splitz2[1].length()-2];
//				cout<<splitz2[1]<<endl;
//				cout<<temp.strand<<endl;
//				cout<<"1"<<endl;
                                vector<string> foundgenes;
//                              cout<<temp.stop<<endl;

//                                      cout<<row<<'\t'<<col<<endl;
  //                                    cout<<Combos[row][col].size()<<endl;
                                for(int j = 0; j < Combos[row][col].size(); j++) {
                                        Combo test=Combos[row][col][j];
    //                                  cout<<test.chr<<endl;
                                        string chrtemp = "chr";
                                        if((temp.chrom).compare(test.chr)==0 && ((temp.start<=test.start && temp.stop >= test.start)||(temp.start <= test.stop && temp.stop >= test.stop)||(temp.start>=test.start && temp.stop <= test.stop)||(temp.start<=test.start&&temp.stop>=test.stop))) {
                                         //       cout<<TFName<<'\t'<<IDtoGeneName[TFName]<<endl;
                                       //         if(IDtoGeneName[TFName].compare("")==0)
                                         //               OutBedFile<<TFName<<'\t'<<TFName<<"\t->\t"<<test.gene<<'\t'<<temp.chrom<<'\t'<<temp.start<<'\t'<<temp.stop<<'\t'<<temp.strand<<'\t'<<splitz[1]<<'\t'<<splitz[2]<<'\t'<<splitz[3]<<'\t'<<splitz[4]<<endl;
                                        //        else {
//							cout<<splitz.size()<<endl;
//							cout<<IDtoGeneName[TFName]<<"\t->\t"<<test.gene<<'\t'<<temp.chrom<<'\t'<<temp.start<<'\t'<<temp.stop<<'\t'<<temp.strand<<'\t'<<splitz[1]<<'\t'<<splitz[2]<<'\t'<<splitz[3]<<'\t'<<splitz[4]<<'\t'<<splitz[5]<<endl;
                                                       OutBedFile<<TFName<<'\t'<<IDtoGeneName[TFName]<<"\t->\t"<<test.gene<<'\t'<<temp.chrom<<'\t'<<temp.start<<'\t'<<temp.stop<<'\t'<<temp.strand<<'\t'<<splitz[1]<<'\t'<<splitz[2]<<'\t'<<splitz[3]<<'\t'<<splitz[4]<<endl;
					//	}
                                                break;
						

                                        }
                                }

                        }
                } else {
                        genomicRegion temp;
                        temp.chrom = splitz[0];
                        istringstream(splitz[1])>>temp.start;
                        istringstream(splitz[2])>>temp.stop;
                        bool found=false;
                        for(int i = 0; i < row1&&!found; i++) {
                                for(int j = 0; j < col1&&!found; j++) {
                                        for(int k = 0; k < Combos[i][j].size(); k++) {
                                                Combo test = Combos[i][j][k];
                                                if(temp.chrom.compare(test.chr)==0 && ((temp.start<=test.start && temp.stop >= test.start)||(temp.start <= test.stop && temp.stop >= test.stop)||(temp.start>=test.start && temp.stop <= test.stop)||(temp.start<=test.start&&temp.stop>=test.stop))) {
                                                        double coverage = -1*(max(test.start,temp.start)-min(test.stop,temp.stop))/(double)(temp.stop-temp.start);
                                                        //cout<<test.start<<'\t'<<test.stop<<'\t'<<temp.start<<'\t'<<temp.stop<<'\t'<<coverage<<endl;
                                                        if(coverage>.5) {

                                                                found = true;
                                                                OutBedFile<<ChIPName<<"\t->\t"<<test.gene<<'\t'<<test.chr<<'\t'<<temp.start<<'\t'<<temp.stop<<'\t'<<i<<'\t'<<j<<'\t'<<endl;
                                                                break;
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
        BedFile.close();
        OutBedFile.close();
}

