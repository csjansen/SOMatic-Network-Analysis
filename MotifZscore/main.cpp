#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <vector>
#include <map>
#include <chrono>
#include <thread>
#include <cmath>

using namespace std;

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

//Move to util header
#define PI 3.1415926

float inverfc(float x) {
        //float u = -2.0/(log(PI*pow(x,2)*log(1.0/x)));
        //float v = log(log(1.0/x))-2+log(PI);
        //float a2 = 1.0/8.0*v;
        //float a3 = -1.0/32.0*(pow(v,2)+6*v-6);
        //float a4 = 1.0/384.0*(4*pow(v,3)+27*pow(v,2)+108*v-300);
        //return pow(u,-1.0/2.0)+a2*pow(u,3.0/2.0)+a3*pow(u,5.0/2.0)+a4*pow(u,7.0/2.0);
//      float t = .5*sqrt(PI)*x;
//      return t+1.0/3*pow(t,3)+7.0/30*pow(t,5)+127.0/630*pow(t,7);
        float w, p;
w = - logf((1.0f-x)*(1.0f+x));
if ( w < 5.000000f ) {
w = w - 2.500000f;
p = 2.81022636e-08f;
p = 3.43273939e-07f + p*w;
p = -3.5233877e-06f + p*w;
p = -4.39150654e-06f + p*w;
p = 0.00021858087f + p*w;
p = -0.00125372503f + p*w;
p = -0.00417768164f + p*w;
p = 0.246640727f + p*w;
p = 1.50140941f + p*w;
}
else {
w = sqrtf(w) - 3.000000f;
p = -0.000200214257f;
p = 0.000100950558f + p*w;
p = 0.00134934322f + p*w;
p = -0.00367342844f + p*w;
p = 0.00573950773f + p*w;
p = -0.0076224613f + p*w;
p = 0.00943887047f + p*w;
p = 1.00167406f + p*w;
p = 2.83297682f + p*w;
}
return p*x;
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

vector <vector<int> > merge(vector <vector<int> > lhs, vector<vector<int> > rhs) {
    vector <vector<int> > result;
    result.insert(result.begin(),lhs.begin(), lhs.end());
    for(int i = 0; i < rhs.size(); i++) {
        bool found = false;
        for(int j = 0; j < lhs.size(); j++) {
            if(rhs[i] == lhs[j]){
                found = true;
                break;
            }
        }
        if(!found) result.push_back(rhs[i]);
    }
    return result;
}
struct Region {
        string chr;
        int start;
        int stop;
        //string qval;
        vector<string> Gene;
        string strand;
};
int main(int argc, char *argv[]) {
        if(argc < 2) {
        cout << "Usage: ./MotifZScore [options] -Metacluster1 <Number of metaclusters in SOM1> -Metacluster2 <Number of metaclusters in SOM2> -FusionBreakup <Location of the fusion breakup folder> -Analysis <Name of Motif Analysis> -ZScoreOutput <Main Output file> -AllOutput <Outputing all connections> -pval <pval cutoff>" <<endl;
        cout << "Options: <default>" <<endl;
        cout << "-MotifSuffix: If fimo.txt is not default, you can add a sufix. <>"<<endl;
	cout << "-Fimo: fimo version if not 5.1.1.  Supports 4.12"<<endl;
        return 0;
    }
        int Metacluster1 = 0;
        int Metacluster2 = 0;
        string FusionBreakup="";
        string Analysis="";
        string ZScoreOutputFileName="";
        string AllOutputFileName="";
        string MotifSuffix = "";
        float ZScore=0;
        float pval=0;
	string Fimo = "5.1.1";
        for(int i = 0; i < argc; i++) {
        string temp = argv[i];
                if(temp.compare("-Metacluster1")==0)
            istringstream(argv[i+1])>>Metacluster1;
                if(temp.compare("-Metacluster2")==0)
            istringstream(argv[i+1])>>Metacluster2;
                if(temp.compare("-FusionBreakup")==0)
            FusionBreakup = argv[i+1];
                if(temp.compare("-Analysis")==0)
            Analysis = argv[i+1];
                if(temp.compare("-ZScoreOutput")==0)
            ZScoreOutputFileName = argv[i+1];
                if(temp.compare("-AllOutput")==0)
            AllOutputFileName = argv[i+1];
                if(temp.compare("-MotifSuffix")==0)
            MotifSuffix = argv[i+1];
                if(temp.compare("-pval")==0)
            istringstream(argv[i+1])>>pval;
		if(temp.compare("-Fimo")==0)
		Fimo = argv[i+1];
        }
        if(pval > 0) {
//              ZScore = sqrt(2)*inverfc(2*pval);
               	pval = pval/(double)(Metacluster1*Metacluster2);
                ZScore = -1*sqrt(2)*inverfc(2*pval-1);
        }
        cout<<ZScore<<endl;
        vector<string> TFList;
        map<string,int> TFDict;
        vector<vector<vector<vector<Region> > > > RegionsForEachTFInEachList;
        //vector<vector<int> > TotalSizeInEachList;
        //vector<vector<int> > TotalRegionsInEachList;
        vector<vector<vector<map<string,int> > > > NumOfMatches;
//        vector<vector<vector<int> > > NumOfMatches2;
        vector<vector<vector<vector<int> > > > nums;

        RegionsForEachTFInEachList.resize(Metacluster1);
       // TotalRegionsInEachList.resize(Metacluster1);
       // TotalSizeInEachList.resize(Metacluster1);
        nums.resize(Metacluster1);
        NumOfMatches.resize(Metacluster1);
        for(int i = 0; i < Metacluster1; i++) {
                RegionsForEachTFInEachList[i].resize(Metacluster2);
         //       TotalRegionsInEachList[i].resize(Metacluster2,0);
         //       TotalSizeInEachList[i].resize(Metacluster2,0);
                nums[i].resize(Metacluster2);
                NumOfMatches[i].resize(Metacluster2);
        }
        int motifNum = 0;
for(int i = 0; i < Metacluster1; i++) {
  //            i=9;
                for(int j = 0; j < Metacluster2; j++) {
//                      j=18;
                        string filename = Analysis+"_"+SSTR(i)+"_"+SSTR(j)+"_fimo/fimo.tsv";
			if(Fimo.compare("4.12")==0) {
				filename = Analysis+"_"+SSTR(i)+"_"+SSTR(j)+"_fimo/fimo.txt";
			}
                        cout<<filename<<endl;
                        ifstream infile(filename.c_str());
                        if(!infile.good())
                                continue;
                        string line;
			bool first = true;
                        while(getline(infile,line)) {
				if(first) {
					first=false;
					continue;
				}
                                if(line[0]=='#')
                                        continue;
                                vector<string> splitz = split(line,'\t');
//                                cout<<line<<endl;
				if(line=="") {
					continue;
				}
				if(splitz.size() < 9) {
					continue;
				} 
                                string TF = splitz[0];

//                              cout<<TF<<endl;
                                map<string,int>::iterator it = TFDict.find(TF);
                                int pos;
                                if(it==TFDict.end()) {
                                        TFList.push_back(TF);
                                        TFDict[TF]=TFList.size()-1;
                                        pos = TFList.size()-1;
//                                      cout<<pos<<endl;
//                                      cout<<RegionsForEachTFInEachList.size()<<endl;
//                                      cout<<RegionsForEachTFInEachList[0].size()<<endl;
                                        for(int k = 0; k < Metacluster1; k++) {
                                                for(int m = 0; m < Metacluster2; m++) {
                                                        vector<Region> temp;
                                                        Region temp3;
                                                        int num = 0;
                                                        vector<int> temp2;
//                                                      cout<<RegionsForEachTFInEachList[k].size()<<'\t'<<m<<endl;
                                                        RegionsForEachTFInEachList[k][m].push_back(temp);
                                                        //cout<<k<<'\t'<<m<<endl;
                                                        map<string,int> tempmap;
                                                        NumOfMatches[k][m].push_back(tempmap);
                                                        nums[k][m].push_back(temp2);
                                                }
                                        }
                                } else {
                                        pos = (int)it->second;
                                }
                                Region temp;
//                                cout<<splitz[1]<<'\t'<<splitz[2]<<'\t'<<splitz[3]<<'\t'<<splitz[4]<<'\t'<<splitz[5]<<endl;
                                temp.chr = splitz[2];
                                istringstream(splitz[3])>>temp.start;
                                istringstream(splitz[4])>>temp.stop;
                                //temp.qval = splitz[8];
                                temp.strand = splitz[5];
//				cout<<temp.strand<<endl;
                                RegionsForEachTFInEachList[i][j][pos].push_back(temp);
                                motifNum++;
                        }
                        infile.close();
        //              break;
                }
        //      break;
        }
	TFDict.clear();
        cout<<"Motif Number: "<<TFList.size()<<endl;
	int NumberOfRegions[Metacluster1][Metacluster2];
//	vector<vector<int>> NumberOfRegions;
	int numNonEmpty = 0;
	for(int i = 0; i < Metacluster1; i++) {
		vector<int> tempRegions;
                for(int j = 0; j < Metacluster2; j++) {
                        string filename = FusionBreakup+"Combo_"+SSTR(i)+"_"+SSTR(j);
                        cout<<filename<<endl;
                        ifstream infile(filename.c_str());
                        if(!infile.good()) {
				cout << "Couldn't find.  Continuing"<<endl;
                                continue;
			}
                        string line;
			vector<Region> Regions;
                        while(getline(infile,line)) {
                                //cout<<line<<endl;
                                vector<string> splitz = split(line,'\t');
				//vector<string> splitz2 = split(splitz[0],':');
				//vector<string> splitz3 = split(splitz2[1],'-');
                                Region temp;
                                temp.chr = splitz[0];
                        //      cout<<temp.chr<<endl;
                                istringstream(splitz[1])>>temp.start;
                                istringstream(splitz[2])>>temp.stop;
				int found = -1;
				for(int k = 0; k < Regions.size(); k++) {
					if(Regions[k].chr==temp.chr && Regions[k].start==temp.start && Regions[k].stop==temp.stop) {
						found = k;
						break;
					}
				}
				if(found==-1) {
                                	temp.Gene.push_back(splitz[3]);
                                	for(int k = 0; k < TFList.size(); k++) {
                                        	NumOfMatches[i][j][k][splitz[3]]=0;
                                	}
					Regions.push_back(temp);
				} else {
					Regions[found].Gene.push_back(splitz[3]);
					for(int k = 0; k < TFList.size(); k++) {
                                                NumOfMatches[i][j][k][splitz[3]]=0;
                                        }
				}
                        }
			NumberOfRegions[i][j]=Regions.size();
//			tempRegions.push_back(Regions.size());
			if(Regions.size()>0) {
				numNonEmpty++;
			}
                        infile.close();
        		for(int k = 0; k < TFList.size(); k++) {
                		cout<<i<<'\t'<<j<<'\t'<<TFList[k]<<'\t'<<k+1<<'/'<<TFList.size()<<endl;
                                for(int n = 0; n < Regions.size(); n++) {
                                        for(int m = 0; m < RegionsForEachTFInEachList[i][j][k].size(); m++) {
			//		cout<<Regions[n].chr<<'\t'<<RegionsForEachTFInEachList[i][j][k][m].chr<<endl;
                                                if(Regions[n].chr.compare(RegionsForEachTFInEachList[i][j][k][m].chr)==0 && Regions[n].start <= RegionsForEachTFInEachList[i][j][k][m].start && Regions[n].stop >= RegionsForEachTFInEachList[i][j][k][m].stop){
				//			cout<<m<<endl;
							//RegionsForEachTFInEachList[i][j][k][m].Gene=Regions[n].Gene;
							for(int t = 0; t < Regions[n].Gene.size(); t++) { 
                                                        	NumOfMatches[i][j][k][Regions[n].Gene[t]]++;
							}
                                                        nums[i][j][k].push_back(m);
                                                }
                                        }
                                }
                        }
                }
//		NumberOfRegions.push_back(tempRegions);
        }
        ofstream outfile(ZScoreOutputFileName.c_str());
        ofstream outfile2(AllOutputFileName.c_str());
        for(int i = 0; i < TFList.size(); i++) {
                double total = 0;
                int totalskipped = 0;
                int totalCounts = 0;
                vector<vector<double> > ValuesMat;
                for(int j = 0; j < Metacluster1; j++) {
                //      j=37;
                        vector<double> ValuesRow;
                        for(int k = 0; k < Metacluster2; k++) {
                		cout<<j<<'\t'<<k<<'\t'<<TFList[i]<<'\t'<<i+1<<'/'<<TFList.size()<<endl;
                //              k=65;
                                //totalCounts++;
                                double total1 = 0;
                                vector<string> keys;
                                for(auto const& imap: NumOfMatches[j][k][i])
                                        keys.push_back(imap.first);
                                int totalcount = 0;
                                for(int n = 0; n < keys.size(); n++) {
                                        if(NumOfMatches[j][k][i][keys[n]]>0) {
                                                if(keys.size() > 0) {
                                                        total1+= 1/(double)keys.size();
                                                        totalcount++;
                                                }
                                        }
                                        //cout<<keys[n]<<'\t'<<NumOfMatches[j][k][i][keys[n]]<<'\t'<<total1<<endl;
                                }
                                //if(totalcount==1) total1=0;
                                ValuesRow.push_back(total1);
                                total+=total1;
                //              break;
                        }
                        ValuesMat.push_back(ValuesRow);
                //      break;
                }
                double average = total/(double)(numNonEmpty);
                double stddevtotal = 0;
                for(int j = 0; j < ValuesMat.size(); j++) {
                        for(int k = 0; k < ValuesMat[j].size(); k++) {
                                stddevtotal+=pow(ValuesMat[j][k]-average,2);
                        }
                }
                double stddev = sqrt(stddevtotal/(float)(Metacluster1*Metacluster2));
                for(int j = 0; j < Metacluster1; j++) {
                        for(int k = 0; k < Metacluster2; k++) {
				double zscore = 0;
				if(stddev > 0)
                                	zscore = (ValuesMat[j][k]-average)/(sqrt(average*(1-average)/(double)(NumberOfRegions[j][k])));
//				if(nums[j][k][i].size()<=1)
//					zscore = 0;
                                if(zscore > ZScore) {
                                        outfile<<TFList[i]<<'\t'<<j<<'\t'<<k<<'\t'<<ValuesMat[j][k]<<'\t'<<zscore;
                                        for(int n = 0; n < nums[j][k][i].size();n++) {
						Region R = RegionsForEachTFInEachList[j][k][i][nums[j][k][i][n]];
                                                outfile<<'\t'<<R.chr<<":"<<R.start<<"-"<<R.stop<<"-"<<R.strand;
                                        }
                                        outfile<<endl;
                                }
                                outfile2<<TFList[i]<<'\t'<<j<<'\t'<<k<<'\t'<<ValuesMat[j][k]<<'\t'<<zscore<<'\t'<<average<<'\t'<<NumberOfRegions[j][k];
                                for(int n = 0; n < nums[j][k][i].size();n++) {
					Region R = RegionsForEachTFInEachList[j][k][i][nums[j][k][i][n]];
                                        outfile2<<'\t'<<R.chr<<":"<<R.start<<"-"<<R.stop<<"-"<<R.strand;
                                }
                                outfile2<<endl;
                        }
                }
        }

        outfile.close();
        outfile2.close();
}

