#ifndef CAVIARMODEL_H
#define CAVIARMODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <armadillo>
#include <sstream>
#include "PostCal.h"

using namespace std;
using namespace arma;

 
class CaviarModel{
	
public:
	double rho;
	double NCP;
	double gamma;
	int snpCount;
	int totalCausalSNP;
	double * sigma;
        double * stat;
        char * configure;
        int * rank;
        bool histFlag;
	PostCal * post;
	string * snpNames;
	string ldFile;
        string zFile;
        string outputFileName;
        string geneMapFile;	
        map<int,double> Weight;
        string weight;
        int nthread;
	CaviarModel(string ldFile, string zFile, string outputFileName, int totalCausalSNP, double NCP,double rho,bool histFlag,double gamma,string  weight, int nthread)           {
             //   cout<<"Creat an object"<<endl;
		int tmpSize = 0;
		this->histFlag = histFlag;
		this->NCP = NCP;
		this->rho = rho;
		this->gamma = gamma;
		this->ldFile = ldFile;
		this->zFile  = zFile;
		this->outputFileName = outputFileName;
		this->totalCausalSNP = totalCausalSNP;
                this->weight = weight;
                this->nthread= nthread;
               // cout<<"It is here1"<<endl;
		fileSize(ldFile, tmpSize);
              //  cout<<"It is here2"<<endl;
	        snpCount = (int)sqrt(tmpSize);
         	sigma     = new double[snpCount * snpCount];
		stat      = new double[snpCount];
		configure = new char[snpCount];
		rank      = new int[snpCount];
		snpNames  = new string [snpCount];
             //   cout<<"It is here3"<<endl;
		importData(ldFile, sigma);
               // cout<<"It is here4"<<endl;
		makeSigmaPositiveSemiDefinite(sigma, snpCount);
//                cout<<"It is here5"<<endl;
	//	for(int i = 0; i < snpCount*snpCount; i++)
	//		cout << sigma[i] << " ";
	//	cout << endl;
//	        cout<<"Parameters are: "<<zFile<<" "<<snpNames<<endl;
		importDataFirstColumn(zFile, snpNames);
  //              cout<<"I'm safe here2"<<endl;
		importDataSecondColumn(zFile, stat);
    //            cout<<"I'm safe here3"<<endl;
		post = new PostCal(sigma, stat, snpCount, totalCausalSNP, snpNames, gamma,nthread);
      //          cout<<"I'm safe here4"<<endl;
          //      extract_weight(weight,Weight);
        //        cout <<"The object was successfully created"<<endl;
	}
        void extract_weight(string weight, map<int,double>& Weight)
          {
            string line = "";
            string fileName=weight;
            ifstream fin(fileName.c_str(), std::ifstream::in);
            int index=-1;
            while (getline(fin, line))
             {
               string snp;
               index++;
               double weight;
               istringstream iss(line);
               iss >> snp;
               iss >> weight;
         //      cout << "Value of SNP is: " << snp << "Value of weight is: " << weight << endl;
              // weight = 1 / (1 + exp(7.2 - weight*1.3));
               weight = 1 / (1 + exp(7.2 - weight*3));
               Weight.insert(map<int, double>::value_type(index, weight));
             } 
          }
	void run() {
        	post->findOptimalSetGreedy(stat, NCP, configure, rank, rho, Weight, nthread);
	}
	void finishUp() {
               // istringstream ss_int_x1("42");
         //       cout<<"Come into finishUp1"<<endl;
//                cout<<"outputFile is: "<<endl;
//		ofstream outputFile;
           //     cout<<"Come into finishUp1_1"<<endl;
               // istringstream ss_int_x2("42");
             //   cout<<"Come into finishUp2"<<endl;
             //   cout<<"outputFileName is: "<<outputFileName<<endl;
                string outFileNameSet = string(outputFileName)+"_set";
             //   cout<<"Come into finishUp2_2"<<endl;
             //   istringstream ss_int_x3("42");
             //   cout<<"Come into finishUp3"<<endl;
                ofstream outputFile;
                outputFile.open(outFileNameSet.c_str(),ios_base::out);
             //   cout<<"Come into finishUp4"<<endl;
             //   istringstream ss_int_x4("42");
             //   cout<<"Come into finishUp5"<<endl;
                for(int i = 0; i < snpCount; i++) {
                        if(configure[i] == '1')
                                outputFile << snpNames[i] << endl;
                }
                post->printPost2File(string(outputFileName)+"_post");
                //output the histogram data to file
                if(histFlag)
                	post->printHist2File(string(outputFileName)+"_hist");
         //       istringstream ss_int_x5("42");
         //       cout<<"I am in finishUp6"<<endl;
	}
        ~CaviarModel() {
	}

};
 
#endif
