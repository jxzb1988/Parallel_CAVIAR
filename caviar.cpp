#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <unistd.h> 

#include "Util.h"
#include "PostCal.h"
#include "TopKSNP.h"
#include "CaviarModel.h"

using namespace std;


int main( int argc, char *argv[]  ){
        //return(1);
	int totalCausalSNP = 2;
	double NCP = 5.7;
	double gamma = 0.01;
	double rho = 0.95;
	bool histFlag = false;
	int oc = 0;	
	string ldFile = "";
	string zFile  = "";
	string outputFileName = "";
	string geneMapFile = "";	
        string weight = "";
        int nthread=1;
	while ((oc = getopt(argc, argv, "vhl:t:o:z:w:g:r:c:w:f:")) != -1) {
		switch (oc) {
			case 'v':
				cout << "version 1.0:" << endl;
			case 'h':
				cout << "Options: " << endl;
  				cout << "-h, --help            		show this help message and exit " << endl;
  				cout << "-o OUTFILE, --out=OUTFILE 	specify the output file" << endl;
  				cout << "-l LDFILE, --ld_file=LDFILE  	the ld input file" << endl;
  				cout << "-z ZFILE, --z_file=ZFILE	the z-score and rsID files" << endl;
  				cout << "-r RHO, --rho-prob=RHO		set $pho$ probability (default 0.95)" << endl;
				cout << "-g GAMMA, --gamma		set $gamma$ the prior of a SNP being causal (default 0.01)" << endl;
				cout << "-c causal			set the maximum number of causal SNPs" << endl;
				cout << "-f 1				to out the probaility of different number of causal SNP" << endl;
                                cout << "-w Weight file, --weight       set the biological annotation to use" << endl;
                                cout << "-t threads to use, --nthread       set the threads to use" << endl;
				exit(0);
			case 'l':
				ldFile = string(optarg);
				break;
			case 'o':
				outputFileName = string(optarg);
				break;
			case 'z':
				zFile = string(optarg);
				break;
			case 'r':
				rho = atof(optarg);
				break;
			case 'c':
				totalCausalSNP = atoi(optarg);
				break;
			case 'g':
				gamma = atof(optarg);
				break;
                        case 'w':
                                weight=string(optarg);
                                break;
                        case 't':
                                nthread=atoi(optarg);
			case 'f':
                                histFlag = true;
                                break;
			case ':':
			case '?':
			default:
				cout << "Strange" << endl;
				break;
		}
	}
	
	//program is running
	cout << "@-------------------------------------------------------------@" << endl;
	cout << "| CAVIAR!		| 	   v1.0         |  28/Sep/2016 | " << endl;
	cout << "|-------------------------------------------------------------|" << endl;
	cout << "| (C) 2016 Farhad Hormozdiari, GNU General Public License, v2 |" << endl;
	cout << "|-------------------------------------------------------------|" << endl;
	cout << "| For documentation, citation & bug-report instructions:      |" << endl;
	cout << "| 		http://genetics.cs.ucla.edu/caviar/            |" << endl;
	cout << "@-------------------------------------------------------------@" << endl;	
        cout <<"Weight file is: "<<weight<<endl;
	CaviarModel caviar(ldFile, zFile, outputFileName, totalCausalSNP, NCP, rho, histFlag, gamma,weight, nthread);
	caviar.run();
	caviar.finishUp();		
	return 0;
}
